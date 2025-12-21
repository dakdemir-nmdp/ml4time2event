#' RuleFit Survival Model (R6)
#'
#' Fits a RuleFit model for survival analysis.
#' Generates rules from an ensemble of trees (rpart) and fits a penalized Cox model
#' (glmnet) on the original features plus the generated rules.
#'
#' @keywords internal
#' @noRd
RulefitSurvivalModel <- R6::R6Class(
  classname = "RulefitSurvivalModel",
  inherit = SurvivalModel,
  public = list(
    model = NULL,
    time_grid = NULL,
    task = NULL,
    varprof = NULL,
    factor_levels = NULL,
    train_matrix = NULL, # Required for prediction (rule matrix structure)
    y_train = NULL, # Required for baseline hazard
    cv_fit = NULL, # The glmnet fit
    rules_list = NULL,
    trees_list = NULL,
    sample_data = NULL, # Keep sample for factor levels consistency

    fit = function(task, time_grid, ...) {
      super$fit(task = task, ...)
      self$task <- task
      self$time_grid <- time_grid

      data <- as.data.frame(task$data)
      expvars <- task$features
      timevar <- task$time_col
      eventvar <- task$event_col

      # Params
      spec <- self$spec
      ntree <- spec$ntree %||% 30
      nsample <- spec$nsample %||% min(300, nrow(data))
      keepvars <- spec$keepvars
      alpha <- spec$alpha %||% 0.5
      maxit <- spec$maxit %||% 2000

      self$varprof <- VariableProfile(data, expvars)

      # Factor levels
      self$factor_levels <- lapply(data[, expvars, drop = FALSE], function(x) {
        if (is.factor(x)) levels(x) else NULL
      })

      # Cuttimes for classification trees
      x_time <- data[[timevar]]
      cuttimes <- stats::quantile(x_time, c(.1, .25, .50, .70, .90), na.rm = TRUE)

      # Formulas
      form_surv <- stats::as.formula(paste("survival::Surv(", timevar, ",", eventvar, ") ~ ."))
      form_class <- stats::as.formula("ClassVar ~ .")

      # 2. Generate Ensemble
      trees <- list()
      for (i in 1:ntree) {
        # Bagging
        # Feature sampling
        n_vars <- length(expvars)
        n_sample_vars <- min(n_vars, sample(1:min(10, n_vars), 1))
        sampcols <- union(keepvars, sample(expvars, n_sample_vars))

        # Row sampling
        samprows <- sample(seq_len(nrow(data)), nsample, replace = TRUE)

        # Model type
        selmodel <- sample(c(1, 1, 1, 1, 2), 1) # Favours survival

        usevars <- c(timevar, eventvar, sampcols)
        datasampl <- data[samprows, usevars, drop = FALSE]

        form <- if (selmodel == 1) form_surv else form_class

        if (selmodel == 2) {
          # Classification tree (event at cut time)
          cut_val <- sample(cuttimes, 1)
          datasampl$ClassVar <- as.factor(
            datasampl[[timevar]] < cut_val & datasampl[[eventvar]] == 1
          )
          datasampl <- datasampl[, !colnames(datasampl) %in% c(timevar, eventvar), drop = FALSE]
        }

        # Fit rpart
        ctrl <- rpart::rpart.control(
          minsplit = stats::rpois(1, 2) + 1,
          minbucket = stats::rpois(1, 20) + 5,
          cp = 0.01 * stats::runif(1),
          maxdepth = stats::rpois(1, 2) + 1,
          usesurrogate = 0,
          xval = 0
        )

        # Capture and store
        fit_tree <- tryCatch(
          {
            rpart::rpart(form, data = datasampl, control = ctrl, model = TRUE)
          },
          error = function(e) NULL
        )

        if (!is.null(fit_tree)) trees[[i]] <- fit_tree
      }
      # Clean nulls
      trees <- trees[!sapply(trees, is.null)]
      self$trees_list <- trees

      # 3. Extract Rules
      rules_list <- lapply(trees, function(x) {
        tryCatch(.extract_rules_from_party(partykit::as.party(x)), error = function(e) list())
      })
      self$rules_list <- rules_list

      # Create Rule Matrix
      rules_mats <- list()
      for (i in seq_along(trees)) {
        tree <- trees[[i]]
        rules <- rules_list[[i]]

        if (length(rules) == 0 || length(unique(names(rules))) <= 1) next

        party_tree <- partykit::as.party(tree)
        # Predict node IDs
        node_preds <- predict(party_tree, newdata = data, type = "node")

        # Create factors with all possible rule nodes as levels
        train_factor <- factor(as.character(node_preds), levels = names(rules))

        mm <- stats::model.matrix(~ -1 + train_factor)
        # Naming convention: i_ruleIdx***NodeID
        colnames(mm) <- paste(paste(i, seq_along(rules), sep = "_"), rules, sep = "***")
        rules_mats[[length(rules_mats) + 1]] <- mm
      }

      # 4. Combine
      train_mat_orig <- stats::model.matrix(~ . - 1, data = data[, expvars, drop = FALSE])

      if (length(rules_mats) > 0) {
        train_mat_rules <- Reduce("cbind", rules_mats)
        final_mat <- cbind(train_mat_orig, train_mat_rules)
      } else {
        final_mat <- train_mat_orig
      }

      # 5. Fit GLMNet
      y_train <- survival::Surv(data[[timevar]], data[[eventvar]] == 1)

      cv_fit <- glmnet::cv.glmnet(
        x = final_mat,
        y = y_train,
        alpha = alpha,
        family = "cox",
        maxit = maxit
      )

      self$cv_fit <- cv_fit
      self$y_train <- y_train
      self$train_matrix <- final_mat
      # Store sample for levels alignment
      self$sample_data <- data[1:min(nrow(data), 10), expvars, drop = FALSE]

      invisible(self)
    },
    predict_survival = function(newdata, times, set = "test", ...) {
      private$ensure_fitted()
      complete_data <- .ensure_prediction_data(newdata, self$task)
      complete_idx <- stats::complete.cases(complete_data[, self$task$features])
      id_values <- complete_data[[self$task$id_col]]

      target_times <- if (is.null(times)) self$time_grid else sort(unique(as.numeric(times)))

      # We need to construct the matrix for newdata
      # 1. Transform factors
      data_test <- complete_data[, self$task$features, drop = FALSE]
      for (v in self$task$features) {
        if (!is.null(self$factor_levels[[v]])) {
          data_test[[v]] <- factor(as.character(data_test[[v]]), levels = self$factor_levels[[v]])
        }
      }

      # 2. Rules Matrix
      rules_mats <- list()
      for (i in seq_along(self$trees_list)) {
        tree <- self$trees_list[[i]]
        rules <- self$rules_list[[i]]
        if (length(rules) == 0 || length(unique(names(rules))) <= 1) next

        node_preds <- predict(partykit::as.party(tree), newdata = data_test, type = "node")
        test_factor <- factor(as.character(node_preds), levels = names(rules))
        mm <- stats::model.matrix(~ -1 + test_factor)
        colnames(mm) <- paste(paste(i, seq_along(rules), sep = "_"), rules, sep = "***")
        rules_mats[[length(rules_mats) + 1]] <- mm
      }

      # 3. Combine
      # Original matrix (aligned via rbind trick usually, or just careful model.matrix)
      # We used model.matrix(~.-1, data[expvars]).
      # To ensure levels match EXACTLY, we can use the sample data trick or explicit factor setting
      # We already set factor levels above explicitly.

      mm_orig_test <- stats::model.matrix(~ . - 1, data = data_test)

      if (length(rules_mats) > 0) {
        test_mat_rules <- Reduce("cbind", rules_mats)
        test_mat <- cbind(mm_orig_test, test_mat_rules)
      } else {
        test_mat <- mm_orig_test
      }

      # 4. Align columns with training matrix
      train_cols <- colnames(self$train_matrix)
      test_cols <- colnames(test_mat)

      missing <- setdiff(train_cols, test_cols)
      if (length(missing) > 0) {
        # Add missing as 0
        zeros <- matrix(0, nrow = nrow(test_mat), ncol = length(missing))
        colnames(zeros) <- missing
        test_mat <- cbind(test_mat, zeros)
      }

      # Reorder and subset
      test_mat <- test_mat[, train_cols, drop = FALSE]

      # 5. Predict
      # Use survfit to get survival curve using training data for baseline
      sf <- survival::survfit(
        self$cv_fit$glmnet.fit,
        s = self$cv_fit$lambda.min,
        x = self$train_matrix,
        y = self$y_train,
        newx = test_mat
      )

      base_times <- sf$time
      surv_probs <- sf$surv # [times, obs] often? Or [obs, times]?
      # glmnet survfit usually returns [times, n_newx]
      # Verified in previous steps: sf$surv from glmnet survfit depends on input.
      # If newx provided, cols are obs.

      if (is.null(dim(surv_probs))) {
        surv_probs <- matrix(surv_probs, ncol = 1)
      }

      # Add 0
      if (!0 %in% base_times) {
        base_times <- c(0, base_times)
        surv_probs <- rbind(rep(1, ncol(surv_probs)), surv_probs)
      }

      preds_mat <- survprobMatInterpolator(surv_probs, base_times, target_times)

      # Flatten
      flat_preds <- as.vector(preds_mat)

      new_survival_prediction(
        id = rep(id_values, each = length(target_times)),
        time = rep(target_times, times = length(id_values)),
        surv = flat_preds,
        model = rep("rulefit", length(flat_preds)),
        ensemble = FALSE,
        set = set
      )
    },
    model_info = function() {
      info <- super$model_info()
      info$label <- "RuleFit Survival"
      info
    },
    required_packages = function() {
      c("rpart", "partykit", "glmnet", "survival")
    }
  ),
  private = list(
    ensure_fitted = function() {
      if (!isTRUE(self$fitted)) rlang::abort("Model not fitted")
    }
  )
)

# Helper for extraction (internal)
.extract_rules_from_party <- function(x) {
  if (is.null(x) || is.null(x$node)) {
    return(list())
  }
  # partykit internal function - risky but standard for rule extraction
  # Check if accesible
  if (exists(".list.rules.party", where = asNamespace("partykit"), mode = "function")) {
    list_rules_party <- utils::getFromNamespace(".list.rules.party", "partykit")
    rule_strings <- list_rules_party(x)
    node_ids <- names(partykit::nodeids(x, terminal = TRUE))
    names(rule_strings) <- node_ids
    return(rule_strings)
  }
  return(list())
}

.register_time_to_event_model(
  engine = "rulefit",
  outcome = "survival",
  constructor = function(spec = list()) {
    RulefitSurvivalModel$new(spec = modifyList(list(engine = "rulefit"), spec))
  },
  packages = c("rpart", "partykit", "glmnet", "survival"),
  tags = c("rules", "tree-based"),
  label = "RuleFit survival"
)
