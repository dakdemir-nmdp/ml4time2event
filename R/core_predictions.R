#' Constructors for tidy prediction frames
#'
#' These helpers standardise the column structure returned by model predictions
#' so that the public API can compose multiple engines without special-casing.
#'
#' @keywords internal
new_survival_prediction <- function(id,
                                    time,
                                    surv,
                                    model,
                                    ensemble = FALSE,
                                    set = "train") {
  validate_prediction_lengths(id, time, surv, model, ensemble, set)
  df <- dplyr::tibble(
    id = id,
    time = time,
    surv = surv,
    model = model,
    ensemble = rep_len(ensemble, length(id)),
    set = rep_len(set, length(id))
  )
  class(df) <- c("t2e_pred_surv", "t2e_pred", class(df))
  df
}

#' @keywords internal
new_risk_prediction <- function(id,
                                risk,
                                model,
                                time = NA_real_,
                                ensemble = FALSE,
                                set = "train") {
  if (length(time) == 1L) {
    time <- rep_len(time, length(id))
  }
  validate_prediction_lengths(id, time, risk, model, ensemble, set)
  df <- dplyr::tibble(
    id = id,
    time = time,
    risk = risk,
    model = model,
    ensemble = rep_len(ensemble, length(id)),
    set = rep_len(set, length(id))
  )
  class(df) <- c("t2e_pred_risk", "t2e_pred", class(df))
  df
}

validate_prediction_lengths <- function(id, time, value, model, ensemble, set) {
  n <- length(id)
  if (!all(lengths(list(time, value, model)) %in% c(1L, n))) {
    rlang::abort("All prediction components must have length 1 or match the number of ids.")
  }
  if (length(ensemble) > 1 && length(ensemble) != n) {
    rlang::abort("`ensemble` must be scalar or length `length(id)`.")
  }
  if (length(set) > 1 && length(set) != n) {
    rlang::abort("`set` must be scalar or length `length(id)`.")
  }
}

#' @keywords internal
new_cif_prediction <- function(id,
                               time,
                               cause,
                               cif,
                               model,
                               ensemble = FALSE,
                               set = "train") {
  n <- length(id)
  if (!all(lengths(list(time, cause, cif, model)) %in% c(1L, n))) {
    rlang::abort("All prediction components must have length 1 or match the number of ids.")
  }
  if (length(ensemble) > 1 && length(ensemble) != n) {
    rlang::abort("`ensemble` must be scalar or length `length(id)`.")
  }
  if (length(set) > 1 && length(set) != n) {
    rlang::abort("`set` must be scalar or length `length(id)`.")
  }
  df <- dplyr::tibble(
    id = id,
    time = time,
    cause = cause,
    cif = cif,
    model = model,
    ensemble = rep_len(ensemble, n),
    set = rep_len(set, n)
  )
  class(df) <- c("t2e_pred_cif", "t2e_pred", class(df))
  df
}
