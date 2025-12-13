#!/usr/bin/env Rscript
# Test whether crisk.pre.bart() actually requires duplicated test data

library(BART)

# Create test data
set.seed(123)
n_train <- 100
n_test <- 20

# Training data
train_data <- data.frame(
    time = rexp(n_train, rate = 0.01),
    event = sample(0:2, n_train, replace = TRUE, prob = c(0.2, 0.4, 0.4)),
    x1 = rnorm(n_train),
    x2 = rnorm(n_train)
)

# Test data
test_data <- data.frame(
    x1 = rnorm(n_test),
    x2 = rnorm(n_test)
)

# Prepare training matrices
x_train <- as.matrix(train_data[, c("x1", "x2")])
times_train <- train_data$time
delta_train <- as.integer(train_data$event)

# Prepare test matrix
x_test <- as.matrix(test_data[, c("x1", "x2")])

cat("=== Test 1: WITHOUT duplication ===\n")
cat("Test data dimensions:", nrow(x_test), "x", ncol(x_test), "\n\n")

# Try WITHOUT duplication
result_no_dup <- tryCatch(
    {
        pre_no_dup <- BART::crisk.pre.bart(
            time = times_train,
            delta = delta_train,
            x.train = BART::bartModelMatrix(x_train),
            x.test = BART::bartModelMatrix(x_test),
            x.train2 = BART::bartModelMatrix(x_train),
            x.test2 = BART::bartModelMatrix(x_test),
            K = 10
        )
        cat("✓ SUCCESS without duplication\n")
        cat("tx.test dimensions:", nrow(pre_no_dup$tx.test), "x", ncol(pre_no_dup$tx.test), "\n")
        cat("tx.test2 dimensions:", nrow(pre_no_dup$tx.test2), "x", ncol(pre_no_dup$tx.test2), "\n\n")
        pre_no_dup
    },
    error = function(e) {
        cat("✗ FAILED without duplication\n")
        cat("Error:", e$message, "\n\n")
        NULL
    }
)

cat("\n=== Test 2: WITH duplication ===\n")
x_test_dup <- rbind(x_test, x_test)
cat("Duplicated test data dimensions:", nrow(x_test_dup), "x", ncol(x_test_dup), "\n\n")

# Try WITH duplication
result_with_dup <- tryCatch(
    {
        pre_with_dup <- BART::crisk.pre.bart(
            time = times_train,
            delta = delta_train,
            x.train = BART::bartModelMatrix(x_train),
            x.test = BART::bartModelMatrix(x_test_dup),
            x.train2 = BART::bartModelMatrix(x_train),
            x.test2 = BART::bartModelMatrix(x_test_dup),
            K = 10
        )
        cat("✓ SUCCESS with duplication\n")
        cat("tx.test dimensions:", nrow(pre_with_dup$tx.test), "x", ncol(pre_with_dup$tx.test), "\n")
        cat("tx.test2 dimensions:", nrow(pre_with_dup$tx.test2), "x", ncol(pre_with_dup$tx.test2), "\n\n")
        pre_with_dup
    },
    error = function(e) {
        cat("✗ FAILED with duplication\n")
        cat("Error:", e$message, "\n\n")
        NULL
    }
)

# Now test actual prediction if preprocessing succeeded
if (!is.null(result_no_dup) || !is.null(result_with_dup)) {
    cat("\n=== Test 3: Fitting BART model ===\n")

    # Fit a BART model
    bart_model <- BART::crisk.bart(
        x.train = x_train,
        times = times_train,
        delta = delta_train,
        x.test = x_train,
        K = 10,
        ntree = 20,
        ndpost = 100,
        nskip = 50,
        keepevery = 5,
        numcut = 2
    )

    cat("Model fitted successfully\n")
    cat("Model K:", bart_model$K, "\n\n")

    # Test prediction WITHOUT duplication
    if (!is.null(result_no_dup)) {
        cat("=== Test 4a: Prediction WITHOUT duplication ===\n")
        pred_no_dup <- tryCatch(
            {
                p <- predict(bart_model, result_no_dup$tx.test, result_no_dup$tx.test2)
                cat("✓ Prediction SUCCESS\n")
                cat("CIF length:", length(p$cif.test.mean), "\n")
                cat("Expected length (N*K):", n_test * bart_model$K, "\n")
                cat("Ratio:", length(p$cif.test.mean) / (n_test * bart_model$K), "\n\n")
                p
            },
            error = function(e) {
                cat("✗ Prediction FAILED\n")
                cat("Error:", e$message, "\n\n")
                NULL
            }
        )
    }

    # Test prediction WITH duplication
    if (!is.null(result_with_dup)) {
        cat("=== Test 4b: Prediction WITH duplication ===\n")
        pred_with_dup <- tryCatch(
            {
                p <- predict(bart_model, result_with_dup$tx.test, result_with_dup$tx.test2)
                cat("✓ Prediction SUCCESS\n")
                cat("CIF length:", length(p$cif.test.mean), "\n")
                cat("Expected length (N*K):", n_test * bart_model$K, "\n")
                cat("Expected with duplication (2*N*K):", 2 * n_test * bart_model$K, "\n")
                cat("Ratio:", length(p$cif.test.mean) / (n_test * bart_model$K), "\n\n")
                p
            },
            error = function(e) {
                cat("✗ Prediction FAILED\n")
                cat("Error:", e$message, "\n\n")
                NULL
            }
        )
    }

    # Compare predictions if both succeeded
    if (!is.null(result_no_dup) && !is.null(result_with_dup) &&
        exists("pred_no_dup") && exists("pred_with_dup") &&
        !is.null(pred_no_dup) && !is.null(pred_with_dup)) {
        cat("\n=== Test 5: Comparing predictions ===\n")

        # Extract first N*K values from both
        len_no_dup <- length(pred_no_dup$cif.test.mean)
        len_with_dup <- length(pred_with_dup$cif.test.mean)
        expected_len <- n_test * bart_model$K

        cif_no_dup <- pred_no_dup$cif.test.mean[1:min(expected_len, len_no_dup)]
        cif_with_dup <- pred_with_dup$cif.test.mean[1:min(expected_len, len_with_dup)]

        if (length(cif_no_dup) == length(cif_with_dup)) {
            max_diff <- max(abs(cif_no_dup - cif_with_dup))
            cat("Maximum difference between predictions:", max_diff, "\n")

            if (max_diff < 1e-10) {
                cat("✓ Predictions are IDENTICAL\n")
                cat("CONCLUSION: Duplication is NOT necessary!\n")
            } else {
                cat("✗ Predictions DIFFER\n")
                cat("CONCLUSION: Duplication may affect results\n")
            }
        } else {
            cat("Cannot compare - different lengths\n")
        }
    }
}

cat("\n=== SUMMARY ===\n")
cat("Without duplication:", ifelse(is.null(result_no_dup), "FAILED", "SUCCESS"), "\n")
cat("With duplication:", ifelse(is.null(result_with_dup), "FAILED", "SUCCESS"), "\n")

if (!is.null(result_no_dup)) {
    cat("\n✓ TEST DATA DUPLICATION IS NOT NECESSARY\n")
    cat("The current implementation can be simplified.\n")
} else {
    cat("\n✗ TEST DATA DUPLICATION APPEARS NECESSARY\n")
    cat("The current implementation should be kept.\n")
}
