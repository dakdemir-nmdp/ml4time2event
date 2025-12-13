# Tests for Statistical Correctness of Metrics
# Focus: IPCW Brier Score and Competing Risks Logic

test_that("BrierScore applies IPCW weights correctly", {
    # Scenario: 3 subjects, 1 censored early
    # Data:
    # 1. T=1, Event=1
    # 2. T=1.5, Event=0 (Censored)
    # 3. T=3, Event=1
    # Eval time t=2

    # Censoring Distribution (Reverse KM):
    # Obs 1: t=1, cens=0 (censored for C)
    # Obs 2: t=1.5, cens=1 (event for C)
    # Obs 3: t=3, cens=0 (censored for C)

    # KM for C:
    # t=1: n=3, d=0. S=1.
    # t=1.5: n=2, d=1. S=1 * (1 - 1/2) = 0.5.
    # t=3: n=1, d=0. S=0.5.

    # Weights for BS at t=2:
    # Obs 1: T=1 <= 2, Event=1. Weight = 1/G(1) = 1/1 = 1.
    # Obs 2: T=1.5 <= 2, Event=0. Weight = 0 (Censored before t).
    # Obs 3: T=3 > 2 (Survivor). Weight = 1/G(2) = 1/0.5 = 2.

    obstimes <- c(1, 1.5, 3)
    obsevents <- c(1, 0, 1)
    pred_times <- c(2)
    # Predictions at t=2 for the 3 subjects
    pred_surv <- matrix(c(0.2, 0.5, 0.8), nrow = 1, ncol = 3) # 1 time x 3 obs

    # Expected BS:
    # Obs 1: (1 - 0.2)^2 = 0.64. Weight 1. Contrib 0.64.
    # Obs 2: Weight 0.
    # Obs 3: Survivor. Target 1 (Surv prob). Y=1.
    # Wait, BS formula: (Y - S)^2.
    # If Y=1 (Survived), score (1 - 0.8)^2 = 0.04.
    # Weight 2. Contrib 2 * 0.04 = 0.08.
    # Total sum: 0.72.
    # Mean: 0.72 / 3 = 0.24.

    # Run function
    bs <- BrierScore(
        predsurv = pred_surv,
        pred_times = pred_times,
        obstimes = obstimes,
        obsevents = obsevents,
        eval_times = 2
    )

    expect_equal(bs$AppErr$model, 0.04, tolerance = 1e-4)
})

test_that("BrierScoreCR includes competing events in risk set", {
    # Scenario: 3 subjects
    # 1. T=1, Status=1 (Event 1)
    # 2. T=2, Status=2 (Competing Event)
    # 3. T=3, Status=0 (Censored)
    # Eval time t=2.5

    # Censoring Distribution:
    # All events (1,2) are "censored" for C.
    # Only 3 is censored at 3 (Event for C).
    # KM for C: G(t)=1 for t < 3.

    obstimes <- c(1, 2, 3)
    obsevents <- c(1, 2, 0)
    pred_times <- c(2.5)
    # Predictions (CIF for cause 1) at t=2.5
    pred_cif <- matrix(c(0.8, 0.4, 0.2), nrow = 1, ncol = 3)

    # SurvObj
    surv_obj <- survival::Surv(obstimes, obsevents, type = "mstate")

    # Weights:
    # Obs 1: T=1 <= 2.5, Status=1 (Event of Interest). Weight 1/1 = 1.
    # Obs 2: T=2 <= 2.5, Status=2 (Competing). Weight 1/1 = 1.
    # Obs 3: T=3 > 2.5. Weight 1/1 = 1.

    # Targets:
    # Obs 1: Target 1. SqError (1 - 0.8)^2 = 0.04.
    # Obs 2: Target 0 (Competing). SqError (0 - 0.4)^2 = 0.16. (Crucial check)
    # Obs 3: Target 0 (Survivor). SqError (0 - 0.2)^2 = 0.04.

    # BS = (0.04 + 0.16 + 0.04)/3 = 0.24 / 3 = 0.08.

    # Run
    bs_cr_val <- ml4time2event:::BrierScoreCR(SurvObj = surv_obj, Predictions = pred_cif, eval_times = 2.5, cause = 1, pred_times = pred_times)

    expect_equal(bs_cr_val, 0.08, tolerance = 1e-4)
})
