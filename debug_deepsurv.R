devtools::load_all()
set.seed(42)

# Test DeepSurv model
test_data3 <- rbind(
   data.frame(time = rexp(20, 0.2), event = rep(1, 20), x1 = 2, x2 = 1),
   data.frame(time = rexp(20, 0.05), event = rep(1, 20), x1 = -2, x2 = -1)
)

cat('Training DeepSurv model...\n')
model_ds <- CRModel_DeepSurv(test_data3, expvars=c('x1', 'x2'), timevar='time', eventvar='event', 
                             size=5, maxit=50, verbose=FALSE)

# Test prediction
test_pred_data <- data.frame(x1 = c(2, -2), x2 = c(1, -1))
pred_ds <- Predict_CRModel_DeepSurv(model_ds, test_pred_data)

cat('DeepSurv CIF results:\n')
cat('High risk patient (x1=2):', pred_ds$CIFs[nrow(pred_ds$CIFs), 1], '\n')
cat('Low risk patient (x1=-2):', pred_ds$CIFs[nrow(pred_ds$CIFs), 2], '\n')

# Check the risk scores directly
newdata_prepared <- test_pred_data
newdata_prepared[, c('x1', 'x2')] <- sweep(newdata_prepared[, c('x1', 'x2'), drop=FALSE], 
                                          2, model_ds$model$scaling_params$mean, "-")
newdata_prepared[, c('x1', 'x2')] <- sweep(newdata_prepared[, c('x1', 'x2'), drop=FALSE], 
                                          2, model_ds$model$scaling_params$sd, "/")

x_new <- stats::model.matrix(~ . - 1, data = newdata_prepared)
log_risk_scores <- forward_pass(x_new, model_ds$model$weights)$output

cat('DeepSurv log risk scores:\n')
print(log_risk_scores)