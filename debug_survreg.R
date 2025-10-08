devtools::load_all()
set.seed(42)

test_data2 <- rbind(
   data.frame(time = rexp(5, 0.2), event = rep(1, 5), x1 = 2, x2 = 1),
   data.frame(time = rexp(5, 0.05), event = rep(1, 5), x1 = -2, x2 = -1)
)

model_sr <- CRModel_SurvReg(test_data2, expvars=c('x1', 'x2'), timevar='time', eventvar='event', verbose=FALSE)

# Check baseline model
cat('Baseline model stored?', !is.null(model_sr$survreg_model$baseline_model), '\n')

# Test prediction step by step  
test_pred_data <- data.frame(x1 = c(2, -2), x2 = c(1, -1))
linear_preds <- predict(model_sr$survreg_model, newdata=test_pred_data, type='linear')
cat('Linear predictors:\n')
print(linear_preds)

# Test survfit
sf <- survival::survfit(
  model_sr$survreg_model$baseline_model,
  newdata = data.frame('score' = linear_preds),
  conf.int = 0.95
)

cat('Survival matrix dimensions:', dim(sf$surv), '\n')
if(!is.null(sf$surv)) {
  cat('Last survival values:\n')
  last_row <- nrow(sf$surv)
  print(sf$surv[last_row, ])
}

# Full prediction using our function
pred_sr <- Predict_CRModel_SurvReg(model_sr, test_pred_data)
cat('Our function CIF results:\n')
print(pred_sr$CIFs[nrow(pred_sr$CIFs), ])