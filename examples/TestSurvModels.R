
#

library(ml4time2event)
library(tidyverse)
library(survival)


dataTrain<-read.csv("examples/TrainingDataFollic.csv", row.names = 1)
dataTest<-read.csv("examples/TestDataFollic.csv", row.names = 1)
dataTrain$factvar<-factor(sample(c("A","B","C","D","E","F","G","H","K"), nrow(dataTrain), replace=TRUE), levels=c("A","B","C","D","E","F","G","H","K"))
dataTest$factvar<-factor(sample(c("A","B","C","D","E","F","G","H","K"), nrow(dataTest), replace=TRUE), levels=c("A","B","C","D","E","F","G","H","K"))

dataTrain$numvar<-rnorm(nrow(dataTrain))
dataTest$numvar<-rnorm(nrow(dataTest))

expvars<-c("age","ch" ,"clinstg","hgb","factvar","numvar")
#expvars<-c("age","ch" ,"clinstg","hgb")

time2eventVars<-c("status","time")
dataTrain$status[dataTrain$status==2]<-0 # converting competing risks to survival
dataTest$status[dataTest$status==2]<-0 # converting competing risks to survival
hist(dataTrain$time)
dataTrain$time<-dataTrain$time*2
dataTest$time<-dataTest$time*2
VariableProfile(dataTrain, expvars)

dataTrain<-dataTrain%>%RemoveMonoVarsData()%>% droplevelsoffactorsData()
dataTrain<-na.omit(dataTrain)
dataTest<-MakeTestDataConfWithTrainData(dataTrain,dataTest)

sum(is.na(dataTest))
dataTest<-na.omit(dataTest)

## deepsurv
deepsurvM<-SurvModel_deepsurv(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars, lr=.1, dropout = 0.1, num_nodes = c(128,64,32), weight_decay = 1e-2, epochs=100L)
PreddeepsurvM<-Predict_SurvModel_deepsurv(modelout=deepsurvM,newdata=dataTest)
str(PreddeepsurvM)
PreddeepsurvM$Probs[1:5,1:5]
matplot(PreddeepsurvM$Times,PreddeepsurvM$Probs, type="l")

# save deepsurv model using savedeepsurv function
savedeepsurv(deepsurvM, "deepsurvmodel")
## load deepsurv model using loaddeepsurv function, we need to try this after restarting R
deepsurvM<-loaddeepsurv("deepsurvmodel")
#use deepsurv model to predict
PreddeepsurvM<-Predict_SurvModel_deepsurv(modelout=deepsurvM,newdata=dataTest)

length(PreddeepsurvM$Times)
dim(PreddeepsurvM$Probs)

PreddeepsurvM$Times[1:2]
PreddeepsurvM$Probs[1:2,1:5]



## survReg
survregM<-SurvModel_SurvReg(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars)
PredsurvregM<-Predict_SurvModel_SurvReg(modelout=survregM,newdata=dataTest)
str(PredsurvregM)
matplot(PredsurvregM$Times,PredsurvregM$Probs, type="l")

cor(PredsurvregM$Probs[-1,10],PreddeepsurvM$Probs[,10])

## survReg
survregM<-SurvModel_SurvReg(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars, dist="weibull")
PredsurvregM<-Predict_SurvModel_SurvReg(modelout=survregM,newdata=dataTest )
str(PredsurvregM)
matplot(PredsurvregM$Times,PredsurvregM$Probs, type="l")




#######Rulefit
RFM<-SurvModel_rulefit(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars, ntree = 500)
PredRFM<-Predict_SurvModel_rulefit(RFM,dataTest)
matplot(PredRFM$Times,PredRFM$Probs, type="l")

length(PredRFM$Times)
dim(PredRFM$Probs)

PredRFM$Times[1:2]
PredRFM$Probs[1:2,1:5]
#######GAM
colnames(dataTrain)
GAMM<-SurvModel_GAM(data=dataTrain, eventvar=time2eventVars[1],
                    timevar=time2eventVars[2], expvars=expvars,
           shrinkTreshold = 10)

PredGAMM<-Predict_SurvModel_GAM(modelout=GAMM, newdata=dataTest)

str(PredGAMM)

matplot(PredGAMM$Times,PredGAMM$Probs, type="l")

length(PredGAMM$Times)
dim(PredGAMM$Probs)

PredGAMM$Times[1:2]
PredGAMM$Probs[1:2,1:5]

####Random Forest
RForestM<-SurvModel_RF(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars, ntree=30, samplesize=300, nsplit=5, trace=TRUE)
PredRForestM<-Predict_SurvModel_RF(modelout=RForestM,newdata=dataTest )
matplot(PredRForestM$Times,PredRForestM$Probs, type="l")

length(PredRForestM$Times)
dim(PredRForestM$Probs)

PredRForestM$Times[1:2]
PredRForestM$Probs[1:2,1:5]

####Cox
str(dataTrain)
CoxM<-SurvModel_Cox(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars)
PredCoxM<-Predict_SurvModel_Cox(modelout=CoxM,newdata=dataTest )
matplot(PredCoxM$Times,PredCoxM$Probs, type="l")

length(PredCoxM$Times)
dim(PredCoxM$Probs)

PredCoxM$Times[1:2]
PredCoxM$Probs[1:2,1:5]

####glmnet
str(dataTrain)
CoxglmnetM<-SurvModel_glmnet(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars)
PredCoxglmnetM<-Predict_SurvModel_glmnet(modelout=CoxglmnetM,newdata=dataTest )
matplot(PredCoxglmnetM$Times,PredCoxglmnetM$Probs, type="l")
length(PredCoxglmnetM$Times)
dim(PredCoxglmnetM$Probs)

PredCoxglmnetM$Times[1:2]
PredCoxglmnetM$Probs[1:2,1:5]


## xgboost
xgboostM<-SurvModel_xgboost(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars)
PredxgboostM<-Predict_SurvModel_xgboost(modelout=xgboostM,newdata=dataTest )
matplot(PredxgboostM$Times,PredxgboostM$Probs, type="l")
length(PredxgboostM$Times)
dim(PredxgboostM$Probs)

PredxgboostM$Times[1:2]
PredxgboostM$Probs[1:2,1:5]


## gbm
gbmM<-SurvModel_gbm(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars)
PredgbmM<-Predict_SurvModel_gbm(modelout=gbmM,newdata=dataTest )
matplot(PredgbmM$Times,PredgbmM$Probs, type="l")
length(PredgbmM$Times)
dim(PredgbmM$Probs)

PredgbmM$Times[1:2]
PredgbmM$Probs[1:2,1:5]



## bart
bartM<-SurvModel_BART(data=dataTrain[sample(nrow(dataTrain),10000, replace=TRUE),], eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars, ntree=5)
PredbartM<-Predict_SurvModel_BART(modelout=bartM,newdata=dataTest)
str(PredbartM)
matplot(PredbartM$Times,PredbartM$Probs, type="l")

### all together

ALLsurvregM<-RunSurvModels(datatrain=dataTrain, ExpVars=expvars, timevar=time2eventVars[2], eventvar=time2eventVars[1], models=c("glmnet","coxph","rulefit","xgboost","deepsurv","gam","gbm","ExpSurvReg","WeibSurvReg","bart"), ntreeRF=5, num_nodes = c(10L,20L,10L))
PredALLsurvregM<-PredictSurvModels(models=ALLsurvregM, newdata=dataTest, newtimes=seq(0,30,length=1000))
dim(PredALLsurvregM$NewProbs)
matplot(c(seq(0,30,length=1000)),PredALLsurvregM$NewProbs, type="l")
PredALLsurvregM$NewProbs[1:5,1:5]

matplot(c(seq(0,30,length=1000)),PredALLsurvregM$ModelPredictions$newprobsrulefit, type="l")
matplot(c(seq(0,30,length=1000)),PredALLsurvregM$ModelPredictions$newprobsgam_Model, type="l")
matplot(c(seq(0,30,length=1000)),PredALLsurvregM$ModelPredictions$newprobsrulefit, type="l")


PredALLsurvregM$ModelPredictions$newprobsrulefit[1:5,1:5]
PredALLsurvregM$ModelPredictions$newprobsglmnet[1:5,1:5]
PredALLsurvregM$ModelPredictions$newprobscph[1:5,1:5]
PredALLsurvregM$ModelPredictions$newprobsxgboost[1:5,1:5]
PredALLsurvregM$ModelPredictions$newprobsdeepsurv_Model[1:5,1:5]

sum(is.na(PredALLsurvregM$NewProbs))


score<-proba2score(PredALLsurvregM$ModelPredictions$newprobsrulefit, Times=c(seq(0,30,length=1000)), ll=0, ul=12, scale=FALSE)
conc<-SurvMetrics::Cindex(survival::Surv(dataTest$time,dataTest$status), score)

conc

##########



#################################
X=scale(model.matrix(~-1+., dataTrain[,expvars]))
times= dataTrain$time
events=dataTrain$status

# install and load the necessary packages
#install.packages(c("keras"))
library(keras)
library(tensorflow)

# initialize the model
model <- keras_model_sequential()

# add layers to the model
model %>%
  layer_dense(units = 32, activation = "relu", input_shape = ncol(X)) %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 1, activation = "linear")

coxph_loss <- function(y_true, y_pred) {
  # Calculate the risk scores from the output of the neural network
  scores <- y_pred

  # Separate the observed times and event indicators
  times <- y_true[,1]
  events <- tf$cast(y_true[,2], 'float32')  # Convert events to float32

  # Sort by time in descending order
  sorted_indices <- tf$argsort(times, direction = 'DESCENDING')
  sorted_times <- tf$gather(times, sorted_indices)
  sorted_events <- tf$gather(events, sorted_indices)
  sorted_scores <- tf$gather(scores, sorted_indices)

  # Calculate the negative partial log likelihood
  risk_set <- tf$math$exp(sorted_scores)
  log_risk_set <- tf$math$log(tf$math$cumsum(risk_set))
  event_risk_scores <- tf$math$reduce_sum(tf$boolean_mask(sorted_scores, sorted_events))
  event_log_risk_set <- tf$math$reduce_sum(tf$boolean_mask(log_risk_set, sorted_events))
  neg_log_likelihood <- event_risk_scores - event_log_risk_set

  return(neg_log_likelihood)
}


# Compile the model with gradient clipping and a lower learning rate
model %>% compile(
  optimizer = optimizer_adam(lr = 0.0001, clipvalue = .1),
  loss = coxph_loss
)

# Train the model
history <- model %>% fit(
  x = X,
  y = cbind(times, events),
  batch_size = nrow(X),
  epochs = 300,
  validation_split = 0.0
)

datasurv=data.frame(times= dataTest$time, events=dataTest$status)
Xtest=scale(model.matrix(~-1+., dataTest[,expvars]))

pred<- model$predict(Xtest)
probs<-score2proba(datasurv, as.matrix(pred))
dim(probs$sf$surv)
prob<-probs$sf$surv
time<-probs$sf$time

matplot(time,prob, type="l")


baseline_hazard <- function(scores, times, events) {
  # sort the data by time
  sorted_indices <- order(times)
  sorted_scores <- scores[sorted_indices]
  sorted_times <- times[sorted_indices]
  sorted_events <- events[sorted_indices]

  # calculate the baseline hazard function
  baseline_hazard <- c()
  for (t in unique(sorted_times)) {
    # calculate the risk set
    risk_set <- which(sorted_times >= t)

    # calculate the number of events at time t
    dN <- sum(sorted_events[risk_set])

    # calculate the sum of the risk scores for the risk set
    risk_set_sum <- sum(exp(sorted_scores[risk_set]))

    # calculate the baseline hazard at time t
    baseline_hazard[t] <- dN / risk_set_sum
  }

  return(baseline_hazard)
}

predX<-model$predict(X)
bh<-baseline_hazard(scores=predX, times=dataTrain$time+.001*rnorm(nrow(dataTrain)), events=dataTrain$status)

survival_probabilities <- function(baseline_hazard, scores, times) {
  # calculate the cumulative baseline hazard function
  cumulative_baseline_hazard <- cumsum(baseline_hazard)

  # calculate the survival probabilities
  survival_prob <- exp(-cumulative_baseline_hazard * exp(scores))

  return(survival_prob)
}

survival_probabilities(bh,predX,0.1)
