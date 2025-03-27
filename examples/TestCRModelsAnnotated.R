#install.packages("/Users/dakdemir/Library/CloudStorage/OneDrive-BeTheMatch/Documents2022_June/ml4time2event_0.1.0.tar.gz", type="source", repos=NULL)
# Load the packages 
library(ml4time2event)
library(tidyverse)
setwd("~/Library/CloudStorage/OneDrive-BeTheMatch/Documents2022_June/ml4time2event/examples/usage_examples")


dataTrain<-read.csv("TrainingDataFollic.csv", row.names = 1)
dataTest<-read.csv("TestDataFollic.csv",row.names = 1)

# Data preparation
dataTrain$factvar<-factor(sample(c("A","B","C","D","E","F","G", "H","K","L"), nrow(dataTrain), replace=TRUE), levels=c("A","B","C","D","E","F","G", "H","K","L"))
dataTest$factvar<-factor(sample(c("A","B","C","D","E","F","G", "H","K","L"), nrow(dataTest), replace=TRUE), levels=c("A","B","C","D","E","F","G", "H","K","L"))
str(dataTrain)
dataTrain$numvar<-rnorm(nrow(dataTrain))
dataTest$numvar<-rnorm(nrow(dataTest))

for (vari in colnames(dataTrain)){
  x<-dataTrain[, vari]
  if (is.character(x)){
    dataTrain[, vari]<-as.factor(dataTrain[, vari])
    dataTest[, vari]<-factor(dataTest[, vari], levels=levels(dataTrain[, vari]))
  }
}
str(dataTrain)

table(dataTrain$factvar)

# Create the model
expvars<-c("age","ch" ,"clinstg","hgb","factvar","numvar" )
#expvars<-c("age","ch" ,"clinstg","hgb")

time2eventVars<-c("status","time")
hist(dataTrain$time)
dataTrain$time<-dataTrain$time*2
dataTest$time<-dataTest$time*2
VariableProfile(dataTrain, expvars)

dataTrain<-dataTrain%>%RemoveMonoVarsData()%>% droplevelsoffactorsData()
dataTrain<-na.omit(dataTrain)
dataTest<-MakeTestDataConfWithTrainData(dataTrain,dataTest)
str(dataTest)
sum(is.na(dataTest))
dataTest<-na.omit(dataTest)

table(dataTest$factvar)
table(dataTrain$factvar)


for (repi in 1:10){
## rulefit
CRModel_rulefitModel<-CRModel_rulefit(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2],
                                 expvars=expvars, ntree=50)
PredCRModel_rulefit<-Predict_CRModel_rulefit(modelout=CRModel_rulefitModel,newdata=dataTest)
matplot(PredCRModel_rulefit$Times,t(PredCRModel_rulefit$CIFs)
        , type="l")
}

#######RF
RandForestM<-CRModel_RF(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars, ntree=10)
PredRandForestM<-Predict_CRModel_RF(RandForestM,dataTest)
str(PredRandForestM)
max(PredRandForestM$CIFs)


####Cox
CoxM<-CRModel_Cox(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars)
PredCoxM<-Predict_CRModel_Cox(modelout=CoxM,newdata=dataTest )
str(PredCoxM)

matplot(PredCoxM$Times,t(PredCoxM$CIFs)
        , type="l")
####FG
str(dataTrain)
CRModel_FineGrayM<-CRModel_FineGray(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars)
PredCRModel_FineGray<-Predict_CRModel_FG(modelout=CRModel_FineGrayM,newdata=dataTest )
str(PredCRModel_FineGray)

matplot(PredCRModel_FineGray$Times,t(PredCRModel_FineGray$CIFs)
        , type="l")

## BART
CRModel_BARTM<-CRModel_BART(data=dataTrain, eventvar=time2eventVars[1], timevar=time2eventVars[2], expvars=expvars, ntree=10)
PredCRModel_BARTM<-Predict_CRModel_BART(modelout=CRModel_BARTM,newdata=dataTest )
str(PredCRModel_BARTM)
matplot(PredCRModel_BARTM$Times,t(PredCRModel_BARTM$CIFs)
        , type="l")



max(PredCRModel_rulefit$CIFs)

### all together
ALLCRModelsM<-RunCRModels(datatrain=dataTrain, ExpVars=expvars, timevar=time2eventVars[2], eventvar=time2eventVars[1], models=c("bart","FG","cox", "rulefit"), ntreeRF=50)
#ALLCRModelsM<-RunCRModels(datatrain=dataTrain, ExpVars=expvars, timevar=time2eventVars[2], eventvar=time2eventVars[1], models=c("rulefit"), ntreeRF=300)

PredALLALLCRModelsM<-PredictCRModels(models=ALLCRModelsM, newdata=dataTest, newtimes=seq(0,30,length=1000))
dim(PredALLALLCRModelsM$NewProbs)
matplot(c(seq(0,30,length=1000)),PredALLALLCRModelsM$NewProbs, type="l")

##########################





PredALLALLCRModelsM$NewProbs[1:5,1:5]
PredALLALLCRModelsM$ModelPredictions$newprobsRF[1:5,1:5]
PredALLALLCRModelsM$ModelPredictions$newprobsRF2[1:5,1:5]
PredALLALLCRModelsM$ModelPredictions$newprobsFG[1:5,1:5]
PredALLALLCRModelsM$ModelPredictions$newprobsBART[1:5,1:5]
PredALLALLCRModelsM$ModelPredictions$newprobsCox[1:5,1:5]




################################################################



matplot(c(seq(0,30,length=1000)),PredALLALLCRModelsM$ModelPredictions$newprobsBART
, type="l")
matplot(c(seq(0,30,length=1000)),PredALLALLCRModelsM$ModelPredictions$newprobsCox
        , type="l")
matplot(c(seq(0,30,length=1000)),PredALLALLCRModelsM$ModelPredictions$newprobsFG
        , type="l")

matplot(c(seq(0,30,length=1000)),PredALLALLCRModelsM$NewProbs
        , type="l")
sum(is.na(PredALLALLCRModelsM$NewProbs))
hist(proba2score(PredALLALLCRModelsM$ModelPredictions$newprobsBART, Times=c(seq(0,30,length=1000)), ll=0, ul=36, scale=TRUE))

###############


