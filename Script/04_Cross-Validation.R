#___________________________________________________________________________________________________________________________________#
###
### 
###           Extract from Hansser et al., 2020 script in supplementary materials
###
#___________________________________________________________________________________________________________________________________#



#==============================
# INTACT GIS Validation
# Author: R. May
# Date: 09.02.2017
#=============================

library(readxl)
library(lme4)
library(MuMIn)
library(effects)
library(lattice)
#library(circular)
#library(stlplus)
library(LMERConvenienceFunctions)
library(AICcmodavg)
library(merTools)
library(lmerTest)
library(cvAUC)

auc.ci <- function(model, data, V=10){
  .cvFolds <- function(Y, V){ #Create CV folds (stratify by outcome)
    Y0 <- split(sample(which(Y==0)), rep(1:V, length=length(which(Y==0))))
    Y1 <- split(sample(which(Y==1)), rep(1:V, length=length(which(Y==1))))
    folds <- vector("list", length=V)
    for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}
    return(folds)
  }
  .doFit <- function(v, folds, data, model){ #Train/test glmer for each fold
    fit <- glm(formula(model), data=data[-folds[[v]],], family=binomial)
    pred <- predict(fit, newdata=data[folds[[v]],] ,allow.new.levels=T)
    pred <- exp(pred)/max(exp(pred), na.rm = TRUE)
    return(pred)
  }
  folds <- .cvFolds(Y=data[,"case_"], V=V) #Create folds
  predictions <- unlist(sapply(seq(V), .doFit, folds=folds, data=data, model=model)) #CV train/predict
  predictions[unlist(folds)] <- predictions #Re-order pred values
  # Get CV AUC and confidence interval
  out <- ci.cvAUC(predictions=predictions, labels=data$case_, folds=folds, confidence=0.95)
  return(out)
}


#Assess model performance:
# Done by hand iteratively (need to fit model without random effect)
source('./Functions/kxv_glmer.r')  # script from Hanssen et al., 2020 supplementary materials

kxvglm(formula(rsf_Tramq3),rsf_Tramq3$frame,k=10,nbin=10,family=binomial)
(auc <- auc.ci(model=rsf_Tramq3,data=rsf_Tramq3$frame,V=10))


