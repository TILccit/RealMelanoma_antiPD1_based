


myCindex <- function(D){
  C= 0.5 *D + 0.5
  C
}



CalculateCharrellFromDxy <- function(validate) {
  ## Test if the object is correct
  stopifnot(class(validate) == "validate")
  
  ## Calculate AUCs from Dxy's
  aucs <- (validate["Dxy", c("index.orig","training","test","optimism","index.corrected")])/2 + 0.5
  
  ## Get n
  n <- validate["Dxy", c("n")]
  
  ## Combine as result
  res <- rbind(validate, C.index = c(aucs, n))
  
  ## Fix optimism
  res["C.index","optimism"] <- res["C.index","optimism"] - 0.5
  
  ## Return results
  res
}




overallAssessment_Apparent <- function(modelList,modelNames, model,modelformula,
                                       formula, data, horizon, time, status,
                                       bootStrapping=FALSE,
                                       bootName=NULL,simple=FALSE){
  # # # # # Runs for one model at a time
  # modelList= list(fullfit)
  # modelNames=c("simple")
  # model = fullfit.rms
  # modelformula =modelFormula(surv_time, surv_status, feats1)
  # formula = modelFormula(surv_time, surv_status, 1)
  # data = data.dev
  # horizon = 35.9999
  
  # model=cox_boot[[1]]
  # modelNames=c("simple")
  # modelformula =modelFormula(time,  status, features)
  # formula=modelFormula(surv_time,  status, 1)
  # data=data
  # horizon=horizon
  # bootStrapping = T
  # bootName = "Apparent"
  # 
  ###########
  ## START 
  modelList=list(model)
  names(modelList) <-modelNames
  # # Apparent Brier and IPA)
  # scores_brier <- riskRegression::Score(modelList,
  #                                       formula = formula, 
  #                                       data = data, 
  #                                       conf.int = TRUE, 
  #                                       times = horizon,
  #                                       cens.method = "ipcw",
  #                                       cens.model = "cox",
  #                                       metrics = "brier",
  #                                       summary = "ipa"
  # )
  # 
  # scores_auc <- riskRegression::Score(modelList,
  #                                     formula = formula, 
  #                                     data = data, 
  #                                     conf.int = TRUE, 
  #                                     times = horizon,
  #                                     cens.method = "ipcw",
  #                                     cens.model = "cox",
  #                                     metrics = "AUC"
  # )
  # 
  # briers <- scores_brier %>% pluck("Brier") %>% pluck("score") %>% dplyr::filter(model %in% modelNames) %>% dplyr::select(Brier,IPA) %>% t() %>% as.data.frame()
  # 
  # auc <- scores_auc %>% pluck("AUC") %>% pluck("score") %>% dplyr::filter(model %in% modelNames) %>% dplyr::select(AUC) %>% t() %>% as.data.frame()
  # 
  ############################
  scores_brier <- riskRegression::Score(modelList,
                                        formula = formula,
                                        data = data,
                                        conf.int = TRUE,
                                        times = horizon,
                                        cens.method = "ipcw",
                                        cens.model = "cox",
                                        # metrics = "brier",
                                        summary = "ipa"
  )
  
  
  briers.i <- scores_brier %>% purrr::pluck("Brier") %>% purrr::pluck("score") %>%
    dplyr::filter(model %in% modelNames) %>%
    dplyr::select(times,Brier,IPA) #%>% as.numeric()# %>% t() %>% as.data.frame()
  
  
  auc.i <- scores_brier %>% purrr::pluck("AUC") %>% purrr::pluck("score") %>%
    dplyr::filter(model %in% modelNames) %>%
    dplyr::select(times,AUC) #%>% as.numeric() # %>% t() %>% as.data.frame()
  
  
  res<- cbind(auc.i[,2],briers.i)
  
  res<-res %>%
    group_split(times,.keep = FALSE)
  names(res) <- horizon
  ###############
  
  
  # Check if simple version that gives out AUC, Brier, IPA, C.index
  if(isTRUE(simple)){
    # # Somers D
    # Dxy <- model$stats["Dxy"][[1]]  %>% as.data.frame()  %>% setNames("V1") %>% set_rownames(c("Dyx"))
    # # C index
    # Cind <-(Dxy$V1/2) + 0.5 %>% as.data.frame() %>% setNames("V1") %>% set_rownames(c("Cindex"))
    # res <-rbind(briers,auc,Cind) %>% setNames("Original")
    # res
    # Somers D
    # Dxy <- model$stats["Dxy"][[1]]  %>% as.data.frame()  %>% setNames("Dyx")# %>% set_rownames(c("Dyx"))
    # C index
    # Cind <-(Dxy$Dyx/2) + 0.5 %>% as.data.frame() %>% setNames("Cindex") #%>% set_rownames(c("Cindex"))
    
    Cind <-concordance(Surv(data[[time]], data[[status]])
                ~ predict(model, newdata = data),
                reverse = TRUE)$concordance %>% as.data.frame() %>% setNames("Cindex")
    Dxy <- (Cind$Cindex - 0.5) *2 %>% as.data.frame()  %>% setNames("Dxy")
    
    res <-sapply(res, function(x) cbind(Dxy,Cind,x) ,simplify = FALSE)
    
    
  }else{
    ## Check if bootstrapping otherwise use other way
    if(!isTRUE(bootStrapping)){
      # Somers D
      # Dxy <- model$stats["Dxy"][[1]]  %>% as.data.frame()  %>% setNames("V1") %>% set_rownames(c("Dyx"))
      # # C index
      # Cind <-(Dxy$V1/2) + 0.5 %>% as.data.frame() %>% setNames("V1") %>% set_rownames(c("Cindex"))
      # # R2
      # R2 <- model$stats["R2"][[1]]  %>% as.data.frame()  %>% setNames("V1") %>% set_rownames(c("R2"))
      # # Gini g index
      # g <-GiniMd(predict(model)) %>% as.data.frame()  %>% setNames("V1") %>% set_rownames(c("g"))
      
      # Dxy <- model$stats["Dxy"][[1]]  %>% as.data.frame()  %>% setNames("Dyx")# %>% set_rownames(c("Dyx"))
      # # C index
      # Cind <-(Dxy$Dyx/2) + 0.5 %>% as.data.frame() %>% setNames("Cindex") #%>% set_rownames(c("Cindex"))
      
      Cind <-concordance(Surv(data[[time]], data[[status]])
                         ~ predict(model, newdata = data),
                         reverse = TRUE)$concordance %>% as.data.frame() %>% setNames("Cindex")
      Dxy <- (Cind$Cindex - 0.5) *2 %>% as.data.frame()  %>% setNames("Dxy")
      # R2
      R2 <- model$stats["R2"][[1]]  %>% as.data.frame()  %>% setNames("R2") #%>% set_rownames(c("R2"))
      # Gini g index
      g <-GiniMd(predict(model)) %>% as.data.frame()  %>% setNames("g") #%>% set_rownames(c("g"))
      # Final results
      # res <-rbind(Dxy,R2,g,briers,auc,Cind) %>% setNames("Original")
      # res
      res <-sapply(res, function(x) cbind(Dxy,Cind,R2,g,x) ,simplify = FALSE)
      
      
    }else{
      # res <-rbind(briers,auc) %>% setNames(bootName)
      # res
      res <-sapply(res, function(x) x %>% set_rownames(bootName) ,simplify = FALSE)
      
    }
  }
  
  
  # gp, g index on the probability scale
  # GiniMd(predict(model, type = fitted))
  # D discrimination index
  # U unrealiability index
  # Q overall quality
  # Cind <-cindex(modelFormula(surv_time, surv_status, feats1),data=data)  %>% pluck("cindex") %>% as.data.frame() %>% setNames("V1") %>% set_rownames(c("Cindex"))
  # # Somers D
  # Dxy <- 2*(Cind$V1-0.5) %>% as.data.frame()  %>% setNames("V1") %>% set_rownames(c("Dyx"))
  # Extract table with Brier, IPA, AUC
  
  # # Update model with new data
  # data$lp <- predict(model, newdata = data)
  # # Harrell C - development
  # harrell_C <- concordance(modelFormula(surv_time, surv_status, "lp"), 
  #                          data, 
  #                          reverse = TRUE)
  
  res
}



overallAssessment_InternalValidation <- function(modelList,modelNames, formula, 
                                                 data, horizon, nBoot=100,mySeed=123){
  #Runs for one model at a time
  # modelList= list(fullfit)
  # modelNames=c("simple")
  # formula = modelFormula(surv_time, surv_status, 1)
  # data = data.dev
  # horizon = 35.9999
  # metric="brier"
  # nBoot=100
  # mySeed=123
  
  ###########
  ## START 
  names(modelList) <-modelNames
  # set.seed(1234)
  # clust <-parallel::makeCluster(5)
  # Apparent Brier and IPA)
  scores_brier <- riskRegression::Score(modelList,
                                        formula = formula, 
                                        data = data, 
                                        conf.int = TRUE, 
                                        times = horizon,
                                        cens.method = "ipcw",
                                        cens.model = "cox",
                                        metrics = "brier",
                                        summary = "ipa",
                                        seed=mySeed,
                                        # parallel = "as.registered", 
                                        split.method="bootcv",B=nBoot,
                                        #M=nrow(data)#,
                                        #ncpus = 5#,
                                        #cl = clust
  )
  # parallel::stopCluster(clust)
  scores_auc <- riskRegression::Score(modelList,
                                      formula = formula, 
                                      data = data, 
                                      conf.int = TRUE, 
                                      times = horizon,
                                      cens.method = "ipcw",
                                      cens.model = "cox",
                                      metrics = "AUC",
                                      seed=mySeed,
                                      # parallel = "as.registered", 
                                      split.method="bootcv",B=nBoot,
                                      #M=nrow(data)
                                      #ncpus = 5#,
                                      #cl = clust
  )
  
  briers <- scores_brier %>% pluck("Brier") %>% pluck("score") %>% dplyr::filter(model %in% modelNames) %>% dplyr::select(Brier,IPA) %>% t() %>% as.data.frame()
  
  auc <- scores_auc %>% pluck("AUC") %>% pluck("score") %>% dplyr::filter(model %in% modelNames) %>% dplyr::select(AUC) %>% t() %>% as.data.frame()
  
  
  # Somers D
  Dxy <- c("")  %>% as.data.frame()  %>% setNames("V1") %>% set_rownames(c("Dyx"))
  # C index
  Cind <-c("") %>% as.data.frame() %>% setNames("V1") %>% set_rownames(c("Cindex"))
  # R2
  R2 <- c("")  %>% as.data.frame()  %>% setNames("V1") %>% set_rownames(c("R2"))
  # Gini g index
  g <-c("") %>% as.data.frame()  %>% setNames("V1") %>% set_rownames(c("g"))
  
  # Extract table with Brier, IPA, AUC
  res <-rbind(Dxy,R2,g,briers,auc,Cind) %>% setNames("Internal_validation")
}





performanceCalculation <- function(model,name,formula,data,horizon){
  
  # model
  # name
  # formula
  # data
  # horizon
  # 
  ###########
  ## START
  # Apparent Brier and IPA)
  modelList = list(model)
  names(modelList)=name
  scores_brier <- riskRegression::Score(modelList,
                                        formula = formula, 
                                        data = data, 
                                        conf.int = TRUE, 
                                        times = horizon,
                                        cens.method = "ipcw",
                                        cens.model = "cox",
                                        # metrics = "brier",
                                        summary = "ipa"
  )
  
  # scores_auc <- riskRegression::Score(modelList,
  #                                     formula = formula, 
  #                                     data = data, 
  #                                     conf.int = TRUE, 
  #                                     times = horizon,
  #                                     # cens.method = "ipcw",
  #                                     cens.model = "cox",
  #                                     metrics = "AUC",
  #                                     weighting = "marginal"
  # )
  
  briers.i <- scores_brier %>% pluck("Brier") %>% pluck("score") %>% 
    dplyr::filter(model %in% name) %>% 
    dplyr::select(Brier,IPA) %>% as.numeric()# %>% t() %>% as.data.frame()
  
  # auc.i <- scores_auc %>% pluck("AUC") %>% pluck("score") %>% 
  #   dplyr::filter(model %in% name) %>% 
  #   dplyr::select(AUC) %>% as.numeric() # %>% t() %>% as.data.frame()
  auc.i <- scores_brier %>% pluck("AUC") %>% pluck("score") %>% 
    dplyr::filter(model %in% name) %>% 
    dplyr::select(AUC) %>% as.numeric() # %>% t() %>% as.data.frame()
  
  
  # res<- c(Dxy.i, Cind.i, R2.i, g.i,briers.i,auc.i)
  # res<- c(R2.i, g.i,briers.i,auc.i)
  res<- c(briers.i,auc.i)
  res
}


performanceCalculationMI <- function(name,formula,data,horizon){
  
  # name=modelName
  # formula=formula(fm %>% as.character())
  # data=data.list[[1]]
  # horizon=horizon
  ##########################################
  ###########
  ## START
  # Calculates measures for MI, fit model for each MI dataset
  environment(formula) <- environment()
  
  # Fit model
  model <-
    survival::coxph(eval(formula), data = data, x = T, y = T)
  
  
  # Brier and IPA)
  modelList = list(model)
  names(modelList)=name
  scores_brier <- riskRegression::Score(modelList,
                                        formula = formula,
                                        data = data,
                                        conf.int = TRUE,
                                        times = horizon,
                                        cens.method = "ipcw",
                                        cens.model = "cox",
                                        # metrics = "brier",
                                        summary = "ipa"
  )
  
  
  
  briers.i <- scores_brier %>% purrr::pluck("Brier") %>% purrr::pluck("score") %>%
    dplyr::filter(model %in% name) %>%
    dplyr::select(Brier,IPA) %>% as.numeric()# %>% t() %>% as.data.frame()
  
  
  
  auc.i <- scores_brier %>% purrr::pluck("AUC") %>% purrr::pluck("score") %>%
    dplyr::filter(model %in% name) %>%
    dplyr::select(AUC) %>% as.numeric() # %>% t() %>% as.data.frame()
  
  
  res<- c(briers.i,auc.i)
  res
}


performanceCalculationMI_horizons <- function(name,formula,data,horizon){
  
  # name=modelName
  # formula=formula(fm %>% as.character())
  # data=data.list[[1]]
  # horizon=horizon
  ##########################################
  ###########
  ## START
  # Calculates measures for MI, fit model for each MI dataset
  environment(formula) <- environment()
  
  # Fit model
  model <-
    survival::coxph(eval(formula), data = data, x = T, y = T)
  
  
  # Brier and IPA)
  modelList = list(model)
  names(modelList)=name
  scores_brier <- riskRegression::Score(modelList,
                                        formula = formula,
                                        data = data,
                                        conf.int = TRUE,
                                        times = horizon,
                                        cens.method = "ipcw",
                                        cens.model = "cox",
                                        # metrics = "brier",
                                        summary = "ipa"
  )
  
  
  briers.i <- scores_brier %>% purrr::pluck("Brier") %>% purrr::pluck("score") %>%
    dplyr::filter(model %in% name) %>%
    dplyr::select(times,Brier,IPA) #%>% as.numeric()# %>% t() %>% as.data.frame()
  
  
  auc.i <- scores_brier %>% purrr::pluck("AUC") %>% purrr::pluck("score") %>%
    dplyr::filter(model %in% name) %>%
    dplyr::select(times,AUC) #%>% as.numeric() # %>% t() %>% as.data.frame()
  
  
  res<- cbind(briers.i,auc.i[,2])
  res
}


performanceCalculationMI_riskPred <- function(name, fm, risks, data, horizon){
  # This function calculated performance measures across all imputations
  # Calculates measures for MI, fit model for each MI dataset
  # name=name
  # formula=formula
  # risks = risk.preds_test[[1]] # Calculated risks on test data, from applying pooled model of train data on test
  # data=data.mi.test[[1]] # original multiply imputed data, stacked- this is one imputed dataset not a list
  # horizon=horizon
  # 
  # risks %>% print
  ###########
  ## START
  environment(fm) <- environment()
  
  # Brier and IPA)
  riskList = list("risk"=risks)
  names(riskList)=name
  # print(names(riskList))
  # print(summary(data))
  scores_brier <- riskRegression::Score(riskList,
                                        formula = fm,
                                        data = data,
                                        conf.int = TRUE,
                                        times = horizon,
                                        cens.method = "ipcw",
                                        cens.model = "cox",
                                        # metrics = "brier",
                                        summary = "ipa"
  )
  
  # scores_auc <- riskRegression::Score(riskList,
  #                                     formula = fm,
  #                                     data = data,
  #                                     conf.int = TRUE,
  #                                     times = horizon,
  #                                     # cens.method = "ipcw",
  #                                     cens.model = "cox",
  #                                     metrics = "AUC",
  #                                     weighting = "marginal"
  # )
  
  briers.i <- scores_brier %>% purrr::pluck("Brier") %>% purrr::pluck("score") %>%
    dplyr::filter(model %in% name) %>%
    dplyr::select(Brier,IPA) %>% as.numeric()# %>% t() %>% as.data.frame()
  
  # auc.i <- scores_auc %>% pluck("AUC") %>% pluck("score") %>%
  #   dplyr::filter(model %in% name) %>%
  #   dplyr::select(AUC) %>% as.numeric() # %>% t() %>% as.data.frame()
  auc.i <- scores_brier %>% purrr::pluck("AUC") %>% purrr::pluck("score") %>%
    dplyr::filter(model %in% name) %>%
    dplyr::select(AUC) %>% as.numeric() # %>% t() %>% as.data.frame()
  
  
  # res<- c(Dxy.i, Cind.i, R2.i, g.i,briers.i,auc.i)
  # res<- c(R2.i, g.i,briers.i,auc.i)
  res<- c(briers.i,auc.i)
  res
}


performanceCalculationMI_riskPred_horizons <- function(name, fm, risks, data, horizon){
  # This function calculated performance measures across all imputations
  # Calculates measures for MI, fit model for each MI dataset
  # name=name
  # formula=formula
  # risks = risk.preds_test[[1]] # Calculated risks on test data, from applying pooled model of train data on test
  # data=data.mi.test[[1]] # original multiply imputed data, stacked- this is one imputed dataset not a list
  # horizon=horizon
  ###########
  ## START
  environment(fm) <- environment()
  
  # Brier and IPA)
  riskList = list("risk"=risks)
  names(riskList)=name
  
  scores_brier <- riskRegression::Score(riskList,
                                        formula = fm,
                                        data = data,
                                        conf.int = TRUE,
                                        times = horizon,
                                        cens.method = "ipcw",
                                        cens.model = "cox",
                                        # metrics = "brier",
                                        summary = "ipa"
  )
  
  
  
  briers.i <- scores_brier %>% purrr::pluck("Brier") %>% purrr::pluck("score") %>%
    dplyr::filter(model %in% name) %>%
    dplyr::select(times,Brier,IPA) #%>% as.numeric()# %>% t() %>% as.data.frame()
  
  
  auc.i <- scores_brier %>% purrr::pluck("AUC") %>% purrr::pluck("score") %>%
    dplyr::filter(model %in% name) %>%
    dplyr::select(times,AUC) #%>% as.numeric() # %>% t() %>% as.data.frame()
  
  
  res<- cbind(briers.i,auc.i[,2])
  res
}


getCoefs <- function(formula,data){
  model <- survival::coxph(formula, data = data, x = T, y = T)
  coefs <- coef(model) #%>% as.data.frame()
  coefs
}



getPredictedRisk <- function(baseCumHazard, time, lp){
  bh<-baseCumHazard[baseCumHazard$time==time,]$hazard
  predRisk <- 1-(exp(-bh)^exp(lp))
  predRisk
}

predictRisk.baseline.coefs <- function(model,horizons){
  # Input is a model fitted to some INITIAL DATA from which we want to extract the baseline and
  # coefficients table and then apply this model to OTHER DATA to generate
  # predicted risks on specific time horizons
  
  ## NEEDS fixing, it should take the baseline hazard of the train fit, and 
  ## the model of the test, in order to extract the test model matrix to be used
  # with the train baseline hazard
  ###########
  ## START
  # Extract coefficients  from fitted model on train data
  B.fit.train <- coef(model) # Equivalent to pooled coefficients from MI models
  # Extract baseline cumulative hazard function from fitted model on train data
  bh.train <- basehaz(model,centered=FALSE) # Equivalent to pooled baseline cum hazard from MI models
  # Interpolate base hazard to get rounded times, all times included in data
  times<- unique(c(floor(bh.train$time),ceiling(bh.train$time)))
  
  bh.train.inter<-data.frame(hazard=approx(bh.train$time, bh.train$hazard,xout=times,f=0, yleft=0,yright=length(bh.train$time),method="constant")$y,
                             time=approx(bh.train$time, bh.train$hazard,xout=times,f=0, yleft=0,yright=length(bh.train$time),method="constant")$x)
  
  # Extract  model matrix
  S.train<- model.matrix(model)
  # Calculate linear predictor
  lp.train <- S.train %*% B.fit.train
  # Calculate predicted risks
  pred.risk.horizons <-sapply(horizons, function(x) getPredictedRisk(baseCumHazard=bh.train.inter, time=x, lp=lp.train)
                              , simplify = FALSE)
  
  pred.risk.horizons <-do.call("cbind",pred.risk.horizons)
  pred.risk.horizons
}




predictRisk.baseline.coefs2 <- function(model,horizons){
  # Input is a model fitted to some INITIAL DATA from which we want to extract the baseline and
  # coefficients table and then apply this model to OTHER DATA to generate
  # predicted risks on specific time horizons
  
  ## NEEDS fixing, it should take the baseline hazard of the train fit, and 
  ## the model of the test, in order to extract the test model matrix to be used
  # with the train baseline hazard
  ###########
  ## START
  # Extract coefficients  from fitted model on train data
  B.fit.train <- coef(model) # Equivalent to pooled coefficients from MI models
  # Extract baseline cumulative hazard function from fitted model on train data
  bh.train <- basehaz(model,centered=FALSE) # Equivalent to pooled baseline cum hazard from MI models
  # Interpolate base hazard to get rounded times, all times included in data
  times<- horizons
  
  bh.train.inter<-data.frame(hazard=approx(bh.train$time, bh.train$hazard,xout=times,f=0, yleft=0,yright=length(bh.train$time),method="constant")$y,
                             time=approx(bh.train$time, bh.train$hazard,xout=times,f=0, yleft=0,yright=length(bh.train$time),method="constant")$x)
  
  # Extract  model matrix
  S.train<- model.matrix(model)
  # Calculate linear predictor
  lp.train <- S.train %*% B.fit.train
  # Calculate predicted risks
  pred.risk.horizons <-sapply(horizons, function(x) getPredictedRisk(baseCumHazard=bh.train.inter, time=x, lp=lp.train)
                              , simplify = FALSE)
  
  pred.risk.horizons <-do.call("cbind",pred.risk.horizons)
  pred.risk.horizons
}

predictRisk.baseline.coefs_onTest <- function(coefs, bh, modelMatrix,horizons){
  # INPUT
  ## coefs of model fitted on train
  ## bazeline cumulative function of model fitted on train
  ## model matrix of test data
  ## horizons
  # coefs = c.table
  # bh = bh.train
  # modelMatrix = mm.test
  # horizons = timesHorizon
  
  # OUT
  ## predictions made using test data
  
  ###########
  ## START
  # Calculate linear predictor
  lp.test <- modelMatrix %*% coefs
  # Calculate predicted risks
  pred.risk.horizons <-sapply(horizons, function(x) getPredictedRisk(bh, time=x, lp=lp.test)
                              , simplify = FALSE)
  
  pred.risk.horizons <-do.call("cbind",pred.risk.horizons)
  pred.risk.horizons
}



## Compute baseline hazard
baseline_strata <- function(model,horizons){
  # Input is a model fitted to some INITIAL DATA from which we want to extract the baseline and
  # coefficients table and then apply this model to OTHER DATA to generate
  # predicted risks on specific time horizons
  
  ## NEEDS fixing, it should take the baseline hazard of the train fit, and 
  ## the model of the test, in order to extract the test model matrix to be used
  # with the train baseline hazard
  ###########
  ## START
  # Extract coefficients  from fitted model on train data
  B.fit.train <- coef(model) # Equivalent to pooled coefficients from MI models
  # Extract baseline cumulative hazard function from fitted model on train data
  bh.train <- basehaz(model,centered=FALSE) # Equivalent to pooled baseline cum hazard from MI models
  # Interpolate base hazard to get rounded times, all times included in data
  times<- horizons
  
  bh.train.inter<-data.frame(hazard=approx(bh.train$time, bh.train$hazard,xout=times,f=0, yleft=0,yright=length(bh.train$time),method="constant")$y,
                             time=approx(bh.train$time, bh.train$hazard,xout=times,f=0, yleft=0,yright=length(bh.train$time),method="constant")$x)
  
  bh.train.inter
}



mypredictCox <- function(
                         coefs.table,
                         model.formula,
                         traindata,
                         lp.train,
                         covMeans=NULL,
                         times,
                         newdata=NULL,
                         centered = TRUE,
                         type=c("cumhazard","survival"),
                         keep.strata = TRUE,
                         keep.times = TRUE,
                         keep.newdata = FALSE,
                         keep.infoVar = FALSE,
                         se = FALSE,
                         band = FALSE,
                         iid = FALSE,
                         confint = (se+band)>0,
                         diag = FALSE,
                         average.iid = FALSE,
                         store.iid = "full"){
  
  #################################################
  # coefs.table=coef(fit.train.strata)
  # model.formula=modelFormula(time="time", event="event", "X1+X2+X7+X9+strata(X5)")
  # traindata=learndat
  # lp.train=fit.train.strata$linear.predictors
  # times=timesHorizon
  # newdata=testdat
  # centered = TRUE
  # type=c("cumhazard","survival")
  # keep.strata = TRUE
  # keep.times = TRUE
  # keep.newdata = FALSE
  # keep.infoVar = FALSE
  # se = FALSE
  # band = FALSE
  # iid = FALSE
  # confint = (se+band)>0
  # diag = FALSE
  # average.iid = FALSE
  # store.iid = "full"
  
  
  # coefs.table=coefs.table
  # model.formula=model.formula
  # traindata=traindata
  # lp.train=lp.train
  # covMeans=covMeans
  # newdata=newdata
  # times=times
  # iid = iid
  # diag = diag
  # average.iid = average.iid
  # store.iid = store.iid
  # type="survival"
  # 
  # se = TRUE
  # band = TRUE
  # iid = FALSE
  # confint = (se+band)>0
  # diag = FALSE
  # average.iid = FALSE
  # store.iid = "full"
  ###########
  ## START
  
  object =survival::coxph(formula=model.formula,data=traindata,x=TRUE,y=T)
  
  
  #################################################
  
  call <- match.call()
  ## centering
  if(!is.null(newdata)){
    if(inherits(centered,"data.frame")){
      df.reference <- centered
      centered2 <- TRUE ## for the linear predictor of the hazard
    }else{
      df.reference <- NULL
      centered2 <- centered ## for the linear predictor of the hazard
    }
    centered <- TRUE ## for the linear predictor of the baseline hazard
  }else{
    centered2 <- FALSE
  }
  
  ## ** Extract elements from object
  if (missing(times)) {
    nTimes <- 0
    times <- numeric(0)
  }else{
    nTimes <- length(times)
  }
  needOrder <- (nTimes[1]>0 && is.unsorted(times))
  if (all(!is.na(times)) && needOrder) {
    order.times <- order(times)
    oorder.times <- order(order.times)
    times.sorted <- sort(times)
  }else{
    if (nTimes==0){
      times.sorted <- numeric(0)
    } else{
      times.sorted <- times
      order.times <- 1:nTimes
      oorder.times <- 1:nTimes
    }
  }
  # Number of patients (total in data which the model was fitted in)
  object.n <- coxN(object)
  #Extract the design matrix used to train a Cox model
  object.modelFrame <- coxModelFrame(object)
  # Extract variable names from a model
  infoVar <- coxVariableName(object, model.frame = object.modelFrame)
  #Extract the type of estimator for the baseline hazard
  object.baseEstimator <- coxBaseEstimator(object)
  
  ## ease access
  is.strata <- infoVar$is.strata
  object.levelStrata <- levels(object.modelFrame$strata) ## levels of the strata variable
  nStrata <- length(object.levelStrata) ## number of strata
  nVar.lp <- length(infoVar$lpvars) ## number of variables in the linear predictor
  
  ## ** normalize model frame
  ## convert strata to numeric
  object.modelFrame[,c("strata.num") := as.numeric(.SD$strata) - 1]
  
  ## linear predictor
  ## if we predict the hazard for newdata then there is no need to center the covariates
  #Compute the linear predictor of a Cox model
  ## HERE CHANGES - LP POOLED WILL BE PROVIDED - must be ordered by TIME
  # object.modelFrame[,c("eXb") := exp(coxLP(object, data = NULL, center = if(is.null(newdata)){centered}else{FALSE}))]
  # object.modelFrame$eXb <- exp(lp.train)
  ## OR
  if(!is.null(covMeans)){
    object.modelFrame[,c("eXb") := exp(myCoxLP.coxph(object=object, 
                                                     data = NULL,
                                                     coef=coefs.table, 
                                                     lp=lp.train,
                                                     center = if(is.null(newdata)){centered}else{FALSE}))]
    
  }else{
    object.modelFrame[,c("eXb") := exp(myCoxLP.coxph_readyCovMeans(object=object, 
                                                                   data = NULL,
                                                                   coef=coefs.table, 
                                                                   lp=lp.train,
                                                                   covmeans=covMeans,
                                                                   center = if(is.null(newdata)){centered}else{FALSE}))]
    
  }
  ## add linear predictor and remove useless columns
  rm.name <- setdiff(names(object.modelFrame),c("start","stop","status","eXb","strata","strata.num"))
  if(length(rm.name)>0){
    object.modelFrame[,c(rm.name) := NULL]
  }
  
  ## sort the data
  object.modelFrame[, c("statusM1") := 1-.SD$status] ## sort by statusM1 such that deaths appear first and then censored events
  object.modelFrame[, c("XXXindexXXX") := 1:.N] ## keep track of the initial positions (useful when calling calcSeCox)
  data.table::setkeyv(object.modelFrame, c("strata.num","stop","start","statusM1"))
  
  ## last event time in each strata
  if(!is.null(attr(times,"etimes.max"))){ ## specified by the user
    etimes.max <- attr(times,"etimes.max")
    attr(times,"etimes.max") <- NULL
    attr(times.sorted,"etimes.max") <- etimes.max
  }else if(is.strata){ ## based on the data
    if(nVar.lp==0){
      iDTtempo <- object.modelFrame[, .SD[which.max(.SD$stop)], by = "strata.num"]
      etimes.max <- iDTtempo[,if(.SD$status==1){1e12}else{.SD$stop}, by = "strata.num"][[2]]
    }else{
      etimes.max <- object.modelFrame[, max(.SD$stop), by = "strata.num"][[2]]
    }
  }else{
    if(nVar.lp==0 && (utils::tail(object.modelFrame$status,1)==1)){ ## no covariates and ends by a death
      etimes.max <- 1e12
    }else{
      etimes.max <- max(object.modelFrame[["stop"]])
    }
  }
  
  ## ** checks
  ## check user imputs 
  if(nTimes[1]>0 && any(is.na(times))){
    stop("Missing (NA) values in argument \'times\' are not allowed.\n")
  }
  type <- tolower(type)
  if(any(type %in% c("lp","hazard","cumhazard","survival") == FALSE)){
    stop("type can only be \"lp\", \"hazard\", \"cumhazard\" or/and \"survival\" \n") 
  }
  if(is.null(newdata) && "lp" %in% type){
    stop("Cannot evaluate the linear predictor when argument \'newdata\' is missing. \n")
  }
  if(length(times)>1 && "lp" %in% type){
    stop("Cannot evaluate the linear predictor when there are multiple timepoints. \n")
  }
  ## predictCox is not compatible with all coxph/cph object (i.e. only handle only simple cox models)
  if(!is.null(object$weights) && !all(object$weights==1)){
    stop("predictCox does not know how to handle Cox models fitted with weights \n")
  }
  if(!is.null(object$naive.var)){
    stop("predictCox does not know how to handle frailty.") 
  }
  if(any(object.modelFrame[["start"]]!=0)){
    warning("The current version of predictCox was not designed to handle left censoring \n",
            "The function may be used on own risks \n") 
  }    
  if(object.baseEstimator == "exact"){
    stop("Prediction with exact handling of ties is not implemented.\n")
  }
  if(!is.null(object$call$tt)){
    stop("predictCox does not know how to handle time varying effects.\n") 
  }
  ## convergence issue
  if(!is.null(coef(object)) && any(is.na(coef(object)))){
    print(coef(object))
    stop("One or several parameters of the regression model have no value, i.e., a value 'NA'.\n")
  }
  ## prediction 
  if (missing(newdata) && (se || iid || average.iid)){
    stop("Argument 'newdata' is missing. Cannot compute standard errors in this case.")
  }
  if(!is.null(newdata)){
    if(nTimes[1]==0 && !identical(as.character(type),"lp")){
      stop("Time points at which to evaluate the predictions are missing \n")
    }
    if(!is.vector(times)){
      stop("Argument \'times\' must be a vector \n")
    }
    name.regressor <- c(infoVar$lpvars.original, infoVar$stratavars.original)
    if(length(name.regressor) > 0 && any(name.regressor %in% names(newdata) == FALSE)){
      stop("Missing variables in argument \'newdata\': \"",
           paste0(setdiff(name.regressor,names(newdata)), collapse = "\" \""),
           "\"\n")
    }
    if(se[1] && ("hazard" %in% type)){
      stop("Standard error cannot be computed for the hazard \n")
    }
    if(band[1] && ("hazard" %in% type)){
      stop("confidence bands cannot be computed for the hazard \n")
    }
  }
  ## diag argument
  if(!is.logical(diag)){
    stop("Argument \'diag\' must be of type logical \n")
  }
  if(diag){
    if(NROW(newdata)!=length(times)){
      stop("When argument \'diag\' is TRUE, the number of rows in \'newdata\' must equal the length of \'times\' \n")
    }
  }
  if(average.iid==TRUE && !is.null(attr(average.iid,"factor"))){
    if(is.null(store.iid) && !is.null(object$iid$store.iid)){
      store.iid <- object$iid$store.iid
    }
    if(store.iid == "full"){
      if(iid){
        stop("Attribute \"factor\" of argument \'average.iid\' not available when \'iid\' is TRUE with argument \'store.iid\' set to \"full\" \n",
             "Consider setting  \'store.iid\' set to \"minimal\" \n")
      }
      if(se){
        stop("Attribute \"factor\" of argument \'average.iid\' not available when \'se\' is TRUE with argument \'store.iid\' set to \"full\" \n",
             "Consider setting  \'store.iid\' set to \"minimal\" \n")
      }
    }
    test.list <- !is.list(attr(average.iid,"factor"))
    if(test.list){
      stop("Attribute \"factor\" of argument \'average.iid\' must be a list \n")
    }
    test.matrix <- any(unlist(lapply(attr(average.iid,"factor"), is.matrix))==FALSE)
    if(test.matrix){
      stop("Attribute \"factor\" of argument \'average.iid\' must be a list of matrices \n")
    }
    for(iFactor in 1:length(attr(average.iid,"factor"))){ ## iFactor <- 1
      ## when only one column and diag = FALSE, use the same weights at all times
      if((diag == FALSE) && (NCOL(attr(average.iid,"factor")[[iFactor]])==1) && (nTimes > 1)){
        attr(average.iid,"factor")[[iFactor]] <- matrix(attr(average.iid,"factor")[[iFactor]][,1],
                                                        nrow = NROW(attr(average.iid,"factor")[[iFactor]]),
                                                        ncol = nTimes, byrow = FALSE)
      }
      ## check dimensions
      if(any(dim(attr(average.iid,"factor")[[iFactor]])!=c(NROW(newdata), diag + (1-diag)*nTimes))){
        stop("Attribute \"factor\" of argument \'average.iid\' must be a list of matrices of size ",new.n,",",diag + (1-diag)*nTimes," \n")
      }
    }
  }
  
  ## ** baseline hazard
  
  ## compute the baseline hazard - I can probably use that, does not need fit object
  Lambda0 <- baseHaz_cpp(starttimes = object.modelFrame$start,
                         stoptimes = object.modelFrame$stop,
                         status = object.modelFrame$status,
                         eXb = object.modelFrame$eXb,
                         strata = object.modelFrame$strata.num,
                         nPatients = object.n,
                         nStrata = nStrata,
                         emaxtimes = etimes.max,
                         predtimes = times.sorted,
                         cause = 1,
                         Efron = (object.baseEstimator == "efron"))
  ## restaure strata levels
  Lambda0$strata <- factor(Lambda0$strata, levels = 0:(nStrata-1), labels = object.levelStrata)
  # print(Lambda0)
  # print(Lambda0)
  ## ** compute cumlative hazard and survival
  if (is.null(newdata)){  
    # print("NO NEW DATA")
    if (!("hazard" %in% type)){
      Lambda0$hazard <- NULL
    } 
    if ("survival" %in% type){  ## must be before cumhazard
      Lambda0$survival = exp(-Lambda0$cumhazard)
    }
    if (!("cumhazard" %in% type)){
      Lambda0$cumhazard <- NULL
    } 
    if (keep.times==FALSE){
      Lambda0$times <- NULL
    }
    if (keep.strata[[1]]==FALSE ||(is.null(call$keep.strata) && !is.strata)){
      Lambda0$strata <- NULL
    }
    if( keep.newdata[1]==TRUE){
      Lambda0$newdata <- object.modelFrame
    }
    # print(Lambda0)
    add.list <- list(lastEventTime = etimes.max,
                     se = FALSE,
                     band = FALSE,                         
                     type = type,
                     nTimes = nTimes,
                     baseline = TRUE,
                     var.lp = infoVar$lpvars.original,
                     var.strata = infoVar$stratavars.original)
    if(keep.infoVar){
      add.list$infoVar <- infoVar
    }
    # print(Lambda0)
    Lambda0[names(add.list)] <- add.list
    class(Lambda0) <- "predictCox"
    return(Lambda0)
  } else {
    # print("NEW DATA")
    out <- list()
    ## *** reformat newdata (compute linear predictor and strata)
    new.n <- NROW(newdata)
    ## newdata <- copy(newdata)
    setDT(newdata)
    ######### ISSUE HERE ##############
    
    ## Need to extract modelmatrix for test data, similar to the model matrix of train
    # to get the linear predictor of test data: model matrix test multiplied by coeficcients of train
    # Model matrix of test
    mm.new <- model.matrix(object,data = newdata)
    Xb <- mm.new%*%coefs.table %>% as.vector()
    # Xb <- coxLP(object, data = newdata, center = FALSE)
    if ("lp" %in% type){
      out$lp <- cbind(Xb)
      lp.iid <- centered2 ## when ask for centered then we need the iid for exporting the correct se
    }else{
      lp.iid <- FALSE
    }
    new.eXb <- exp(Xb)
    # Define the strata for a new dataset, Define the strata in a dataset to match those of a stratified Cox model
    # NOW ITS OKAY BECAUSE IT RUNS WITH NEW DATA
    new.strata <- coxStrata(object, data = newdata, 
                            sterms = infoVar$strata.sterms, 
                            strata.vars = infoVar$stratavars, 
                            strata.levels = infoVar$strata.levels)
    
    new.levelStrata <- levels(droplevels(new.strata))
    
    ## *** subject specific hazard
    if (is.strata==FALSE && !identical(as.character(type),"lp")){
      if(diag){
        if(needOrder){
          iTimes <- prodlim::sindex(jump.times = Lambda0$times, eval.times = times.sorted[oorder.times])
        }else{
          iTimes <- prodlim::sindex(jump.times = Lambda0$times, eval.times = times.sorted)
        }                  
      }
      
      if ("hazard" %in% type){
        if(diag){
          out$hazard <- cbind(new.eXb * Lambda0$hazard[iTimes])
        }else{
          out$hazard <- (new.eXb %o% Lambda0$hazard)
          if (needOrder) out$hazard <- out$hazard[,oorder.times,drop=0L]
        }
      }
      if ("cumhazard" %in% type || "survival" %in% type){
        if(diag){
          cumhazard <- cbind(new.eXb * Lambda0$cumhazard[iTimes])
        }else{
          cumhazard <- new.eXb %o% Lambda0$cumhazard
          if (needOrder){cumhazard <- cumhazard[,oorder.times,drop=0L]}
        }
        if ("cumhazard" %in% type){
          out$cumhazard <- cumhazard
        }
        if ("survival" %in% type){
          out$survival <- exp(-cumhazard)
        }
      }              
      
    }else if(!identical(as.character(type),"lp")){ 
      ## initialization
      if ("hazard" %in% type){
        out$hazard <- matrix(0, nrow = new.n, ncol = nTimes*(1-diag)+diag)
      }
      if ("cumhazard" %in% type){
        out$cumhazard <- matrix(NA, nrow = new.n, ncol = nTimes*(1-diag)+diag)                
      }
      if ("survival" %in% type){
        out$survival <- matrix(NA, nrow = new.n, ncol = nTimes*(1-diag)+diag)
      }
      
      ## loop across strata
      for(S in new.levelStrata){ ## S <- 1
        ####
        # S=new.levelStrata[1]
        ####
        id.S <- which(Lambda0$strata==S)
        newid.S <- which(new.strata==S)
        
        if(diag){
          if(needOrder){
            iSTimes <- prodlim::sindex(jump.times = Lambda0$times[id.S], eval.times = times.sorted[oorder.times[newid.S]])
          }else{
            iSTimes <- prodlim::sindex(jump.times = Lambda0$times[id.S], eval.times = times.sorted[newid.S])
          }                  
        }
        
        if ("hazard" %in% type){
          if(diag){
            out$hazard[newid.S] <- new.eXb[newid.S] * Lambda0$hazard[id.S][iSTimes]
          }else{
            out$hazard[newid.S,] <- new.eXb[newid.S] %o% Lambda0$hazard[id.S]
            if (needOrder){
              out$hazard[newid.S,] <- out$hazard[newid.S,oorder.times,drop=0L]
            }
          }
        }
        if ("cumhazard" %in% type || "survival" %in% type){
          if(diag){
            cumhazard.S <-  cbind(new.eXb[newid.S] * Lambda0$cumhazard[id.S][iSTimes])
          }else{
            cumhazard.S <-  new.eXb[newid.S] %o% Lambda0$cumhazard[id.S]
            if (needOrder){
              cumhazard.S <- cumhazard.S[,oorder.times,drop=0L]
            }
          }
          
          if ("cumhazard" %in% type){
            out$cumhazard[newid.S,] <- cumhazard.S
          }
          if ("survival" %in% type){
            out$survival[newid.S,] <- exp(-cumhazard.S)
          }
        }
      }
    }
    
    #################################
    # Take each new case and perform the standard operation
    # with the t matrix to get the pred. err.
    # predErr <- apply(mm.new, 1, function(X) sqrt(t(X) %*% fitObj$var %*% X))
    # results <- data.frame(out$survival, lwr=out$survival-predErr*1.96, upr=out$survival+predErr*1.96)
    #################################
    
    if(se[[1]] || band[[1]] || iid[[1]] || average.iid[[1]]){
      # print("CALCULATING SE")
      if({nVar.lp} > 0){
        ## get the (new) design matrix
        new.LPdata <- model.matrix(object, data = newdata)
        if(NROW(new.LPdata)!=NROW(newdata)){
          stop("NROW of the design matrix and newdata differ. \n",
               "Maybe because newdata contains NA values \n")
        }
        if(any(sort(colnames(new.LPdata))!=sort(names(coef(object))))){
          stop("Names of the design matrix and model parameters differ. \n",
               "Possible error in model.matrix due to special operator in the formula. \n")
        }
      }else{
        new.LPdata <- matrix(0, ncol = 1, nrow = new.n)
      }
      
      ## restaure original ordering
      data.table::setkeyv(object.modelFrame,"XXXindexXXX")
      if(diag){
        Lambda0$oorder.times <- oorder.times
      }else{
        Lambda0$oorder.times <- 1:nTimes
      }
      
      ## Computation of the influence function and/or the standard error
      export <- c("iid"[(iid+band+lp.iid)>0],"se"[(se+band)>0],"average.iid"[average.iid==TRUE])
      if(!is.null(attr(average.iid,"factor"))){
        if(diag){
          attr(export,"factor") <- attr(average.iid,"factor")
        }else{
          ## re-order columns according to times
          attr(export,"factor") <- lapply(attr(average.iid,"factor"), function(iF){
            iF[,order.times,drop=FALSE]
          })
        }
      }
      if(diag){
        times2 <- times
      }else{
        times2 <- times.sorted
      }
      attr(times2,"etimes.max") <- attr(times.sorted,"etimes.max")
      # Computation of standard errors for predictions
      outSE <- riskRegression:::calcSeCox(object,
                         times = times2,
                         nTimes = nTimes,
                         type = type,
                         diag = diag,
                         Lambda0 = Lambda0,
                         object.n = object.n,
                         object.time = object.modelFrame$stop,
                         object.eXb = object.modelFrame$eXb,
                         object.strata =  object.modelFrame$strata, 
                         nStrata = nStrata,
                         new.n = new.n,
                         new.eXb = new.eXb,
                         new.LPdata = new.LPdata,
                         new.strata = new.strata,
                         new.survival = if(diag){out$survival}else{out$survival[,order.times,drop=FALSE]},
                         nVar.lp = nVar.lp, 
                         export = export,
                         store.iid = store.iid)
      
      ## restaure orginal time ordering
      if((iid+band)>0){
        if ("lp" %in% type){
          out$lp.iid <- outSE$lp.iid
        }
        if ("hazard" %in% type){
          if (needOrder[1] && (diag[1] == FALSE))
            out$hazard.iid <- outSE$hazard.iid[,oorder.times,,drop=0L]
          else
            out$hazard.iid <- outSE$hazard.iid
        }
        if ("cumhazard" %in% type){
          if (needOrder[1] && (diag[1] == FALSE))
            out$cumhazard.iid <- outSE$cumhazard.iid[,oorder.times,,drop=0L]
          else
            out$cumhazard.iid <- outSE$cumhazard.iid
        }
        if ("survival" %in% type){
          if (needOrder[1] && (diag[1] == FALSE))
            out$survival.iid <- outSE$survival.iid[,oorder.times,,drop=0L]
          else
            out$survival.iid <- outSE$survival.iid
        }
      }
      if(average.iid == TRUE){
        if("lp" %in% type){
          out$lp.average.iid <- outSE$lp.average.iid
        }
        if ("hazard" %in% type){
          if (needOrder && (diag[1] == FALSE)){
            if(is.list(outSE$hazard.average.iid)){
              out$hazard.average.iid <- lapply(outSE$hazard.average.iid, function(iIID){iIID[,oorder.times,drop=0L]})
            }else{
              out$hazard.average.iid <- outSE$hazard.average.iid[,oorder.times,drop=0L]
            }
          }else{
            out$hazard.average.iid <- outSE$hazard.average.iid
          }
        }
        if ("cumhazard" %in% type){
          if (needOrder && (diag[1] == FALSE)){
            if(is.list(outSE$cumhazard.average.iid)){
              out$cumhazard.average.iid <- lapply(outSE$cumhazard.average.iid, function(iIID){iIID[,oorder.times,drop=0L]})
            }else{
              out$cumhazard.average.iid <- outSE$cumhazard.average.iid[,oorder.times,drop=0L]
            }
          }else{
            out$cumhazard.average.iid <- outSE$cumhazard.average.iid
          }
        }
        if ("survival" %in% type){
          if (needOrder && (diag[1] == FALSE)){
            if(is.list(outSE$survival.average.iid)){
              out$survival.average.iid <- lapply(outSE$survival.average.iid, function(iIID){iIID[,oorder.times,drop=0L]})
            }else{
              out$survival.average.iid <- outSE$survival.average.iid[,oorder.times,drop=0L]
            }
          }else {
            out$survival.average.iid <- outSE$survival.average.iid
          }
        }
        
        if(is.list(attr(average.iid,"factor"))){
          if("hazard" %in% type){
            names(out$hazard.average.iid) <- names(attr(average.iid,"factor"))
          }
          if("cumhazard" %in% type){
            names(out$cumhazard.average.iid) <- names(attr(average.iid,"factor"))
          }
          if("survival" %in% type){
            names(out$survival.average.iid) <- names(attr(average.iid,"factor"))
          }
        }
        
      }
      if((se+band)>0){
        if("lp" %in% type){
          out$lp.se <- outSE$lp.se
        }
        if ("cumhazard" %in% type){
          if (needOrder && (diag[1] == FALSE)){
            out$cumhazard.se <- outSE$cumhazard.se[,oorder.times,drop=0L]
          }else{
            out$cumhazard.se <- outSE$cumhazard.se
          }          
        }
        if ("survival" %in% type){
          if (needOrder && (diag[1] == FALSE)){
            out$survival.se <- outSE$survival.se[,oorder.times,drop=0L]
          } else{
            out$survival.se <- outSE$survival.se
          }          
        }
      }      
    }
    
    ## ** substract reference
    if("lp" %in% type && centered2){
      if(is.null(df.reference)){
        data <- try(eval(object$call$data), silent = TRUE)
        if(inherits(x=data,what="try-error")){
          stop("Could not evaluate the dataset used to fit the model to define a reference level. \n",
               "Set argument \'centered\' to FALSE or to a data.frame definining the reference level. \n")
        }
        var.original <- infoVar$lpvars.original
        
        ls.ref <- lapply(var.original, function(iVar){
          if(is.numeric(data[[iVar]])){
            return(unname(mean(data[[iVar]])))
          }else if(is.factor(data[[iVar]])){
            return(unname(factor(levels(data[[iVar]])[1],levels(data[[iVar]]))))
          }else if(is.character(data[[iVar]])){
            return(unname(sort(unique(data[[iVar]])[1])))
          }
        })
        df.reference <- as.data.frame(setNames(ls.ref,var.original))
      }
      
      ls.args <- as.list(call)[-1]
      ls.args$newdata <- df.reference
      ls.args$centered <- FALSE
      ls.args$type <- "lp"
      ls.args$se <- (se+band>0)
      ls.args$iid <- (iid+se+band>0)
      ls.args$band <- FALSE
      ls.args$confint <- FALSE
      outRef <- do.call(predictCox, args = ls.args)
      
      out$lp <- out$lp - as.double(outRef$lp)
      if(band[1] || se[1]){
        out$lp.se <- cbind(sqrt(colSums(colCenter_cpp(outSE$lp.iid, outRef$lp.iid)^2))) ## use outSe$lp.iid instead of out$lp.iid for when not exporting the iid
      }
      if(band[1] || iid[1]){
        out$lp.iid <- colCenter_cpp(out$lp.iid, outRef$lp.iid)
      }
      
    }
    
    ## ** add information to the predictions
    add.list <- list(lastEventTime = etimes.max,
                     se = se,
                     band = band,
                     type = type,
                     diag = diag,
                     nTimes = nTimes,
                     baseline = FALSE,
                     var.lp = infoVar$lpvars.original,
                     var.strata = infoVar$stratavars.original)
    if (keep.times==TRUE){
      add.list$times <- times
    }
    if (is.strata[1] && keep.strata[1]==TRUE){
      add.list$strata <- new.strata
    }
    
    if( keep.infoVar){
      add.list$infoVar <- infoVar
    }
    all.covars <- c(infoVar$stratavars.original,infoVar$lpvars.original)
    if( keep.newdata[1]==TRUE && length(all.covars)>0){
      if (data.table::is.data.table(newdata))
        add.list$newdata <- newdata[, all.covars, with = FALSE]
      else
        add.list$newdata <- newdata[, all.covars, drop=FALSE]
    }
    out[names(add.list)] <- add.list
    class(out) <- "predictCox"
    
    ## ** confidence intervals/bands
    if(confint){
      out <- stats::confint(out)
    }
    if(band[1] && se[1]==FALSE){
      out[paste0(type,".se")] <- NULL
    }
    if(band[1] && iid[1]==FALSE){
      out[paste0(type,".iid")] <- NULL
    }
    
    # return(list(preds=out,
    #             bazeCum=Lambda0))
    
    #################################
    # # Take each new case and perform the standard operation
    # # with the t matrix to get the pred. err.
    # predErr <- apply(mm.new, 1, function(X) sqrt(t(X) %*% fitObj$var %*% X))
    # results <- data.frame(out$survival, lwr=out$survival-predErr*1.96, upr=out$survival+predErr*1.96)
    #################################
    return(out)
  }
  
}



mypredictRisk.coxph <- function(coefs.table,
                                model.formula,
                                traindata,
                                lp.train,
                                covMeans=NULL,
                                newdata,
                                times,
                                product.limit = FALSE,
                                diag = FALSE,
                                iid = FALSE,
                                average.iid = FALSE,
                                centered = TRUE,
                                type=c("cumhazard","survival"),
                                keep.strata = TRUE,
                                keep.times = TRUE,
                                keep.newdata = FALSE,
                                keep.infoVar = FALSE,
                                se = FALSE,
                                band = FALSE,
                                confint = (se+band)>0,
                                store.iid = "full", fullRes=FALSE){
  
  ##################################
  # coefs.table=coefs.app_pooled
  # model.formula=form
  # traindata=traindata
  # lp.train=lp.app_pooled
  # covMeans=covmeans.app_pooled
  # newdata=x
  # times=horizon
  # product.limit = FALSE
  # diag = FALSE
  # iid = FALSE
  # average.iid = FALSE
  # centered = TRUE
  # type=c("cumhazard","survival")
  # keep.strata = TRUE
  # keep.times = TRUE
  # keep.newdata = FALSE
  # keep.infoVar = FALSE
  # se = FALSE
  # band = FALSE
  # confint = (se+band)>0
  # store.iid = "full"
  
  # coefs.table=c.table
  # model.formula=fm
  # traindata=data
  # lp.train= lp #+ sum(cox_model$means*coef(cox_model)),
  # covMeans=covmeans
  # newdata=testdata
  # times=timePoint
  # 
  # product.limit = FALSE
  # diag = FALSE
  # iid = FALSE
  # average.iid = FALSE
  # centered = TRUE
  # type=c("cumhazard","survival")
  # keep.strata = TRUE
  # keep.times = TRUE
  # keep.newdata = FALSE
  # keep.infoVar = FALSE
  # se = TRUE
  # band = TRUE
  # confint = TRUE#(se+band)>0
  # store.iid = "full"
  ##########
  ## START
  
  # dots <- list(...)
  # type <- dots$type ## hidden argument for ate
  # store.iid <- dots$store ## hidden argument for ate
  
  if(product.limit){
    outPred <- predictCoxPL(object=object,
                            newdata=newdata,
                            times=times,
                            iid = iid,
                            average.iid = average.iid,
                            keep.times=FALSE,
                            diag = diag,
                            store.iid = store.iid,
                            type="survival")
  }else{
    outPred <- mypredictCox(coefs.table=coefs.table,
                            model.formula=model.formula,
                            traindata=traindata,
                            lp.train=lp.train,
                            covMeans=covMeans,
                            newdata=newdata,
                            times=times,
                            iid = iid,
                            diag = diag,
                            average.iid = average.iid,
                            store.iid = store.iid,
                            type="survival",
                            se = se,
                            band = band,
                            confint = confint)
  }
  if(identical(type,"survival")){
    out <- outPred$survival
  }else{
    out <- 1-outPred$survival
  }
  if(iid){
    if(identical(type,"survival")){
      attr(out,"iid") <- outPred$survival.iid
    }else{
      attr(out,"iid") <- -outPred$survival.iid
    }
  }
  if(average.iid){
    if(identical(type,"survival")){
      attr(out,"average.iid") <- outPred$survival.average.iid
    }else{
      if(is.list(outPred$survival.average.iid)){
        attr(out,"average.iid") <- lapply(outPred$survival.average.iid, function(iIID){-iIID})
      }else{
        attr(out,"average.iid") <- -outPred$survival.average.iid
      }
    }
  }
  if(isTRUE(fullRes)){
    out<- outPred
  }
  return(out)
}


## Calculate linear predictor not based on fitted object but rather already the coefs table
myCoxLP.coxph <- function(object, data,coef,lp, center){
  
  # ####
  # object = object
  # data =NULL
  # center = if(is.null(newdata)){centered}else{FALSE}
  # coef=coefs.table
  # lp=lp.train
  #####
  ###########
  ## START
  n.varLP <- length(coef)
  
  if(is.null(data)){ ## training dataset
    # print("TRAIN")
    #Xb <- object$linear.predictors
    
    Xb <- lp
    
    if(center[[1]] == FALSE && n.varLP != 0){
      Xb <- Xb + sum(coxCenter(object)*coef)# ISSUE if we predict the hazard for newdata then there is no need to center the covariates
    }
    
  }else{ ## new dataset
    # print("NEWDATA")
    if(n.varLP>0){
      is.strata <- attr(object$terms, "special")$strata
      
      
      if(length(is.strata)>0){
        object.strata <- object[["strata"]]
        object[["strata"]] <- NULL # solve a bug in survival:::predict.coxph when fitting the model with x = TRUE
        
        Xb <- try(rowSums(stats::predict(object, newdata = as.data.frame(data), #### HERE IS THE ISSUE=make predictions based on lambda calculated above
                                         type = "terms")), silent = TRUE)
        if(inherits(x=Xb,what="try-error")){ ## Fix an error when the dataset used to fit the object is removed from the global environment
          ## survival:::predict.coxph search for it and read (at least) the status variable
          txt <- paste0("survival::predict.coxph returns the following error:\n",
                        as.character(Xb),
                        "It seems that the dataset used to fit the model is no more compatible with the model,\n",
                        "probably because it has been modified afterwards.\n",
                        "coxLP.coxph will try to reconstruct the original dataset and continue the execution.\n")
          warning(txt)
          ## So avoid an error, the following code re-create the original dataset
          object[["strata"]] <- object.strata
          object$call$data <- reconstructData(object)
          object[["strata"]] <- NULL
          
          Xb <- rowSums(stats::predict(object, newdata = as.data.frame(data),
                                       type = "terms"))
        }
      }else{
        Xb <- stats::predict(object, newdata = as.data.frame(data), type = "lp")
      }
      if(center == FALSE){
        Xb <- Xb + sum(coxCenter(object)*coef)
      }
      
    }else{
      Xb <- rep(0, NROW(data))
    }
  }
  
  return(unname(Xb))
}




myCoxLP.coxph_readyCovMeans <- function(object, data,coef,lp, covmeans, center){
  
  # ####
  # object = object
  # data =NULL
  # center = if(is.null(newdata)){centered}else{FALSE}
  # coef=coefs.table
  # lp=lp.train
  # object = object
  # data =NULL
  # center = if(is.null(newdata)){centered}else{FALSE}
  # coef=coefs.table
  # lp=lp.train
  # covmeans= myCoxLP.coxph_readyCovMeans
  # 
  
  ###########
  ## START
  n.varLP <- length(coef)
  
  if(is.null(data)){ ## training dataset
    # print("TRAIN")
    #Xb <- object$linear.predictors
    
    Xb <- lp
    
    if(center[[1]] == FALSE && n.varLP != 0){
      Xb <- Xb + sum(covmeans*coef)# ISSUE if we predict the hazard for newdata then there is no need to center the covariates
    }
    
  }else{ ## new dataset
    # print("NEWDATA")
    if(n.varLP>0){
      is.strata <- attr(object$terms, "special")$strata
      
      
      if(length(is.strata)>0){
        object.strata <- object[["strata"]]
        object[["strata"]] <- NULL # solve a bug in survival:::predict.coxph when fitting the model with x = TRUE
        
        Xb <- try(rowSums(stats::predict(object, newdata = as.data.frame(data), #### HERE IS THE ISSUE=make predictions based on lambda calculated above
                                         type = "terms")), silent = TRUE)
        if(inherits(x=Xb,what="try-error")){ ## Fix an error when the dataset used to fit the object is removed from the global environment
          ## survival:::predict.coxph search for it and read (at least) the status variable
          txt <- paste0("survival::predict.coxph returns the following error:\n",
                        as.character(Xb),
                        "It seems that the dataset used to fit the model is no more compatible with the model,\n",
                        "probably because it has been modified afterwards.\n",
                        "coxLP.coxph will try to reconstruct the original dataset and continue the execution.\n")
          warning(txt)
          ## So avoid an error, the following code re-create the original dataset
          object[["strata"]] <- object.strata
          object$call$data <- reconstructData(object)
          object[["strata"]] <- NULL
          
          Xb <- rowSums(stats::predict(object, newdata = as.data.frame(data),
                                       type = "terms"))
        }
      }else{
        Xb <- stats::predict(object, newdata = as.data.frame(data), type = "lp")
      }
      if(center == FALSE){
        Xb <- Xb + sum(covmeans*coef)
      }
      
    }else{
      Xb <- rep(0, NROW(data))
    }
  }
  
  return(unname(Xb))
}



lp.coxph <- function(object){
  return(object$linear.predictors)
}



#' Pooling performance measures over multiply imputed datasets
#'
#' \code{pool_performance_internal} Pooling performance measures
#'
#' @param data Data frame with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded 
#'   from the dataset.
#' @param formula A formula object to specify the model as normally used by glm.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes 
#'   the imputed datasets.
#'
#' @keywords internal 
#'  
#' @export
pool_performance_internal.cox <- function(data, 
                                          formula, 
                                          nimp, 
                                          impvar,
                                          horizon,
                                          time,
                                          status,
                                          modelName)
{
  # # # # #####################
  # data=boot_data_compl
  # nimp = nimp_mice
  # impvar=impVar
  # formula = fm
  # horizon=horizon
  # time=Y_time
  # status=Y_status
  # modelName=modelName
  # # # # # #####################
  ###########
  ## START
  
  coef_f <- se_f <- pred.group <- obs.group <- list()
  
  if(is_empty(formula))
    stop("\n", "Model not specified in formula object")
  environment(formula) <- environment()
  
  fm <-
    formula
  
  # perf_stats <-
  #   matrix(NA, nimp, 3)
  
  # Extract performance measures for EACH MI dataset
  ## Foreach package
  # tic()
  data.list <-data %>% as.data.frame() %>% miceafter::df2milist(impvar = impvar) # ,.packages=c("survival","rms","dplyr","purrr"),.export=c("data.list","fm")
  
  # # #Following fails when strata are included
  perf_stats<-foreach(i = seq_along(data.list),.combine = 'rbind',.export = c("performanceCalculationMI_horizons","strata",
                                                                              "fm","rcs")
  ) %do% {#,.combine = 'rbind'
    
    g <-performanceCalculationMI_horizons(name=modelName,
                                          formula=formula(fm %>% as.character()),
                                          data=data.list[[i]],
                                          horizon=horizon)
    # g
  }
  
  # Create list of performance measures, each list object is for specific time
  perf_stats<-perf_stats %>%
    group_split(times)
  names(perf_stats) <- horizon
  
  # Get coefficients
  coef_f <-sapply(data.list, function(x) getCoefs(fm,x), simplify = FALSE)
  
  # C.index
  C.index <-data %>% miceafter::df2milist(impvar = impvar) %>%
    with(miceafter::cindex(coxph(formula=as.formula(as.character(fm))))) %>%
    miceafter::pool_cindex()
  # Dxy
  Dxy <-(C.index - 0.5)*2
  # End pooling performance measures in multiply imputed data
  
  # Pooling over measure values across MI datasets
  # Time dependent ROC (just average)
  # Go through each time and pool values
  # roc_res <-
  #   pool_auc(perf_stats[, 3], sqrt(pROC::var(perf_stats[, 3])),
  #            nimp = nimp, log_auc = FALSE)
  
  roc_res <-sapply(perf_stats, function(x) pool_auc(x[, "AUC"], sqrt(pROC::var(x[, "AUC"])),
                                                    nimp = nimp, log_auc = FALSE), simplify=FALSE)
  
  
  
  # Pool Brier
  # brier_pool <-
  #   mean(perf_stats[, 1])
  
  brier_pool <-sapply(perf_stats, function(x) mean(x[, "Brier"] %>% pull(),na.rm=T), simplify=FALSE)
  
  
  # Pool IPA
  # ipa_pool <-
  #   mean(perf_stats[, 2])
  
  ipa_pool <-sapply(perf_stats, function(x) mean(x[, "IPA"] %>% pull(),na.rm=T), simplify=FALSE)
  #
  
  # Colmeans of predictors in multiply imputed datasets
  coef_pooled <-
    colMeans(do.call("rbind", coef_f))
  
  
  # Pooling Performance measures and coefficients
  pobjperform <- list(coef_pooled=coef_pooled,
                      # R2_pooled=rsq.n,
                      Dxy_pooled=Dxy,
                      C.index_pooled= C.index,
                      timeROC_pooled=roc_res,
                      Brier_pooled = brier_pool,
                      IPA_pooled=ipa_pool,
                      nimp=nimp)
  # Pooled info in each bootstrap sample
  pobjperform
}


my.pool_RR_covMeans <- function(impList,estimate,impN,n,type=c("lp","means"),impvar){
  # impList = pool_covCenter.list
  # estimate = "LDH_level"
  # impN = 5
  # n = nrow(data %>% dplyr::filter(.imp==1))
  # type="means"
  # impvar=impvar
  # 
  
  ###
  Q <- array(dim = impN)
  U <- array(dim = impN)
  
  for (i in 1:impN){
    print(paste0("Imputation:",i,"| Parameter:", estimate))
    print(impList[[i]][[estimate]])
    # i=1
    Q[i] <- impList[[i]][[estimate]]
    if(type=="means"){
      U[i] <- var(impList[[i]][[estimate]])/n
      
    }else{
      U[i] <- var(impList[[i]][[estimate]])
      
    }
  }
  B <-  var(Q)
  mean_of_means <- mean(Q) # OVerall estimate Q over imputations
  total_variance_of_means <- mean(U) + (1 + 1/impN) * B #  Total variance T for overall MI estimate
  # CI calculation
}

stderror <- function(x) sd(x)/sqrt(length(x))


my.pool_RR_covMeans <- function (est, se, conf.level = 0.95, n, k){
  m <- length(est) # number of imputations
  mean_est <- mean(est) # overall estimate Q
  var_w <- mean(se^2)/n # Within imputation variance U, for means it must be divided by number of observations
  var_b <- var(est) # Between imputation variance B
  var_T <- var_w + (1 + (1/m)) * var_b # Total variance for overall MI estimate
  se_total <- sqrt(var_T)
  r <- (1 + 1/m) * (var_b/var_w)
  v_old <- (m - 1) * (1 + (1/r))^2
  lambda <- (var_b + (var_b/m))/var_T
  v_obs <- (((n - k) + 1)/((n - k) + 3)) * (n - k) * (1 - lambda)
  v_adj <- (v_old * v_obs)/(v_old + v_obs)
  alpha <- 1 - (1 - conf.level)/2
  t_stats <- mean_est/se_total
  p_val <- 2 * pt(-abs(t_stats), df = v_adj)
  t <- qt(alpha, v_adj)
  ci_upper <- mean_est + (t * se_total)
  ci_lower <- mean_est - (t * se_total)
  res <- round(c(mean_est, ci_lower, ci_upper, p_val), 7)
  names(res) <- c("Estimate", "95% CI L", "95% CI U", 
                  "P-val")
  # return(res)
  mean_est
}


cindex_fromRiskPrediction <- function(formula,data){
  
  # data = data.mi.test[[1]]
  # risks = risk.preds_test[[1]][,1]
  # time=Y_time
  # status=Y_status
  
  # cox.formula=modelFormula(time, status, risks)
  # cpred<-survConcordance(as.formula(paste0("Surv(",time,",",status,")~","risks")), data)
  
  # data=data.risks[[1]]
  # mf_call <- match.call()
  # fit <- eval(mf_call[[2L]], parent.frame())
  risks <- data$risks
  time <- data$time
  status <-data$status
  data <- data$data
  
  cfit <-concordance(coxph(as.formula(paste0("Surv(",time,",",status,")~","risks")), data))
  
  cest <- cfit$concordance
  cse <- sqrt(cfit$var)
  n <- cfit$n
  dfcom <- n - 1
  output <- matrix(c(cest, cse, dfcom), 1, 3)
  colnames(output) <- c("c-index", "se", "dfcom")
  return(output)
  # return(mf_call)
}



pool_performance_nestedBootMI.cox <- function(data, 
                                              form, 
                                              nimp, 
                                              impvar,
                                              horizon,
                                              name,
                                              coefs.app_pooled,
                                              traindata,
                                              lp.app_pooled,
                                              covmeans.app_pooled,
                                              time,
                                              status){
  
  # data=pobj$data
  # form=fm
  # nimp=nimp_mice
  # impvar=impVar
  # horizon=horizon
  # name=modelName
  # coefs.app_pooled = lp_app_pooled
  # traindata=boot_train
  # lp.app_pooled=pool_model_cox$pool_lp_final  %>% magrittr::extract2(1)
  # covmeans.app_pooled=pool_model_cox$pool_covCenter_final %>% magrittr::extract2(1)
  # time=Y_time
  # status = Y_status
  ###########
  
  # fit.pooled2<-fit.pooled
  # fit.pooled$linear.predictors<- fit.pooled$linear.predictors- sum(coxCenter(fit.pooled)*fit.pooled$coefficients)
  # 
  ## START
  
  if(is_empty(form))
    stop("\n", "Model not specified in formula object")
  
  environment(form) <- environment()
  # List of imputed datasets
  data.mi.test<- data %>% df2milist(impvar = impvar)
  # Extract risk predictions for all imputed datasets
  risk.preds_test <- sapply(data.mi.test, 
                            function(x) mypredictRisk.coxph(coefs.table=coefs.app_pooled,
                                                            model.formula=form,
                                                            traindata=traindata,
                                                            lp.train=lp.app_pooled,
                                                            covMeans=covmeans.app_pooled,
                                                            newdata=x,
                                                            times=horizon,
                                                            product.limit = FALSE,
                                                            diag = FALSE,
                                                            iid = FALSE,
                                                            average.iid = FALSE,
                                                            centered = TRUE,
                                                            type=c("cumhazard","survival"),
                                                            keep.strata = TRUE,
                                                            keep.times = TRUE,
                                                            keep.newdata = FALSE,
                                                            keep.infoVar = FALSE,
                                                            se = FALSE,
                                                            band = FALSE,
                                                            confint = FALSE,
                                                            store.iid = "full"),
                            simplify = FALSE)
  
  
 
  
  #########################
  # Calculate performance measures over imputed datasets using predicted risks
  # Extract performance measures for EACH MI dataset
  # loop over number of imputed datasets
  
  ## Foreach package
  # tic()
  perf_stats<-foreach(i = 1:nimp,.combine = 'rbind',.export = "performanceCalculationMI_riskPred_horizons") %do% {
    g <-performanceCalculationMI_riskPred_horizons(name=name,
                                                   fm=form,
                                                   risks = risk.preds_test[[i]], # Calculated risks on test data, from applying pooled model of train data on test
                                                   data=data.mi.test[[i]], # original multiply imputed data, stacked- this is one imputed dataset not a list
                                                   horizon=horizon)
    # g
  }
  # toc()
  
  
  # Create list of performance measures, each list object is for specific time
  perf_stats<-perf_stats %>%
    group_split(times)
  names(perf_stats) <- horizon
  
  
  # C.index-is calculated for one time point, usually its the maximum time in the data-As I tested, time does NOT MATTER, so I can give as input one list
  # of predicted risks at ANY TIMEPOINT
  # Make a list of data and risks
  data.risks <- sapply(1:nimp, function(x) list(data=data.mi.test[[x]],
                                                risks=risk.preds_test[[x]][,1],
                                                time=time,
                                                status=status
  ), simplify = FALSE)
  
  
  ## Calculate Cindex over imputed datasets - Cannot use with so do it manually
  statistics <- as.list(seq_len(length(data.risks)))
  # output statistical analysis
  for (i in seq_along(statistics)) {
    # i=1
    df_m <- data.risks[[i]]
    statistics[[i]] <- eval(expr = cindex_fromRiskPrediction(df_m,formula =as.formula(paste0("Surv(",time,",",status,")~","risks")) ),
                            envir = df_m, enclos = parent.frame())
    if (is.expression(statistics[[i]])) {
      statistics[[i]] <- eval(expr = statistics[[i]],
                              envir = df_m, enclos = parent.frame())
    }
  }
  
  # return the repeated statistical analyses
  obj <- list(call = call, statistics = statistics)
  class(obj) <- c("mistats")
  
  C.index <-pool_cindex(obj)
  
  # Dxy 
  Dxy <-(C.index - 0.5)*2
  # End pooling performance measures in multiply imputed data
  
  # Pooling over measure values across MI datasets
  # Time dependent ROC (just average)
  # roc_res <-
  #   pool_auc(perf_stats[, 3], sqrt(pROC::var(perf_stats[, 3])),
  #            nimp = nimp, log_auc = FALSE)
  roc_res <-sapply(perf_stats, function(x) pool_auc(x[, "AUC"], sqrt(pROC::var(x[, "AUC"])),
                                                    nimp = nimp, log_auc = FALSE), simplify=FALSE)
  
  # Pool Brier
  # brier_pool <-
  #   mean(perf_stats[, 1])
  brier_pool <-sapply(perf_stats, function(x) mean(x[, "Brier"] %>% pull(),na.rm=T), simplify=FALSE)
  # Pool IPA
  # ipa_pool <-
  #   mean(perf_stats[, 2])
  ipa_pool <-sapply(perf_stats, function(x) mean(x[, "IPA"] %>% pull(),na.rm=T), simplify=FALSE)
  
  # # Colmeans of predictors in multiply imputed datasets
  # coef_pooled <-
  #   colMeans(do.call("rbind", coef_f))
  
  
  # Pooling Performance measures and coefficients
  pobjperform <- list(#coef_pooled=coef_pooled,
    # R2_pooled=rsq.n,
    Dxy_pooled=Dxy,
    C.index_pooled= C.index,
    timeROC_pooled=roc_res,
    Brier_pooled = brier_pool, 
    IPA_pooled=ipa_pool,
    nimp=nimp)
  # Pooled info in each bootstrap sample
  pobjperform
  
  
}



# ## Moderate calibration - fixed time point
# ## Bootstrap confidence intervals
# ## for the cal------


#
# Bootstrap calibration measures
numsum_boot <- function(split) {
  # split <- boot_straps$splits[[1]]
  # split <- as.data.frame(split)
  pred <- riskRegression::predictRisk(fitObj,
                                      newdata = analysis(split),
                                      times = horizon)



  pred.cll <- log(-log(1 - pred))

  # Estimate actual risk - basic model
  # d <- datadist(as.data.frame(split))
  # options(datadist=d)
  vcal <- rms::cph(formula=modelFormula(time,status,"rcs(pred.cll, 4)"),
                   x = T,
                   y = T,
                   surv = T,
                   data = analysis(split)
  )


  # Save objects needed
  cox.surv <- rms::survest(vcal,
                           times = horizon,
                           newdata = analysis(split))
  db_cal_boot <- data.frame(
    "obs" = 1 - cox.surv$surv,

    "lower" = 1 - cox.surv$upper,

    "upper" = 1 - cox.surv$lower,

    "pred" = pred
  )

  # data.frame(hazard=approx(db_cal_boot$obs, db_cal_boot$pred,xout=times,f=0, yleft=0,yright=length(bh.train$time),method="constant")$y,
  #
  #            time=approx(bh.train$time, bh.train$hazard,xout=times,f=0, yleft=0,yright=length(bh.train$time),method="constant")$x)

  absdiff_boot <- abs(db_cal_boot$obs - db_cal_boot$pred)


  res_cal_boot <- data.frame(
    "ICI" = mean(absdiff_boot),
    "E50" = quantile(absdiff_boot, probs = .5),
    "E90" = quantile(absdiff_boot, probs = .9)
  )
  res_cal_boot

  # db_cal_boot
  # db_cal_boot <- db_cal_boot[order(db_cal_boot$pred), ]
  # 
  # par(xaxs = "i", yaxs = "i", las = 1)
  # plot(
  #   db_cal_boot$pred,
  #   db_cal_boot$obs,
  #   type = "l",
  #   lty = 1,
  #   xlim = c(0, 1),
  #   ylim = c(0, 1),
  #   lwd = 2,
  #   xlab = "Predicted risk from developed model",
  #   ylab = "Predicted risk from refitted model", bty = "n"
  # )
  # lines(db_cal_boot$pred,
  #       db_cal_boot$lower,
  #       type = "l",
  #       lty = 2,
  #       lwd = 2)
  # lines(db_cal_boot$pred,
  #       db_cal_boot$upper,
  #       type = "l",
  #       lty = 2,
  #       lwd = 2)
  # abline(0, 1, lwd = 2, lty = 2, col = 2)
  # legend("bottomright",
  #        c("Ideal calibration",
  #          "Calibration curve based on secondary Cox model",
  #          "95% confidence interval"),
  #        col = c(2, 1, 1),
  #        lty = c(2, 1, 2),
  #        lwd = c(2, 2, 2),
  #        bty = "n",
  #        cex = 0.85)
  # title("Basic model - validation data ")


}




getSurvProb <- function(data, testdata, fitObj, time, status,manualCalc, timePoint, features, predName){
  
  time.un <- as.name(time)
  testdata <- testdata %>% mutate(time.un = timePoint)
  
  c.table=coef(fitObj)
  lp = fitObj$linear.predictors
  covmeans = fitObj$means
  fm =modelFormula(time=time, event=status, features)
  
  
  if(isTRUE(manualCalc)){
    predRisk <- predictRisk(fitObj, newdata =testdata ,times = timePoint) %>% as.vector() %>% round(4)
  }else{
    predRisk <- mypredictRisk.coxph(coefs.table=c.table,
                                    model.formula=fm,
                                    traindata=data,
                                    lp.train= lp ,#+ sum(cox_model$means*coef(cox_model)),
                                    covMeans=covmeans,
                                    newdata=testdata,
                                    times=timePoint) %>% as.vector() %>% round(4)
  }
  
  survProb <- 1-predRisk
  
  testdata$survival <- round(survProb,2)*100
  res = list(data=testdata,
             prob=survProb)
  res
  
}


getSurvProb_nice <- function(data, testdata, fitObj, time, status,manualCalc, timePoint, features, predName, strataVar="First_line_treatments"){
  # print("HERE")
  ####################
  # data=somData
  # testdata=somData[605,]%>% dplyr::select(-First_line_treatments)
  # fitObj=fit.pooled2
  # time=surv_time
  # status=surv_status
  # manualCalc=F
  # timePoint=myHorizons[c(1,2,4,6)]
  # features=feats.trans.strata.f
  # predName="MI_model"
  # strataVar= "First_line_treatments"
  ####################
  if(!strataVar %in% colnames(testdata)){
    show=length(data[[strataVar]] %>% levels())
    if(show==2){
      testdata <- rbind(testdata,testdata)
      
    }else if(show==3){
      testdata <- rbind(testdata,testdata,testdata)
      
    }else if(show==4){
      testdata <- rbind(testdata,testdata,testdata,testdata)
      
    }
    testdata[[strataVar]]<- factor(levels(data[[strataVar]]), levels = levels(data[[strataVar]]))
    
  }else{
    show=1
  }
  
  if(nrow(testdata)<=6){
    testdata<- rbind(testdata,data)
  }
  
  
  # 
  # time.un <- as.name(time)
  # testdata <- testdata %>% mutate(time.un = timePoint)
  # 
  c.table=coef(fitObj)
  lp = fitObj$linear.predictors
  covmeans = fitObj$means
  fm =modelFormula(time=time, event=status, features)

  # print(testdata)
  if(isTRUE(manualCalc)){
    predRisk <- mypredictRisk.coxph(coefs.table=c.table,
                                    model.formula=fm,
                                    traindata=data,
                                    lp.train= lp ,#+ sum(cox_model$means*coef(cox_model)),
                                    covMeans=covmeans,
                                    newdata=testdata,
                                    times=timePoint, type="survival", se=T, confint=T, band=T,fullRes=T) #%>% as.vector() %>% round(4)
    # predRisk <- cbind(
    #                   predRisk$strata %>% as.character,
    #                   predRisk$times,
    #                   predRisk$survival %>% as.vector(),
    #                   predRisk$survival.se %>% as.vector(),
    #                   predRisk$survival.lower %>% as.vector(),
    #                   predRisk$survival.upper %>% as.vector()) %>% as.data.frame()
    # names(predRisk) <- c("strata","times","survival", "se", "lower","upper")
    predRisk <- print(predRisk) %>% as.data.frame()
    
  }else{
    print("CORRECT")
    # predRisk2 <- riskRegression::predictRisk(fitObj, newdata =testdata ,times = timePoint, type="survival",se=T, band=T,confint=T)# %>% as.vector() %>% round(4)
    predRisk <- riskRegression::predictCox(fitObj, newdata =testdata ,times = timePoint, se=T, confint=T,type="survival")
    predRisk <- print(predRisk) %>% as.data.frame()

  }
  
  #################################
  # Take each new case and perform the standard operation
  # with the t matrix to get the pred. err.
  # mm.new <- model.matrix(fitObj,data = testdata)
  # predErr <- apply(mm.new, 1, function(X) sqrt(t(X) %*% fitObj$var %*% X)) 
  # predErr.surv <- 1-exp(-predErr)
  # results <- data.frame(as.vector(predRisk), lwr=as.vector(predRisk)-predErr.surv*1.96, upr=as.vector(predRisk)+predErr.surv*1.96)
  #################################
  
  ######################
  
  testdata<- testdata[1:show,]
  
  predRisk <- predRisk %>% dplyr::filter(observation %in% 1:show)
  # predRisk <- predRisk[1:show,]
  # survProb <- 1-predRisk %>% as.data.frame()
  # survProb.list <- as.list(survProb)
  # survProb.list <- as.list(predRisk)
  # names(survProb.list) <-timePoint  %>% as.character()
  # testdata.final <-sapply(1:length(survProb.list) , function(x) testdata %>%
  #                           add_column(time=rep(names(survProb.list)[x],nrow(testdata))) %>%
  #                           add_column(survival=round(survProb.list[[x]],2)*100), simplify=FALSE) %>% ldply()
  # 
  testdata.final <- testdata
  if(length(timePoint)>1){
    for(i in 2:length(timePoint)){
      testdata.final <- rbind(testdata.final,testdata)
    }
  }
  
  testdata.final<- cbind(testdata.final,predRisk %>% dplyr::select(strata, times, survival, survival.lower,survival.upper,survival.se))
  testdata.final$survival <- round(testdata.final$survival,2)*100
  testdata.final$survival.lower<- round(testdata.final$survival.lower,2)*100
  testdata.final$survival.upper <- round(testdata.final$survival.upper,2)*100
  

  # print(head(testdata.final))
  res = list(data=testdata.final,
             times=timePoint)
  res
  
}



getSurvProb_nice2 <- function(data, testdata, fitObj, time, status,manualCalc, timePoint, features, predName, strataVar="First_line_treatments"){
  # print("HERE")
  ####################
  # data=somData
  # testdata=somData[605,]%>% dplyr::select(-First_line_treatments)
  # fitObj=fit.pooled2
  # time=surv_time
  # status=surv_status
  # manualCalc=F
  # timePoint=myHorizons[c(1,2,4,6)]
  # features=feats.trans.strata.f
  # predName="MI_model"
  # strataVar= "First_line_treatments"
  ####################
  if(!strataVar %in% colnames(testdata)){
    show=length(data[[strataVar]] %>% levels())
    if(show==2){
      testdata <- rbind(testdata,testdata)
      
    }else if(show==3){
      testdata <- rbind(testdata,testdata,testdata)
      
    }else if(show==4){
      testdata <- rbind(testdata,testdata,testdata,testdata)
      
    }
    testdata[[strataVar]]<- factor(levels(data[[strataVar]]), levels = levels(data[[strataVar]]))
    
  }else{
    show=1
  }
  
  # 
  c.table=coef(fitObj)
  lp = fitObj$linear.predictors
  covmeans = fitObj$means
  fm =modelFormula(time=time, event=status, features)
  
  # print(testdata)
  if(isTRUE(manualCalc)){
    predRisk <- mypredictRisk.coxph(coefs.table=c.table,
                                    model.formula=fm,
                                    traindata=data,
                                    lp.train= lp ,#+ sum(cox_model$means*coef(cox_model)),
                                    covMeans=covmeans,
                                    newdata=testdata,
                                    times=timePoint, type="survival", se=T, confint=T, band=T,fullRes=T) #%>% as.vector() %>% round(4)
    
    predRisk <- print(predRisk) %>% as.data.frame()
    
  }else{
    print("CORRECT")
    # predRisk2 <- riskRegression::predictRisk(fitObj, newdata =testdata ,times = timePoint, type="survival",se=T, band=T,confint=T)# %>% as.vector() %>% round(4)
    predRisk <- riskRegression::predictCox(fitObj, newdata =testdata ,times = timePoint, se=T, confint=T,type="survival")
    predRisk <- print(predRisk) %>% as.data.frame()
    
  }
  
  #################################
  # Take each new case and perform the standard operation
  # with the t matrix to get the pred. err.
  # mm.new <- model.matrix(fitObj,data = testdata)
  # predErr <- apply(mm.new, 1, function(X) sqrt(t(X) %*% fitObj$var %*% X)) 
  # predErr.surv <- 1-exp(-predErr)
  # results <- data.frame(as.vector(predRisk), lwr=as.vector(predRisk)-predErr.surv*1.96, upr=as.vector(predRisk)+predErr.surv*1.96)
  #################################
  
  ######################
  
  testdata<- testdata[1:show,]
  
  predRisk <- predRisk %>% dplyr::filter(observation %in% 1:show)
  
  testdata.final <- testdata
  if(length(timePoint)>1){
    for(i in 2:length(timePoint)){
      testdata.final <- rbind(testdata.final,testdata)
    }
  }
  
  testdata.final<- cbind(testdata.final,predRisk %>% dplyr::select(strata, times, survival, survival.lower,survival.upper,survival.se))
  testdata.final$survival <- round(testdata.final$survival,2)*100
  testdata.final$survival.lower<- round(testdata.final$survival.lower,2)*100
  testdata.final$survival.upper <- round(testdata.final$survival.upper,2)*100
  
  
  # print(head(testdata.final))
  res = list(data=testdata.final,
             times=timePoint)
  res
  
}



getHazard <- function(data,fitObj, time, status, features, timePoint){
  # data=finalModel_noInter$si_data
  # fitObj = finalModel_noInter$bw
  # time=surv_time
  # status=surv_status
  # features <-finalModel_noInter$feats_bw.c
  # timePoint <-c(12,24,36)
  ################################
  c.table=coef(fitObj)
  lp = fitObj$linear.predictors
  covmeans = fitObj$means
  fm =modelFormula(time=time, event=status, features)
  
  baseHaz <-mypredictCox(coefs.table=c.table,
                         model.formula=fm,
                         traindata=data,
                         lp.train=lp,
                         covMeans=covmeans,
                         newdata=NULL,
                         times=timePoint,type=c("hazard","cumhazard","survival"))
  baseHaz
}

makeClinExample <- function(data, test, model, time, status, calcPred, horizon, features, Modname, strataVar, selectStrat,colsRemove,finalColsName,nRows){
  
  # data=somData
  # test=patient_ex[1,] %>% dplyr::select(-First_line_treatments)
  # model=finalModel_noInter$bw
  # time=surv_time
  # status=surv_status
  # calcPred=TRUE
  # horizon=myHorizons[c(2)]
  # features=finalModel_noInter$feats_bw.c
  # Modname="MI_model"
  # strataVar="First_line_treatments"
  # selectStrat="BRAFi"
  # colsRemove=c("Patient_ID", "Disease_stage","Performance_status.v1","Gender","Comorbidities", "Autoimmunity" , "Other_malignancies","BRAF_status",
  #               "LDH_level.log","strata","times","First_line_treatments",
  #               "Steroids","Performance_status","Stage_Brain_mets","Disease_stage")
  # finalColsName=c("Age", "LDH level", "Estimated 1-Year Survival, %")
  # nRows=1
  ###
  strata.un <- as.name(strataVar)
  
  example.Surv.bw<- sapply(1:nRows, function(x)getSurvProb_nice2(data=data, testdata=test[x,],
                                                             fitObj=model, time=time, status=status,
                                                             manualCalc=calcPred,
                                                             timePoint=horizon,#,4,6
                                                             features=features,
                                                             predName=Modname,
                                                             strataVar=strataVar), simplify = FALSE)
  
  example.Surv.bw.dat<- sapply(example.Surv.bw, function(x) x %>% extract2(1), simplify=FALSE)
  
  example.Surv.bw.merged <- ldply(example.Surv.bw.dat)
  # Only anti-Pd-01
  
  example.Surv.bw.merged.f <- example.Surv.bw.merged %>% dplyr::filter(!!strata.un==selectStrat) %>% droplevels()
  # Female patient that had reveived Steroids, has no other malignancies, with PS=0/1
  example.Surv.bw.merged.f.present <-example.Surv.bw.merged.f %>% dplyr::select(-!!colsRemove)
  
  example.Surv.bw.merged.f.present <- example.Surv.bw.merged.f.present %>% 
              add_column(estim.surv=paste0(example.Surv.bw.merged.f.present$survival," (",example.Surv.bw.merged.f.present$survival.lower, 
                                           "-",example.Surv.bw.merged.f.present$survival.upper, ")"),.before="survival")
  # 
  example.Surv.bw.merged.f.present <- example.Surv.bw.merged.f.present %>% dplyr::select(-survival,-survival.lower, -survival.upper,-survival.se) %>% setNames(finalColsName)
  return(example.Surv.bw.merged.f.present)
  
}