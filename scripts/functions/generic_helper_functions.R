`%notin%` <- Negate(`%in%`)

detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

fmt_pvalue_with_stars <- function(x) {
  dplyr::case_when(
    x < 0.001 ~ paste0(style_pvalue(x), " ***"),
    x < 0.01 ~ paste0(style_pvalue(x), " **"),
    x < 0.05 ~ paste0(style_pvalue(x), " *"),
    TRUE ~ style_pvalue(x, digits = 3)
  )
}


modelFormula <- function(time, event, features){
  # time = surv_time
  # event = surv_status
  
  formula.f = as.formula(paste0("Surv(",time,",",event,")","~ ", features))
  formula.f
}


# Function that makes formula including interactions-provide vector of interactions
# for rcs provide vector of features, and then vector of knots
makeComplexFormula <- function(features, rcsFeats=NULL, knots=NULL, interactions=NULL,strataFeats=NULL,rms=FALSE,knotsPos=NULL){
  # features=c(categ_feats[-10],numeric_feats.ext[3]) # SHOULD NOT INCLUDE THE RCS FEATS/STRATA FEATS
  # rcsFeats=c(numeric_feats.ext[1],numeric_feats.ext[3])
  # knots=c(3,4)
  # interactions=list(c(categ_feats[1],categ_feats[3]),
  #                   # c(categ_feats[1],categ_feats[10]),
  #                   # c(categ_feats[5],rcsFeats[1]),
  #                   c(numeric_feats.ext[1],categ_feats[10]))
  # strataFeats=c(categ_feats[10])
  # features=c(categ_feats[-c(9)],numeric_feats.ext[3]) 
  # rcsFeats=c(numeric_feats.ext[1])
  # knots=c(3)
  # interactions=list(
  #   c(categ_feats[1],categ_feats[2]),
  #   c(categ_feats[1],categ_feats[3]),
  #   c(categ_feats[3],categ_feats[5]),
  #   c(categ_feats[5],numeric_feats.ext[3]),
  #   c(categ_feats[6],numeric_feats.ext[1]))
  # strataFeats=c(categ_feats[9])
  # knotsPos=list(Knots)
  # rms=F
  #############
  ## START
  if(any(rcsFeats %in% features) | any(strataFeats %in% features)){
    stop("Cannot include same feature in general features and rcs/strata feats!")
  }
  
  featsFinal<- c(features)
  # Begin with rcs
  if(!is.null(rcsFeats)){
    rcs_cov =c()
    
    for (i in 1:length(rcsFeats)){
      # i=1
      # print(i)
      if(!is.null(knotsPos)){
        rcs_cov = append(rcs_cov,
                         paste0("rcspline.eval(",rcsFeats[i],",knots = c(",knotsPos[[i]] %>% as.character() %>% paste(collapse = ","),"), inclx =TRUE)"))
        # rcs_cov = append(rcs_cov, 
        #                 paste0("rcspline.eval(",rcsFeats[i],",knots = c(",knotsPos[[i]] %>% as.character() %>% paste(collapse = ","),"))"))
      }else{
        rcs_cov = append(rcs_cov, 
                         paste0("rcs(",rcsFeats[i],",",knots[i],")"))
      }
      
    }
    
    featsFinal<- append(featsFinal,rcs_cov)
  }
  
  # Continue with strata
  if(!is.null(strataFeats)){
    strata_cov =c()
    
    for (i in 1:length(strataFeats)){
      # i=1
      # print(i)
      strata_cov = append(strata_cov, 
                          paste0("strata(",strataFeats[i],")"))
    }
    
    featsFinal<- append(featsFinal,strata_cov)
    
  }
  
  # Continue with interactions
  if(!is.null(interactions)){
    inter_cov =c()
    
    for (i in 1:length(interactions)){
       # i=5
      # print(i)
      ## Need to maked sure that if I have strata/rcs this need to be considered
      ## If strata vars are included in interactions
      
      #which(strataFeats %in% interactions[[2]])
      if(!is_empty(which(interactions[[i]]  %in% strataFeats))){
        str.idx <-which(interactions[[i]]  %in% strataFeats)
        interactions[[i]][str.idx] <- paste0("strata(",interactions[[i]][str.idx],")")
        
      }
      
      if(!is_empty(which(interactions[[i]]  %in% rcsFeats))){
        rcs.idx <-which(interactions[[i]]  %in% rcsFeats)
        rcs.knots.idx <-which(rcsFeats %in% interactions[[i]])
        
        if(!is.null(knotsPos)){
          # rcs_cov = append(rcs_cov,
          #                  paste0("rcspline.eval(",rcsFeats[i],",knots = c(",knotsPos[[i]] %>% as.character() %>% paste(collapse = ","),"), inclx =TRUE)"))
          # 
          interactions[[i]][rcs.idx] <- paste0("rcspline.eval(",interactions[[i]][rcs.idx],",","knots = c(",knotsPos[[rcs.knots.idx]] %>% as.character() %>% paste(collapse = ","),"), inclx =TRUE)")
                                               
        }else{
          interactions[[i]][rcs.idx] <- paste0("rcs(",interactions[[i]][rcs.idx],",",
                                               knots[rcs.knots.idx],")")
        }
        
        
        
      }
      
      
      inter_cov = append(inter_cov, 
                         paste(interactions[[i]], collapse = "*"))
    }
    
    featsFinal<- append(featsFinal,inter_cov)
    
  }
  
  form = paste(featsFinal, collapse = " + ")
  if(isTRUE(rms)){
    form <- str_replace_all(form,"strata","strat")
  }
  form
  
}


myclean_P <- function(variable) {
  variable <- gsub(paste(c("factor", "[()]", "rcs", "strata","strat",
                           ",", " ", c(3:7)), collapse = "|"), 
                   "", variable)
  return(variable)
}

## POST PROCESSING
fixScoresDF <- function(DF){
  # DF <- resCox_MI_vs_complete.noiter[[1]]$stats_val$`6`
  # DF <- resCox_MI_vs_complete.noiter[[1]]$stats_val$`6`
  # DF2 <- resCox_MI_vs_complete.noiter[[2]]$stats_val$`6`

  ########
  DF.t <- DF %>% dplyr::select(Index,Corrected) %>% column_to_rownames("Index") %>% data.table::transpose()
  rownames(DF.t) <- DF %>% dplyr::select(Index,Corrected) %>% column_to_rownames("Index") %>% colnames()
  colnames(DF.t) <- DF %>% dplyr::select(Index,Corrected) %>% column_to_rownames("Index") %>% rownames()
  # DF.t$model <- rep(modelNames,nrow(DF.t))
  DF.t
  
  
}


postProcessing_scores<- function(bootMIresObj,modelName){
  # bootMIresObj <-resCox_MI_vs_complete.noiter[[1]]
  # modelName="simple"
  ################
  
  scores.DF <-sapply(bootMIresObj$stats_val, function(x) fixScoresDF(x) , simplify=FALSE)
  indexes.DF <-ldply(scores.DF,.id="times")
  indexes.DF$model <- modelName
  indexes.DF
  
}


####################
#' @importFrom rms rcs asis pol catg scored matrx strat
#' @importFrom stats poly
#' @importFrom survival strata
trans.base2rms <- function(call.new){
  call.new=gsub('I\\(','asis(',call.new)
  call.new=gsub('ns\\(','rcs(',call.new)
  call.new=gsub('poly\\(','pol(',call.new)
  call.new=gsub('factor\\(','catg(',call.new)
  call.new=gsub('ordered\\(','scored(',call.new)
  call.new=gsub('matrix\\(','matrx(',call.new)
  call.new=gsub('strata\\(','strat(',call.new)
  call.new
}
trans.rms2base <- function(call.new){
  if (grepl('lsp\\(',call.new)) stop('We can not trans lsp() in rms to base function')
  call.new=gsub('asis\\(','I(',call.new)
  call.new=gsub('rcs\\(','ns(',call.new)
  call.new=gsub('pol\\(','poly(',call.new)
  call.new=gsub('catg\\(','factor(',call.new)
  call.new=gsub('scored\\(','ordered(',call.new)
  call.new=gsub('matrx\\(','matrix(',call.new)
  call.new=gsub('strat\\(','strata(',call.new)
  call.new
}


coxph2cph.v2 <- function(fit,dataIn){
  #   fit <- coxph(Surv(mpg,vs)~am+gear+strata(carb),data=mtcars)
  # coxph2cph(fit)
  #   
  if (class(fit) != 'coxph') stop('fit must be coxph results()')
  fit=update(fit,x=TRUE,y=TRUE,model=TRUE)
  call=paste0(deparse(fit$call),collapse = '')
  call.new=sub('coxph','cph',call)
  call.new=trans.base2rms(call.new)
  data.name=fit$call$data
  if (!is.null(data.name)){
    fit$model$timeggg=as.numeric(fit$model[,1])[1:nrow(fit$model)]
    fit$model$eventggg=as.numeric(fit$model[,1])[-c(1:nrow(fit$model))]
    colnames(fit$model)[(ncol(fit$model)-1):ncol(fit$model)]=        strsplit(do::Replace0(call,c('.*formula = Surv\\(','\\) ~.*')),', ')[[1]]
    fit$model=fit$model[,-1]
    # data=paste0(deparse(data.name),'=','fit$model')
    data=paste0(deparse(data.name),'=','dataIn')
    
    eval(parse(text = data))
    eval(parse(text=call.new), envir=.GlobalEnv)
    # eval(call.new)
    
  }else{
    stop('data must be given in formula')
  }
}



my.in_order <- function(x) {
  # droping unobserved value if needed
  forcats::fct_inorder(as.character(x))
}


# Make fits of two way interactions

analyseCox.interTwoWay <- function(data,featsVector, rms=T) {
  # data=data.complete
  # featsVector=featsCombos_final3[[11]]
  # print(featsVector)
  #rms=T
  cat <- NULL
  for(i in 1:length(featsVector)){
    i=1
    if(!is.numeric(data[[myclean_P(featsVector[i])]])){# it's categorical
      # print("CAT factor")
      # data[[featsVector[i]]] <- as.numeric(data[[featsVector[i]]])-1
      cat<- i
    }
  }
  
  
  featsVector <- c(setdiff(featsVector,featsVector[cat]),featsVector[cat])
  
  if(isTRUE(rms)){
    fit <-cph(formula = as.formula(paste0("Surv(",surv_time,",",surv_status,")","~ ",paste(featsVector, collapse = "*"))), data = data)# %>% coef()

  }else{
    fit <-coxph(formula = as.formula(paste0("Surv(",surv_time,",",surv_status,")","~ ",paste(featsVector, collapse = "*"))), data = data, model=TRUE)# %>% coef()
    
  }
  
  fit
  
  list(cat=cat, fit=fit)
}



getInteractions <- function(data, featsVector,typePlot="2D"){
  # data=real_data
  # featsVector=featsCombos[[1]]
  # typePlot="2D"
  
  res <-analyseCox.interTwoWay(data,featsVector)
  #Now plot
  featsVector <- c(setdiff(featsVector,featsVector[res$cat]),featsVector[res$cat])
  
  p <-mbrPlotInteraction(res$fit, vars=featsVector,plotype=typePlot)
  p <- recordPlot()
  # p <- p %>% as.grob()
  p
}


myplot_model <- function(model,data, featsVector) {
  # data=real_data
  # featsVector=featsCombos[[1]]
  cat <- NULL
  for(i in 1:length(featsVector)){
    if(!is.numeric(data[[featsVector[i]]])){# it's categorical
      data[[featsVector[i]]] <- as.numeric(data[[featsVector[i]]])-1
      cat<- i
    }
  }
  
  
  featsVector <- c(setdiff(featsVector,featsVector[cat]),featsVector[cat])
  
  plot_model(model, type = "emm", terms = featsVector, show.values = T, colors =pal_uchicago("default")(9) ) %>% print()
}


explorePH.inter.trans <- function(data, time, status, features, inters){
  # data=data.complete
  # time=surv_time
  # status=surv_status
  # features=feats1
  # inters=featsInters1
  
  # Fit model with rms
  dd <- datadist ( data )
  options ( datadist = dd)
  units ( data[[time]] ) <- 'Month'
  
  f <- cph(modelFormula(time,status,features), data = data, x=T, y=T)
  # print (f , coefs = FALSE )
  # Extract:
  LR.chi.2 <- f$stats[[3]]
  LR.p <-f$stats[[5]]
  dfs <- f$stats[[4]]
  R.2 <- f$stats[[8]]
  Dxy <- f$stats[[13]]
  
  # Calculate
  AIC.chi.2 = (LR.chi.2-2*dfs)
  AIC <- AIC(f)
  rough.Shrinkage=(LR.chi.2-dfs)/LR.chi.2
  
  model.stats <- c(LR.chi.2,LR.p,dfs,R.2,Dxy,AIC.chi.2,AIC,rough.Shrinkage)
  # Run anova: test for the significance of main effects
  anov.res <- anova(f)
  # aov1 <- aov(modelFormula(time,status,features), data = data)
  
  # Test PH
  phtest <- cox.zph (f , transform = 'identity')
  global.ph <- phtest$table["GLOBAL","p"]
  phtest.sign <- phtest$table %>% as.data.frame() %>% dplyr::filter(p< 0.05) %>% rownames() %>% paste(collapse=", " )
  ph.plots <- ggcoxzphFixed(phtest, resid=T, se=T)
  ### Interactions
  f.i <- cph(modelFormula(time,status,inters), data = data, x=T, y=T)
  # print (f.i , coefs = FALSE )
  # Extract:
  LR.chi.2.i <- f.i$stats[[3]]
  LR.p.i <-f.i$stats[[5]]
  dfs.i <- f.i$stats[[4]]
  R.2.i <- f.i$stats[[8]]
  Dxy.i <- f.i$stats[[13]]
  
  # Calculate
  AIC.chi.2.i = (LR.chi.2-2*dfs.i)
  AIC.i <- AIC(f.i)
  rough.Shrinkage.i=(LR.chi.2.i-dfs.i)/LR.chi.2.i
  
  model.stats.i <- c(LR.chi.2.i,LR.p.i,dfs.i,R.2.i,Dxy.i,AIC.chi.2.i,AIC.i,rough.Shrinkage.i)
  
  # Run anova: test for the significance of main effects
  anov.res.i <- anova(f.i)
  res.anova.i <- anov.res.i %>% 
    as.data.frame() %>% 
    dplyr::filter(P<0.05) %>% 
    rownames_to_column(var="term") %>% 
    pull(term) #%>% 
    # str_replace("\\...Factor.Higher.Order.Factors.","") %>% 
    # str_replace_all("\\...",":")
  
  x <- as.data.frame(res.anova.i)
  colnames(x) <- "term"
  x$term <-x$term %>% str_replace_all("...Factor.Higher.Order.Factors.","")
  x <- x %>% dplyr::filter(!grepl("*X.Nonlinear*|X.All.Interactions*|TOTAL*",term))
  x$term <-x$term %>% str_replace("\\.\\.+",":")
  x <- x %>% dplyr::filter(grepl("[A-Za-z0-9]+:[a-zA-Z]",term))
  # x %>% pull(term)
  
  # Test PH
  phtest.i <- cox.zph (f.i , transform = 'identity')
  global.ph.i <- phtest.i$table["GLOBAL","p"]
  phtest.sign.i <- phtest.i$table %>% as.data.frame() %>% dplyr::filter(p< 0.05) %>% rownames() %>% paste(collapse=", " )
  ph.plots.i <- ggcoxzphFixed(phtest.i, resid=T, se=T)
  
  # Gather results
  model.stat.all =rbind(model.stats,model.stats.i)
  colnames(model.stat.all)=c("LR.chi.2","LR.p","dfs","R.2","Dxy","AIC.chi.2","AIC","rough.Shrinkage")
  rownames(model.stat.all)=c("No interactions","With interactions")
  res<- list(model_stats= model.stat.all,
             anova=list(nointer=anov.res,
                        inter=anov.res.i),
             phtest=list(nointer=phtest,
                         inter=phtest.i),
             phtest.sign=list(nointer=phtest.sign,
                              inter=phtest.sign.i),
             ph.plots=list(nointer=ph.plots,
                           inter=ph.plots.i),
             anova.inter.res=x %>% pull(term))
  
}


getPooledFitObj <- function(dataToFit, dataMI, time, status, features, pooledObj){
  # dataToFit =orig.imp.list.cox.s
  # dataMI = orig.imp.list.cox
  # time = surv_time
  # status = surv_status
  # features = feats.trans.strata
  # pooledObj = pool_coxr.cox.nointer
  # 
  ## START
  # Fit cox on all MIs and pool
  fit.mi <- with(dataMI,coxph(formula = as.formula(modelFormula(time,status,features) %>% as.character)))
  fit.mi.pooled <- pool(fit.mi)   
  
  # out.res <-
  #   suppressWarnings(summary(pool(fit)))
  # Create fit cox obj to replace values with pooled
  fit.pool<- coxph(modelFormula(time,status,features ), data = dataToFit, x=T, y=T, model=T)
  # fit.pooled <- fit.strata.nointer.nobw
  # Add final pooled coefs to the model
  coef_names <-pooledObj %>%
    pluck("RR_model_final") %>%
    magrittr::extract2(1) %>%
    pull(term) %>% as.character()
  
  fit.pool$coefficients <-pooledObj %>%
    pluck("RR_model_final") %>%
    magrittr::extract2(1) %>%
    pull(estimate)  %>% as.vector() %>%
    setNames(coef_names)
  
  fit.pool$linear.predictors  <- pooledObj$pool_lp_final %>%  magrittr::extract2(1) %>% setNames(NULL)
  
  fit.pool$means <- pooledObj$pool_covCenter_final %>%
    magrittr::extract2(1)
  
  comb <- (MIcombine(fit.mi))
  # pooled.coef <- comb$coefficients
  pooled.vars<- comb$variance
  colnames(pooled.vars)<-NULL
  rownames(pooled.vars) <- NULL
  
  # Average Residuals
  # res= sapply(fit.mi,residuals)
  # res.avg <- apply(res,1,mean)
  # # fit.pool$residuals <- res.avg
  # plot(fit.pool$residuals,res.avg)
  # res2 = sapply(fit.mi,residuals, simplify = FALSE)
  # res.pooled =(MIcombine(res2,pooled.vars))
  # plot(fit.pool$residuals,res.pooled$coefficients)
  # 
  resids<-MIextract(fit.mi,fun=residuals)
  resids.pooled=MIcombine(resids,pooled.vars)
  # plot(fit.pool$residuals,resids.pooled$coefficients)
  fit.pool$residuals <-resids.pooled$coefficients
  
  # # betas<-MIextract(fit.mi,fun=coef)
  # vcov_imp<-MIextract(fit.mi, fun=vcov)
  # vcov_ary <- lm2a(vcov_imp)
  # # uncertainty due to sample size.
  # sample_size <- apply(X = vcov_ary, MARGIN = c(1,2), FUN = mean)
  # # MIextract(fit.mi,fun=logLik)
  # #uncertainty due to missing data. Here we quantify how much the regression coefficients from each imputed dataset differ from the pooled value
  # # Last we have the uncertainty due to finite simulation
  #################
  fit.pool$var<- pooled.vars
  
  return(fit.pool)
}


getPooledFitObj.cph <- function(dataToFit, dataMI, time, status, features, pooledObj){
  # dataToFit=orig.imp.list.cox.s
  # dataMI = orig.imp.list.cox
  # time = surv_time
  # status = surv_status
  # features = feats.trans.strata.f
  # pooledObj = pool_coxr.cox.nointer
  
  ## START
  # Fit cox on all MIs and pool
  fit.mi <- with(dataMI,cph(formula = as.formula(modelFormula(time,status,features %>% str_replace_all("strata","strat")) %>% as.character )))
  fit.mi.pooled <- pool(fit.mi)                  
  # Create fit cox obj to replace values with pooled
  fit.pool<- coxph(modelFormula(time,status,features ), data = dataToFit, x=T, y=T, model=T)
  # fit.pooled <- fit.strata.nointer.nobw
  # Add final pooled coefs to the model
  coef_names <-pooledObj %>%
    pluck("RR_model_final") %>%
    magrittr::extract2(1) %>%
    pull(term) %>% as.character()
  
  fit.pool$coefficients <-pooledObj %>%
    pluck("RR_model_final") %>%
    magrittr::extract2(1) %>%
    pull(estimate)  %>% as.vector() %>%
    setNames(coef_names)
  
  fit.pool$linear.predictors  <- pooledObj$pool_lp_final %>%  magrittr::extract2(1) %>% setNames(NULL)
  
  fit.pool$means <- pooledObj$pool_covCenter_final %>%
    magrittr::extract2(1)
  
  comb <- (MIcombine(fit.mi))
  # pooled.coef <- comb$coefficients
  pooled.vars<- comb$variance
  colnames(pooled.vars)<-NULL
  rownames(pooled.vars) <- NULL
  # # betas<-MIextract(fit.mi,fun=coef)
  # vcov_imp<-MIextract(fit.mi, fun=vcov)
  # vcov_ary <- lm2a(vcov_imp)
  # # uncertainty due to sample size.
  # sample_size <- apply(X = vcov_ary, MARGIN = c(1,2), FUN = mean)
  # # MIextract(fit.mi,fun=logLik)
  # #uncertainty due to missing data. Here we quantify how much the regression coefficients from each imputed dataset differ from the pooled value
  # # Last we have the uncertainty due to finite simulation
  #################
  fit.pool$var<- pooled.vars
  
  return(fit.pool)
}

rcsKnotAIC_Cox <- function(DF, endopoint, event, predictor, nk=3,knots=NULL, nBoots=5){
  # DF <- data.original.ext.cox %>% na.omit()
  # endopoint <- OS_time
  # event <- OS_status
  # predictor <- "Age_first_Treatment"
  # nk=3
  # knots=knots_ldh_man
  # 
  predictor.un <- as.name(predictor)
  plotTitle = paste0(predictor, ": ", nk," knots")
  
  dd <- datadist(DF)
  options(datadist = dd)
  
  if(is.null(knots)){
    myFit <- rms::cph(as.formula(paste0("Surv(",endopoint,", ",event,") ~ ","rcs(",predictor,",",nk,")")),
                      data = DF, x = T, y = T, surv = T)
  }else{
    print("Specify knot positions")
    print(as.formula(paste0("Surv(",endopoint,", ",event,") ~ ","rcs(",predictor,",c(",paste(knots,collapse = ","),"))")))
    print(paste0("Number of knots:",nk))
    myFit <- rms::cph(as.formula(paste0("Surv(",endopoint,", ",event,") ~ ","rcs(",predictor,",c(",paste(knots,collapse = ","),"))")),
                      data = DF, x = T, y = T, surv = T)
  }
  
  
  # myBootFit <- bootcov(myFit, B=9999, pr=F, coef.reps=T, loglik=F)
  my.valid <- validate(myFit, method="boot", B=nBoots)
  # my.calib <- calibrate(myFit, method="boot", B=1000)
  
  if(is.null(knots)){
    knotsPos <-rcspline.eval(DF[[predictor]] %>% na.omit(), nk=nk, inclx=TRUE, knots.only = T)
    
  }else{
    knotsPos <-rcspline.eval(DF[[predictor]] %>% na.omit(), nk=nk,knots=knots, inclx=TRUE, knots.only = T)
    
  }
  
  
  print(knotsPos)
  p <- rms::Predict(myFit) %>% as.data.frame() %>% ggplot(aes(!!predictor.un,yhat)) +
    # geom_jitter(size=0.9) + 
    geom_line(color="red",size=.9) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + 
    ggtitle(plotTitle) + labs(y="log Relative Hazard", x=predictor)
  
  for( i in 1:length(knotsPos)){
    # print(i)
    p <- p + geom_vline(xintercept = knotsPos[i], linetype="dashed", 
                        color = "black", size=.8)
  }
  # p
  # AIC of model
  getAIC <- extractAIC(myFit)
  res = list(plot = p,
             AIC = getAIC[2],
             fit=myFit,
             bootVal=my.valid,
             knots=knotsPos)
  res
}