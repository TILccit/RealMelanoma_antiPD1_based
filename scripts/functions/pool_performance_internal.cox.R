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
  #####################
  # data=boot_data_compl
  # nimp = nimp_mice
  # impvar=".imp"
  # formula = fm
  # horizon=45
  # time <- Y_time 
  # status <- Y_status
  # modelName="simple"
  #####################
  
  
  coef_f <- se_f <- pred.group <- obs.group <- list()
  
  if(is_empty(formula))
    stop("\n", "Model not specified in formula object")
  
  fm <-
    formula
  
  # perf_stats <-
  #   matrix(NA, nimp, 3)
  
  # Extract performance measures for EACH MI dataset
  set.seed(1234)
  clust <-parallel::makeCluster(5)
  parallel::clusterExport(clust,list("performanceCalculationMI","modelName","fm","horizon","coxph","Surv",
                                     "%>%","pluck","rcs"), envir = environment())
  all <- parallel::parSapply(clust, data %>% df2milist(impvar = impvar), function(x) performanceCalculationMI(name=modelName,formula=fm,data=x,horizon=horizon),simplify = FALSE)
  parallel::stopCluster(clust)
  
  perf_stats = ldply(all)
  
  # Get coefficients 
  coef_f <-sapply(data %>% df2milist(impvar = impvar), function(x) getCoefs(fm,x), simplify = FALSE)

  # C.index
  C.index <-data %>% df2milist(impvar = impvar) %>%
    with(cindex(coxph(formula=as.formula(as.character(fm))))) %>%
    pool_cindex()
  # Dxy 
  Dxy <-(C.index - 0.5)*2
  # End pooling performance measures in multiply imputed data
  
  # Pooling over measure values across MI datasets
  # Time dependent ROC (just average)
  roc_res <-
    pool_auc(perf_stats[, 3], sqrt(pROC::var(perf_stats[, 3])),
             nimp = nimp, log_auc = FALSE)
 
  # Pool Brier
  brier_pool <-
    mean(perf_stats[, 1])
  # Pool IPA
  ipa_pool <-
    mean(perf_stats[, 2])
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




pool_performance_internal.cox.OLD <- function(data, 
                                          formula, 
                                          nimp, 
                                          impvar)
{
  # #####################
  # data=boot_data_compl
  # nimp = nimp_mice
  # impvar=".imp"
  # formula = fm
  # horizon=45
  # time <- Y_time 
  # status <- Y_status
  # modelName="simple"
  # #####################
  
  
  coef_f <- se_f <- pred.group <- obs.group <- list()
  
  if(is_empty(formula))
    stop("\n", "Model not specified in formula object")
  
  fm <-
    formula
  
  perf_stats <-
    matrix(NA, nimp, 3)
  
  # lp_mi <- matrix(NA, nrow(data[data[impvar] == 1, ]), nimp)
  
  # # Pool performance measures over imputed datasets
  # system.time({ all<- map(
  #   boot_data_compl %>% df2milist(impvar = ".imp"),
  #   ~ performanceCalculationMI(name=modelName,formula=fm,data=.,horizon=horizon)
  # ) 
  # })
  # 
  # #  user  system elapsed 
  # #21.80    4.09   25.92  
  # 
  # system.time({ all<- sapply(
  #   boot_data_compl %>% df2milist(impvar = ".imp"),
  #   function(x) performanceCalculationMI(name=modelName,formula=fm,data=x,horizon=horizon)
  # ) 
  # })
  
  # user  system elapsed
  # 18.36    4.31   22.65
  
  
  set.seed(1234)
  clust <-parallel::makeCluster(5)
  parallel::clusterExport(clust,list("performanceCalculationMI","modelName","fm","horizon","coxph","Surv",
                                     "%>%","pluck","rcs"), envir = environment())
  all <- parallel::parSapply(clust, boot_data_compl %>% df2milist(impvar = ".imp"), function(x) performanceCalculationMI(name=modelName,formula=fm,data=x,horizon=horizon),simplify = FALSE)
  parallel::stopCluster(clust)
  
  perf_stats = ldply(all)
  
  
  
  coefs <-sapply(boot_data_compl %>% df2milist(impvar = ".imp"), function(x) getCoefs(fm,x), simplify = FALSE)
  
  # user  system elapsed 
  # 0.13    0.01   12.13 
  
  # system.time({
  # for (i in 1:nimp) {
  #   ##############################
  #   ## For each imputed dataset ##
  #   ##############################
  #   # i=1
  #   #######################
  #   # Extract imputed dataset
  #   data_compl <- data[data[impvar] == i, ]
  #   # Fit Cox model on imputed dataset- here use rms : MAKE A FUNCTION FOR THESE METRICS
  #   # f <-
  #   #   cph(fm, data = data_compl, x = T, y = T)
  #   # 
  #   f <-
  #     coxph(fm, data = data_compl, x = T, y = T)
  #   
  #   perf_stats[i, ] <-
  #     performanceCalculation(model=f,name=modelName,formula=fm,data=data_compl,horizon=horizon)
  #   # Brier, IPA, AUC # R2, g NOT FOR NOW, if you want to include them, fit the model with cph() of rms package
  #   
  #   # # 
  #   # coef_f[[i]] <-
  #   #   coef(f)
  #   # 
  # }
  # 
  # })  
  
  # user  system elapsed 
  # 20.38    4.19   24.57
  # C.index
  C.index <-boot_data_compl %>% df2milist(impvar = ".imp") %>%
    with(cindex(coxph(formula=as.formula(as.character(fm))))) %>%
    pool_cindex()
  # Dxy 
  Dxy <-(C.index - 0.5)*2
  # End pooling performance measures in multiply imputed data
  
  
  # Time dependent ROC (just average)
  roc_res <-
    pool_auc(perf_stats[, 3], sqrt(pROC::var(perf_stats[, 3])),
             nimp = nimp, log_auc = TRUE)
  # # Pooling R square
  # # Fisher z Transformation
  # rsq.n <-
  #   tanh(mean(atanh(perf_stats[, 1])))
  # # Pool gini
  # g_pool <-
  #   mean(perf_stats[, 2])
  # Pool Brier
  brier_pool <-
    mean(perf_stats[, 3])
  # Pool IPA
  ipa_pool <-
    mean(perf_stats[, 4])
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