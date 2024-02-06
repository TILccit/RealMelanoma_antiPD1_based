sample_boot <- function(db, B) {
  # db = db
  # B = 10
  db_boot <- matrix(NA, nrow = nrow(db) * B, ncol = ncol(db))
  sample_row <- list()
  for (j in 1:B) {
    sample_row[[j]] <- sample(nrow(db), size = nrow(db), replace = TRUE)
  }
  sample_row <- unlist(sample_row)
  db_boot <- db[sample_row, ]
  db_boot$id_boot <- sort(rep(1:B, nrow(db)))
  db_boot <- db_boot %>%
    group_by(id_boot) %>%
    nest() %>%
    dplyr::rename(boot_data = data)
  return(db_boot)
}


boot_Validation_Cox <- function(modelNames, time,status,features, data, horizon,readyBoot=NULL,nBoot, seed){
  
  # Model is constructed in each bootstrap sample and performance is evaluated in both bootstrap (assessment, analysiss) and original cohort. 
  # The difference between the train and test of the bootstrap is called optimism. 
  # By substracting the optimism from the apparent performance
  # You get the corrected indexes, expected model performance for future patients
  # similar to derivation (development) cohort.
  # 
  # provide functionality to give already bootstrapped data, so that it can accomotate for example BOOT-MI (pooled), or BOOT-SI
  # # # # # # # Runs for one model at a time
  # modelNames=c("simple")
  # time =surv_time
  # status = surv_status
  # features = feats1
  # data = data.all$original
  # horizon = 35.9999
  # nBoot = 10
  # seed=123
  # readyBoot=NULL
  ## START  
  # # Fit model on original dataset
  # fit.orig <- cph(formula=modelFormula(time, status, features),
  #            data = data,surv=T,x=T,y=T
  # )
  
  ###########
  ## START
  # Create # bootstrap from original data with missings included
  if(is.null(readyBoot)){
    set.seed(seed)
    data_boot <- bootstraps(data,
                            times = nBoot,
                            apparent = FALSE)
  }else{
    data.boot <- readyBoot
  }
  
  
  # Fit models on each bootstrap analusis set, that has the same size as the original data
  # fit <- cph(formula=modelFormula(time, status, features),
  #                  data = analysis(data),surv=T,x=T,y=T
  # )
  
  cox_boot = map(
    data_boot$splits,
    ~ cph(modelFormula(time, status, features), data = analysis(.), x = T, y = T)
  )
  
  ## Now I want to extract performance measures
  # Evaluate performance on resampled data
  # Evaluate performance on original data
  
  # C index
  #compute the concordance over the restricted range ymin <= y <= ymax. (For survival data this is a time range.)
  Harrell_C_orig =
    map_dbl(cox_boot,
            ~ concordance(Surv(data[[time]], data[[status]])
                          ~ predict(.x, newdata = data),
                          reverse = TRUE,ymax =horizon)$concordance
    )
  
  Harrell_C_boot =
    map2_dbl(
      data_boot$splits, cox_boot,
      ~ concordance(Surv(analysis(.x)[[time]],analysis(.x)[[status]])
                    ~ predict(.y, newdata = analysis(.x)),
                    reverse = TRUE,ymax =horizon)$concordance
    )
  Harrell_C_boot.test =
    map2_dbl(
      data_boot$splits, cox_boot,
      ~ concordance(Surv(assessment(.x)[[time]],assessment(.x)[[status]])
                    ~ predict(.y, newdata = assessment(.x)),
                    reverse = TRUE,ymax =horizon)$concordance
    )
  
  # Somers D 
  Dxy_orig =
    map_dbl(Harrell_C_orig,
            ~ 2*(.x-0.5)
    )
  
  Dxy_boot =
    map_dbl(
      Harrell_C_boot,
      ~ 2*(.x-0.5)
    )
  Dxy_boot.test =
    map_dbl(Harrell_C_boot.test,
            ~ 2*(.x-0.5)
    )
  
  # Brier, IPA, AUC
  Brier.Auc_orig =
    map(cox_boot,
        ~ overallAssessment_Apparent(model=.x,modelNames=c("simple"), 
                                     modelformula =modelFormula(time,  status, features),
                                     formula=modelFormula(surv_time,  status, 1), 
                                     data=data, horizon=horizon, bootStrapping = T,bootName = "Apparent")
    )
  
  
  
  Brier.Auc_boot =
    map2(
      data_boot$splits, cox_boot,
      ~ overallAssessment_Apparent(model=.y,modelNames=c("simple"), 
                                   modelformula =modelFormula(time,  status, features),
                                   formula=modelFormula(surv_time,  status, 1), 
                                   data=analysis(.x), horizon=horizon, bootStrapping = T,bootName = "Train")
    )
  
  Brier.Auc_boot.test =
    map2(
      data_boot$splits, cox_boot,
      ~ overallAssessment_Apparent(model=.y,modelNames=c("simple"), 
                                   modelformula =modelFormula(time,  status, features),
                                   formula=modelFormula(surv_time,  status, 1), 
                                   data=assessment(.x), horizon=horizon, bootStrapping = T,bootName = "Test")
    )
  
  ## Take the differences between bootstrapped calculated indexes bootstraps and on original data to calculate optimism for each index
  # # Optimism version 2 is substracting boot test from train
  # Harrell_C_diff2 = map2_dbl(
  #   Harrell_C_boot, Harrell_C_boot.test,
  #   function(a, b) {
  #     a - b
  #   }
  # )
  # 
  # Dxy_diff2 = map2_dbl(
  #   Dxy_boot, Dxy_boot.test,
  #   function(a, b) {
  #     a - b
  #   }
  # )
  # 
  # Brier.Auc_diff2 = map2(
  #   Brier.Auc_boot, Brier.Auc_boot.test,
  #   function(a, b) {
  #     a - b
  #   }
  # ) %>% map(., function(x){x %>% as.data.frame() %>% setNames("optimism")} ) %>% bind_cols
  
  # Optimism 1 is substracting boot apparent from train  -this seems more similar to results of validate()
  
  Harrell_C_diff = map2_dbl(
    Harrell_C_boot, Harrell_C_orig,
    function(a, b) {
      a - b
    }
  )
  
  Dxy_diff = map2_dbl(
    Dxy_boot, Dxy_orig,
    function(a, b) {
      a - b
    }
  )
  
  Brier.Auc_diff = map2(
    Brier.Auc_boot, Brier.Auc_orig,
    function(a, b) {
      a - b
    }
  ) %>% map(., function(x){x %>% as.data.frame() %>% setNames("optimism")} ) %>% bind_cols
  
  # Take average of differences (optimism) for each index, and substract it from original index estimate to get bias corrected estimate
  # First create DF for original
  res <- overallAssessment_Apparent(modelNames=c("simple"), 
                                    model = fit.orig, modelformula =modelFormula(time,  status, features),
                                    formula=modelFormula(time, status, 1), 
                                    data=data, horizon=horizon) %>% setNames("index.orig")
  res$training <- c(mean(Dxy_boot),-0,-0,Brier.Auc_boot %>% bind_cols() %>% rowMeans() %>% as.vector(),mean(Harrell_C_boot)) # training is the assessment set of each bootstrap
  res$test <- c(mean(Dxy_orig),-0,-0,Brier.Auc_orig %>% bind_cols() %>% rowMeans() %>% as.vector(),mean(Harrell_C_orig)) # test here is the fitted bootstrapped model on original data
  
  res$optimism <- c(mean(Dxy_diff),-0,-0,Brier.Auc_diff %>% bind_cols() %>% rowMeans() %>% as.vector(),mean(Harrell_C_diff))
  res$index.corrected <- res$index.orig - res$optimism
  
  # res$test2 <- c(mean(Dxy_boot.test),-0,-0,Brier.Auc_boot.test %>% bind_cols() %>% rowMeans() %>% as.vector(),mean(Harrell_C_boot.test))
  # res$optimism2 <- c(mean(Dxy_diff2),-0,-0,Brier.Auc_diff2 %>% bind_cols() %>% rowMeans() %>% as.vector(),mean(Harrell_C_diff2))
  # res$index.corrected2 <- res$index.orig - res$optimism2
  #############################
  
  
  res
}





#' Calculate optimism-corrected bootstrap internal validation for AUC, Brier and Scaled Brier Score
#' @param db data to calculate the optimism-corrected bootstrap
#' @param B number of bootstrap sample (default 10)
#' @param time follow-up time
#' @param status indicator variable (0=censored, 1=event)
#' @param formula_model formula for the model (Cox model)
#' @param pred.time time horizon as predictor
#' @param formula_ipcw formula to calculate inverse probability censoring weighting
#'
#' @return
#'
#' @author Daniele Giardiello
#'
#' @examples


# pacman::p_load(
#   rio,
#   survival,
#   rms,
#   pec,
#   tidyverse,
#   timeROC,
#   riskRegression
# )

## Simplified version including only Brier,IPA,AUC,Cindex
boot_internal_validation <- function(modelNames, time,status,features, data, 
                                     horizon,readyBoot=NULL,nBoot, seed,forPlot=FALSE){
  
  # Model is constructed in each bootstrap sample and performance is evaluated in both bootstrap (assessment, analysiss) and original cohort. 
  # The difference between the train and test of the bootstrap is called optimism. 
  # By substracting the optimism from the apparent performance
  # You get the corrected indexes, expected model performance for future patients
  # similar to derivation (development) cohort.
  # 
  # provide functionality to give already bootstrapped data, so that it can accomotate for example BOOT-MI (pooled), or BOOT-SI
  # # # # # # Runs for one model at a time
  # modelNames=c("simple")
  # time =surv_time
  # status = surv_status
  # features = str_replace(feats.trans.strata,"strata","strat")
  # data = data.complete
  # horizon = myHorizons
  # readyBoot=NULL
  # nBoot = 5
  # seed=123
  # # 
  ###########
  ## START  
  # Fit model on original dataset
  fit.orig <- cph(formula=modelFormula(time, status, features),
                  data = data,x=T,y=T#surv=T,
  )
  # Create bootstraps
  if(is.null(readyBoot)){
    set.seed(seed)
    data_boot <- bootstraps(data,
                            times = nBoot,
                            apparent = FALSE)
  }else{
    data.boot <- readyBoot
  }
  
  
  # Fit models on each bootstrap analysis set, that has the same size as the original data
  
  cox_boot = map(
    data_boot$splits,
    ~ cph(modelFormula(time, status, features), data = analysis(.), x = T, y = T)
  )
  
  ## Now I want to extract performance measures
  
  # Brier, IPA, AUC
  ## TEST
  Brier.Auc_orig =
    map(cox_boot,
        ~ overallAssessment_Apparent(model=.x,modelNames=c("simple"), 
                                     modelformula =modelFormula(time,  status, features),
                                     formula=modelFormula(surv_time,  status, 1), 
                                     data=data, horizon=horizon,time=time, status=status, bootStrapping = T,bootName = "Apparent",simple = TRUE)
    )
  names(Brier.Auc_orig) <- paste0("Boot ",1:nBoot)
  # Over bootstraps
  perf_measures_TEST <- sapply(Brier.Auc_orig, function(x) ldply(x,.id = "times"), simplify=FALSE)
  perf_measures_TEST <- ldply(perf_measures_TEST,.id = "Boot")
  perf_measures_TEST<-perf_measures_TEST %>%
    group_split(times,.keep = FALSE)
  names(perf_measures_TEST) <- horizon
  perf_measures_TEST <- sapply(perf_measures_TEST, function(x) x %>% setNames(paste0(colnames(x),"_test")), simplify=FALSE)
  perf_measures_TEST <- sapply(perf_measures_TEST, function(x) x %>% column_to_rownames(var ="Boot_test" ), simplify=FALSE)
  
  ## APPARENT
  Brier.Auc_boot =
    map2(
      data_boot$splits, cox_boot,
      ~ overallAssessment_Apparent(model=.y,modelNames=c("simple"), 
                                   modelformula =modelFormula(time,  status, features),
                                   formula=modelFormula(surv_time,  status, 1), 
                                   data=analysis(.x), horizon=horizon,time=time, status=status, bootStrapping = T,bootName = "Train",simple = TRUE)
    )
  
  names(Brier.Auc_boot) <- paste0("Boot ",1:nBoot)
  
  perf_measures_APP <- sapply(Brier.Auc_boot, function(x) ldply(x,.id = "times"), simplify=FALSE)
  perf_measures_APP <- ldply(perf_measures_APP,.id = "Boot")
  perf_measures_APP<-perf_measures_APP %>%
    group_split(times,.keep = FALSE)
  names(perf_measures_APP) <- horizon
  perf_measures_APP <- sapply(perf_measures_APP, function(x)  x %>% setNames(paste0(colnames(x),"_app")), simplify=FALSE)
  perf_measures_APP <- sapply(perf_measures_APP, function(x) x %>% column_to_rownames(var ="Boot_app" ), simplify=FALSE)
  
  # Merge over bootstraps the apparent and test
  res_boot <-sapply(as.character(horizon), function(x) cbind(perf_measures_APP[[x]]["Dxy_app"],
                                                             perf_measures_TEST[[x]]["Dxy_test"],
                                                             perf_measures_APP[[x]]["Cindex_app"],
                                                             perf_measures_TEST[[x]]["Cindex_test"],
                                                             perf_measures_APP[[x]]["AUC_app"],
                                                             perf_measures_TEST[[x]]["AUC_test"],
                                                             perf_measures_APP[[x]]["Brier_app"],
                                                             perf_measures_TEST[[x]]["Brier_test"],
                                                             perf_measures_APP[[x]]["IPA_app"],
                                                             perf_measures_TEST[[x]]["IPA_test"]), simplify = FALSE)
  
  ## Take the differences between bootstrapped calculated indexes bootstraps (APPARENT) and on original data (TEST) 
  ## to calculate optimism for each index
  # Calculate optimism for EACH HORIZON
  
  Dxy_optimism <- sapply(horizon %>% as.character, function(x) res_boot[[x]]$Dxy_app -  res_boot[[x]]$Dxy_test, simplify=FALSE)
  Cindex_optimism <- sapply(horizon %>% as.character, function(x) res_boot[[x]]$Cindex_app -  res_boot[[x]]$Cindex_test, simplify=FALSE)
  AUC_optimism <- sapply(horizon %>% as.character, function(x) res_boot[[x]]$AUC_app -  res_boot[[x]]$AUC_test, simplify=FALSE)
  Brier_optimism <- sapply(horizon %>% as.character, function(x) res_boot[[x]]$Brier_app -  res_boot[[x]]$Brier_test, simplify=FALSE)
  IPA_optimism <- sapply(horizon %>% as.character, function(x) res_boot[[x]]$IPA_app -  res_boot[[x]]$IPA_test, simplify=FALSE)
  
  # Average optimism over B bootstrap samples
  message("\n", "Average optimism over B bootstrap samples","\n")
  res_boot_m <-sapply(horizon %>% as.character, function(x) colMeans(data.frame(res_boot[[x]],
                                                                                "Dxy_optimism"=Dxy_optimism[[x]], 
                                                                                "Cindex_optimism"=Cindex_optimism[[x]], 
                                                                                "AUC_optimism"=AUC_optimism[[x]],
                                                                                "Brier_optimism"=Brier_optimism[[x]],
                                                                                "IPA_optimism"=IPA_optimism[[x]]),na.rm = TRUE), simplify=FALSE)
  
  
  
  ## ORIGINAL
  # Calculate performance measures for original
  res_orig <- overallAssessment_Apparent(modelNames=c("simple"), 
                                         model = fit.orig, modelformula =modelFormula(time,  status, features),
                                         formula=modelFormula(time, status, 1), 
                                         data=data, horizon=horizon,time=time, status=status, simple=TRUE) #%>% setNames("Original")
  
  names(res_orig) <- horizon
  
  #########################
  # Correct indexes
  
  Dxy_corr <-sapply(horizon %>% as.character(), function(x) res_orig[[x]]$Dxy - res_boot_m[[x]]["Dxy_optimism"] %>% as.numeric(), simplify=FALSE)
  Dxy_val <-sapply(horizon %>% as.character(), function(x) c(res_orig[[x]]$Dxy, res_boot_m[[x]]["Dxy_app"],
                                                             res_boot_m[[x]]["Dxy_test"], res_boot_m[[x]]["Dxy_optimism"], Dxy_corr[[x]]), simplify=FALSE)
  
  Cindex_corr <-sapply(horizon %>% as.character(), function(x) res_orig[[x]]$Cindex - res_boot_m[[x]]["Cindex_optimism"]%>% as.numeric(), simplify=FALSE)
  Cindex_val <-sapply(horizon %>% as.character(), function(x) c(res_orig[[x]]$Cindex, res_boot_m[[x]]["Cindex_app"],
                                                                res_boot_m[[x]]["Cindex_test"], res_boot_m[[x]]["Cindex_optimism"], Cindex_corr[[x]]), simplify=FALSE)
  
  
  AUC_corr <-sapply(horizon %>% as.character(), function(x) res_orig[[x]]$AUC - res_boot_m[[x]]["AUC_optimism"]%>% as.numeric(), simplify=FALSE)
  AUC_val <-sapply(horizon %>% as.character(), function(x) c(res_orig[[x]]$AUC, res_boot_m[[x]]["AUC_app"],
                                                             res_boot_m[[x]]["AUC_test"], res_boot_m[[x]]["AUC_optimism"], AUC_corr[[x]]), simplify=FALSE)
  
  Brier_corr <-sapply(horizon %>% as.character(), function(x) res_orig[[x]]$Brier - res_boot_m[[x]]["Brier_optimism"]%>% as.numeric(), simplify=FALSE)
  Brier_val <-sapply(horizon %>% as.character(), function(x) c(res_orig[[x]]$Brier, res_boot_m[[x]]["Brier_app"],
                                                               res_boot_m[[x]]["Brier_test"], res_boot_m[[x]]["Brier_optimism"], Brier_corr[[x]]), simplify=FALSE)
  
  
  IPA_corr <-sapply(horizon %>% as.character(), function(x) res_orig[[x]]$IPA - res_boot_m[[x]]["IPA_optimism"]%>% as.numeric(), simplify=FALSE)
  IPA_val <-sapply(horizon %>% as.character(), function(x) c(res_orig[[x]]$IPA, res_boot_m[[x]]["IPA_app"],
                                                             res_boot_m[[x]]["IPA_test"], res_boot_m[[x]]["IPA_optimism"], IPA_corr[[x]]), simplify=FALSE)
  
  pobjval <-sapply(horizon %>% as.character(), function(x) as.data.frame(matrix(c(Dxy_val[[x]], Cindex_val[[x]], AUC_val[[x]], 
                                                                                  Brier_val[[x]], IPA_val[[x]]), 5, 5, byrow = T)), simplify=FALSE)
  
  
  
  
  pobjval <-lapply(pobjval, setNames, c("Orig", "Apparent", "Test", "Optimism", "Corrected"))
  pobjval <-sapply(pobjval, function(x) {rownames(x) <- c("Dxy", "Cindex", "AUC", "Brier","IPA"); x}, simplify=FALSE)
  
  pobjval <- sapply(pobjval, function(x) x %>% as.data.frame() %>% rownames_to_column(var="Index") ,simplify = FALSE)
  
  pobjval <- list(stats_val = pobjval, #intercept_test = res_boot_m[7],
                  res_boot = res_boot,
                  model_orig=modelFormula(time, status, features))
  return(pobjval)

}



bootstrap_cv <- function(db, B = 10,
                         time,
                         status,
                         formula_model,
                         pred.time,
                         formula_ipcw) {
  
  
  # db = data.dev
  # B = 10
  # time = surv_time
  # status = surv_status
  # formula_model = modelFormula(surv_time, surv_status, feats1)
  # formula_ipcw = modelFormula(surv_time, surv_status, 1)
  # pred.time = 35.9999
  # # 
  ################
  # frm_model <- as.formula(formula_model)
  # frm_ipcw <- as.formula(formula_ipcw)
  ###########
  ## START
  frm_model <- formula_model
  frm_ipcw <- formula_ipcw
  
  db$id <- 1:nrow(db)
  
  
  # Duplicate data
  db_ext <- db %>% slice(rep(row_number(), B))
  db_ext$.rep <- with(db_ext, ave(seq_along(id), id, FUN = seq_along)) # add an index identifying the replications
  
  db_tbl <- db_ext %>%
    group_by(.rep) %>%
    nest() %>%
    dplyr::rename(
      orig_data = data,
      id_boot = .rep
    )
  
  # Create bootstrap samples
  
  # Join original data and the bootstrap data in a nested tibble
  a <- sample_boot(db, B)
  # a
  # # 
  
  # ##### USING THE BOOTSTRAPS FUNCTION ### - 
  #now we get the same results as bootstrapping with overallAssessment--WHAT IS THE DIFFERENCE?- WRONG
  # 
  # data_boot <- bootstraps(db,
  #                         times = B,
  #                         apparent = FALSE)
  # 
  # a <-   sapply(data_boot$splits, function(x) as_tibble(x[[1]]), simplify = FALSE)
  # names(a) <- seq(1,B,1)
  # a <- ldply(a,.id="id_boot")
  # a$id_boot <- as.numeric(a$id_boot)
  # a <- a %>% dplyr::select(-id)
  # a<- a  %>%
  #   group_by(id_boot) %>%
  #   nest() %>%
  #   dplyr::rename(boot_data = data)
  
  b <- a %>% left_join(db_tbl)
  
  # Create optimism-corrected performance measures
  b <- b %>% mutate(
    cox_boot = map(
      boot_data,
      ~ coxph(frm_model, data = ., x = T, y = T)
    ),
    
    cox_apparent = map(
      orig_data,
      ~ coxph(frm_model, data = ., x = T, y = T)
    ),
    
    # Model fitted on apparent data, index calculated on apparent data
    Harrell_C_app =
      map2_dbl(
        orig_data, cox_apparent,
        ~ concordance(Surv(.x[[time]], .x[[status]])
                      ~ predict(.y, newdata = .x),
                      reverse = TRUE)$concordance
      ),
    
    # Model fitted on bootstrapped data, index calculated on apparent data
    Harrell_C_orig =
      map2_dbl(
        orig_data, cox_boot,
        ~ concordance(Surv(.x[[time]], .x[[status]])
                      ~ predict(.y, newdata = .x),
                      reverse = TRUE)$concordance
      ),
    # Model fitted on bootstrapped data, index calculated on bootstrapped data
    Harrell_C_boot =
      map2_dbl(
        boot_data, cox_boot,
        ~ concordance(Surv(.x[[time]],.x[[status]])
                      ~ predict(.y, newdata = .x),
                      reverse = TRUE)$concordance
      ),
    
    Harrell_C_diff = map2_dbl(
      Harrell_C_boot, Harrell_C_orig,
      function(a, b) {
        a - b
      }
    ),
    
    Uno_C_app =
      map2_dbl(
        orig_data, cox_apparent,
        ~ concordance(Surv(.x[[time]], .x[[status]])
                      ~ predict(.y, newdata = .x),
                      reverse = TRUE,
                      timewt = "n/G2")$concordance
      ),
    
    Uno_C_orig =
      map2_dbl(
        orig_data, cox_boot,
        ~ concordance(Surv(.x[[time]], .x[[status]])
                      ~ predict(.y, newdata = .x),
                      reverse = TRUE,
                      timewt = "n/G2")$concordance
      ),
    
    Uno_C_boot =
      map2_dbl(
        boot_data, cox_boot,
        ~ concordance(Surv(.x[[time]],.x[[status]])
                      ~ predict(.y, newdata = .x),
                      reverse = TRUE,
                      timewt = "n/G2")$concordance
      ),
    
    Uno_C_diff = map2_dbl(
      Uno_C_boot, Uno_C_orig,
      function(a, b) {
        a - b
      }
    ),
    
    AUC_app = map2_dbl(
      orig_data, cox_apparent,
      ~ timeROC(
        T = .x[[time]], delta = .x[[status]],
        marker = predict(.y, newdata = .x),
        cause = 1, weighting = "marginal", times = pred.time,
        iid = FALSE
      )$AUC[[2]]
    ),
    
    AUC_orig = map2_dbl(
      orig_data, cox_boot,
      ~ timeROC(
        T = .x[[time]], delta = .x[[status]],
        marker = predict(.y, newdata = .x),
        cause = 1, weighting = "marginal", times = pred.time,
        iid = FALSE
      )$AUC[[2]]
    ),
    
    AUC_boot = map2_dbl(
      boot_data, cox_boot,
      ~ timeROC(
        T = .x[[time]], delta = .x[[status]],
        marker = predict(.y, newdata = .x),
        cause = 1, weighting = "marginal", times = pred.time,
        iid = FALSE
      )$AUC[[2]]
    ),
    
    AUC_diff = map2_dbl(
      AUC_boot, AUC_orig,
      function(a, b) {
        a - b
      }
    ),
    
    # Briers on apparent data from model fitted on apparent data
    Score_app = map2(
      orig_data, cox_apparent,
      ~ Score(list("Cox" = .y),
              formula = frm_ipcw,
              data = .x, times = pred.time,
              cens.method = "ipcw",
              cens.model = "cox",
              metrics = "brier",
              summary = "ipa"
      )$Brier[[1]]
    ),
    
    Brier_app = map_dbl(Score_app, ~ .x$Brier[[2]]),
    IPA_app = map_dbl(Score_app, ~ .x$IPA[[2]]),
    
    # Briers on apparent data from model fitted on bootstrapped data
    Score_orig = map2(
      orig_data, cox_boot,
      ~ Score(list("Cox" = .y),
              formula = frm_ipcw,
              data = .x, times = pred.time,
              cens.method = "ipcw",
              cens.model = "cox",
              metrics = "brier",
              summary = "ipa"
      )$Brier[[1]]
    ),
    
    Brier_orig = map_dbl(Score_orig, ~ .x$Brier[[2]]),
    IPA_orig = map_dbl(Score_orig, ~ .x$IPA[[2]]),
    # Briers on bootstrapped data from model fitted on bootstrapped data
    Score_boot = map2(
      boot_data, cox_boot,
      ~ Score(list("Cox" = .y),
              formula = frm_ipcw,
              data = .x, times = pred.time,
              cens.method = "ipcw",
              cens.model = "cox",
              metrics = "brier",
              summary = "ipa"
      )$Brier[[1]]
    ),
    
    Brier_boot = map_dbl(Score_boot, ~ .x$Brier[[2]]),
    IPA_boot = map_dbl(Score_boot, ~ .x$IPA[[2]]),
    Brier_diff = map2_dbl(
      Brier_boot, Brier_orig,
      function(a, b) {
        a - b
      }
    ),
    
    IPA_diff = map2_dbl(
      IPA_boot, IPA_orig,
      function(a, b) {
        a - b
      }
    )
  )
  
  # b
  # Generate output
  # Apparent
  AUC_apparent <- b$AUC_app[1]
  Brier_apparent <- b$Brier_app[1]
  IPA_apparent <- b$IPA_app[1]
  Harrell_C_apparent <- b$Harrell_C_app[1]
  Uno_C_apparent <- b$Uno_C_app[1]
  
  res.app <- c(Brier_apparent,IPA_apparent,AUC_apparent, Harrell_C_apparent, Uno_C_apparent) %>% as.data.frame() %>% setNames("Apparent") %>% set_rownames(c("Brier","IPA","AUC","Cindex","Uno_C"))
  # Internal validation
  
  AUC_corrected <- b$AUC_app[1] - mean(b$AUC_diff)
  Brier_corrected <- b$Brier_app[1] - mean(b$Brier_diff)
  IPA_corrected <- b$IPA_app[1] - mean(b$IPA_diff)
  Harrell_C_corrected <- b$Harrell_C_app[1] - mean(b$Harrell_C_diff)
  Uno_C_corrected <- b$Uno_C_app[1] - mean(b$Uno_C_diff)
  res.int <- c(Brier_corrected,IPA_corrected,AUC_corrected, Harrell_C_corrected, Uno_C_corrected) %>% as.data.frame() %>% setNames("Internal_validation") %>% set_rownames(c("Brier","IPA","AUC","Cindex","Uno_C"))
  res<- cbind(res.app,res.int)
  
  # res <- c("AUC corrected" = AUC_corrected,
  #          "Brier corrected" = Brier_corrected,
  #          "IPA corrected" = IPA_corrected,
  #          "Harrell C corrected" = Harrell_C_corrected,
  #          "Uno C corrected" = Uno_C_corrected)
  return(res)
}


#' Bootstrap validation in Multiply Imputed datasets
#'
#' \code{boot_MI} Bootstrapping followed by Multiple Imputation for internal validation.
#'  Called by function \code{psfmi_perform}.
#'
#' @param pobj An object of class \code{pmods} (pooled models), produced by a previous
#'  call to \code{psfmi_lr}.
#' @param data_orig dataframe of original dataset that contains missing data.
#' @param nboot The number of bootstrap resamples, default is 10.
#' @param nimp_mice Numerical scalar. Number of multiple imputation runs.
#' @param p.crit A numerical scalar. P-value selection criterium used for backward
#'  or forward selection during validation. When set at 1, validation is done
#'  without variable selection.
#' @param direction The direction of predictor selection, "BW" is for backward selection and
#'  "FW" for forward selection.
#' @param miceImp Wrapper function around the \code{mice} function.
#' @param ...  Arguments as predictorMatrix, seed, maxit, etc that can be adjusted for
#'  the \code{mice} function.
#'
#' @details 
#'  This function bootstraps from the incomplete dataset and applies MI in each 
#'  bootstrap sample. The model that is selected by the \code{psfmi_lr} function is
#'  validated. When p.crit != 1, internal validation is conducted with variable selection.
#'  The performance measures in the multiply imputed bootstrap samples are tested in the 
#'  original multiply imputed datasets (pooled) to determine the optimism.  
#'
#' @seealso \code{\link{psfmi_perform}}
#' @author Martijn Heymans, 2020
#' @keywords internal 
#' 
#' @export
boot_MI.cox <- function(pobj, 
                        data_orig, 
                        nboot = 10, 
                        nimp_mice, 
                        p.crit, 
                        direction, 
                        miceImp, 
                        predictorMatrix,
                        seed,
                        horizon,
                        modelName,
                        miceMethod,
                        smcfcsMethod,
                        smcfcsFormula,
                        keep.predictors=NULL,
                        strata.predictors=NULL,
                        ...)
{
  # # ## Performance metrics for one time horizon, not multiple ones
  # # # # # # # ####
  # pobj= pobj
  # data_orig = data_orig
  # nboot = nboot
  # nimp_mice = nimp_mice
  # p.crit = p.crit
  # direction = direction
  # miceImp = miceImp
  # predictorMatrix = predictorMatrix
  # seed=seed
  # horizon=horizon
  # modelName=modelName
  # miceMethod=miceMethod
  # smcfcsMethod=smcfcsMethod
  # smcfcsFormula=smcfcsFormula
  # keep.predictors=keep.predictors
  # strata.predictors=strata.predictors
  # # # # # ####
  
  # bootstrap from original data with missings included
  boot_data_orig <-
    bootstraps(data_orig, times = nboot)
  
  boot_seq <-
    as.list(1:nboot)
  
  # Name variable for imputations
  impVar<- pobj$impvar
  impVar.un <- as.name(impVar)
  # use bootstrap data as input
  # Start bootstrap run
  # tic()
  opt_boot <-
    mapply(function(x, y) {
      # PROCESS BELOW IS FOR EACH BOOTSTRAP SAMPLE
      # # # # ####
      # x = boot_data_orig$splits[[1]]
      # y=boot_seq[[1]]
      # # # #####
      message("\n", "Boot ", y)
      
      # print(paste0("\n", "Boot ", y))
      
      
      #print("\n", "Boot ", y)
      x <- as.data.frame(x)
      
      # Multiply Impute bootstrap sample with mice
      # boot_data_imp <-
      #   miceImp(data=x, m=nimp_mice, ...)
      
      if(!is.null(miceMethod)){
        
        if(isTRUE( smcfcsMethod)){
          # If I want to run smcfcsMethod the methods for all covariates need to be defined!--MAYBE CREATE ERROR>?
          
          boot_mi_smc.pl=smcfcs.parallel(
            smcfcs_func="smcfcs",
            seed = 3456,
            n_cores = 10,
            originaldata = x,
            smtype="coxph",
            smformula=smcfcsFormula,
            method=miceMethod,m=nimp_mice,predictorMatrix = predictorMatrix,noisy=F
          )
          
          
          # test<-smcfcs(x, smtype="coxph",smformula=smcfcsFormula,
          #        method=miceMethod,m=nimp_mice,predictorMatrix = predictorMatrix)
          # 
          # 
          boot_mi_smc.pl <- imputationList(boot_mi_smc.pl$impDatasets)
          boot_mi_smc.pl <- boot_mi_smc.pl$imputations
          names(boot_mi_smc.pl) <- paste0("imp",1:10)
          # Add id to each 
          
          boot_data_compl <- sapply(boot_mi_smc.pl, function(x) add_column(.id=seq.int(nrow(x)),.data=x, .before=1), simplify=FALSE)
          # Imputed datasets stacked vertically
          boot_data_compl <- ldply(boot_data_compl,.id=".imp")
          boot_data_compl$.imp <- as.character(boot_data_compl$.imp)
          boot_data_compl$.imp <- str_replace(boot_data_compl$.imp,"imp","")
          boot_data_compl$.imp <- as.numeric(boot_data_compl$.imp)
          
        }else{
          boot_data_imp <-
            miceImp(DF=x, m=nimp_mice, seed=seed,predictorMatrix=as.matrix(predictorMatrix),meth=miceMethod, print=FALSE)
          # Imputed datasets stacked vertically
          boot_data_compl <-
            complete(boot_data_imp, action = "long", include = FALSE)
        }
        
        
      }else{
        boot_data_imp <-
          miceImp(DF=x, m=nimp_mice, seed=seed,predictorMatrix=as.matrix(predictorMatrix), print=FALSE)
        # Imputed datasets stacked vertically
        boot_data_compl <-
          complete(boot_data_imp, action = "long", include = FALSE)
      }
      
      
      
      
      
      # Change impvar inside stacked imputed df column name (it is the first),
      # so that it is the same as in the trained original multiply imputed datasets
      colnames(boot_data_compl)[1] <-impVar
      
      #######################################
      # Extract time, status features and model formula
      Y_boot.time <-pobj$time
      Y_boot.status <-pobj$status
      
      fm_boot <-
        modelFormula(Y_boot.time,Y_boot.status, paste(pobj$predictors_final, collapse = "+"))
      
      # Pool model in completed (MI stacked) bootstrap sample
      # if(direction=="BW"){
      #   message("\n", "Pooling model in bootstrap sample (MI stacked) ","\n")
      # }else{
      #   message("\n", "Pooling model in bootstrap sample (MI stacked), using backward selection! ","\n")
      # }
      message("\n", "Pooling model in bootstrap sample (MI stacked) ","\n")
      pool_model_cox <-
        psfmi_coxr.Fixed(formula = fm_boot, data =  boot_data_compl, nimp=nimp_mice, impvar = impVar,
                         p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                         method = pobj$method, direction = direction,
                         strata.predictors=pobj$strata.predictors) #  NOT CORRECT FOR BW
      
      # With added pooled LPs and covCenter
      
      
      # Choose predictors according to variable selection method (if chosen)
      if(p.crit!=1){
        if(direction=="BW")
          predictors_selected <-
            ifelse(pool_model_cox$predictors_out[nrow(pool_model_cox$predictors_out), ], 0, 1)
        #names_predictors_selected <- pool_model_cox$predictors_initial
        if(direction=="FW")
          predictors_selected <-
            pool_model_cox$predictors_in[nrow(pool_model_cox$predictors_in), ]
      } else {
        predictors_selected <-
          ifelse(pool_model_cox$predictors_out[nrow(pool_model_cox$predictors_out), ], 0, 1)
      }
      
      # Extract apparent pooled estimates LP (LP in each bootstrap sample)
      Y_time <-pobj$time
      Y_status <-pobj$status
      
      if(is_empty(pool_model_cox$predictors_final)) {
        pool_model_cox$predictors_final <- 1
        fm <-
          modelFormula(Y_time,Y_status, paste(pool_model_cox$predictors_final, collapse = "+"))
        
        lp_app_pooled <- 1
      } else {
        fm <-
          modelFormula(Y_time,Y_status, paste(pool_model_cox$predictors_final, collapse = "+"))
        # lp_app_pooled <-
        #   pool_model_cox$RR_model[[1]][, 2]
        ## EXTRACT POOLED MEASURES
        
        coef_names <-pool_model_cox %>% 
          pluck("RR_model") %>% 
          magrittr::extract2(1) %>%
          pull(term) %>% as.character()
        
        lp_app_pooled <-pool_model_cox %>% 
          pluck("RR_model") %>% 
          magrittr::extract2(1) %>% 
          pull(estimate) %>%
          setNames(coef_names)
        
        # dataSim <- boot_data_compl %>% dplyr::filter(!!impVar.un==1) %>% droplevels()
        # fit.pooled <- coxph(formula = as.formula(fm %>% as.character), data = dataSim, x=T,y=T)
        # fit.mi <- with(boot_data_imp %>% miceadds::datlist2Amelia() %>% pluck("imputations"),coxph(formula = as.formula(fm %>% as.character), x=T,y=T))
        # # fit.mi <- with(boot_data_compl %>% split(f=impVar),coxph(formula = as.formula(fm %>% as.character), x=T,y=T))
        # comb <- (MIcombine(fit.mi))
        # # pooled.coef <- comb$coefficients
        # pooled.vars<- comb$variance
        # colnames(pooled.vars)<-NULL
        # rownames(pooled.vars) <- NULL
        # pooled.coef <- pool_model_cox %>% 
        #   pluck("RR_model") %>% 
        #   magrittr::extract2(1) %>%
        #   pull(estimate)  %>%
        #   setNames(coef_names)
        # fit.pooled$coefficients<- pooled.coef
        # fit.pooled$linear.predictors <- pool_model_cox$pool_lp_final %>% 
        #   magrittr::extract2(1) %>% setNames(NULL)
        # fit.pooled$means <- pool_model_cox$pool_covCenter_final %>% 
        #   magrittr::extract2(1) 
        # fit.pooled$var <-pooled.vars
        
        
        
        #
        
        if(p.crit<1)
          coef_names <-pool_model_cox %>% 
          pluck("RR_model_final") %>% 
          magrittr::extract2(1) %>%
          pull(term) %>% as.character()
        
        lp_app_pooled <-pool_model_cox %>% 
          pluck("RR_model_final") %>% 
          magrittr::extract2(1) %>% 
          pull(estimate) %>%
          setNames(coef_names)
        
          # dataSim <- boot_data_compl %>% dplyr::filter(!!impVar.un==1) %>% droplevels()
          # fit.pooled <- coxph(formula = as.formula(pool_model_cox$formula_final %>% as.character), data = dataSim, x=T,y=T)
          # fit.mi <- with(boot_data_imp %>% miceadds::datlist2Amelia() %>% pluck("imputations"),coxph(formula = as.formula(pool_model_cox$formula_final %>% as.character), x=T,y=T))
          # comb <- (MIcombine(fit.mi))
          # # pooled.coef <- comb$coefficients
          # pooled.vars<- comb$variance
          # colnames(pooled.vars)<-NULL
          # rownames(pooled.vars) <- NULL
          # pooled.coef <- pool_model_cox %>% 
          #   pluck("RR_model_final") %>% 
          #   magrittr::extract2(1) %>%
          #   pull(estimate)  %>%
          #   setNames(coef_names)
          # fit.pooled$coefficients<- pooled.coef
          # fit.pooled$linear.predictors <- pool_model_cox$pool_lp_final
          # fit.pooled$means <- pool_model_cox$pool_covCenter_final
          # fit.pooled$var <-pooled.vars
        
      }
      
      # Obtain apparent pooled performance measures
      # tic()
      message("\n", "Obtain apparent pooled performance measures ","\n") #  TRAIN
      
      perform_app <-
        pool_performance_internal.cox(data=boot_data_compl, nimp = nimp_mice,
                                      impvar=impVar, formula = fm,
                                      horizon=horizon,
                                      time=Y_time,
                                      status=Y_status,
                                      modelName=modelName)
      
      # FROM ABOVE I GET WARNINGS AND ERRORS
      # Warning messages:
      #   1: In Score.list(modelList, formula = formula, data = data,  ... :
      #                      Cannot (not yet) estimate standard errors for AUC with Cox IPCW.
      #                    Therefore, force cens.model to be marginal.
      #                    2: In coxLP.coxph(object, data = newdata, center = FALSE) :
      #                      survival::predict.coxph returns the following error:
      #                      Error in predict.coxph(object, newdata = as.data.frame(data), type = "terms") : 
      #                      Data is not the same size as it was in the original fit
      #                    It seems that the dataset used to fit the model is no more compatible with the model,
      #                    probably because it has been modified afterwards.
      #                    coxLP.coxph will try to reconstruct the original dataset and continue the execution.
      
      # toc()
      
      # Test apparent LP in original Multiply Imputed data
      coef_slope_test <-
        list()
      perform_test <-
        matrix(NA, pobj$nimp, 5)
      #
      # Define train data for input in predicting risk in test data
      # They will be used to fit a model and help extract design matrix
      # that will be used in calculating baseline cum hazard
      # For each bootstrap we have multiple train data( because of MI)
      # so since -AS I UNDERSTAND- the train data are only used for the
      # design matrix, I will choose one imputation as train (in each bootstrap)
      
      
      boot_train <- boot_data_compl %>% dplyr::filter(!!impVar.un==1) %>% droplevels() 
      # This is used for predictions in the next step as train data
      # and it affects the model frame (matrix) and eventually the cum hazard and out predictions and performance measures
      # Here with the above I use the SAME one MI dataset as train (from the bootstrap)
      # Alternatively I could much MI data from train with MI from test,
      # So for each risk prediction in each original MI dataset, to use as train, one MI data from the train bootstrap, different each time.
      # ##################
      # fit1 <- cph(formula=str_replace(fm,"strata","strat") %>% as.formula,data=boot_data_compl %>% dplyr::filter(!!impVar.un==1) %>% droplevels(),x=TRUE, y=TRUE, surv=TRUE, time.inc=6)
      # fit2 <- cph(formula=str_replace(fm,"strata","strat") %>% as.formula,data=boot_data_compl %>% dplyr::filter(!!impVar.un==2) %>% droplevels(),x=TRUE, y=TRUE, surv=TRUE, time.inc=6)
      # 
      # #################
      
      ## Apply pooled model in original multiply imputed datasets and  obtain pooled AUC, Brier etc
      ### Following runs for all nested imputations within each bootstrap
      # tic()
      message("\n", "Apply pooled model(train) in original mutliply imputed datasets, get pooled performance measures","\n")
      # Here I get error on the 5th bootstrapm on this function (I guess>?)
      perform_test <-pool_performance_nestedBootMI.cox(data=pobj$data,
                                                       form=fm,
                                                       nimp=nimp_mice,
                                                       impvar=impVar,
                                                       horizon=horizon,
                                                       name=modelName,
                                                       coefs.app_pooled = lp_app_pooled,
                                                       traindata=boot_train,
                                                       lp.app_pooled=pool_model_cox$pool_lp_final  %>% magrittr::extract2(1),
                                                       covmeans.app_pooled=pool_model_cox$pool_covCenter_final %>% magrittr::extract2(1),
                                                       time=Y_time,
                                                       status = Y_status)
      # toc()
      
      # End pooling in original multiply imputed test data
      
      # Dxy
      Dxy_app <-perform_app$Dxy_pooled[1]
      Dxy_test <-perform_test$Dxy_pooled[1]
      # Cindex
      Cindex_app <-perform_app$C.index_pooled[1]
      Cindex_test <-perform_test$C.index_pooled[1]
      # timeROC
      # AUC_app <-perform_app$timeROC_pooled[2]
      # AUC_test <-perform_test$timeROC_pooled[2]
      AUC_app <-sapply(perform_app$timeROC_pooled, function(x) x[2],simplify = FALSE)
      AUC_test <-sapply(perform_test$timeROC_pooled, function(x) x[2],simplify = FALSE)
      # Brier
      Brier_app <-perform_app$Brier_pooled
      Brier_test <-perform_test$Brier_pooled
      # IPA
      IPA_app <-perform_app$IPA_pooled
      IPA_test <-perform_test$IPA_pooled
      
      
      perf.measures<-sapply(horizon %>% as.character(), function(x) c(Dxy_app, Dxy_test, 
                                                                      Cindex_app, Cindex_test,
                                                                      AUC_app[[x]], AUC_test[[x]],
                                                                      Brier_app[[x]], Brier_test[[x]],
                                                                      IPA_app[[x]],IPA_test[[x]]), simplify=FALSE)
      # Gather
      # opt_perform <-
      #   list(c(Dxy_app, Dxy_test, 
      #          Cindex_app, Cindex_test,
      #          AUC_app, AUC_test,
      #          Brier_app, Brier_test,
      #          IPA_app,IPA_test),
      #        predictors_selected)
      
      # list of vectors with performance measures at each timepoint
      opt_perform <-
        list(perf.measures,
             predictors_selected)
      
      
      return(opt_perform)
      
    }, x = boot_data_orig$splits, y=boot_seq, SIMPLIFY = FALSE)
  
  message("\n", "Completed calculating measures over bootstraps","\n")
  # toc()
  predictors_selected <-
    data.frame(do.call("rbind", lapply(opt_boot, function(x) x[[2]])))
  # If there are strata, they are not included in the predictors_selected, but are present in predictors_final
  # I have to find a way to workaround that.- Below works only for ONE STRATA PREDICTOR
  
  if(!is_empty(pobj$strata.predictors)){
    colnames(predictors_selected) <-
      setdiff(pobj$predictors_final, paste0("strata(",pobj$strata.predictors,")"))
  }else{
    colnames(predictors_selected) <-
      pobj$predictors_final
    
  }
  
  # colnames(predictors_selected) <-
  #   pobj$predictors_final
  row.names(predictors_selected) <-
    paste("Boot", 1:nboot)
  
  res_boot <-
    sapply(horizon %>% as.character, function(y) data.frame(do.call("rbind", lapply(opt_boot, function(x) x[[1]][[y]]))), simplify=FALSE)
  
  # Fix colnames in perf measures tables
  perfNames<- c("Dxy_app", "Dxy_test", 
                "Cindex_app", "Cindex_test",
                "AUC_app", "AUC_test",
                "Brier_app", "Brier_test",
                "IPA_app","IPA_test")
  
  res_boot <-lapply(res_boot, setNames, perfNames)
  
  res_boot <-sapply(res_boot, function(x) {rownames(x) <- paste("Boot", 1:nboot); x}, simplify=FALSE)
  # colnames(res_boot) <-
  #   c("Dxy_app", "Dxy_test", 
  #     "Cindex_app", "Cindex_test",
  #     "AUC_app", "AUC_test",
  #     "Brier_app", "Brier_test",
  #     "IPA_app","IPA_test")
  # row.names(res_boot) <-
  #   paste("Boot", 1:nboot)
  
  # Calculate optimism over perfermance measures for all bootstraps
  # Dxy_optimism <-
  #   res_boot$Dxy_app - res_boot$Dxy_test
  # Cindex_optimism <-
  #   res_boot$Cindex_app - res_boot$Cindex_test
  # AUC_optimism <-
  #   res_boot$AUC_app - res_boot$AUC_test
  # Brier_optimism <-
  #   res_boot$Brier_app - res_boot$Brier_test
  # IPA_optimism <-
  #   res_boot$IPA_app - res_boot$IPA_test
  
  Dxy_optimism <- sapply(horizon %>% as.character, function(x) res_boot[[x]]$Dxy_app -  res_boot[[x]]$Dxy_test, simplify=FALSE)
  Cindex_optimism <- sapply(horizon %>% as.character, function(x) res_boot[[x]]$Cindex_app -  res_boot[[x]]$Cindex_test, simplify=FALSE)
  AUC_optimism <- sapply(horizon %>% as.character, function(x) res_boot[[x]]$AUC_app -  res_boot[[x]]$AUC_test, simplify=FALSE)
  Brier_optimism <- sapply(horizon %>% as.character, function(x) res_boot[[x]]$Brier_app -  res_boot[[x]]$Brier_test, simplify=FALSE)
  IPA_optimism <- sapply(horizon %>% as.character, function(x) res_boot[[x]]$IPA_app -  res_boot[[x]]$IPA_test, simplify=FALSE)
  
  # Average optimism over B bootstrap samples
  message("\n", "Average optimism over B bootstrap samples","\n")
  res_boot_m <-sapply(horizon %>% as.character, function(x) colMeans(data.frame(res_boot[[x]],
                                                                                "Dxy_optimism"=Dxy_optimism[[x]], 
                                                                                "Cindex_optimism"=Cindex_optimism[[x]], 
                                                                                "AUC_optimism"=AUC_optimism[[x]],
                                                                                "Brier_optimism"=Brier_optimism[[x]],
                                                                                "IPA_optimism"=IPA_optimism[[x]]),na.rm = TRUE), simplify=FALSE)
  
  # res_boot_m <-
  #   colMeans(data.frame(res_boot,
  #                       Dxy_optimism, Cindex_optimism, AUC_optimism,Brier_optimism,IPA_optimism))
  # 
  # Perform original model in multiply imputed original data
  # if p.crit != 1 selection takes place during validation
  if(p.crit != 1) {
    Y_orig_time <- pobj$time
    Y_orig_status <- pobj$status
    if(is_empty(pobj$predictors_initial)) stop("You can not validate an empty model")
    # fm_orig <-
    #   as.formula(paste(Y_orig, paste(pobj$predictors_initial, collapse = "+")))
    # 
    
    # if(!is_empty(pobj$strata.predictors)){
    #   fm_orig <-
    #     modelFormula(Y_orig_time,Y_orig_status, paste(c(pobj$predictors_final, paste0("strata(",pobj$strata.predictors,")")), collapse = "+"))
    # }else{
    #   fm_orig <-
    #     modelFormula(Y_orig_time,Y_orig_status, paste(pobj$predictors_final, collapse = "+"))
    #   
    # }
    fm_orig <-
      modelFormula(Y_orig_time,Y_orig_status, paste(pobj$predictors_final, collapse = "+"))
    
    message("\n", "Pool model in original multiply imputed data, extract pooled original performance measures || Variable selection!","\n")
    
    pool_coxr_orig <-
      psfmi_coxr.Fixed(formula = fm_orig, data =  pobj$data, nimp=pobj$nimp,
                       impvar = pobj$impvar, p.crit = p.crit,
                       keep.predictors = pobj$keep.predictors,
                       method = pobj$method, direction = direction,
                       strata.predictors=pobj$strata.predictors)
    
    # psfmi_coxr.Fixed(data=mydata, nimp=5, impvar="Impnr", 
    #                 formula = modelFormula(surv_time, surv_status, feats1), p.crit=1,
    #                 method="D1", direction = "BW")
    
    if(is_empty(pool_coxr_orig$predictors_final)) {
      pool_coxr_orig$predictors_final <- 1
      # fm_orig <-
      #   as.formula(paste(Y_orig, paste(pool_coxr_orig$predictors_final, collapse = "+")))
      
      # if(!is_empty(pobj$strata.predictors)){
      #   fm_orig <-
      #     modelFormula(Y_orig_time,Y_orig_status, paste(c(pobj$predictors_final, paste0("strata(",pobj$strata.predictors,")")), collapse = "+"))
      # }else{
      #   fm_orig <-
      #     modelFormula(Y_orig_time,Y_orig_status, paste(pobj$predictors_final, collapse = "+"))
      #   
      # }
      fm_orig <-
        modelFormula(Y_orig_time,Y_orig_status, paste(pobj$predictors_final, collapse = "+"))
      
    } else {
      # # fm_orig <-
      # #   as.formula(paste(Y_orig, paste(pool_coxr_orig$predictors_final, collapse = "+")))
      # fm_orig <-
      #   modelFormula(Y_orig_time,Y_orig_status, paste(pobj$predictors_final, collapse = "+"))
      
      # if(!is_empty(pobj$strata.predictors)){
      #   fm_orig <-
      #     modelFormula(Y_orig_time,Y_orig_status, paste(c(pobj$predictors_final, paste0("strata(",pobj$strata.predictors,")")), collapse = "+"))
      # }else{
      #   fm_orig <-
      #     modelFormula(Y_orig_time,Y_orig_status, paste(pobj$predictors_final, collapse = "+"))
      #   
      # }
      fm_orig <-
        modelFormula(Y_orig_time,Y_orig_status, paste(pobj$predictors_final, collapse = "+"))
      
    }
    message("\n", "Extract pooled original performance measures","\n")
    
    perform_mi_orig <-
      pool_performance_internal.cox(data=pobj$data, nimp = pobj$nimp,
                                    impvar=pobj$impvar, formula = fm_orig,
                                    horizon=horizon,
                                    time=Y_time,
                                    status=Y_status,
                                    modelName=modelName)
    
    
    
  } else {
    # Y_orig <- c(paste(pobj$Outcome, paste("~")))
    Y_orig_time <- pobj$time
    Y_orig_status <- pobj$status
    
    if(is_empty(pobj$predictors_final)) {
      pobj$predictors_final <- 1
      # # fm_orig <-
      # #   as.formula(paste(Y_orig, paste(pobj$predictors_final, collapse = "+")))
      # fm_orig <-
      #   modelFormula(Y_orig_time,Y_orig_status, paste(pobj$predictors_final, collapse = "+"))
      
      # if(!is_empty(pobj$strata.predictors)){
      #   fm_orig <-
      #     modelFormula(Y_orig_time,Y_orig_status, paste(c(pobj$predictors_final, paste0("strata(",pobj$strata.predictors,")")), collapse = "+"))
      # }else{
      #   fm_orig <-
      #     modelFormula(Y_orig_time,Y_orig_status, paste(pobj$predictors_final, collapse = "+"))
      #   
      # }
      # 
      fm_orig <-
        modelFormula(Y_orig_time,Y_orig_status, paste(pobj$predictors_final, collapse = "+"))
      
    } else {
      # # fm_orig <-
      # #   as.formula(paste(Y_orig, paste(pobj$predictors_final, collapse = "+")))
      # fm_orig <-
      #   modelFormula(Y_orig_time,Y_orig_status, paste(pobj$predictors_final, collapse = "+"))
      
      # if(!is_empty(pobj$strata.predictors)){
      #   fm_orig <-
      #     modelFormula(Y_orig_time,Y_orig_status, paste(c(pobj$predictors_final, paste0("strata(",pobj$strata.predictors,")")), collapse = "+"))
      # }else{
      #   fm_orig <-
      #     modelFormula(Y_orig_time,Y_orig_status, paste(pobj$predictors_final, collapse = "+"))
      #   
      # }
      # 
      fm_orig <-
        modelFormula(Y_orig_time,Y_orig_status, paste(pobj$predictors_final, collapse = "+"))
    }
    message("\n", "Pool model in original multiply imputed data, extract pooled original performance measures","\n")
    perform_mi_orig <-
      pool_performance_internal.cox(data=pobj$data, nimp = pobj$nimp,
                                    impvar=pobj$impvar, formula = fm_orig,
                                    horizon=horizon,
                                    time=Y_time,
                                    status=Y_status,
                                    modelName=modelName)
  }
  
  # Extract performance measures
  # Dxy
  Dxy_orig <-perform_mi_orig$Dxy_pooled[1]
  # Cindex
  Cindex_orig <-perform_mi_orig$C.index_pooled[1]
  # timeROC
  AUC_orig <-sapply(perform_mi_orig$timeROC_pooled, function(x) x[2],simplify = FALSE)
  # Brier
  Brier_orig <-perform_mi_orig$Brier_pooled
  # IPA
  IPA_orig <-perform_mi_orig$IPA_pooled
  
  
  
  # Correct indexes
  # Dxy_corr <-
  #   Dxy_orig - res_boot_m["Dxy_optimism"]
  # 
  # Dxy_val <-
  #   c(Dxy_orig, res_boot_m["Dxy_app"],
  #     res_boot_m["Dxy_test"], res_boot_m["Dxy_optimism"], Dxy_corr)
  # 
  Dxy_corr <-sapply(horizon %>% as.character(), function(x) Dxy_orig - res_boot_m[[x]]["Dxy_optimism"], simplify=FALSE)
  Dxy_val <-sapply(horizon %>% as.character(), function(x) c(Dxy_orig, res_boot_m[[x]]["Dxy_app"],
                                                             res_boot_m[[x]]["Dxy_test"], res_boot_m[[x]]["Dxy_optimism"], Dxy_corr[[x]]), simplify=FALSE)
  
  # Cindex_corr <-
  #   Cindex_orig - res_boot_m["Cindex_optimism"]
  # Cindex_val <-
  #   c(Cindex_orig, res_boot_m["Cindex_app"],
  #     res_boot_m["Cindex_test"], res_boot_m["Cindex_optimism"], Cindex_corr)
  
  Cindex_corr <-sapply(horizon %>% as.character(), function(x) Cindex_orig - res_boot_m[[x]]["Cindex_optimism"], simplify=FALSE)
  Cindex_val <-sapply(horizon %>% as.character(), function(x) c(Cindex_orig, res_boot_m[[x]]["Cindex_app"],
                                                                res_boot_m[[x]]["Cindex_test"], res_boot_m[[x]]["Cindex_optimism"], Cindex_corr[[x]]), simplify=FALSE)
  
  
  # AUC_corr <-
  #   AUC_orig - res_boot_m["AUC_optimism"]
  # AUC_val <-
  #   c(AUC_orig, res_boot_m["AUC_app"],
  #     res_boot_m["AUC_test"], res_boot_m["AUC_optimism"], AUC_corr)
  
  
  AUC_corr <-sapply(horizon %>% as.character(), function(x) AUC_orig[[x]] - res_boot_m[[x]]["AUC_optimism"], simplify=FALSE)
  AUC_val <-sapply(horizon %>% as.character(), function(x) c(AUC_orig[[x]], res_boot_m[[x]]["AUC_app"],
                                                             res_boot_m[[x]]["AUC_test"], res_boot_m[[x]]["AUC_optimism"], AUC_corr[[x]]), simplify=FALSE)
  
  
  # Brier_corr <-
  #   Brier_orig - res_boot_m["Brier_optimism"]
  # Brier_val <-
  #   c(Brier_orig, res_boot_m["Brier_app"],
  #     res_boot_m["Brier_test"], res_boot_m["Brier_optimism"], Brier_corr)
  # 
  
  Brier_corr <-sapply(horizon %>% as.character(), function(x) Brier_orig[[x]] - res_boot_m[[x]]["Brier_optimism"], simplify=FALSE)
  Brier_val <-sapply(horizon %>% as.character(), function(x) c(Brier_orig[[x]], res_boot_m[[x]]["Brier_app"],
                                                               res_boot_m[[x]]["Brier_test"], res_boot_m[[x]]["Brier_optimism"], Brier_corr[[x]]), simplify=FALSE)
  
  
  # IPA_corr <-
  #   IPA_orig - res_boot_m["IPA_optimism"]
  # IPA_val <-
  #   c(IPA_orig, res_boot_m["IPA_app"],
  #     res_boot_m["IPA_test"], res_boot_m["IPA_optimism"], IPA_corr)
  
  IPA_corr <-sapply(horizon %>% as.character(), function(x) IPA_orig[[x]] - res_boot_m[[x]]["IPA_optimism"], simplify=FALSE)
  IPA_val <-sapply(horizon %>% as.character(), function(x) c(IPA_orig[[x]], res_boot_m[[x]]["IPA_app"],
                                                             res_boot_m[[x]]["IPA_test"], res_boot_m[[x]]["IPA_optimism"], IPA_corr[[x]]), simplify=FALSE)
  
  
  
  # pobjval <-
  #   as.data.frame(matrix(c(Dxy_val, Cindex_val, AUC_val, Brier_val, IPA_val), 5, 5, byrow = T))
  # 
  
  pobjval <-sapply(horizon %>% as.character(), function(x) as.data.frame(matrix(c(Dxy_val[[x]], Cindex_val[[x]], AUC_val[[x]], 
                                                                                  Brier_val[[x]], IPA_val[[x]]), 5, 5, byrow = T)), simplify=FALSE)
  
  
  
  
  pobjval <-lapply(pobjval, setNames, c("Orig", "Apparent", "Test", "Optimism", "Corrected"))
  pobjval <-sapply(pobjval, function(x) {rownames(x) <- c("Dxy", "Cindex", "AUC", "Brier","IPA"); x}, simplify=FALSE)
  
  pobjval <- sapply(pobjval, function(x) x %>% as.data.frame() %>% rownames_to_column(var="Index") ,simplify = FALSE)
  
  # pobjval <- pobjval %>% as.data.frame() %>% rownames_to_column(var="Index")
  # colnames(pobjval) <-
  #   c("Orig", "Apparent", "Test", "Optimism", "Corrected")
  # row.names(pobjval) <-
  #   c("Dxy", "Cindex", "AUC", "Brier","IPA")
  
  ## HERE IF YOU HAVE BW, DO STABILITY ANALYSIS
  if(p.crit!=1){
    if (direction == "BW"){
      bif <- predictors_selected
      bif_pat <- bif %>% group_by_all %>% dplyr::count()
      bif_pat_sort <- data.frame(bif_pat %>% arrange(desc(n)) %>% 
                                   select(1:ncol(bif_pat)))
      bif_pat_perc <- round((bif_pat_sort$n/nboot) * 100, 0)
      bif_pat_sort <- data.frame(bif_pat_sort, bif_pat_perc)
      
      colnames(bif_pat_sort) <- c(colnames(bif), 
                                  "freq", "bif_pat_perc")
      
      # rownames(bif) <- paste("boot", 1:nboot)
      bif_total <- colSums(bif)
      bif_perc <- round((bif_total/nboot) * 100, 3)
      stabobj <- list(bif = bif, bif_total = bif_total, bif_perc = bif_perc, 
                      model_stab = bif_pat_sort)
      
      
      pobjval <- list(stats_val = pobjval, #intercept_test = res_boot_m[7],
                      res_boot = res_boot, predictors_selected = predictors_selected, strata.predictors=pobj$strata.predictors,
                      model_orig=fm_orig,
                      stability=stabobj)
    }
    
  }else{
    pobjval <- list(stats_val = pobjval, #intercept_test = res_boot_m[7],
                    res_boot = res_boot, predictors_selected = predictors_selected, strata.predictors=pobj$strata.predictors,
                    model_orig=fm_orig)
  }
  
  
  
  ###########################
  
  
  return(pobjval)
}



validation_calibration_rms<- function(data, time, status, features, set_seed, nBoots,timePoint,set.units="Month", BW=FALSE, rule="p", sls=.05){
  
  # data=orig.imp.list.cox.s
  # time=surv_time
  # status=surv_status
  # features=feats.trans.strata %>% as.character %>% str_replace_all("strata","strat")
  # set_seed=2371
  # nBoots=myBoot
  # timePoint=6
  # set.units="Month"
  # BW=TRUE
  # rule="p"
  # sls=.05
  # 
  units(data) <- set.units
  d <- datadist(data)
  options(datadist=d)   
  timePoint <- as.numeric(timePoint)
  
  f <- cph(formula = as.formula(modelFormula(time,status,features) %>% as.character),data=data,x=T, y=T,surv=T,time.inc=timePoint)
  f$units <- set.units
  # ggplot ( Predict (f ) , sepdiscrete = 'vertical', nlevels =4 ,
  # vnames = 'names')
  if(isTRUE(BW)){
    print("running BW with fastbw")
    fit.bw <- rms::fastbw(f, rule=rule,sls=sls)
    # compare with initial features
    init_feats <- str_split(features," \\+ ") %>% magrittr::extract2(1)
    final_feats <- setdiff(init_feats,init_feats[fit.bw$factors.deleted])
    final.feat.fm <- makeComplexFormula(features=final_feats, 
                                        rcsFeats=NULL, 
                                        knots=NULL, 
                                        interactions=NULL,
                                        strataFeats=NULL)
    f.bw <- cph(formula = as.formula(modelFormula(time,status,final.feat.fm) %>% as.character),data=data,x=T, y=T,surv=T,time.inc=timePoint)
    f.bw$units <- set.units
    # VALIDATION
    set.seed (set_seed)
    v <-validate(f , B =nBoots,u=timePoint, bw=BW,rule=rule,sls=sls)
    # CALIBRATION
    set.seed (set_seed)
    cal <-calibrate(f ,method="boot", B =nBoots , u =timePoint ,cmethod='hare',legend = TRUE, digits = 3, subtitles = T,cex.subtitles=0.5, bw=BW,rule=rule,sls=sls)#, maxdim =3
  }else{
    print("running complete model")
    # VALIDATION
    set.seed (set_seed)
    v <-validate(f , B =nBoots,u=timePoint, bw=BW,rule=rule,sls=sls)
    # CALIBRATION
    set.seed (set_seed)
    cal <-calibrate(f ,method="boot", B =nBoots , u =timePoint ,cmethod='hare',legend = TRUE, digits = 3, subtitles = T,cex.subtitles=0.5, bw=BW,rule=rule,sls=sls
                    )#, maxdim =3 # pred=seq(0,1, by=0.02)
  }
  
  # plot ( cal)
  set.seed (set_seed)
  p <-myCalPlot(cal,par.corrected = list(col="black"), lwd=0.2, cex=0.6)
  # p<- plot( cal,par.corrected = list(col="Indian red"))
  p <- recordPlot()
  if(isTRUE(BW)){
    res = list(fit=f,
               fit.bw=f.bw,
               val=v,
               cal=cal,
               cPlot=p)
  }else{
    res = list(fit=f,
               val=v,
               cal=cal,
               cPlot=p)
  }
  
  res
}



