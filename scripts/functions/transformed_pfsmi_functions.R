#' Pooling and Predictor selection function for backward or forward selection of
#' Cox regression models across multiply imputed data.
#'
#' \code{psfmi_coxr} Pooling and backward or forward selection of Cox regression
#' prediction models in multiply imputed data using selection methods D1, D2 and MPR.
#'
#' @param data Data frame with stacked multiple imputed datasets.
#'   The original dataset that contains missing values must be excluded from the
#'   dataset. The imputed datasets must be distinguished by an imputation variable,
#'   specified under impvar, and starting by 1.
#' @param formula A formula object to specify the model as normally used by coxph.
#'   See under "Details" and "Examples" how these can be specified. If a formula
#'   object is used set predictors, cat.predictors, spline.predictors or int.predictors
#'   at the default value of NULL.
#' @param nimp A numerical scalar. Number of imputed datasets. Default is 5.
#' @param impvar A character vector. Name of the variable that distinguishes the
#' imputed datasets.
#' @param time Survival time.
#' @param status The status variable, normally 0=censoring, 1=event.
#' @param predictors Character vector with the names of the predictor variables.
#'   At least one predictor variable has to be defined. Give predictors unique names
#'   and do not use predictor name combinations with numbers as, age2, gnder10, etc.
#' @param cat.predictors A single string or a vector of strings to define the
#' categorical variables. Default is NULL categorical predictors.
#' @param spline.predictors A single string or a vector of strings to define the
#' (restricted cubic) spline variables. Default is NULL spline predictors. See details.
#' @param int.predictors A single string or a vector of strings with the names of the variables that form
#'   an interaction pair, separated by a ":" symbol.
#' @param keep.predictors A single string or a vector of strings including the variables that are forced
#'   in the model during predictor selection. Categorical and interaction variables are allowed.
#' @param nknots A numerical vector that defines the number of knots for each spline predictor separately.
#' @param p.crit A numerical scalar. P-value selection criterion. A value of 1
#'   provides the pooled model without selection.
#' @param method A character vector to indicate the pooling method for p-values to pool the
#'   total model or used during predictor selection. This can be "RR", D1", "D2", or "MPR".
#'   See details for more information. Default is "RR".
#' @param direction The direction of predictor selection, "BW" means backward selection and "FW"
#'   means forward selection.
#'
#' @details The basic pooling procedure to derive pooled coefficients, standard errors, 95
#'  confidence intervals and p-values is Rubin's Rules (RR). However, RR is only possible when
#'  the model included continuous or dichotomous variables. Specific procedures are
#'  available when the model also included categorical (> 2 categories) or restricted cubic spline
#'  variables. These pooling methods are: "D1" is pooling of the total covariance matrix,
#'  "D2" is pooling of Chi-square values and "MPR" is pooling of median p-values (MPR rule).
#'  Spline regression coefficients are defined by using the rcs function for restricted cubic
#'  splines of the rms package. A minimum number of 3 knots as defined under knots is required.
#'
#'  A typical formula object has the form \code{Surv(time, status) ~ terms}. Categorical variables has to
#'  be defined as \code{Surv(time, status) ~ factor(variable)}, restricted cubic spline variables as
#'  \code{Surv(time, status) ~ rcs(variable, 3)}. Interaction terms can be defined as
#'  \code{Surv(time, status) ~ variable1*variable2} or \code{Surv(time, status) ~ variable1 + variable2 + 
#'  variable1:variable2}. All variables in the terms part have to be separated by a "+". If a formula
#'  object is used set predictors, cat.predictors, spline.predictors or int.predictors
#'  at the default value of NULL.
#'
#'
#' @return An object of class \code{pmods} (multiply imputed models) from
#'  which the following objects can be extracted: 
#'  \itemize{
#'  \item  \code{data} imputed datasets
#'  \item  \code{RR_model} pooled model at each selection step
#'  \item  \code{RR_model_final} final selected pooled model
#'  \item  \code{multiparm} pooled p-values at each step according to pooling method
#'  \item  \code{multiparm_final} pooled p-values at final step according to pooling method
#'  \item  \code{multiparm_out} (only when direction = "FW") pooled p-values of removed predictors 
#'  \item  \code{formula_step} formula object at each step
#'  \item  \code{formula_final} formula object at final step
#'  \item  \code{formula_initial} formula object at final step
#'  \item  \code{predictors_in} predictors included at each selection step
#'  \item  \code{predictors_out} predictors excluded at each step
#'  \item  \code{impvar} name of variable used to distinguish imputed datasets
#'  \item  \code{nimp} number of imputed datasets
#'  \item  \code{status} name of the status variable
#'  \item  \code{time} name of the time variable
#'  \item  \code{method} selection method
#'  \item  \code{p.crit} p-value selection criterium
#'  \item  \code{call} function call
#'  \item  \code{model_type} type of regression model used
#'  \item  \code{direction} direction of predictor selection
#'  \item  \code{predictors_final} names of predictors in final selection step
#'  \item  \code{predictors_initial} names of predictors in start model
#'  \item  \code{keep.predictors} names of predictors that were forced in the model    
#' }
#'
#' @references Eekhout I, van de Wiel MA, Heymans MW. Methods for significance testing of categorical
#'   covariates in logistic regression models after multiple imputation: power and applicability
#'   analysis. BMC Med Res Methodol. 2017;17(1):129.
#' @references Enders CK (2010). Applied missing data analysis. New York: The Guilford Press.
#' @references van de Wiel MA, Berkhof J, van Wieringen WN. Testing the prediction error difference between
#'   2 predictors. Biostatistics. 2009;10:550-60.
#' @references Marshall A, Altman DG, Holder RL, Royston P. Combining estimates of interest in prognostic
#'   modelling studies after multiple imputation: current practice and guidelines. BMC Med Res Methodol.
#'   2009;9:57.
#' @references Van Buuren S. (2018). Flexible Imputation of Missing Data. 2nd Edition. Chapman & Hall/CRC
#'   Interdisciplinary Statistics. Boca Raton.
#' @references EW. Steyerberg (2019). Clinical Prediction MOdels. A Practical Approach
#'  to Development, Validation, and Updating (2nd edition). Springer Nature Switzerland AG.
#'
#' @references http://missingdatasolutions.rbind.io/
#'
#' @section Vignettes:
#'   https://mwheymans.github.io/psfmi/articles/psfmi_CoxModels.html
#'   
#' @author Martijn Heymans, 2020
#' 
#' @examples
#'  pool_coxr <- psfmi_coxr(formula = Surv(Time, Status) ~ Pain + Tampascale +
#'                        Radiation + Radiation*Pain + Age + Duration + Previous,
#'                      data=lbpmicox, p.crit = 0.05, direction="BW", nimp=5, impvar="Impnr",
#'                      keep.predictors = "Radiation*Pain", method="D1")
#'                      
#'  pool_coxr$RR_model_final
#'  
#' @export
psfmi_coxr.Fixed <- function(data,
                             formula = NULL,
                             nimp=5,
                             impvar=NULL,
                             status = NULL,
                             time = NULL,
                             predictors=NULL,
                             cat.predictors=NULL,
                             spline.predictors=NULL,
                             int.predictors=NULL,
                             keep.predictors=NULL,
                             nknots=NULL,
                             strata.predictors=NULL,
                             p.crit=1,
                             method="RR",
                             direction=NULL)
{
  
  #########################
  # data=mydata
  # nimp=2
  # impvar="Impnr" 
  # formula = modelFormula(surv_time, surv_status, feats1)
  # p.crit=1
  # method="D1"
  # direction = "BW"
  # status = NULL
  # time = NULL
  # predictors=NULL
  # cat.predictors=NULL
  # spline.predictors=NULL
  # int.predictors=NULL
  # keep.predictors=NULL
  # nknots=NULL
  
  # data=mydata
  # nimp=2
  # impvar="Impnr"
  # formula = modelFormula(surv_time, surv_status, feats1)
  # p.crit=1
  # method="D1"
  # direction = "BW"
  # status = NULL
  # time = NULL
  # predictors=NULL
  # cat.predictors=NULL
  # spline.predictors=NULL
  # int.predictors=NULL
  # keep.predictors=NULL
  # nknots=NULL
  
  # 
  # data=orig.imp.final.cox
  # nimp=10
  # impvar="Impnr"
  # formula = modelFormula(surv_time, surv_status, feats.trans.strata )
  # p.crit=1
  # method="D1"
  # direction = "BW"
  # strata.predictors=strataVars
  # status = NULL
  # time = NULL
  # predictors=NULL
  # cat.predictors=NULL
  # spline.predictors=NULL
  # int.predictors=NULL
  # keep.predictors=NULL
  # nknots=NULL
  # 
  #########################
  
  
  call <- match.call()
  
  
  
  
  if(method=="D3")
    stop("\n", "Method D3 not available for survival data")
  if(is_empty(formula)) {
    if(!all(data[status]==1 | data[status]==0))
      stop("Status should be a 0 - 1 variable")
    
    keep_temp <- keep.predictors
    
    P <-
      predictors
    cat.P <-
      cat.predictors
    int.P <-
      gsub(":", "*", int.predictors)
    s.P <-
      spline.predictors
    strata.P <-  strata.predictors
  } else{
    form <-
      terms(formula)
    form_vars <-
      attr(form, "term.labels")
    if(is_empty(form_vars))
      stop("\n", "No predictors defined, model is empty")
    time <-
      as.character(attr(form, "variables")[[2]][[2]])
    status <-
      as.character(attr(form, "variables")[[2]][[3]])
    int.P <-
      form_vars[grepl(paste(c("[*]", ":"), collapse = "|"), form_vars)]
    int.P_temp <-
      unique(unlist(str_split(int.P, paste(c("[*]", ":"), collapse = "|"))))
    form_vars <-
      form_vars[!grepl(paste(c("[*]", ":"), collapse = "|"), form_vars)]
    form_vars <-
      unique(c(form_vars, int.P_temp))
    cat.P <-
      form_vars[grepl("factor", form_vars)]
    form_vars <-
      form_vars[!grepl("factor", form_vars)]
    s.P <-
      form_vars[grepl("rcs", form_vars)]
    nknots <-
      c(readr::parse_number(s.P))
    form_vars <-
      form_vars[!grepl("rcs", form_vars)]
    int.P <-
      gsub(":", "*", myclean_P(int.P))
    cat.P <- clean_P(cat.P)
    s.P <- clean_P(s.P)
    # Get variables stratified
    strata.P <-
      form_vars[grepl("strata", form_vars)]
    form_vars <-
      form_vars[!grepl("strata", form_vars)]
    strata.P <-  myclean_P(strata.P)
    ##
    P <- form_vars
  }
  
  keep.P <-
    gsub(":", "*", keep.predictors)
  keep.P <-
    sapply(as.list(keep.P), clean_P)
  
  P.check <-
    c(P, cat.P, s.P,strata.P)
  
  # Check data input
  if(p.crit!=1){
    if(is_empty(direction))
      stop("Specify FW or BW for forward or backward predictor selection")
  }
  if(!(is.data.frame(data)))
    stop("Data should be a data frame")
  data <- data.frame(as_tibble(data))
  # data <- mutate_if(data, is.factor, ~ as.numeric(as.character(.x)))
  if(!all(data[status]==1 | data[status]==0))
    stop("Status should be a 0 - 1 variable")
  if ((nvar <- ncol(data)) < 2)
    stop("Data should contain at least two columns")
  if(is_empty(impvar))
    stop("Imputation variable is not defined")
  if(is_empty(method)) method="RR"
  if(all(!is_empty(cat.P) | !is_empty(s.P)) & method=="RR")
    stop("Categorical or spline variables in model, define selection method: D1, D2 or MPR")
  if (sort(unique(data[, impvar]))[1] == 0)
    stop("Original dataset should not be included")
  if(is_empty(nimp))
    stop("Number of imputed datasets is not defined, use nimp!")
  if (nimp < 2) {
    stop("\n", "Number of imputed datasets must be > 1", "\n\n")
  }
  if (p.crit > 1)
    stop("\n", "P-value criterium > 1", "\n")
  if (any(nknots<3))
    stop("\n", "Number of knots must be > 2", "\n")
  if (length(nknots) != length(s.P))
    stop("\n", "Number of knots not specified for every spline variable", "\n")
  if (!is_empty(cat.P)) {
    if(any(cat.P%in%P)){
      cat.P.double <- cat.P[cat.P%in%P]
      stop("\n", "Categorical variable(s) -", cat.P.double,
           "- also defined as Predictor", "\n\n")
    }
  }
  if (!is_empty(s.P)){
    if(any(s.P%in%P)){
      s.P.double <- s.P[s.P%in%P]
      stop("\n", "Do not include Spline variable(s) -", s.P.double,
           "- in predictors", "\n\n")
    }
  }
  
  if (!is_empty(strata.P)){
    if(any(strata.P%in%P)){
      strata.P.double <- strata.P[strata.P%in%P]
      stop("\n", "Do not include strata variable(s) -", strata.P.double,
           "- in predictors", "\n\n")
    }
  }
  
  if(any(duplicated(P))){
    stop("\n", "Predictor(s) - ", c(P[duplicated(P)]),
         " - defined more than once", "\n\n")
  }
  # Check if al variables are available in dataset
  if(any(!P.check %in% names(data))) {
    P.mis <- P.check[!P.check %in% names(data)]
    stop("\n", "Predictor(s) - ", P.mis,
         "- not available in dataset", "\n\n")
  }
  if(!is_empty(int.P)) {
    int.P.check <- lapply(int.P[grep("[*]", int.P)],
                          function(x) { unlist(strsplit(x, split="[*]")) })
    int.P.check <- unique(unlist(int.P.check))
    if(any(!int.P.check %in% P.check))
      stop("\n", "Not all interaction terms defined as
        Predictor or Categorical Predictor", "\n\n")
  }
  # First predictors, second categorical
  # predictors and last interactions
  P <-
    c(P, cat.P, s.P, int.P,strata.P)
  if (is_empty(P))
    stop("\n", "No predictors to select, model is empty", "\n\n")
  
  if (!is_empty(keep.P)) {
    for(i in 1:length(keep.P)){
      if(grepl("[*]", keep.P[i])) {
        keep.P.spl <- unlist(strsplit(keep.P[i], split="[*]"))
        if(length(P[Reduce("&", lapply(keep.P.spl, grepl, P))])==0)
          stop("Interaction term in keep.predictors not defined
            as int.predictors, incorrect")
        keep.P[i] <- P[Reduce("&", lapply(keep.P.spl, grepl, P))]
      }
    }
  }
  
  if (!is_empty(cat.P)) {
    if(length(cat.P)==1){
      P <-
        gsub(cat.P,
             replacement=paste0("factor(", cat.P, ")"), P)
      if(!is_empty(keep.P)){
        keep.P <-
          gsub(cat.P,
               replacement=paste0("factor(", cat.P, ")"), keep.P)
      }
    } else {
      for(i in 1:length(cat.P)) {
        P <-
          gsub(cat.P[i],
               replacement=paste0("factor(", cat.P[i], ")"), P)
        if(!is_empty(keep.P)){
          keep.P <-
            gsub(cat.P[i],
                 replacement=paste0("factor(", cat.P[i], ")"), keep.P)
        }
      }
    }
  }
  if (!is_empty(s.P)) {
    if(length(s.P)==1){
      P <-
        gsub(s.P,
             replacement=paste0("rcs(", s.P, ",", nknots, ")"), P)
      if(!is_empty(keep.P)){
        keep.P <-
          gsub(s.P,
               replacement=paste0("rcs(", s.P, ",", nknots, ")"), keep.P)
      }
    } else {
      for(i in 1:length(s.P)) {
        P <- gsub(s.P[i],
                  replacement=paste0("rcs(", s.P[i], ",", nknots[i], ")"), P)
        if(!is_empty(keep.P)){
          keep.P <-
            gsub(s.P[i],
                 replacement=paste0("rcs(", s.P[i], ",", nknots[i], ")"), keep.P)
        }
      }
    }
  }
  # Same as for splines but for strata
  if (!is_empty(strata.P)) {
    if(length(strata.P)==1){
      P <-stringi::stri_replace_all_regex(P,
                                          pattern=strata.P,
                                          replacement=paste0("strata(", strata.P, ")"),
                                          vectorize=FALSE)
      
      if(!is_empty(keep.P)){
        keep.P <-stringi::stri_replace_all_regex(keep.P,
                                                 pattern=strata.P,
                                                 replacement=paste0("strata(", strata.P, ")"),
                                                 vectorize=FALSE)
          # gsub(strata.P,
          #      replacement=paste0("strata(", strata.P, ")"), keep.P)
        
        
        
      }
    } else {
      for(i in 1:length(strata.P)) {
        P <- gsub(strata.P[i],
                  replacement=paste0("strata(", strata.P[i], ")"), P)
        if(!is_empty(keep.P)){
          keep.P <-
            gsub(strata.P[i],
                 replacement=paste0("strata(", strata.P[i], ")"), keep.P)
        }
      }
    }
  }
  #####
  levels.cat.P <- lapply(cat.P, function(x) {
    nr.levels.cat.P <- length(table(data[data[impvar] == 1, ][, x]))
    if (nr.levels.cat.P < 3) {
      stop("\n", "Categorical variable(s) only 2 levels,
        do not define as categorical", "\n\n")
    }
  })
  
  if(any(!keep.P %in% P))
    stop("\n", "Variables to keep not defined as Predictor", "\n\n")
  
  if(p.crit==1){
    # pobjpool <-
    #   my.psfmi_coxr_bw(data = data, nimp=nimp, impvar = impvar, status = status, time = time,
    #                 P = P, p.crit = p.crit, method = method, keep.P = keep.P)
    
    pobjpool <-my.psfmi_coxr_bw(data = data, nimp=nimp, impvar = impvar, status = status, time = time,
                                P = P, p.crit = p.crit, method = method, keep.P = keep.P,strata.P=strata.P)
    class(pobjpool) <-
      "pmods"
    return(pobjpool)
  }
  if(direction=="FW"){
    
    pobjfw <-
      psfmi_coxr_fw(data = data, nimp = nimp, impvar = impvar, status = status, time = time, p.crit = p.crit,
                    P = P, keep.P = keep.P, method = method)
    class(pobjfw) <-
      "pmods"
    return(pobjfw)
  }
  if(direction=="BW"){
    # pobjbw <-
    #   my.psfmi_coxr_bw(data = data, nimp=nimp, impvar = impvar, status = status, time = time,
    #                 P = P, p.crit = p.crit, method = method, keep.P = keep.P)
    
    pobjbw <-
      my.psfmi_coxr_bw(data = data, nimp=nimp, impvar = impvar, status = status, time = time,
                       P = P, p.crit = p.crit, method = method, keep.P = keep.P,strata.P=strata.P)
    class(pobjbw) <-
      "pmods"
    return(pobjbw)
  }
}

my.psfmi_coxr_bw <- function(data, nimp, impvar, status, time, P, p.crit, method, keep.P,strata.P)
{
  
  
  # # # # # #####################
  # data = data
  # nimp=nimp
  # impvar = impvar
  # status = status
  # time = time
  # P = P
  # p.crit = p.crit
  # method = method
  # keep.P = keep.P
  # strata.P=strata.P
  # 
  # 
  # #####################
  # message("CORRECT FUNCTION")
  
  call <- match.call()
  
  impvar.un <- as.name(impvar)
  ndata <- nrow(data %>% dplyr::filter(!!impvar.un==1))
  RR.model <- P_rm_step <- fm_step <- imp.dt <- multiparm <- pool_lp <- pool_covCenter<- list()
  
  # If there are strata variables, remove them from the predictors toloop through -ISSUE WHEN STRATA IS PART OF INTERACTION
  # if(!is_empty(strata.P)){
  #   P.strata.toRemove <-P[grepl("strata", P)]
  #   P <- setdiff(P,P.strata.toRemove)
  # }
  
  ## Strat in interactions
  P.inter=P[grepl("\\*", P)]
  
  if(!is_empty(strata.P)){
    P.strata.toRemove <-P[grepl("strata", P)]
    P <- c(setdiff(P,P.strata.toRemove),P.inter) %>% unique()
  }
  
  
  P_orig <-
    P
  keep.P <-
    sapply(as.list(keep.P), clean_P)
  
  P_orig_temp <-
    myclean_P(P)
  
  if(!is_empty(keep.P))
    if(length(P_orig)==1)
      if(P_orig_temp == keep.P)
        stop("\n", "No need to define keep.predictors. Exclude keep.predictors and set p.crit = 1","\n")
  
  # Loop k, to pool models in multiply imputed datasets
  for (k in 1:(length(P)+1)) {
    # # # ##
    # k=1
    # # # # print(paste0("STEP:",k))
    # 
    # #
    # set regression formula fm
    Y <-
      c(paste0("Surv(", time, ",", status, ")~"))
    
    if(!is_empty(strata.P)){
      fm <-
        as.formula(paste(Y, paste(c(P, paste0("strata(",strata.P,")")), collapse = "+")))
      
    }else{
      fm <-
        as.formula(paste(Y, paste(P, collapse = "+")))
      
    }
    
    
    
    
    # Extract df of freedom for MPR
    if(method=="MPR" | method=="RR"){
      
      chi.LR <-
        data.frame(matrix(0, length(attr(terms(fm), "term.labels")), nimp))
      chi.p <-
        data.frame(matrix(0, length(attr(terms(fm), "term.labels")), nimp))
      
      fit <- list()
      for (i in 1:nimp) {
        imp.dt[[i]] <- data[data[impvar] == i, ]
        fit[[i]] <- coxph(fm, data = imp.dt[[i]])
        if(length(attr(terms(fm), "term.labels")) == 1){
          chi.LR[, i] <- car::Anova(fit[[i]])[-1, 2]
          chi.p[, i] <- car::Anova(fit[[i]])[-1, 4]
        } else {
          chi.LR[, i] <- car::Anova(fit[[i]])[, 1]
          chi.p[, i] <- car::Anova(fit[[i]])[, 3]
        }
      }
      
      # Rubin's Rules
      out.res <-
        suppressWarnings(summary(pool(fit)))
      HR <-
        exp(out.res$estimate)
      lower.EXP <-
        exp(out.res$estimate - (qt(0.975, out.res$df)*out.res$std.error))
      upper.EXP <-
        exp(out.res$estimate + (qt(0.975, out.res$df)*out.res$std.error))
      model.res <-
        data.frame(cbind(out.res, HR, lower.EXP, upper.EXP))
      RR.model[[k]] <-
        model.res
      names(RR.model)[[k]] <-
        paste("Step", k)#
    }
    
    # D1 and D2 pooling methods
    if(method=="D1" | method == "D2") {
      
      pool.p.val <-
        matrix(0, length(P), 2)
      P_test <-
        myclean_P(P)
      
      # ####################################
      # # To not repeat again and again with predictors looping, extract the full model info here
      fit1 <- list()
      
      
      if(!is_empty(strata.P)){
        
        form1 <-
          as.formula(paste(Y, paste(c(P, paste0("strata(",strata.P,")")), collapse = "+")))
        
      }else{
        form1 <-
          as.formula(paste(Y, paste(P, collapse = "+")))
        
      }
      
      
      
      for (i in 1:nimp) {
        imp.dt[[i]] <-
          data[data[impvar] == i, ]
        fit1[[i]] <-
          coxph(form1, data = imp.dt[[i]])
        # fit0[[i]] <-
        #   coxph(form0, data = imp.dt[[i]])
      }
      out.res1 <-
        suppressWarnings(summary(pool(fit1)))
      HR <-
        exp(out.res1$estimate)
      lower.EXP <-
        exp(out.res1$estimate - (qt(0.975, out.res1$df)*out.res1$std.error))
      upper.EXP <-
        exp(out.res1$estimate + (qt(0.975, out.res1$df)*out.res1$std.error))
      model.res1 <-
        data.frame(cbind(out.res1, HR, lower.EXP, upper.EXP))
      RR.model[[k]] <-
        model.res1
      names(RR.model)[[k]] <-
        paste("Step", k)
      
      # Get pooled linear predictors with Rubin's Rules
      pool_lp.list <- sapply(fit1, function(x) x$linear.predictors,simplify = FALSE)
      pool_lp.df <-ldply(pool_lp.list)
      pool_lp[[k]]<-apply(pool_lp.df,2, function(x) pool_RR(est=x,
                                                            se=stderror(x),
                                                            n=ndata,
                                                            k=1) %>% pluck("Estimate"))
      names(pool_lp)[[k]] <-
        paste("Step", k)
      # Get pooled covariate means
      pool_covCenter.list <- sapply(fit1, function(x) coxCenter(x),simplify = FALSE)
      pool_covCenter.df <-ldply(pool_covCenter.list)
      pool_covCenter[[k]] <-apply(pool_covCenter.df,2, function(x) my.pool_RR_covMeans(est=x,
                                                                                       se=stderror(x),
                                                                                       n=ndata,
                                                                                       k=1)
      )
      names(pool_covCenter)[[k]] <-
        paste("Step", k)
      # ####################################
      
      for (j in 1:length(P)) {
        # ## Here is backward selection, removes each predictor sequentially from full model
        # # ###
        # j=1
        # # # ###
        cov.nam0 <-
          P[-j]
        if (length(P) == 1) {
          cov.nam0 <-
            "1"
        }
        Y <-
          c(paste0("Surv(", time, ",", status, ")~"))
        
        
        if(!is_empty(strata.P)){
          form0 <-
            as.formula(paste(Y, paste(c(cov.nam0, paste0("strata(",strata.P,")")), collapse = "+")))
          
        }else{
          form0 <-
            as.formula(paste(Y, paste(cov.nam0, collapse = "+")))
          
        }
        
        
        
        if(any(grepl(P_test[j], P_test[-j]))){
          cov.nam0 <-
            P[-grep(P_test[j], P_test)]
          
          
          if(!is_empty(strata.P)){
            
            form0 <-
              as.formula(paste(Y, paste(c(cov.nam0, paste0("strata(",strata.P,")")), collapse = "+")))
          }else{
            form0 <-
              as.formula(paste(Y, paste(c(cov.nam0), collapse = "+")))
            
          }
          
          
        }
        
        if(method =="D1" | method =="D2")
          # fit1 <- fit0 <- imp.dt <- list()
          fit0 <- imp.dt <- list()
        for (i in 1:nimp) {
          imp.dt[[i]] <-
            data[data[impvar] == i, ]
          # fit1[[i]] <-
          #   coxph(form1, data = imp.dt[[i]])
          fit0[[i]] <-
            coxph(form0, data = imp.dt[[i]])
        }
        
        test_P <-
          mitml::testModels(fit1, fit0, method = method)
        pvalue <-
          test_P$test[4]
        fstat <-
          test_P$test[1]
        pool.p.val[j, ] <-
          c(pvalue, fstat)
        
        
      }
      p.pool <-
        data.frame(pool.p.val)
      row.names(p.pool) <-
        P
      names(p.pool) <-
        c(paste("p-values", method), "F-statistic")
    }
    
    # MPR Pooling
    if(method=="MPR") {
      p.pool <-
        data.frame(apply(chi.p, 1 , median))
      rownames(p.pool) <-
        P
      names(p.pool) <-
        c("p-value MPR")
    }
    
    # RR Pooling
    if(method=="RR") {
      p.pool <-
        data.frame(RR.model[[k]][, 6],
                   row.names=P)
      names(p.pool) <-
        c("p-value RR")
    }
    
    # Extract regression formula's
    fm_step[[k]] <-
      formula(fm)
    names(fm_step)[[k]] <-
      paste("Step", k)
    # Extract multiparameter pooling
    multiparm[[k]] <-
      p.pool
    names(multiparm)[[k]] <-
      paste("Step", k)
    
    # Clean variable names for selection
    P_temp <-
      clean_P(row.names(p.pool))
    p.pool_temp <-
      p.pool
    row.names(p.pool_temp) <-
      P_temp
    # detect interaction terms and exclude variables
    # that are part of interaction
    if(any(grepl("[*]", P))){
      P_int_id <-
        P_temp[grep("[*]", P_temp)]
      # Detect main effects of interaction terms
      P_main <-
        unique(unlist(str_split(P_int_id, "[*]")))
      # Exclude variables to keep and main effect from selection
      remove_id <-
        !row.names(p.pool_temp) %in% unique(c(P_main))
      p.pool <-
        p.pool[remove_id, , FALSE]
      p.pool_temp <-
        p.pool_temp[remove_id, , FALSE]
    }
    if(!is_empty(keep.P)){
      remove_P_keep <-
        !row.names(p.pool_temp) %in% keep.P
      p.pool <-
        p.pool[remove_P_keep, , FALSE]
    }
    
    if(p.crit==1){
      # message("\n", "No backward selection chosen, break out of loop", "\n")
      break()
    }
    
    if(nrow(p.pool)==0)
      break()
    
    # print("Continuing to select variables")
    # Select variables
    P_temp <-
      row.names(p.pool)
    P_excl <-
      which(p.pool[, 1] == max(p.pool[, 1]))
    if(length(P_excl) > 1) {
      P_excl <-
        P_excl[1]
    }
    
    if (p.pool[, 1][P_excl] < p.crit) {
      message("\n", "Selection correctly terminated, ",
              "\n", "No more variables removed from the model", "\n")
      (break)()
    }
    
    P_out <-
      P_temp[P_excl]
    
    if(p.pool[, 1][P_excl] > p.crit) {
      message("Removed at Step ", k,
              " is - ", P_out)
      
    }
    
    P_rm_step[[k]] <-
      P_out
    # Variable excluded on each step
    P <-
      P[!P %in% P_out]
    
    if(is_empty(P)){
      print("Empty variable vector")
      fm_step[[k+1]] <-
        as.formula(paste(Y, 1))
      names(fm_step)[[k+1]] <-
        paste("Step", k+1)
      multiparm[[k+1]] <-
        0
      names(multiparm)[[k+1]] <-
        paste("Step", k+1)
      fit <- list()
      for (i in 1:nimp) {
        print("Fit models for each imputed dataset")
        imp.dt[[i]] <-
          data[data[impvar] == i, ]
        fit[[i]] <-
          coxph(fm_step[[k+1]], data = imp.dt[[i]])
      }
      
      if(is_empty(P)){
        RR.model[[k+1]] <-
          fit[[1]]
        names(RR.model)[[k+1]] <-
          paste("Step", k+1)
        break()
      }
      
      print("Create output results of pooled model")
      # Rubin's Rules
      out.res <-
        suppressWarnings(summary(pool(fit)))
      HR <-
        exp(out.res$estimate)
      lower.EXP <-
        exp(out.res$estimate - (qt(0.975, out.res$df)*out.res$std.error))
      upper.EXP <-
        exp(out.res$estimate + (qt(0.975, out.res$df)*out.res$std.error))
      model.res <-
        data.frame(cbind(out.res, HR, lower.EXP, upper.EXP))
      RR.model[[k+1]] <-
        model.res
      names(RR.model)[[k+1]] <-
        paste("Step", k+1)
      
      # Get pooled linear predictors
     
      pool_lp.list <- sapply(fit, function(x) x$linear.predictors,simplify = FALSE)
      pool_lp.df <-ldply(pool_lp.list)
      pool_lp[[k+1]]<-apply(pool_lp.df,2, function(x) pool_RR(est=x,
                                                              se=stderror(x),
                                                              n=ndata, 
                                                              k=1) %>% pluck("Estimate"))
      names(pool_lp)[[k+1]] <-
        paste("Step", k+1)
      # Get pooled covariate means
      pool_covCenter.list <- sapply(fit, function(x) coxCenter(x),simplify = FALSE)
      pool_covCenter.df <-ldply(pool_covCenter.list)
      
      pool_covCenter[[k+1]] <-apply(pool_covCenter.df,2, function(x) my.pool_RR_covMeans(est=x,
                                                                                         se=stderror(x),
                                                                                         n=ndata, 
                                                                                         k=1)
      )
      names(pool_covCenter)[[k+1]] <-
        paste("Step", k+1)
      
      break()
    }
  } # End k loop

  P_rm_step_final <-
    P_rm_step
  if(is_empty(P)) {
    P_rm_step <-
      lapply(P_rm_step, clean_P)
  } else {
    P_rm_step <-
      lapply(P_rm_step[-k], clean_P)
  }
  
  # Extract selected models
  outOrder_step <-
    P_orig_temp
  if(!is_empty(P_rm_step)){
    P_remove <-
      data.frame(do.call("rbind",
                         lapply(P_rm_step, function(x) {
                           outOrder_step %in% x })))
    names(P_remove) <-
      P_orig
    row.names(P_remove) <-
      paste("Step", 1:length(P_rm_step))
    if(nrow(P_remove)!=1) {
      P_remove <-
        apply(P_remove, 2, function(x) ifelse(x, 1, 0))
      P_remove_final <-
        colSums(P_remove)
      P_remove <-
        rbind(P_remove, P_remove_final)
      row.names(P_remove)[nrow(P_remove)] <-
        "Removed"
    } else {
      P_remove <-
        matrix(apply(P_remove, 2, function(x) ifelse(x, 1, 0)), 1, length(P_orig))
      dimnames(P_remove) <-
        list("Removed", P_orig)
    }
  } else {
    P_remove <-
      matrix(rep(0, length(P_orig)), 1, length(P_orig))
    dimnames(P_remove) <-
      list("Removed", P_orig)
  }
  
  P_included <-
    as_tibble(names(P_remove[nrow(P_remove), ][P_remove[nrow(P_remove), ] == 0] ))
  predictors_final <-
    names(P_remove[nrow(P_remove), ][P_remove[nrow(P_remove), ] == 0])
  if(length(P_orig)==1 & !is_empty(P))
    P_included <- predictors_final <- P
  
  
  # If there are predictors and NO VARIABLE SELECTION
  if(is_empty(P) | p.crit==1){
    RR_model_step <-
      RR.model
    multiparm_step <-
      multiparm
    fm_step_total <-
      fm_step
    # My extras
    pool_lp_step <-
      pool_lp
    pool_covCenter_step <-
      pool_covCenter
    ####
    if(p.crit!=1) names(RR_model_step) <- names(multiparm_step) <- names(fm_step_total) <-
      paste("Step", 1:(k+1), "- removal -", c(unlist(P_rm_step_final), "ended"))
    if(p.crit==1)  names(RR_model_step) <- names(multiparm_step) <- names(fm_step_total) <-
      paste("Step", 1, "- no variables removed -")
    RR_model_final <-
      RR_model_step[k+1]
    multiparm_final <-
      multiparm_step[k+1]
    fm_step_final <-
      fm_step_total[k+1]
    if(p.crit==1) {
      Y_initial <-
        c(paste0("Surv(", time, ",", status, ")~"))
      
      
      if(!is_empty(strata.P)){
        
        formula_initial <-
          as.formula(paste(Y_initial, paste(c(P_orig, paste0("strata(",strata.P,")")), collapse = "+")))
      }else{
        formula_initial <-
          as.formula(paste(Y_initial, paste(P_orig, collapse = "+")))
        
      }
      
      
      fm_step_final <- 
        formula_initial
      RR_model_final <-
        RR_model_step[k]
      multiparm_final <-
        multiparm_step[k]
      # My extras
      pool_lp_final <-
        pool_lp[k]
      pool_covCenter_final <-
        pool_covCenter[k]
      ####
    }
  }
  
  # If there are predictors and VARIABLE SELECTION
  if(!is_empty(P) & p.crit !=1){
    # Steps
    RR_model_step <-
      RR.model
    multiparm_step <-
      multiparm
    fm_step_total <-
      fm_step
    
    # My extras
    pool_lp_step <-
      pool_lp
    pool_covCenter_step <-
      pool_covCenter
    ####
    
    names(RR_model_step) <-
      names(multiparm_step) <- names(fm_step_total) <- names(pool_lp_step) <- names(pool_covCenter_step) <-
      paste("Step", 1:k, "- removal -", c(unlist(P_rm_step_final), "ended"))
    # Final
    RR_model_final <-
      RR.model[k]
    multiparm_final <-
      multiparm[k]
    fm_step_final <-
      fm_step[k]
    # My extras
    pool_lp_final <-
      pool_lp[k]
    pool_covCenter_final <-
      pool_covCenter[k]
    ####
  }
  
  Y_initial <-
    c(paste0("Surv(", time, ",", status, ")~"))
  
  if(!is_empty(strata.P)){
    formula_initial <-
      as.formula(paste(Y_initial,  paste(c(P_orig, paste0("strata(",strata.P,")")), collapse = "+")))
    
  }else{
    formula_initial <-
      as.formula(paste(Y_initial, paste(P_orig, collapse = "+")))
    
  }
  
  ## final fix if we have strata
  if(!is_empty(strata.P)){
    predictors_final = c(predictors_final,paste0("strata(",strata.P,")"))
    predictors_initial = c(P_orig,paste0("strata(",strata.P,")"))
    
  }else{
    predictors_initial <- P_orig
  }
  
  
  bw <-
    list(data = data, 
         RR_model = RR_model_step, RR_model_final = RR_model_final,
         multiparm = multiparm_step, multiparm_final = multiparm_final, 
         pool_lp = pool_lp_step, pool_lp_final = pool_lp_final,
         pool_covCenter = pool_covCenter_step,pool_covCenter_final=pool_covCenter_final,
         formula_step = fm_step_total, formula_final = fm_step_final,
         formula_initial = formula_initial, 
         predictors_in = P_included, predictors_out = P_remove, 
         impvar = impvar, nimp = nimp, status = status, time = time,
         method = method, p.crit = p.crit, call = call, 
         model_type = "survival", direction = "BW",
         predictors_final = predictors_final, predictors_initial = predictors_initial, 
         keep.predictors = keep.P, strata.predictors=strata.P)
  return(bw)
}




#' Internal validation and performance of logistic prediction models across Multiply Imputed datasets
#'
#' \code{psfmi_validate} Evaluate Performance of logistic regression models selected with
#'  the \code{psfmi_lr} function of the \code{psfmi} package by using cross-validation
#'  or bootstrapping.
#'
#' @param pobj An object of class \code{pmods} (pooled models), produced by a previous
#'  call to \code{psfmi_lr}.
#' @param val_method Method for internal validation. MI_boot for first Multiple Imputation and
#'  than bootstrapping in each imputed dataset and boot_MI for first bootstrapping and than
#'  multiple imputation in each bootstrap sample, and cv_MI, cv_MI_RR and
#'  MI_cv_naive for the combinations of cross-validation and multiple imputation.
#'  To use cv_MI, cv_MI_RR and boot_MI, data_orig has to be specified. See details for more information.
#' @param data_orig dataframe of original dataset that contains missing data for methods
#'  cv_MI, cv_MI_RR and boot_MI.
#' @param int_val If TRUE internal validation is conducted using bootstrapping or cross-validation.
#'  Default is TRUE. If FALSE only apparent performance measures are calculated.
#' @param nboot The number of bootstrap resamples, default is 10. Used for methods boot_MI and MI_boot.
#' @param folds The number of folds, default is 3. Used for methods cv_MI, cv_MI_RR and MI_cv_naive.
#' @param nimp_cv Numerical scalar. Number of (multiple) imputation runs for method cv_MI.
#' @param nimp_mice Numerical scalar. Number of imputed datasets for method cv_MI_RR and boot_MI.
#'  When not defined, the number of multiply imputed datasets is used of the
#'  previous call to the function \code{psfmi_lr}.
#' @param p.crit A numerical scalar. P-value selection criterium used for backward or forward
#'  selection during validation. When set at 1, pooling and internal validation is done without
#'  backward selection.
#' @param BW Only used for methods cv_MI, cv_MI_RR and MI_cv_naive. If TRUE backward selection is
#'  conducted within cross-validation. Default is FALSE.
#' @param direction Can be used together with val_methods boot_MI and MI_boot. The direction of
#'  predictor selection, "BW" is for backward selection and "FW" for forward selection.
#' @param cv_naive_appt Can be used in combination with val_method MI_cv_naive. Default is TRUE for
#'  showing the cross-validation apparent (train) and test results. Set to FALSE to only give test results.
#' @param cal.plot If TRUE a calibration plot is generated. Default is FALSE. Can be used in combination
#'  with int_val = FALSE.
#' @param plot.method If "mean" one calibration plot is generated, first taking the 
#'   mean of the linear predictor across the multiply imputed datasets (default), if 
#'   "individual" the calibration plot of each imputed dataset is plotted, 
#'   if "overlay" calibration plots from each imputed datasets are plotted in one figure.       
#' @param groups_cal A numerical scalar. Number of groups used on the calibration plot and. 
#'  for the Hosmer and Lemeshow test. Default is 10. If the range of predicted probabilities. 
#'  is low, less than 10 groups can be chosen, but not < 3. 
#' @param miceImp Wrapper function around the \code{mice} function.
#' @param ...  Arguments as predictorMatrix, seed, maxit, etc that can be adjusted for
#'  the \code{mice} function. To be used in combination with validation methods cv_MI,
#'  cv_MI_RR and MI_boot. For method cv_MI the number of imputed datasets is fixed at 1 and cannot
#'  be changed.
#'
#' @details For internal validation five methods can be used, cv_MI, cv_MI_RR, MI_cv_naive,
#'  MI_boot and boot_MI. Method cv_MI uses imputation within each cross-validation fold definition.
#'  By repeating this in several imputation runs, multiply imputed datasets are generated. Method
#'  cv_MI_RR uses multiple imputation within the cross-validation definition. MI_cv_naive, applies
#'  cross-validation within each imputed dataset. MI_boot draws for each bootstrap step the same
#'  cases in all imputed datasets. With boot_MI first bootstrap samples are drawn from the original
#'  dataset with missing values and than multiple imputation is applied. For multiple imputation
#'  the \code{mice} function from the \code{mice} package is used. It is recommended to use a minumum
#'  of 100 imputation runs for method cv_MI or 100 bootstrap samples for method boot_MI or MI_boot.
#'  Methods cv_MI, cv_MI_RR and MI_cv_naive can be combined with backward selection during
#'  cross-validation and with methods boot_MI and MI_boot, backward and forward selection can
#'  be used. For methods cv_MI and cv_MI_RR the outcome in the original dataset has to be complete.
#'
#'@return A \code{psfmi_perform} object from which the following objects can be extracted: \code{res_boot},
#'  result of pooled performance (in multiply imputed datasets) at each bootstrap step of ROC app (pooled
#'  ROC), ROC test (pooled ROC after bootstrap model is applied in original multiply imputed datasets),
#'  same for R2 app (Nagelkerke's R2), R2 test, Scaled Brier app and Scaled Brier test. Information is also provided
#'  about testing the Calibration slope at each bootstrap step as interc test and Slope test.
#'  The performance measures are pooled by a call to the function \code{pool_performance}. Another
#'  object that can be extracted is \code{intval}, with information of the AUC, R2, Scaled Brier score and
#'  Calibration slope averaged over the bootstrap samples, in terms of: Orig (original datasets),
#'  Apparent (models applied in bootstrap samples), Test (bootstrap models are applied in original datasets),
#'  Optimism (difference between apparent and test) and Corrected (original corrected for optimism).
#'
#' @references Heymans MW, van Buuren S, Knol DL, van Mechelen W, de Vet HC. Variable selection under
#'  multiple imputation using the bootstrap in a prognostic study. BMC Med Res Methodol. 2007(13);7:33.
#' @references F. Harrell. Regression Modeling Strategies. With Applications to
#'  Linear Models, Logistic and Ordinal Regression, and Survival Analysis (2nd edition). Springer,
#'  New York, NY, 2015.
#' @references Van Buuren S. (2018). Flexible Imputation of Missing Data. 2nd Edition. Chapman &
#'  Hall/CRC Interdisciplinary Statistics. Boca Raton.
#' @references Harel, O. (2009). The estimation of R2 and adjusted R2 in
#'  incomplete data sets using multiple imputation. Journal of Applied Statistics,
#'  36(10), 1109-1118.
#' @references Musoro JZ, Zwinderman AH, Puhan MA, ter Riet G, Geskus RB. Validation of prediction
#'  models based on lasso regression with multiply imputed data. BMC Med Res Methodol. 2014;14:116.
#' @references Wahl S, Boulesteix AL, Zierer A, Thorand B, van de Wiel MA. Assessment of
#'  predictive performance in incomplete data by combining internal validation and multiple
#'  imputation. BMC Med Res Methodol. 2016;16(1):144.
#' @references EW. Steyerberg (2019). Clinical Prediction MOdels. A Practical Approach
#'  to Development, Validation, and Updating (2nd edition). Springer Nature Switzerland AG.
#'
#' @references http://missingdatasolutions.rbind.io/
#'
#' @section Vignettes:
#' 
#' \itemize{
#' \item  \href{https://mwheymans.github.io/psfmi/articles/cv_MI.html}{MI and Cross-validation - Method cv_MI}
#' \item  \href{https://mwheymans.github.io/psfmi/articles/cv_MI_RR.html}{MI and Cross-validation - Method cv_MI_RR}
#' \item  \href{https://mwheymans.github.io/psfmi/articles/MI_cv_naive.html}{MI and Cross-validation - Method MI_cv_naive} 
#' \item  \href{https://mwheymans.github.io/psfmi/articles/boot_MI.html}{MI and Bootstrapping - Method boot_MI}
#' \item  \href{https://mwheymans.github.io/psfmi/articles/MI_boot.html}{MI and Bootstrapping - Method MI_boot}
#' }
#' 
#' @author Martijn Heymans, 2020
#'
#' @examples
#' pool_lr <- psfmi_lr(data=lbpmilr, formula = Chronic ~ Pain + JobDemands + rcs(Tampascale, 3) +
#'            factor(Satisfaction) + Smoking, p.crit = 1, direction="FW",
#'            nimp=5, impvar="Impnr", method="D1")
#'            
#' pool_lr$RR_model
#'
#' res_perf <- psfmi_validate(pool_lr, val_method = "cv_MI", data_orig = lbp_orig, folds=3,
#'             nimp_cv = 2, p.crit=0.05, BW=TRUE, miceImp = miceImp, printFlag = FALSE)
#'             
#' res_perf
#'
#'\dontrun{
#'  set.seed(200)
#'   res_val <- psfmi_validate(pobj, val_method = "boot_MI", data_orig = lbp_orig, nboot = 5,
#'   p.crit=0.05, BW=TRUE, miceImp = miceImp, nimp_mice = 5, printFlag = FALSE, direction = "FW")
#'   
#'   res_val$stats_val
#'}
#'
#' @export
psfmi_validate.cox <- function(pobj, 
                               val_method = NULL, 
                               data_orig = NULL, 
                               int_val = TRUE, 
                               nboot = 10,
                               folds=3, 
                               nimp_cv = 5, 
                               nimp_mice = 5, 
                               p.crit = 1, 
                               BW = FALSE,
                               direction = NULL, 
                               cv_naive_appt=FALSE,
                               cal.plot=FALSE, 
                               plot.method="mean", 
                               groups_cal=5, 
                               miceImp, 
                               predictorMatrix,
                               seed,
                               horizon,
                               modelName,
                               miceMethod = NULL,
                               smcfcsMethod=FALSE,
                               smcfcsFormula=NULL,
                               keep.predictors=NULL,
                               strata.predictors=NULL,
                               ...)
{
  # # ##############################
  # pobj = pool_coxr.cox.nointer
  # val_method = "boot_MI"
  # data_orig = data.original.ext.cox
  # int_val = TRUE
  # nboot = myBoot
  # folds=NULL
  # nimp_cv = NULL
  # nimp_mice = 10
  # p.crit = 1
  # BW = FALSE
  # direction = NULL
  # cv_naive_appt=FALSE
  # cal.plot=FALSE
  # plot.method="mean"
  # groups_cal=5
  # miceImp = collectPooledMice
  # predictorMatrix =pred.cox
  # seed=3456
  # horizon=myHorizons
  # modelName="simple"
  # miceMethod=meth.cox
  # smcfcsMethod=FALSE
  # smcfcsFormula=NULL
  # keep.predictors=NULL
  # strata.predictors=NULL
  
  
  # pobj = pool_coxr
  # val_method = "boot_MI"
  # data_orig = data.original.ext
  # int_val = TRUE
  # nboot = 5
  # folds=NULL
  # nimp_cv = NULL
  # nimp_mice = 10
  # p.crit = 1
  # BW = FALSE
  # direction = NULL
  # cv_naive_appt=FALSE
  # cal.plot=FALSE
  # plot.method="mean"
  # groups_cal=5
  # miceImp = collectPooledMice4#collectPooledMice2
  # predictorMatrix =predMat
  # seed=3456
  # horizon=c(6,12,24,36)
  # modelName="simple"
  # miceMethod = NULL
  # smcfcsMethod=FALSE
  # smcfcsFormula=NULL
  # 
  
  # pobj = pool_coxr.cox.nointer
  # val_method = "boot_MI"
  # data_orig = data.original.ext.cox
  # int_val = TRUE
  # nboot = 1000
  # folds=NULL 
  # nimp_cv = NULL
  # nimp_mice = 10
  # p.crit = 1
  # BW = FALSE
  # direction = NULL
  # cv_naive_appt=FALSE
  # cal.plot=FALSE
  # plot.method="mean" 
  # groups_cal=5
  # miceImp = collectPooledMice
  # predictorMatrix =pred.cox
  # seed=3456
  # horizon=myHorizons
  # modelName="simple"
  # miceMethod=meth.cox
  # smcfcsMethod=FALSE
  # smcfcsFormula=NULL
  ##############################
  #General Settings
  # print("INTERNAL VALIDATION")
  call <- match.call()
  
  if(!inherits(pobj, "pmods"))
    stop("\n", "Object should be of type pmods", "\n")
  # if(pobj$model_type!="binomial")
  #   stop("\n", "Methods only available for models of type binomial", "\n")
  if(is_empty(pobj$predictors_final))
    stop("\n", "Model is empty. Cannot validate empty model", "\n")
  if(is.null(val_method))
    stop("\n", "Validation method not defined, choose boot_MI, MI_boot, cv_MI, cv_MI_RR or MI_cv_naive")
  
  # Check if you define smcfcs for MI you also need to provide methods and formuka
  if(isTRUE(smcfcsMethod)){
    if(is.null(smcfcsFormula) | is.null(miceMethod) | is_empty(miceMethod)){
      stop("\n", "Cannot run smcfcs MI without formula and methods, please check to provide both!")
      
    }
  }
  
  
  if(!int_val){ ### NEED TO CHECK THAT
    message("\n", "No validation - Apparent performance", "\n")
    Y <- c(paste(pobj$Outcome, paste("~")))
    if(is_empty(pobj$predictors_final)) {
      pobj$predictors_final <- 1
      fm <- as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))
    } else {
      fm <- as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))
    }
    
    perform_mi_orig <- pool_performance(data=pobj$data, nimp = pobj$nimp,
                                        impvar=pobj$impvar,
                                        formula = fm,
                                        cal.plot=cal.plot,
                                        plot.method=plot.method,
                                        groups_cal=groups_cal)
    return(perform_mi_orig)
  }
  if(int_val){
    
    # Specific for cv_MI methods
    if(val_method=="cv_MI" | val_method=="cv_MI_RR") {
      if(BW==FALSE & p.crit!=1)
        stop("\n", "If BW == FALSE, p.crit must be 1", "\n")
      if(BW & p.crit==1)
        stop("\n", "If BW == TRUE, p.crit must be < 1","\n")
      if(is_empty(data_orig))
        stop("\n", "data_orig not defined", "\n")
      if(any(is.na(pobj$Outcome)))
        stop("\n", "Outcome variable contains missing data,
             cv_MI or cv_MI_RR can only be used when outcome is complete", "\n")
      if(val_method=="cv_MI"){
        if(nimp_cv==10){
          message("\n", "Recommended to increase nimp_cv to 100", "\n")
        }
        pobjcv <- cv_MI(pobj, data_orig = data_orig, nimp_cv = nimp_cv,
                        folds = folds, p.crit = p.crit, BW=BW, miceImp = miceImp, ...)
        return(pobjcv)
      }
      
      if(val_method=="cv_MI_RR"){
        if(is_empty(nimp_mice))
          nimp_mice <- pobj$nimp
        pobjcv <- cv_MI_RR(pobj, data_orig = data_orig, nimp_mice = nimp_mice,
                           p.crit = p.crit, BW=BW, folds = folds, miceImp = miceImp, ...)
        return(pobjcv)
      }
    }
    
    if(val_method=="MI_cv_naive"){
      pobjcv <-
        MI_cv_naive(pobj, folds = folds, BW = BW, p.crit = p.crit, cv_naive_appt = cv_naive_appt)
      return(pobjcv)
    }
    
    # Specific for boot_MI and MI_boot
    if(val_method=="boot_MI" | val_method=="MI_boot") {
      message("\n", "Running BOOT-MI approach for the following time horizons: ", paste(horizon, collapse = ", "), "\n")
      if(val_method=="boot_MI"){
        # Part boot_MI
        if(is_empty(data_orig))
          stop("data_orig not defined")
        if(is_empty(nimp_mice))
          nimp_mice <- pobj$nimp
        if (nimp_mice < 2)
          stop("\n", "Number of imputed datasets for MI_boot must be > 1", "\n")
        # pobjboot <-
        #   boot_MI.cox(pobj, data_orig = data_orig, nboot = nboot,
        #           nimp_mice = nimp_mice, p.crit = p.crit, direction = direction, 
        #           miceImp = miceImp,seed=seed,
        #           horizon=horizon,
        #           modelName=modelName, ...)
        
        
        pobjboot <-
          boot_MI.cox(pobj= pobj,
                      data_orig = data_orig,
                      nboot = nboot,
                      nimp_mice = nimp_mice,
                      p.crit = p.crit,
                      direction = direction,
                      miceImp = miceImp,
                      predictorMatrix = predictorMatrix,
                      seed=seed,
                      horizon=horizon,
                      modelName=modelName,
                      miceMethod=miceMethod,
                      smcfcsMethod=smcfcsMethod,
                      smcfcsFormula=smcfcsFormula,
                      keep.predictors=keep.predictors,
                      strata.predictors=strata.predictors)
        
        
        if(p.crit==1)
          message("\n", "p.crit = 1, validation is done without variable selection", "\n")
        return(pobjboot)
      }
      if(val_method=="MI_boot"){
        pobjboot <-
          MI_boot(pobj, p.crit = p.crit, nboot = nboot, direction = direction)
        if(p.crit==1)
          message("\n", "p.crit = 1, validation is done without variable selection", "\n")
        return(pobjboot)
      }
    }
    cal.plot==FALSE
    message("\n", "Calibration plot not made, only when int.val = FALSE", "\n")
  }
}


collectPooledMice <- function(DF,...){
  # DF = dataDF.fix
  # numbM = 4
  # impMethod = defMethod1
  # myseed = 3456
  # myMat = mat.all
  imp <- mice(DF,...)
  imp
}



#' Function to evaluate bootstrap predictor and model stability in multiply imputed datasets.
#'
#' \code{psfmi_stab} Stability analysis of predictors and prediction models selected with
#'  the \code{psfmi_lr}, \code{psfmi_coxr} or \code{psfmi_mm} functions of the \code{psfmi} package.
#'
#' @param pobj An object of class \code{pmods} (pooled models), produced by a previous call to
#'  \code{psfmi_lr}, \code{psfmi_coxr} or \code{psfmi_mm}.
#' @param boot_method A single string to define the bootstrap method. Use "single" after a call to
#'  \code{psfmi_lr} and \code{psfmi_coxr} and "cluster" after a call to \code{psfmi_mm}.
#' @param nboot A numerical scalar. Number of bootstrap samples to evaluate the stability. Default is 20.
#' @param p.crit A numerical scalar. Used as P-value selection criterium during bootstrap model selection.
#' @param start_model If TRUE the bootstrap evaluation takes place from the start model of object pobj, if
#'  FALSE the final model is used for the evaluation.
#' @param direction The direction of predictor selection, "BW" for backward selection and "FW"
#'   for forward selection.
#'#'
#' @details The function evaluates predictor selection frequency in stratified or cluster bootstrap samples.
#'  The stratification factor is the variable that separates the imputed datasets. The same bootstrap cases
#'  are drawn in each bootstrap sample. It uses as input an object of class \code{pmods} as a result of a
#'  previous call to the \code{psfmi_lr}, \code{psfmi_coxr} or \code{psfmi_mm} functions.
#'  In combination with the \code{psfmi_mm} function a cluster bootstrap method is used where bootstrapping
#'  is used on the level of the clusters only (and not also within the clusters).
#'
#'@return A \code{psfmi_stab} object from which the following objects can be extracted: bootstrap
#'  inclusion (selection) frequency of each predictor \code{bif}, total number each predictor is
#'  included in the bootstrap samples as \code{bif_total}, percentage a predictor is selected
#'  in each bootstrap sample as \code{bif_perc} and number of times a prediction model is selected in
#'  the bootstrap samples as \code{model_stab}.
#'
#' @references Heymans MW, van Buuren S. et al. Variable selection under multiple imputation using the bootstrap
#'   in a prognostic study. BMC Med Res Methodol. 2007;13:7-33.
#' @references Eekhout I, van de Wiel MA, Heymans MW. Methods for significance testing of categorical
#'   covariates in logistic regression models after multiple imputation: power and applicability
#'   analysis. BMC Med Res Methodol. 2017;17(1):129.
#' @references Sauerbrei W, Schumacher M. A bootstrap resampling procedure for model building:
#'   application to the Cox regression model. Stat Med. 1992;11:2093-109.
#' @references Royston P, Sauerbrei W (2008) Multivariable model-building - a pragmatic approach to
#'   regression analysis based on fractional polynomials for modelling continuous variables. (2008).
#'   Chapter 8, Model Stability. Wiley, Chichester
#' @references Heinze G, Wallisch C, Dunkler D. Variable selection - A review and
#'  recommendations for the practicing statistician. Biom J. 2018;60(3):431-449.
#'
#' @references http://missingdatasolutions.rbind.io/
#'
#' @section Vignettes:
#'  https://mwheymans.github.io/psfmi/articles/psfmi_StabilityAnalysis.html
#'
#' @examples
#'  pool_lr <- psfmi_coxr(formula = Surv(Time, Status) ~ Pain + factor(Satisfaction) + 
#'    rcs(Tampascale,3) + Radiation + Radiation*factor(Satisfaction) + Age + Duration + 
#'    Previous + Radiation*rcs(Tampascale, 3), data=lbpmicox, p.crit = 0.157, direction="FW",
#'    nimp=5, impvar="Impnr", keep.predictors = NULL, method="D1")
#'
#'  pool_lr$RR_Model
#'  pool_lr$multiparm
#'
#' \dontrun{
#'  stab_res <- psfmi_stab(pool_lr, direction="FW", start_model = TRUE,
#'      boot_method = "single", nboot=20, p.crit=0.05)
#'  stab_res$bif
#'  stab_res$bif_perc
#'  stab_res$model_stab
#'}
#'
#' @export
my.psfmi_stab <- function(pobj, boot_method=NULL, nboot=20,
                       p.crit = 0.05, start_model = TRUE, direction = NULL)
{
  # pobj=pool_lr
  # direction="BW"
  # start_model = TRUE
  # boot_method = "single"
  # nboot=20
  # p.crit=0.05
  
  if(!inherits(pobj, "pmods"))
    stop("\n", "Object should be of type pmods", "\n")
  if(is.null(boot_method))
    stop("\n", "boot_method is not defined, choose single or cluster", "\n")
  if(p.crit==1)
    stop("\n", "To determine Model Stability p.crit must be < 1", "\n\n")
  
  if(boot_method=="single" & pobj$model_type=="binomial" | pobj$model_type=="survival"){
    if(is_empty(direction))
      stop( "\n", "Specify FW or BW for forward or backward predictor selection", "\n")
    if(start_model == FALSE & is_empty(pobj$predictors_final))
      stop( "\n", "Final model is empty. You cannot determine the stability of an empty model", "\n\n")
    if(pobj$model_type=="binomial"){
      if(pobj$direction=="FW" & start_model == FALSE){
        Y <-
          c(paste(pobj$Outcome, paste("~")))
        pobj$formula_final <-
          as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))
      }
    }
    if(pobj$model_type=="survival"){
      if(pobj$direction=="FW" & start_model == FALSE){
        Y <-
          c(paste0("Surv(", pobj$time, ",", pobj$status, ")~"))
        pobj$formula_final <-
          as.formula(paste(Y, paste(pobj$predictors_final, collapse = "+")))
      }
    }
  }
  call <- match.call()
  nboot <- nboot
  data <- pobj$data
  boot_seq <- as.list(1:nboot)
  
  if (inherits(pobj, "pmods") & boot_method == "single"){
    
    boot_data <- bootstraps(data, strata = pobj$impvar, times = nboot)
    boot_pred_pat <- mapply(function(x, y) {
      #######################
      x<-boot_data$splits[[1]]
      y<-boot_seq[[1]]
      #######################
      message("\n", "Boot ", y)
      x <- as.data.frame(x)
      
      if(pobj$model_type=="binomial") {
        
        if(start_model){
          psfmi_boot <- psfmi_lr(formula = pobj$formula_initial, data = x, nimp=pobj$nimp, impvar = pobj$impvar,
                                 p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                                 method = pobj$method, direction = direction)
          # print(psfmi_boot)
          
        } else {
          psfmi_boot <- psfmi_lr(formula = pobj$formula_final, data = x, nimp=pobj$nimp, impvar = pobj$impvar,
                                 p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                                 method = pobj$method, direction = direction)
        }
      }
      if(pobj$model_type=="survival") {
        
        if(start_model){
          psfmi_boot <- psfmi_coxr.Fixed(formula = pobj$formula_initial, data = x, nimp=pobj$nimp, impvar = pobj$impvar,
                                   p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                                   method = pobj$method, direction = direction)
          
        } else {
          psfmi_boot <- psfmi_coxr.Fixed(formula = pobj$formula_final, data = x, nimp=pobj$nimp, impvar = pobj$impvar,
                                   p.crit = p.crit, keep.predictors = pobj$keep.predictors,
                                   method = pobj$method, direction = direction)
        }
      }
      
      if(direction=="BW")
        boot_predictors_in <- ifelse(psfmi_boot$predictors_out[nrow(psfmi_boot$predictors_out), ], 0, 1)
      if(direction=="FW")
        boot_predictors_in <- psfmi_boot$predictors_in[nrow(psfmi_boot$predictors_in), ]
      return(boot_predictors_in)
    }, x = boot_data$splits, y=boot_seq, SIMPLIFY = FALSE)
    
    bif <- data.frame(do.call("rbind", boot_pred_pat))
    
    if(!start_model)
      colnames(bif) <- pobj$predictors_final
    if(start_model)
      colnames(bif) <- pobj$predictors_initial
    
    # Group selected models
    bif_pat <- bif %>%
      group_by_all %>%
      count
    
    # desc order of selected models
    bif_pat_sort <- data.frame(bif_pat %>%
                                 arrange(desc(n)) %>%
                                 select(1:ncol(bif_pat)))
    bif_pat_perc <- round((bif_pat_sort$n / nboot) * 100, 0)
    bif_pat_sort <- data.frame(bif_pat_sort, bif_pat_perc)
    
    if(!start_model) {
      colnames(bif_pat_sort) <- c(pobj$predictors_final, "freq", "bif_pat_perc")
    } else {
      colnames(bif_pat_sort) <- c(pobj$predictors_initial, "freq", "bif_pat_perc")
    }
    rownames(bif) <- paste("boot", 1:nboot)
    
    bif_total <- colSums(bif)
    bif_perc <- round((bif_total / nboot) * 100, 3)
    
    stabobj <- list(bif = bif, bif_total = bif_total, bif_perc = bif_perc,
                    model_stab = bif_pat_sort, call = call)
    return(stabobj)
  }
  
  if(inherits(pobj, "pmods") & boot_method == "cluster")
  {
    if(is_empty(pobj$clusvar))
      stop("\n", "No cluster variable defined, use single bootstrapping method")
    cluster_bootmi <- function(data, impvar, clusvar, nimp, nboot)
    {
      boot_imp <- list()
      for(i in 1:nimp){
        
        # Extract each imputed dataset
        imp_data <- subset(data, data[impvar] == i)
        nest_data <- imp_data %>% nest(data=c(-which(names(imp_data) == clusvar)))
        nest_data
        
        bs_data <- bootstraps(nest_data, times = nboot)
        bs_data
        
        df_boot <-lapply(bs_data$splits,
                         function(x) {
                           data_boot <- as_tibble(x) %>%
                             unnest(cols = c(clusvar))
                           n_row <- map_int(data_boot$data, ~ sapply(.x[1], NROW))
                           df_booti <- as_tibble(data.frame(clus_id=rep(data_boot[[clusvar]], n_row),
                                                            do.call("rbind", data_boot$data)))
                           names(df_booti)[1] <- clusvar
                           return(df_booti)
                         })
        boot_imp[[i]] <- df_boot
      }
      
      # binds first elements of lists
      dfs_bootmi <- transpose(boot_imp) %>% map(bind_rows)
      
      bootnr <- (1:nboot)
      nrow_boot <- map_int(dfs_bootmi, ~ sapply(.x[1], NROW))
      dfs_bootmi_tot <- mapply(function(x, y, z) {
        idboot <- rep(x, y)
        data.frame(idboot, z)},
        x=bootnr, y=nrow_boot, z=dfs_bootmi, SIMPLIFY = FALSE)
      return(dfs_bootmi_tot)
    }
    
    boot_clus_res <- cluster_bootmi(data=pobj$data, impvar = pobj$impvar,
                                    clusvar = pobj$clusvar, nimp = pobj$nimp, nboot = nboot)
    
    boot_pred_pat <- mapply(function(x, y) {
      message("\n", "Boot ", y)
      x <- as.data.frame(x)
      psfmi_boot <- psfmi_mm(data=x, nimp=pobj$nimp, impvar = pobj$impvar, random.eff = pobj$random.eff,
                             Outcome = pobj$Outcome, predictors = pobj$predictors, family = pobj$family,
                             p.crit = pobj$p.crit, cat.predictors = pobj$cat.predictors,
                             spline.predictors = pobj$spline.predictors, clusvar = pobj$clusvar,
                             int.predictors = pobj$int.predictors, keep.predictors = pobj$keep.predictors,
                             nknots = pobj$nknots, method = pobj$method, print.method = pobj$print.method)
      
      boot_predictors_in <- psfmi_boot$predictors_in[nrow(psfmi_boot$predictors_in), ]
      return(boot_predictors_in)
    }, x = boot_clus_res, y=boot_seq, SIMPLIFY = FALSE)
    
    
    bif <- data.frame(do.call("rbind", boot_pred_pat))
    
    if(!start_model)
      colnames(bif) <- pobj$predictors_final
    if(start_model)
      colnames(bif) <- pobj$predictors_initial
    
    names_temp <- colnames(bif)
    # Group selected models
    bif_pat <- bif %>%
      group_by_all() #%>%
    n <- count(bif_pat)$n
    bif_pat <- count(bif_pat)
    bif_pat <- data.frame(bif_pat, n)[, -length(colnames(bif_pat))]
    colnames(bif_pat) <- c(names_temp, "n")
    #names_temp <- names(bif_pat)
    
    bif_pat_sort <- arrange(bif_pat, desc(n))
    #bif_pat_sort <- data.frame(bif_pat %>%
    #                             arrange(desc(n)) %>%
    #                             select(1:ncol(bif_pat)))
    bif_pat_perc <- round((bif_pat_sort$n / nboot) * 100, 0)
    bif_pat_sort <- data.frame(bif_pat_sort, bif_pat_perc)
    
    bif_pat_perc <- round((bif_pat_sort$n / nboot) * 100, 0)
    bif_pat_sort <- data.frame(bif_pat_sort, bif_pat_perc)
    colnames(bif_pat_sort) <- c(names(pobj$predictors_in), "freq", "bif_pat_perc")
    
    rownames(bif) <- paste("boot", 1:nboot)
    
    bif_total <- colSums(bif)
    bif_perc <- round((bif_total / nboot) * 100, 3)
    
    stabobj <- list(bif = bif, bif_total = bif_total, bif_perc = bif_perc,
                    model_stab = bif_pat_sort, call = call)
    return(stabobj)
  }
}



my.psfmi_stab2 <- function (pobj, boot_method = NULL, nboot = 20, p.crit = 0.05, 
          start_model = TRUE, direction = NULL) 
{
  
  # pobj=pool_lr
  # direction="BW"
  # start_model = TRUE
  # boot_method = "single"
  # nboot=5
  # p.crit=0.05
  
  
  if (!inherits(pobj, "pmods")) 
    stop("\n", "Object should be of type pmods", 
         "\n")
  if (is.null(boot_method)) 
    stop("\n", "boot_method is not defined, choose single or cluster", 
         "\n")
  if (p.crit == 1) 
    stop("\n", "To determine Model Stability p.crit must be < 1", 
         "\n\n")
  if (boot_method == "single" & pobj$model_type == "binomial" | 
      pobj$model_type == "survival") {
    if (is_empty(direction)) 
      stop("\n", "Specify FW or BW for forward or backward predictor selection", 
           "\n")
    if (start_model == FALSE & is_empty(pobj$predictors_final)) 
      stop("\n", "Final model is empty. You cannot determine the stability of an empty model", 
           "\n\n")
    if (pobj$model_type == "binomial") {
      if (pobj$direction == "FW" & start_model == 
          FALSE) {
        Y <- c(paste(pobj$Outcome, paste("~")))
        pobj$formula_final <- as.formula(paste(Y, paste(pobj$predictors_final, 
                                                        collapse = "+")))
      }
    }
    if (pobj$model_type == "survival") {
      if (pobj$direction == "FW" & start_model == 
          FALSE) {
        Y <- c(paste0("Surv(", pobj$time, ",", 
                      pobj$status, ")~"))
        pobj$formula_final <- as.formula(paste(Y, paste(pobj$predictors_final, 
                                                        collapse = "+")))
      }
    }
  }
  call <- match.call()
  nboot <- nboot
  data <- pobj$data
  boot_seq <- as.list(1:nboot)
  if (inherits(pobj, "pmods") & boot_method == "single") {
    boot_data <- bootstraps(data, strata = pobj$impvar, times = nboot)
    boot_pred_pat <- mapply(function(x, y) {
      
      #######################
      # x<-boot_data$splits[[1]]
      # y<-boot_seq[[1]]
      # #######################
      message("\n", "Boot ", y)
      x <- as.data.frame(x)
      if (pobj$model_type == "binomial") {
        if (start_model) {
          psfmi_boot <- psfmi_lr(formula = pobj$formula_initial, 
                                 data = x, nimp = pobj$nimp, impvar = pobj$impvar, 
                                 p.crit = p.crit, keep.predictors = pobj$keep.predictors, 
                                 method = pobj$method, direction = direction)
        }
        else {
          psfmi_boot <- psfmi_lr(formula = pobj$formula_final, 
                                 data = x, nimp = pobj$nimp, impvar = pobj$impvar, 
                                 p.crit = p.crit, keep.predictors = pobj$keep.predictors, 
                                 method = pobj$method, direction = direction)
        }
      }
      if (pobj$model_type == "survival") {
        if (start_model) {
          psfmi_boot <- psfmi_coxr(formula = pobj$formula_initial, 
                                   data = x, nimp = pobj$nimp, impvar = pobj$impvar, 
                                   p.crit = p.crit, keep.predictors = pobj$keep.predictors, 
                                   method = pobj$method, direction = direction)
        }
        else {
          psfmi_boot <- psfmi_coxr(formula = pobj$formula_final, 
                                   data = x, nimp = pobj$nimp, impvar = pobj$impvar, 
                                   p.crit = p.crit, keep.predictors = pobj$keep.predictors, 
                                   method = pobj$method, direction = direction)
        }
      }
      if (direction == "BW") 
        boot_predictors_in <- ifelse(psfmi_boot$predictors_out[nrow(psfmi_boot$predictors_out), 
        ], 0, 1)
      if (direction == "FW") 
        boot_predictors_in <- psfmi_boot$predictors_in[nrow(psfmi_boot$predictors_in), 
        ]
      return(boot_predictors_in)
    }, x = boot_data$splits, y = boot_seq, SIMPLIFY = FALSE)
    bif <- data.frame(do.call("rbind", boot_pred_pat))
    if (!start_model) 
      colnames(bif) <- pobj$predictors_final
    if (start_model) 
      colnames(bif) <- pobj$predictors_initial
    bif_pat <- bif %>% group_by_all %>% dplyr::count()
    bif_pat_sort <- data.frame(bif_pat %>% arrange(desc(n)) %>% 
                                 select(1:ncol(bif_pat)))
    bif_pat_perc <- round((bif_pat_sort$n/nboot) * 100, 0)
    bif_pat_sort <- data.frame(bif_pat_sort, bif_pat_perc)
    if (!start_model) {
      colnames(bif_pat_sort) <- c(pobj$predictors_final, 
                                  "freq", "bif_pat_perc")
    }
    else {
      colnames(bif_pat_sort) <- c(pobj$predictors_initial, 
                                  "freq", "bif_pat_perc")
    }
    rownames(bif) <- paste("boot", 1:nboot)
    bif_total <- colSums(bif)
    bif_perc <- round((bif_total/nboot) * 100, 3)
    stabobj <- list(bif = bif, bif_total = bif_total, bif_perc = bif_perc, 
                    model_stab = bif_pat_sort, call = call)
    return(stabobj)
  }
  if (inherits(pobj, "pmods") & boot_method == "cluster") {
    if (is_empty(pobj$clusvar)) 
      stop("\n", "No cluster variable defined, use single bootstrapping method")
    cluster_bootmi <- function(data, impvar, clusvar, nimp, 
                               nboot) {
      boot_imp <- list()
      for (i in 1:nimp) {
        imp_data <- subset(data, data[impvar] == i)
        nest_data <- imp_data %>% nest(data = c(-which(names(imp_data) == 
                                                         clusvar)))
        nest_data
        bs_data <- bootstraps(nest_data, times = nboot)
        bs_data
        df_boot <- lapply(bs_data$splits, function(x) {
          data_boot <- as_tibble(x) %>% unnest(cols = c(clusvar))
          n_row <- map_int(data_boot$data, ~sapply(.x[1], 
                                                   NROW))
          df_booti <- as_tibble(data.frame(clus_id = rep(data_boot[[clusvar]], 
                                                         n_row), do.call("rbind", data_boot$data)))
          names(df_booti)[1] <- clusvar
          return(df_booti)
        })
        boot_imp[[i]] <- df_boot
      }
      dfs_bootmi <- transpose(boot_imp) %>% map(bind_rows)
      bootnr <- (1:nboot)
      nrow_boot <- map_int(dfs_bootmi, ~sapply(.x[1], NROW))
      dfs_bootmi_tot <- mapply(function(x, y, z) {
        idboot <- rep(x, y)
        data.frame(idboot, z)
      }, x = bootnr, y = nrow_boot, z = dfs_bootmi, SIMPLIFY = FALSE)
      return(dfs_bootmi_tot)
    }
    boot_clus_res <- cluster_bootmi(data = pobj$data, impvar = pobj$impvar, 
                                    clusvar = pobj$clusvar, nimp = pobj$nimp, nboot = nboot)
    boot_pred_pat <- mapply(function(x, y) {
      message("\n", "Boot ", y)
      x <- as.data.frame(x)
      psfmi_boot <- psfmi_mm(data = x, nimp = pobj$nimp, 
                             impvar = pobj$impvar, random.eff = pobj$random.eff, 
                             Outcome = pobj$Outcome, predictors = pobj$predictors, 
                             family = pobj$family, p.crit = pobj$p.crit, cat.predictors = pobj$cat.predictors, 
                             spline.predictors = pobj$spline.predictors, clusvar = pobj$clusvar, 
                             int.predictors = pobj$int.predictors, keep.predictors = pobj$keep.predictors, 
                             nknots = pobj$nknots, method = pobj$method, print.method = pobj$print.method)
      boot_predictors_in <- psfmi_boot$predictors_in[nrow(psfmi_boot$predictors_in), 
      ]
      return(boot_predictors_in)
    }, x = boot_clus_res, y = boot_seq, SIMPLIFY = FALSE)
    bif <- data.frame(do.call("rbind", boot_pred_pat))
    if (!start_model) 
      colnames(bif) <- pobj$predictors_final
    if (start_model) 
      colnames(bif) <- pobj$predictors_initial
    names_temp <- colnames(bif)
    bif_pat <- bif %>% group_by_all()
    n <- count(bif_pat)$n
    bif_pat <- count(bif_pat)
    bif_pat <- data.frame(bif_pat, n)[, -length(colnames(bif_pat))]
    colnames(bif_pat) <- c(names_temp, "n")
    bif_pat_sort <- arrange(bif_pat, desc(n))
    bif_pat_perc <- round((bif_pat_sort$n/nboot) * 100, 0)
    bif_pat_sort <- data.frame(bif_pat_sort, bif_pat_perc)
    bif_pat_perc <- round((bif_pat_sort$n/nboot) * 100, 0)
    bif_pat_sort <- data.frame(bif_pat_sort, bif_pat_perc)
    colnames(bif_pat_sort) <- c(names(pobj$predictors_in), 
                                "freq", "bif_pat_perc")
    rownames(bif) <- paste("boot", 1:nboot)
    bif_total <- colSums(bif)
    bif_perc <- round((bif_total/nboot) * 100, 3)
    stabobj <- list(bif = bif, bif_total = bif_total, bif_perc = bif_perc, 
                    model_stab = bif_pat_sort, call = call)
    return(stabobj)
  }
}