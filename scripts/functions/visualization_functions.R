plotMeasures_byTime <- function(modelList,modelNames, formula, data, horizons,metric="brier",timesLabel="Months",filterModels=NULL){
  # Applies to AUC, Brier and IPA
  # Can be used to compare competing models
  # Allows to filter out models from list based on their names, by default keeps the null model for brier, so needs to be manually specified  "Null model"
  
  # modelList= list(fullfit)
  # modelNames=c("simple")
  # formula = modelFormula(surv_time, surv_status, 1)
  # data = data.dev
  # horizons = c(6,12,18,24,30,36,42)
  # metric="ipa"
  # timesLabel="Months"
  # filterModels="Null model"
  ## START  
  # modelList=list(model)
  names(modelList) <-modelNames
  # Apparent Brier and IPA)
  scores_brier <- riskRegression::Score(modelList,
                                        formula = formula, 
                                        data = data, 
                                        conf.int = TRUE, 
                                        times = horizons,
                                        cens.method = "ipcw",
                                        cens.model = "cox",
                                        metrics = "brier",
                                        summary = "ipa"
  )
  
  
  scores_auc <- riskRegression::Score(modelList,
                                      formula = formula, 
                                      data = data, 
                                      conf.int = TRUE, 
                                      times = horizons,
                                      cens.method = "ipcw",
                                      cens.model = "cox",
                                      metrics = "AUC"
  )
  
  
  
  if(metric=="auc2"){
    # Default method for plotting AUC from riskRegression package
    p <- autoplot(scores_auc)
    p
  }else{
    
    if(metric=="brier"){
      scores <- scores_brier$Brier$score
      yVar="Brier"
      yVar.un<-as.name(yVar)
      yLabel<-"Brier score"
      scores$Brier <- scores$Brier*100
      scores$lower <- scores$lower*100
      scores$upper <- scores$upper*100
      
    }else if(metric=="ipa"){
      scores <- scores_brier$Brier$score
      yVar="IPA"
      yVar.un<-as.name(yVar)
      yLabel<-"IPA score"
      scores$IPA <- scores$IPA*100
      scores$lower <- scores$lower*100
      scores$upper <- scores$upper*100
    }else if(metric=="auc"){
      scores <- scores_auc$AUC$score
      yVar="AUC"
      yVar.un<-as.name(yVar)
      yLabel<-"AUC"
      scores$AUC <- scores$AUC*100
      scores$lower <- scores$lower*100
      scores$upper <- scores$upper*100
    }
    
    if(!is.null(filterModels)){
      scores <- scores %>% dplyr::filter(model %notin% filterModels)
    }
    
    p <- scores %>%
      ggplot( aes(x=times, y=!!yVar.un, group=model, color=model)) +
      geom_line(size=1) +
      scale_color_manual(values = ccf_palette("main")) + 
      theme_ipsum() +
      ylab(yLabel) + 
      xlab(timesLabel) + 
      labs(color="Model") + 
      scale_x_continuous(limits=c(horizons[1], horizons[length(horizons)]),breaks=horizons)
    
    
    if(metric=="brier"){
      p <- p + scale_y_continuous(labels = scales::percent_format(scale = 1),limits=c(0,100),breaks = c(0,10,20,30,40,50,60,70,80,90,100))
      
    }else if(metric=="ipa"){
      p <- p + scale_y_continuous(labels = scales::percent_format(scale = 1),limits=c(0,100),breaks = c(0,10,20,30,40,50,60,70,80,90,100))
    }else if(metric=="auc"){
      p <- p + scale_y_continuous(labels = scales::percent_format(scale = 1),limits=c(50,100),breaks = c(50,60,70,80,90,100))
      
    }
    
    
    
    p <- p + theme_bw(base_family = "Palatino", base_size = 15) +
      theme(panel.border = element_rect(colour = "black", linewidth=0.2),
            panel.background = element_rect(fill = "white"),
            panel.grid=element_blank(),
            axis.text = element_text(size = 10),
            axis.title=element_text(size = 14, face="bold"),
            text=element_text(family="Palatino Linotype"),
            legend.text=element_text(size = 8),
            legend.title=element_text(size = 10),
            legend.justification = 'center', 
            legend.position = 'bottom', 
            legend.box = 'horizontal')
    
    p
  }
}



plotMeasures_byTimeBoot <- function(scores, horizons,metric="brier",timesLabel="Months",palleteName="main",filterModels=NULL, scale=TRUE,vlines=NULL){
  # Applies to AUC, Brier, IPA, C.index
  # Can be used to compare competing models
  # Allows to filter out models from list based on their names, by default keeps the null model for brier, so needs to be manually specified  "Null model"
  # Input is a dataframe with calculated scores over time points, can include multiple models, model names in column called model - input ideally is after
  # validating each model with bootstraps to get optimism corrected indexes
  
  # scores = indexes.time.all
  # horizons = c(6,12,18,24,30,36,42)
  # metric="c"
  # timesLabel="Months"
  # filterModels="Null model"
  
  
  # scores = indexes.simple
  # horizons = c(6,12,24,36)
  # metric="brier"
  # timesLabel="Months"
  # filterModels=NULL
  
  # START
  scores$times <- as.numeric(as.character(scores$times))
  
  
  if(metric=="brier"){
    yVar="Brier"
    yVar.un<-as.name(yVar)
    yLabel<-"Brier score"
    scores$Brier <- scores$Brier*100
    
    
  }else if(metric=="ipa"){
    yVar="IPA"
    yVar.un<-as.name(yVar)
    yLabel<-"IPA score"
    scores$IPA <- scores$IPA*100
    
  }else if(metric=="auc"){
    yVar="AUC"
    yVar.un<-as.name(yVar)
    yLabel<-"AUC"
    # scores$AUC <- scores$AUC*100
    scores$AUC <- scores$AUC
    
  }else if(metric=="c"){
    yVar="Cindex"
    yVar.un<-as.name(yVar)
    yLabel<-"Harrell's C.index"
    
  }
  
  if(!is.null(filterModels)){
    scores <- scores %>% dplyr::filter(model %notin% filterModels)
  }
  
  p <- scores %>%
    ggplot( aes(x=times, y=!!yVar.un, group=model, color=model)) +
    geom_line(size=.8) +
    scale_color_manual(values = ccf_palette(name=palleteName,type = "discrete") %>% rev()) + 
    theme_ipsum() +
    ylab(yLabel) + 
    xlab(timesLabel) + 
    labs(color="Model") + 
    scale_x_continuous(limits=c(horizons[1], horizons[length(horizons)]),breaks=horizons)
  
  
  if(isTRUE(scale)){
    if(metric=="brier"){
      p <- p + scale_y_continuous(labels = scales::percent_format(scale = 1),limits=c(0,100),breaks = c(0,10,20,30,40,50,60,70,80,90,100))
      
    }else if(metric=="ipa"){
      p <- p + scale_y_continuous(labels = scales::percent_format(scale = 1),limits=c(0,50),breaks = c(0,10,20,30,40,50))
    }else if(metric=="auc"){
      # p <- p + scale_y_continuous(labels = scales::percent_format(scale = 1),limits=c(50,100),breaks = c(50,60,70,80,90,100))
      p <- p + scale_y_continuous(limits=c(.50,1),breaks = c(.50,.60,.70,.80,.90,1))
      
    }else if(metric=="c"){
      p <- p + scale_y_continuous(limits=c(.5,1),breaks = seq(.5,1,by=0.1))
      
    }
    
  }
  
  if(!is.null(vlines)){
    p <- p + geom_vline(xintercept = vlines, linetype="dotted", 
                          color = "black", size=.7)
    p <- p + geom_hline(yintercept = .8, linetype="dotted", 
                        color = "grey", size=.7)
  }
  
  
  p <- p + theme_bw(base_family = "Palatino", base_size = 15) +
    theme(panel.border = element_rect(colour = "black", linewidth=0.2),
          panel.background = element_rect(fill = "white"),
          panel.grid=element_blank(),
          axis.text = element_text(size = 12),
          axis.title=element_text(size = 14, face="bold"),
          text=element_text(family="Palatino"),
          legend.text=element_text(size = 11),
          legend.title=element_text(size = 12),
          legend.justification = 'center', 
          legend.position = 'bottom', 
          legend.box = 'horizontal')
  
  p
  
}



getLegendData2 <- function(object,
                           models,
                           times,
                           brier.in.legend=TRUE,
                           ipa.in.legend=TRUE,
                           format.brier,
                           auc.in.legend=TRUE,
                           format.auc,
                           drop.null.model=TRUE,
                           scale=100,
                           digits=1,
                           ...){
  model=AUC=lower=upper=Brier=NULL
  if (missing(models)) {
    models <- names(object$models)
  }
  if(!is.null(object$null.model)){
    if (drop.null.model==TRUE) models <- models[models!=object$null.model]
  }
  maxlen <- max(nchar(as.character(models)))
  legend.text.models <- sprintf(paste0("%",maxlen,"s"),models)
  if (missing(format.auc))
    format.auc <- paste0("%1.",digits,"f [%1.",digits,"f;%1.",digits,"f]")
  if (missing(format.brier))
    format.brier <- paste0("%1.",digits,"f [%1.",digits,"f;%1.",digits,"f]")
  if (is.null(object$null.model) || drop.null.model[[1]]==TRUE){
    keep.null.model <- FALSE
  }else{
    if (brier.in.legend==TRUE){
      keep.null.model=TRUE
    }else{
      keep.null.model=match(object$null.model,models,nomatch=FALSE)
    }
  }
  if (auc.in.legend==TRUE){
    if (is.null(object$AUC)){
      warning("Cannot show AUC as it is not stored in object. Set metrics='auc' in the call of Score.")
      legend.text.auc <- NULL
    }else{
      auc.data <- object$AUC$score[(model%in%models)]
      if (object$response.type!="binary"){
        if (missing(times)){
          tp <- max(auc.data[["times"]])
          if (length(unique(auc.data$times))>1)
            warning("Time point not specified, use max of the available times: ",tp)
        } else{ ## can only do one time point
          tp <- times[[1]]
          if (!(tp%in%unique(auc.data$times)))
            stop(paste0("Requested time ",times[[1]]," is not in object"))
        }
        auc.data <- auc.data[times==tp]
      }else tp <- NULL
      if (keep.null.model==FALSE){
        if(!is.null(object$null.model)){
          auc.data <- auc.data[model!=object$null.model]
        }
      }
      ## user's order
      auc.data[,model:=factor(model,levels=models)]
      setkey(auc.data,model)
      legend.text.auc <- auc.data[,sprintf(fmt=format.auc,scale*AUC,scale*lower,scale*upper)]
    }
  }else{
    legend.text.auc <- NULL
  }
  if (brier.in.legend==TRUE){
    if (is.null(object$Brier)){
      warning("Cannot show Brier score as it is not stored in object. Set metrics='brier' in the call of Score.")
      legend.text.brier <- NULL
    }else{
      brier.data <- object$Brier$score[(model%in%models)]
      if (object$response.type!="binary"){
        if (missing(times)){
          tp <- max(brier.data[["times"]])
          if (length(unique(brier.data$times))>1)
            warning("Time point not specified, use max of the available times: ",tp)
        } else{ ## can only do one time point
          tp <- times[[1]]
          if (!(tp%in%unique(brier.data$times)))
            stop(paste0("Requested time ",times[[1]]," is not in object"))
        }
        brier.data <- brier.data[times==tp]
      }else tp <- NULL
      if (!is.null(object$null.model) && keep.null.model[[1]]==FALSE){
        brier.data <- brier.data[model!=object$null.model]
      }
      brier.data[,model:=factor(model,levels=models)]
      setkey(brier.data,model)
      legend.text.brier <- brier.data[,sprintf(fmt=format.brier,scale*Brier,scale*lower,scale*upper)]
    }
  }else{
    legend.text.brier <- NULL
  }
  
  if (ipa.in.legend==TRUE){
    # print("ADDING IPA")
    if (is.null(object$Brier$score$IPA)){
      warning("Cannot show IPA score as it is not stored in object. Set metrics='brier' and summary='ipa' in the call of Score.")
      legend.text.ipa <- NULL
    }else{
      # ipa.data <- object$Brier$score[(model%in% c(models))]
      # print(ipa.data)
      # print(models)
      ipa.data <- object$Brier$score[(model%in% levels(models))]
      # print(ipa.data)
      ipa.data$IPA_lower <- 1-(ipa.data$upper/ipa.data$upper[1])
      ipa.data$IPA_upper <- 1-(ipa.data$lower/ipa.data$lower[1])
      ipa.data <- ipa.data[(model%in% c(models))]
      # print(ipa.data)
      
      if (object$response.type!="binary"){
        if (missing(times)){
          tp <- max(ipa.data[["times"]])
          if (length(unique(ipa.data$times))>1)
            warning("Time point not specified, use max of the available times: ",tp)
        } else{ ## can only do one time point
          tp <- times[[1]]
          if (!(tp%in%unique(ipa.data$times)))
            stop(paste0("Requested time ",times[[1]]," is not in object"))
        }
        ipa.data <- ipa.data[times==tp]
      }else tp <- NULL
      if (!is.null(object$null.model) && keep.null.model[[1]]==FALSE){
        ipa.data <- ipa.data[model!=object$null.model]
      }
      ipa.data[,model:=factor(model,levels=models)]
      setkey(brier.data,model)
      legend.text.ipa <- ipa.data[,sprintf(fmt=format.brier,scale*IPA,scale*IPA_lower,scale*IPA_upper)]
    }
  }else{
    legend.text.ipa <- NULL
  }
  
  out <- cbind(legend.text.models,
               AUC=legend.text.auc,
               Brier=legend.text.brier,
               IPA=legend.text.ipa)
  out
}


## plotROC.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun 23 2016 (10:27) 
## Version: 
## last-updated: Mar  9 2022 (14:33) 
##           By: Thomas Alexander Gerds
##     Update #: 184
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Plot ROC curve 
##'
##' @title Plot ROC curves
#' @param x Object obtained with function \code{Score}
#' @param models Choice of models to plot
#' @param times Time point(s) specifying the prediction horizon
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#' @param col line color
#' @param lwd line width
#' @param lty line style
#' @param cex point size
#' @param pch point style
#' @param legend logical. If \code{1L} draw a legend with the values
#'     of AUC.
#' @param auc.in.legend Logical. If \code{TRUE} add AUC to legend.
#' @param brier.in.legend Logical. If \code{TRUE} add Brier score to
#'     legend.
#' @param add logical. If \code{1L} add lines to an existing plot.
#' @param ... Used for additional control of the subroutines: plot,
#'     axis, lines, legend, addtable2plot. See \code{\link{SmartControl}}.
##' @examples
##' ## binary
##' set.seed(18)
##' if (require("randomForest",quietly=TRUE)){
##' library(randomForest)
##' library(prodlim)
##' bdl <- sampleData(40,outcome="binary")
##' bdt <- sampleData(58,outcome="binary")
##' bdl[,y:=factor(Y)]
##' bdt[,y:=factor(Y)]
##' fb1 <- glm(y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=bdl,family="binomial")
##' fb2 <- randomForest(y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=bdl)
##' xb <- Score(list("glm"=fb1,"rf"=fb2),y~1,data=bdt,
##'             plots="roc",metrics=c("auc","brier"))
##' plotROC(xb,brier.in.legend=1L)
##' 
##' # with cross-validation
##' \dontrun{
##' xb3 <- Score(list("glm"=fb1,"rf"=fb2),y~1,data=bdl,
##'             plots="roc",B=3,split.method="bootcv",
##'             metrics=c("auc"))
##' }
##' }
##' ## survival
##' set.seed(18)
##' library(survival)
##' sdl <- sampleData(40,outcome="survival")
##' sdt <- sampleData(58,outcome="survival")
##' fs1 <- coxph(Surv(time,event)~X3+X5+X6+X7+X8+X10,data=sdl,x=TRUE)
##' fs2 <- coxph(Surv(time,event)~X1+X2+X9,data=sdl,x=TRUE)
##' xs <- Score(list(model1=fs1,model2=fs2),Hist(time,event)~1,data=sdt,
##'             times=5,plots="roc",metrics="auc")
##' plotROC(xs)
##' ## competing risks
##' data(Melanoma)
##' f1 <- CSC(Hist(time,status)~age+sex+epicel+ulcer,data=Melanoma)
##' f2 <- CSC(Hist(time,status)~age+sex+logthick+epicel+ulcer,data=Melanoma)
##' x <- Score(list(model1=f1,model2=f2),Hist(time,status)~1,data=Melanoma,
##'             cause=1,times=5*365.25,plots="roc",metrics="auc")
##' plotROC(x)
#' @export
myplotROC <- function(x,
                      models,
                      times,
                      xlab="1-Specificity",
                      ylab="Sensitivity",
                      col,
                      lwd,
                      lty=1,
                      cex=1,
                      pch=1,
                      legend=TRUE,
                      auc.in.legend=TRUE,
                      brier.in.legend=FALSE,
                      ipa.in.legend=FALSE,
                      add=FALSE,
                      ...){
  if (is.null(x$ROC))
    stop("Object has no information for ROC curves.\nYou should call the function \"riskRegression::Score\" with plots=\"ROC\".")
  model=FPR=TPR=NULL
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  pframe <- x$ROC$plotframe
  if (x$response.type!="binary"){
    if (missing(times)){
      tp <- max(pframe[["times"]])
      if (length(unique(pframe$times))>1)
        warning("Time point not specified, use max of the available times: ",tp)
    } else{ ## can only do one time point
      tp <- times[[1]]
      if (!(tp%in%unique(pframe$times)))
        stop(paste0("Requested time ",times[[1]]," is not in object"))
    }
    pframe <- pframe[times==tp]
  }else tp <- NULL
  if (!missing(models)){
    if (!all(models%in% pframe[,unique(model)]))
      stop(paste0("Cannot identify model names.\nRequested models: ",paste(models,collapse=", "),"\n",
                  "Available models: ",paste(pframe[,unique(model)],collapse=", ")))
    pframe <- pframe[model%in%models]
    ## user's order
    pframe[,model:=factor(model,levels=models)]
  }else{
    if (length(x$null.model)>0){
      pframe <- pframe[model!=x$null.model]
    }
    models=unique(pframe$model)
  }
  setkey(pframe,model)
  mm <- unique(pframe$model)
  lenmm <- length(mm)
  if(missing(col)) col <- rep(cbbPalette,length.out=lenmm)
  names(col) <- mm
  if(missing(lwd)) lwd <- 3
  lwd <- rep(lwd,length.out=lenmm)
  names(lwd) <- mm
  pch <- rep(pch,length.out=lenmm)
  names(pch) <- mm
  if(missing(lwd)) lty <- 1
  lty <- rep(lty,length.out=lenmm)
  names(lty) <- mm
  lines.DefaultArgs <- list(pch=pch,cex=cex,lwd=lwd,type="l",col=col,lty=lty)
  # {{{ legend
  if (is.character(legend[[1]])|| legend[[1]]==TRUE){
    legend.data <- getLegendData2(object=x,
                                  models=models,
                                  times=tp,
                                  auc.in.legend=auc.in.legend,
                                  brier.in.legend=brier.in.legend,
                                  ipa.in.legend=ipa.in.legend,
                                  drop.null.model=TRUE)
    # legend.data
    if (is.character(legend))
      legend.text <- legend
    else
      legend.text <- unlist(legend.data[,1])
    nrows.legend <- NROW(legend.data)
    if (nrows.legend==1){
      legend.lwd <- NA
    }else{
      legend.lwd <- lwd
    }
    legend.DefaultArgs <- list(legend=legend.text,
                               lwd=legend.lwd,
                               col=col,
                               ncol=1,
                               lty=lty,
                               cex=cex,
                               bty="n",
                               y.intersp=1,
                               x="bottom",
                               title="")
    if (NCOL(legend.data)>1){
      addtable2plot.DefaultArgs <- list(yjust=1.18, cex=cex, table=legend.data[,-1,drop=FALSE])
    }else{
      addtable2plot.DefaultArgs <- NULL
    }
  }else{
    legend.DefaultArgs <- NULL
    addtable2plot.DefaultArgs <- NULL
  }
  # }}}
  plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = c(0,1),xlim = c(0,1),ylab=ylab,xlab=xlab)
  axis1.DefaultArgs <- list(side=1,las=1,at=seq(0,1,.25))
  axis2.DefaultArgs <- list(side=2,las=1,at=seq(0,1,.25))
  control <- prodlim::SmartControl(call= list(...),
                                   keys=c("plot","lines","legend","addtable2plot","axis1","axis2"),
                                   ignore=NULL,
                                   ignore.case=TRUE,
                                   defaults=list("plot"=plot.DefaultArgs,
                                                 "lines"=lines.DefaultArgs,
                                                 "legend"=legend.DefaultArgs,
                                                 "addtable2plot"=addtable2plot.DefaultArgs,
                                                 "axis1"=axis1.DefaultArgs,
                                                 "axis2"=axis2.DefaultArgs),
                                   forced=list("plot"=list(axes=FALSE),
                                               "axis1"=list(side=1),"axis2"=list(side=2)),
                                   verbose=TRUE)
  if (add==0L) {do.call("plot",control$plot)
    if (is.null(control$axis1$labels)){
      control$axis1$labels <- paste(100*control$axis1$at,"%")
    }
    if (is.null(control$axis2$labels)){
      control$axis2$labels <- paste(100*control$axis2$at,"%")
    }
    do.call("axis",control$axis1)
    do.call("axis",control$axis2)
  }
  if (is.character(legend[[1]]) || legend[[1]]==TRUE){
    legend.coords <- do.call("legend",control$legend)
    if (!is.null(addtable2plot.DefaultArgs)){
      if (is.null(control$addtable2plot[["x"]]))
        control$addtable2plot[["x"]] <- legend.coords$rect$left+legend.coords$rect$w
      ## strange error: $y does not work as yjust is matched, thus use [["y"]]
      if (is.null(control$addtable2plot[["y"]]))
        control$addtable2plot[["y"]] <- legend.coords$rect$top-legend.coords$rect$h
      do.call(plotrix::addtable2plot,control$addtable2plot)
    }
  }
  abline(a=0,b=1,col="gray77",lwd=3)
  pframe[,{thisline <- control$line;
  thisline$col=thisline$col[[as.character(model[1])]];
  thisline$lwd=thisline$lwd[[as.character(model[1])]];
  thisline$lty=thisline$lty[[as.character(model[1])]];
  thisline$pch=thisline$pch[[as.character(model[1])]];
  thisline$x=c(0,FPR,1);
  thisline$y=c(0,TPR,1);
  do.call("lines",thisline)},by=model]
  invisible(x)
}

#----------------------------------------------------------------------
### plotROC.R ends here



myTileHeatPlot <- function(DF, xVar, yVar,xLabel, yLabel, labelVar=NULL, colorVar,facetVar =NULL){
  # DF <- res_w.treat# %>% dplyr::filter(Model=="TCR Richness") %>% droplevels()
  # xVar <-"times"
  # yVar <-"Variable"
  # colorVar <-"IPA.drop"
  # labelVar <- "IPA.drop"
  # xLabel = "Months"
  # yLabel = "Variable"
  
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  colorVar.un <- as.name(colorVar)
  p <- ggplot(DF, aes(x=!!xVar.un, y=!!yVar.un, fill=!!colorVar.un)) + viridis::scale_fill_viridis(option = "D",
                                                                                                   direction=1,
                                                                                                   alpha=0.5,
                                                                                                   label = function(x) sprintf("%.2f", x)) + #, + 
    # limits=c(0.44,.89)
    
    
    geom_tile(alpha=0.7,color="gray")
  if(!is.null(labelVar)){
    labelVar.un <- as.name(labelVar)
    
    p <- p +geom_text(aes(label = round(!!labelVar.un,2)), color="black", size=4, fontface="bold")
    
  }
  
  if(!is.null(facetVar)){
    p <- p + facet_wrap(facetVar,nrow=5,scales='free')
    
  }
  
  p <- p + theme_bw(base_family = "Palatino", base_size = 12)
  p <- p + theme(legend.position = "bottom",
                 panel.border = element_blank(),#element_rect(colour = "black", size=0.2),
                 # remove the horizontal grid lines
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.text.x = element_text(size=12,face="bold"),
                 axis.text.y = element_text(size = 12,face="bold"),#,angle = 90,vjust=0.9
                 axis.ticks.x = element_blank(),
                 axis.ticks.y = element_blank(),
                 # axis.title.x = element_text(size=14, face="bold"),
                 # axis.title.y = element_text(size=14, face="bold"),
                 legend.key = element_rect(size = 2),
                 legend.key.height=unit(0.6,"line"),
                 legend.key.width=unit(2,"line"),
                 legend.text = element_text(size = 11),
                 # legend.background = element_rect(colour = "black"),
                 legend.direction='horizontal',legend.box='horizontal',
                 legend.title = element_text(colour="black", size=12, face="plain"))
  p<-p + labs(x = "", y = "")
  p
}


myTileHeatPlot.2 <- function(DF, xVar, yVar,xLabel, yLabel, labelVar=NULL, colorVar,colorCont=FALSE,facetVar =NULL, removeLegend=TRUE){
  # DF <- res_w.treat# %>% dplyr::filter(Model=="TCR Richness") %>% droplevels()
  # xVar <-"times"
  # yVar <-"Variable"
  # colorVar <-"IPA.drop"
  # labelVar <- "IPA.drop"
  # xLabel = "Months"
  # yLabel = "Variable"
  
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  colorVar.un <- as.name(colorVar)
  
  if(isTRUE(colorCont)){
    p <- ggplot(DF, aes(x=!!xVar.un, y=!!yVar.un, fill=!!colorVar.un)) + viridis::scale_fill_viridis(option = "D",
                                  direction=1,
                                  alpha=0.5,
                                  label = function(x) sprintf("%.2f", x)) + geom_tile(alpha=0.7,color="gray")
  }else{
    p <- ggplot(DF, aes(x=!!xVar.un, y=!!yVar.un, fill=!!colorVar.un)) + scale_fill_manual(values=c("white","black")) + #, + #
      geom_tile(alpha=0.7,color="gray")
    
  }
  if(!is.null(labelVar)){
    labelVar.un <- as.name(labelVar)
    
    p <- p +geom_text(aes(label = round(!!labelVar.un,2)), color="black", size=4, fontface="bold")
    
  }
  
  if(!is.null(facetVar)){
    p <- p + facet_wrap(facetVar,nrow=5,scales='free')
    
  }
  
  p <-p+ scale_x_discrete(position = "top") 
  p <- p + ylab("Boot models")
  p <- p + theme_bw(base_family = "Palatino", base_size = 12)
  p <- p + theme(legend.position = "bottom",
                 panel.border = element_blank(),#element_rect(colour = "black", size=0.2),
                 # remove the horizontal grid lines
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.text.x = element_text(size=12,face="bold"),
                 axis.text.y = element_blank(),#element_text(size = 12,face="bold"),#,angle = 90,vjust=0.9
                 axis.ticks.x = element_blank(),
                 axis.ticks.y = element_blank(),
                 # axis.title.x = element_text(size=14, face="bold"),
                 # axis.title.y = element_text(size=14, face="bold"),
                 legend.key = element_rect(size = 2),
                 legend.key.height=unit(0.5,"line"),
                 legend.key.width=unit(1.5,"line"),
                 legend.text = element_text(size = 11),
                 legend.background = element_rect(colour = "black"),
                 legend.direction='horizontal',legend.box='horizontal',
                 legend.title = element_blank())#element_text(colour="black", size=11, face="plain")
  p<-p + labs(x = "", y = "")
  
  if(isTRUE(removeLegend))
    p <- p + theme(legend.position = "none")
  p
}



myCalPlot <-function (x, xlab, ylab, subtitles = TRUE, conf.int = TRUE, cex.subtitles = 0.75, 
                      riskdist = TRUE, add = FALSE, scat1d.opts = list(nhistSpike = 200), 
                      par.corrected = NULL, ...) 
{
  # x=cal
  # subtitles = TRUE
  # conf.int = TRUE
  # cex.subtitles = 0.75
  # riskdist = TRUE
  # add = FALSE
  # scat1d.opts = list(nhistSpike = 200)
  # par.corrected = NULL
  # xlab = NULL
  # ylab = NULL
  ##########
  at <- attributes(x)
  u <- at$u
  units <- at$units
  if (length(par.corrected) && !is.list(par.corrected)) 
    stop("par.corrected must be a list")
  z <- list(col = "black", lty = 1, lwd = 1.6, pch = 4)
  if (!length(par.corrected)) 
    par.corrected <- z
  else for (n in setdiff(names(z), names(par.corrected))) par.corrected[[n]] <- z[[n]]
  predicted <- at$predicted
  if ("KM" %in% colnames(x)) {
    type <- "stratified"
    pred <- x[, "mean.predicted"]
    cal <- x[, "KM"]
    cal.corrected <- x[, "KM.corrected"]
    se <- x[, "std.err"]
  }else {
    type <- "smooth"
    pred <- x[, "pred"]
    cal <- x[, "calibrated"]
    cal.corrected <- x[, "calibrated.corrected"]
    se <- NULL
  }
  un <- if (u == 1) 
    paste(units, "s", sep = "")
  else units
  if (missing(xlab)) 
    xlab <- paste("Predicted ", format(u), units, "Survival")
  if (missing(ylab)) 
    ylab <- paste("Fraction Surviving ", format(u), 
                  " ", un, sep = "")
  if (length(se) && conf.int) {
    ciupper <- function(surv, d) ifelse(surv == 0, 0, pmin(1, 
                                                           surv * exp(d)))
    cilower <- function(surv, d) ifelse(surv == 0, 0, surv * 
                                          exp(-d))
    errbar(pred, cal, cilower(cal, 1.959964 * se), ciupper(cal, 
                                                           1.959964 * se), xlab = xlab, ylab = ylab, type = "b", 
           add = add, ...)
  }
  else if (add) 
    lines(pred, cal, type = if (type == "smooth") 
      "l"
      else "b")
  else plot(pred, cal, xlab = xlab, ylab = ylab, type = if (type == 
                                                            "smooth") 
    "o"
    else "b",cex=0.3, col="gray",...)
  err <- NULL
  if (riskdist && length(predicted)) {
    do.call("scat1d", c(list(x = predicted), scat1d.opts))
    if (type == "smooth") {
      s <- !is.na(pred + cal.corrected)
      err <- predicted - approxExtrap(pred[s], cal.corrected[s], 
                                      xout = predicted, ties = mean)$y
    }
  }
  if (subtitles && !add) {
    if (type == "smooth") {
      Col <- par.corrected$col
      substring(Col, 1, 1) <- toupper(substring(Col, 1, 
                                                1))
      title(sub = sprintf("Gray: observed  Orange: ideal\n%s : optimism corrected", 
                          Col), adj = 0, cex.sub = cex.subtitles)
      w <- if (length(err)) 
        paste("B=", at$B, " based on ", at$what, 
              "\nMean |error|=", round(mean(abs(err)), 
                                       3), "  0.9 Quantile=", round(quantile(abs(err), 
                                                                             0.9, na.rm = TRUE), 3), sep = "")
      else paste("B=", at$B, "\nBased on ", 
                 at$what, sep = "")
      title(sub = w, adj = 1, cex.sub = cex.subtitles)
    }
    else {
      title(sub = paste("n=", at$n, " d=", 
                        at$d, " p=", at$p, ", ", at$m, " subjects per group\nGray: ideal", 
                        sep = ""), adj = 0, cex.sub = cex.subtitles)
      title(sub = paste("X - resampling optimism added, B=", 
                        at$B, "\nBased on ", at$what, sep = ""), 
            adj = 1, cex.sub = cex.subtitles)
    }
  }
  # abline(0, 1, col = gray(0.9))
  abline(0, 1, lwd = 1.8, lty = 2, col = "#FFA500")
  if (type == "stratified") 
    points(pred, cal.corrected, pch = par.corrected$pch, 
           col = par.corrected$col)
  else lines(pred, cal.corrected, col = par.corrected$col, 
             lty = par.corrected$lty, lwd = par.corrected$lwd)
  invisible()
}


calibrationGGplot_TEST <- function(calObj,optCorColor){
  # calObj =res.cal.strat.nointer_RED$mice.bw$`36`$cal
  # optCorColor = "orange"
  ### START HERE
  
  at    <- attributes(calObj)
  u     <- at$u
  units <- at$units
  un <- if(u==1) paste(units, 's', sep='') else units
  xlab <- paste("Predicted ",format(u),units,"Survival")
  ylab <- paste("Fraction Surviving ",format(u)," ",un,
                sep="")
  
  predicted <- at$predicted #  for histogram
  pred <- calObj[,'pred']
  cal  <- calObj[,'calibrated']
  cal.corrected <- calObj[,'calibrated.corrected']
  
  x.max=round(max(cal.corrected,cal,pred),1)
  y.max=round(max(cal.corrected,cal,pred),1)
  print(c(x.max, y.max))
  # if(x.max >1){
  #   x.max= 1
  # }
  # 
  # if(y.max >1){
  #   y.max=1
  # }
  # se <- NULL
  # s <- !is.na(pred + cal.corrected)
  # err <- predicted -
  #   approxExtrap(pred[s], cal.corrected[s], xout=predicted, ties=mean)$y
  # 
  # sub.1 =paste('B=', at$B, ' based on ', at$what,
  #              '\nMean |error|=', round(mean(abs(err)), 3),
  #              '  0.9 Quantile=',
  #              round(quantile(abs(err), .9, na.rm=TRUE), 3),
  #              sep='')
  # sub.2 = sprintf('Black: optimism corrected\nGray: observed %s : ideal',
  #                 optCorColor)
  # 
  # # Create dataframe for plotting
  # df = data.frame(pred=pred, cal=cal, cal.corrected=cal.corrected)
  # predicted.df <- as.data.frame(predicted)
  # predicted.df$y <-y.max
  # h = predicted.df %>%
  #   mutate(breaks = cut(predicted, breaks=seq(0,x.max,0.005), labels=seq(0.0005,x.max,0.005),
  #                       include.lowest=TRUE),
  #          breaks = as.numeric(as.character(breaks))) %>%
  #   dplyr::group_by(y, breaks) %>% 
  #   dplyr::summarise(n = n()) %>%
  #   mutate(pct = x.max - n/sum(n)) 
  # h$breaks <- h$breaks %>% as.character() %>% as.numeric()
  # # Start plotting
  # p <- ggplot(data=df, aes(x=pred))
  # # Add observed line
  # p <- p + geom_line(aes(y = cal), color = "gray",size=.75) 
  # # Add optimism corrected line
  # p <- p +  geom_line(aes(y = cal.corrected), color="black",size=.75) 
  # # Add ideal line
  # p <- p + geom_abline(intercept = 0, slope = 1, color="#FFA500", linetype="dashed",size=.9) 
  # p <- p + ylab(ylab)
  # p <- p + xlab(xlab)
  # p <- p + labs(caption = c(sub.1, sub.2))
  # # Add risk density as rug plot on top
  # p <- p + geom_rug(data=predicted.df,aes(x=predicted), sides = "t",length = unit(0.01, "npc"))
  # # p <- p + ggMarginal(data = predicted.df, x=predicted, y=1,type="histogram")
  # # p + geom_segment(dat=predicted.df, aes(x=predicted, xend=predicted, y=1, yend=1.02), size=0.2, colour="grey30")
  # p <- p + geom_segment(data=h, size=.6, show.legend=FALSE,
  #                       aes(x=breaks, xend=breaks, y=y, yend=pct, colour="black"), colour="black")
  # p <- p + xlim(0,x.max) + ylim(0,y.max)
  # p <- p + theme_bw(base_family = "Palatino", base_size = 13) #+ labs_pubr(base_size = 14, base_family = "Palatino Linotype")
  # 
  # p <- p + theme(legend.position = "right",
  #                legend.box = "horizontal",
  #                legend.direction = "vertical",
  #                legend.title = element_blank(),
  #                legend.text = element_text(size=12),
  #                legend.key=element_blank(),
  #                axis.text = element_text(size=11),
  #                axis.title = element_text(size=13),
  #                panel.grid.major = element_blank(),
  #                panel.grid.minor = element_blank(),
  #                panel.background = element_blank(),
  #                plot.caption = element_text(hjust=c(1, 0)))
  # 
  # 
  # p
}


calibrationGGplot <- function(calObj,optCorColor){
  # calObj =res.cal.strat.nointer_RED$mice.bw$`36`$cal
  # optCorColor = "orange"
  ### START HERE
  
  at    <- attributes(calObj)
  u     <- at$u
  units <- at$units
  un <- if(u==1) paste(units, 's', sep='') else units
  xlab <- paste("Predicted ",format(u),units,"Survival")
  ylab <- paste("Fraction Surviving ",format(u)," ",un,
                sep="")
  
  predicted <- at$predicted #  for histogram
  pred <- calObj[,'pred']
  cal  <- calObj[,'calibrated']
  cal.corrected <- calObj[,'calibrated.corrected']
  
  x.max=round(max(cal.corrected,cal,pred),1)
  y.max=round(max(cal.corrected,cal,pred),1)
  if(x.max >1){
    x.max= 1
  }
  
  if(y.max >1){
    y.max=1
  }
  se <- NULL
  s <- !is.na(pred + cal.corrected)
  err <- predicted -
    approxExtrap(pred[s], cal.corrected[s], xout=predicted, ties=mean)$y
  
  sub.1 =paste('B=', at$B, ' based on ', at$what,
               '\nMean |error|=', round(mean(abs(err)), 3),
               '  0.9 Quantile=',
               round(quantile(abs(err), .9, na.rm=TRUE), 3),
               sep='')
  sub.2 = sprintf('Black: optimism corrected\nGray: observed %s : ideal',
                  optCorColor)
  
  # Create dataframe for plotting
  df = data.frame(pred=pred, cal=cal, cal.corrected=cal.corrected)
  predicted.df <- as.data.frame(predicted)
  predicted.df$y <-y.max
  h = predicted.df %>%
    mutate(breaks = cut(predicted, breaks=seq(0,x.max,0.005), labels=seq(0.0005,x.max,0.005),
                        include.lowest=TRUE),
           breaks = as.numeric(as.character(breaks))) %>%
    dplyr::group_by(y, breaks) %>% 
    dplyr::summarise(n = n()) %>%
    mutate(pct = x.max - n/sum(n)) 
  h$breaks <- h$breaks %>% as.character() %>% as.numeric()
  # Start plotting
  p <- ggplot(data=df, aes(x=pred))
  # Add observed line
  p <- p + geom_line(aes(y = cal), color = "gray",size=.75) 
  # Add optimism corrected line
  p <- p +  geom_line(aes(y = cal.corrected), color="black",size=.75) 
  # Add ideal line
  p <- p + geom_abline(intercept = 0, slope = 1, color="#FFA500", linetype="dashed",size=.9) 
  p <- p + ylab(ylab)
  p <- p + xlab(xlab)
  p <- p + labs(caption = c(sub.1, sub.2))
  # Add risk density as rug plot on top
  p <- p + geom_rug(data=predicted.df,aes(x=predicted), sides = "t",length = unit(0.01, "npc"))
  # p <- p + ggMarginal(data = predicted.df, x=predicted, y=1,type="histogram")
  # p + geom_segment(dat=predicted.df, aes(x=predicted, xend=predicted, y=1, yend=1.02), size=0.2, colour="grey30")
  p <- p + geom_segment(data=h, size=.6, show.legend=FALSE,
                 aes(x=breaks, xend=breaks, y=y, yend=pct, colour="black"), colour="black")
  p <- p + xlim(0,x.max) + ylim(0,y.max)
  p <- p + theme_bw(base_family = "Palatino", base_size = 13) #+ labs_pubr(base_size = 14, base_family = "Palatino Linotype")
  
  p <- p + theme(legend.position = "right",
          legend.box = "horizontal",
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(size=12),
          legend.key=element_blank(),
          axis.text = element_text(size=11),
          axis.title = element_text(size=13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.caption = element_text(hjust=c(1, 0)))
  
  
  p
}


get_ggcoefsData_pooled <- function(fitObj,pooledObj, significance=.05){
  # fitObj <- fit.pooled2
  # pooledObj <- pool_coxr.cox.nointer %>%
  #           pluck("RR_model_final") %>%
  #           magrittr::extract2(1)

  significance=0.05
  ## START
  data.res <- ggcoef_model(fitObj,
                           exponentiate = TRUE, # Equivalent to HR values
                           add_reference_rows = FALSE, # this fits the output of pooled model coefs
                           # colour = NULL, 
                           # stripped_rows = FALSE, 
                           return_data = T) #%>% as.data.frame()
  
  # dim(data.res)
  # dim(pooledObj)
  # keep.init <- data.res
  # Replace pooled HR, confint,. p values, statistics stderror in data.res from pooled
  # data.res$estimate <- pooledObj$HR # no need, it's the same
  data.res$std.error <- pooledObj$std.error
  data.res$statistic <- pooledObj$statistic
  data.res$p.value <- pooledObj$p.value
  data.res$conf.low <- pooledObj$lower.EXP
  data.res$conf.high <- pooledObj$upper.EXP
  # Now fix significance columns in data.res.they need to be updated
  significance_labels <- paste(c("p \u2264", "p >"), significance)
  data.res$significance <- factor(
    !is.na(data.res$p.value) & data.res$p.value <= significance,
    levels = c(TRUE, FALSE),
    labels = significance_labels)
  
  data.res$signif_stars <- signif_stars(data.res$p.value, point = NULL)
  data.res$p_value_label <- ifelse(
    is.na(data.res$p.value),
    "",
    scales::pvalue(data.res$p.value, add_p = TRUE)
  )
  
  # Fix add_to_label: p value found plus the stars of significance 
  data.res$add_to_label <- paste0(data.res$p_value_label, data.res$signif_stars)
  # Fix label (only levels)
  mylabels<- sapply(1:nrow(data.res), function(x) ifelse(data.res$term[x]==data.res$variable[x],data.res$term[x],str_remove(data.res$term[x], data.res$variable[x])), simplify=FALSE) %>% ldply() %>% pull(V1)
  data.res$label <-mylabels #%>% as.factor()
  
  
  data.res$label <- forcats::fct_inorder(factor(paste0(data.res$label, 
                                                       ifelse(data.res$add_to_label == "", "", paste0(" (", 
                                                                                                      data.res$add_to_label, ")")))))
  
  data.res$label_light <- forcats::fct_inorder(factor(paste0(data.res$label_light, 
                                                             ifelse(data.res$add_to_label == "", "", paste0(" (", 
                                                                                                            data.res$add_to_label, ")")))))
  
  
  data.res$label_light <- dplyr::if_else(
    as.character(data.res$label) == as.character(data.res$var_label) &
      ((!grepl("^nmatrix", data.res$var_class)) | is.na(data.res$var_class)),
    "",
    as.character(data.res$label)
  ) %>%
    my.in_order()
  data.res
}



final_model_MI_poolCoefs <- function(origMIstacked, impSel, impvar, time, status, features,pooledObj){
  
  # origMIstacked = orig.imp.final.cox
  # impSel = 1
  # impvar = "Impnr"
  # time = surv_time
  # status = surv_status
  # features = feats.trans.strata
  # pooledObj = pool_coxr.cox.nointer
  
  ######
  # Here I construct a fit oject but replace with pooled estimates from MI pooling process
  
  impvar.un <- as.name(impvar)
  somData<-origMIstacked %>% dplyr::filter(!!impvar.un==impSel)
  
  fit.pooled<- coxph(modelFormula(time,status,features ), data = somData)
  # fit.pooled <- fit.strata.nointer.nobw
  # Add final pooled coefs to the model
  coef_names <-pooledObj %>% 
    pluck("RR_model_final") %>% 
    magrittr::extract2(1) %>%
    pull(term) %>% as.character()
  
  fit.pooled$coefficients <-pooledObj %>% 
    pluck("RR_model_final") %>%
    magrittr::extract2(1) %>%
    pull(estimate)  %>% as.vector() %>%
    setNames(coef_names)
  
  fit.pooled$linear.predictors  <- pooledObj$pool_lp_final %>%  magrittr::extract2(1) %>% setNames(NULL)
  
  fit.pooled$means <- pooledObj$pool_covCenter_final %>%
    magrittr::extract2(1)
  #fit.pooled
  coefs.data<- get_ggcoefsData_pooled(fitObj=fit.pooled,pooledObj=pooledObj %>% 
                                        pluck("RR_model_final") %>%
                                        magrittr::extract2(1), significance=.05)
  model.pooled <- list(fit=fit.pooled,
                       coefs=coefs.data)
  model.pooled
  
}


final_model_MI_pooled.Fit <- function(origMIstacked, impSel, impvar, time, status, features,pooledObj){
  
  # origMIstacked = orig.imp.final.cox
  # impSel = 1
  # impvar = "Impnr"
  # time = surv_time
  # status = surv_status
  # features = feats.trans.strata
  # pooledObj = pool_coxr.cox.nointer

  ######
  # Here I construct a fit oject but replace with pooled estimates from MI pooling process
  
  impvar.un <- as.name(impvar)
  somData<-origMIstacked %>% dplyr::filter(!!impvar.un==impSel)
  
  fit.pooled<- coxph(modelFormula(time,status,features ), data = somData)
  # fit.pooled <- fit.strata.nointer.nobw
  # Add final pooled coefs to the model
  coef_names <-pooledObj %>% 
    pluck("RR_model_final") %>% 
    magrittr::extract2(1) %>%
    pull(term) %>% as.character()
  
  fit.pooled$coefficients <-pooledObj %>% 
    pluck("RR_model_final") %>%
    magrittr::extract2(1) %>%
    pull(estimate)  %>% as.vector() %>%
    setNames(coef_names)
  
  fit.pooled$linear.predictors  <- pooledObj$pool_lp_final %>%  magrittr::extract2(1) %>% setNames(NULL)
  
  fit.pooled$means <- pooledObj$pool_covCenter_final %>%
    magrittr::extract2(1)
  # #fit.pooled
  # coefs.data<- get_ggcoefsData_pooled(fitObj=fit.pooled,pooledObj=pooledObj %>% 
  #                                       pluck("RR_model_final") %>%
  #                                       magrittr::extract2(1), significance=.05)
  # model.pooled <- list(fit=fit.pooled,
  #                      coefs=coefs.data)
  # model.pooled
  return(fit.pooled)
}










#####################





### For model comparing
compare_dataCoefs<- function(data){
  # data <- coefs.models
  # type = "dodged"#c("dodged", "faceted")
  # exponentiate = TRUE
  data <- dplyr::bind_rows(data, .id = "model")
  coefficients_label <- attr(data, "coefficients_label")
  data$model <- forcats::fct_inorder(data$model)
  # data <- data %>% broom.helpers::tidy_select_variables(
  #                 model = fit.strata.nointer.nobw.COMPL) %>% broom.helpers::tidy_detach_model()
  # data <- data %>% tidyr::complete(.data$model, tidyr::nesting(!!sym("variable"), 
  #                                                              !!sym("var_label"), !!sym("var_class"), !!sym("var_type"), 
  #                                                              !!sym("contrasts"), !!sym("reference_row"), 
  #                                                              !!sym("label"), !!sym("label_light")))
  attr(data, "coefficients_label") <- coefficients_label
  data
}



mybarplot <- function(df,xVar,yVar,xLabel,yLabel,facetVar,bootFreq=0.5){
  # Basic barplot
  p<-ggplot(data=df, aes_string(x=xVar, y=yVar)) +
    geom_bar(stat="identity")
  if(!is.null(facetVar)){
    p <- p + facet_wrap(facetVar,ncol=2,scales='free')
    
  }
  p <-p + geom_hline(yintercept = bootFreq, colour="indian red", linetype="dashed")
  p <- p + ylab(yLabel)
  p <- p + xlab(xLabel)
  p <- p + ylim(0,100)
  # Horizontal bar plot
  p <-p + coord_flip()
  #p
  p <- p + theme_ipsum(base_family = "Palatino", base_size = 10)
  p <- p + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                 axis.title.x = element_text(size=10, face="bold"),
                 axis.text.y = element_text( hjust = 1, vjust = 0,size = 10),
                 axis.title.y = element_text(size=11, face="bold"),
                 legend.title = element_text(size=11,color="black"),
                 legend.text = element_text(size = 10),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.position="top") 
  p
}


#'Graphical Test of Proportional Hazards with ggplot2
#'@description Displays a graph of the scaled Schoenfeld residuals, along with a
#'  smooth curve using \pkg{ggplot2}. Wrapper around \link{plot.cox.zph}.
#'@param fit an object of class \link{cox.zph}.
#'@param resid	a logical value, if TRUE the residuals are included on the plot,
#'  as well as the smooth fit.
#'@param se a logical value, if TRUE, confidence bands at two standard errors
#'  will be added.
#'@param df	the degrees of freedom for the fitted natural spline, df=2 leads to
#'  a linear fit.
#'@param nsmo	number of points used to plot the fitted spline.
#'@param var the set of variables for which plots are desired. By default, plots
#'  are produced in turn for each variable of a model.
#'@param point.col,point.size,point.shape,point.alpha color, size, shape and visibility to be used for points.
#'@param caption the caption of the final \link{grob} (\code{bottom} in \link{arrangeGrob})
#'@param ggtheme function, ggplot2 theme name.
#'  Allowed values include ggplot2 official themes: see \code{\link[ggplot2]{theme}}.
#'@param ... further arguments passed to either the print() function or to the \code{\link[ggpubr]{ggpar}} function for customizing the plot (see Details section).
#'@details \strong{Customizing the plots}: The plot can be easily
#'  customized using additional arguments to be passed to the function ggpar().
#'  Read ?ggpubr::ggpar. These arguments include
#'  \emph{font.main,font.submain,font.caption,font.x,font.y,font.tickslab,font.legend}:
#'  a vector of length 3 indicating respectively the size (e.g.: 14), the style
#'  (e.g.: "plain", "bold", "italic", "bold.italic") and the color (e.g.: "red")
#'  of main title, subtitle, caption, xlab and ylab and axis tick labels,
#'  respectively. For example \emph{font.x = c(14, "bold", "red")}.  Use font.x
#'  = 14, to change only font size; or use font.x = "bold", to change only font
#'  face.
#'@return Returns an object of class \code{ggcoxzph} which is a list of ggplots.
#'
#'@author Marcin Kosinski , \email{m.p.kosinski@@gmail.com}
#'
#'@examples
#'
#' library(survival)
#' fit <- coxph(Surv(futime, fustat) ~ age + ecog.ps + rx, data=ovarian)
#' cox.zph.fit <- cox.zph(fit)
#' # plot all variables
#' ggcoxzph(cox.zph.fit)
#' # plot all variables in specified order
#' ggcoxzph(cox.zph.fit, var = c("ecog.ps", "rx", "age"), font.main = 12)
#' # plot specified variables in specified order
#' ggcoxzph(cox.zph.fit, var = c("ecog.ps", "rx"), font.main = 12, caption = "Caption goes here")
#'
#'@describeIn ggcoxzph Graphical Test of Proportional Hazards using ggplot2.
#'@export
ggcoxzphFixed <- function (fit, resid = TRUE, se = TRUE, df = 4, nsmo = 40, var,
                           point.col = "red", point.size = 1, point.shape = 19, point.alpha = 1,
                           caption = NULL,
                           ggtheme = theme_survminer(), ...){
  
  x <- fit
  if(!methods::is(x, "cox.zph"))
    stop("Can't handle an object of class ", class(x))
  
  xx <- x$x
  yy <- x$y
  d <- nrow(yy)
  df <- max(df)
  nvar <- ncol(yy)
  pred.x <- seq(from = min(xx), to = max(xx), length = nsmo)
  temp <- c(pred.x, xx)
  lmat <- splines::ns(temp, df = df, intercept = TRUE)
  pmat <- lmat[1:nsmo, ]
  xmat <- lmat[-(1:nsmo), ]
  qmat <- qr(xmat)
  if (qmat$rank < df)
    stop("Spline fit is singular, try a smaller degrees of freedom")
  if (se) {
    bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
    xtx <- bk %*% t(bk)
    # seval <- d * ((pmat %*% xtx) * pmat) %*% rep(1, df)
    seval <- ((pmat %*% xtx) * pmat) %*% rep(1, df)
    
  }
  ylab <- paste("Beta(t) for", dimnames(yy)[[2]])
  if (missing(var))
    var <- 1:nvar
  else {
    if (is.character(var))
      var <- match(var, dimnames(yy)[[2]])
    if (any(is.na(var)) || max(var) > nvar || min(var) <
        1)
      stop("Invalid variable requested")
  }
  if (x$transform == "log") {
    xx <- exp(xx)
    pred.x <- exp(pred.x)
  }
  else if (x$transform != "identity") {
    xtime <- as.numeric(dimnames(yy)[[1]])
    indx <- !duplicated(xx)
    apr1 <- approx(xx[indx], xtime[indx], seq(min(xx), max(xx),
                                              length = 17)[2 * (1:8)])
    temp <- signif(apr1$y, 2)
    apr2 <- approx(xtime[indx], xx[indx], temp)
    xaxisval <- apr2$y
    xaxislab <- rep("", 8)
    for (i in 1:8) xaxislab[i] <- format(temp[i])
  }
  plots <- list()
  lapply(var, function(i) {
    invisible(round(x$table[i, 3],4) -> pval)
    ggplot() + labs(title = paste0('Schoenfeld Individual Test p: ', pval)) + ggtheme -> gplot
    y <- yy[, i]
    yhat <- as.vector(pmat %*% qr.coef(qmat, y))
    if (resid)
      yr <- range(yhat, y)
    else yr <- range(yhat)
    if (se) {
      temp <- as.vector(2 * sqrt(x$var[i, i] * seval))
      yup <- yhat + temp
      ylow <- yhat - temp
      yr <- range(yr, yup, ylow)
    }
    if (x$transform == "identity") {
      gplot + geom_line(aes(x=pred.x, y=yhat)) +
        xlab("Time") +
        ylab(ylab[i]) +
        ylim(yr) -> gplot
    } else if (x$transform == "log") {
      gplot + geom_line(aes(x=log(pred.x), y=yhat)) +
        xlab("Time") +
        ylab(ylab[i]) +
        ylim(yr)  -> gplot
    } else {
      gplot + geom_line(aes(x=pred.x, y=yhat)) +
        xlab("Time") +
        ylab(ylab[i]) +
        scale_x_continuous(breaks = xaxisval,
                           labels = xaxislab) +
        ylim(yr)-> gplot
    }
    
    if (resid)
      gplot <- gplot + geom_point(aes(x = xx, y =y),
                                  col = point.col, shape = point.shape, size = point.size, alpha = point.alpha)
    
    if (se) {
      gplot <- gplot + geom_line(aes(x=pred.x, y=yup), lty = "dashed") +
        geom_line(aes( x = pred.x, y = ylow), lty = "dashed")
    }
    
    ggpubr::ggpar(gplot, ...)
    
    
  }) -> plots
  names(plots) <- var
  class(plots) <- c("ggcoxzph", "ggsurv", "list")
  
  if("GLOBAL" %in% rownames(x$table)) # case of multivariate Cox
    global_p <- x$table["GLOBAL", 3]
  else global_p <- NULL # Univariate Cox
  attr(plots, "global_pval") <- global_p
  attr(plots, "caption") <- caption
  plots
  
}

#' @param x an object of class ggcoxzph
#' @param newpage open a new page. See \code{\link{grid.arrange}}.
#' @method print ggcoxzph
#' @rdname ggcoxzph
#' @export
print.ggcoxzph <- function(x, ..., newpage = TRUE){
  if(!inherits(x, "ggcoxzph"))
    stop("An object of class ggcoxzph is required.")
  plots <- x
  pval <- attr(x, "global_pval")
  grobs <- widths <- list()
  for (i in 1:length(plots)) {
    grobs[[i]] <- ggplotGrob(plots[[i]])
    widths[[i]] <- grobs[[i]]$widths[2:5]
  }
  maxwidth <- do.call(grid::unit.pmax, widths)
  for (i in 1:length(grobs)) {
    grobs[[i]]$widths[2:5] <- as.list(maxwidth)
  }
  
  if(!is.null(pval)) main <- paste0("Global Schoenfeld Test p: ", signif(pval, 4), "\n")
  else main <- NULL
  
  caption <- attr(plots, "caption")
  
  do.call(gridExtra::grid.arrange, c(grobs, top = main, bottom = caption, newpage = newpage))
}


myInteraction_plot <- function(data, fittedObj, nameObj,int){
  # data=data
  # fittedObj <-fits[[12]]
  # nameObj <- names(fits[12])
  # int=12
  # Get predictors
  factors <- str_split(nameObj,"\\: ") %>% magrittr::extract2(1)
  # print(factors)
  
  if(is_empty(fittedObj$call)){
    print("No interactions evaluated, no model fitted: X matrix deemed to be singular")
    df <- data.frame()
    p <- ggplot(df) + geom_point()
    p <- p+ggtitle(paste0("Inter:",int, "|",factors[1],":",factors[2]))
    p <-p + labs(y = "log Relative Hazard", x=factors[1])
    p <- p + theme_bw(base_family = "Palatino", base_size = 12)
    # p <- NULL
  }else{
    dd <- datadist(data)
    options(datadist=dd)  
    
    # Run predict
    pred_df <- eval(parse(text=paste0("Predict(fittedObj,",myclean_P(factors[1]),",",myclean_P(factors[2]),")")
    ))
    factors <- myclean_P(factors)
    # Get data for plotting
    dataPlot <-as.data.frame(pred_df)
    
    ## Here I need to check if both xVar, colorVar are categorical, numerical or one cat one numerical
    cat <- c(NULL,NULL)
    for(i in 1:length(factors)){
      # i=2
      if(!is.numeric(pred_df[[factors[i]]])){# it's categorical
        cat[i]<- i
      }
    }
    
    if(length(cat)==2){# both categorical
      # Both categorical
      xVar=factors[1] # fisrt categorical or one numerical goes here
      yVar="yhat"
      colorVar = factors[2]
      
      xVar.un <- as.name(xVar)
      yVar.un <- as.name(yVar)
      colorVar.un <- as.name(colorVar)
      
      p <- ggplot(dataPlot, aes(x=!!xVar.un, y=!!yVar.un, color=!!colorVar.un, group=!!colorVar.un))
      p <- p + geom_point() + geom_line()
      p <- p + geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                             position=position_dodge(0.05))
      
      
    }else if(length(cat)==1){# one only categorical
      xVar=setdiff(factors,factors[cat])
      yVar="yhat"
      colorVar = factors[cat]
      xVar.un <- as.name(xVar)
      yVar.un <- as.name(yVar)
      colorVar.un <- as.name(colorVar)
      
      p <- ggplot(dataPlot, aes(x=!!xVar.un, y=!!yVar.un, color=!!colorVar.un, group=!!colorVar.un))
      p <- p + geom_line()
      p <- p + geom_ribbon(aes(ymin=lower, ymax=upper,color=!!colorVar.un, fill=!!colorVar.un), alpha=0.1,linetype="dashed")
      
    }else{# both numeric
      xVar=factors[1] # fisrt categorical or one numerical goes here
      yVar="yhat"
      colorVar = factors[2]
      
      xVar.un <- as.name(xVar)
      yVar.un <- as.name(yVar)
      colorVar.un <- as.name(colorVar)
      
      Mean.cont <- mean(dataPlot[[colorVar]], na.rm=T)
      SD = sd(dataPlot[[colorVar]])
      dataPlot$Cont.cat <- ifelse(dataPlot[[colorVar]]>(Mean.cont + SD), 'High',
                                  ifelse(dataPlot[[colorVar]]<(Mean.cont - SD), 'Low','Mid'))
      
      dataPlot.fin <-dataPlot[!duplicated(dataPlot%>% 
                                            select(!!xVar.un,Cont.cat)),]
      
      colnames(dataPlot.fin)[6] <- paste0(colorVar,".cat")
      colorVar = paste0(colorVar,".cat")
      colorVar.un <- as.name(colorVar)
      
      
      p <- ggplot(dataPlot.fin, aes(x=!!xVar.un, y=!!yVar.un, color=!!colorVar.un, group=!!colorVar.un))
      p <- p + geom_line()
      p <- p + geom_ribbon(aes(ymin=lower, ymax=upper,color=!!colorVar.un), alpha=0.1,linetype="dashed", fill="gray")
      
      
      
    }
    p <- p+ggtitle(paste0("Inter:",int))
    p<-p + labs(y = "log Relative Hazard")
    p <- p + theme_bw(base_family = "Palatino", base_size = 12)
  }
  
  p 
  
}


plotExploreInteractions <- function(features,data, problematic=NULL,toRemove=NULL){
  # features = feats1
  # data=data.complete
  # problematic = 
  # toRemove = 
  
  
  feat.split <- str_split(features," \\+ ") %>% magrittr::extract2(1)
  featsCombos <-split(t(combn(feat.split ,2)), f = 1:ncol(combn(feat.split ,2)))
  names(featsCombos) <- sapply(featsCombos, function(x) paste(x,collapse = ": "), simplify=FALSE) %>% unlist() %>% setNames(NULL)
  featsCombos_final <- featsCombos
  # # Remove problematic, when you have perfect classification
  # featsCombos_final <- featsCombos[-c(11,13,22,25,29,32,35,38,43,47,50)]
  # 
  # ## Remove first line treatment interactions - problematic for some reason
  idx_strata <-grep("strata|strat", names(featsCombos_final))
  if(!is_empty(idx_strata)){
    featsCombos_final <- featsCombos_final[-idx_strata]
    
  }
  
  # Make fits of two way interactions
  fits <-sapply(featsCombos_final, function(x) analyseCox.interTwoWay(data, x, rms=TRUE), simplify = FALSE)
  fits <-sapply(fits, function(x) x$fit, simplify = FALSE)
  
  InterPlots <-sapply(1:length(fits), function(x) myInteraction_plot(data=data,
                                                                     fittedObj =fits[[x]],
                                                                     nameObj =names(fits[x]),int=x), simplify = FALSE)
  
  resPlots=list(p1=ggarrange(plotlist = InterPlots[1:11],nrow=4, ncol=3),
                p2=ggarrange(plotlist = InterPlots[12:22],nrow=4, ncol=3),
                p3=ggarrange(plotlist = InterPlots[23:33],nrow=4, ncol=3),
                p4=ggarrange(plotlist = InterPlots[34:44],nrow=4, ncol=3),
                p5=ggarrange(plotlist = InterPlots[45:55],nrow=4, ncol=3),
                p6=ggarrange(plotlist = InterPlots[56:66],nrow=4, ncol=3))
  
  res=list(plots=resPlots,
           combos=featsCombos_final)
  res
}

myDCAggPlot <- function(stdcaObj,predictors){
  # stdcaObj=netB.strata.nointer.bw
  # predictors="model"
  ## START HERE
  xstart=0.01
  xstop=0.99
  ymin=-0.05
  nb = stdcaObj$dcaPlot.data$net.benefit
  
  pred.n=length(predictors)
  
  #getting maximum net benefit
  ymax=max(nb[names(nb)!="threshold"],na.rm = TRUE)
  
  xlab="Threshold probability"
  ylab="Net benefit"
  
 
  nb.m <- melt(nb, id.vars = "threshold")
  levels(nb.m$variable)<- c("All","None",predictors)
  y.max <- round(nb.m$value%>% max,1)
  # nb.m$variable <- ifelse(nb.m$variable=="all","All",
  #                         ifelse(nb.m$variable=="none","None",nb.m$variable))
  p <- ggplot(data=nb.m, aes(x=threshold, y=value))
  # Add treat all
  p <- p + geom_line(aes(color = variable, linetype = variable), size=.75)
  p <- p + scale_color_manual(values = c("grey", "black", "black"))#
  p <- p + scale_linetype_manual(values=c("solid", "solid", "dashed"))
  # Add ideal line
  # p <- p + geom_hline(yintercept = 0,color="black",size=1.3) 
  p <- p + ylab(ylab)
  p <- p + xlab(xlab)
  # p <- p + labs(caption = c(sub.1, sub.2))
  p <- p + scale_y_continuous(breaks=seq(0,y.max,0.1))
  # p <- p + ylim(ymin,y.max)
  p <- p + theme_bw(base_family = "Palatino", base_size = 13) #+ labs_pubr(base_size = 14, base_family = "Palatino Linotype")
  
  p <- p + theme(legend.justification = c(1, 1), 
                 legend.position = c(1, 1),
                 legend.box = "horizontal",
                 legend.direction = "vertical",
                 legend.title = element_blank(),
                 legend.text = element_text(size=11),
                 legend.background = element_rect(colour = "black"),
                 # legend.key=element_blank(),
                 axis.text = element_text(size=12),
                 axis.title = element_text(size=13),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 plot.caption = element_text(hjust=c(1, 0)))
  
  
  p
}
  


myDCAggPlot_mult <- function(stdcaObj,predictors,ymin=-0.05){
  # stdcaObj=net_both
  # predictors=c("Parsimonious_model","Complete_model")
  ## START HERE
  xstart=0.01
  xstop=0.99
  
  nb = stdcaObj
  
  pred.n=length(predictors)
  
  #getting maximum net benefit
  ymax=max(nb[names(nb)!="threshold"],na.rm = TRUE)
  
  xlab="Threshold probability"
  ylab="Net benefit"
  
  
  nb.m <- melt(nb, id.vars = "threshold")
  levels(nb.m$variable)<- c("All survive","None survive",predictors)
  y.max <- round(nb.m$value%>% max,1)
  # nb.m$variable <- ifelse(nb.m$variable=="all","All",
  #                         ifelse(nb.m$variable=="none","None",nb.m$variable))
  p <- ggplot(data=nb.m, aes(x=threshold, y=value))
  # Add treat all
  p <- p + geom_line(aes(color = variable, linetype = variable), size=.75)
  p <- p + scale_color_manual(values = c("grey", "black", "indian red", "black"))#
  p <- p + scale_linetype_manual(values=c("solid", "solid", "dashed", "dashed"))
  # Add ideal line
  # p <- p + geom_hline(yintercept = 0,color="black",size=1.3) 
  p <- p + ylab(ylab)
  p <- p + xlab(xlab)
  # p <- p + labs(caption = c(sub.1, sub.2))
  p <- p + scale_y_continuous(breaks=seq(0,y.max,0.1))
 p <- p + ylim(ymin,y.max)
  p <- p + theme_bw(base_family = "Palatino", base_size = 13) #+ labs_pubr(base_size = 14, base_family = "Palatino Linotype")
  
  p <- p + theme(legend.justification = c(1, 1), 
                 legend.position = c(1, 1),
                 legend.box = "horizontal",
                 legend.direction = "vertical",
                 legend.title = element_blank(),
                 legend.text = element_text(size=11),
                 legend.background = element_rect(colour = "black"),
                 # legend.key=element_blank(),
                 axis.text = element_text(size=12),
                 axis.title = element_text(size=13),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 plot.caption = element_text(hjust=c(1, 0)))
  
  
  p
}


myDCA <- function(data, testdata, fitObj, time, status,manualCalc, timePoint, features, predName){
  
  # data=somData
  # testdata=somData
  # fitObj=fit.pooled
  # time=surv_time
  # status=surv_status
  # manualCalc=T
  # timePoint=42
  # features=feats.trans.strata
  # predName="MI_model"
  
  time.un <- as.name(time)
  testdata <- testdata %>% mutate(time.un = timePoint)
  
  c.table=coef(fitObj)
  lp = fitObj$linear.predictors
  covmeans = fitObj$means
  fm =modelFormula(time=time, event=status, features) #%>% as.character()
  environment(formula) <- environment()
  
  
  if(isTRUE(manualCalc)){
    predRisk <- mypredictRisk.coxph(coefs.table=c.table,
                                    model.formula=fm,
                                    traindata=data,
                                    lp.train= lp ,#+ sum(cox_model$means*coef(cox_model)),
                                    covMeans=covmeans,
                                    newdata=testdata,
                                    times=timePoint) %>% as.vector() %>% round(4)
    
    
  }else{
    predRisk <- riskRegression::predictRisk(fitObj, newdata =testdata ,times = timePoint) %>% as.vector() %>% round(4)
    
    # predRisk2 <- riskRegression::predictRisk(fitObj, newdata =testdata ,times = timePoint, type="survival",se=T, band=T,confint=T)# %>% as.vector() %>% round(4)
    # predRisk <- riskRegression::predictCox(fitObj, newdata =testdata ,times = timePoint, se=T, confint=T,type="survival")
  }
  
  testdata[[predName]] <- predRisk
  
  dcaPlot <- dca(modelFormula(time=time, event=status, predName),
                 data = testdata,
                 time = timePoint,
                 thresholds = 1:50 / 100) %>%
    plot(smooth = TRUE)
  
  testdata[[status]] = ifelse(testdata[[status]]==TRUE, 1, 0)
  
  
  dcaPlot.bnw <- stdca(data=testdata%>% as.data.frame(), 
                       outcome=status, 
                       ttoutcome=time, timepoint=timePoint, 
                       predictors=predName, probability=T, xstop=.7)
  
  # print(dcaPlot.bnw)
  res <- list(tesdata.out=testdata %>% as.data.frame(),
              dcaPlot=dcaPlot,
              dcaPlot.bnw = dcaPlot.bnw$plot,
              dcaPlot.data=dcaPlot.bnw)
  return(res)
}

