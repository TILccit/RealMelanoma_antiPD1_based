#' Get unused and missing packages for RMD or R scripts - does not work very good
#' 
#' Returns age, decimal or not, from single value or vector of strings
check_packages <- function(path,filename,typeFile=".RMD", unused=T, missing=T){
  ####
  # path='D:/BIOINF/PROJECTS/2_ANALYSES/18_GIT_PAPERS/Long_term_antiPD1_based_prognostic_model/scripts/functions/'
  # filename="data_preprocessing_functions"
  # typeFile=".R"
  # Find which packages do used functions belong to ----
  packageload <- p_loaded(all=T) %>% sort() #, character.only=T
  packageload <- search()[!(search() %in% packageload)] |> grep(pattern = "package:", value = TRUE) |> gsub(pattern = "package:", replacement = "")
  
  # R or RMD
  if(typeFile==".RMD"){
    knitr::purl(paste0(path,filename, typeFile), output=paste0(path,filename,".R"))
  }
  # Find functions and packages used
  used.functions <- NCmisc::list.functions.in.file(filename = paste0(path,filename,".R"), alphabetic = FALSE) #|> print()
  
  
  # Find which loaded packages are not used ----
  n <- used.functions |> names()
  used.packages <- unlist(regmatches(n, gregexpr("(?<=package:)[^,\\)\\\"]+", n, perl = TRUE))) %>% unique()
  unused.packages <- packageload[!(packageload %in% used.packages)] %>% sort()#|> print()
  # Find missing packages
  missing.packages <- used.packages[!(used.packages %in% packageload)] %>% sort()#|> print()
  #
  out <- list(used=used.packages,
              unused=unused.packages,
              missing=missing.packages)
  out
}



#' Get age
#' 
#' Returns age, decimal or not, from single value or vector of strings
#' or dates, compared to a reference date defaulting to now. Note that
#' default is NOT the rounded value of decimal age.
#' @param from_date vector or single value of dates or characters
#' @param to_date date when age is to be computed
#' @param dec return decimal age or not
#' @examples
#' get_age("2000-01-01")
#' get_age(lubridate::as_date("2000-01-01"))
#' get_age("2000-01-01","2015-06-15")
#' get_age("2000-01-01",dec = TRUE)
#' get_age(c("2000-01-01","2003-04-12"))
#' get_age(c("2000-01-01","2003-04-12"),dec = TRUE)
get_age <- function(from_date,to_date = lubridate::now(),dec = FALSE,decPoints=1){
  if(is.character(from_date)) from_date <- lubridate::as_date(from_date)
  if(is.character(to_date))   to_date   <- lubridate::as_date(to_date)
  if (dec) { age <- lubridate::interval(start = from_date, end = to_date)/(lubridate::days(365)+lubridate::hours(6))
  age <- round(age,decPoints)
  } else   { age <- lubridate::year(lubridate::as.period(lubridate::interval(start = from_date, end = to_date)))}
  age
}

#' Plot for missing data
#' 
ggplot_missing <- function(x){

x %>% 
  is.na %>%
  melt %>%
  ggplot(data = .,
         aes(x = Var2,
             y = Var1)) +
  geom_raster(aes(fill = value)) +
  scale_fill_grey(name = "",
                  labels = c("Present","Missing")) +
  theme_minimal() + 
  theme(axis.text.x  = element_text(angle=45, vjust=0.5)) + 
  labs(x = "Variables in Dataset",
       y = "Rows / observations")
}

#' Barplots for nominal variables
#' 
barplot_nominal <- function(DF,variable,angle, yLabel){
  df<-as.data.frame(table(DF[[variable]]))
  p <- ggplot(df, aes(x=Var1, y=Freq, fill=Var1))
  p <- p + geom_bar(stat="identity", color="black")
  p <- p+scale_fill_brewer(palette="Dark2")
  p <- p + theme_ipsum(base_family = "Palatino", base_size = 10)
  p <- p + theme(axis.text.x = element_text(angle = angle, hjust = 1, vjust = 1,size = 10),
                 axis.title.x = element_text(size=10, face="bold"),
                 axis.text.y = element_text( hjust = 1, vjust = 1,size = 10),
                 axis.title.y = element_text(size=11, face="bold"),
                 legend.title = element_text(size=11,color="black"),
                 legend.text = element_text(size = 10),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"))
  p <- p + labs(fill=variable,
                title = paste0(""),
                x = variable,
                y = yLabel
  )
  p
}

#' Barplots for nominal variables with facets
#' 
barplot_nominal_fac <- function(df,variable,angle,facetVariable=NULL){
  
  # df <- dataDF
  # variable <-"Age categories"
  # angle <- 45
  # facetVariable="PS = 0"
  # print(facetVariable)
  variable.un <- as.name(variable)
  facetVariable.un <- as.name(facetVariable)
  df <- df %>% select(!!variable.un, !!facetVariable.un)
  df<-as.data.frame(table(df))
  colnames(df) <- c(variable, facetVariable,"Freq")
  variable <- rlang::sym(variable)
  facetVariable <- rlang::sym(facetVariable)
  countVar <-rlang::sym("Freq")
  # p <- ggplot(df, aes_string(x=str_replace(variable," ","."), y="Freq", fill=str_replace(facetVariable," ",".")))
  p <- ggplot(df, aes(x=!!variable, y=!!countVar, fill=!!facetVariable))
  
  p <- p + geom_bar(stat="identity", color="black")
  # p <- p+scale_fill_brewer(palette="Dark2")
  p <- p+scale_fill_grey()
  # if(!is.null(facetVariable)){
  #   p <- p + facet_wrap(facetVariable,ncol=2,scales='free')
  # }
  p <- p + theme_ipsum(base_family = "Palatino", base_size = 10)
  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 10),
                 axis.title.x = element_text(size=10, face="bold"),
                 axis.text.y = element_text( hjust = 1, vjust = 1,size = 10),
                 axis.title.y = element_text(size=11, face="bold"),
                 legend.title = element_text(size=11,color="black"),
                 legend.text = element_text(size = 10),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"))
  # p <- p + labs(fill=variable,
  #               title = paste0(""),
  #               x = variable,
  #               y = yLabel
  # )
  print(p)
}

#' Pie chart, allowing for filtering variables and grouping
#' 
myPieChart <- function(DF, groupVar,filterVar=NULL,filterSel=NULL){
  # DF <- dataDF
  # groupVar <-"firstLine"
  # filterVar <-"eligible"
  # filterSel <-"eligible"
  # facetVar <- "Year_of_diagnosis"
  # 
  # filterVar.un <- as.name(filterVar)
  groupVar.un <- as.name(groupVar)
  # facetVar.un <- as.name(facetVar)
  # # Filter
  # data.n <- DF %>% dplyr::filter(!!filterVar.un==filterSel) %>% droplevels()
  data.n <- DF
  # Data transformation
  
  df <- data.n %>% 
    dplyr::group_by(!!groupVar.un) %>% # Variable to be transformed
    count() %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(perc = `n` / sum(`n`)) %>% 
    dplyr::arrange(perc) %>%
    dplyr::mutate(labels = scales::percent(perc))
  
  # myPalette <- brewer.pal(5, "Set2") 
  # myPalette <- pal_uchicago("default")(9)
  myPalette <- ccf_cols(c("surf_crest","silver_sand","ccf_green",  "maroon_flush"))
  p <-ggplot(df, aes(x = "", y = perc, fill = !!groupVar.un)) +
    geom_col() +
    geom_label(aes(label = labels),
               position = position_stack(vjust = 0.5),
               show.legend = FALSE, size=8) +
    # coord_polar(theta = "y", start=0) + scale_fill_brewer(palette="Set2") + theme_ipsum(base_family = "Palatino", base_size = 10) + theme_void()
    coord_polar(theta = "y", start=0) + scale_fill_manual(values=myPalette) + 
    theme_ipsum(base_family = "Palatino", base_size = 15) + 
    theme_void() + 
    theme(legend.text = element_text(size = 13),
          legend.title = element_blank(),
          legend.position="bottom")
  
  # + theme(
  #     # axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
  #     axis.title.x = element_blank(),
  #     # axis.text.y = element_text( hjust = 1, vjust = 1,size = 10),
  #     axis.title.y = element_blank(),
  #     # legend.title = element_text(size=11,color="black"),
  #     # legend.text = element_text(size = 10),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     panel.background = element_blank(),
  #     # axis.line = element_line(colour = "black"),
  #     legend.position="bottom")
  p <- p + ggtitle(filterSel)
  # if(!is.null(facetVar)){
  #   p <- p + facet_wrap(facetVar,ncol=plotCols,scales='free')
  # }
  p
}

#' Pie chart, allowing for filtering variables and grouping, setting legend title and coloring
#' 
myPieChart2 <- function(DF, groupVar,filterVar,filterSel, legendTitle,fontSize = 1,fillPal=NULL, repel=TRUE){
  # DF <- dataDF
  # groupVar <-"eligible"
  # filterVar <-"Year_of_diagnosis"
  # filterSel <-"2012"
  # facetVar <- "Year_of_diagnosis"
  # 
  filterVar.un <- as.name(filterVar)
  groupVar.un <- as.name(groupVar)
  # facetVar.un <- as.name(facetVar)
  # Filter
  data.n <- DF %>% dplyr::filter(!!filterVar.un==filterSel) %>% droplevels()
  # Data transformation
  
  df <- data.n %>% 
    dplyr::group_by(!!groupVar.un) %>% # Variable to be transformed
    count() %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(perc = `n` / sum(`n`)) %>% 
    dplyr::arrange(perc) %>%
    dplyr::mutate(labels = scales::percent(perc))
  
  df <- df %>% mutate(perc.n = perc*100,
                      csum = rev(cumsum(rev(perc.n))), 
                      pos = perc.n/2 + lead(csum, 1),
                      pos = if_else(is.na(pos), perc.n/2, pos))
  # myPalette <- brewer.pal(5, "Set2") 
  
  
  # p <- ggplot(df, aes(x = "" , y = perc.n, fill = fct_inorder(!!groupVar.un))) +
  #   geom_col(width = 1, color = 1) +
  #   coord_polar(theta = "y") +
  #   scale_fill_brewer(palette = "Pastel1") +
  #   geom_label_repel(data = df,
  #                    aes(y = pos, label = labels),#paste0(round(perc.n,2), "%"))
  #                    size = 4, nudge_x = .6*fontSize, show.legend = FALSE) +
  #   guides(fill = guide_legend(title = legendTitle)) +
  #   theme_void()
  # p <- p + ggtitle(filterSel)
  # # if(!is.null(facetVar)){
  # #   p <- p + facet_wrap(facetVar,ncol=plotCols,scales='free')
  # # }
  # p
  # 
  # p <- ggplot(df, aes(x = "" , y = perc.n, fill = fct_inorder(!!groupVar.un))) +
  #   geom_col(width = 1, color = 1) + 
  #   geom_text(aes(label = round(perc.n,2)),
  #             position = position_stack(vjust = 0.5)) +
  #   coord_polar(theta = "y") +
  #   guides(fill = guide_legend(title = legendTitle)) +
  #   scale_y_continuous(breaks = df$pos, labels = df[[groupVar]]) +
  #   scale_fill_brewer(palette = "Pastel1") +
  #   theme(axis.ticks = element_blank(),
  #         axis.title = element_blank(),
  #         axis.text = element_text(size = 15),
  #         legend.position = "none", # Removes the legend
  #         panel.background = element_rect(fill = "white"))
  # p
  p <- ggplot(df, aes(x = "" , y = perc.n, fill = fct_inorder(!!groupVar.un))) +
    geom_col(width = 1, color = 1) + 
    geom_text(aes(label = round(perc.n,2)),
              position = position_stack(vjust = 0.5),label.size = 3.5,fontface="bold",size=6*fontSize) +
    coord_polar(theta = "y") +
    guides(fill = guide_legend(title = legendTitle)) 
  if(isTRUE(repel)){
    p <- p +
      geom_label_repel(data = df,
                       aes(y = pos, label = df[[groupVar]]),#paste0(round(perc.n,2), "%"))
                       size = 4*fontSize, nudge_x = 1,nudge_y = 0.3, show.legend = FALSE)
  } 
  p <-p + ggtitle(filterSel) 
  # scale_y_continuous(breaks = df$pos, labels = df[[groupVar]]) +
  if(!is.null(fillPal)){
    p <- p + scale_fill_manual(values = fillPal)
  }else{
    scale_fill_brewer(palette = "Pastel2")
  }
  p <- p + theme_ipsum(base_family = "Palatino", base_size = 25) + 
    theme(plot.title = element_text(size = 45*fontSize, face = "bold"), legend.text =element_text(size = 45*fontSize) ) +
    theme_void(base_family = "Palatino", base_size = 15*fontSize)
  # theme(axis.ticks = element_blank(),
  #       axis.title = element_blank(),
  #       axis.text = element_text(size = 15),
  #       legend.position = "none", # Removes the legend
  #       panel.background = element_rect(fill = "white"))
  
  # p <- p 
  p
}

#' Table styling function,add stars to p values
#' 
fmt_pvalue_with_stars <- function(x) {
  dplyr::case_when(
    x < 0.001 ~ paste0(style_pvalue(x), " ***"),
    x < 0.01 ~ paste0(style_pvalue(x), " **"),
    x < 0.05 ~ paste0(style_pvalue(x), " *"),
    TRUE ~ style_pvalue(x, digits = 3)
  )
}

#' Extract Kaplan Meier plot, survival probabilities and risk tables
#'
getKMplot2<- function(DF,timeColumn,eventColumn, byColumn, xLabel, yLabel, timeUnit="months", 
                      timeBreak=12, addPHR=FALSE,riskTab=FALSE,addMedian=FALSE, tbH=1.5,scaleX="m_y",xmax=122,tbX=120, addLegend=FALSE, LegLabel=NULL,legpos=NULL){
  # 
  # DF = dataDF.16.18
  # timeColumn = surv_time
  # eventColumn = surv_status
  # byColumn = "BRAF_YEAR"
  # xLabel = xLabel
  # yLabel = yLabel
  # timeUnit="months"
  # timeBreak=12
  # addPHR=FALSE
  # riskTab=TRUE
  # addMedian=TRUE
  # tbH=1.5
  # scaleX="m_y"
  # xmax=132
  # tbX=139
  # addLegend=FALSE
  # LegLabel=NULL
  # legpos = c(.85,0.4)
  # print("START")
  timeColumn.un = as.name(timeColumn)
  eventColumn.un = as.name(eventColumn)
  byColumn.un = as.name(byColumn)
  
  dataObj <- deparse(substitute(DF))
  survFit.obJ<-eval(parse(text=paste0("survfit(Surv(",timeColumn,", ",eventColumn,") ~ ",byColumn,", data = ",dataObj,")")))
  survpval <-eval(parse(text=paste0("surv_pvalue(survFit.obJ",", data = ",dataObj,")")))
  
  # survFit.obJ <- SurvFun(fun.time=timeColumn, fun.event=eventColumn, grouping = byColumn, fun.dat=DF)
  # Single plot
  if(addLegend){
    survPlot <- ggsurvplot(
      survFit.obJ,
      data = DF,
      pval = FALSE,
      pval.method = FALSE,
      xscale=scaleX,
      break.time.by=timeBreak,
      xlim=c(0,xmax),
      legend.title=LegLabel,
      legend.labs = survFit.obJ$strata %>% names() %>% str_replace_all(paste0(byColumn,"="),""),
      conf.int = FALSE,
      xlab = xLabel,
      ylab = yLabel,
      # ggtheme = theme_minimal(base_family = "Palatino", base_size = 13),
      ggtheme = theme_pubclean(base_family = "Palatino", base_size = 13),
      surv.median.line = "hv",
      size = 0.75,
      # palette = c("grey","black"),
      palette = grey_range(DF[[byColumn]] %>% levels() %>% length()),#"uchicago",
      font.main=16,
      font.legend = 14,
      # Add risk table
      risk.table = TRUE,
      tables.height = .2,
      tables.theme = theme_cleantable(base_size = 11)
    )
  }else{
    # print("HERE")
    survPlot <- ggsurvplot(
      survFit.obJ,
      data = DF,
      pval = FALSE,
      pval.method = FALSE,
      xscale=scaleX,
      break.time.by=timeBreak,
      xlim=c(0,xmax),
      # legend.title=analyzeSurvObj$by,
      legend.labs = survFit.obJ$strata %>% names() %>% str_replace_all(paste0(byColumn,"="),""),
      conf.int = FALSE,
      xlab = xLabel,
      ylab = yLabel,
      # ggtheme = theme_minimal(base_family = "Palatino", base_size = 13),
      ggtheme = theme_pubclean(base_family = "Palatino", base_size = 13),
      surv.median.line = "hv",
      size = 0.75,
      # palette = c("grey","black"),
      palette = grey_range(DF[[byColumn]] %>% levels() %>% length()),#"uchicago",
      font.main=16,
      font.legend = 14,
      # Add risk table
      risk.table = TRUE,
      tables.height = .2,
      tables.theme = theme_cleantable(base_size = 11)
    )
  }
  
  if(addLegend){
    
    if(!is_empty(legpos)){
      # print("LEGPOS")
      kmPlot <- survPlot$plot + theme_pubclean(base_family = "Palatino", base_size = 13) +
        theme(legend.position = legpos,
              legend.box = "horizontal",
              legend.direction = "vertical",
              legend.title = element_text(size=12),
              legend.text = element_text(size=12),
              axis.text.x = element_text(size=11),
              legend.key=element_blank())
    }else{
      
      kmPlot <- survPlot$plot + theme_pubclean(base_family = "Palatino", base_size = 13) +
        theme(legend.position = "right",
              legend.box = "horizontal",
              legend.direction = "vertical",
              legend.title = element_text(size=12),
              legend.text = element_text(size=12),
              axis.text.x = element_text(size=11),
              legend.key=element_blank())
    }
    
  }else{
    if(!is_empty(legpos)){
      # print("LEGPOS")
      kmPlot <- survPlot$plot + theme_pubclean(base_family = "Palatino", base_size = 13) +
        theme(legend.position = legpos,
              legend.box = "horizontal",
              legend.direction = "vertical",
              legend.title = element_blank(),
              legend.text = element_text(size=12),
              axis.text.x = element_text(size=11),
              legend.key=element_blank())
    }else{
      
      kmPlot <- survPlot$plot + theme_pubclean(base_family = "Palatino", base_size = 13) +
        theme(legend.position = "right",
              legend.box = "horizontal",
              legend.direction = "vertical",
              legend.title = element_blank(),
              legend.text = element_text(size=12),
              axis.text.x = element_text(size=11),
              legend.key=element_blank())
    }
    
  }
  
  
  riskPlot <-survPlot$table + ggplot2::theme_bw(base_size = 11)+ theme(
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=11),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )
  riskPlot <-riskPlot + geom_label(size=3.5,label.size = NA, fill = "white")
  if(isTRUE(addPHR)){
    if(length(analyzeSurvObj$factorLabels)>2){
      tb <- survivalDFs$hazardRatios %>% as.data.frame() %>% dplyr::filter(!str_detect(label, paste0(analyzeSurvObj$factorLabels[1]," vs."))) %>%
        dplyr::mutate(HR=paste0("HR = ",HR," (",Lower.CI,"-",Upper.CI,") ")) %>%
        dplyr::mutate(p=paste0("p = ",p)) %>%
        dplyr::select(label, HR, p) %>% as_tibble()
      # colnames(tb ) <- c("","","")
      # tb <- tb %>% as_tibble()
      # Filter out rows that start with reference level
      data.tb <- tibble(x = max(DF[[timeColumn]], na.rm = T), #-(max(DF[[timeColumn]], na.rm = T)/48)
                        y = c(.99), tb = list(tb))
      kmPlot.f <- kmPlot +
        geom_table(data = data.tb, aes(x, y, label = tb),
                   table.theme = ttheme_minimal,
                   size = 3, colour = "black")
      
      
    }else{
      # Make annotations
      annotations <- data.frame(
        # xpos = c(-Inf),
        # ypos =  c(-Inf),
        xpos = max(DF[[timeColumn]], na.rm = T)-(max(DF[[timeColumn]], na.rm = T)/24),
        # xpos = max(DF[[timeColumn]], na.rm = T),
        ypos = c(.9) ,
        annotateText = c(paste0("HR = ",survivalDFs$hazardRatios[1,]$HR,"(",survivalDFs$hazardRatios[1,]$Lower.CI,"-",survivalDFs$hazardRatios[1,]$Upper.CI,")\n",survivalDFs$cohortMetadata$log.rank.p)),
        hjustvar = c(1) ,
        vjustvar = c(0))
      
      kmPlot.f <- kmPlot+
        ggplot2::geom_label(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,
                                                    vjust=vjustvar,label=annotateText),family="Palatino",size=3.5)
    }
    
  }else{
    kmPlot.f <- kmPlot
  }
  
  if(addMedian){
    medianSurv <- tbl_survfit(
      DF,
      y = Surv(!!timeColumn.un,!!eventColumn.un),
      include = c(!!byColumn.un),
      probs = 0.5,
      label_header = "**Median Survival**",
    )
    
    medianSurv.df <-  medianSurv$table_body %>%
      filter(label %in% DF[[byColumn]]) %>%
      dplyr::select(label, stat_1) %>%
      as.data.frame() %>%
      column_to_rownames(var="label") %>%
      t() %>%
      as_tibble() %>%
      as.data.frame()
    
    # if(addLegend){
    #   rownames(medianSurv.df) <- paste0(LegLabel,"\n","Median OS\n(",timeUnit,")")
    #   
    # }else{
    #   rownames(medianSurv.df) <- paste0("Median OS\n(",timeUnit,")")
    #   
    # }
    rownames(medianSurv.df) <- paste0("Median OS\n(",timeUnit,")")
    medianSurv.df[1,] <-medianSurv.df[1,] %>% str_replace_all(" to ",",")
    medianSurv.df[1,] <-medianSurv.df[1,] %>% str_replace_all(" ","\n")
    medianSurv.df <- medianSurv.df %>% rownames_to_column(var=" ")
    
    data.tb <- tibble(x = tbX + 1.2, #-(max(DF[[timeColumn]], na.rm = T)/48)
                      y = c(.99), tb = list(medianSurv.df))
    kmPlot.f <- kmPlot +
      geom_table(data = data.tb, aes(x, y, label = list(medianSurv.df)),
                 table.theme = ttheme_minimal,
                 size = 3.8, colour = "black")
    
  }
  
  
  
  if(isTRUE(riskTab)){
    kmPlotFull <- kmPlot.f /riskPlot + patchwork::plot_layout(ncol=1,nrow=2,heights = c(10-tbH,tbH)) #+ theme(legend.position = "bottom")
  }else{
    kmPlotFull <- kmPlot.f
    
  }
  
  # Surv probs
  if(scaleX=="m_y"){
    surv_prob <- tbl_survfit(
      survFit.obJ,
      times = c(12,24,36,48,60,120),
      label_header = "**{time/12} yr**",statistic="{estimate}", label=NULL
    )|>
      modify_header(label = "")
  }else{
    surv_prob <- tbl_survfit(
      survFit.obJ,
      times = c(1,2,3,4,5,10),
      label_header = "**{time} yr**",statistic="{estimate}", label=NULL
    )|>
      modify_header(label = "")
  }
  
  # surv_prob$table_body$label[1] <- 
  
  res <- list(km=kmPlotFull,
              s_table=surv_prob,
              risk=survPlot$table)
  res
}