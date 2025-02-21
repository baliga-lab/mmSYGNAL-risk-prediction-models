## Functions for drug therapy analysis

## Get subtype specific network constrained activity
## subby: genetic subtype
## pheno1: ia12 phenotype file. Samples are ordered by Guan risk top=high, bottom =low
## act: data.frame of constrained activity with first columns being drug target symbol.
##          Remaining columns are ia12 sample IDs that match pheno$sample. table can be one 
##          drug target or multiple targets
## quant: Number of quantiles to generate. default is decile (10 values)
##
## Return quantiles of activity for each drug target
getSubtypeQuantile <- function(subby, pheno1, act, quant=num_quantile) {
  
  ## make sure subtype is in pheno file
  if(!(subby %in% colnames(pheno1))) {
    stop(paste("Subtype", subby, "not found in pheno file"))
  }
  
  pheno1 %>% 
    dplyr::filter(!!as.name(subby)==1) ->
    sub_samps
  
  actor <- act[, sub_samps$sample]
  decs <- t(apply(actor, 1, getPieceWiseMean, n=quant))
  rownames(decs) <- unlist(act[,1])
  return(decs)
} ## end getSubtypeQuantile

## Get subtype specific mmSYGNAL risk probability
## subby: genetic subtype
## pheno1: ia12 mmSYGNAL risk predictions. Risk prediction values for each subtype a
##         subject exhibits
## quant: Number of quantiles to generate. default is decile (10 values)
##
## Return quantiles of activity for each risk subtype
getSubtypeRiskQuantile <- function(subby, pheno1, quant=num_quantile) {
  
  ## change the annoying t414 subtype to match
  pheno1 %<>% dplyr::mutate(subtype=plyr::mapvalues(subtype, from="t(4;14)", 
                                                    to="WHSC1")) 
  
  ## make sure subtype is in pheno file
  if(!(subby %in% pheno1$subtype)) {
    stop(paste("Subtype", subby, "not found in pheno file"))
  }
  
  pheno1 %>% 
    dplyr::filter(subtype==subby) %>%
    dplyr::arrange(-pred) ->
    sub_samps
  
  decs <-  getPieceWiseMean(sub_samps$pred, n=quant)
  return(decs)
} ## end getSubtypeRiskQuantile


## calculate deciles of a row of gene expression
## exp - vector of gene expression
## n - number of quantiles 
## returns running mean of lengths ~ length(exp)/n
getPieceWiseMean <- function(exp, n=10) {
  splits <- splitIndices(length(exp), n)
  sapply(splits, function(x, nums=exp) return(mean(nums[x])))
} ## end getPieceWiseMean


## Generate heatmap of constrained network activity
## target: Name of target, either Drug
## pheno: ia12 phenotype file. Samples are ordered by Guan risk top=high, bottom =low
## activity: data.frame of constrained activity with first column being either Drug 
##           or DrugTarget. It must match the type that is defined in 'target'
## num_quantile: number of quantiles to use in plot. default is 10 (deciles)
## sub_pheno: vector of subtype labels that map to pheno
## sub_anno: vector of subtype labels for printing
generateDrugTherapyHeatMap <- function(target, pheno, activity, num_quantile=10,
                                       sub_pheno=subtypes, sub_anno=sub_labels) {
  
  colnames(activity)[1] <- "Target"
  activity %>%
    dplyr::filter(Target %in% target) %>%
    dplyr::distinct() ->
    mini
  
  ## collect subtype specific activity
  all_quans1 <- lapply(sub_pheno, getSubtypeQuantile, pheno1=pheno, act=mini, quant=num_quantile)
  all_quans <-do.call(cbind, all_quans1)
  colnames(all_quans) <- 1:ncol(all_quans)
  
  ## create subtype column annotation
  col_anno <- data.frame(subtype=rep(sub_anno, each=num_quantile))
  rownames(col_anno) <- colnames(all_quans)
  
  ## create gaps between subtype groups
  gaps <- seq(from=num_quantile, to=ncol(all_quans), by=num_quantile)
  ##setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.8, name="vp", just=c("right","top"))), action="prepend")
  
  pheatmap(all_quans, cluster_rows=FALSE, cluster_cols=FALSE, 
           annotation_col=col_anno, 
           annotation_names_col=FALSE,
           gaps_col=gaps,
           show_colnames=FALSE,
           cellheight=10,
           cellwidth=5,
           fontsize=7)
  ##setHook("grid.newpage", NULL, "replace")
  ##grid.text("Actual clinical outcome: high to low",x=0, y=0.2, gp=gpar(fontsize=12))
  
}

## Generate heatmap of constrained network activity
## target: Name of target, either Drug or DrugTarget
## pheno: ia12 phenotype file. Samples are ordered by Guan risk top=high, bottom =low
## activity: data.frame of constrained activity with first column being either Drug 
##           or DrugTarget. It must match the type that is defined in 'target'
## num_quantile: number of quantiles to use in plot. default is 10 (deciles)
## sub_pheno: vector of subtype labels that map to pheno
## sub_anno: vector of subtype labels for printing
## heat_colors: vector of colors for subtype annotation
## c_height: height of cells - play around with it to make nice plot
generateDrugTherapyComplex <- function(target, pheno, activity, num_quantile=10,
                                       sub_pheno=subtypes, sub_anno=sub_labels,
                                       heat_colors=sub_colors, c_height=5) {
  
  colnames(activity)[1] <- "Target"
  activity %>%
    dplyr::filter(Target %in% target) %>%
    dplyr::distinct() ->
    mini
  
  ## collect subtype specific activity
  all_quans1 <- lapply(as.character(sub_pheno), getSubtypeQuantile, pheno1=pheno, act=mini, quant=num_quantile)
  all_quans <-do.call(cbind, all_quans1)
  colnames(all_quans) <- 1:ncol(all_quans)
  
  ## plot deciles
  dec_sub_labs <- rep(sub_anno, each=num_quantile)
  sub_color <- factor(dec_sub_labs,
                      levels=unique(sub_anno),
                      labels=heat_colors)
  names(sub_color) <- dec_sub_labs
  pfs <- HeatmapAnnotation(subtype=dec_sub_labs, subtype=sub_color)
  hm_bin <- Heatmap(all_quans, 
                    column_title="subjects:", 
                    row_title="drug target",
                    cluster_rows=FALSE, 
                    cluster_columns=FALSE,
                    show_column_names=FALSE, 
                    top_annotation=pfs, 
                    column_split=paste(dec_sub_labs, sep="_"),
                    row_split=c(1:nrow(all_quans)),
                    row_gap=unit(0.2, "mm"),
                    height = unit(c_height, "mm"),
                    heatmap_legend_param=list(direction="horizontal",
                                              title="network activity") ) 
  draw(hm_bin, heatmap_legend_side="top", annotation_legend_side="right")    
}


## Generate survival curve for a regulon
## regulon: regulon label
## pheno: ia12 phenotype file. 
## activity: data.frame of regulon activity with first column being regulon label and all
##           other columns are samples. 
## druggy: Data.frame of master file of drug/regulon/program/mmSYGNAL mappings
## patient_regulon: data.frame of regulon activity. first column is regulon
##                  label and second is patient activity
generateKMPlot <- function(target, pheno, activity, druggy,
                           patient_regulon) {
  
  ## get patient activity for this regulon
  colnames(patient_regulon)[2] <- "activity"
  patient_regulon %>%
    dplyr::filter(regulon==target) %>%
    dplyr::select(activity) %>%
    unlist() -> 
    reggie 
  
  reg_activity <- "neutral"
  if(reggie == 1) reg_activity <- "over"
  if(reggie == 0) reg_activity <- "under"
  
  reg_lab <- paste("patient activity:", reg_activity)
  
  ## match color of regulon activity to KM plot
  reg_color <- "green"
  if(reg_activity=="over") {reg_color <- "red"}
  if(reg_activity=="under") {reg_color <- "blue"}
  
  ## collect drug targets for this regulon
  druggy %>%
    dplyr::filter(Regulon==target) %>%
    dplyr::select(TargetSymbol) %>%
    dplyr::distinct() %>%
    unlist() %>%
    sort() %>%
    paste(collapse=", ")-> 
    all_targs
  
  colnames(activity)[1] <- "Target"
  
  ## order pheno and activity file
  risk_sub <- data.frame(network=unlist(activity[activity$Target==target, pheno$sample]),
                         PFS=pheno$PFS, 
                         relapse=pheno$D_PFS_FLAG)
  
  ## map network activity to over, under, neutral activity
  risk_sub$activity <- plyr::mapvalues(risk_sub$network, 
                                       from=c(-1,0,1), 
                                       to=c("under", "neutral", "over"))
  
  ## RMST values: only compare over and under active curves
  risk_sub %>%
    dplyr::filter(network != 0) ->
    dat_in
  
  ## rescale for rmst2
  dat_in$rmst <- plyr::mapvalues(dat_in$network, c(-1,1), c(1,0)) 
  
  res <- rmst2(time=dat_in$PFS, status=dat_in$relapse, 
               arm=dat_in$rmst) 
  ret <- c(RMST=signif(res$unadjusted.result[1, 1], 3), 
           pvalue=signif(res$unadjusted.result[1,4],3))
  rmst_lab <- paste0("RMST: ", ret[1], "\npvalue: ", ret[2])
  
  ## Cox HR ratio
  cox_res <- summary(coxph(Surv(PFS, relapse) ~ network, data=risk_sub))
  cox_ret <- c(HR=signif(cox_res$coefficients[2], 3), 
               pvalue=signif(cox_res$coefficients[5], 3))
  cox_lab <- paste0("Cox HR: ", cox_ret[1], "\npvalue: ", cox_ret[2])
  
  ## survival object
  cox_ia <- survfit(Surv(time=PFS, event=relapse)~ activity,
                    data=risk_sub)
  
  ## generate KM plot
  labs <- names(cox_ia$strata)
  labs <- str_remove(labs, "activity=")
  tab <- cox_ia$strata
  
  plot_ind <- factor(c(rep(labs[1], tab[1]), rep(labs[2], tab[2]), rep(labs[3], tab[3])))
  
  surv_df <- data.frame(SYGNAL_risk=plot_ind, 
                        PFS=cox_ia$surv, 
                        time=cox_ia$time)
  
  ## get med PFS (50% survival)
  med_pfs <- summary(cox_ia)$table[,7]
  names(med_pfs) <- labs
  pfs_tab <- data.frame(pfs=med_pfs, y1=rep(0.5, length(med_pfs)),
                        y2=rep(0,length(med_pfs)),
                        SYGNAL_risk=names(med_pfs))
  
  ggplot(surv_df, aes(x=time, y=PFS, color=SYGNAL_risk)) + 
    geom_point() + geom_line() +
    scale_color_manual(values=c("green", "red", "blue")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1),
          legend.background=element_rect(fill="transparent")) + 
    xlab("progression free survival (months)") + ylab("survival probability") +
    ggtitle(paste("regulon:", target, "\ntargets:", all_targs)) +
    geom_hline(yintercept = 0.5, linetype="dotted") +
    geom_segment(data=pfs_tab, mapping=aes(x=pfs, xend=pfs, y=y1, yend=y2),
                 linetype="dashed") + labs(color="mmSYGNAL activity") + 
    annotate('text', x=11, y=0.30, label=reg_lab, color=reg_color, 
             fontface = 'bold', size=4) +
    annotate('text', x=10, y=0.20, label=rmst_lab, color="black", 
             fontface = 'bold', size=4)  +
    annotate('text', x=10, y=0.05, label=cox_lab, color="black", 
             fontface = 'bold', size=4)  +
    theme(aspect.ratio=4/4) +
    theme(plot.title = element_textbox_simple()) 
  
}

## target: Name of Drug 
## pheno: ia12 phenotype file. Samples are ordered by Guan risk top=high, bottom =low
## activity: data.frame of constrained activity with first column being
##           Drug.  It must match the type that is defined in 'target'
## patient_score: table of drug therapy info from mapping patient regulon
##                and program activity to our master drug therapy info file
## sub_pheno: vector of subtype labels that map to pheno
## sub_anno: vector of subtype labels for printing
plotDrugActivityHistogram <- function(target, pheno, activity,
                                      patient_score,
                                      sub_pheno=subtypes,
                                      sub_anno=sub_labels) {
  
  ## extract target constrained regulon activity
  patient_score %>%
    dplyr::filter(Drug==target) %>%
    dplyr::select(DrugConstrainedRegulonActivity) %>%
    unlist() ->
    p_reg_score
  
  colnames(activity)[1] <- "Target" 
  drug_act <- activity[activity$Target==target, pheno$sample]
  
  drug_hist <- NULL
  for(sub in 1:length(sub_pheno)) {
    pheno %>%
      dplyr::filter(!!as.name(sub_pheno[sub])==1) %>%
      dplyr::select(sample) %>%
      unlist() ->
      sub_samps 
    
    drug_hist <- rbind(drug_hist, 
                       data.frame(subtype=sub_anno[sub], 
                                  activity=unlist(drug_act[1, sub_samps])))
  } ## end for sub
  
  ggplot(drug_hist, aes(x=activity, fill=subtype)) + geom_histogram() +
    ggtitle(target) + geom_vline(xintercept=p_reg_score, linetype="dotted") +
    xlab("drug constrained regulon activity")
}


## Generate plot for drugs constrained network activity
## target: drug label
## pheno: ia12 phenotype file. 
## activity: data.frame of regulon constrained network activity with first column being
## regulon label and all other columns are samples. 
generateDrugPseudoKMPlot <- function(target, pheno, activity) {
  
  colnames(activity)[1] <- "Target"
  
  over_cut <- 0.333
  under_cut <- -0.333
  
  ## order pheno and activity file
  risk_sub <- data.frame(network=unlist(activity[activity$Target==target, pheno$sample]),
                         PFS=pheno$PFS, 
                         PFS_FLAG=pheno$D_PFS_FLAG,
                         relapse=plyr::mapvalues(pheno$D_PFS_FLAG, from=c(0,1), 
                                                 to=c(FALSE, TRUE)),
                         Clinical_Outcome=pheno$GuanScore)
  
  ## map network activity to over, under, neutral activity
  risk_sub$constrained_regulon="neutral"
  risk_sub$constrained_regulon[risk_sub$network >= over_cut] <- "over"
  risk_sub$constrained_regulon[risk_sub$network <= under_cut] <- "under"
  risk_sub$constrained_network <- as.numeric(plyr::mapvalues(risk_sub$constrained_regulon,
                                                             from=c("under", "neutral", "over"),
                                                             to=c(-1,0,1)))
  
  ## RMST values: only compare over and under active curves
  risk_sub %>%
    dplyr::filter(constrained_regulon != "neutral") ->
    dat_in
  
  ## rescale for rmst2
  dat_in$rmst <- plyr::mapvalues(dat_in$constrained_regulon, c("under","over"), c(1,0)) 
  
  res <- rmst2(time=dat_in$PFS, status=dat_in$relapse, 
               arm=dat_in$rmst) 
  ret <- c(RMST=signif(res$unadjusted.result[1, 1], 3), 
           pvalue=signif(res$unadjusted.result[1,4],3))
  rmst_lab <- paste0("RMST: ", ret[1], "\npvalue: ", ret[2])
  
  ## Cox HR ratio
  cox_res <- summary(coxph(Surv(PFS, relapse) ~ constrained_network, data=risk_sub))
  cox_ret <- c(HR=signif(cox_res$coefficients[2], 3), 
               pvalue=signif(cox_res$coefficients[5], 3))
  cox_lab <- paste0("Cox HR: ", cox_ret[1], "\npvalue: ", cox_ret[2])
  
  
  ## survival object
  cox_ia <- survfit(Surv(time=PFS, event=relapse)~ constrained_regulon, data=risk_sub)
  
  ## generate KM plot
  labs <- names(cox_ia$strata)
  labs <- str_remove(labs, "constrained_regulon=")
  tab <- cox_ia$strata
  
  plot_ind <- factor(c(rep(labs[1], tab[1]), rep(labs[2], tab[2]), rep(labs[3], tab[3])))
  
  surv_df <- data.frame(SYGNAL_risk=plot_ind, 
                        PFS=cox_ia$surv, 
                        time=cox_ia$time)
  
  ## get med PFS (50% survival)
  med_pfs <- summary(cox_ia)$table[,7]
  names(med_pfs) <- labs
  pfs_tab <- data.frame(pfs=med_pfs, y1=rep(0.5, length(med_pfs)),
                        y2=rep(0,length(med_pfs)),
                        SYGNAL_risk=names(med_pfs))
  
  ggplot(surv_df, aes(x=time, y=PFS, color=SYGNAL_risk)) + geom_point() + geom_line() +
    scale_color_manual(values=c("green", "red", "blue")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1),
          legend.background=element_rect(fill="transparent")) + 
    xlab("progression free survival (months)") + ylab("survival probability") +
    ggtitle(paste0(target, ": regulon constrained activity")) +
    geom_hline(yintercept = 0.5, linetype="dotted") +
    geom_segment(data=pfs_tab, mapping=aes(x=pfs, xend=pfs, y=y1, yend=y2),
                 linetype="dashed") + labs(color="mmSYGNAL activity") + 
    annotate('text', x=10, y=0.30, label=rmst_lab, color="black", 
             fontface = 'bold', size=4)  +
    annotate('text', x=10, y=0.10, label=cox_lab, color="black", 
             fontface = 'bold', size=4)  +
    theme(aspect.ratio=4/4)
  
  
} ## end generateDrugPseudoKMPlot 

## Generate plot for drugs constrained network activity
## target: drug label
## pheno: ia12 phenotype file. 
## activity: data.frame of regulon constrained network activity with first column being
## regulon label and all other columns are samples. 
generateDrugOutcomePlot <- function(target, pheno, activity) {
  
  colnames(activity)[1] <- "Target"
  
  ## order pheno and activity file
  risk_sub <- data.frame(network=unlist(activity[activity$Target==target, pheno$sample]),
                         PFS=pheno$PFS, 
                         PFS_FLAG=pheno$D_PFS_FLAG,
                         relapse=plyr::mapvalues(pheno$D_PFS_FLAG, from=c(0,1), 
                                                 to=c("no", "yes")),
                         Clinical_Outcome=1-pheno$GuanScore)

  
  ## Cox HR ratio
  cox_res <- summary(coxph(Surv(PFS, PFS_FLAG) ~ network, data=risk_sub))
  cox_ret <- c(HR=signif(cox_res$coefficients[2], 3), 
               pvalue=signif(cox_res$coefficients[5], 3))
  cox_lab <- paste0("Cox HR: ", cox_ret[1], "\npvalue: ", cox_ret[2])
  
  
  ggplot(risk_sub, aes(x=Clinical_Outcome, y=network)) + 
    geom_point(aes(color=relapse)) + 
    xlab("Clinical Outcome: High -> Low disease progression") +
    ylab("drug constrained regulon activity") +
    geom_smooth(method="lm") + 
    annotate('text', x=0.15, y=-0.60, label=cox_lab, color="black", 
             fontface = 'bold', size=4) +
    ggtitle(target) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
}

## pheno: ia12 mmSYGNAL risk file. Columns are sample labels, mmSYGNAL risk prediction 
##        and subtype that a patient exhibits 
## patient risk: data.frame of risk probabilitesf for patient. Column 1 is subtype and
## column 2 is mmSYGNAL risk probability
## patient_lab: patient ID
## num_quantile: number of percentiles to be calculated: default is 10 (decile)
## sub_pheno: vector of subtype labels that map to pheno
## sub_anno: vector of subtype labels for printing
## doPlot: if TRUE print plotpheno
## best_score: value for best model risk prediction based on A,B,C quality subtypes
plotRiskHistogram <- function(pheno, patient_risk, patient_lab="", num_quantile=10,
                              sub_pheno=subtypes, sub_anno=sub_labels,
                              doPlot=TRUE,
                              best_score=NULL) {
  
  pheno$subtype <- plyr::mapvalues(pheno$subtype, 
                                   from=c("all", "del13", "del1p", "amp1q",
                                          "FGFR3", "t(4;14)"),
                                   to=c("all_patients", "del(13)", "del(1p)", "amp(1q)",
                                        "FGFR3", "t(4;14)"))
  
  library(RColorBrewer)
  colorrs <-brewer.pal(6, "Accent")
  names(colorrs) <- c("all_patients", "del(13)", "del(1p)", "amp(1q)", "FGFR3", "t(4;14)")
  
  patient_risk$subtype[patient_risk$subtype=="agnostic"] <- "all_patients"
  patient_risk %>% 
    dplyr::filter(subtype %in% unique(pheno$subtype)) %>%
    dplyr::mutate(y1=-10) %>%
    dplyr::mutate(y2=0) %>%
    dplyr::mutate(color=colorrs[subtype]) %>%
    dplyr::mutate(size=0.7)->
    risk_mat 
  
  if(!is.null(best_score)) {
     risk_mat <- rbind(risk_mat, c("best", best_score, -10, 0, "black", 0.5))
  }
  
  g1 <- ggplot(pheno, aes(x=pred)) + 
    scale_fill_manual(values=colorrs) +
    geom_histogram(aes(fill=subtype), bins=30) +
    geom_rect(data=NULL, aes(xmin=0, xmax=0.5, ymin=-5,ymax=0),
              fill="blue", alpha=0.1) +
    geom_rect(data=NULL, aes(xmin=0.5, xmax=0.6, ymin=-5,ymax=0),
              fill="green", alpha=0.1) +
    geom_rect(data=NULL, aes(xmin=0.6, xmax=1, ymin=-5,ymax=0),
              fill="red", alpha=0.1) +
    ggtitle(paste(patient_lab, "mmSYGNAL risk probability")) + 
    xlab("mmSYGNAL risk probability") +
    labs(color=patient_lab) +
    scale_color_manual(labels=risk_mat$subtype, values=risk_mat$color) +
    geom_vline(data=risk_mat, aes(xintercept=as.numeric(risk), color=subtype), 
               linetype="solid", size=c(rep(1.5, nrow(risk_mat)-1), 0.5)) 
  
  if(doPlot) {
     print(g1)
  }
  
  return(g1)
} ## end 


## target: either drug or drug target label. Must be found in first column of activity
## pheno: ia12 mmSYGNAL risk file. Columns are sample labels, mmSYGNAL risk prediction 
##        and subtype that a patient exhibits 
## patient risk: data.frame of risk probabilites for patient. Column 1 is subtype and
## column 2 is mmSYGNAL risk probability
## activity: data.frame of drug or drug target constrained regulon activity with 
## first column being regulon label and all other columns are samples. 
## num_quantile: number of percentiles to be calculated: default is 10 (decile)
## sub_pheno: vector of subtype labels that map to pheno
## sub_anno: vector of subtype labels for printing
## cell_color: color to outline cell that matches patient risk probability prediction
## heat_colors: vector of colors for subtype annotation
## c_height: height of cells - play around with it to make nice plot
## doHook: Boolean: if True then use hack to print row title
## y_point: Sets y axis of row title
plotRiskHeatmap <- function(target, pheno, patient_risk, activity, num_quantile=10,
                            sub_pheno=subtypes, sub_anno=sub_labels,
                            cell_color="black",
                            heat_colors=sub_colors,
                            c_height=5,
                            doHook=FALSE,
                            y_point=0.3) { 
  
  ## set doHook to FALSE as I put the info into plot title. Keeping the
  ## parameters in function call to not break previous function calls
  doHook <- FALSE
  
  ## rename target column (either drug or drug target)
  colnames(activity)[1] <- "Target" 
  
  ## extract activity for target
  activity %>%
    dplyr::filter(Target %in% target) %>%
    dplyr::distinct() ->
    targ_act
  
  ## relabel the t(4;14) subtype
  pheno$subtype[pheno$subtype=="t(4;14)"] <- "WHSC1"
  
  pheno %<>% dplyr::arrange(subtype, -pred) 
  
  ## get risk and activity quantiles
  all_quans <- NULL
  risk_quans <- list()
  for(sub in sub_pheno) {
    pheno %>%
      dplyr::filter(subtype==sub) ->
      sub_dat
    
    ## order target activity by ia12 patient mmSYGNAL risk
    sub_act <- targ_act[sub_dat$sample]
    
    ## get quantiles
    ##x <- getPieceWiseMean(unlist(sub_act[5,]))
    act_quan <- t(apply(sub_act, 1, getPieceWiseMean, n=num_quantile))
    risk_quan <- getPieceWiseMean(sub_dat$pred, n=num_quantile)
    
    all_quans <- cbind(all_quans, act_quan)
    risk_quans[[sub]] <- risk_quan
  }
  
  colnames(all_quans) <- 1:ncol(all_quans) 
  rownames(all_quans) <- targ_act$Target
  
  ## create subtype column annotation
  col_anno <- data.frame(subtype=rep(sub_anno, each=num_quantile))
  rownames(col_anno) <- colnames(all_quans)
  
  ## create gaps between subtype groups
  gaps <- seq(from=num_quantile, to=ncol(all_quans), by=num_quantile)
  
  ## get cell number for each risk prediction by subtype
  borders <- NULL
  for(sub in 1:length(sub_pheno)) {
    quans <- risk_quans[[sub_pheno[sub]]]
    border1 <- rep("NA", length(quans))
    
    patient_risk %>%
      dplyr::filter(subtype==sub_anno[sub]) ->
      tmp
    
    ## determine location in quantile if subject has risk prediction for this subtype
    if(nrow(tmp) >0) {
      if(tmp$risk > max(quans)) {
        cell_num <- 1
      } else {
        cell_num <- max(which(quans > tmp$risk))
      }
      border1[cell_num] <- cell_color
    }
    borders <- c(borders, border1)
  } ## end for sub
  
  all_borders <- NULL
  for(i in 1:nrow(targ_act)) {
    all_borders <- rbind(all_borders, borders)
  }
  
  
  ## horrible hack to add row title
  if(doHook) {
    setHook("grid.newpage", function()    pushViewport(viewport(x=1,y=1,width=0.9, height=0.8, name="vp", just=c("right","top"))), action="prepend")
  }
  pheatmap(all_quans, cluster_rows=FALSE, cluster_cols=FALSE, 
           annotation_col=col_anno, 
           annotation_names_col=FALSE,
           legend_breaks=seq(from=0, to=1, by=0.25),
           gaps_col=gaps,
           border_color=all_borders,
           show_colnames=FALSE,
           cellheight=10,
           cellwidth=5,
           main="Constrained drug activity ordered from high to low mmSYGNAL risk",
           fontsize=7)
  if(doHook) {
    setHook("grid.newpage", NULL, "replace")
    grid.text("mmSYGNAL risk probability: ordered high (1) to low (0)", 
              x=0.4, y=y_point, gp=gpar(fontsize=12))
  }
  
  
  
  if(FALSE) {
    
    heatmap3(as.matrix(all_quans), Rowv=NA, Colv=NA, balanceColor=TRUE,
             ColSideColors=as.character(sub_color),
             cexRow=0.75, cexCol=0.75)
    
    ## plot deciles
    dec_sub_labs <- rep(sub_anno, each=num_quantile)
    anno_labs <- as.factor(dec_sub_labs)
    ##anno_labs[-seq(from=1, to=length(dec_sub_labs), by=10)] <- ""
    sub_color <- factor(dec_sub_labs,
                        levels=unique(sub_anno),
                        labels=heat_colors)
    names(sub_color) <- dec_sub_labs
    pfs <- HeatmapAnnotation(subtype=dec_sub_labs, subtype=sub_color, 
                             show_annotation_name=FALSE,
                             labels_gp=gpar(col = "black", fontsize = 5),
                             annotation_legend_param=list(labels_gp = gpar(fontsize = 10)))
    pfs <- HeatmapAnnotation(subtype=anno_block(labels=sub_anno))
    
    col_fun = colorRamp2(c(1, 0, -1), c("red", "white", "blue"))
    lgd = Legend(col_fun = col_fun, title = "constrained regulon activity",
                 direction = "horizontal", )
    hm_bin <- Heatmap(all_quans, 
                      column_title="", 
                      row_title="",
                      cluster_rows=FALSE, 
                      cluster_columns=FALSE,
                      show_column_names=FALSE, 
                      ##top_annotation=pfs, 
                      column_split=paste(dec_sub_labs, sep="_"),
                      ##column_km = 3, column_title_gp = gpar(fill = c("red", "blue", "green"), font = 1:3),
                      column_title_gp = gpar(col = "black", 
                                             fontsize = rep(10, 6)),
                      ##column_split=rep(1:6, each=10),
                      ##row_split=c(1:nrow(all_quans)),
                      row_gap=unit(0.2, "mm"),
                      height = unit(c_height, "mm"),
                      ##width= unit(c_height, "mm"))
                      show_heatmap_legend = FALSE)
    ##heatmap_legend_param=list(direction="horizontal",
    ##                          title="constrained regulon activity") ) 
    draw(hm_bin, heatmap_legend_list=lgd, heatmap_legend_side="top")    
  }
  
  
}

## pheno: ia12 mmSYGNAL risk file. Columns are sample labels, mmSYGNAL risk prediction 
##        and subtype that a patient exhibits 
## patient risk: data.frame of risk probabilitesf for patient. Column 1 is subtype and
## column 2 is mmSYGNAL risk probability
## num_quantile: number of percentiles to be calculated: default is 10 (decile)
## sub_pheno: vector of subtype labels that map to pheno
## sub_anno: vector of subtype labels for printing
## cell_color: color to outline cell that matches patient risk probability prediction
## heat_colors: vector of colors for subtype annotation
## c_height: height of cells - play around with it to make nice plot
plotRiskHeatmapComplex <- function(pheno, patient_risk, num_quantile=10,
                                   sub_pheno=subtypes, sub_anno=sub_labels,
                                   cell_color="black",
                                   heat_colors=sub_colors, c_height=5) { 
  
  ## collect subtype specific activity
  all_quans1 <- lapply(sub_pheno, getSubtypeRiskQuantile, pheno1=pheno,
                       quant=num_quantile)
  names(all_quans1) <- sub_anno
  all_quans <-t(data.frame(c(unlist(all_quans1))))
  colnames(all_quans) <- 1:ncol(all_quans)
  rownames(all_quans) <- NULL
  
  ## create subtype column annotation
  col_anno <- data.frame(subtype=rep(sub_anno, each=num_quantile))
  rownames(col_anno) <- colnames(all_quans)
  
  ## create gaps between subtype groups
  gaps <- seq(from=num_quantile, to=ncol(all_quans), by=num_quantile)
  
  ## get cell number for each risk prediction by subtype
  borders <- NULL
  for(sub in sub_anno) {
    quans <- all_quans1[[sub]]
    border1 <- rep("NA", length(quans))
    
    patient_risk %>%
      dplyr::filter(subtype==sub) ->
      tmp
    
    ## determine location in quantile if subject has risk prediction for this subtype
    if(nrow(tmp) >0) {
      cell_num <- max(which(quans > tmp$risk))
      border1[cell_num] <- cell_color
    }
    borders <- c(borders, border1)
  } ## end for sub
  
  
  
  ## plot deciles
  dec_sub_labs <- rep(sub_anno, each=num_quantile)
  anno_labs <- as.factor(dec_sub_labs)
  anno_labs[-seq(from=1, to=length(dec_sub_labs), by=10)] <- ""
  sub_color <- factor(dec_sub_labs,
                      levels=unique(sub_anno),
                      labels=heat_colors)
  names(sub_color) <- dec_sub_labs
  pfs <- HeatmapAnnotation(subtype=dec_sub_labs, subtype=sub_color,
                           text = anno_text(anno_labs, which="column", rot=0, just="center", gp=gpar(size=8)))
  hm_bin <- Heatmap(all_quans, 
                    column_title="", 
                    row_title="",
                    cluster_rows=FALSE, 
                    cluster_columns=FALSE,
                    show_column_names=FALSE, 
                    top_annotation=pfs, 
                    column_split=paste(dec_sub_labs, sep="_"),
                    row_split=c(1:nrow(all_quans)),
                    row_gap=unit(0.2, "mm"),
                    height = unit(c_height, "mm"),
                    heatmap_legend_param=list(direction="horizontal",
                                              title="mmSYGNAL risk probability") ) 
  draw(hm_bin, heatmap_legend_side="top", annotation_legend_side="right")    
}
