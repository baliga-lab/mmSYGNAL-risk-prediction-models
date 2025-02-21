library(data.table)

#' Given a drug response score return whether the score is a responder (score is
#' greater than or equal cut_up), non-responder (score is less than or equal 
#' cut_down), or neutral (score is between cut_up and cut_down).
#' 
#' @param x drug response score (either scalar or vector)
#' @param cut_up if x is greater than or equal this value: responder
#' @param cut_down x is less than or equal this value: non-responder
#' 
#' @return character scalar or vector labelling each element responder, nonresponder,
#' or neutral
#' @export
getDrugResponse <- function(x, cut_up, cut_down) {
   ret <- rep("neutral", length(x))
   ret[x >= cut_up] <- "responder"
   ret[x <= cut_down] <- "nonresponder"
   return(ret)
} ## end getDrugResponse



## calculate mean of vector and remove NAs
meanNA <- function(x) mean(x, na.rm=TRUE)

## get vector of regulons for drug that match regulon activity
## mat - data.table of drug therapy info (getDrugTherapyActivity)
## drug - name of drug
## activity - value of activity (-1, 1)
getRegulonsByActivity <- function(drug, mat, value) {
  mat %>% 
    dplyr::filter(Drug==drug) %>%
    dplyr::filter(reg_activity==value) ->
    tmp
  regs <- ifelse(nrow(tmp) == 0, NA, paste(sort(unique(tmp$Regulon[tmp$reg_activity==value])), 
                                           collapse=","))
  
  return(c(Drug=drug, regulons=regs))
}

## get vector of programs for drug that match regulon activity
## mat - data.table of drug therapy info (getDrugTherapyActivity)
## drug - name of drug
## activity - value of activity (-1, 1)
getProgramsByActivity <- function(drug, mat, value) {
  mat %>% 
    dplyr::filter(Drug==drug) %>%
    dplyr::filter(prog_activity==value) ->
    tmp
  regs <- ifelse(nrow(tmp) == 0, NA, paste(sort(unique(tmp$Program[tmp$prog_activity==value])), 
                                           collapse=","))
  
  return(c(Drug=drug, programs=regs))
}


## Generate drug target specific program activity  
## dug_mat - data.frame of drug target info (see drug_information_processing.Rmd
##           output/drugs/drug_map_me_to_patient.csv"
## program - data.frame of program activity for one subject col1: Program label
##                                                          col2: program activity
## regulon - data.frame of program activity for one subject col1: regulon label
##                                                          col2: regulon activity
## doOverActive - if TRUE, calculate constrained regulon activity with only active regulons
##                and programs else use both under and over activity
## fdr_cut - fdr cutoff for program disease relevance
getDrugTargetTherapyActivity <- function(drug_mat, program, regulon, doOverActivity=TRUE, 
                                         fdr_cut=0.2) {
  
  ## rename inputs for joining
  colnames(program) <- c("program", "prog_activity")
  colnames(regulon) <- c("regulon", "reg_activity")
  
  ## convert under active to neutral
  if(doOverActivity) {
       program$program[program$program==-1] <- 0
       regulon$regulon[regulon$regulon==-1] <- 0
  }
  
  ## join regulon and program activity to drug info
  drug_mat %>%
    dplyr::left_join(regulon,by=c("Regulon"="regulon")) %>%
    dplyr::left_join(program, by=c("Program"="program")) %>%
    data.table() ->
    drug_rel
  
  ##drug_rel %>%
  ##  dplyr::filter(Drug=="ACAMPROSATE") ->
  ##  tmp
  
  druggies <- unique(drug_rel$Drug)
  
  
  if(FALSE) {
    
    ## get list of all drug targets for each drug
    drug_rel[, DrugTargetAll:=paste0(sort(unique(TargetSymbol)), collapse=","),
             by=.(Drug, TargetSymbol)]
    
    ## get list of all regulons for each drug 
    drug_rel[, DrugRegulonAll:=paste0(sort(unique(Regulon)), collapse=","),
             by=.(Drug, TargetSymbol)]
    
  ## get list of regulons that are over active
  act_regs <- data.frame(t(sapply(druggies, getRegulonsByActivity, mat=drug_rel, value=1)))
  drug_rel %>%
    dplyr::left_join(act_regs, by=c("Drug"="Drug")) %>%
    dplyr::rename(DrugRegulonOverActive=regulons) ->
    drug_rel
  
  ## get list of regulons that are under active
  act_regs <- data.frame(t(sapply(druggies, getRegulonsByActivity, mat=drug_rel, value=-1)))
  drug_rel %>%
    dplyr::left_join(act_regs, by=c("Drug"="Drug")) %>%
    dplyr::rename(DrugRegulonUnderActive=regulons) ->
    drug_rel 
  
  ## get list of programs that are over active
  act_progs <- data.frame(t(sapply(druggies, getProgramsByActivity, mat=drug_rel, value=1)))
  drug_rel %>%
    dplyr::left_join(act_progs, by=c("Drug"="Drug")) %>%
    dplyr::rename(DrugProgramOverActive=programs) ->
    drug_rel
  
  ## get list of programs that are under active
  act_progs <- data.frame(t(sapply(druggies, getProgramsByActivity, mat=drug_rel, value=-1)))
  drug_rel %>%
    dplyr::left_join(act_progs, by=c("Drug"="Drug")) %>%
    dplyr::rename(DrugProgramUnderActive=programs) ->
    drug_rel 
  
  ## get list of all programs for each drug 
  drug_rel[, DrugProgramAll:=paste0(sort(unique(Program)), collapse=","),
           by=.(Drug)]
  
  } ## end if FALSE
  
 
  ## TODO: check sign of mean of RSMT program direction
  ## calculate regulon disease relevant statistics
  drug_rel[, DrugMeanRMSTRegulon:=meanNA(RMST_UnderMinusOver_reg), by=.(Drug, TargetSymbol)]
  drug_rel[, DrugRMSTRegulon_MeanDirection:=meanNA(RMSTRegulonDirection), by=.(Drug, TargetSymbol)]
  drug_rel[, DrugRMSTRegulonMeanPvalue:=meanNA(RMST_reg_pval), by=.(Drug, TargetSymbol)]
  
  ## calculate program disease relevant statistics
  drug_rel[, DrugRMSTProgram_MeanDirection:=meanNA(RMSTProgramDirection), by=.(Drug, TargetSymbol)]
  drug_rel[, DrugMeanRMSTProgram:=meanNA(RMST_UnderMinusOver_prog), by=.(Drug, TargetSymbol)]
  drug_rel[, DrugRMSTProgramMeanPvalue:=meanNA(RMST_prog_pval), by=.(Drug, TargetSymbol)]
  
  ## calculate drug constrained regulon activity
  drug_rel[, DrugConstrainedRegulonActivity:=meanNA(reg_activity),
           by=.(Drug, TargetSymbol)]
  
  ## calculate drug constrained program activity
  drug_rel[, DrugConstrainedProgramActivity:=meanNA(prog_activity),
           by=.(Drug, TargetSymbol)] 
  
  drug_rel %>%
    dplyr::select(-Regulon, -Program, -RegulonRegulator, -RegulatorSymbol, 
                  -RegulatorTarget, -cox_HR_reg, -cox_HR_reg_pvalue, -cox_HR_reg_fdr,
                  -RMST_UnderMinusOver_reg, -RMST_reg_pval, -RMST_reg_fdr, -cox_HR_prog,
                  -cox_HR_prog_pvalue, -cox_HR_prog_fdr, -RMST_UnderMinusOver_prog,
                  -RMST_prog_pval, -RMST_prog_fdr, -Cox_RegulonDirection,
                  -CoxProgramDirection, -RMSTRegulonDirection, -RMSTProgramDirection,
                  -drug_direction, -RegulatorRegulon_Spearman_R, -RegulonMutationAll,
                  -reg_activity, -prog_activity) %>%
    distinct() %>%
    dplyr::select(everything(),DrugConstrainedRegulonActivity,DrugConstrainedProgramActivity)  ->
    drug_sum

  return(drug_sum)
}


## Generate drug specific program activity - this is called in 
## drug_information_processing_updated.csv to generate IA12 DCNA
## dug_mat - data.frame of drug target info (see drug_information_processing.Rmd
##           output/drugs/drug_map_me_to_patient.csv"
## program - data.frame of program activity for one subject col1: Program label
##                                                          col2: program activity
## regulon - data.frame of program activity for one subject col1: regulon label
##                                                          col2: regulon activity
## fdr_cut - fdr cutoff for program disease relevance
getDrugTherapyActivity <- function(drug_mat, program, regulon, fdr_cut=0.2) {
  
  ## rename inputs for joining
  colnames(program) <- c("program", "prog_activity")
  colnames(regulon) <- c("regulon", "reg_activity")
  
  ## join regulon and program activity to drug info
  drug_mat %>%
    dplyr::left_join(regulon,by=c("Regulon"="regulon")) %>%
    dplyr::left_join(program, by=c("Program"="program")) %>%
    data.table() ->
    drug_rel
  
  ##drug_rel %>%
  ##  dplyr::filter(Drug=="ACAMPROSATE") ->
  ##  tmp
  
  druggies <- unique(drug_rel$Drug)
  
  ## get list of all drug targets for each drug
  drug_rel[, DrugTargetAll:=paste0(sort(unique(TargetSymbol)), collapse=","),
           by=.(Drug)]
  
  ## get list of all regulons for each drug 
  drug_rel[, DrugRegulonAll:=paste0(sort(unique(Regulon)), collapse=","),
           by=.(Drug)]
  
  ## get list of regulons that are over active
  act_regs <- data.frame(t(sapply(druggies, getRegulonsByActivity, mat=drug_rel, value=1)))
  drug_rel %>%
    dplyr::left_join(act_regs, by=c("Drug"="Drug")) %>%
    dplyr::rename(DrugRegulonOverActive=regulons) ->
    drug_rel
    
  ## get list of regulons that are under active
  act_regs <- data.frame(t(sapply(druggies, getRegulonsByActivity, mat=drug_rel, value=-1)))
  drug_rel %>%
    dplyr::left_join(act_regs, by=c("Drug"="Drug")) %>%
    dplyr::rename(DrugRegulonUnderActive=regulons) ->
    drug_rel 
  
  ## get list of all programs for each drug 
  drug_rel[, DrugProgramAll:=paste0(sort(unique(Program)), collapse=","),
           by=.(Drug)]
  
  ## get list of programs that are over active
  act_progs <- data.frame(t(sapply(druggies, getProgramsByActivity, mat=drug_rel, value=1)))
  drug_rel %>%
    dplyr::left_join(act_progs, by=c("Drug"="Drug")) %>%
    dplyr::rename(DrugProgramOverActive=programs) ->
    drug_rel
  
  ## get list of programs that are under active
  act_progs <- data.frame(t(sapply(druggies, getProgramsByActivity, mat=drug_rel, value=-1)))
  drug_rel %>%
    dplyr::left_join(act_progs, by=c("Drug"="Drug")) %>%
    dplyr::rename(DrugProgramUnderActive=programs) ->
    drug_rel 
  
  ## TODO: check sign of mean of RSMT program direction
  ## calculate regulon disease relevant statistics
  drug_rel[, DrugMeanRMSTRegulon:=meanNA(RMST_UnderMinusOver_reg), by=.(Drug)]
  drug_rel[, DrugRMSTRegulon_MeanDirection:=meanNA(RMSTRegulonDirection), by=.(Drug)]
  drug_rel[, DrugRMSTRegulonMeanPvalue:=meanNA(RMST_reg_pval), by=.(Drug)]
  
  ## calculate program disease relevant statistics
  drug_rel[, DrugRMSTProgram_MeanDirection:=meanNA(RMSTProgramDirection), by=.(Drug)]
  drug_rel[, DrugMeanRMSTProgram:=meanNA(RMST_UnderMinusOver_prog), by=.(Drug)]
  drug_rel[, DrugRMSTProgramMeanPvalue:=meanNA(RMST_prog_pval), by=.(Drug)]
  
  ## calculate drug constrained regulon activity
  drug_rel[, DrugConstrainedRegulonActivity:=meanNA(reg_activity),
           by=.(Drug)]
  
  ## calculate drug constrained program activity
  drug_rel[, DrugConstrainedProgramActivity:=meanNA(prog_activity),
           by=.(Drug)] 
  
  drug_rel %>%
    dplyr::select(-Regulon, -Program, -RegulonRegulator, -Regulon, -RegulatorSymbol, 
                  -RegulatorTarget, -cox_HR_reg, -cox_HR_reg_pvalue, -cox_HR_reg_fdr,
                  -RMST_UnderMinusOver_reg, -RMST_reg_pval, -RMST_reg_fdr, -cox_HR_prog,
                  -cox_HR_prog_pvalue, -cox_HR_prog_fdr, -RMST_UnderMinusOver_prog,
                  -RMST_prog_pval, -RMST_prog_fdr, -Cox_RegulonDirection,
                  -CoxProgramDirection, -RMSTRegulonDirection, -RMSTProgramDirection,
                  -drug_direction, -RegulatorRegulon_Spearman_R, -RegulonMutationAll,
                  -reg_activity, -prog_activity, -TargetSymbol, -DrugTarget) %>%
    distinct() ->
    drug_tmp
  
    ## reorder
    drug_tmp %>%
    dplyr::select(Drug, TargetClass, TargetDisease,                     
                  DrugMechOfAction, effect, DrugTargetAll,                      
                  ot, soc, rr,
                  ot_phase, Cancer_PhaseIV, MM_Phase,
                  DrugRegulonAll, DrugProgramAll, 
                  DrugMeanRMSTProgram, DrugRMSTProgram_MeanDirection, DrugRMSTProgramMeanPvalue,
                  DrugMeanRMSTRegulon, DrugRMSTRegulon_MeanDirection, DrugRMSTRegulonMeanPvalue, 
                  DrugProgramOverActive, DrugProgramUnderActive,
                  DrugConstrainedProgramActivity, 
                  DrugRegulonOverActive, DrugRegulonUnderActive,
                  DrugConstrainedRegulonActivity) %>%
    dplyr::arrange(-DrugConstrainedRegulonActivity) ->
    drug_sum
  
  ##drug_sum %>%
  ##  dplyr::filter(Drug=="DACTOLISIB") ->
  ##  tmp1
  
  return(drug_sum)
}


## Generate drug specific program activity  
## dug_mat - data.frame of drug target info (see drug_information_processing.Rmd
##           output/drugs/drug_map_me_to_patient.csv"
## program - data.frame of program activity for one subject col1: Program label
##                                                          col2: program activity
## regulon - data.frame of regulon activity for one subject col1: regulon label
##                                                          col2: regulon activity
## doRegulon - if TRUE, calculate constrained regulon activity else calculate constrained
##              program activity
getProgramDrugActivity <- function(drug_mat, program, regulon, doRegulon=TRUE) {
  
  ## rename inputs for joining
  colnames(program) <- c("program", "prog_activity")
  colnames(regulon) <- c("regulon", "reg_activity")
  
  ## join regulon and program activity to drug info
  drug_mat %>%
    dplyr::left_join(regulon,by=c("Regulon"="regulon")) %>%
    dplyr::left_join(program, by=c("Program"="program")) ->
    drug_causal_complete
  
  
  
  ## make sure direction of program/regulon and drug type are the same.
  if(doRegulon) {
     ## filter for drug relevance
     drug_causal_complete %>%
        dplyr::filter(RMST_reg_fdr <= fdr_cut) %>%
        data.table() ->
        drug_rel
     drug_rel %<>% dplyr::filter(RMSTRegulonDirection == drug_rel$drug_direction) 
  } else {
     ## filter for drug relevance
     drug_causal_complete %>%
       dplyr::filter(RMST_prog_fdr <= fdr_cut) %>%
       data.table() ->
       drug_rel
     ##sum(drug_rel$RMSTProgramDirection == drug_rel$drug_direction, na.rm=FALSE)
     drug_rel %<>% dplyr::filter(RMSTProgramDirection == drug_rel$drug_direction) 
  }
  
  
  druggies <- unique(drug_rel$drug)
  
  ## get list of all drug targets for each drug
  drug_rel[, DrugTargetAll:=paste0(unique(TargetSymbol), collapse=","),
           by=.(Drug)]
  
  ## get list of all regulons for each drug 
  drug_rel[, DrugRegulonAll:=paste0(unique(Regulon), collapse=","),
           by=.(Drug)]
  
  ## get list of all programs for each drug target
  drug_rel[, DrugProgramAll:=paste0(unique(Program), collapse=","),
           by=.(Drug)]
  
  ## TODO: check sign of mean of RSMT program direction
  ## calculate regulon and program disease relevant statistics
  if(doRegulon) {
     drug_rel[, DrugMeanRMSTRegulon:=meanNA(RMST_UnderMinusOver_reg), by=.(Drug)]
     drug_rel[, DrugRMSTRegulon_MeanDirection:=meanNA(RMSTRegulonDirection), by=.(Drug)]
     drug_rel[, DrugRMSTRegulonMeanPvalue:=meanNA(RMST_reg_pval), by=.(Drug)]
     
     ## calculate drug constrained regulon activity
     drug_rel[, DrugConstrainedRegulonActivity:=meanNA(reg_activity),
              by=.(Drug)]
  } else {
     drug_rel[, DrugRMSTProgram_MeanDirection:=meanNA(RMSTProgramDirection), by=.(Drug)]
     drug_rel[, DrugMeanRMSTProgram:=meanNA(RMST_UnderMinusOver_prog), by=.(Drug)]
     drug_rel[, DrugRMSTProgramMeanPvalue:=meanNA(RMST_prog_pval), by=.(Drug)]
     
     ## calculate drug constrained program activity
     drug_rel[, DrugConstrainedProgramActivity:=meanNA(prog_activity),
              by=.(Drug)] 
  }
  
  
  drug_rel %>%
    dplyr::select(-Regulon, -Program, -RegulonRegulator, -Regulon, -RegulatorSymbol, 
                  -RegulatorTarget, -cox_HR_reg, -cox_HR_reg_pvalue, -cox_HR_reg_fdr,
                  -RMST_UnderMinusOver_reg, -RMST_reg_pval, -RMST_reg_fdr, -cox_HR_prog,
                  -cox_HR_prog_pvalue, -cox_HR_prog_fdr, -RMST_UnderMinusOver_prog,
                  -RMST_prog_pval, -RMST_prog_fdr, -Cox_RegulonDirection,
                  -CoxProgramDirection, -RMSTRegulonDirection, -RMSTProgramDirection,
                  -drug_direction, -RegulatorRegulon_Spearman_R, -RegulonMutationAll,
                  -reg_activity, -prog_activity, -TargetSymbol, -DrugTarget) %>%
    distinct() ->
    drug_sum
  
  return(drug_sum)
}
