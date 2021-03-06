#' Computes threshold values for QC-ART Scores
#' 
#' This function will compute threshold values for QC-ART scores computed using the `qcart` function.
#' 
#' @param scores  the QC-ART scores returned by `qcart`
#' @param baseline  a vector of logical values (TRUE/FALSE) indicating which observations were used as the baseline
#' @param time_stamps  the vector of time stamps assocaited with the QC-ART scores
#' @param type  character specifying the type of threhsold to be returned: "static", "dynamic" or "both"
#' @param percentile what percentile of the model should be use as a cut-off, e.g. 0.9, 0.95, or 0.99
#' 
#' @return  a list containing
#' \tabular{ll}{
#' Results  \tab a `data.frame` with columns for the QC-ART scores, baseline indicators, time stamps and threshold values for each qc-art score\cr
#'  \tab \cr
#' Static_Model  \tab the 'lm' object that defines the final static model fit to the QC-ART scores \cr
#'  \tab \cr
#' Dynamic_Model  \tab the 'lm' object that defines the final static model fit to the QC-ART scores \cr
#'  }
#' @export
#' @examples 
#' #Compute the QC-ART scores.  For step-by-step instructions, see the documentation
#' #for the `qcart` function
#' library(lubridate)
#' library(QCART)
#' library(dplyr)
#' library(ggplot2)
#' data("amidan_et_al")  
#' amidan$Acq_Time_Start <- mdy_hm(amidan$Acq_Time_Start)
#' chosen_inst <- "VOrbiETD04"
#' vorb_04 <- filter(amidan,Instrument==chosen_inst)
#' vorb_04 <- arrange(vorb_04,Acq_Time_Start)
#' vorb_04 <- vorb_04[-c(1:6),]
#' chosen_vars <- c("P_2C","MS1_Count","MS2_Count","MS1_2B","MS2_1","MS2_2","MS2_3","MS2_4A","MS2_4B","RT_MS_Q1","RT_MS_Q4","RT_MSMS_Q1","RT_MSMS_Q4","XIC_WideFrac")
#' bline_ob <- 1:50
#' vorb_04$Baseline <- FALSE
#' vorb_04$Baseline[bline_ob] <- TRUE
#' vorb_04$QC_Art <- qcart(all_data = vorb_04, baseline = bline_ob, variables = chosen_vars)
#' 
#' #Compute a static threshold for the QC-ART scores.  Because there is no significant relationship between
#' #the QC-ART scores and time since the analysis began, the threshold value is the same for all time points.
#' s_threshold <- compute_threshold(scores = vorb_04$QC_Art, baseline = vorb_04$Baseline, time_stamps = vorb_04$Acq_Time_Start, type='static')
#' sthresh_df <- s_threshold$Results
#' qplot(Time_Stamp,Scores,data=sthresh_df,colour=Baseline)+xlab("Date")+ylab("QC-ART Score")+
#'   geom_line(aes(Time_Stamp,Static_Threshold),colour=1)
#'   
#'   
#' #Compute a dynamic threshold for the QC-ART scores.
#' d_threshold <- compute_threshold(scores = vorb_04$QC_Art, baseline = vorb_04$Baseline, time_stamps = vorb_04$Acq_Time_Start, type='dynamic')
#' dthresh_df <- d_threshold$Results
#' qplot(Time_Stamp,Scores,data=dthresh_df,colour=Baseline)+xlab("Date")+ylab("QC-ART Score")+
#'   geom_line(aes(Time_Stamp,Dynamic_Threshold),colour=1)
#' 
#' #Both dynamic and static thresholds can be computed in one call to `compute_threshold` by setting the 
#' #`type` argument to "both" 
#' 
#' both_thresholds <- compute_threshold(scores = vorb_04$QC_Art, baseline = vorb_04$Baseline, time_stamps = vorb_04$Acq_Time_Start, type='both')
#' all_thresh_df <- both_thresholds$Results
#' qplot(Time_Stamp,Scores,data=all_thresh_df,colour=Baseline)+xlab("Date")+ylab("QC-ART Score")+
#'   geom_line(aes(Time_Stamp,Static_Threshold),colour=1,linetype=1)+
#'   geom_line(aes(Time_Stamp,Dynamic_Threshold),colour=1,linetype=2)
#' 

compute_threshold <- function(scores, baseline, time_stamps, type = 'static', perc = 0.95){
  
  #Check that the scores, baseline and time_stamps are all the same length
  if(length(scores)!=length(baseline) || length(scores)!=length(time_stamps)){
    stop("The length of scores, baseline and time_stamps all need to match")
  }
  
  #Check `type` argument was correctly specified
  type <- match.arg(tolower(type),c("static","dynamic","both"))
  if(type=="both"){
    type <- c("static","dynamic")
  }
  
  if("dynamic"%in%type){
    if(!require(dlm)){
      warning("The 'dlm' package is missing and must be installed to compute dynamic thresholds, therefore only static threhsolds will be returned.")
      type <- "static"
    }
  }
  
  
  #Order the scores and baseline by time_stamps
  qc_data_frame <- data.frame(Scores=scores, Baseline=baseline, Time_Stamp=time_stamps)
  qc_data_frame <- qc_data_frame%>%arrange(Time_Stamp)%>%mutate(Time=difftime(Time_Stamp,min(Time_Stamp),units='hours'))
  qc_data_frame$Time <- as.numeric(qc_data_frame$Time)
  
  #filter down to baseline data
  qc_baseline <- filter(qc_data_frame, Baseline)
  
  if("static"%in%type){
    ##----------------------------------------##
    ##----- Static threshold computation  ----##
    ##----------------------------------------##
    
    #Try linear model in time
    lin_mod <- lm(log(Scores)~Time,data=qc_baseline)
    
    #If the linear model doesn't fit the data well, just use the intercept model
    if(anova(lin_mod)[1,5]>0.05){
      lin_mod <- lm(log(Scores)~1,data=qc_baseline)
      lnorm_mean <- unname(lin_mod$coefficients[1])
      lnorm_sd <- sqrt(anova(lin_mod)[1,3])
    }else{
      lnorm_mean <- lin_mod$coefficients[1]+lin_mod$coefficients[2]*qc_data_frame$Time
      lnorm_sd <- sqrt(anova(lin_mod)[2,3])*sqrt(1+1/nrow(qc_baseline))
    }
    
    #The threshold is the 
    qc_data_frame$Static_Threshold <- qlnorm(perc, meanlog = lnorm_mean, sdlog = lnorm_sd)
    static_model <- lin_mod
  }else{
    static_model <- NULL
  }
  
  if("dynamic"%in%type){
    
    ##----------------------------------------##
    ##---- Dynamic threshold computation  ----##
    ##----------------------------------------##
    #Define the DLM model
    buildCutpt <- function(theta,xvar){
      dlmModReg(X=xvar, dV=exp(theta[1]), dW=exp(theta[2:3]),addInt=TRUE) 
    }
    npars <- 3
    
    #Estimate parameters based on baseline data only
    fit <- dlmMLE(log(qc_baseline$Scores),parm=rep(0,npars),buildCutpt, xvar=matrix(qc_baseline$Time,ncol=1))
    
    #Now apply the model to the full dataset
    modCutpt <- buildCutpt(fit$par, xvar=qc_data_frame$Time)
    
    #filter estiamte for full data set
    filtQC <- dlmFilter(log(qc_data_frame$Scores),modCutpt)
    filtPred <- filtQC$f
    filtSD <- residuals(filtQC,sd=TRUE)$sd
    
    qc_data_frame$Dynamic_Threshold <- qlnorm(perc,meanlog=filtPred,sdlog=filtSD)

    dynamic_model <- filtQC
    
  }else{
    dynamic_model <- NULL
  }
  
  return(list(Results=qc_data_frame, Static_Model = static_model, Dynamic_Model = dynamic_model))
}

# ##------ Testing -------##
# library(dlm)
# scores <- vorb_04$QC_Art
# baseline <- vorb_04$Baseline
# time_stamps <- vorb_04$Acq_Time_Start
# 
# qc_data_frame <- compute_threshold(scores = scores, baseline = baseline, time_stamps = time_stamps, type='static')
# qc_baseline <- filter(qc_data_frame,Baseline)
# 
# buildCutpt <- function(theta,xvar){
#   dlmModReg(X=xvar, dV=exp(theta[1]), dW=exp(theta[2:3]),addInt=TRUE) #+ dlmModARMA(ar = c(theta[4]))
# }
# npars <- 3
# 
# #Estimate parameters based on baseline data only
# fit <- dlmMLE(log(qc_baseline$Scores),parm=rep(0,npars),buildCutpt, xvar=matrix(qc_baseline$Time,ncol=1))
# 
# #Now apply the model to the full dataset
# modCutpt <- buildCutpt(fit$par, xvar=qc_data_frame$Time)
# 
# #filter estiamte for full data set
# filtQC <- dlmFilter(log(qc_data_frame$Scores),modCutpt)
# filtPred <- filtQC$f
# filtSD <- residuals(filtQC,sd=TRUE)$sd
# 
# qc_data_frame$Dynamic_Threshold <- pmin(qlnorm(0.9,meanlog=filtPred,sdlog=filtSD),200)
# 
# 
# qplot(Time_Stamp,Scores,data=qc_data_frame,colour=Baseline)+xlab("Date")+ylab("QC-ART Score")+
#   geom_line(aes(Time_Stamp,Dynamic_Threshold),colour=1)




