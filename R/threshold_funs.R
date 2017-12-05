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
#' @return a data.frame containing the QC-ART scores, baseline indicators, time stamps and threshold values for each qc-art score
#' @export
#' 

compute_threshod <- function(scores, baseline, time_stamps, type = 'static', perc = 0.95){
  
  #Check that the scores, baseline and time_stamps are all the same length
  if(length(scores)!=length(baseline) || length(scores)!=length(time_stamps)){
    stop("The length of scores, baseline and time_stamps all need to match")
  }
  
  #Check `type` argument was correctly specified
  type <- match.arg(tolower(type),c("static","dynamic","both"))
  if(type=="both"){
    type <- c("static","dynamic")
  }
  
  
  #Order the scores and baseline by time_stamps
  qc_data_frame <- data.frame(Scores=scores, Baseline=baseline, T_Stamps=time_stamps)
  qc_data_frame <- qc_data_frame%>%arrange(T_Stamps)%>%mutate(Time=difftime(T_Stamps,min(T_Stamps),units='hours'))
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
      lnorm_sd <- sqrt(anova(lin_mod)[2,3])
    }
    
    qc_data_frame$Static_Threshold <- qlnorm(perc, meanlog = lnorm_mean, sdlog = lnorm_sd)
  }
  
  if("dynamic"%in%type){
    
    ##----------------------------------------##
    ##---- Dynamic threshold computation  ----##
    ##----------------------------------------##
    
    warning("Dynamic thresholds aren't currently available.")    

  }
  
  return(qc_data_frame)
}


#scores <- vorb_04$QC_Art
#baseline <- vorb_04$Baseline
#time_stamps <- vorb_04$Acq_Time_Start

#data_with_thresh <- compute_threshod(scores = scores, baseline = baseline, time_stamps = time_stamps, type='static')

#qplot(Time,Scores,data=data_with_thresh,colour=Baseline)+xlab("Date")+ylab("QC-ART Score")+
# geom_line(aes(Time,Static_Threshold),colour=1)

