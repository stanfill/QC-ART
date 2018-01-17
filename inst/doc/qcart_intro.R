## ----load----------------------------------------------------------------
library(QCART)
data("amidan_et_al")

## ----forease, include=FALSE----------------------------------------------
library(lubridate)
amidan$Acq_Time_Start <- mdy_hm(amidan$Acq_Time_Start)
chosen_inst <- "VOrbiETD04"

## ----subset,message=FALSE------------------------------------------------
library(dplyr)
chosen_inst <- "VOrbiETD04"
vorb_04 <- filter(amidan,Instrument==chosen_inst)
vorb_04 <- vorb_04[-c(1:6),]

## ----plot,messages=FALSE,warning=FALSE,fig.subcap=c("P 2C","MS1_Count"),fig.width=4,fig.height=4,fig.align="center"----
library(lubridate)
amidan$Acq_Time_Start <- mdy_hm(amidan$Acq_Time_Start)
library(ggplot2)
qplot(Acq_Time_Start,P_2C,data=vorb_04,xlab="Date")
qplot(Acq_Time_Start,MS1_Count,data=vorb_04,xlab="Date")

## ----vars----------------------------------------------------------------
chosen_vars <- c("P_2C","MS1_Count","MS2_Count","MS1_2B","MS2_1","MS2_2","MS2_3","MS2_4A","MS2_4B","RT_MS_Q1","RT_MS_Q4","RT_MSMS_Q1","RT_MSMS_Q4","XIC_WideFrac")

## ----run_qcart,fig.align="center",fig.width=6,fig.height=6---------------
vorb_04 <- arrange(vorb_04,Acq_Time_Start)
bline_ob <- 1:50
vorb_04$Baseline <- FALSE
vorb_04$Baseline[bline_ob] <- TRUE
vorb_04$QC_Art <- qcart(all_data = vorb_04, baseline = bline_ob, variables = chosen_vars)
qplot(Acq_Time_Start,QC_Art,data=vorb_04,colour=Baseline)+xlab("Date")+ylab("QC-ART Score")

## ----sthresh, fig.width=8,fig.height=5,fig.align="center",message=FALSE----
s_threshold <- compute_threshold(scores = vorb_04$QC_Art, baseline = vorb_04$Baseline, time_stamps = vorb_04$Acq_Time_Start, type='static')
sthresh_df <- s_threshold$Results

qplot(Time_Stamp,Scores,data=sthresh_df,colour=Baseline)+xlab("Date")+ylab("QC-ART Score")+
  geom_line(aes(Time_Stamp,Static_Threshold),colour=1)

## ----dthresh, fig.width=8,fig.height=5,fig.align="center",message=FALSE----
d_threshold <- compute_threshold(scores = vorb_04$QC_Art, baseline = vorb_04$Baseline, time_stamps = vorb_04$Acq_Time_Start, type = 'dynamic')
dthresh_df <- d_threshold$Results

qplot(Time_Stamp,Scores,data=dthresh_df,colour=Baseline)+xlab("Date")+ylab("QC-ART Score")+
  geom_line(aes(Time_Stamp,Dynamic_Threshold),colour=1)

## ------------------------------------------------------------------------
library(dlm)

#Convert time stamps into time since cohort study began
vorb_04 <- vorb_04%>%arrange(Acq_Time_Start)%>%
  mutate(Time=difftime(Acq_Time_Start,min(Acq_Time_Start),units='hours'))

buildCutpt <- function(theta,xvar){
  dlmModReg(X=xvar, dV=exp(theta[1]), dW=exp(theta[2:3]),addInt=TRUE) 
}

#Estimate parameters based on baseline data only
fit <- dlmMLE(log(vorb_04$QC_Art),parm=rep(0,3),buildCutpt, xvar=matrix(vorb_04$Time,ncol=1))
    
#Now apply the model to the full dataset
modCutpt <- buildCutpt(fit$par, xvar=vorb_04$Time)
    
#Filterd estiamte for full data set
filtQC <- dlmFilter(log(vorb_04$QC_Art),modCutpt)
filtPred <- filtQC$f
filtSD <- residuals(filtQC,sd=TRUE)$sd
    
Dynamic_Threshold <- qlnorm(0.95,meanlog=filtPred,sdlog=filtSD)

## ----fig.align='center'--------------------------------------------------
resids <- residuals(filtQC,sd=FALSE)
qplot(qnorm(ppoints(length(resids))),sort(resids))+geom_abline(aes(intercept=0,slope=1))+
  xlab("Theoretical Quantiles")+ylab("Sample Quantiles")

## ---- fig.align='center',fig.width=5,fig.height=10-----------------------
tsdiag(filtQC,gof.lag = 25)

## ---- fig.align='center',fig.width=5,fig.height=10-----------------------

buildCutpt_MA2 <- function(theta,xvar){
    dlmModReg(X=xvar, dV=1e-7, dW=exp(theta[1:2]),addInt=TRUE)+
      dlmModARMA(ma=exp(theta[3:4]),sigma2=exp(theta[4]))
}

#Estimate parameters based on baseline data only
fit_MA2 <- dlmMLE(log(vorb_04$QC_Art),parm=rep(0,5),buildCutpt_MA2, xvar=matrix(vorb_04$Time,ncol=1))
    
#Now apply the model to the full dataset
modCutpt_MA2 <- buildCutpt_MA2(fit_MA2$par, xvar=vorb_04$Time)
    
#Filterd estiamte for full data set
filtQC_MA2 <- dlmFilter(log(vorb_04$QC_Art),modCutpt_MA2)

tsdiag(filtQC_MA2,gof.lag = 25)


## ----new_dthresh, fig.width=8,fig.height=5,fig.align="center",message=FALSE----

filtPred_MA2 <- filtQC_MA2$f
filtSD_MA2 <- residuals(filtQC_MA2,sd=TRUE)$sd
    
dthresh_df$NewThreshold <- qlnorm(0.95,meanlog=filtPred_MA2,sdlog=filtSD_MA2)

qplot(Time_Stamp,Scores,data=dthresh_df,colour=Baseline)+xlab("Date")+ylab("QC-ART Score")+
  geom_line(aes(Time_Stamp,NewThreshold),colour=1)+ylim(c(0,250))


