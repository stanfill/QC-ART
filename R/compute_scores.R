#' Computes QC-ART scores
#' 
#' @param all_data  a data frame of all of the data to be analyzed
#' @param baseline  vector identifying which rows of `all_data` are the baseline observations
#' @param variables  vector of numbers or column names identifying which columns are the variables to use to compute scores
#' @param prop  how many latent variables should be retained?  proportion of variability explianed by the latent variables
#' 
#' @return vector of QC-ART scores corresponding to the rows of `all_data`
#' @export
#' 
#' @examples
#' 
#' library(lubridate)
#' library(QCART)
#' library(dplyr)
#' library(ggplot2)
#' 
#' #Load the Amidan et al. (2014) dataset
#' data("amidan_et_al")  
#' 
#' #Use `lubridate` to  make wokring with time stamps easier
#' amidan$Acq_Time_Start <- mdy_hm(amidan$Acq_Time_Start)
#' 
#' #Choose one particular instrument to analyze and arrange the rows by time
#' chosen_inst <- "VOrbiETD04"
#' vorb_04 <- filter(amidan,Instrument==chosen_inst)
#' vorb_04 <- arrange(vorb_04,Acq_Time_Start)
#' 
#' #Remove the first 6 rows that occured prior to an instrument cleaning
#' vorb_04 <- vorb_04[-c(1:6),]
#' 
#' #Choose the variables that will be used to compute QC-ART scores
#' chosen_vars <- c("P_2C","MS1_Count","MS2_Count","MS1_2B","MS2_1","MS2_2","MS2_3","MS2_4A","MS2_4B","RT_MS_Q1","RT_MS_Q4","RT_MSMS_Q1","RT_MSMS_Q4","XIC_WideFrac")
#' 
#' #Define which observations will be used as the baseline
#' bline_ob <- 1:50
#' vorb_04$Baseline <- FALSE
#' vorb_04$Baseline[bline_ob] <- TRUE
#' 
#' #Compute the scores for that instrument, baseline set and variables
#' vorb_04$QC_Art <- qcart(all_data = vorb_04, baseline = bline_ob, variables = chosen_vars)
#' 
#' #Plot the scores over time
#' qplot(Acq_Time_Start,QC_Art,data=vorb_04,colour=Baseline)+xlab("Date")+ylab("QC-ART Score")
#' 

qcart <- function(all_data, baseline, variables, prop = 0.95){
  scores <- rep(NA,nrow(all_data))
  
  for(i in 1:nrow(all_data)){
    both <- apply_sign2_mod(all_data[i,],baseobs = all_data[baseline,], vars = variables, explvar = prop)
    scores[i] <- both$Modified
  }
  return(scores)
}

apply_sign2_mod <- function(x,baseobs,vars, mads=FALSE,...){
  fullx <- rbind(x,baseobs)
  res <- sign2_mod(fullx[,vars], return_mads = mads, keep_all=TRUE,...)
  if(mads){
    return(res)
  }
  return(list(Modified=res$x.dist_mod[1], Original=res$x.dist_orig[1]))
}

#This is a modified version of mvoutlier's sign2 function
sign2_mod <- function (x, makeplot = FALSE, explvar = 0.95, qcrit = 0.975, return_mads=FALSE, keep_all=TRUE,...) {
  p = ncol(x)
  n = nrow(x)
  x.mad = apply(x, 2, mad)
  
  
  #I handle mad==0 by using the mean as center instead of median
  if(any(x.mad == 0)){
    redo <- which(x.mad==0)
    if(keep_all){
      for(ii in redo){
        x.mad[ii] <- round(mad(x[,ii],center=mean(x[,ii])),5)
        #x.sc[,ii] <- (x[,ii]-mean(x[,ii]))/x.mad[ii]
        #If that doesn't work, use the mean absolute deviation
        if(x.mad[ii]==0){
          x.mad[ii] <- mean(abs(x[,ii]-mean(x[,ii])))
          #x.sc[,ii] <- (x[,ii]-mean(x[,ii]))/x.mad[ii]
        }
      }
    }else{
      x <- x[,-redo]
      x.mad <- apply(x,2,mad)
    }
  } 
  
  
  #This is how a mad of zero was handeled before
  if (any(x.mad == 0)) 
    stop("More than 50% equal values in one or more variables!")
  
  x.sc <- scale(x, apply(x, 2, median), x.mad)
  
  xs <- x.sc/sqrt(apply(x.sc^2, 1, sum))
  if(return_mads){
    return(apply(xs,2,sd))
  }
  
  svdxs <- svd(xs)
  xs.evec <- svdxs$v
  xs.pc <- x.sc %*% xs.evec
  
  xs.pcscal <- apply(xs.pc, 2, mad)^2
  xs.pcorder <- order(xs.pcscal, decreasing = TRUE)
  p1 <- (1:p)[(cumsum(xs.pcscal[xs.pcorder])/sum(xs.pcscal) > explvar)][1]
  x.pc <- x.sc %*% xs.evec[, xs.pcorder[1:p1]]
  xpc.sc <- scale(x.pc, apply(x.pc, 2, median), apply(x.pc,2, mad))
  xpc.norm <- sqrt(apply(xpc.sc^2, 1, sum))
  xpc.out <- xpc.norm/median(xpc.norm)
  x.dist <- xpc.out * sqrt(qchisq(0.5, p1))
  
  p1new <- min(which(cumsum(svdxs$d)/sum(svdxs$d)>explvar))
  x.pcnew <- x.pc#x.sc %*% xs.evec[, c(1:p1new)]
  xpc.scnew <- scale(x.pcnew, apply(x.pcnew, 2, median), apply(x.pcnew,2, mad))
  #wts <- (svdxs$d[1:p1new])/sum(svdxs$d[1:p1new])
  #for(i in 1:ncol(xpc.scnew)){
  #  xpc.scnew[,i] <- xpc.scnew[,i]*wts[i]
  #}
  xpc.normnew <- sqrt(rowSums(xpc.scnew^2))
  
  
  #const <- sqrt(qchisq(qcrit, p1))
  #wfinal01 <- rep(0, n)
  #wfinal01[x.dist < const] <- 1
  if (makeplot) {
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 2))
    on.exit(par(op))
    plot(x.dist, xlab = "Index", ylab = "Distance", ...)
    abline(h = const)
    plot(wfinal01, xlab = "Index", ylab = "Final 0/1 weight", 
         ylim = c(0, 1), ...)
  }
  list(x.dist_orig = x.dist, x.dist_mod = xpc.normnew)
} 

