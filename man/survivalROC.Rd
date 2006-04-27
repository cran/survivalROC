\name{survivalROC}
\alias{survivalROC}
\title{Time-dependent ROC curve estimation from censored survival data}
\description{
  This function creates time-dependent ROC curve from censored survival
  data using the Kaplan-Meier (KM) or Nearest Neighbor Estimation (NNE)
  method of Heagerty, Lumley and Pepe, 2000}
\usage{
survivalROC(Stime, status, marker, entry = NULL, predict.time, cut.values =
NULL, method = "NNE", lambda = NULL, span = NULL, window =
"symmetric")
}
\arguments{
  \item{Stime}{Event time or censoring time for subjects}
  \item{status}{Indicator of status, 1 if death or event, 0 otherwise }
  \item{marker}{Predictor or marker value}
  \item{entry}{Entry time for the subjects}
  \item{predict.time}{Time point of the ROC curve}
  \item{cut.values}{marker values to use as a cut-off for
    calculation of sensitivity and specificity} 
  \item{method}{Method for fitting joint distribution of (x,t), either
    of KM or NNE, the default method is NNE}
  \item{lambda}{smoothing parameter for NNE}
  \item{span}{Span for the NNE, need either lambda or span for NNE}
  \item{window}{window for NNE, either of symmetric or asymmetric}
}
\details{ Suppose we have censored survival data along with a baseline
  marker value and  we want to see how well the marker predicts the 
  survival time for the subjects in the dataset. In particular,
  suppose we have survival times in days and  we want to see how well the
  marker predicts the one-year survival (predict.time=365 days). This
  function roc.KM.calc(), returns the unique marker values, TP (True
  Positive), FP (False Positive),  Kaplan-Meier survival
  estimate  corresponding to the  time point of interest (predict.time)
  and AUC (Area Under (ROC) Curve) at the time point of interest. 
}
\value{Returns a list of the following items:
\item{cut.values}{unique marker values for calculation of TP and FP}
  \item{TP}{True Positive corresponding to the cut offs in x}
  \item{FP}{False Positive corresponding to the cut offs in x}
  \item{predict.time}{time point of interest}
  \item{Survival}{Kaplan-Meier survival estimate at predict.time}
  \item{AUC}{Area Under (ROC) Curve at time predict.time}
}

\references{Heagerty, P.J., Lumley, T., Pepe, M. S. (2000)
  Time-dependent ROC Curves for Censored Survival Data and a Diagnostic
  Marker \emph{Biometrics}, \bold{56}, 337 -- 344} 
\author{Patrick J. Heagerty }
\examples{
data(mayo)
nobs <- NROW(mayo)
cutoff <- 365
  ## MAYOSCORE 4, METHOD = NNE
  Mayo4.1= survivalROC(Stime=mayo$time,  
    status=mayo$censor,      
    marker = mayo$mayoscore4,     
    predict.time = cutoff,span = 0.25*nobs^(-0.20) )
  plot(Mayo4.1$FP, Mayo4.1$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
  xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.1$AUC,3)), 
  ylab="TP",main="Mayoscore 4, Method = NNE \n  Year = 1")
  abline(0,1)

  ## MAYOSCORE 4, METHOD = KM
  Mayo4.2= survivalROC(Stime=mayo$time,  
    status=mayo$censor,      
    marker = mayo$mayoscore4,     
    predict.time =  cutoff, method="KM")
  plot(Mayo4.2$FP, Mayo4.2$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
  xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.2$AUC,3)), 
  ylab="TP",main="Mayoscore 4, Method = KM \n Year = 1")
  abline(0,1)
  }
\keyword{survival}
