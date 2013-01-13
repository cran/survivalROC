\name{survivalROC.C}
\alias{survivalROC.C}
\title{Time-dependent ROC curve estimation from censored survival data}
\description{
  This function creates time-dependent ROC curve from censored survival
  data using the Nearest Neighbor Estimation (NNE) method of Heagerty,
  Lumley and Pepe, 2000}
\usage{
survivalROC.C(Stime,status,marker,predict.time,span)
}
\arguments{
  \item{Stime}{Event time or censoring time for subjects}
  \item{status}{Indicator of status, 1 if death or event, 0 otherwise }
  \item{marker}{Predictor or marker value}
  \item{predict.time}{Time point of the ROC curve}
  \item{span}{Span for the NNE}  
}

\details{ Suppose we have censored survival data along with a baseline
  marker value and
  we want to see how well the marker predicts the survival time for the
  subjects in the dataset. In particular, suppose we have survival times in
  days and  we want to see how well the
  marker predicts the one-year survival (PredictTime=365 days). This
  function returns the unique marker values, sensitivity (True positive
  or TP), (1-specificity) (False positive or FP) and Kaplan-Meier
  survival estimate  corresponding to the
  time point of interest (PredictTime). The (FP,TP) values then can be
  used to construct ROC curve at the time point of interest.
}
\value{Returns a list of the following items:
  \item{cut.values}{unique marker values for calculation of TP and FP}
  \item{TP}{TP corresponding to the cut off in marker}
  \item{FP}{FP corresponding to the cut off in marker}
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
Staltscore4 <- NULL
Mayo.fit4 <- survivalROC.C( Stime = mayo$time,  
      status = mayo$censor,      
      marker = mayo$mayoscore4,     
      predict.time = cutoff,      
      span = 0.25*nobs^(-0.20))
Staltscore4 <- Mayo.fit4$Survival
plot(Mayo.fit4$FP, Mayo.fit4$TP, type = "l",
xlim = c(0,1), ylim = c(0,1),
xlab = paste( "FP \n AUC =",round(Mayo.fit4$AUC,3)),
ylab = "TP",main = "Year = 1" )
abline(0,1)
}
\keyword{survival}
