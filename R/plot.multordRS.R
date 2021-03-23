#' Plot function for MultOrdRS
#' 
#' Plot function for a \code{MultOrdRS} object. Plots show coefficients of the explanatory variables, both with repect to location and response styles.
#' The coefficient pairs are displayed as stars, where the rays represent (1-alpha) confidence intervals.
#' 
#' @usage \method{plot}{MultOrdRS}(x, alpha = 0.05, CIfactor = 0.9, \dots)
#' @param x \code{MultOrdRS} object
#' @param alpha Specifies the confidence level 1-alpha of the confidence interval. 
#' @param CIfactor Argument that helps to control the appearance (the width) of the stars that represent the confidence intervals of both 
#' parameters (location and response style) corresponding to one covariate.
#' @param ... Further plot arguments.
#' @return No return value, called for side effects
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}\cr
#' \url{https://orcid.org/0000-0002-0392-1580}
#' @seealso \code{\link{multordRS}}, \code{\link{ctrl.multordRS}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2021): Multivariate Ordinal Random Effects Models Including Subject and Group Specific Response Style Effects, 
#' \emph{Statistical Modelling}, \url{https://journals.sagepub.com/doi/10.1177/1471082X20978034}
#' @examples
#' \donttest{
#' data(tenseness)
#' 
#' ## create a small subset of the data to speed up calculations
#' set.seed(1860)
#' tenseness <- tenseness[sample(1:nrow(tenseness), 300),]
#' 
#' ## scale all metric variables to get comparable parameter estimates
#' tenseness$Age <- scale(tenseness$Age)
#' tenseness$Income <- scale(tenseness$Income)
#' 
#' ## two formulas, one without and one with explanatory variables (gender and age)
#' f.tense0 <- as.formula(paste("cbind(",paste(names(tenseness)[1:4],collapse=","),") ~ 1"))
#' f.tense1 <- as.formula(paste("cbind(",paste(names(tenseness)[1:4],collapse=","),") ~ Gender + Age"))
#' 
#' 
#' 
#' ####
#' ## Adjacent Categories Models
#' ####
#' 
#' ## Multivariate adjacent categories model, without response style, without explanatory variables
#' m.tense0 <- multordRS(f.tense0, data = tenseness, control = ctrl.multordRS(RS = FALSE))
#' m.tense0
#' 
#' 
#' ## Multivariate adjacent categories model, with response style as a random effect, 
#' ## without explanatory variables
#' m.tense1 <- multordRS(f.tense0, data = tenseness)
#' m.tense1
#' 
#' ## Multivariate adjacent categories model, with response style as a random effect, 
#' ## without explanatory variables for response style BUT for location
#' m.tense2 <- multordRS(f.tense1, data = tenseness, control = ctrl.multordRS(XforRS = FALSE))
#' m.tense2
#' 
#' ## Multivariate adjacent categories model, with response style as a random effect, with 
#' ## explanatory variables for location AND response style
#' m.tense3 <- multordRS(f.tense1, data = tenseness)
#' m.tense3
#' 
#' plot(m.tense3)
#' 
#' 
#' 
#' ####
#' ## Cumulative Models
#' ####
#' 
#' ## Multivariate cumulative model, without response style, without explanatory variables
#' m.tense0.cumul <- multordRS(f.tense0, data = tenseness, 
#'   control = ctrl.multordRS(RS = FALSE), model = "cumulative")
#' m.tense0.cumul
#' 
#' ## Multivariate cumulative model, with response style as a random effect, 
#' ## without explanatory variables
#' m.tense1.cumul <- multordRS(f.tense0, data = tenseness, model = "cumulative")
#' m.tense1.cumul
#' 
#' ## Multivariate cumulative model, with response style as a random effect, 
#' ## without explanatory variables for response style BUT for location
#' m.tense2.cumul <- multordRS(f.tense1, data = tenseness, 
#'   control = ctrl.multordRS(XforRS = FALSE), model = "cumulative")
#' m.tense2.cumul
#' 
#' ## Multivariate cumulative model, with response style as a random effect, with 
#' ## explanatory variables for location AND response style
#' m.tense3.cumul <- multordRS(f.tense1, data = tenseness, model = "cumulative")
#' m.tense3.cumul
#' 
#' plot(m.tense3.cumul)
#' 
#' #################################################################
#' ## Examples from Schauberger and Tutz (2020) 
#' ## Data from the German Longitudinal Election Study (GLES) 2017
#' #################################################################
#' 
#' ####
#' ## Source: German Longitudinal Election Study 2017 
#' ## Rossteutscher et al. 2017, https://doi.org/10.4232/1.12927
#' ####
#' 
#' ## load GLES data
#' data(GLES17)
#' 
#' ## scale data
#' GLES17[,7:11] <- scale(GLES17[,7:11])
#' 
#' ## define formula
#' f.GLES <- as.formula(cbind(RefugeeCrisis, ClimateChange, Terrorism, 
#'                        Globalization, Turkey, NuclearEnergy) ~ 
#'                        Age + Gender + Unemployment + EastWest + Abitur)
#' 
#' ## fit adjacent categories model without and with response style parameters
#' m.GLES0 <- multordRS(f.GLES, data = GLES17, control =  ctrl.multordRS(RS = FALSE, cores = 6))
#' m.GLES <- multordRS(f.GLES, data = GLES17, control =  ctrl.multordRS(cores = 6))
#' 
#' m.GLES0
#' m.GLES
#' 
#' plot(m.GLES, main = "Adjacent categories model")
#' 
#' 
#' ## fit cumulative model without and with response style parameters (takes pretty long!!!)
#' m.GLES20 <- multordRS(f.GLES, data = GLES17,  model="cumul", 
#' control = ctrl.multordRS(opt.method = "nlminb", cores = 6, RS = FALSE))
#' 
#' m.GLES2 <- multordRS(f.GLES, data = GLES17,  model="cumul", 
#' control = ctrl.multordRS(opt.method = "nlminb", cores = 6))
#' 
#' m.GLES20
#' m.GLES2
#' 
#' plot(m.GLES2, main = "Cumulative model")
#' 
#'}
plot.MultOrdRS <- function(x, alpha = 0.05, CIfactor = 0.9, ...){
  
  quant <- qnorm(1-alpha/2)
  
  if(is.na(x$beta.X[1])|is.na(x$beta.XRS[1])){
    stop("Plotting is only possible if covariates are used both for 
         the location and response style effect!")
  }
  
  betaX <- x$beta.X
  betaX.KI <- exp(cbind(betaX-quant*x$se.X,betaX+quant*x$se.X))
  betaX <- exp(betaX)
  
  betaXRS <- x$beta.XRS
  betaXRS.KI <- exp(cbind(betaXRS-quant*x$se.XRS,betaXRS+quant*x$se.XRS))
  betaXRS <- exp(betaXRS)
  
  
  
  plot(betaX,betaXRS,pch=16,xlim=range(c(1,betaX.KI)),ylim=range(c(1,betaXRS.KI)),
       xlab=expression(exp(gamma)),ylab=expression(exp(alpha)), ...)
  
  p.X <- length(betaX)
  
  label.x <- label.y <- c()

  for(i in 1:p.X){
    
    x <- c(betaX.KI[i,1],betaX.KI[i,1]+(betaX[i]-betaX.KI[i,1])*(CIfactor),betaX[i],betaX[i]+(betaX[i]-betaX.KI[i,1])*(1-CIfactor),
           betaX.KI[i,2],betaX[i]+(betaX[i]-betaX.KI[i,1])*(1-CIfactor),betaX[i],betaX.KI[i,1]+(betaX[i]-betaX.KI[i,1])*(CIfactor),
           betaX.KI[i,1])
    
    y <- c(betaXRS[i],betaXRS.KI[i,1]+(betaXRS[i]-betaXRS.KI[i,1])*(CIfactor),betaXRS.KI[i,1],betaXRS.KI[i,1]+(betaXRS[i]-betaXRS.KI[i,1])*(CIfactor),betaXRS[i],
           betaXRS[i]+(betaXRS[i]-betaXRS.KI[i,1])*(1-CIfactor),betaXRS.KI[i,2],betaXRS[i]+(betaXRS[i]-betaXRS.KI[i,1])*(1-CIfactor),betaXRS[i])
    
    polygon(x,y,col=grey(0.9))
    label.x <- c(label.x,x[6])
    label.y <- c(label.y,y[6])
  }
  points(betaX,betaXRS,pch=16)
  abline(h=1,lty=2,lwd=2,col="gray")
  abline(v=1,lty=2,lwd=2,col="gray")
  
  text(label.x,label.y,labels=names(betaX),adj=c(-0.1,-0.1))
}
