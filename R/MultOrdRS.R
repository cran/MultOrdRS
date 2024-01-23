


#' Control Function for multordRS
#'
#' Control function for multordRS, a model for multivariate ordinal responses including response styles
#'
#'
#' @param RS Logical value indicating whether response style should be modelled.
#' @param thresholds.acat Type of parametrization used for thresholds: \code{thresholds = "full"} implies
#' separate estimates of threshold values for each response variable; \code{thresholds = "shift"} implies
#' equal threshold parameter across all response variables modified by shift parameters for
#' each response variable; \code{thresholds = "minimal"} implies equal threshold parameter across all response variables. This option only applies
#' for adjacent categories models (\code{model = "acat"} and is not implemented for cumulative models.)
#' @param XforRS Logical value indicating whether also covariate effects on the
#' response style should be considered. Only relevant if \code{RS = TRUE}.
#' @param opt.method Specifies optimization algorithm used by \code{\link{optim}}, either
#' \code{L-BFGS-B} or \code{nlminb}.
#' @param Q Number of nodes to be used (per dimension) in Gauss-Hermite-Quadrature. If \code{RS = TRUE},
#' Gauss-Hermite-Quadrature is two-dimensional.
#' @param cores Number of cores to be used in parallelized computation.
#' @param lambda Tuning parameter for internal ridge penalty. It is supposed to be set to a small value
#' to stabilize estimates.
#' @return Returns list of control parameters used in \code{\link{multordRS}}.
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}\cr
#' \url{https://orcid.org/0000-0002-0392-1580}
#' @seealso \code{\link{multordRS}} \code{\link{MultOrdRS-package}} \code{\link{plot.MultOrdRS}}
#' @keywords multivariate ordinal response style adjacent categories cumulative
#' @references Schauberger, Gunther and Tutz, Gerhard (2021): Multivariate Ordinal Random Effects Models Including Subject and Group Specific Response Style Effects, 
#' \emph{Statistical Modelling}, \doi{10.1177/1471082X20978034}
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
#' ## Multivariate adjacent categories model, without response style, 
#' ## without explanatory variables
#' m.tense0 <- multordRS(f.tense0, data = tenseness, control = ctrl.multordRS(RS = FALSE))
#' m.tense0
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
#' ################################################################
#' ## Examples from Schauberger and Tutz (2020) on 
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
ctrl.multordRS <-
  function(RS = TRUE,
           thresholds.acat = c("full", "shift", "minimal"),
           XforRS = TRUE,
           opt.method = c("L-BFGS-B", "nlminb"),
           Q = 10,
           cores = 5,
           lambda = 1e-2) {
    
    thresholds <- match.arg(thresholds.acat)
    opt.method = match.arg(opt.method)
    rs.scaled  <-  TRUE
    
    ret.list <-
      list(
        RS = RS,
        rs.scaled = rs.scaled,
        thresholds = thresholds,
        XforRS = XforRS,
        opt.method = opt.method,
        Q = Q,
        cores = cores,
        lambda = lambda
      )
    ret.list
  }



## response function for acat
resp.fun.acat <- function(eta) {
  q <- length(eta)
  eta.help <- matrix(rep(c(0, eta), each = q + 1), ncol = q + 1)
  eta.help[upper.tri(eta.help)] <- 0
  pi <- cumprod(c(1, exp(eta[-q]))) / sum(apply(exp(eta.help), 1, prod))
  pi
}

## create responses for acat from ordinal values
create.resp.acat <- function(Y) {
  c(t(model.matrix( ~ 0 + Y)[, -length(levels(Y))]))
}

## create cumulative response vector
create.resp.cumul <- function(Y) {
  cum.resp <- c()
  q <- length(levels(Y)) - 1
  for (i in 1:length(Y)) {
    cum.resp <- c(cum.resp, as.numeric(as.numeric(Y[i]) <= 1:q))
  }
  cum.resp
}

#' Model Multivariate Ordinal Responses Including Response Styles
#'
#' A model for multivariate ordinal responses. The response is modelled
#' using a mixed model approach that is also capable of the inclusion
#' of response style effects of the respondents.
#'
#' @param formula Formula containing the (multivariate) ordinal response on the left side and the explanatory variables on the right side.
#' @param data Data frame containing the ordinal responses as well as the explanatory variables from the \code{formula}.
#' @param control Control argument for \code{multord()} function. For details see \code{\link{ctrl.multordRS}}.
#' @param se Should standard errors be calculated for the regression coefficients? Default is \code{TRUE}.
#' @param model Specifies, which type of model is used, either the (multivariate) adjacent categories model (\code{model = "acat"}) or the (multivariate) cumulative model (\code{model = "cumulative"}).
#' @return
#' \item{beta.thresh}{Matrix containing all threshold parameters for the respective model.}
#' \item{beta.shift}{Vector containing all shift parameters. Only relevant if \code{model = "acat"} and \code{thresholds.acat = "shift"}.}
#' \item{beta.X}{Vector containing parameter estimates for the location effects of the explanatory variables.}
#' \item{beta.XRS}{Vector containing parameter estimates for the response style effects of the explanatory variables.}
#' \item{Sigma}{Estimate of the variance (or covariance matrix) of the random effects.  The estimate is a matrix if person-specific random response style effects are considered in the model (i.e. if \code{RS = TRUE}).}
#' \item{Y}{Matrix containing the explanatory variables.}
#' \item{X}{Data frame containing the multivariate ordinal response, one row per obeservation, one column per response variable.}
#' \item{se.thresh}{Matrix containing all standard errors of the threshold parameters for the respective model.}
#' \item{se.shift}{Vector containing all standard errors of the shift parameters. Only relevant if \code{model = "acat"} and \code{thresholds.acat = "shift"}.}
#' \item{se.X}{Vector containing standard errors of the parameter estimates for the location effects of the explanatory variables.}
#' \item{se.XRS}{Vector containing standard errors of the parameter estimates for the response style effects of the explanatory variables.}
#' \item{coef.vec}{Complete vector of all parameter estimates (for internal use).}
#' \item{se.vec}{Complete vector of all standard errors (for internal use).}
#' \item{design.values}{Some values of the design matrix (for internal use).}
#' \item{loglik}{(Marginal) Log Likelihood}
#' \item{call}{Function call}
#' \item{df}{Degrees of freedom}
#' \item{control}{Control argument from function call}
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}\cr
#' \url{https://orcid.org/0000-0002-0392-1580}
#' @seealso \code{\link{ctrl.multordRS}} \code{\link{MultOrdRS-package}} \code{\link{plot.MultOrdRS}}
#' @keywords multivariate ordinal response style adjacent categories
#' @references Schauberger, Gunther and Tutz, Gerhard (2021): Multivariate Ordinal Random Effects Models Including Subject and Group Specific Response Style Effects, 
#' \emph{Statistical Modelling}, \doi{10.1177/1471082X20978034}
#' @examples
#' 
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
#' m.tense0 <- multordRS(f.tense0, data = tenseness, control = ctrl.multordRS(RS = FALSE, cores = 2))
#' m.tense0
#' 
#' \donttest{
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
#' m.tense0.cumul <- multordRS(f.tense0, data = tenseness, control = ctrl.multordRS(RS = FALSE), 
#'   model = "cumulative")
#' m.tense0.cumul
#' 
#' ## Multivariate cumulative model, with response style as a random effect, 
#' ## without explanatory variables
#' m.tense1.cumul <- multordRS(f.tense0, data = tenseness, model = "cumulative")
#' m.tense1.cumul
#' 
#' ## Multivariate cumulative model, with response style as a random effect, 
#' ## without explanatory variables for response style BUT for location
#' m.tense2.cumul <- multordRS(f.tense1, data = tenseness, control = ctrl.multordRS(XforRS = FALSE), 
#'   model = "cumulative")
#' m.tense2.cumul
#' 
#' ## Multivariate cumulative model, with response style as a random effect, with 
#' ## explanatory variables for location AND response style
#' m.tense3.cumul <- multordRS(f.tense1, data = tenseness, model = "cumulative")
#' m.tense3.cumul
#' 
#' plot(m.tense3.cumul)
#' 
#' ################################################################
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
multordRS <-
  function(formula,
           data = NULL,
           control = ctrl.multordRS(),
           se = TRUE,
           model = c("acat", "cumulative")) {
    if (is.null(data)) {
      X <- model.matrix(formula)
      Y <- model.response(model.frame(formula))
    } else{
      X <- model.matrix(formula, data = data)
      Y <- model.response(model.frame(formula, data = data))
    }
    
    
    if (colnames(X)[1] == "(Intercept)") {
      X <- X[, -1, drop = FALSE]
    }
    
    if (length(unique(sapply(apply(Y, 2, unique), length))) != 1) {
      stop("All response variable must have equal numbers of response categories!")
    }
    

    ## get basic ordinal model
    model <- match.arg(model)
    
    ## initalize response vector and other parameters
    y.vec <- as.factor(c(t(Y)))
    k <- length(levels(y.vec))
    q <- k - 1
    n <- nrow(Y)
    I <- ncol(Y)
    
    if (model == "acat") {
      ## create design for threshold parameters
      if (control$thresholds == "full") {
        design.thresh <- diag(q * I)
        p.thresh <- q * I
        p.shift <- 0
        start.thresh <- rep(1:q, I) * 0.1
      }
      if (control$thresholds == "shift") {
        design.thresh <-
          cbind(matrix(rep(diag(q), I), ncol = q, byrow = TRUE),
                rbind(matrix(rep(
                  diag(I - 1), each = q
                ), ncol = I - 1),
                matrix(
                  0, ncol = I - 1, nrow = q
                )))
        p.thresh <- q
        p.shift <- I - 1
        start.thresh <- (1:q) * 0.1
      }
      if (control$thresholds == "minimal") {
        design.thresh <- matrix(rep(diag(q), I), ncol = q, byrow = TRUE)
        p.thresh <- q
        p.shift <- 0
        start.thresh <- (1:q) * 0.1
      }
    } else{
      p.thresh <- q * I
      p.shift = 0
      start.thresh <- rep(0.1, q * I)
    }
    
    ## check if response style is intended to be modelled
    ## also, create start values and upper and lower bounds for estimators
    p.rnd <- 1
    start.rnd <- 1
    l.bound.rnd <- 5e-3
    u.bound.rnd <- Inf
    if (control$RS) {
      p.rnd <- 3
      start.rnd <- c(1, 0, 0.02)
      l.bound.rnd <- c(5e-3,-1 + 5e-3, 5e-3)
      u.bound.rnd <- c(Inf, 1 - 5e-3, Inf)
    }
    
    ## check if covariates are specified
    if (is.null(X)) {
      X <- matrix(0, 0, 0)
    }
    p.X <- ncol(X)
    
    ## check if covariates are intended to be used also for response style
    p.XRS <- 0
    if (control$XforRS & control$RS) {
      p.XRS <- p.X
    }
    
    ## initialize starting values
    p.fix <- p.thresh + p.shift + p.X + p.XRS
    p.all <- p.fix + p.rnd
    
    
    ## create correct response
    if (model == "acat") {
      response <- create.resp.acat(y.vec)
    }
    if (model == "cumulative") {
      response <- create.resp.cumul(y.vec)
    }
    
    
    ## get nodes and weights for Gauss-Hermite quadrature
    her_poly <- gauss.quad(control$Q, "hermite")
    nodes <- her_poly$nodes
    weights <- her_poly$weights * dnorm(nodes) * exp(nodes ^ 2)
    
    if (model == "acat") {
      results.estim <-
        estimate.acat(
          control,
          start.thresh,
          p.all,
          p.shift,
          p.X,
          p.fix,
          p.XRS,
          p.thresh,
          p.rnd,
          q,
          I,
          n,
          response,
          X,
          weights,
          nodes,
          design.thresh,
          l.bound.rnd,
          u.bound.rnd,
          start.rnd,
          se
        )
    } else{
      results.estim <-
        estimate.cumul(
          control,
          start.thresh,
          p.all,
          p.shift,
          p.X,
          p.fix,
          p.XRS,
          p.thresh,
          p.rnd,
          q,
          I,
          n,
          response,
          X,
          weights,
          nodes,
          l.bound.rnd,
          u.bound.rnd,
          start.rnd,
          se
        )
    }
    
    coefs <- results.estim$coefs
    if(se){
      se.vec <- results.estim$se.vec
    }else{
      se.vec <- se.thresh <- se.shift <- se.X <- se.XRS <- NA
    }
    loglik <- results.estim$loglik
    
    ########################
    ## extract results and prepare return
    
    beta.thresh <- coefs[1:p.thresh]
    thresh.orig <- matrix(beta.thresh, nrow = I)
    if(se){
      se.thresh <- se.vec[1:p.thresh]
    }
    
    beta.shift <- coefs[(p.thresh + 1):(p.thresh + p.shift)]
    if(se){
    se.shift <- se.vec[(p.thresh + 1):(p.thresh + p.shift)]
    }
    
    beta.X <- coefs[(p.thresh + p.shift + 1):(p.thresh + p.shift + p.X)]
    if(se){
    se.X <- se.vec[(p.thresh + p.shift + 1):(p.thresh + p.shift + p.X)]
    }
    
    beta.XRS <-
      coefs[(p.thresh + p.shift + p.X + 1):(p.thresh + p.shift + p.X + p.XRS)]
    if(se){
    se.XRS <-
      se.vec[(p.thresh + p.shift + p.X + 1):(p.thresh + p.shift + p.X + p.XRS)]
    }
    
    beta.rnd <-
      coefs[(p.thresh + p.shift + p.X + p.XRS + 1):(p.thresh + p.shift + p.X +
                                                      p.XRS + p.rnd)]
    
    
    names.vec <- c()
    
    if (model == "acat") {
      if (p.thresh == q * I) {
        beta.thresh <- matrix(beta.thresh, byrow = TRUE, ncol = q)
        if(se){
        se.thresh <- matrix(se.thresh, byrow = TRUE, ncol = q)
        }
        colnames(beta.thresh) <-paste0("Thresh.", 1:q)
        if(se){
          colnames(se.thresh) <- paste0("Thresh.", 1:q)
        }
        
        rownames(beta.thresh) <- colnames(Y)
        
        if(se){
          rownames(se.thresh) <- colnames(Y)
        }
        
        names.vec <-
          c(names.vec, paste(rep(colnames(Y), each = q), rep(1:q, I), sep = ":"))
      } else{
        names(beta.thresh) <- paste0("Thresh.", 1:q)
        if(se){
          names(se.thresh) <- paste0("Thresh.", 1:q)
        }
        if(se){
          names.vec <- c(names.vec, names(se.thresh))
        }
      }
    } else{
      beta.thresh <- matrix(beta.thresh, ncol = q)
      if(se){
        se.thresh <- matrix(se.thresh, ncol = q)
      }
      colnames(beta.thresh) <- c("gamma", paste0("delta", 2:q))
      if(se){
        colnames(se.thresh) <- c("gamma", paste0("delta", 2:q))
      }
      
      rownames(beta.thresh) <- colnames(Y)
      
      if(se){
        rownames(se.thresh) <- colnames(Y)
      }
      
      names.vec <-
        c(names.vec, paste(rep(colnames(Y), q), rep(c("gamma", paste0("delta", 2:q)), each = I), sep = ":"))
    }
    
    
    
    if (p.shift > 0) {
      names(beta.shift) <- names(se.shift) <- head(colnames(Y), I - 1)
      names.vec <- c(names.vec, names(se.shift))
    } else{
      beta.shift <- se.shift <- NA
    }
    
    if (p.X > 0) {
      names(beta.X) <- colnames(X)
      names.vec <- c(names.vec, names(beta.X))
      if(se){
        names(se.X) <- colnames(X)
      }
    } else{
      beta.X <- se.X <- NA
    }
    
    if (p.XRS > 0) {
      names(beta.XRS) <- colnames(X)
      names.vec <- c(names.vec, names(beta.XRS))
      if(se){
        names(se.XRS) <- colnames(X)
      }
    } else{
      beta.XRS <- se.XRS <- NA
    }
    

    if (p.rnd == 3) {
      Sigma <-
        matrix(c(
          beta.rnd[1],
          beta.rnd[2] * sqrt(beta.rnd[1]) * sqrt(beta.rnd[3]),
          beta.rnd[2] * sqrt(beta.rnd[1]) * sqrt(beta.rnd[3]),
          beta.rnd[3]
        ),
        ncol = 2)
      
      colnames(Sigma) <- rownames(Sigma) <- c("Intercept", "RespStyle")
      names.vec <-
        c(names.vec,  c("Intercept", "Correlation", "RespStyle"))
    } else{
      Sigma <- beta.rnd
      names(Sigma) <- "Intercept"
      names.vec <- c(names.vec,  "Intercept")
    }
    
    
    names(coefs) <- names.vec
    
    if(se){
      names(se.vec) <- names.vec
    }
    
    fun.call <- match.call()
    design.values <- list(
      k = k,
      n = n,
      I = I,
      p.thresh = p.thresh,
      p.shift = p.shift,
      p.X = p.X,
      p.XRS = p.XRS,
      p.rnd = p.rnd,
      data.name = fun.call$data,
      se = se
    )
    
    
    
    ret.list <-
      list(
        beta.thresh = beta.thresh,
        beta.shift = beta.shift,
        beta.X = beta.X,
        beta.XRS = beta.XRS,
        Sigma = Sigma,
        Y = Y,
        X = X,
        control = control,
        se.thresh = se.thresh,
        se.shift = se.shift,
        se.X = se.X,
        se.XRS = se.XRS,
        coef.vec = coefs,
        se.vec = se.vec,
        design.values = design.values,
        loglik = -loglik,
        call = fun.call,
        df = length(coefs),
        control = control,
        thresh.orig = thresh.orig
      )
    
    class(ret.list) <- "MultOrdRS"
    
    return(ret.list)
  }
