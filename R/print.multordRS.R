print.MultOrdRS <- function(x, ...){
  cat("Output of multordRS estimation:","\n")
  
  cat("---","\n")
  
  cat("\n")

  cat("Call:\n")
  
  print(x$call,...)
  
  
  p.rnd <- x$design.values$p.rnd
  p.thresh <- x$design.values$p.thresh + x$design.values$p.shift
  p.X <- x$design.values$p.X
  p.XRS <- x$design.values$p.XRS
  p <- x$df
  p.fix <- p - p.rnd
  
  coef.mat <- cbind(x$coef.vec, x$se.vec)
  coef.mat <- cbind(coef.mat, coef.mat[,1]/coef.mat[,2])
  coef.mat <- cbind(coef.mat, 2*pnorm(-abs(coef.mat[, 3])))
  colnames(coef.mat) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
  
  coef.thresh <- coef.mat[1:p.thresh,]

  

  cat("\n")

  if(x$design.values$p.shift>0){
    cat("Threshold and shift coefficients:\n")
    print(x$beta.thresh, ...)
    print(x$beta.shift, ...)
    # printCoefmat(coef.thresh, ...)
  }else{
    cat("Threshold coefficients:\n")
    # printCoefmat(coef.thresh, ...)
    print(x$beta.thresh, ...)
  }
  
  cat("\n")
  
  if(p.X>0){
    coef.X <- coef.mat[(p.thresh+1):(p.thresh+p.X),]
    
    cat("Location effects:\n")
    printCoefmat(coef.X, ...)
    
    cat("\n")
  }
  
  if(p.XRS>0){
    coef.XRS <- coef.mat[(p.thresh+p.X+1):(p.thresh+p.X+p.XRS),]
  
  cat("Response style effects:\n")
  printCoefmat(coef.XRS, ...)
  
  cat("\n")
  }
  
  cat("Estimates of random effects (co-)variance Sigma","\n")
  print(x$Sigma, ...)
  
  cat("\n")
  
  cat("Log-Likelihood:",x$loglik, "with", x$df,"degrees of freedom\n")
  
  invisible(coef.mat)
}





