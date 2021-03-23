estimate.acat <- function(control, start.thresh, p.all, p.shift, p.X, p.fix, p.XRS, p.thresh, p.rnd, q, I, n, response, X, weights, 
              nodes, design.thresh, l.bound.rnd, u.bound.rnd, start.rnd, se){

  
if(control$RS){
  alpha.start <- c(start.thresh,rep(0,p.shift+p.X),1)
  l.bound <- c(rep(-Inf,p.fix-p.XRS),0)
  u.bound <- c(rep(Inf,p.fix-p.XRS),Inf)
  m.opt <- try(optim(par = alpha.start, fn = loglikMO_noRS, gr = scoreMO_noRS,
                     Q = control$Q,  q = q, I = I, n = n, Y = response, X = X, pall = length(alpha.start),
                     pX = p.X, pthresh = p.thresh, pshift = p.shift, 
                     GHweights = weights, GHnodes = nodes, 
                     dthresh = design.thresh, cores = control$cores, lambda = control$lambda,
                     lower = l.bound, upper = u.bound, method="L-BFGS-B"), silent = TRUE)
  
  if(!inherits(m.opt,"try-error") ){
    alpha.start <- c(m.opt$par[1:(p.fix-p.XRS)],rep(0,p.XRS),m.opt$par[length(alpha.start)],0,0.02)
  }else{
    alpha.start <- c(start.thresh, rep(0,p.shift+p.X+p.XRS), start.rnd)
  }
}else{
  alpha.start <- c(start.thresh, rep(0,p.shift+p.X+p.XRS), start.rnd)
}

if(control$RS){
  weights <- weights %*% t(weights)
  control$XforRS <- FALSE
}


## upper and lower bound for all parameters
l.bound <- c(rep(-Inf,p.fix),l.bound.rnd)
u.bound <- c(rep(Inf,p.fix),u.bound.rnd)

## initialize parscale equal for all parameters
par.scale <- rep(1,p.all)

converged <- FALSE
tries.lbfgs <- tries.nlminb <- 0
########################



while(!converged){ 
  if(control$opt.method == "L-BFGS-B"){
    if(control$RS){
      
      m.opt <- try(optim(par = alpha.start, fn = loglikMO, gr = scoreMO,
                         Q = control$Q,  q = q, I = I, n = n, Y = response, X = X, pall = p.all,
                         pX = p.X, pXRS = p.XRS, pthresh = p.thresh, pshift = p.shift, prnd = p.rnd,
                         GHweights = weights, GHnodes = nodes, 
                         scaled = as.numeric(control$rs.scaled), dthresh = design.thresh,
                         cores = control$cores,  lambda = control$lambda,
                         lower = l.bound, upper = u.bound, method="L-BFGS-B",
                         control = list(parscale = par.scale)), silent = TRUE)
    }else{
      m.opt <- try(optim(par = alpha.start, fn = loglikMO_noRS, gr = scoreMO_noRS,
                         Q = control$Q,  q = q, I = I, n = n, Y = response, X = X, pall = p.all,
                         pX = p.X, pthresh = p.thresh, pshift = p.shift, 
                         GHweights = weights, GHnodes = nodes, 
                         dthresh = design.thresh, cores = control$cores, lambda = control$lambda,
                         lower = l.bound, upper = u.bound, method="L-BFGS-B",
                         control = list(parscale = par.scale)), silent = TRUE)
      
    }
    
    if( inherits(m.opt,"try-error") ){
      par.scale <- par.scale*0.1
    }else{
      converged <- TRUE
      loglik <- m.opt$value
    }
    tries.lbfgs <- tries.lbfgs + 1
    
    
    if(inherits(m.opt,"try-error")  & tries.lbfgs == 2){
      par.scale <- rep(1,p.all)
      control$opt.method <- "nlminb"
      warning("Optimization method is automatically switched to nlminb because L-BFGS-B did not converge!")
    }
  }
  
  if(control$opt.method == "nlminb"){
    if(control$RS){
      m.opt <- try(nlminb(start = alpha.start, objective = loglikMO, gradient = scoreMO, 
                          Q = control$Q,  q = q, I = I, n = n, Y = response, X = X, pall = p.all,
                          pX = p.X, pXRS = p.XRS,pthresh = p.thresh, pshift = p.shift,  prnd = p.rnd,
                          GHweights = weights, GHnodes = nodes, 
                          scaled = as.numeric(control$rs.scaled), dthresh = design.thresh,
                          cores = control$cores, lambda = control$lambda, 
                          lower = l.bound, upper = u.bound, scale= par.scale), 
                   silent = TRUE)
    }else{
      m.opt <- try(nlminb(start = alpha.start, objective = loglikMO_noRS, gradient = scoreMO_noRS, 
                          Q = control$Q,  q = q, I = I, n = n, Y = response, X = X, pall = p.all,
                          pX = p.X, pthresh = p.thresh, pshift = p.shift,
                          GHweights = weights, GHnodes = nodes, 
                          dthresh = design.thresh, cores = control$cores,  
                          lambda = control$lambda, lower = l.bound, upper = u.bound, scale= par.scale), 
                   silent = TRUE)
    }
    
    
    if(inherits(m.opt,"try-error")){
      par.scale <- par.scale*0.1
    }else{
      converged <- TRUE
      loglik <- m.opt$objective
    }
    tries.nlminb <- tries.nlminb + 1
    
    if(inherits(m.opt,"try-error") & tries.nlminb == 3){
      stop("Optimization did not converge!")
    }
  }
}

coefs <- m.opt$par

if(se){
  if(control$RS){
    hess <- optimHess(coefs, loglikMO, scoreMO, 
                      Q = control$Q,  q = q, I = I, n = n, Y = response, X = X, pall = p.all,
                      pX = p.X, pXRS = p.XRS,pthresh = p.thresh, pshift = p.shift,  prnd = p.rnd,
                      GHweights = weights, GHnodes = nodes, 
                      scaled = as.numeric(control$rs.scaled), dthresh = design.thresh,
                      cores = control$cores, lambda = control$lambda )
    
  }else{
    hess <- optimHess(coefs, loglikMO_noRS, scoreMO_noRS, 
                      Q = control$Q,  q = q, I = I, n = n, Y = response, X = X, pall = p.all,
                      pX = p.X, pthresh = p.thresh, pshift = p.shift,
                      GHweights = weights, GHnodes = nodes, 
                      dthresh = design.thresh, cores = control$cores, 
                      lambda = control$lambda)
    
  }
  
  hess.inv <- solve(hess)
  se.vec <- sqrt(diag(hess.inv))
}else{
  se.vec <- NULL
}

return(list(coefs = coefs, se.vec = se.vec, loglik = loglik))

}