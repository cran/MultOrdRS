
// [[Rcpp::depends(RcppArmadillo)]]
// #define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include <vector>
#include <algorithm>
#include <assert.h>
#include <omp.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(openmp)]]

// square root of a matrix
mat matsqrt2(mat A){
  vec eigval;
  mat eigvec;
  
  eig_sym(eigval, eigvec, A);
  
  colvec d = (eigval+abs(eigval))*0.5;
  colvec d2 = sqrt(d);
  mat B = eigvec*diagmat(d2)*trans(eigvec);
  
  return B;
}


// response function for acat
vec respFun_acat(vec eta){
  
  int q = eta.n_rows;
  mat eta_help1 = ones(q+1)*trans(join_cols(zeros(1),eta));
  eta_help1 = eta_help1 % trimatl(ones(eta_help1.n_rows,eta_help1.n_cols));
  vec pi_vec = ones(q)/as_scalar(sum(prod(exp(eta_help1),1)));
  for(int k=1; k<q ;k++){
    pi_vec(k) =pi_vec(k-1)*exp(eta(k-1));
  }
  pi_vec = (pi_vec-0.5)*0.9999999+0.5;
  return pi_vec;
}

// response function for cumulative
vec respFun_cumul(vec eta){
  vec pi_vec = exp(eta)/(1+exp(eta));
  return pi_vec;
}

vec respFun_cumul2(vec eta){
  eta(find(eta>20)) = ones(size(find(eta>20)))*20;
  eta(find(eta<-20)) = -ones(size(find(eta<-20)))*20;
  vec pi_vec = exp(eta)/(1+exp(eta));
  return pi_vec;
}

// create inverse sigma for acat
mat createSigmaInv_acat(vec mu){
  mat Sigma = diagmat(mu) - mu * trans(mu);
  mat SigmaInv;
  try{
    SigmaInv = inv(Sigma);
  }
  catch(...){
    SigmaInv = pinv(Sigma);
  }
  return SigmaInv;
}

// create inverse sigma for cumulative
mat createSigmaInv_cumul(vec mu){
  mat Sigma = mu*trans(1-mu);
  
  Sigma = symmatu(Sigma);
  // Sigma = Sigma + diagmat(ones(Sigma.n_cols)*0.00001);
  
  mat SigmaInv = inv(Sigma); 
  return SigmaInv;
}


// create derivative matrix for acat
mat createD_acat(vec mu){
  int q = mu.n_rows;

  mat D2 = zeros(q,q) - diagmat(1/mu);

  D2(span(1,q-1),span(0,q-2)) = D2(span(1,q-1),span(0,q-2)) + diagmat(1/mu(span(1,q-1)));
  D2(span::all,q-1) = -ones(q,1)/(1-as_scalar(sum(mu)));
  D2(q-1,q-1) = as_scalar(-(1-sum(mu(span(0,q-2))))/((1-sum(mu))*mu(q-1)));
  // D2 = D2 + diagmat(ones(D2.n_cols)*0.00001);
  mat D;
  try{
    D = inv(D2);
  }
  catch(...){
    D = pinv(D2);
  }
  
  return D;
}


// create derivative matrix for cumulative
mat createD_cumul(vec eta){
  vec D = exp(eta)/(1+exp(eta))%(1-exp(eta)/(1+exp(eta)));
  
  return diagmat(D);
}
  
  // [[Rcpp::export]]
  arma::vec scoreMO(arma::vec alpha,
           arma::vec Y,
           arma::mat X,
           int Q,
           int q,
           int n,
           int I,
           int pall,
           int pX,
           int pXRS,
           int pthresh,
           int pshift,
           int prnd,
           arma::mat GHweights,
           arma::vec GHnodes,
           int scaled,
           arma::mat dthresh,
           double cores,
           double lambda) { 
             

      // initialize score vector
      vec s = zeros(pall);
    
    vec P2 = 2*alpha*lambda;

      // initialize cij matrix for all persons and all knots, 
      // will be needed for normalization per person afterwards
      // in contrast, for implementation of PCM this was calculated BEFORE looping through persons
      // should be faster this way
      mat cij_mat = zeros(n,Q*Q);
      
      // initialize matrix containing score contributions per person and per parameter
      // has to be kept because one has to use cij_mat before one can sum over persons
      mat help_mat = zeros(pall,n);

      // initialize design for one person, without random effects
      mat Z = join_rows(dthresh,zeros(q*I,prnd+pXRS+pX));
      mat etai = Z*alpha;
      // current linear predictor without random effects, per item and category
      etai = reshape(etai,q,I);

      // create design matrix for both random effects
      vec resp_style_help = zeros(q);

      for(int r=0; r<q ;r++){
        resp_style_help(r) = double(q-1)/2 - r;
      }
      
      if(scaled==0){
        resp_style_help = sign(resp_style_help);
      }
      mat design_rnd = join_rows(ones(q), resp_style_help);

      // create sigma matrix from current parameters
      double co_var = alpha(pall-2)*sqrt(alpha(pall-1))*sqrt(alpha(pall-3));
      mat sigma = zeros(2,2);        
      sigma(0,0) = alpha(pall-3);
      sigma(1,0) = co_var;
      sigma(0,1) = co_var;
      sigma(1,1) = alpha(pall-1);

      
      mat sigma12 = chol(sigma);

      vec etaTheta = zeros(n);
      vec etaGamma = zeros(n);
      if(pX>0){
        etaTheta = X * alpha(span(pthresh+pshift,pthresh+pshift+pX-1));
      }
      if(pXRS>0){
        etaGamma = X * alpha(span(pthresh+pshift+pX,pthresh+pshift+pX+pXRS-1));
      }
      mat etaBoth = join_rows(etaTheta,etaGamma);
      
  // initialize number of different threads
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif

#pragma omp parallel for
    // loop through persons
    for(int i=0; i<n ;i++){ 

      // initialize a running vector over all Q*Q knots
      int pos_knot = 0;
      
      int pos = i*q*I;
      // get response of person i
      vec yi = Y(span(pos,pos+q*I-1));

      mat Xihelp;
      if(pX>0){
      Xihelp = X(i,span::all);
      if(q>1){
        for(int f=1; f<q ;f++){ 
          Xihelp = join_cols(Xihelp,X(i,span::all));
      }
      }

      if(pXRS>0){
      mat Xihelp2 = Xihelp;
      
        for(int f=0; f<pX ;f++){ 
          Xihelp2(span::all,f) = Xihelp2(span::all,f)%resp_style_help;
        }
      Xihelp = join_rows(Xihelp,Xihelp2);
      }
      }

      // initialize design for person i
      mat Zij = Z;

      
      // loop through all knots, for both random effects
      for(int j=0; j<Q ;j++){  
        for(int jj=0; jj<Q ;jj++){
          // initialize derivatives, inverse sigma and mu for all items and current knots
          mat D_i = zeros(q*I,q*I);
          mat SigmaInv_i = zeros(q*I,q*I);
          vec mu_i = zeros(q*I);

          // initialize current random effects or knots
          vec rnd_act = zeros(2);
          rnd_act(0) = GHnodes(j);
          rnd_act(1) = GHnodes(jj);
          
          rnd_act = sigma12*rnd_act + trans(etaBoth(i,span::all));

          // initialize current cij value with weight and prob of current knot-combination
          double cij_help = GHweights(j,jj);
            
          // loop through all items  
          for(int k=0; k<I ;k++){  
            
            // response of person i for current item k
            vec yi_k = yi(span(k*q,q+k*q-1));
            
              yi_k = join_cols(yi_k,1-sum(yi_k,0));

            // create current eta and mu
            vec eta_k = etai(span::all,k) + design_rnd*rnd_act;
            vec mu_k = respFun_acat(eta_k);
            
            mu_i(span(k*q,q+k*q-1)) = mu_k;

            // create colmuns for random effects in design matrix, 
            // one column per parameter, therefore three columns
            mat design_act2 = zeros(q,3);
            
            design_act2(span::all,0) = (design_rnd(span::all,0)*GHnodes(j))/(2*sqrt(sigma(0,0)));
            design_act2(span::all,1) = GHnodes(j)*sqrt(sigma(1,1))*design_rnd(span::all,1) - (design_rnd(span::all,1)*GHnodes(jj)*sigma(1,1)*alpha(pall-2))/sqrt(sigma(1,1)-sigma(1,1)*pow(alpha(pall-2),2));
            design_act2(span::all,2) = design_rnd(span::all,1)*GHnodes(j)*alpha(pall-2)/(2*sqrt(sigma(1,1))) + (GHnodes(jj)*(1-pow(alpha(pall-2),2))*design_rnd(span::all,1))/(2*sqrt(sigma(1,1)-sigma(1,1)*pow(alpha(pall-2),2)));
            
            
            // update the respective part of design matrix
            if(pX>0){
              Zij(span(k*q,q+k*q-1),span(pthresh+pshift,pall-1)) = join_rows(Xihelp, design_act2);
            }else{
              Zij(span(k*q,q+k*q-1),span(pthresh+pshift,pall-1)) = design_act2;
            }

            // update derivative and inverse sigma for current item
              D_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createD_acat(mu_k);
              SigmaInv_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createSigmaInv_acat(mu_k);


            // update current cij by probability of current response
              mu_k = join_cols(mu_k,1-sum(mu_k,0));

          
            cij_help = cij_help * prod(mu_k % yi_k - (yi_k-1));
          }
          // after looping through all items, update cij_mat
          cij_mat(i,pos_knot) = cij_help;
       
          // update all contribution of current knot-combination and person i to column of person i
          help_mat(span::all,i) = help_mat(span::all,i) + (trans(Zij)*D_i*SigmaInv_i*(yi-mu_i))*cij_help;
          
          pos_knot = pos_knot + 1;
        }}
    }
    
      // normalization value per person
      vec cij_norm = sum(cij_mat,1);
      
      // normalize row of each parameter by cij_norm
      for(int e=0; e<pall ;e++){  
        help_mat(e,span::all) = help_mat(e,span::all)%trans(1/cij_norm);
      }

      // sum score contributions over all persons per covariate
      s = -sum(help_mat,1)+P2;
      
      return s;
           }
           
 // [[Rcpp::export]]
  double loglikMO(arma::vec alpha,
           arma::vec Y,
           arma::mat X,
           int Q,
           int q,
           int n,
           int I,
           int pall,
           int pX,
           int pXRS,
           int pthresh,
           int pshift,
           int prnd,
           arma::mat GHweights,
           arma::vec GHnodes,
           int scaled,
           arma::mat dthresh,
           int cores,
           double lambda) { 

      // initialize loglikelihood       
      double f = 0;   

// current value of L2 penalty
double P2 = accu(alpha%alpha)*lambda;

      // initialize design for one person, without random effects
      mat Z = join_rows(dthresh,zeros(q*I,prnd+pX+pXRS));
      mat etai = Z*alpha;
      // current linear predictor without random effects, per item and category
      etai = reshape(etai,q,I);

      // create design matrix for both random effects
      vec resp_style_help = zeros(q);

        for(int r=0; r<q ;r++){
          resp_style_help(r) = double(q-1)/2 - r;
        }

      if(scaled==0){
        resp_style_help = sign(resp_style_help);
      }
      mat design_rnd = join_rows(ones(q), resp_style_help);

      // create sigma matrix from current parameters
      double co_var = alpha(pall-2)*sqrt(alpha(pall-1))*sqrt(alpha(pall-3));
      mat sigma = zeros(2,2);        
      sigma(0,0) = alpha(pall-3);
      sigma(1,0) = co_var;
      sigma(0,1) = co_var;
      sigma(1,1) = alpha(pall-1);
      
    
    mat sigma12 = chol(sigma);
    
      vec etaTheta = zeros(n);
      vec etaGamma = zeros(n);
      if(pX>0){
        etaTheta = X * alpha(span(pthresh+pshift,pthresh+pshift+pX-1));
      }
      if(pXRS>0){
        etaGamma = X * alpha(span(pthresh+pshift+pX,pthresh+pshift+pX+pXRS-1));
      }
      mat etaBoth = join_rows(etaTheta,etaGamma);
      
  // initialize number of different threads
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif

      // loop through all persons
       #pragma omp parallel for reduction(-:f)
      for(int i=0; i<n ;i++){
int pos = i*q*I;
        // get response of person i
        vec yi = Y(span(pos,pos+q*I-1));
 
        mat prods_i = ones(Q,Q);
          // loop through all knots, for both random effects
          for(int j=0; j<Q ;j++){ 
            for(int jj=0; jj<Q ;jj++){ 
              // initialize current random effects or knots
              vec rnd_act = zeros(2);
              rnd_act(0) = GHnodes(j);
              rnd_act(1) = GHnodes(jj);

              rnd_act = sigma12*rnd_act + trans(etaBoth(i,span::all));

              
              // loop through all items  
              for(int k=0; k<I ;k++){  
                // response of person i for current item k
                vec yi_k = yi(span(k*q,q+k*q-1));
                
                yi_k = join_cols(yi_k,1-sum(yi_k,0));
         
                // get eta and mu of person i for current item k
                vec eta_k = etai(span::all,k) + design_rnd*rnd_act;

                vec mu_k = respFun_acat(eta_k);
                mu_k = join_cols(mu_k,1-sum(mu_k,0));

                // create prob for current item, update prob matrix for respective knots
                vec help_pow = mu_k % yi_k - (yi_k-1);

                prods_i(j,jj) = prods_i(j,jj)*prod(help_pow);
              }
            }
          }
         f  -= log(accu(prods_i%GHweights));
         // accumulate all likelihood contributions, weights by respective weights and probs of knots
      }

  return (f+P2);
}


// [[Rcpp::export]]
double loglikMO_noRS(arma::vec alpha,
                arma::vec Y,
                arma::mat X,
                int Q,
                int q,
                int n,
                int I,
                int pall,
                int pX,
                int pthresh,
                int pshift,
                arma::vec GHweights,
                arma::vec GHnodes,
                arma::mat dthresh,
                int cores,
                double lambda) { 
  
  // initialize loglikelihood       
  double f = 0;   
  
  // current value of L2 penalty
  double P2 = accu(alpha%alpha)*lambda;
  
  // initialize design for one person, without random effects
  mat Z = join_rows(dthresh,zeros(q*I,1+pX));
  mat etai = Z*alpha;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  
  vec etaTheta = zeros(n);
  if(pX>0){
    etaTheta = X * alpha(span(pthresh+pshift,pthresh+pshift+pX-1));
  }

  
  // initialize number of different threads
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif
  
  // loop through all persons
#pragma omp parallel for reduction(-:f)
  for(int i=0; i<n ;i++){
    int pos = i*q*I;
    // get response of person i
    vec yi = Y(span(pos,pos+q*I-1));
    
    mat prods_i = ones(Q);
    // loop through all knots, for both random effects
    for(int j=0; j<Q ;j++){ 
        
        // loop through all items  
        for(int k=0; k<I ;k++){  
          // response of person i for current item k
          vec yi_k = yi(span(k*q,q+k*q-1));
            yi_k = join_cols(yi_k,1-sum(yi_k,0));


          // get eta and mu of person i for current item k
          vec eta_k = etai(span::all,k) + sqrt(alpha(pall-1))*GHnodes(j) + etaTheta(i);
          vec mu_k = respFun_acat(eta_k);
            mu_k = join_cols(mu_k,1-sum(mu_k,0));

          // create prob for current item, update prob matrix for respective knots
          vec help_pow = mu_k % yi_k - (yi_k-1);

          prods_i(j) = prods_i(j)*prod(help_pow);
        }
      }
    f  -= log(accu(prods_i%GHweights));
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  
  return (f+P2);
}


// [[Rcpp::export]]
arma::vec scoreMO_noRS(arma::vec alpha,
                  arma::vec Y,
                  arma::mat X,
                  int Q,
                  int q,
                  int n,
                  int I,
                  int pall,
                  int pX,
                  int pthresh,
                  int pshift,
                  arma::vec GHweights,
                  arma::vec GHnodes,
                  arma::mat dthresh,
                  double cores,
                  double lambda) { 
  
  
  // initialize score vector
  vec s = zeros(pall);
  
  vec P2 = 2*alpha*lambda;
  
  // initialize cij matrix for all persons and all knots, 
  // will be needed for normalization per person afterwards
  // in contrast, for implementation of PCM this was calculated BEFORE looping through persons
  // should be faster this way
  mat cij_mat = zeros(n,Q);
  
  // initialize matrix containing score contributions per person and per parameter
  // has to be kept because one has to use cij_mat before one can sum over persons
  mat help_mat = zeros(pall,n);
  
  // initialize design for one person, without random effects
  mat Z = join_rows(dthresh,zeros(q*I,1+pX));
  mat etai = Z*alpha;
  // current linear predictor without random effects, per item and category
  etai = reshape(etai,q,I);
  
  
  vec etaTheta = zeros(n);
  if(pX>0){
    etaTheta = X * alpha(span(pthresh+pshift,pthresh+pshift+pX-1));
  }
  
  // initialize number of different threads
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif
  
  
#pragma omp parallel for
  // loop through persons
  for(int i=0; i<n ;i++){ 
    
    // initialize a running vector over all Q*Q knots
    int pos_knot = 0;
    
    int pos = i*q*I;
    // get response of person i
    vec yi = Y(span(pos,pos+q*I-1));
    
    mat Xihelp;
    if(pX>0){
      Xihelp = X(i,span::all);
      if(q>1){
        for(int f=1; f<q ;f++){ 
          Xihelp = join_cols(Xihelp,X(i,span::all));
        }
      }
    }
    
    // initialize design for person i
    mat Zij = Z;
    
    
    // loop through all knots, for both random effects
    for(int j=0; j<Q ;j++){  

        // initialize derivatives, inverse sigma and mu for all items and current knots
        mat D_i = zeros(q*I,q*I);
        mat SigmaInv_i = zeros(q*I,q*I);
        vec mu_i = zeros(q*I);
        
        // initialize current cij value with weight and prob of current knot-combination
        double cij_help = GHweights(j);
        
        // loop through all items  
        for(int k=0; k<I ;k++){  
          
          // response of person i for current item k
          vec yi_k = yi(span(k*q,q+k*q-1));
            yi_k = join_cols(yi_k,1-sum(yi_k,0));

          
          // create current eta and mu
          vec eta_k = etai(span::all,k) + sqrt(alpha(pall-1))*GHnodes(j) + etaTheta(i);

          vec mu_k = respFun_acat(eta_k);

          mu_i(span(k*q,q+k*q-1)) = mu_k;
          
          // create colmuns for random effects in design matrix, 
          // one column per parameter, therefore three columns
          vec design_act2 = (ones(q)*GHnodes(j))/(2*sqrt(alpha(pall-1)));
      
          
          // update the respective part of design matrix
          if(pX>0){
            Zij(span(k*q,q+k*q-1),span(pthresh+pshift,pall-1)) = join_rows(Xihelp, design_act2);
          }else{
            Zij(span(k*q,q+k*q-1),span(pthresh+pshift,pall-1)) = design_act2;
          }
          
          // update derivative and inverse sigma for current item
            D_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createD_acat(mu_k);
            SigmaInv_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createSigmaInv_acat(mu_k);

          // update current cij by probability of current response
            mu_k = join_cols(mu_k,1-sum(mu_k,0));

          
          cij_help = cij_help * prod(mu_k % yi_k - (yi_k-1));
        }
        // after looping through all items, update cij_mat
        cij_mat(i,pos_knot) = cij_help;
        
        // update all contribution of current knot-combination and person i to column of person i
        help_mat(span::all,i) = help_mat(span::all,i) + (trans(Zij)*D_i*SigmaInv_i*(yi-mu_i))*cij_help;
        
        pos_knot = pos_knot + 1;
      }
  }
  
  // normalization value per person
  vec cij_norm = sum(cij_mat,1);
  
  // normalize row of each parameter by cij_norm
  for(int e=0; e<pall ;e++){  
    help_mat(e,span::all) = help_mat(e,span::all)%trans(1/cij_norm);
  }
  
  // sum score contributions over all persons per covariate
  s = -sum(help_mat,1)+P2;
  
  return s;
}

double loglikMO_cumul2(arma::vec alpha,
                       arma::vec Y,
                       arma::mat X,
                       int Q,
                       int q,
                       int n,
                       int I,
                       int pall,
                       int pX,
                       int pXRS,
                       int pthresh,
                       int pshift,
                       int prnd,
                       arma::mat GHweights,
                       arma::vec GHnodes,
                       int scaled,
                       int cores,
                       double lambda) { 
  
  // initialize loglikelihood       
  
  vec f = zeros(n);
  
  
  // current value of L2 penalty
  double P2 = accu(alpha%alpha)*lambda;
  
  // create sigma matrix from current parameters
  double co_var = alpha(pall-2)*sqrt(alpha(pall-1))*sqrt(alpha(pall-3));
  mat sigma = zeros(2,2);        
  sigma(0,0) = alpha(pall-3);
  sigma(1,0) = co_var;
  sigma(0,1) = co_var;
  sigma(1,1) = alpha(pall-1);
  
  
  mat sigma12 = chol(sigma);
  
  vec etaTheta = zeros(n);
  vec etaGamma = zeros(n);
  if(pX>0){
    etaTheta = X * alpha(span(pthresh+pshift,pthresh+pshift+pX-1));
  }
  if(pXRS>0){
    etaGamma = X * alpha(span(pthresh+pshift+pX,pthresh+pshift+pX+pXRS-1));
  }
  
  vec intercepts = alpha(span(0,I-1));
  
  // initialize design for all thresholds
  mat threshs = zeros(q,I);
  mat expvals = trans(reshape(exp(alpha(span(I,pthresh-1))),I,q-1));
  int m;
  
  // falls Rest Null ist k ungerade
  if(q%2==0){
    m = (q+2)/2;
    
    for(int i=0; i < q; i++){
      
      if(i==(m-2)){
        threshs(i,span::all) = -expvals(m-2,span::all)/2;
      }
      if(i==(m-1)){
        threshs(i,span::all) = expvals(m-2,span::all)/2;
      }
      if(i<(m-2)){
        threshs(i,span::all) =  -expvals(m-2,span::all)/2 - sum(expvals(span(i, m-3),span::all), 0);
      }
      
      if(i>(m-1)){
        threshs(i,span::all) = expvals(m-2,span::all)/2 + sum(expvals(span(m-1,i-1),span::all), 0);
      }
    }
  }else{
    m = (q+1)/2;
    for(int i=0; i < q; i++){
      // threshs(i,span::all) = intercepts(span::all);
      if(i<(m-1)){
        threshs(i,span::all) = -sum(expvals(span(i, m-2),span::all), 0);
      }
      if(i>(m-1)){
        threshs(i,span::all) =  sum(expvals(span(m-1,i-1),span::all), 0);
      }
      
    }
  }
  
  // initialize number of different threads
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif
  
  // loop through all persons
  // #pragma omp parallel for reduction(-:f)
  
  int i; int pos; vec yi; mat prods_i; int j; int jj; vec rnd_act; int k; vec yi_k; 
  vec eta_k; vec mu_k; vec help_pow;
  
  
#ifdef _OPENMP
#pragma omp parallel for private(i, pos, yi, prods_i, j, jj, rnd_act, k,  yi_k, eta_k, mu_k, help_pow) shared(f)
#endif
  for(i=0; i<n ;i++){
    pos = i*q*I;
    // get response of person i
    yi = Y(span(pos,pos+q*I-1));
    
    prods_i = ones(Q,Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 
      for(jj=0; jj<Q ;jj++){ 
        // initialize current random effects or knots
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        
        rnd_act = sigma12*rnd_act;
        
        rnd_act(1) = exp(rnd_act(1) + etaGamma(i)); 
        rnd_act(0) = rnd_act(0) + etaTheta(i); 
        
        // loop through all items  
        for(k=0; k<I ;k++){  
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          
          yi_k = diff(join_cols(join_cols(zeros<vec>(1),yi_k),ones<vec>(1)));
          
          // get eta and mu of person i for current item k
          
          eta_k = intercepts(k) + threshs(span::all,k)*rnd_act(1) + rnd_act(0);
          
          
          mu_k = respFun_cumul2(eta_k);
          
          mu_k = diff(join_cols(join_cols(zeros<vec>(1),mu_k),ones<vec>(1)));
          mu_k = (mu_k-0.5)*0.9999999+0.5;
          
          // create prob for current item, update prob matrix for respective knots
          help_pow = mu_k % yi_k - (yi_k-1);
          
          prods_i(j,jj) = prods_i(j,jj)*prod(help_pow);
        }
        
      }
    }
    
    f(i)  = -log(accu(prods_i%GHweights));
    
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  double fx = sum(f) + P2;

  return (fx);
}

// [[Rcpp::export]]
double loglikMO_cumul(arma::vec alpha,
                arma::vec Y,
                arma::mat X,
                int Q,
                int q,
                int n,
                int I,
                int pall,
                int pX,
                int pXRS,
                int pthresh,
                int pshift,
                int prnd,
                arma::mat GHweights,
                arma::vec GHnodes,
                int scaled,
                int cores,
                double lambda) { 
  
  // initialize loglikelihood       
  
  vec f = zeros(n);

  
  // current value of L2 penalty
  double P2 = accu(alpha%alpha)*lambda;
  
  // create sigma matrix from current parameters
  double co_var = alpha(pall-2)*sqrt(alpha(pall-1))*sqrt(alpha(pall-3));
  mat sigma = zeros(2,2);        
  sigma(0,0) = alpha(pall-3);
  sigma(1,0) = co_var;
  sigma(0,1) = co_var;
  sigma(1,1) = alpha(pall-1);
  
  
  mat sigma12 = chol(sigma);
  
  vec etaTheta = zeros(n);
  vec etaGamma = zeros(n);
  if(pX>0){
    etaTheta = X * alpha(span(pthresh+pshift,pthresh+pshift+pX-1));
  }
  if(pXRS>0){
    etaGamma = X * alpha(span(pthresh+pshift+pX,pthresh+pshift+pX+pXRS-1));
  }
  
  vec intercepts = alpha(span(0,I-1));
  
  // initialize design for all thresholds
  mat threshs = zeros(q,I);
  mat expvals = trans(reshape(exp(alpha(span(I,pthresh-1))),I,q-1));
  int m;
  
  // falls Rest Null ist k ungerade
  if(q%2==0){
    m = (q+2)/2;

    for(int i=0; i < q; i++){

      if(i==(m-2)){
        threshs(i,span::all) = -expvals(m-2,span::all)/2;
      }
      if(i==(m-1)){
        threshs(i,span::all) = expvals(m-2,span::all)/2;
      }
      if(i<(m-2)){
        threshs(i,span::all) =  -expvals(m-2,span::all)/2 - sum(expvals(span(i, m-3),span::all), 0);
      }
      
      if(i>(m-1)){
        threshs(i,span::all) = expvals(m-2,span::all)/2 + sum(expvals(span(m-1,i-1),span::all), 0);
      }
    }
  }else{
    m = (q+1)/2;
    for(int i=0; i < q; i++){
      // threshs(i,span::all) = intercepts(span::all);
      if(i<(m-1)){
        threshs(i,span::all) = -sum(expvals(span(i, m-2),span::all), 0);
      }
      if(i>(m-1)){
        threshs(i,span::all) =  sum(expvals(span(m-1,i-1),span::all), 0);
      }
      
    }
  }

  // initialize number of different threads
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif
  
  // loop through all persons
// #pragma omp parallel for reduction(-:f)

int i; int pos; vec yi; mat prods_i; int j; int jj; vec rnd_act; int k; vec yi_k; 
vec eta_k; vec mu_k; vec help_pow;


#ifdef _OPENMP
#pragma omp parallel for private(i, pos, yi, prods_i, j, jj, rnd_act, k,  yi_k, eta_k, mu_k, help_pow) shared(f)
#endif
  for(i=0; i<n ;i++){
    pos = i*q*I;
    // get response of person i
    yi = Y(span(pos,pos+q*I-1));
    
    prods_i = ones(Q,Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 
      for(jj=0; jj<Q ;jj++){ 
        // initialize current random effects or knots
        rnd_act = zeros(2);
        rnd_act(0) = GHnodes(j);
        rnd_act(1) = GHnodes(jj);
        
        rnd_act = sigma12*rnd_act;
        
        rnd_act(1) = exp(rnd_act(1) + etaGamma(i)); 
        rnd_act(0) = rnd_act(0) + etaTheta(i); 
        
        // loop through all items  
        for(k=0; k<I ;k++){  
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          
          yi_k = diff(join_cols(join_cols(zeros<vec>(1),yi_k),ones<vec>(1)));
          
          // get eta and mu of person i for current item k

          eta_k = intercepts(k) + threshs(span::all,k)*rnd_act(1) + rnd_act(0);
          

          mu_k = respFun_cumul(eta_k);
          
          mu_k = diff(join_cols(join_cols(zeros<vec>(1),mu_k),ones<vec>(1)));
          mu_k = (mu_k-0.5)*0.9999999+0.5;

          // create prob for current item, update prob matrix for respective knots
          help_pow = mu_k % yi_k - (yi_k-1);
          
          prods_i(j,jj) = prods_i(j,jj)*prod(help_pow);
        }
        
      }
    }
    
    f(i)  = -log(accu(prods_i%GHweights));
    
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  double fx = sum(f) + P2;

  if(std::isnan(fx)){
    fx = loglikMO_cumul2(alpha,
                         Y,
                         X,
                         Q,
                         q,
                         n,
                         I,
                         pall,
                         pX,
                         pXRS,
                         pthresh,
                         pshift,
                         prnd,
                         GHweights,
                         GHnodes,
                         scaled,
                         cores,
                         lambda);
  }
      
  return (fx);
}




// [[Rcpp::export]]
double loglikMO_cumul_noRS(arma::vec alpha,
                      arma::vec Y,
                      arma::mat X,
                      int Q,
                      int q,
                      int n,
                      int I,
                      int pall,
                      int pX,
                      int pthresh,
                      int pshift,
                      int prnd,
                      arma::vec GHweights,
                      arma::vec GHnodes,
                      int scaled,
                      int cores,
                      double lambda) { 
  
  // initialize loglikelihood       
  
  vec f = zeros(n);
  
  
  // current value of L2 penalty
  double P2 = accu(alpha%alpha)*lambda;
  
  vec etaTheta = zeros(n);
  if(pX>0){
    etaTheta = X * alpha(span(pthresh+pshift,pthresh+pshift+pX-1));
  }
  
  vec intercepts = alpha(span(0,I-1));
  
  // initialize design for all thresholds
  mat threshs = zeros(q,I);
  mat expvals = trans(reshape(exp(alpha(span(I,pthresh-1))),I,q-1));
  int m;
  
  // falls Rest Null ist k ungerade
  if(q%2==0){
    m = (q+2)/2;
    
    for(int i=0; i < q; i++){
      
      if(i==(m-2)){
        threshs(i,span::all) = -expvals(m-2,span::all)/2;
      }
      if(i==(m-1)){
        threshs(i,span::all) = expvals(m-2,span::all)/2;
      }
      if(i<(m-2)){
        threshs(i,span::all) =  -expvals(m-2,span::all)/2 - sum(expvals(span(i, m-3),span::all), 0);
      }
      
      if(i>(m-1)){
        threshs(i,span::all) = expvals(m-2,span::all)/2 + sum(expvals(span(m-1,i-1),span::all), 0);
      }
    }
  }else{
    m = (q+1)/2;
    for(int i=0; i < q; i++){
      // threshs(i,span::all) = intercepts(span::all);
      if(i<(m-1)){
        threshs(i,span::all) = -sum(expvals(span(i, m-2),span::all), 0);
      }
      if(i>(m-1)){
        threshs(i,span::all) =  sum(expvals(span(m-1,i-1),span::all), 0);
      }
      
    }
  }

  // initialize number of different threads
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif
  
  // loop through all persons
  // #pragma omp parallel for reduction(-:f)
  
  int i; int pos; vec yi; mat prods_i; int j; int k; vec yi_k; 
  vec eta_k; vec mu_k; vec help_pow;
  
  
#ifdef _OPENMP
#pragma omp parallel for private(i, pos, yi, prods_i, j, k,  yi_k, eta_k, mu_k, help_pow) shared(f)
#endif
  for(i=0; i<n ;i++){
    pos = i*q*I;
    // get response of person i
    yi = Y(span(pos,pos+q*I-1));
    
    prods_i = ones(Q);
    // loop through all knots, for both random effects
    for(j=0; j<Q ;j++){ 

        
        // loop through all items  
        for(k=0; k<I ;k++){  
          // response of person i for current item k
          yi_k = yi(span(k*q,q+k*q-1));
          
          yi_k = diff(join_cols(join_cols(zeros<vec>(1),yi_k),ones<vec>(1)));
          
          // get eta and mu of person i for current item k
          
          eta_k = intercepts(k) + threshs(span::all,k) + etaTheta(i) + GHnodes(j)*sqrt(alpha(pall-1));
          
          mu_k = respFun_cumul(eta_k);
          mu_k = diff(join_cols(join_cols(zeros<vec>(1),mu_k),ones<vec>(1)));
          // create prob for current item, update prob matrix for respective knots
          help_pow = mu_k % yi_k - (yi_k-1);
          
          prods_i(j) = prods_i(j)*prod(help_pow);
      }
    }
    f(i)  = -log(accu(prods_i%GHweights));
    
    // accumulate all likelihood contributions, weights by respective weights and probs of knots
  }
  double fx = sum(f) + P2;
  
  return (fx);
}
