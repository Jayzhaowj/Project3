// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
Rcpp::List forward_filter_backward_smooth(arma::mat yt, arma::mat F1, arma::mat F2, 
                                          int n_t, int I, int m,  int type, 
                                          double delta1, double delta2){
  // some constants
  int sign = 1;
  arma::mat I_n(I, I, arma::fill::eye);
  // initial states
  arma::colvec mk_0(I, arma::fill::zeros);
  arma::mat Ck_0(I, I, arma::fill::eye);
  arma::mat Ck_s_0(I, I, arma::fill::eye);
  double n_0 = 1;
  double d_0 = 1;
  arma::vec nt(n_t, arma::fill::zeros);
  arma::vec dt(n_t, arma::fill::zeros);
  arma::vec St(n_t, arma::fill::zeros);
  
  // parcor states
  arma::mat at(I, n_t, arma::fill::zeros);
  arma::mat mt(I, n_t, arma::fill::zeros);
  arma::mat ft(I, n_t, arma::fill::zeros);
  arma::mat et(I, n_t, arma::fill::zeros);
  
  arma::cube Rt(I, I, n_t, arma::fill::zeros);
  arma::cube Rt_s(I, I, n_t, arma::fill::zeros);
  arma::cube Ct(I, I, n_t, arma::fill::zeros);
  arma::cube Ct_s(I, I, n_t, arma::fill::zeros);
  arma::cube Ut_s(I, I, n_t, arma::fill::zeros);
  arma::cube Qt(I, I, n_t, arma::fill::zeros);
  arma::cube Qt_s(I, I, n_t, arma::fill::zeros);
  arma::cube inv_Qt_s(I, I, n_t, arma::fill::zeros);
  
  //structure level
  arma::mat akt(I, n_t, arma::fill::zeros);
  arma::mat mkt(I, n_t, arma::fill::zeros);
  arma::cube Rkt(I, I, n_t, arma::fill::zeros);
  arma::cube Rkt_s(I, I, n_t, arma::fill::zeros);
  arma::cube Ckt(I, I, n_t, arma::fill::zeros);
  arma::cube Ckt_s(I, I, n_t, arma::fill::zeros);
  arma::cube V2t(I, I, n_t, arma::fill::zeros);
  arma::cube V2t_s(I, I, n_t, arma::fill::zeros);
  arma::cube Ukt_s(I, I, n_t, arma::fill::zeros);
  arma::cube F1t(I, I, n_t, arma::fill::zeros);
  
  // smooth part
  arma::mat mnt(I, n_t, arma::fill::zeros);
  arma::cube Cnt(I, I, n_t, arma::fill::zeros);
  arma::mat mnkt(I, n_t, arma::fill::zeros);
  arma::cube Cnkt(I, I, n_t, arma::fill::zeros);
  arma::cube Ant(I, I, n_t, arma::fill::zeros);
  arma::cube Ankt(I, I, n_t, arma::fill::zeros);
  arma::mat resid(I, n_t, arma::fill::zeros);
  
  // asign the lower bound and upper bound
  int ubound = 0;
  int lbound = 0;
  if(type == 1){
    ubound = n_t;
    lbound = m;
    sign = 1;
  }else{
    ubound = n_t - m;
    lbound = 0;
    sign = -1;
  }
  
  
  // filtering conditional on V:
  for(int i = lbound; i < ubound; i++){
    
    // prior distribution update
    
    F1t.slice(i) = arma::diagmat(F1.col(i - sign*m)); 
    if(i == lbound){
      akt.col(i) = mk_0;
      Rkt_s.slice(i) = Ck_s_0/delta1;
    }else{
      akt.col(i) = mkt.col(i-1);
      Rkt_s.slice(i) = Ckt_s.slice(i-1)/delta1;
    }
    
    at.col(i) = F2*akt.col(i);
    Rt_s.slice(i) = F2*Rkt_s.slice(i)*arma::trans(F2)/delta2;
    V2t_s.slice(i) = (1 - delta2)/delta2 * F2 * Rkt_s.slice(i) * arma::trans(F2);
    
    // predictive distribution update
    ft.col(i) = F1t.slice(i) * at.col(i);
    et.col(i) = yt.col(i) - ft.col(i);
    if(i == lbound){
      Qt_s.slice(i) = F1t.slice(i) * Rt_s.slice(i) * arma::trans(F1t.slice(i)) + I_n;
      Qt_s.slice(i) = 0.5 * Qt_s.slice(i) + 0.5 * arma::trans(Qt_s.slice(i));
    }else{
      Qt_s.slice(i) = F1t.slice(i) * Rt_s.slice(i) * arma::trans(F1t.slice(i)) + I_n;
      Qt_s.slice(i) = 0.5 * Qt_s.slice(i) + 0.5 * arma::trans(Qt_s.slice(i));
    }
    
    // posterior distribution
    // evolution equation
    Ukt_s.slice(i) = Rkt_s.slice(i) * arma::trans(F1t.slice(i)*F2);
    inv_Qt_s.slice(i) = arma::inv_sympd(Qt_s.slice(i));
    //arma::mat inv_Qt = arma::inv(Qt.slice(i));
    //Rprintf("computation is completed! \n");

    
    mkt.col(i) = akt.col(i) + Ukt_s.slice(i) * inv_Qt_s.slice(i)*et.col(i);
    Ckt_s.slice(i) = Rkt_s.slice(i) - Ukt_s.slice(i)*inv_Qt_s.slice(i)*arma::trans(Ukt_s.slice(i));
    
    // Structural equation
    Ut_s.slice(i) = Rt_s.slice(i) * arma::trans(F1t.slice(i));
    mt.col(i) = at.col(i) + Ut_s.slice(i) * inv_Qt_s.slice(i) * et.col(i);
    Ct_s.slice(i) = Rt_s.slice(i) - Ut_s.slice(i) * inv_Qt_s.slice(i) * arma::trans(Ut_s.slice(i));
  }
  // for the precision V^-1
  for(int i = lbound; i < ubound; i++){
    if(i == lbound){
      nt(i) = n_0;
      dt(i) = d_0;
      
    }else{
      nt(i) = nt(i-1) + I;
      dt(i) = dt(i-1) + arma::as_scalar(arma::trans(et.col(i)) * inv_Qt_s.slice(i) * et.col(i));
    }
    St(i) = dt(i)/nt(i);
  }
  
  // filtering unconditional on V
  for(int i = lbound; i < ubound; i++){
    if(i == lbound){
      Qt.slice(i) = d_0/n_0*Qt_s.slice(i);
      
      Rt.slice(i) = d_0/n_0*Rt_s.slice(i);
      Rkt.slice(i) = d_0/n_0*Rkt_s.slice(i);
      
      // V2t.slice(i) = d_0/n_0*V2t_s.slice(i);
      // 
      // Ct.slice(i) = d_0/n_0*Ct_s.slice(i);
      // Ckt.slice(i) = d_0/n_0*Ckt_s.slice(i);
    }else{
      Qt.slice(i) = dt(i-1)/nt(i-1)*Qt_s.slice(i);
      
      Rt.slice(i) = dt(i-1)/nt(i-1)*Rt_s.slice(i);
      Rkt.slice(i) = dt(i-1)/nt(i-1)*Rkt_s.slice(i);
      
      // V2t.slice(i) = dt(i)/nt(i)*V2t_s.slice(i);
      // 
      // Ct.slice(i) = dt(i)/nt(i)*Ct_s.slice(i);
      // Ckt.slice(i) = dt(i)/nt(i)*Ckt_s.slice(i);
    }
    V2t.slice(i) = dt(i)/nt(i)*V2t_s.slice(i);
    
    Ct.slice(i) = dt(i)/nt(i)*Ct_s.slice(i);
    Ckt.slice(i) = dt(i)/nt(i)*Ckt_s.slice(i);
  }
  
  // smooth part
  mnt.col(ubound-1) = mt.col(ubound-1);
  Cnt.slice(ubound-1) = Ct.slice(ubound-1);
  mnkt.col(ubound-1) = mkt.col(ubound-1);
  Cnkt.slice(ubound-1) = Ckt.slice(ubound-1);
  resid.col(ubound-1) = yt.col(ubound-1) - F1t.slice(ubound-1) * mnt.col(ubound-1);

  for(int i = (ubound - 2); i > lbound - 1; i--){
    //F1t.slice(i) = arma::diagmat(F1.col(i - sign*m)); 
    arma::mat V02t_s = dt(ubound-1)/nt(ubound-1) * I_n + V2t.slice(i);
    arma::mat inv_V02t_s = arma::inv_sympd(0.5*V02t_s + 0.5*arma::trans(V02t_s));
 
    Ant.slice(i) = F2*Ckt.slice(i)*arma::trans((I_n - V2t.slice(i)*arma::trans(F1t.slice(i))*inv_V02t_s*F1t.slice(i))*F2);
    Ankt.slice(i) = Ckt.slice(i);

    arma::mat inv_Rtp1 = arma::inv_sympd(0.5*Rt.slice(i+1) + 0.5*arma::trans(Rt.slice(i+1)));
    arma::mat inv_Rktp1 = arma::inv_sympd(0.5*Rkt.slice(i+1) + 0.5*arma::trans(Rkt.slice(i+1)));

    mnt.col(i) = mt.col(i) + arma::trans(Ant.slice(i))*inv_Rtp1*(mnt.col(i+1) - at.col(i+1));
    Cnt.slice(i) = (Ct.slice(i) - arma::trans(Ant.slice(i))*inv_Rtp1*(Rt.slice(i+1) - Cnt.slice(i+1))*inv_Rtp1*Ant.slice(i));
    
    mnkt.col(i) = mkt.col(i) + arma::trans(Ankt.slice(i))*inv_Rktp1*(mnkt.col(i+1) - akt.col(i+1));
    Cnkt.slice(i) = (Ckt.slice(i) - arma::trans(Ankt.slice(i))*inv_Rktp1*(Rkt.slice(i+1) - Cnkt.slice(i+1))*inv_Rktp1*Ankt.slice(i+1));
    
    resid.col(i) = yt.col(i) - F1t.slice(i) * mnt.col(i);

  }
  return Rcpp::List::create(Rcpp::Named("mnt") = mnt, 
                            Rcpp::Named("mnkt") = mnkt,
                            Rcpp::Named("Cnt") = Cnt,
                            Rcpp::Named("Cnkt") = Cnkt,
                            Rcpp::Named("mt") = mt,
                            Rcpp::Named("Ct") = Ct,
                            //Rcpp::Named("Ant") = Ant,
                            //Rcpp::Named("Ankt") = Ankt,
                            Rcpp::Named("residuals") = resid,
                            Rcpp::Named("sigma2t") = dt/nt, 
                            Rcpp::Named("ft") = ft,
                            Rcpp::Named("Qt") = Qt,
                            Rcpp::Named("nt") = nt);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
