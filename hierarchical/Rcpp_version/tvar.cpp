// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <RcppDist.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
Rcpp::List ffbs_tvar(arma::vec yt, arma::mat FF,
                double n0, double S0,
                double delta){
  int n_t = FF.n_rows;
  int P = FF.n_cols;
  arma::vec nt(n_t, arma::fill::zeros);
  arma::vec dt(n_t, arma::fill::zeros);
  arma::vec St(n_t, arma::fill::zeros);

  // TVAR coefficients filtering part
  arma::mat at(P, n_t, arma::fill::zeros);
  arma::mat mt(P, n_t, arma::fill::zeros);
  arma::vec ft(n_t, arma::fill::zeros);
  arma::vec et(n_t, arma::fill::zeros);
  arma::cube Rt(P, P, n_t, arma::fill::zeros);
  arma::cube Ct(P, P, n_t, arma::fill::zeros);
  arma::vec Qt(n_t, arma::fill::zeros);
  arma::mat At;
  double ll = 0.0;

  // TVAR coefficients smoothing part
  arma::mat mnt(P, n_t, arma::fill::zeros);
  arma::cube Cnt(P, P, n_t, arma::fill::zeros);
  arma::mat Bt(P, P, arma::fill::zeros);
  // initial states
  arma::vec m0(P, arma::fill::zeros);
  arma::mat C0(P, P, arma::fill::eye);
  double d0 = n0 * S0;

  for(int i = 0; i < n_t; i++){

    if(i == 0){
      at.col(i) = m0;
      Rt.slice(i) = S0 * C0 / delta;
      Rt.slice(i) = 0.5*Rt.slice(i) + 0.5*arma::trans(Rt.slice(i));
      Qt(i) = arma::as_scalar(S0 + FF.row(i) * Rt.slice(i) * arma::trans(FF.row(i)));

    }else{
      at.col(i) = mt.col(i-1);
      Rt.slice(i) = Ct.slice(i-1) / delta;
      Rt.slice(i) = 0.5*Rt.slice(i) + 0.5*arma::trans(Rt.slice(i));
      Qt(i) = arma::as_scalar(St(i-1) + FF.row(i) * Rt.slice(i) * arma::trans(FF.row(i)));
    }


    ft(i) = arma::as_scalar(FF.row(i) * at.col(i));
    et(i) = yt(i) - ft(i);
    ll += arma::log_normpdf(yt(i), ft(i), std::sqrt(Qt(i)));
    if(i == 0){
      nt(i) = n0 + 1.0;
      dt(i) = d0 + S0 * std::pow(et(i), 2) / Qt(i);
    }else{
      nt(i) = nt(i-1) + 1.0;
      dt(i) = dt(i-1) + St(i-1) * std::pow(et(i), 2) / Qt(i);
    }
    St(i) = dt(i) / nt(i);
 
    At = Rt.slice(i) * arma::trans(FF.row(i)) / Qt(i);

    mt.col(i) = at.col(i) + At * et(i);
    if(i == 0){
      Ct.slice(i) = (St(i)/S0)*(Rt.slice(i) - At * arma::trans(At) * Qt(i));
    }else{
      Ct.slice(i) = (St(i)/St(i-1))*(Rt.slice(i) - At * arma::trans(At) * Qt(i));
    }
    Ct.slice(i) = 0.5*Ct.slice(i) + 0.5*arma::trans(Ct.slice(i));
    //Rcout << "Ct is" << std::endl << Ct.slice(i) << std::endl;
  }

  // smooth part
  mnt.col(n_t-1) = mt.col(n_t-1);
  Cnt.slice(n_t-1) = Ct.slice(n_t-1);
  for(int i = n_t-2; i > 0; i--){
    //Rcout << "Rt is" << std::endl << Rt.slice(i+1) << std::endl;

    arma::mat inv_Rtp1 = arma::inv_sympd(Rt.slice(i+1));
    Bt = Ct(i) * inv_Rtp1;
    mnt.col(i) = mt.col(i) + Bt*(mnt.col(i+1) - at.col(i+1));
    Cnt.slice(i) = (Ct.slice(i) + Bt*(Cnt.slice(i+1) - Rt.slice(i+1))*arma::trans(Bt));
    Cnt.slice(i) = 0.5*Cnt.slice(i) + 0.5*arma::trans(Cnt.slice(i));
  }
  for(int i = 0; i < n_t; i++){
    Cnt.slice(i) = St(n_t-1)/St(i) * Cnt.slice(i);
  }
  return Rcpp::List::create(Rcpp::Named("mt") = mt,
                            Rcpp::Named("Ct") = Ct,
                            Rcpp::Named("mnt") = mnt,
                            Rcpp::Named("Cnt") = Cnt,
                            Rcpp::Named("et") = et,
                            Rcpp::Named("St") = St,
                            Rcpp::Named("ll") = ll);

}


// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::export]]

Rcpp::List tvar(arma::vec yt, arma::mat FF,
                double n0, double S0,
                arma::vec delta){
  
  int delta_n = delta.n_elem;
  double delta_min;
  Rcpp::List result_opt = ffbs_tvar(yt, FF, n0, S0, delta(0));
  double ll_max = result_opt["ll"];
  for(int i = 1; i < delta_n; i++){
    Rcpp::List result_new = ffbs_tvar(yt, FF, n0, S0, delta(i));
    double ll_new = result_new["ll"];
    if(ll_max < ll_new){
      result_opt = result_new;
      ll_max = ll_new;
      delta_min = delta(i);
    }
  }
  return Rcpp::List::create(Rcpp::Named("mt") = result_opt["mt"],
                            Rcpp::Named("Ct") = result_opt["Ct"],
                            Rcpp::Named("mnt") = result_opt["mnt"], 
                            Rcpp::Named("Cnt") = result_opt["Cnt"],
                            Rcpp::Named("St") = result_opt["St"],
                            Rcpp::Named("ll_max") = ll_max,
                            Rcpp::Named("delta_min") = delta_min);
}

