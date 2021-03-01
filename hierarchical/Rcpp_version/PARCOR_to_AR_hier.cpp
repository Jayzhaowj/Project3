

#include <RcppArmadillo.h>
#include <Rcpp.h>


// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List PARCOR_to_AR(arma::mat phi_forward, arma::mat phi_backward,
                        arma::cube akm_prev, arma::cube dkm_prev, 
                        int cur_level){
  int n_t = phi_forward.n_cols;
  int n_I = phi_forward.n_rows;
  arma::cube akm_cur(n_I, n_t, cur_level + 1);
  arma::cube dkm_cur(n_I, n_t, cur_level + 1);
  if (cur_level == 0){
    akm_cur.slice(0) = phi_forward;
    dkm_cur.slice(0) = phi_backward;
  }
  else
  {
    akm_cur.slice(cur_level) = phi_forward;
    dkm_cur.slice(cur_level) = phi_backward;
    for(int i = 0; i < cur_level; i++){
      arma::mat akm_temp(n_I, n_t);
      arma::mat dkm_temp(n_I, n_t);
      arma::mat akm_i_prev = akm_prev.slice(i);
      arma::mat dkm_i_prev = dkm_prev.slice(i);
      arma::mat akmm_i_prev = akm_prev.slice(cur_level - i - 1);
      arma::mat dkmm_i_prev = dkm_prev.slice(cur_level - i - 1);
      //for(int j = PP; j < (n_t - PP); j++){
      for(int j = 0; j < n_t; j++){
        akm_temp.col(j) = akm_i_prev.col(j) - phi_forward.col(j) % dkmm_i_prev.col(j);
        dkm_temp.col(j) = dkm_i_prev.col(j) - phi_backward.col(j) % akmm_i_prev.col(j);
      }
      akm_cur.slice(i) = akm_temp;
      dkm_cur.slice(i) = dkm_temp;
    }
  }
  return Rcpp::List::create(Rcpp::Named("forward") = akm_cur, Rcpp::Named("backward") = dkm_cur);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List PAR_to_AR_fun(arma::cube phi_fwd, arma::cube phi_bwd){
  int n_t = phi_fwd.n_cols;
  int P = phi_fwd.n_slices;
  int n_I = phi_fwd.n_rows;
  arma::cube akm_prev(n_I, n_t, P);
  arma::cube dkm_prev(n_I, n_t, P);
  Rcpp::List ar_coef_prev;
  Rcpp::List ar_coef(P);
  for(int i = 0; i < P; i++){
    if(i == 0){
      ar_coef(i) = PARCOR_to_AR(phi_fwd.slice(i), phi_bwd.slice(i), akm_prev, dkm_prev, i);
    }else{
      ar_coef_prev = ar_coef(i - 1);
      akm_prev = Rcpp::as<arma::cube>(ar_coef_prev["forward"]);
      dkm_prev = Rcpp::as<arma::cube>(ar_coef_prev["backward"]);
      ar_coef(i) = PARCOR_to_AR(phi_fwd.slice(i), phi_bwd.slice(i), akm_prev, dkm_prev, i);
    }
  }
  return ar_coef;
}


