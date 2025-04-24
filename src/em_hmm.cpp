#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// double non0_min(vec &p);
vec my_pava(vec &values, vec &weight, bool decreasing);
arma::vec replis(arma::vec pi, arma::mat A, arma::vec f1, arma::vec f2);
arma::vec fdr(arma::vec repLIS);

// [[Rcpp::export]]
SEXP em_hmm(SEXP pa_in, SEXP pb_in, SEXP pi0a_in, SEXP pi0b_in) {
  try{
    const int maxIter = 200;
    const double tol = 1e-3;
    const double pvalue_cutoff = 1e-15;
    const double f_cutoff = 1e-15;
    
    vec pa = as<arma::vec>(pa_in), pb = as<arma::vec>(pb_in);
    const double pi0_pa = Rcpp::as<double>(pi0a_in), pi0_pb = Rcpp::as<double>(pi0b_in);
    
    int J = pa.size();
    double min_a = pa.elem(find(pa>0)).min(), min_b = pb.elem(find(pb>0)).min();
    pa.elem(find(pa<=0)).fill(pvalue_cutoff < min_a ? pvalue_cutoff : min_a);
    pb.elem(find(pb<=0)).fill(pvalue_cutoff < min_b ? pvalue_cutoff : min_b);
    
    vec p1 = sort(pa), p2 = sort(pb);
    uvec ix1 = sort_index(pa), ix2 = sort_index(pb);
    
    vec p1_diff(J), p2_diff(J);
    p1_diff(0) = p1(0);
    p2_diff(0) = p2(0);
    for(int i = 1; i<J; ++i){
      p1_diff(i) = p1(i) - p1(i-1);
      p2_diff(i) = p2(i) - p2(i-1);
    }
    
    // Initialization
    vec pi(4);
    pi(0) = pi0_pa * pi0_pb; 
    pi(1) = pi0_pa * (1 - pi0_pb);
    pi(2) = (1 - pi0_pa) * pi0_pb; 
    pi(3) = (1 - pi0_pa) * (1 - pi0_pb);
    
    mat A = {{0.9, 0.04, 0.04, 0.02},
    {0.28, 0.3, 0.14, 0.28},
    {0.28, 0.14, 0.3, 0.28},
    {0.14, 0.28, 0.28, 0.3}};
    
    vec f0 = ones(J,1), f1(J), f2(J);
    f1 = 1 - pa;
    f2 = 1 - pb;
    
    mat f = join_rows(f0, f2, f1, f1%f2);
    
    vec loglik(maxIter);
    loglik(0) = -datum::inf;
    
    mat alpha = zeros(J, 4), beta = zeros(J, 4), gamma = zeros(J, 4);
    cube xi = zeros(4, 4, J-1);
    
    vec xi_kl(J-1), Q1(J), Q2(J), y1(J), y2(J), res1(J), res2(J);
    mat xi_k(4,J-1);
    mat f_jp1(4,4);
    double sub_loglik, loglik_delta;
    
    std::cout << "EM begins:" << std::endl;
    
    for (int i = 1; i < maxIter; i++){
      // E-step
      // calculate the forward, backward and posterior probabilities based on current HMM
      alpha.row(0) = pi.t() % f.row(0);
      alpha.row(0) = alpha.row(0)/sum(alpha.row(0));
      
      for (int j = 1; j < J; j++){
        alpha(j,0) = sum(alpha.row(j-1).t() % A.col(0)) * f(j,0);
        alpha(j,1) = sum(alpha.row(j-1).t() % A.col(1)) * f(j,1);
        alpha(j,2) = sum(alpha.row(j-1).t() % A.col(2)) * f(j,2);
        alpha(j,3) = sum(alpha.row(j-1).t() % A.col(3)) * f(j,3);
        
        alpha.row(j) = alpha.row(j)/sum(alpha.row(j));
      }
      
      beta.row(J-1).fill(0.25);
      for(int j = J-2; j >=0; j--){
        beta(j,0) = sum(beta.row(j+1) % A.row(0) % f.row(j+1));
        beta(j,1) = sum(beta.row(j+1) % A.row(1) % f.row(j+1));
        beta(j,2) = sum(beta.row(j+1) % A.row(2) % f.row(j+1));
        beta(j,3) = sum(beta.row(j+1) % A.row(3) % f.row(j+1));
        
        beta.row(j) = beta.row(j)/sum(beta.row(j));
      }
      
      for(int j = 0; j < J; j++){
        gamma.row(j) = alpha.row(j) % beta.row(j) / sum(alpha.row(j) % beta.row(j));
      }
      
      for(int j = 0; j < J-1; j++){
        f_jp1 = repmat(f.row(j+1), 4, 1);
        xi.slice(j) = alpha.row(j).t() * beta.row(j+1) % A % f_jp1/
          accu(alpha.row(j).t() * beta.row(j+1) % A % f_jp1);
      }
      
      // M-step
      // update the parameters pi and A
      pi = gamma.row(0).t();
      for(int k = 0; k < 4; k++){
        for(int l = 0; l < 4; l++){
          xi_kl = xi.tube(k,l);
          xi_k = xi.row(k);
          A(k,l) = sum(xi_kl)/accu(xi_k);
        }
      }
      
      // update f1 and f2
      Q1 = gamma.col(2) + gamma.col(3);
      Q2 = gamma.col(1) + gamma.col(3);
      Q1 = Q1(ix1);
      Q2 = Q2(ix2);
      
      y1 = - p1_diff * sum(Q1) / Q1;
      y2 = - p2_diff * sum(Q2) / Q2;
      
      y1.elem(find_nonfinite(y1)).fill(y1.elem(find_finite(y1)).min());
      y2.elem(find_nonfinite(y2)).fill(y2.elem(find_finite(y2)).min());
      
      res1 = my_pava(y1, Q1, true);
      res2 = my_pava(y2, Q2, true);
      
      f1 = -1 / res1;
      f1 = f1 / sum(f1 % p1_diff);
      f1(ix1) = f1;
      f1.elem(find_nan(f1)).fill(f1.min());
      
      f2 = -1 / res2;
      f2 = f2 / sum(f2 % p2_diff);
      f2(ix2) = f2;
      f2.elem(find_nan(f2)).fill(f2.min());
      
      double min_f1 = f1.elem(find(f1>0)).min(), min_f2 = f2.elem(find(f2>0)).min();
      f1.elem(find(f1<=0)).fill(f_cutoff < min_f1 ? f_cutoff : min_f1);
      f2.elem(find(f2<=0)).fill(f_cutoff < min_f2 ? f_cutoff : min_f2);
      
      f = join_rows(f0, f2, f1, f1%f2);
      
      
      // calculate the updated log-likelihood
      sub_loglik = 0;
      for(int j = 0; j < J - 1; j++){
        sub_loglik = sub_loglik + accu(log(A) % xi.slice(j));
      }
      
      loglik(i) = sum(log(pi).t() % gamma.row(1)) + sub_loglik + accu(gamma % log(f));
      loglik_delta = abs((loglik(i) - loglik(i-1))/loglik(i-1));
      
      std::cout<<i<<". "<< loglik(i) << ", delta = " << loglik_delta << std::endl;
      
      if(loglik_delta < tol){
        break;
      }
    }
    
    // === Update pi using dominant eigenvector of A^T ===
    arma::mat At = A.t();
    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_gen(eigval, eigvec, At);
    arma::uword max_index = index_max(abs(eigval));
    pi = real(eigvec.col(max_index));
    pi = pi / sum(pi);  // Normalize
    
    vec repLIS = replis(pi, A, f1, f2);
    vec radj = fdr(repLIS);
    
    return Rcpp::List::create(Rcpp::Named("repLIS") = repLIS,
                              Rcpp::Named("fdr") = radj,
                              Rcpp::Named("loglik") = loglik,
                              Rcpp::Named("pi") = pi.t(),
                              Rcpp::Named("A") = A,
                              Rcpp::Named("f1") = f1,
                              Rcpp::Named("f2") = f2);
  } catch( std::exception &ex ) {
    forward_exception_to_r(ex);
    return Rcpp::List::create();
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)..." );
    return Rcpp::List::create();
  }
}

// double non0_min(vec &p){
//   double _min = std::numeric_limits<double>::max();
//   
//   p.for_each([&_min](double &val) { if(val > 0 && val < _min) _min = val; });
//   
//   return _min;
// }

double na_rm(vec &p){
  double _min = std::numeric_limits<double>::max();
  
  p.for_each([&_min](double &val) { if(val > 0 && val < _min) _min = val; });
  
  return _min;
}

// [[Rcpp::export]]
arma::vec replis(arma::vec pi, arma::mat A, arma::vec f1, arma::vec f2){
  int J = f1.size();
  vec f0 = ones(J,1);
  mat f = join_rows(f0, f2, f1, f1%f2);
  
  mat alpha = zeros(J, 4), beta = zeros(J, 4), f_jp1(4,4);
  
  alpha.row(0) = pi.t() % f.row(0);
  alpha.row(0) = alpha.row(0)/sum(alpha.row(0));
  
  for (int j = 1; j < J; j++){
    alpha(j,0) = sum(alpha.row(j-1).t() % A.col(0)) * f(j,0);
    alpha(j,1) = sum(alpha.row(j-1).t() % A.col(1)) * f(j,1);
    alpha(j,2) = sum(alpha.row(j-1).t() % A.col(2)) * f(j,2);
    alpha(j,3) = sum(alpha.row(j-1).t() % A.col(3)) * f(j,3);
    
    alpha.row(j) = alpha.row(j)/sum(alpha.row(j));
  }
  
  beta.row(J-1).fill(0.25);
  for(int j = J-2; j >=0; j--){
    beta(j,0) = sum(beta.row(j+1) % A.row(0) % f.row(j+1));
    beta(j,1) = sum(beta.row(j+1) % A.row(1) % f.row(j+1));
    beta(j,2) = sum(beta.row(j+1) % A.row(2) % f.row(j+1));
    beta(j,3) = sum(beta.row(j+1) % A.row(3) % f.row(j+1));
    
    beta.row(j) = beta.row(j)/sum(beta.row(j));
  }
  
  vec repLIS(J);
  for(int j = 0; j < J; j++){
    repLIS(j) = sum(alpha.row(j).head(3) % beta.row(j).head(3))/
      sum(alpha.row(j) % beta.row(j));
  }
  
  return repLIS;
}

// [[Rcpp::export]]
arma::vec fdr(arma::vec repLIS){
  int J = repLIS.size();
  
  vec ordered_lis = sort(repLIS), s = linspace(1,J,J);
  uvec ix_lis = sort_index(repLIS);
  
  vec radj = cumsum(ordered_lis)/s;
  radj(ix_lis) = radj;
  
  return radj;
}

// inline bool compare(double x, double y, Ordering ordering)
// {
//   return ordering == Ordering::Increasing ? x > y : x < y;
// }

vec my_pava(vec &values, vec &weight, bool decreasing)
{
  if(decreasing){
    values = reverse(values);
    weight = reverse(weight);
  }
  vec w(values.size(), fill::zeros);
  vec x(values.size(), fill::zeros);
  x[0] = values[0];
  w[0] = weight[0];
  unsigned j = 0;
  vec s(values.size(), fill::zeros);
  
  for (unsigned i = 1; i < values.size(); i++) {
    j += 1;
    x[j] = values[i];
    w[j] = weight[i];
    while (j > 0 && x[j - 1]>x[j]) {
      x[j - 1] = (w[j] * x[j] + w[j - 1] * x[j - 1]) / (w[j] + w[j - 1]);
      w[j - 1] += w[j];
      j -= 1;
    }
    s[j + 1] = i + 1;
  }
  
  vec ww(values.size(), fill::zeros);
  vec xx(values.size(), fill::zeros);
  for (unsigned k = 0; k < j + 1; k++) {
    for (unsigned i = s[k]; i < s[k + 1]; i++) {
      ww[i] = w[k];
      xx[i] = x[k];
    }
  }
  
  if(decreasing){
    xx = reverse(xx);
  }
  
  return xx;
}