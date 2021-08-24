#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
arma::vec basicsolution1cpt(arma::vec tm, double kR, double k10, double v1, double init) {
  arma::vec cons = (kR/k10*(1.0-exp(-tm*k10)) + init*v1 * exp(-tm*k10)) / v1;
  return cons;
}


// [[Rcpp::export]]
arma::mat pksol1cpt(arma::vec& tms, arma::vec& pars, arma::vec& begin,
                    arma::vec& end, arma::vec& infs, double init) {

  double k_10 = pars[0];
  double v_1  = pars[1];

  // append values 0 and infinity to begin/end sequences
  // tms_all contains all begin and end times, supplemented with 0 and infinity
  // arma::vec vec_0(1, arma::fill::zeros);
  arma::vec vec_inf(1, arma::fill::zeros);
  vec_inf.fill(std::numeric_limits<int>::max());
  // arma::vec m1 = arma::join_cols(vec_0, begin); //
  // arma::vec m2 = arma::join_cols(m1, end);
  arma::vec m2 = arma::join_cols(begin, end);
  arma::vec tms_all = arma::join_cols(m2, vec_inf);

  // sort and subset to unique begin/end times
  arma::vec tms_all_sort = arma::sort(tms_all);
  arma::uvec unique_order = arma::find_unique(tms_all_sort);
  arma::vec prd = tms_all_sort.elem(unique_order);

  // Initialize matrices for starting and evaluated concentrations
  arma::vec cons(tms.size(), arma::fill::zeros);
  arma::vec inits(prd.size()+1, arma::fill::zeros);

  // Set initial concentrations
  inits(0) = init;

  int P = prd.size()-1;
  // Create vector of infusion rates
  arma::vec infi(prd.size(), arma::fill::zeros);
  for(int i = 0; i < P; i++) {
    for(int j = 0; j < begin.size(); j++){
      if(prd[i] >= begin[j] && prd[i] < end[j])
        infi[i] += infs[j];
    }
  }

  arma::vec tmsi;
  for(int ii = 0; ii < P; ii++) {

    double eps = 1e-6;
    // Identify indices and times to evaluate
    arma::uvec ix_tms;
    if(ii == 0){
      // if first evaluation, include zero if specified
      ix_tms = find(tms >= (prd(ii)-eps) && tms <= (prd(ii+1)+eps));
    } else{
      // Include upper bound to use for initial values
      ix_tms = find(tms > prd(ii) && tms <= (prd(ii+1)+eps));
    }

    arma::vec tmsix = tms.elem(ix_tms); // times specified
    arma::vec val_prd(1, arma::fill::zeros); // time corresponding to end of period (prd(ii+1))
    val_prd.fill(prd(ii+1));
    arma::vec tmsi = arma::join_cols(tmsix, val_prd); // joint set of times

    // Evaluate piece-wise PK solution
    arma::vec consi = basicsolution1cpt(tmsi - prd[ii], infi[ii], k_10, v_1, inits(ii));

    // Update initial values
    inits(ii+1) = consi(consi.n_rows-1);

    if(ix_tms.size() > 0){
      // Store values
      cons.subvec(min(ix_tms),max(ix_tms)) = consi.subvec(0,consi.n_rows-2); // Store concentrations in rows
    }
    // std::cout << "cons: " << cons;

  }

  // Replace any negative values
  arma::uvec neg_ids = find(cons < 0.0);
  cons.elem(neg_ids).fill(0.0);

  return std::move(cons);
}





// [[Rcpp::export]]
arma::mat basicsolution2cpt(arma::vec tm, double kR, double k10, double k12, double k21, double v1,
                            double v2, arma::vec c0) {

  double k20 = 0.0;
  double E1 = k10+k12;
  double E2 = k21+k20;
  double E1E2 = E1+E2;

  double lambda1 = 0.5*((E1E2)+sqrt(pow(E1E2,2.0)-4.0*(E1*E2-k12*k21)));
  double lambda2 = 0.5*((E1E2)-sqrt(pow(E1E2,2.0)-4.0*(E1*E2-k12*k21)));

  double A1last = c0[0]*v1;
  double A2last = c0[1]*v2;
  double Doserate = kR;

  arma::vec A1term1 = (((A1last*E2+Doserate+A2last*k21)-A1last*lambda1)*exp(-tm*lambda1)-((A1last*E2+Doserate+A2last*k21)-A1last*lambda2)*exp(-tm*lambda2))/(lambda2-lambda1);
  arma::vec A1term2 = Doserate*E2*(1/(lambda1*lambda2)+exp(-tm*lambda1)/(lambda1*(lambda1-lambda2))-exp(-tm*lambda2)/(lambda2*(lambda1-lambda2)));
  arma::vec a_1 = A1term1+A1term2;

  arma::vec A2term1 = (((A2last*E1+A1last*k12)-A2last*lambda1)*exp(-tm*lambda1)-((A2last*E1+A1last*k12)-A2last*lambda2)*exp(-tm*lambda2))/(lambda2-lambda1);
  arma::vec A2term2 = Doserate*k12*(1/(lambda1*lambda2)+exp(-tm*lambda1)/(lambda1*(lambda1-lambda2))-exp(-tm*lambda2)/(lambda2*(lambda1-lambda2)));
  arma::vec a_2 = A2term1+A2term2;

  int nobs = tm.size();

  arma::mat cons(2,nobs);
  cons(0,arma::span(0,nobs-1)) = a_1.t()/v1;
  cons(1,arma::span(0,nobs-1)) = a_2.t()/v2;

  return cons;
}


// [[Rcpp::export]]
arma::mat pksol2cpt(arma::vec& tms, arma::vec& pars, arma::vec& begin,
                    arma::vec& end, arma::vec& infs, arma::vec& init) {

  double k_10 = pars[0];
  double k_12 = pars[1];
  double k_21 = pars[2];
  double v_1  = pars[3];
  double v_2  = pars[4];

  // append values 0 and infinity to begin/end sequences
  // arma::vec vec_0(1, arma::fill::zeros);
  arma::vec vec_inf(1, arma::fill::zeros);
  vec_inf.fill(std::numeric_limits<int>::max());
  // arma::vec m1 = arma::join_cols(vec_0, begin);
  // arma::vec m2 = arma::join_cols(m1, end);
  arma::vec m2 = arma::join_cols(begin, end);
  arma::vec tms_all = arma::join_cols(m2, vec_inf);

  // sort and subset to unique times
  arma::vec tms_all_sort = arma::sort(tms_all);
  arma::uvec unique_order = arma::find_unique(tms_all_sort);
  arma::vec prd = tms_all_sort.elem(unique_order);

  // Initialize matrices for starting and evaluated concentrations
  arma::mat cons(2,tms.size(), arma::fill::zeros);
  arma::mat inits(2,prd.size()+1, arma::fill::zeros);

  // Set initial concentrations
  inits(arma::span(0,1),0) = init;

  int P = prd.size()-1;
  // Create vector of infusion rates
  arma::vec infi(prd.size(), arma::fill::zeros);
  for(int i = 0; i < P; i++) {
    for(int j = 0; j < begin.size(); j++){
      if(prd[i] >= begin[j] && prd[i] < end[j])
        infi[i] += infs[j];
    }
  }

  arma::vec tmsi;
  for(int ii = 0; ii < P; ii++) {

    double eps = 1e-6;
    // Identify indices and times to evaluate
    arma::uvec ix_tms;
    if(ii == 0){
      // if first evaluation, include zero if specified
      ix_tms = find(tms >= (prd(ii)-eps) && tms <= (prd(ii+1)+eps));
    } else{
      // Include upper bound to use for initial values
      ix_tms = find(tms > prd(ii) && tms <= (prd(ii+1)+eps));
    }

    arma::vec tmsix = tms.elem(ix_tms); // times specified
    arma::vec val_prd(1, arma::fill::zeros); // time corresponding to end of period (prd(ii+1))
    val_prd.fill(prd(ii+1));
    arma::vec tmsi = arma::join_cols(tmsix, val_prd); // joint set of times

    // Evaluate piecewise PK solution
    arma::mat consi = basicsolution2cpt(tmsi - prd[ii], infi[ii], k_10, k_12, k_21, v_1, v_2, inits(arma::span(0,1),ii));

    // Update initial values
    inits(arma::span(0,1),ii+1) = consi(arma::span(0,1),consi.n_cols-1);

    if(ix_tms.size() > 0){
      // Store values
      cons.submat(0, min(ix_tms), 1, max(ix_tms)) = consi.submat(0,0,1,consi.n_cols-2); // Store concentrations in rows
    }
    // std::cout << "cons: " << cons;

  }

  // Replace any negative values
  arma::umat neg_ids = find(cons < 0.0);
  cons.elem(neg_ids).fill(0.0);

  return cons;
}



// [[Rcpp::export]]
arma::mat basicsolution3cpt(arma::vec tm, double kR, double k10, double k12, double k21, double k13, double k31, double v1,
                            double v2, double v3, arma::vec c0) {

  double k20 = 0.0;
  double k30 = 0.0;

  double E1 = k10+k12+k13;
  double E2 = k21+k20;
  double E3 = k31+k30;
  double a = E1+E2+E3;
  double b = E1*E2+E3*(E1+E2)-k12*k21-k13*k31;
  double c = E1*E2*E3-E3*k12*k21-E2*k13*k31;
  double m = (3.0*b - pow(a, 2.0))/3.0;
  double n = (2.0*pow(a, 3.0) - 9.0*a*b + 27.0*c)/27.0;
  double Q = (pow(n,2.0))/4.0 + (pow(m,3.0))/27.0;
  double Qmx = max(Rcpp::NumericVector::create(0.0,-1.0*Q));
  double alpha = sqrt(Qmx);
  double beta  = -1.0*n/2.0;
  double gamma =  sqrt(pow(beta,2.0)+pow(alpha,2.0));
  double theta =  atan2(alpha,beta);

  double lambda1 = a/3.0 + pow(gamma,1.0/3.0)*(cos(theta/3.0) + sqrt(3.0)*sin(theta/3.0));
  double lambda2 = a/3.0 + pow(gamma,1.0/3.0)*(cos(theta/3.0) - sqrt(3.0)*sin(theta/3.0));
  double lambda3 = a/3.0 -(2.0*pow(gamma,1.0/3.0)*cos(theta/3.0));

  double A1last = c0[0]*v1;
  double A2last = c0[1]*v2;
  double A3last = c0[2]*v3;
  double Doserate = kR;

  double B = A2last*k21+A3last*k31;
  double C = E3*A2last*k21+E2*A3last*k31;
  double I = A1last*k12*E3-A2last*k13*k31+A3last*k12*k31;
  double J = A1last*k13*E2+A2last*k13*k21-A3last*k12*k21;

  arma::vec A1term1 = A1last*(exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
  arma::vec A1term2 = exp(-tm*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2));
  arma::vec A1term3 = Doserate*((E2*E3)/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)));
  arma::vec a_1 = A1term1+A1term2+A1term3;

  arma::vec A2term1 = A2last*(exp(-tm*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
  arma::vec A2term2 = exp(-tm*lambda1)*(I-A1last*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(A1last*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(A1last*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2));
  arma::vec A2term3 = Doserate*k12*(E3/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)));
  arma::vec a_2 = A2term1+A2term2+A2term3;

  arma::vec A3term1 = A3last*(exp(-tm*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
  arma::vec A3term2 = exp(-tm*lambda1)*(J-A1last*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(A1last*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(A1last*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2));
  arma::vec A3term3 = Doserate*k13*(E2/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E2-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E2-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E2-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)));
  arma::vec a_3 = A3term1+A3term2+A3term3;

  int nobs = tm.size();

  arma::mat cons(3,nobs);
  cons(0,arma::span(0,nobs-1)) = a_1.t()/v1;
  cons(1,arma::span(0,nobs-1)) = a_2.t()/v2;
  cons(2,arma::span(0,nobs-1)) = a_3.t()/v3;

  return cons;
}


// [[Rcpp::export]]
arma::mat pksol3cpt(arma::vec& tms, arma::vec& pars, arma::vec& begin,
                    arma::vec& end, arma::vec& infs, arma::vec& init) {

  double k_10 = pars[0];
  double k_12 = pars[1];
  double k_21 = pars[2];
  double k_13 = pars[3];
  double k_31 = pars[4];
  double v_1  = pars[5];
  double v_2  = pars[6];
  double v_3  = pars[7];

  // append values 0 and infinity to begin/end sequences
  // arma::vec vec_0(1, arma::fill::zeros);
  arma::vec vec_inf(1, arma::fill::zeros);
  vec_inf.fill(std::numeric_limits<int>::max());
  // arma::vec m1 = arma::join_cols(vec_0, begin);
  // arma::vec m2 = arma::join_cols(m1, end);
  arma::vec m2 = arma::join_cols(begin, end);
  arma::vec tms_all = arma::join_cols(m2, vec_inf);

  // sort and subset to unique times
  arma::vec tms_all_sort = arma::sort(tms_all);
  arma::uvec unique_order = arma::find_unique(tms_all_sort);
  arma::vec prd = tms_all_sort.elem(unique_order);

  // Initialize matrices for starting and evaluated concentrations
  arma::mat cons(3,tms.size(), arma::fill::zeros);
  arma::mat inits(3,prd.size()+1, arma::fill::zeros);

  // Set initial concentrations
  inits(arma::span(0,2),0) = init;

  int P = prd.size()-1;
  // Create vector of infusion rates
  arma::vec infi(prd.size(), arma::fill::zeros);
  for(int i = 0; i < P; i++) {
    for(int j = 0; j < begin.size(); j++){
      if(prd[i] >= begin[j] && prd[i] < end[j])
        infi[i] += infs[j];
    }
  }

  arma::vec tmsi;
  for(int ii = 0; ii < P; ii++) {

    double eps = 1e-6;
    // Identify indices and times to evaluate
    arma::uvec ix_tms;
    if(ii == 0){
      // if first evaluation, include zero if specified
      ix_tms = find(tms >= (prd(ii)-eps) && tms <= (prd(ii+1)+eps));
    } else{
      // Include upper bound to use for initial values
      ix_tms = find(tms > prd(ii) && tms <= (prd(ii+1)+eps));
    }

    arma::vec tmsix = tms.elem(ix_tms); // times specified
    arma::vec val_prd(1, arma::fill::zeros); // time corresponding to end of period (prd(ii+1))
    val_prd.fill(prd(ii+1));
    arma::vec tmsi = arma::join_cols(tmsix, val_prd); // joint set of times

    // Evaluate piecewise PK solution
    arma::mat consi = basicsolution3cpt(tmsi - prd[ii], infi[ii], k_10, k_12, k_21,k_13, k_31, v_1, v_2, v_3, inits(arma::span(0,2),ii));

    // Update initial values
    inits(arma::span(0,2),ii+1) = consi(arma::span(0,2),consi.n_cols-1);

    if(ix_tms.size() > 0){
      // Store values
      cons.submat(0, min(ix_tms), 2, max(ix_tms)) = consi.submat(0,0,2,consi.n_cols-2); // Store concentrations in rows
    }
    // std::cout << "cons: " << cons;

  }

  // Replace any negative values
  arma::umat neg_ids = find(cons < 0.0);
  cons.elem(neg_ids).fill(0.0);

  return cons;
}


// [[Rcpp::export]]
arma::mat basicsolution3cptm(arma::vec tm, double kR, double k10, double k12, double k21, double k13, double k31, double v1,
                                           double v2, double v3, double ke0, arma::vec c0) {

  double kme = ke0;
  double km = ke0 / 100000.0;
  double v4 = v1 / 100000.0;
  double k20 = 0.0;
  double k30 = 0.0;
  double E1 = k10+k12+k13+km;
  double E2 = k21+k20;
  double E3 = k31+k30;
  double a = E1+E2+E3;
  double b = E1*E2+E3*(E1+E2)-k12*k21-k13*k31;
  double c = E1*E2*E3-E3*k12*k21-E2*k13*k31;
  double m = (3.0*b - pow(a, 2.0))/3.0;
  double n = (2.0*pow(a, 3.0) - 9.0*a*b + 27.0*c)/27.0;
  double Q = (pow(n,2.0))/4.0 + (pow(m,3.0))/27.0;
  double Qmx = max(Rcpp::NumericVector::create(0.0,-1.0*Q));
  double alpha = sqrt(Qmx);
  double beta  = -1.0*n/2.0;
  double gamma =  sqrt(pow(beta,2.0)+pow(alpha,2.0));
  double theta =  atan2(alpha,beta);

  double lambda1 = a/3.0 + pow(gamma,1.0/3.0)*(cos(theta/3.0) + sqrt(3.0)*sin(theta/3.0));
  double lambda2 = a/3.0 + pow(gamma,1.0/3.0)*(cos(theta/3.0) - sqrt(3.0)*sin(theta/3.0));
  double lambda3 = a/3.0 -(2.0*pow(gamma,1.0/3.0)*cos(theta/3.0));

  double A1last = c0[0]*v1;
  double A2last = c0[1]*v2;
  double A3last = c0[2]*v3;
  double Amlast = c0[3]*v4;
  double Doserate = kR;

  double B = A2last*k21+A3last*k31;
  double C = E3*A2last*k21+E2*A3last*k31;
  double I = A1last*k12*E3-A2last*k13*k31+A3last*k12*k31;
  double J = A1last*k13*E2+A2last*k13*k21-A3last*k12*k21;

  arma::vec A1term1 = A1last*(exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
  arma::vec A1term2 = exp(-tm*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2));
  arma::vec A1term3 = Doserate*((E2*E3)/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)));
  arma::vec a_1 = A1term1+A1term2+A1term3;

  arma::vec A2term1 = A2last*(exp(-tm*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
  arma::vec A2term2 = exp(-tm*lambda1)*(I-A1last*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(A1last*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(A1last*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2));
  arma::vec A2term3 = Doserate*k12*(E3/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)));
  arma::vec a_2 = A2term1+A2term2+A2term3;

  arma::vec A3term1 = A3last*(exp(-tm*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-tm*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
  arma::vec A3term2 = exp(-tm*lambda1)*(J-A1last*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-tm*lambda2)*(A1last*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-tm*lambda3)*(A1last*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2));
  arma::vec A3term3 = Doserate*k13*(E2/(lambda1*lambda2*lambda3)-exp(-tm*lambda1)*(E2-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E2-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E2-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)));
  arma::vec a_3 = A3term1+A3term2+A3term3;

  arma::vec Amterm1 = Amlast*exp(-tm*kme) +km*A1last*(exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(kme-lambda1))+exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/((kme-lambda2)*(lambda1-lambda2)*(lambda3-lambda2))+exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/((kme-lambda3)*(lambda1-lambda3)*(lambda2-lambda3))+exp(-tm*kme)*(E2-kme)*(E3-kme)/((lambda1-kme)*(lambda2-kme)*(lambda3-kme)));
  arma::vec Amterm2 = km*(exp(-tm*lambda1)*(B*lambda1-C)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-kme))+exp(-tm*lambda2)*(C-B*lambda2)/((lambda1-lambda2)*(lambda2-lambda3)*(lambda2-kme))+exp(-tm*lambda3)*(C-B*lambda3)/((lambda1-lambda3)*(lambda3-lambda2)*(lambda3-kme))-exp(-tm*kme)*(B*kme-C)/((lambda1-kme)*(kme-lambda2)*(kme-lambda3)));
  arma::vec Amterm3 = km*Doserate*((E2*E3)/(lambda1*lambda2*lambda3*kme)-exp(-tm*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(kme-lambda1)*(lambda2-lambda1)*(lambda3-lambda1))-exp(-tm*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(kme-lambda2)*(lambda1-lambda2)*(lambda3-lambda2))-exp(-tm*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(kme-lambda3)*(lambda1-lambda3)*(lambda2-lambda3))-exp(-tm*kme)*(E2-kme)*(E3-kme)/(kme*(lambda1-kme)*(lambda2-kme)*(lambda3-kme)));
  arma::vec a_4 = Amterm1+Amterm2+Amterm3;

  int nobs = tm.size();

  arma::mat cons(4,nobs);
  cons(0,arma::span(0,nobs-1)) = a_1.t()/v1;
  cons(1,arma::span(0,nobs-1)) = a_2.t()/v2;
  cons(2,arma::span(0,nobs-1)) = a_3.t()/v3;
  cons(3,arma::span(0,nobs-1)) = a_4.t()/v4;

  return cons;
}


// [[Rcpp::export]]
arma::mat pksol3cptm(arma::vec& tms, arma::vec& pars, arma::vec& begin,
                                     arma::vec& end, arma::vec& infs, arma::vec& init) {

  double k_10 = pars[0];
  double k_12 = pars[1];
  double k_21 = pars[2];
  double k_13 = pars[3];
  double k_31 = pars[4];
  double v_1  = pars[5];
  double v_2  = pars[6];
  double v_3  = pars[7];
  double k_e0 = pars[8];

  // append values 0 and infinity to begin/end sequences
  // arma::vec vec_0(1, arma::fill::zeros);
  arma::vec vec_inf(1, arma::fill::zeros);
  vec_inf.fill(std::numeric_limits<int>::max());
  // arma::vec m1 = arma::join_cols(vec_0, begin);
  // arma::vec m2 = arma::join_cols(m1, end);
  arma::vec m2 = arma::join_cols(begin, end);
  arma::vec tms_all = arma::join_cols(m2, vec_inf);

  // sort and subset to unique times
  arma::vec tms_all_sort = arma::sort(tms_all);
  arma::uvec unique_order = arma::find_unique(tms_all_sort);
  arma::vec prd = tms_all_sort.elem(unique_order);

  // Initialize matrices for starting and evaluated concentrations
  arma::mat cons(4,tms.size(), arma::fill::zeros);
  arma::mat inits(4,prd.size()+1, arma::fill::zeros);

  // Set initial concentrations
  inits(arma::span(0,3),0) = init;

  int P = prd.size()-1;
  // Create vector of infusion rates
  arma::vec infi(prd.size(), arma::fill::zeros);
  for(int i = 0; i < P; i++) {
    for(int j = 0; j < begin.size(); j++){
      if(prd[i] >= begin[j] && prd[i] < end[j])
        infi[i] += infs[j];
    }
  }

  arma::vec tmsi;
  for(int ii = 0; ii < P; ii++) {

    double eps = 1e-6;
    // Identify indices and times to evaluate
    arma::uvec ix_tms;
    if(ii == 0){
      // if first evaluation, include zero if specified
      ix_tms = find(tms >= (prd(ii)-eps) && tms <= (prd(ii+1)+eps));
    } else{
      // Include upper bound to use for initial values
      ix_tms = find(tms > prd(ii) && tms <= (prd(ii+1)+eps));
    }

    arma::vec tmsix = tms.elem(ix_tms); // times specified
    arma::vec val_prd(1, arma::fill::zeros); // time corresponding to end of period (prd(ii+1))
    val_prd.fill(prd(ii+1));
    arma::vec tmsi = arma::join_cols(tmsix, val_prd); // joint set of times

    // Evaluate piecewise PK solution
    arma::mat consi = basicsolution3cptm(tmsi - prd[ii], infi[ii], k_10, k_12, k_21,k_13, k_31, v_1, v_2, v_3, k_e0, inits(arma::span(0,3),ii));

    // Update initial values
    inits(arma::span(0,3),ii+1) = consi(arma::span(0,3),consi.n_cols-1);

    if(ix_tms.size() > 0){
      // Store values
      cons.submat(0, min(ix_tms), 3, max(ix_tms)) = consi.submat(0,0,3,consi.n_cols-2); // Store concentrations in rows
    }
    // std::cout << "cons: " << cons;

  }

  // Replace any negative values
  arma::umat neg_ids = find(cons < 0.0);
  cons.elem(neg_ids).fill(0.0);

  return cons;
}



