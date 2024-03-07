// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

#include <omp.h>
//[[Rcpp::plugins(openmp)]]

using namespace std;
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
double rinvGauss(double nu, double lambda)
{
  double y = R::rchisq(1);
  double x;
  double u = 4 * lambda / (nu * y);
  if(u > 1e-11)
  {
    x = nu * (1 - 2 / (1 + sqrt(1 + u)));
  } else {
    x = lambda / y;
  }
  double z = R::runif(0,1);
  if(z <= nu / (nu + x))
    return(x);
  else
    return(pow(nu, 2)/x);
}

double propose_lognormal(double x, double stepsize = 2)
{
  x *= exp(R::rnorm(0, stepsize));
  return(x);
}

double lpexp(double x, double y)
{
  if(x > y)
  {
    return(x + log1p(exp(y-x)));
  } else {
    return(y + log1p(exp(x-y)));
  }
}

double lmexp(double x, double y)
{
  return(x + log1p(-exp(y-x)));
}

// [[Rcpp::export]]
void h2D2_pretrain(S4 h2D2, 
                   const Eigen::SparseMatrix<double> W,
                   const Eigen::VectorXd NbetaHat,
                   const int pre_mcmc_n = 200,
                   const int pre_use = 100,
                   const double pre_p = 0.05,
                   const int pre_maxiter = 10,
                   const double stepsize = 2, 
                   const unsigned int seed = 428)
{
  Function set_seed("set.seed");
  set_seed(seed);
  
  // Data
  const int M = h2D2.slot("M");
  const double N = h2D2.slot("N");
  
  Environment methods = Environment::namespace_env("methods");
  Function AS = methods["as"];
  const SparseMatrix<double> R = as<SparseMatrix<double>>(
    AS(h2D2.slot("R"), _["Class"] = "dgCMatrix"));
  const SparseMatrix<double> LD_pairs = h2D2.slot("LD_pairs");
  
  // Hyper parameters
  const NumericVector a = h2D2.slot("a");
  const double asum = sum(a);
  const double sqrt2 = sqrt(2);
  double bmax = asum * (1-1e-6) * 1e6;
  double bmin = asum * (1-0.1) * 10;
  double b = asum * (1-1e-4) * 1e4; // initial b
  double h2_est = 1e-4;
  
  // Initial samples
  int i = 0;
  int ii;
  
  double h2, log_res;
  VectorXd beta, log_sigma2, psi, keep_h2_beta;
  
  beta = VectorXd::Zero(M);
  log_sigma2 = VectorXd::Constant(M,log(1e-4/M));
  psi = VectorXd::Constant(M,1);
  h2 = 1e-4;
  
  keep_h2_beta = VectorXd::Zero(pre_use);
  
  // mcmc
  int j, j1, k;
  int j2 = 0;
  log_res = log(1 - h2);
  double sigma2_j, log_sigma2_j_new, sigma2_j_new, tsigma2_j, tsigma2_j_new;
  double u_j, nu_j, acc_prob, rss, ch, r, u;
  double W_12, p_21;
  double p_12 = 0;
  NumericVector order_sigma2;
  
  NumericVector z = h2D2.slot("z");
  Function order("order");
  NumericVector order_SNP = order(abs(z), _["decreasing"] = true);
  order_SNP = order_SNP - 1;
  
  Environment stats = Environment::namespace_env("stats");
  Function t_test = stats["t.test"];
  List test;
  double pval, mean;
  
  bool dir = true;
  
  while(i < pre_maxiter)
  {
    ii = 0;
    while(ii < pre_mcmc_n)
    {
      for(k = 0; k < M; ++k)
      {
        j = order_SNP(k);
        
        //1. sigma2
        sigma2_j = exp(log_sigma2(j));
        log_res = lpexp(log_res, log_sigma2(j));
        
        log_sigma2_j_new = propose_lognormal(log_sigma2(j) - log_res, stepsize) + log_res;
        acc_prob = (log_sigma2_j_new - log_res) / (log_sigma2(j) - log_res);
        
        sigma2_j_new = exp(log_sigma2_j_new);
        
        tsigma2_j = 1 / (N + 1 / (sigma2_j * psi(j)));
        tsigma2_j_new = 1 / (N + 1 / (sigma2_j_new * psi(j)));
        
        u_j = NbetaHat(j);
        for(SparseMatrix<double>::InnerIterator it(W,j); it; ++it)
        {
          u_j -= it.value() * beta(it.row());
        }
        
        acc_prob *= sqrt((N * sigma2_j * psi(j) + 1) /
          (N * sigma2_j_new * psi(j) + 1)) *
            exp(pow(u_j, 2) * (tsigma2_j_new - tsigma2_j) / 2) *
            exp(a(j) * (log_sigma2_j_new - log_sigma2(j))) *
            exp((b-1) * (lmexp(log_res, log_sigma2_j_new) - 
            lmexp(log_res, log_sigma2(j))));
        if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
        {
          log_sigma2(j) = log_sigma2_j_new;
          sigma2_j = sigma2_j_new;
          tsigma2_j = tsigma2_j_new;
        }
        log_res = lmexp(log_res, log_sigma2(j));
        
        //2. beta
        beta(j) = R::rnorm(tsigma2_j * u_j, sqrt(tsigma2_j));
        
        //3. psi
        if(beta(j) == 0)
        {
          psi(j) = 1;
        } else {
          nu_j = sqrt2 * exp(log_sigma2(j)/2) / abs(beta(j));
          psi(j) = 1 / rinvGauss(nu_j, 2);
        }
      }
      
      if(ii >= pre_mcmc_n - pre_use)
      {
        keep_h2_beta(ii + pre_use - pre_mcmc_n) = 0;
        for(j = 0; j < M; ++j)
        {
          for(SparseMatrix<double>::InnerIterator it(R,j); it; ++it)
          {
            k = it.row();
            if(k > j)
            {
              keep_h2_beta(ii + pre_use - pre_mcmc_n) += 
                2 * it.value() * beta(k) * beta(j);
            }
          }
          keep_h2_beta(ii + pre_use - pre_mcmc_n) += pow(beta(j), 2);
        }
      }
      
      ++ii;
      
      if(((ii + 1) % 10 == 0) && (ii < pre_mcmc_n))
      {
        //Choose the index for switching.
        order_sigma2 = order(log_sigma2, _["decreasing"] = true);
        order_sigma2 = order_sigma2 - 1;
        k = 10 + R::rgeom(0.3);
        if(k > M) k = M - 1;
        
        for(j=0; j<k; ++j)
        {
          j1 = order_sigma2(j); //index of SNP 1
          if(LD_pairs.col(j1).nonZeros() >= 1)
          {
            r = R::runif(0,1);
            u = 0;
            for(SparseMatrix<double>::InnerIterator it(LD_pairs,j1); it; ++it)
            {
              u += it.value();
              if(u > r)
              {
                j2 = it.row();
                p_12 = it.value();
                break;
              }
            }
            p_21 = LD_pairs.coeff(j1,j2);
            
            rss = 0;
            W_12 = W.coeff(j1,j2);
            if(W_12 > 0)
            {
              for(SparseMatrix<double>::InnerIterator it(W,j1); it; ++it)
              {
                rss += it.value() * beta(it.row());
              }
              for(SparseMatrix<double>::InnerIterator it(W,j2); it; ++it)
              {
                rss -= it.value() * beta(it.row());
              }
              rss += W_12 * (beta(j1) - beta(j2));
              rss *= (beta(j1) - beta(j2));
              
              rss += (NbetaHat(j1) - NbetaHat(j2)) * (beta(j2) - beta(j1));
            } else {
              for(SparseMatrix<double>::InnerIterator it(W,j1); it; ++it)
              {
                rss += it.value() * beta(it.row());
              }
              for(SparseMatrix<double>::InnerIterator it(W,j2); it; ++it)
              {
                rss += it.value() * beta(it.row());
              }
              rss -= W_12 * (beta(j1) + beta(j2));
              rss *= (beta(j1) + beta(j2));
              
              rss += -(NbetaHat(j1) + NbetaHat(j2)) * (beta(j1) + beta(j2));
            }
            
            acc_prob = exp(rss) *
              exp((a(j1) - a(j2)) * (log_sigma2(j2) - log_sigma2(j1))) *
              p_21 / p_12;
            if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
            {
              if(W_12 > 0)
              {
                ch = beta(j1);
                beta(j1) = beta(j2);
                beta(j2) = ch;
              } else {
                ch = beta(j1);
                beta(j1) = -beta(j2);
                beta(j2) = -ch;
              }
              
              ch = log_sigma2(j1);
              log_sigma2(j1) = log_sigma2(j2);
              log_sigma2(j2) = ch;
              
              ch = psi(j1);
              psi(j1) = psi(j2);
              psi(j2) = ch;
            }
          }
        }
        
        if(ii >= pre_mcmc_n - pre_use)
        {
          keep_h2_beta(ii + pre_use - pre_mcmc_n) = 0;
          for(j = 0; j < M; ++j)
          {
            for(SparseMatrix<double>::InnerIterator it(R,j); it; ++it)
            {
              k = it.row();
              if(k > j)
              {
                keep_h2_beta(ii + pre_use - pre_mcmc_n) += 
                  2 * it.value() * beta(k) * beta(j);
              }
            }
            keep_h2_beta(ii + pre_use - pre_mcmc_n) += pow(beta(j), 2);
          }
        }
        
        ++ii;
      }
    }
    
    test = t_test(_["x"] = keep_h2_beta, _["mu"] = h2_est);
    pval = test["p.value"];
    mean = test["estimate"];
    if(i == 0) dir = mean > h2_est;
    
    cout << "Pretrain step " << i << ", b=" << b << ", h2_est=" <<
      h2_est << ", mean=" << mean << ", p=" << pval << 
        "." << endl;
    
    if((pval >= pre_p) || (mean > h2_est) != dir)
    {
      break;
    } else {
      b = asum * (1 - mean) / mean;
      h2_est = mean;
    }
    
    ++i;
  }
  
  cout << "End pretrain. b=" << b << "." << endl;
  if(b > bmax) b = bmax;
  if(b < bmin) b = bmin;
  
  h2D2.slot("b") = b;
}

// [[Rcpp::export]]
void h2D2_sampling(S4 h2D2, 
                   const Eigen::SparseMatrix<double> W,
                   const Eigen::VectorXd NbetaHat,
                   const int mcmc_n = 100, 
                   const int thin = 1, 
                   const double stepsize = 2, 
                   const unsigned int seed = 428)
{
  Function set_seed("set.seed");
  set_seed(seed);
  
  // Data
  const int M = h2D2.slot("M");
  const double N = h2D2.slot("N");
  
  Environment methods = Environment::namespace_env("methods");
  Function AS = methods["as"];
  const SparseMatrix<double> R = as<SparseMatrix<double>>(
    AS(h2D2.slot("R"), _["Class"] = "dgCMatrix"));
  const SparseMatrix<double> LD_pairs = h2D2.slot("LD_pairs");
  
  // Hyper parameters
  const NumericVector a = h2D2.slot("a");
  const double b = h2D2.slot("b");
  const double sqrt2 = sqrt(2);
  
  // Initial samples
  int i = 0;
  int ii;
  List mcmc_samples = h2D2.slot("mcmc_samples");
  
  int n_samples, n_burnin;
  double h2, log_res;
  VectorXd beta, log_sigma2, psi;
  VectorXd keep_h2, keep_h2_beta;
  MatrixXd keep_beta, keep_log_sigma2, keep_psi;
  
  if (mcmc_samples.size() == 0)
  {
    n_samples = 0;
    n_burnin = 0;
    
    beta = VectorXd::Zero(M);
    log_sigma2 = VectorXd::Constant(M,log(1e-4/M));
    psi = VectorXd::Constant(M,1);
    h2 = 1e-4;
    
    keep_log_sigma2 = MatrixXd::Zero(mcmc_n, M);
    keep_log_sigma2.row(0) = log_sigma2;
    
    keep_beta = MatrixXd::Zero(mcmc_n, M);
    keep_beta.row(0) = beta;
    
    keep_psi = MatrixXd::Zero(mcmc_n, M);
    keep_psi.row(0) = psi;
    
    keep_h2 = VectorXd::Zero(mcmc_n);
    keep_h2(0) = h2;
    
    keep_h2_beta = VectorXd::Zero(mcmc_n);
    keep_h2_beta(0) = 0;
    
    ++i;
  } else {
    n_samples = mcmc_samples["n_samples"];
    n_burnin = mcmc_samples["n_burnin"];
    
    keep_log_sigma2 = mcmc_samples["log_sigma2"];
    log_sigma2 = keep_log_sigma2.row(n_samples - 1);
    keep_log_sigma2.conservativeResize(n_samples + mcmc_n, M);
    
    keep_beta = mcmc_samples["beta"];
    beta = keep_beta.row(n_samples - 1);
    keep_beta.conservativeResize(n_samples + mcmc_n, M);
    
    keep_psi = mcmc_samples["psi"];
    psi = keep_psi.row(n_samples - 1);
    keep_psi.conservativeResize(n_samples + mcmc_n, M);
    
    keep_h2 = mcmc_samples["h2"];
    h2 = keep_h2(n_samples - 1);
    keep_h2.conservativeResize(n_samples + mcmc_n);
    
    keep_h2_beta = mcmc_samples["h2_beta"];
    keep_h2_beta.conservativeResize(n_samples + mcmc_n);
  }
  
  // mcmc
  int j, j1, k;
  int j2 = 0;
  log_res = log(1 - h2);
  double sigma2_j, log_sigma2_j_new, sigma2_j_new, tsigma2_j, tsigma2_j_new;
  double u_j, nu_j, acc_prob, rss, ch, r, u;
  double W_12, p_21;
  double p_12 = 0;
  NumericVector order_sigma2;
  
  NumericVector z = h2D2.slot("z");
  Function order("order");
  NumericVector order_SNP = order(abs(z), _["decreasing"] = true);
  order_SNP = order_SNP - 1;
  
  while(i < mcmc_n)
  {
    for(ii = 0; ii < thin; ++ii)
    {
      for(k = 0; k < M; ++k)
      {
        j = order_SNP(k);
        
        //1. sigma2
        sigma2_j = exp(log_sigma2(j));
        log_res = lpexp(log_res, log_sigma2(j));
        
        log_sigma2_j_new = propose_lognormal(log_sigma2(j) - log_res, stepsize) + log_res;
        acc_prob = (log_sigma2_j_new - log_res) / (log_sigma2(j) - log_res);
        
        sigma2_j_new = exp(log_sigma2_j_new);
        
        tsigma2_j = 1 / (N + 1 / (sigma2_j * psi(j)));
        tsigma2_j_new = 1 / (N + 1 / (sigma2_j_new * psi(j)));
        
        u_j = NbetaHat(j);
        for(SparseMatrix<double>::InnerIterator it(W,j); it; ++it)
        {
          u_j -= it.value() * beta(it.row());
        }
        
        acc_prob *= sqrt((N * sigma2_j * psi(j) + 1) /
          (N * sigma2_j_new * psi(j) + 1)) *
            exp(pow(u_j, 2) * (tsigma2_j_new - tsigma2_j) / 2) *
            exp(a(j) * (log_sigma2_j_new - log_sigma2(j))) *
            exp((b-1) * (lmexp(log_res, log_sigma2_j_new) - 
            lmexp(log_res, log_sigma2(j))));
        if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
        {
          log_sigma2(j) = log_sigma2_j_new;
          // sigma2_j = sigma2_j_new;
          tsigma2_j = tsigma2_j_new;
        }
        log_res = lmexp(log_res, log_sigma2(j));
        
        //2. beta
        beta(j) = R::rnorm(tsigma2_j * u_j, sqrt(tsigma2_j));
        
        //3. psi
        if(beta(j) == 0)
        {
          psi(j) = 1;
        } else {
          nu_j = sqrt2 * exp(log_sigma2(j)/2) / abs(beta(j));
          psi(j) = 1 / rinvGauss(nu_j, 2);
        }
      }
    }
    h2 = 1 - exp(log_res);
    
    keep_log_sigma2.row(n_samples + i) = log_sigma2;
    keep_beta.row(n_samples + i) = beta;
    keep_psi.row(n_samples + i) = psi;
    keep_h2(n_samples + i) = h2;
    keep_h2_beta(n_samples + i) = 0;
    for(j = 0; j < M; ++j)
    {
      for(SparseMatrix<double>::InnerIterator it(R,j); it; ++it)
      {
        k = it.row();
        if(k > j)
        {
          keep_h2_beta(n_samples + i) += 2 * it.value() * beta(k) * beta(j);
        }
      }
      keep_h2_beta(n_samples + i) += pow(beta(j), 2);
    }
    
    ++i;
    
    // 5. Switch 2 variables after every 10 steps.
    if(((i + 1) % 10 == 0) && (i < mcmc_n))
    {
      //Choose the index for switching.
      order_sigma2 = order(log_sigma2, _["decreasing"] = true);
      order_sigma2 = order_sigma2 - 1;
      k = 10 + R::rgeom(0.3);
      if(k > M) k = M - 1;
      
      for(j=0; j<k; ++j)
      {
        j1 = order_sigma2(j); //index of SNP 1
        if(LD_pairs.col(j1).nonZeros() >= 1)
        {
          r = R::runif(0,1);
          u = 0;
          for(SparseMatrix<double>::InnerIterator it(LD_pairs,j1); it; ++it)
          {
            u += it.value();
            if(u > r)
            {
              j2 = it.row();
              p_12 = it.value();
              break;
            }
          }
          p_21 = LD_pairs.coeff(j1,j2);
          
          rss = 0;
          W_12 = W.coeff(j1,j2);
          if(W_12 > 0)
          {
            for(SparseMatrix<double>::InnerIterator it(W,j1); it; ++it)
            {
              rss += it.value() * beta(it.row());
            }
            for(SparseMatrix<double>::InnerIterator it(W,j2); it; ++it)
            {
              rss -= it.value() * beta(it.row());
            }
            rss += W_12 * (beta(j1) - beta(j2));
            rss *= (beta(j1) - beta(j2));
            
            rss += (NbetaHat(j1) - NbetaHat(j2)) * (beta(j2) - beta(j1));
          } else {
            for(SparseMatrix<double>::InnerIterator it(W,j1); it; ++it)
            {
              rss += it.value() * beta(it.row());
            }
            for(SparseMatrix<double>::InnerIterator it(W,j2); it; ++it)
            {
              rss += it.value() * beta(it.row());
            }
            rss -= W_12 * (beta(j1) + beta(j2));
            rss *= (beta(j1) + beta(j2));
            
            rss += -(NbetaHat(j1) + NbetaHat(j2)) * (beta(j1) + beta(j2));
          }
          
          acc_prob = exp(rss) *
            exp((a(j1) - a(j2)) * (log_sigma2(j2) - log_sigma2(j1))) *
            p_21 / p_12;
          if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
          {
            if(W_12 > 0)
            {
              ch = beta(j1);
              beta(j1) = beta(j2);
              beta(j2) = ch;
            } else {
              ch = beta(j1);
              beta(j1) = -beta(j2);
              beta(j2) = -ch;
            }
            
            ch = log_sigma2(j1);
            log_sigma2(j1) = log_sigma2(j2);
            log_sigma2(j2) = ch;
            
            ch = psi(j1);
            psi(j1) = psi(j2);
            psi(j2) = ch;
          }
        }
      }
      
      keep_log_sigma2.row(n_samples + i) = log_sigma2;
      keep_beta.row(n_samples + i) = beta;
      keep_psi.row(n_samples + i) = psi;
      keep_h2(n_samples + i) = h2;
      keep_h2_beta(n_samples + i) = 0;
      for(j = 0; j < M; ++j)
      {
        for(SparseMatrix<double>::InnerIterator it(R,j); it; ++it)
        {
          k = it.row();
          if(k > j)
          {
            keep_h2_beta(n_samples + i) += 2 * it.value() * beta(k) * beta(j);
          }
        }
        keep_h2_beta(n_samples + i) += pow(beta(j), 2);
      }
      
      ++i;
    }
  }
  
  mcmc_samples = List::create(Named("n_samples") = keep_beta.rows(),
                              Named("n_burnin") = n_burnin,
                              Named("log_sigma2") = keep_log_sigma2,
                              Named("beta") = keep_beta,
                              Named("psi") = keep_psi,
                              Named("h2") = keep_h2,
                              Named("h2_beta") = keep_h2_beta);
  h2D2.slot("mcmc_samples") = mcmc_samples;
}