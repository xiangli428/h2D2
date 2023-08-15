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

// [[Rcpp::export]]
double propose_lognormal(double x, double stepsize = 2)
{
  x *= exp(R::rnorm(0, stepsize));
  return(x);
}

// [[Rcpp::export]]
void h2D2_sampling(S4 h2D2, 
                   List L,
                   int mcmc_n = 100, 
                   int thin = 1, 
                   double stepsize = 2, 
                   double scale2 = 1, 
                   unsigned int seed = 428)
{
  Function set_seed("set.seed");
  set_seed(seed);
  
  // Data
  const int M = h2D2.slot("M");
  
  // const VectorXd betaHat = h2D2.slot("betaHat");
  // const VectorXd S = h2D2.slot("S");
  // VectorXd S_2 = VectorXd::Zero(M);
  // VectorXd S_2betaHat = VectorXd::Zero(M);
  // int j;
  // for(j = 0; j < M; ++j)
  // {
  //   S_2(j) = 1 / pow(S(j), 2);
  //   S_2betaHat(j) = S_2(j) * betaHat(j);
  // }
  
  // const VectorXd S_2 = h2D2.slot("S_2");
  // const VectorXd S_2betaHat = h2D2.slot("S_2betaHat");
  
  const VectorXd S_2 = L[1];
  const VectorXd S_2betaHat = L[2];
  
  Environment methods = Environment::namespace_env("methods");
  Function AS = methods["as"];
  const SparseMatrix<double> R = as<SparseMatrix<double>>(
    AS(h2D2.slot("R"), _["Class"] = "dgCMatrix"));
  const SparseMatrix<double> W = L[0];
  // const SparseMatrix<double> W = as<SparseMatrix<double>>(
  //   AS(h2D2.slot("W"), _["Class"] = "dgCMatrix"));
  const SparseMatrix<double> LD_pairs = h2D2.slot("LD_pairs");
  
  // Hyper parameters
  const VectorXd a = h2D2.slot("a");
  const double b = h2D2.slot("b");
  const double scale = sqrt(scale2);
  
  // Initial samples
  int i = 0;
  int ii;
  List mcmc_samples = h2D2.slot("mcmc_samples");
  
  int n_samples, n_burnin;
  double h2, res, log_res;
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
  res = 1 - h2;
  double sigma2_j, log_sigma2_j_new, sigma2_j_new, tsigma2_j, tsigma2_j_new;
  double u_j, nu_j, acc_prob, rss, ch, r, u;
  double W_12, u_1, u_2, sigma2_1, sigma2_2, p_21;
  double p_12 = 0;
  double log_sigma2_1_new, log_sigma2_2_new, sigma2_1_new, sigma2_2_new;
  double beta_1_new, beta_2_new, psi_1_new, psi_2_new;
  NumericVector order_sigma2;
  
  NumericVector z = h2D2.slot("z");
  Function order("order");
  NumericVector order_SNP = order(abs(z), _["decreasing"] = true);
  order_SNP = order_SNP - 1;
  NumericVector rank_SNP (M);
  for(k = 0; k < M; ++k)
  {
    rank_SNP[order_SNP[k]] = k;
  }
  
  while(i < mcmc_n)
  {
    for(ii = 0; ii < thin; ++ii)
    {
      for(k = 0; k < M; ++k)
      {
        j = order_SNP(k);
        
        //1. sigma2
        sigma2_j = exp(log_sigma2(j));
        res += sigma2_j;
        log_res = log(res);
        
        log_sigma2_j_new = propose_lognormal(log_sigma2(j) - log_res, stepsize) + log_res;
        acc_prob = (log_sigma2_j_new - log_res) / (log_sigma2(j) - log_res);
        
        sigma2_j_new = exp(log_sigma2_j_new);
        
        tsigma2_j = 1 / (S_2(j) + scale2 / (sigma2_j * psi(j)));
        tsigma2_j_new = 1 / (S_2(j) + scale2 / (sigma2_j_new * psi(j)));
        
        u_j = S_2betaHat(j);
        for(SparseMatrix<double>::InnerIterator it(W,j); it; ++it)
        {
          u_j -= it.value() * beta(it.row());
        }
        
        acc_prob *= sqrt((S_2(j) * sigma2_j * psi(j) + scale2) /
          (S_2(j) * sigma2_j_new * psi(j) + scale2)) *
            exp(pow(u_j, 2) * (tsigma2_j_new - tsigma2_j) / 2) *
            exp(a(j) * (log_sigma2_j_new - log_sigma2(j))) *
            pow((res - sigma2_j_new) / (res - sigma2_j), b - 1);
        if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
        {
          log_sigma2(j) = log_sigma2_j_new;
          sigma2_j = sigma2_j_new;
          tsigma2_j = tsigma2_j_new;
        }
        res -= sigma2_j;
        
        //2. beta
        beta(j) = R::rnorm(tsigma2_j * u_j, sqrt(tsigma2_j));
        
        //3. psi
        if(beta(j) == 0)
        {
          psi(j) = 1;
        } else {
          nu_j = exp(log_sigma2(j)/2) / (abs(beta(j)) * scale);
          psi(j) = 1 / rinvGauss(nu_j, 1);
        }
      }
    }
    h2 = 1 - res;
    
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
    
    // 5. Switch 2 variables after every ? steps.
    if(((i + 1) % 10 == 0) && (i < mcmc_n))
    {
      if((i + 1) % 1000 == 0)
      {
        for(k = 0; k < M; ++k)
        {
          j1 = order_SNP(k); //index of SNP 1
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
                break;
              }
            }
            
            if(rank_SNP(j1) > rank_SNP(j2))
            {
              W_12 = W.coeff(j1,j2);
              
              // Joint sampling log_sigma2(j1) and log_sigma2(j2)
              sigma2_1 = exp(log_sigma2(j1));
              sigma2_2 = exp(log_sigma2(j2));
              
              res += sigma2_1 + sigma2_2;
              
              log_sigma2_1_new = propose_lognormal(log_sigma2(j1), stepsize);
              log_sigma2_2_new = propose_lognormal(log_sigma2(j2), stepsize);
              
              sigma2_1_new = exp(log_sigma2_1_new);
              sigma2_2_new = exp(log_sigma2_2_new);
              
              if(res - sigma2_1_new - sigma2_2_new > 0)
              {
                u_1 = S_2betaHat(j1) + W_12 * beta(j2);
                for(SparseMatrix<double>::InnerIterator it(W,j1); it; ++it)
                {
                  u_1 -= it.value() * beta(it.row());
                }
                
                u_2 = S_2betaHat(j2) + W_12 * beta(j1);
                for(SparseMatrix<double>::InnerIterator it(W,j2); it; ++it)
                {
                  u_2 -= it.value() * beta(it.row());
                }
                
                psi_1_new = R::rexp(1);
                psi_2_new = R::rexp(1);
                
                beta_1_new = R::rnorm(0, sqrt(psi_1_new * sigma2_1_new / scale2));
                beta_2_new = R::rnorm(0, sqrt(psi_2_new * sigma2_2_new / scale2));
                
                acc_prob = u_1 * (beta_1_new - beta(j1)) + 
                  u_2 * (beta_2_new - beta(j2)) +
                  S_2(j1) * (pow(beta(j1),2) - pow(beta_1_new,2)) / 2 +
                  S_2(j2) * (pow(beta(j2),2) - pow(beta_2_new,2)) / 2 +
                  W_12 * (beta(j1) * beta(j2) - beta_1_new * beta_2_new) +
                  a(j1) * (log_sigma2_1_new - log_sigma2(j1)) +
                  a(j2) * (log_sigma2_2_new - log_sigma2(j2)) +
                  (b-1) * log((res - sigma2_1_new - sigma2_2_new) /
                    (res - sigma2_1 - sigma2_2)) +
                      log((log_sigma2_1_new / log_sigma2(j1)) *
                      (log_sigma2_2_new / log_sigma2(j2)));
                acc_prob = exp(acc_prob);
                
                if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
                {
                  log_sigma2(j1) = log_sigma2_1_new;
                  log_sigma2(j2) = log_sigma2_2_new;
                  
                  sigma2_1 = sigma2_1_new;
                  sigma2_2 = sigma2_2_new;
                  
                  psi(j1) = psi_1_new;
                  psi(j2) = psi_2_new;
                  
                  beta(j1) = beta_1_new;
                  beta(j2) = beta_2_new;
                }
              }
              
              res -= (sigma2_1 + sigma2_2);
            }
          }
        }
      } else {
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
              
              rss += (S_2betaHat(j1) - S_2betaHat(j2)) * (beta(j2) - beta(j1)) +
                (S_2(j1) - S_2(j2)) * (pow(beta(j1),2) - pow(beta(j2),2)) / 2;
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
              
              rss += -(S_2betaHat(j1) + S_2betaHat(j2)) * (beta(j1) + beta(j2)) +
                (S_2(j1) - S_2(j2)) * (pow(beta(j1),2) - pow(beta(j2),2)) / 2;
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
      }
      
      h2 = 1 - res;
      
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