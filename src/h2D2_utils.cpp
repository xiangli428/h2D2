// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

#include <omp.h>
//[[Rcpp::plugins(openmp)]]

using namespace std;
using namespace Rcpp;
using namespace Eigen;

//' The inverse Gaussian distribution
//' 
//' @description Random generation for the inverse Gaussian distribution with 
//' parameters \\code{mu} and \\code{lambda}.
//'
//' @usage
//' rinvGauss(mu, lambda)
//' 
//' @param mu The mean parameter. Must be positive.
//' @param lambda The shape parameter. Must be positive.
//' 
//' @details
//' The inverse Gaussian distribution with parameters $mean = \\mu$ and 
//' $shape = \\lambda$ has density
//' $$
//' f(x) = \\sqrt\{\\frac\{\\lambda\}\{2 \\pi x^3\}\} 
//' \\exp \{ -\\frac\{\\lambda (x-\\mu)^2\}\{2 \\mu^2 x\} \}
//' $$
//' 
//' @return A random number from inverse Gaussian distribution.
// [[Rcpp::export]]
double rinvGauss(double mu, double lambda)
{
  double y = R::rchisq(1);
  double x;
  double u = 4 * lambda / (mu * y);
  if(u > 1e-11)
  {
    x = mu * (1 - 2 / (1 + sqrt(1 + u)));
  } else {
    x = lambda / y;
  }
  double z = R::runif(0,1);
  if(z <= mu / (mu + x))
    return(x);
  else
    return(pow(mu, 2)/x);
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
List h2D2_pretrain(S4 h2D2,
                   List sample,
                   const Eigen::VectorXd dW,
                   const Eigen::SparseMatrix<double> W,
                   const Eigen::VectorXd mu,
                   const Eigen::SparseMatrix<double> LD_pairs,
                   const unsigned int n_chain,
                   const NumericVector temp,
                   const unsigned int pre_mcmc_n = 400,
                   const unsigned int pre_burn_in = 200,
                   const double pre_p = 0.05,
                   const unsigned int pre_maxiter = 10,
                   const unsigned int pre_miniter = 2,
                   const double stepsize = 2,
                   const unsigned int seed = 428)
{
  Function set_seed("set.seed");
  set_seed(seed);

  // Data
  const int M = h2D2.slot("M");

  Environment methods = Environment::namespace_env("methods");
  Function AS = methods["as"];
  const SparseMatrix<double> R = as<SparseMatrix<double>>(
    AS(h2D2.slot("R"), _["Class"] = "generalMatrix"));

  // Hyper parameters
  const NumericVector a = h2D2.slot("a");
  const double sqrt2 = sqrt(2);
  const double asum = sum(a);
  double b = asum * (1-1e-4) * 1e4; // initial b

  // Initial samples
  int i = 0; // index for mcmc iterations
  int ii; // index for thinning
  int j, j1, j2 = 0, k; // index for variant
  int n; // index for chain

  MatrixXd beta(M, n_chain), tau(M, n_chain), psi(M, n_chain);
  VectorXd log_res = VectorXd::Constant(n_chain,1e-4);
  VectorXd ch_v(M);
  double h2_beta;

  beta = sample[0];
  tau = sample[1];
  psi = sample[2];

  NumericVector keep_h2_beta(pre_mcmc_n - pre_burn_in);

  // mcmc
  double sigma2_j, tau_j_new, sigma2_j_new, tsigma2_j, tsigma2_j_new;
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
  double pval, estimate;

  double h2_est = 1e-4;
  bool dir = true;

  while(i < pre_maxiter)
  {
    ii = 0;
    while(ii < pre_mcmc_n)
    {
      for(k = 0; k < M; ++k)
      {
        j = order_SNP(k);

        for(n = 0; n < n_chain; ++n)
        {
          //1. sigma2
          sigma2_j = exp(tau(j,n));
          log_res(n) = lpexp(log_res(n), tau(j,n));

          tau_j_new = propose_lognormal(tau(j,n) - log_res(n), stepsize) +
            log_res(n);
          acc_prob = (tau_j_new - log_res(n)) / (tau(j,n) - log_res(n));

          sigma2_j_new = exp(tau_j_new);

          tsigma2_j = 1 / (dW(j) / temp(n) + 1 / (sigma2_j * psi(j,n)));
          tsigma2_j_new = 1 / (dW(j) / temp(n) + 1 / (sigma2_j_new * psi(j,n)));

          u_j = mu(j);
          for(SparseMatrix<double>::InnerIterator it(W,j); it; ++it)
          {
            u_j -= it.value() * beta(it.row(),n);
          }
          u_j /= temp(n);

          acc_prob *= sqrt((dW(j) * sigma2_j * psi(j,n) + temp(n)) /
            (dW(j) * sigma2_j_new * psi(j,n) + temp(n))) *
              exp(pow(u_j, 2) * (tsigma2_j_new - tsigma2_j) / 2) *
              exp(a(j) * (tau_j_new - tau(j,n))) *
              exp((b-1) * (lmexp(log_res(n), tau_j_new) -
              lmexp(log_res(n), tau(j,n))));
          if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
          {
            tau(j,n) = tau_j_new;
            sigma2_j = sigma2_j_new;
            tsigma2_j = tsigma2_j_new;
          }
          log_res(n) = lmexp(log_res(n), tau(j,n));

          //2. beta
          beta(j,n) = R::rnorm(tsigma2_j * u_j, sqrt(tsigma2_j));

          //3. psi
          if(beta(j,n) == 0)
          {
            psi(j,n) = 1;
          } else {
            nu_j = sqrt2 * exp(tau(j,n)/2) / abs(beta(j,n));
            psi(j,n) = 1 / rinvGauss(nu_j, 2);
          }
        }
      }

      h2_beta = 0;
      for(j = 0; j < M; ++j)
      {
        for(SparseMatrix<double>::InnerIterator it(R,j); it; ++it)
        {
          k = it.row();
          if(k > j)
          {
            h2_beta += 2 * it.value() * beta(k,0) * beta(j,0);
          }
        }
        h2_beta += pow(beta(j,0), 2);
      }

      if(ii >= pre_burn_in)
      {
        keep_h2_beta(ii - pre_burn_in) = h2_beta;
      }

      ++ii;

      if(((ii + 1) % 5 == 0) && (ii < pre_mcmc_n))
      {
        if(((ii + 1) % 20 != 0))
        {
          for(n = 0; n < n_chain; ++n)
          {
            //Choose the index for switching.
            order_sigma2 = order(tau.col(n), _["decreasing"] = true);
            order_sigma2 = order_sigma2 - 1;
            k = 10 + R::rgeom(0.3);
            if(k > M) k = M - 1;

            for(j = 0; j < k; ++j)
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
                    rss += it.value() * beta(it.row(),n);
                  }
                  for(SparseMatrix<double>::InnerIterator it(W,j2); it; ++it)
                  {
                    rss -= it.value() * beta(it.row(),n);
                  }
                  rss += W_12 * (beta(j1,n) - beta(j2,n));
                  rss *= (beta(j1,n) - beta(j2,n));

                  rss += (mu(j1) - mu(j2)) * (beta(j2,n) - beta(j1,n)) +
                    (dW(j1) - dW(j2)) * (pow(beta(j1,n),2) - pow(beta(j2,n),2)) / 2;
                } else {
                  for(SparseMatrix<double>::InnerIterator it(W,j1); it; ++it)
                  {
                    rss += it.value() * beta(it.row(),n);
                  }
                  for(SparseMatrix<double>::InnerIterator it(W,j2); it; ++it)
                  {
                    rss += it.value() * beta(it.row(),n);
                  }
                  rss -= W_12 * (beta(j1,n) + beta(j2,n));
                  rss *= (beta(j1,n) + beta(j2,n));

                  rss += -(mu(j1) + mu(j2)) * (beta(j1,n) + beta(j2,n)) +
                    (dW(j1) - dW(j2)) * (pow(beta(j1,n),2) - pow(beta(j2,n),2)) / 2;
                }

                acc_prob = exp(rss / temp(n)) *
                  exp((a(j1) - a(j2)) * (tau(j2) - tau(j1))) *
                  p_21 / p_12;
                if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
                {
                  if(W_12 > 0)
                  {
                    ch = beta(j1,n);
                    beta(j1,n) = beta(j2,n);
                    beta(j2,n) = ch;
                  } else {
                    ch = beta(j1,n);
                    beta(j1,n) = -beta(j2,n);
                    beta(j2,n) = -ch;
                  }

                  ch = tau(j1,n);
                  tau(j1,n) = tau(j2,n);
                  tau(j2,n) = ch;

                  ch = psi(j1,n);
                  psi(j1,n) = psi(j2,n);
                  psi(j2,n) = ch;
                }
              }
            }
          }
        } else {
          for(n = n_chain-2; n >= 0; --n)
          {
            rss = mu.dot(beta.col(n)) - beta.col(n).dot(W * beta.col(n)) / 2;
            for(j = 0; j < M; ++j)
            {
              rss -= dW(j) * pow(beta(j,n),2) / 2;
            }
            acc_prob = rss * (1 / temp(n+1) - 1 / temp(n));

            rss = mu.dot(beta.col(n+1)) - beta.col(n+1).dot(W * beta.col(n+1)) / 2;
            for(j = 0; j < M; ++j)
            {
              rss -= dW(j) * pow(beta(j,n+1),2) / 2;
            }
            acc_prob -= rss * (1 / temp(n+1) - 1 / temp(n));
            acc_prob = exp(acc_prob);

            if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
            {
              ch_v = beta.col(n);
              beta.col(n) = beta.col(n+1);
              beta.col(n+1) = ch_v;

              ch_v = tau.col(n);
              tau.col(n) = tau.col(n+1);
              tau.col(n+1) = ch_v;

              ch_v = psi.col(n);
              psi.col(n) = psi.col(n+1);
              psi.col(n+1) = ch_v;

              ch = log_res(n);
              log_res(n) = log_res(n+1);
              log_res(n+1) = ch;
            }
          }
        }

        h2_beta = 0;
        for(j = 0; j < M; ++j)
        {
          for(SparseMatrix<double>::InnerIterator it(R,j); it; ++it)
          {
            k = it.row();
            if(k > j)
            {
              h2_beta += 2 * it.value() * beta(k,0) * beta(j,0);
            }
          }
          h2_beta += pow(beta(j,0), 2);
        }

        if(ii >= pre_burn_in)
        {
          keep_h2_beta(ii - pre_burn_in) = h2_beta;
        }

        ++ii;
      }
    }

    estimate = mean(keep_h2_beta);
    if(estimate >= 1e-6 && estimate < 1)
    {
      test = t_test(_["x"] = keep_h2_beta, _["mu"] = h2_est);
      pval = test["p.value"];
      estimate = test["estimate"];
      if(i == pre_miniter) dir = estimate > h2_est;

      cout << "Pretrain step " << i << ", b=" << b << ", h2_est=" <<
        h2_est << ", mean=" << estimate << ", p=" << pval <<
          "." << endl;

      if(((pval >= pre_p) || (estimate > h2_est) != dir) && i >= pre_miniter)
      {
        break;
      } else {
        b = asum * (1 - estimate) / estimate;
        h2_est = estimate;
      }
    } else {
      estimate = 1e-6;
      b = asum * (1 - estimate) / estimate;
      h2_est = estimate;
    }

    ++i;
  }

  cout << "End pretrain. b=" << b << "." << endl;

  h2D2.slot("b") = b;
  
  List last_sample = List::create(beta, tau, psi);
  return(last_sample);
}

// [[Rcpp::export]]
List h2D2_sampling(S4 h2D2, 
                   List sample,
                   const Eigen::VectorXd dW,
                   const Eigen::SparseMatrix<double> W,
                   const Eigen::VectorXd mu,
                   const Eigen::SparseMatrix<double> LD_pairs,
                   const unsigned int n_chain,
                   const NumericVector temp,
                   const unsigned int mcmc_n = 100, 
                   const unsigned int thin = 1, 
                   const double stepsize = 2, 
                   const unsigned int seed = 428)
{
  Function set_seed("set.seed");
  set_seed(seed);
  
  // Data
  const int M = h2D2.slot("M");
  
  Environment methods = Environment::namespace_env("methods");
  Function AS = methods["as"];
  const SparseMatrix<double> R = as<SparseMatrix<double>>(
    AS(h2D2.slot("R"), _["Class"] = "generalMatrix"));
  
  // Hyper parameters
  const NumericVector a = h2D2.slot("a");
  const double b = h2D2.slot("b");
  const double sqrt2 = sqrt(2);
  
  // Initial samples
  int i = 0; // index for mcmc iterations
  int ii; // index for thinning
  int j, j1, j2 = 0, k; // index for variant
  int n; // index for chain
  
  MatrixXd beta(M, n_chain), tau(M, n_chain), psi(M, n_chain);
  VectorXd log_res(n_chain);
  VectorXd ch_v(M);
  
  beta = sample[0];
  tau = sample[1];
  psi = sample[2];
  
  for(n = 0; n < n_chain; ++n)
  {
    log_res(n) = 0;
    for(j = 0; j < M; ++j)
    {
      log_res(n) = lmexp(log_res(n), tau(j, n));
    }
  }
  
  MatrixXd keep_beta = MatrixXd::Zero(M,mcmc_n);
  MatrixXd keep_tau = MatrixXd::Zero(M,mcmc_n);
  MatrixXd keep_psi = MatrixXd::Zero(M,mcmc_n);
  VectorXd keep_h2_beta = VectorXd::Zero(mcmc_n);
  VectorXd keep_h2 = VectorXd::Zero(mcmc_n);
  
  // mcmc
  double sigma2_j, tau_j_new, sigma2_j_new, tsigma2_j, tsigma2_j_new;
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
        
        for(n = 0; n < n_chain; ++n)
        {
          //1. sigma2
          sigma2_j = exp(tau(j,n));
          log_res(n) = lpexp(log_res(n), tau(j,n));
          
          tau_j_new = propose_lognormal(tau(j,n) - log_res(n), stepsize) + 
            log_res(n);
          acc_prob = (tau_j_new - log_res(n)) / (tau(j,n) - log_res(n));
          
          sigma2_j_new = exp(tau_j_new);
          
          tsigma2_j = 1 / (dW(j) / temp(n) + 1 / (sigma2_j * psi(j,n)));
          tsigma2_j_new = 1 / (dW(j) / temp(n) + 1 / (sigma2_j_new * psi(j,n)));
          
          u_j = mu(j);
          for(SparseMatrix<double>::InnerIterator it(W,j); it; ++it)
          {
            u_j -= it.value() * beta(it.row(),n);
          }
          u_j /= temp(n);
          
          acc_prob *= sqrt((dW(j) * sigma2_j * psi(j,n) + temp(n)) /
            (dW(j) * sigma2_j_new * psi(j,n) + temp(n))) *
              exp(pow(u_j, 2) * (tsigma2_j_new - tsigma2_j) / 2) *
              exp(a(j) * (tau_j_new - tau(j,n))) *
              exp((b-1) * (lmexp(log_res(n), tau_j_new) - 
              lmexp(log_res(n), tau(j,n))));
          if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
          {
            tau(j,n) = tau_j_new;
            sigma2_j = sigma2_j_new;
            tsigma2_j = tsigma2_j_new;
          }
          log_res(n) = lmexp(log_res(n), tau(j,n));
          
          //2. beta
          beta(j,n) = R::rnorm(tsigma2_j * u_j, sqrt(tsigma2_j));
          
          //3. psi
          if(beta(j,n) == 0)
          {
            psi(j,n) = 1;
          } else {
            nu_j = sqrt2 * exp(tau(j,n)/2) / abs(beta(j,n));
            psi(j,n) = 1 / rinvGauss(nu_j, 2);
          }
        }
      }
    }
    
    keep_beta.col(i) = beta.col(0);
    keep_tau.col(i) = tau.col(0);
    keep_psi.col(i) = psi.col(0);
    keep_h2(i) = 1 - exp(log_res(0));
    
    keep_h2_beta(i) = 0;
    for(j = 0; j < M; ++j)
    {
      for(SparseMatrix<double>::InnerIterator it(R,j); it; ++it)
      {
        k = it.row();
        if(k > j)
        {
          keep_h2_beta(i) += 2 * it.value() * beta(k,0) * beta(j,0);
        }
      }
      keep_h2_beta(i) += pow(beta(j,0), 2);
    }
    
    ++i;
    
    if(((i + 1) % 10 == 0) && (i < mcmc_n))
    {
      if(((i + 1) % 50 != 0))
      {
        for(n = 0; n < n_chain; ++n)
        {
          //Choose the index for switching.
          order_sigma2 = order(tau.col(n), _["decreasing"] = true);
          order_sigma2 = order_sigma2 - 1;
          k = 10 + R::rgeom(0.3);
          if(k > M) k = M - 1;
          
          for(j = 0; j < k; ++j)
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
                  rss += it.value() * beta(it.row(),n);
                }
                for(SparseMatrix<double>::InnerIterator it(W,j2); it; ++it)
                {
                  rss -= it.value() * beta(it.row(),n);
                }
                rss += W_12 * (beta(j1,n) - beta(j2,n));
                rss *= (beta(j1,n) - beta(j2,n));
                
                rss += (mu(j1) - mu(j2)) * (beta(j2,n) - beta(j1,n)) +
                  (dW(j1) - dW(j2)) * (pow(beta(j1,n),2) - pow(beta(j2,n),2)) / 2;
              } else {
                for(SparseMatrix<double>::InnerIterator it(W,j1); it; ++it)
                {
                  rss += it.value() * beta(it.row(),n);
                }
                for(SparseMatrix<double>::InnerIterator it(W,j2); it; ++it)
                {
                  rss += it.value() * beta(it.row(),n);
                }
                rss -= W_12 * (beta(j1,n) + beta(j2,n));
                rss *= (beta(j1,n) + beta(j2,n));
                
                rss += -(mu(j1) + mu(j2)) * (beta(j1,n) + beta(j2,n)) +
                  (dW(j1) - dW(j2)) * (pow(beta(j1,n),2) - pow(beta(j2,n),2)) / 2;
              }
              
              acc_prob = exp(rss / temp(n)) *
                exp((a(j1) - a(j2)) * (tau(j2) - tau(j1))) *
                p_21 / p_12;
              if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
              {
                if(W_12 > 0)
                {
                  ch = beta(j1,n);
                  beta(j1,n) = beta(j2,n);
                  beta(j2,n) = ch;
                } else {
                  ch = beta(j1,n);
                  beta(j1,n) = -beta(j2,n);
                  beta(j2,n) = -ch;
                }
                
                ch = tau(j1,n);
                tau(j1,n) = tau(j2,n);
                tau(j2,n) = ch;
                
                ch = psi(j1,n);
                psi(j1,n) = psi(j2,n);
                psi(j2,n) = ch;
              }
            }
          }
        }
      } else {
        for(n = n_chain-2; n >= 0; --n)
        {
          rss = mu.dot(beta.col(n)) - beta.col(n).dot(W * beta.col(n)) / 2;
          for(j = 0; j < M; ++j)
          {
            rss -= dW(j) * pow(beta(j,n),2) / 2;
          }
          acc_prob = rss * (1 / temp(n+1) - 1 / temp(n));
          
          rss = mu.dot(beta.col(n+1)) - beta.col(n+1).dot(W * beta.col(n+1)) / 2;
          for(j = 0; j < M; ++j)
          {
            rss -= dW(j) * pow(beta(j,n+1),2) / 2;
          }
          acc_prob -= rss * (1 / temp(n+1) - 1 / temp(n));
          acc_prob = exp(acc_prob);
          
          if((acc_prob > 1) || (R::runif(0,1) < acc_prob))
          {
            ch_v = beta.col(n);
            beta.col(n) = beta.col(n+1);
            beta.col(n+1) = ch_v;
            
            ch_v = tau.col(n);
            tau.col(n) = tau.col(n+1);
            tau.col(n+1) = ch_v;
            
            ch_v = psi.col(n);
            psi.col(n) = psi.col(n+1);
            psi.col(n+1) = ch_v;
            
            ch = log_res(n);
            log_res(n) = log_res(n+1);
            log_res(n+1) = ch;
          }
        }
      }
      
      keep_beta.col(i) = beta.col(0);
      keep_tau.col(i) = tau.col(0);
      keep_psi.col(i) = psi.col(0);
      keep_h2(i) = 1 - exp(log_res(0));
      
      keep_h2_beta(i) = 0;
      for(j = 0; j < M; ++j)
      {
        for(SparseMatrix<double>::InnerIterator it(R,j); it; ++it)
        {
          k = it.row();
          if(k > j)
          {
            keep_h2_beta(i) += 2 * it.value() * beta(k,0) * beta(j,0);
          }
        }
        keep_h2_beta(i) += pow(beta(j,0), 2);
      }
      
      ++i;
    }
  }
  
  List last_sample = List::create(beta, tau, psi);
  
  return(List::create(keep_beta, keep_tau, keep_psi, keep_h2_beta, keep_h2,
                      last_sample));
}