#include <Rcpp.h>
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
List generate_censored_data(int r, double pf, double beta, double eta)
{
  int n;
  n = floor(r/pf);
  
  NumericVector y_obs (n);
  
  double x = 0, y, t_c;
  
  t_c = eta*pow(-log(1-pf), 1/beta);
  NumericVector u_array = runif(n, 0, 1);
  
  int i;
  for(i=0; i<n; i++)
  {
    y = 1 - (1-x)*pow(1-u_array[i], 1.0/(n-i));

    if(y > pf)
      break;
    
    y_obs[i] = y;
    x = y;
  }
  
  if(i < 2)
  {
    return generate_censored_data(r, pf, beta, eta);
  }
  
  NumericVector y_out (i);
  for(int j=0; j<i; j++)
  {
    y_out[j] = eta*pow(-log(1-y_obs[j]),1/beta);
  }
  
  List L = List::create(_["Number_of_Failures"] = i, _["Censor_Time"] = t_c, _["Failure_Times"] = y_out, _["Total_Number"] = n);
  
  return L;
}

// [[Rcpp::export]]
List bootstrap_sample(NumericVector mles, double t_c, int n)
{
  double beta = mles[0], eta = mles[1];
  NumericVector boot_sample = rweibull(n, beta, eta), censored = {};
  int i, r = 0;
  for (i=0; i<boot_sample.length(); i++)
  {
    if (boot_sample[i] <= t_c)
    {
      r++;
      censored.push_back(boot_sample[i]);
    }
  }

  if (r<2)
  {
    return bootstrap_sample(mles, t_c, n);
  }

  List L = List::create(_["Number_of_Failures"] = r, _["Censor_Time"] = t_c, _["Failure_Times"] = censored, _["Total_Number"] = n);
  return L;
}

// [[Rcpp::export]]
NumericVector generate_y_n(int r, int n, double t_c, double t_w, double beta_used, double eta_used, int number_of_y_n)
{
  double p = (R::pweibull(t_w, beta_used, eta_used, true, false)-R::pweibull(t_c, beta_used, eta_used, true, false))/R::pweibull(t_c, beta_used, eta_used, false, false);
  if (std::isnan(p))
  {
    p = 1;
  }
  NumericVector y_n = Rcpp::rbinom(number_of_y_n, n-r, p);

  return y_n;
}

// [[Rcpp::export]]
double compute_p(double t_c, double t_w, double beta, double eta)
{
  double p;
  p = (R::pweibull(t_w, beta, eta, true, false)-R::pweibull(t_c, beta, eta, true, false))/R::pweibull(t_c, beta, eta, false, false);
  if (std::isnan(p))
  {
    p = 1;
  }
  return p;
}

// [[Rcpp::export]]
double pred_dist(int y, NumericVector p, int n)
{
  int B = p.length();
  
  double pred = 0;
  for(int i=0; i<B; i++)
  {
    pred += R::pbinom(y, n, p[i], true, false);
  }
  
  pred = pred/B;
  return pred;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
