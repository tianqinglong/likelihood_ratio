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

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
generate_censored_data(5, 0.2, 2, 1)
*/
