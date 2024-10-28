/* The random number generator used by all samplers. */

/*
Set up random number generation such that:
  (a) Results are reprodicible;
  (b) In R, random number generation carries on from where it left over previously;
  (c) By default, the same stream is used for all random number generation;
  (d) Multiple independent streams can be used for parallel random number generation.
  
This is acheived by defining the RNG in its own namespace, which is compiled together with all of the random number generating functions.
*/

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "pcg_random.hpp"

using uniform_distribution = boost::random::uniform_real_distribution<double>;
using normal_distribution = boost::random::normal_distribution<double>;

namespace RNG
{
    pcg32 rng(12345);
    uniform_distribution runif(0., 1.); // Boost function
    normal_distribution rnorm(0., 1.);  // Ziggurat method, much faster than R's inversion.

    Eigen::VectorXd GaussianVector(const int &d)
    {
        Eigen::VectorXd result(d);
        for (int i = 0; i < d; i++)
        {
            result[i] = RNG::rnorm(RNG::rng);
        }
        return result;
    }

} // namespace RNG


// Re-Seed RNG from R
// 
//' @export
// [[Rcpp::export(rng = false)]]
void SetSeed_cpp(const int &seed, const int &stream = 0)
{
    pcg32 rng_refresh(seed, stream);

    RNG::rng = rng_refresh;
}
