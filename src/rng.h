/* This header links to the random number generator, which is intialized exactly once for all samplers when the package is loaded. */

#ifndef RNG_H
#define RNG_H

// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "pcg_random.hpp"

using uniform_distribution = boost::random::uniform_real_distribution<double>;
using normal_distribution = boost::random::normal_distribution<double>;

namespace RNG
{
    // Declare rng
    extern pcg32 rng;

    // Declare distributions
    extern uniform_distribution runif;
    extern normal_distribution rnorm;

} // namesepace RNG

void SetSeed_cpp (int seed, int stream);

#endif /* RNG_H */
