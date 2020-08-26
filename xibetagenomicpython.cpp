#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <list>
#include <forward_list>
#include <assert.h>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <gsl/gsl_matrix.h>
#include "xibetagenomiconstants.hpp"
#include "xibetagenomicpython.hpp"

const unsigned samplesize = 180;
const unsigned lgroups = 20;
const unsigned Imin = 1;
const unsigned Imax = 1;
const unsigned Jmin = 15;
const unsigned Jmax = samplesize-1;
const unsigned Kmin = Imin;
const unsigned Kmax = Jmax;
const unsigned trials = 1;

/*
in bash:
g++ -Wall -fPIC -shared -o <thisfile>.so <thisfile>.cpp -lm -lgsl -lgslcblas
*/
