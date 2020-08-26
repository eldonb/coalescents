#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <utility>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <list>
#include <string>
#include <fstream>
#include <forward_list>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/uniform_01.hpp>
#include "clambdabetapython.hpp"

  

using namespace std;
using namespace boost;
using namespace boost::random;
using namespace math;

/* compiling:
g++ -Wall -c -fPIC clambdabetapython.cpp -o clambdabetapython.o 
g++ -shared -o clambdabetapython.so clambdabetapython.so
------------------
in bash:
g++ -Wall -fPIC -shared -o clambdabetapython.so clambdabetapython.cpp
in python:
import clambdabetalib
import numpy as np
from scipy import stats
m=np.loadtxt('egresskra0', usecols=range(0,3))
k=stats.gaussian_kde(m.T, bw_method=0.01)
given datapoint
gagnapunktur=[.16, 0.5, 236]
print( k(gagnapunktur))
[0.00118435]
*/
