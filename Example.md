Example file running simulation code:


File clambdabetapython.cpp  runs the Lambda-variant of  Beta(2-a,a)-coalescent, i.e.\ the one derived by Schweinsberg (2003). 
It also runs  the Kingman-coalescent, and  a time-changed variant corresponding to exponential population growth. 
Requirements:  A C++ compiler (e.g. g++ version 8.3.0), and  Python version  3.7.3, and GSL (GNU Scientific Library) version 2.5.
Note that the code is called from a python script, hence one compiles with
g++ -Wall -fPIC -shared -o clambdabetapython.so clambdabetapython.cpp
