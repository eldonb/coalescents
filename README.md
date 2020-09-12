# coalescents
Scripts for simulating coalescents and analysing the output. 

The scripts are part of a manuscript currently undergoing peer-review.  Citation information will appear in due course. 


# Example

The  example file {codexample.sfs}  contains  examples for compiling and  running the code and   generating  simulated datapoints and comparing them to  data.

Depending on the number of  points in the grid of parameter values,  the simulations can take a while, easily several days.  

# 1. System requirements
The main requirements are  a  C++ compiler (e.g. g++ version 8.3.0),  Python version 3.7.3,   GSL (GNU Scientific Library) version 2.5, and the Boost C++ library. 

The code compiles successfully on a Linux Debian 4.19.0-10-amd64 x86-64 architecture. 
No non-standard hardware is required.

# 2. Installation guide
The file {codexample.sfs} contains instructions for compiling the code.  Compiling takes only few seconds.

# 3. Demo
The file {codexample.sfs} contains instructions for compiling and  running the code to generate  simulated datapoints, and  
for further processing the simulated datapoints, e.g. for comparison with actual data.

The expected output is described in the file {codexample.sfs}.

Expected run time varies depending on the specific setup, e.g.  the number of  grid points in the grid of parameter values, and the 
number of available CPUs.  Running time can, therefore, take from a few minutes to several days. 

# 4. Instructions for use
The file {codexample.sfs} contains instructions for compiling and running the scripts, and for reproducing the output. 
However,  the running time for the  results in the manuscript was many days on a powerful computer with  64 CPUs.
