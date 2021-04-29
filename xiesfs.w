\pdfoutput=1
\documentclass[10pt]{cweb}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}%
\usepackage{marvosym}%
\usepackage{bm}
\usepackage{upgreek}
%%\usepackage[super,round]{natbib}
\usepackage{natbib}
\usepackage[all]{xy}%
\usepackage{color}%
\usepackage{rotating}%
\usepackage{a4wide,fullpage}%
\usepackage{setspace}%
\usepackage{enumerate}%
\usepackage{dsfont}%
\usepackage{textcomp}
\usepackage{environ}%
\usepackage{eufrak}%
\usepackage{hyperref}%
\setstretch{1.5}%
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\} }  } }
\newcommand{\g}{\,\boldsymbol{|}\,}%
\newcommand{\EE}[1]{\mathds{E}\left[ #1 \right]}%
\newcommand{\im}{\ensuremath{\imath} }%
\newcommand{\jm}{\ensuremath{\jmath} }%
\newcommand{\be}{\begin{equation}}%
\newcommand{\ee}{\end{equation}}%
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\hp}{\ensuremath{\theta_1}}%
\newcommand{\hs}{\ensuremath{\theta_2}}%
\newcommand{\D}{\ensuremath{\mathbb{D}}}%
\newcommand{\F}{\ensuremath{\mathbb{F}} }%
\newcommand{\G}{\ensuremath{\mathbb{G}} }%
%%
\newcommand{\bi}[1]{\textcolor{blue}{\it #1}}%
%%
\NewEnviron{esplit}[1]{%
\begin{equation}
\label{#1}
\begin{split}
  \BODY
\end{split}\end{equation}
}
%%
\title{Exact computation of the normalised expected site-frequency spectrum associated with  $\Xi$-coalescents}%%
\author{ Bjarki Eldon\\ Leibniz Institute at MfN Berlin\\ \href{mailto:eldon@@mfn-berlin.de}{\texttt{\Letter: my.name(\emph{you know the symbol})mfn-berlin.de}} }%
\date{\today }%
\begin{document}
\maketitle
\begin{abstract}%
Exact computation of normalised (with expected total length of the genealogy)  expected  branch lengths of a non-recombining one-locus genealogy associated with $\Xi$-coalescents \citep{Blath2016}.     This CWEB 
       \citep{knuth1994cweb} document (the {\tt .w} file) can be compiled with {\tt cweave} to
       generate a {\tt .tex} file, and with {\tt ctangle} to generate
       a {\tt .c}  \citep{kernighan1988c}  file.    Please cite \cite{Blath2016}. 



\end{abstract}%

\tableofcontents

@q generate shasum: shasum -a 512224 file.tgz > file.tgz.sha1 @>
@q check shasum: shasum -a 512224 -c file.tgz.sha1 @>


@* {\bf Copyright}. 

Copyright {\copyright} {\the\year}  Bjarki Eldon \newline

This document and any source code it contains  is distributed under the terms of the GNU General Public Licence (version $\ge 3$).  You
should have received a copy of the licence along with this file (see file COPYING).  


    The source codes  described in this document  are  free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This document and the code it contains   is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.


@* {\bf intro}. 
\label{intro}


Define $\one{A} := 1$ if $A$ holds, and zero otherwise $(\one{1 < 0} = 0, \one{1 > 0} = 1)$ .    Let $b\ge 2$ denote the current number of active blocks, and  $k_1, \ldots, k_r \ge 2$ the merger sizes of $r$ mergers, $2\le k_1 + \cdots + k_r \le b$, and $s = b - k_1 - \cdots - k_r$.   In a diploid population $1 \le r \le 4$.    Define, for some $0<\psi \le 1$, 
\begin{esplit}{Cconst}
C_{b; k_1, \ldots, k_r;s } &=  \frac{4}{\psi^2} \sum_{\ell = 0}^{s \wedge (4-r)} \binom{s}{\ell} (4)_{r+\ell} (1-\psi)^{s-\ell }\left( \frac{\psi}{4} \right) ^{k_1 + \cdots + k_r + \ell}  \\
\end{esplit}
   The generator of the model in  \cite{BBE13} is given by
\begin{equation}\label{G}
G = \begin{cases}%
1 + c\frac{\psi^2}{4}C_{b;2;b-2} & \text{pairwise merger} \\
  c\frac{\psi^2}{4} C_{b;k_1, \ldots, k_r; s } & \text{at least three-merger} \\
\end{cases}%
\end{equation}%
 where $ C_{b;k_1, \ldots, k_r; s } $ is from Eq \eqref{Cconst}.   
For given merger sizes $k_1, \ldots, k_r$ let  $\ell_j := \sum_{i=1}^r\one{k_i = j} $.  The number of possible such  $(k_1, \ldots, k_r)$ mergers is
\be\label{N}
    \mathcal{N}(b; k_1, \ldots, k_r ) = \binom{b}{k_1 \ldots k_r\, s} \frac{1}{ \prod_{j=2}^b\ell_j!  }
\ee
 \cite{S00} for $b\ge 2$ and $2 \le k_1 + \cdots + k_r \le b$ and $s = b - k_1 - \cdots - k_r$.   Observe that $ \mathcal{N}(b; 2) = b(b-1)/2$ for any unsigned $b\ge 2$. 

%%The total rate of a $(k_1, \ldots, k_r)$-merger is then
%%\be
%%\lambda_{b;k_1, \ldots, k_r} =    \mathcal{N}(b; k_1, \ldots, k_r ) G(b;k_1, \ldots, k_r)
%%\ee
%%where $ G(b;k_1, \ldots, k_r) $ is from  Eq \eqref{G}.


The  rate in a general  $\Xi$-coalescent is
\be\label{xilambdabk}
    \lambda_{b;k_1, \ldots, k_r} =  \int_\Delta  \sum_{\ell = 0}^s \sum_{i_1\neq \ldots \neq i_{r+\ell}} x_1^{k_1}\cdots x_{i_r }^{k_r } \left(1 - \sum_{j=1}^\infty \right)^{s-\ell}\frac{1}{\sum_{j=1}^\infty x_j^2 }  \Xi(dx)  +  \one{r=1,k_1 = 2}
\ee
where $\Xi$ is a measure on $\Delta$ \cite{S00} with no atom at zero.   From the form in Eq \eqref{xilambdabk} the total coalescence rate given $b$ blocks can be shown to be
\be\label{lambdab}
  \lambda_b =  \int_\Delta \left(1 -  \sum_{\ell=0}^b \sum_{i_1\neq \ldots \neq i_\ell} \binom{b}{\ell}x_{i_1} \cdots x_{i_\ell}\left(1 - \sum_{j=1}^\infty x_j  \right)^{b-\ell}  \right) \frac{1}{\sum_{j=1}^\infty x_j^2 } \Xi(dx)  +  \binom{b}{2} 
\ee
\cite{S00}.   Since the $\Xi$-coalescent with generator $G$ as in Eq \eqref{Cconst} corresponds to a $\Xi$-coalescent with rate as the integral part in Eq \eqref{xilambdabk}  with $\Delta = (\psi/4, \psi/4, \psi/4, \psi/4, 0, \ldots)$ one can obtain the total rate of the multiple-merger part  from Eq \eqref{xilambdab} to be 
\begin{esplit}{xilambdab}
\lambda_b^{(\Xi)} & = \frac{4 }{\psi^2 } \left(1 - \sum_{\ell = 0 }^{b\wedge 4} \binom{b}{\ell} (4)_\ell \left( \frac{\psi}{4} \right)^\ell (1 - \psi)^{b-\ell } \right) 
\end{esplit}
%%
To add the Kingman-part, we  observe that the  rate at which one observes  large families  is $c\psi^2/4$, and the rate at which we observe small families is one.  Normalising the  full measure, $F$ say, with $F([0,1])$,  we see that the full process is a mixture process  behaving like a  $\Xi$-coalescent on $(\psi/4, \psi/4, \psi/4, \psi/4, 0,0, \ldots)$ with  probability  $(c\psi^2/4)/(1 + c\psi^2/4) $, and like  a Kingman-coalescent with probability     $1/(1 + c\psi^2/4) $.  The total coalescence rate is then, with $\lambda_b^{(\Xi)}$ from Eq \eqref{xilambdab}, 
\begin{esplit}{lambdab}
\lambda_b & = \frac{c\psi^2/4}{1 +  c\psi^2/4} \lambda_b^{(\Xi)}  +   \frac{1}{1 +  c\psi^2/4} \binom{b}{2}.
\end{esplit}
One can check that $\lambda_b$ in Eq \eqref{lambdab} gives the same rate as summing over  all the individual rates  $\lambda_{b;k_1, \ldots, k_r}$ when $  \mathcal{N}(b; k_1, \ldots, k_r ) $ is as in Eq \eqref{N} and  $C_{b; k_1, \ldots, k_r;s }$  as in Eq \eqref{Cconst} and 
\be
      \lambda_{b;k_1, \ldots, k_r} =    \mathcal{N}(b; k_1, \ldots, k_r ) \left( \frac{c\psi^2/4}{1 +  c\psi^2/4}C_{b; k_1, \ldots, k_r;s } +     \frac{ \one{r=1, k_1 = 2} }{1 +  c\psi^2/4}  \right),
\ee
$$
   \lambda_b =  \sum_{\substack{ b \ge  k_1 \ge  \ldots \ge k_r \ge 2 \\ 2 \le k_1 + \cdots + k_r \le b \\ 1 \le r \le 4  } } \lambda_{b;k_1, \ldots, k_r}.
$$
Observe that $ \lambda_{2} = \lambda_{2;2} = 1$.


Having an explicit formula for $\lambda_b$ allows us to sample the
time $T$ between events by drawing a random uniform $U$ on the unit
interval, and through the inverse transform $T = -(\log(1 -
U))/\lambda_b$.  We then draw another independent uniform $U$, 
if $U \le 1/(1 + c\psi^2/4)$ then we pick two blocks uniformly at random without
replacement and merge them.   If $U >   1/(1 + c\psi^2/4)  $    we    toss a  five-sided dice with sides $0,1,\ldots,4$, and   resp probabilities $1-\psi, \psi/4, \ldots, \psi/4$; the blocks with same  non-zero side are merged.  

@* {\bf Code}. 


Refer to file "xiesfsw.h" for the modules doing the computation.
Module \bi{Sconst} computes the number of mergers $\mathcal{N}(b;k_1,
\ldots, k_r)$ Eq\eqref{N} given merger sizes $k_1, \ldots,
k_r$. Module \bi{Cconst} computes $C_{b;k_1, \ldots, k_r }$ Eq
\eqref{Cconst}.  Given merger sizes $k_1, \ldots, k_r$, the
coalescence rate can be computed as $\lambda_{b;k_1, \ldots, k_r} =
\mathcal{N}(b;k_1, \ldots, k_r) C_{b;k_1, \ldots, k_r } $.  We then
need to enumerate all mergers, ordered in non-increasing order. 



@*1 {\bf The main function}. 



@C

@<Includes@>@#

int main( int argc, char * argv[] )
{@#

 /* run as {\tt ./a.out <sample size> <m> <cp> <psi> <alpha> }, where $m < 1$ returns results for the  $\Xi$-Beta-coalescent, $m> 1$ for the $\Xi$-Dirac coalescent,  $c$ and $psi$ are $c$ and $\psi$ of the $\Xi$-Dirac, and  $alpha$ is $\alpha$ of the $\Xi$-Beta coal */ 

 @q  static void printvarphi( int lauf,  double process, double cparameter, double psiparameter, double alphaparameter ) @>
 @q n, alpha, K, C @>
  printvarphi( atoi(argv[1]), atof(argv[2])/100.0, atof(argv[3]), atof(argv[4]) ); 

  gsl_rng_free( rngtype ) ; 
  
  return GSL_SUCCESS ;

/* end $main$ function */
}


@*1 {\bf Includes}. 

@<Includes@>=@#
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_combination.h>
#include "xiesfs.h"



@* {\bf Funding}. 

Funded by the \textcolor{blue}{\bf DFG}  {\bf SPP Priority Programme  1819}: \emph{Rapid Evolutionary Adaptation}.



@* {\bf References}. 


\bibliographystyle{plain}
\bibliography{refs.bib}




@
\end{document}%