

 /* obtain a seed out of thin air for the random number engine */
  std::random_device randomseed;
  /* Standard mersenne twister  random number engine seeded with rng() */
  std::mt19937_64 rng(randomseed());
 /* construct a uniform on the unit interval; call with uni01( rng ) */
std::uniform_real_distribution<double> uni01(0.0, 1.0);
/* construct a beta distribution; if boost beta is available */
/*    boost::random::beta_distribution<double> randombeta(1.0, 1.0); */
/* constuct a binomial */
std::binomial_distribution<unsigned> randombinom(1, .5) ;
std::gamma_distribution<double> rgamma(2.0, 1.0); 



/* returns a random gamma(a) variate */
static double samplegamma( const double a)
{
  rgamma.reset();  
  rgamma = std::gamma_distribution<double>(a, 1.0);
  return( rgamma( rng ) ); 
}



/* returns a random beta(a,b) variate */
static double bg(const double a, const double b)
{

  double x = samplegamma(a);
  return( x/(x + samplegamma(b) ) );
}


static double xreject( const double a, const double u, const double M)
{
/* need boost library */
 double x;

 /* no random beta in boost library 
 randombeta.reset();
 randombeta =  boost::random::beta_distribution<double>(2.-a, a);
 x = randombeta(rng); 
 */
 x = bg( 2.0 - a, a); 
  /* x =  boost::math::ibeta_inv(2.-a, a, u); */

 while( (x > M) || ( x <= 0.) ){
   /*   x = randombeta(rng); } */
   x = bg( 2.0 - a, a); }
 /* x =   boost::math::ibeta_inv(2.-a, a, u) ;} */

  assert( x > 0. );
  assert( x <= M);
  return( x);
}


static double itime( const double n, const double a,   const double M,  const double b,   double* x, double* Sjp)
{
  /* $n$ is  number of  blocks; $a$ is $\alpha$; $M=M$, $b$ is exponential growth parameter */
double timi = 0. ;
double ct = n*(n - 1.)/2. ;

/*
  timi = timi + ( -log(1. - uni01(rng) )/(ct) ); 
 x[0] = xreject( a, uni01(rng), M); 
*/

if ( a < 2. ){ 
do{
  timi = timi + ( -log(1. - uni01(rng) )/(ct) ); 
  x[0] = xreject( a, uni01(rng), M); }
  while ( uni01(rng) >  (((1. - pow(1.-x[0], n) - (n*x[0]*pow(1.-x[0], n-1.)))/(x[0]*x[0]))/(ct) ) );
}
else{
   timi =  ( -log(1. - uni01(rng) )/( ct ) );} 

  double Sj = (a < 2 ? 0. : (b > 0. ? log( exp(b*Sjp[0]) - (b*log(uni01(rng)) / ( ct )))/b : Sjp[0] + timi ) ); 

  double svar =   (a < 2. ? timi : Sj - Sjp[0]) ;
  Sjp[0] = Sj ;
  return( svar ); 
}

static unsigned int samplek( const double n,  const double x )
{

  unsigned k;
  randombinom.reset(); 
   randombinom = std::binomial_distribution<unsigned>(n-2., x) ;
   k = 2 + (x > 0.0 ? randombinom(rng) : 0) ;
   while( uni01(rng) > 2./( (double)(k*(k - 1)) ) ){
     k = 2 + (x > 0.0 ? randombinom( rng ) : 0); }
 
 assert( k > 1);
 assert( k <= (unsigned)n);
 
  return(k);
}



static bool inRange( const unsigned x, const unsigned l, const unsigned h)
{

  return(  (x - l) <=  (h-l));  
  /*  return( x >= Imin ? (x <= Imax ? 0 : ( x >= Jmin ? (x <= Jmax ? 1 : 3) : 3)) : 3);  */

}


static unsigned group( const unsigned x)
{

  /* return( inRange( x, Imin, Imax) ? 0 : (inRange(x, Jmin, Jmax) ? 1 : 3) ); */
  return (1);

}


/*
extern "C" unsigned lambdaijk( const unsigned n, const double a, const double K, const double b, const double theta, const unsigned Imin, const unsigned Imax, const unsigned Jmin, const unsigned Jmax,  const unsigned Kmin, const unsigned Kmax,  const unsigned trials, const unsigned numer )
g++ -Wall -fPIC -shared -o clambdabetapython.so clambdabetapython.cpp
import numpy as np
np.loadtxt('egresskra123', usecols=range(0,3))
*/

/* change between a and b depending on process */
extern "C" void lambdaijk(const double a, const double theta, double * ebi)
{
  /*   NumericMatrix& res */
 /* std::cout << n << ' ' << a << ' ' << K << ' ' << b << '\n' ; */
  /* set dimenions of matrix `res` 
  res.attr("dim") = Dimension(trials, 3) ;
  */
  const double b = 0;

  /* south-coast data is 79 fish; thistilfj data is 90 fish */
  const unsigned n = 79;
  const unsigned numer=123; 
  const double K = 0;
  /* const double theta = 0; */
  const unsigned Imin = 1;
  const unsigned Imax = 1;
  const unsigned Jmin = 15;
  const unsigned Jmax = n-1;
  const unsigned Kmin = Imin;
  const unsigned Kmax = Jmax;
  const unsigned trials = 12;
  /*   const unsigned numer = 2323; */

  
  std::string resskra  = "egresskra";
  resskra = resskra + std::to_string(numer);
  
  std::vector<unsigned> v;
  /*  std::vector<double> ebi; */
  std::list<unsigned> tre;
  std::list<unsigned>::iterator it;

  unsigned r, size, kblocks, g;
  double   timi, z;
  /*
   const double   a = atof(argv[2]);
   const double K = atof(argv[3]);
   const double b = atof(argv[4]); 
   const double theta = atof(argv[5]);
   */
   const double M = (a < 2. ? (K > 0. ? K/(K + 1. + (pow(2,1.-a)/(a-1.))) : 1.) : 1.);
   double *x = (double *)calloc(1, sizeof(double)); 
   double *Sjp = (double *)calloc(1, sizeof(double)); 

   /* generating a poisson random number generator */
    std::poisson_distribution<unsigned> rpois(1.0) ;
    
   /* generating a beta random number generator using boost */

  /* generating an exponential dist */

  /* initialise spectrum */
    /* ebi.assign( (unsigned)n  + 1, 0. ); 
    ebi.assign( 4, 0. ); 
 r = 0; 
    */
    /*
       # pragma omp parallel for
 for( r =0; r < trials; ++r ){
    */
    assert( n > 1 );
    tre.assign( n, 1.);
    ebi[0] = 0.0;
    ebi[1] = 0.0;
    ebi[2] = 0.0;
    ebi[3] = 0.0;
   /* ebi.assign( 4, 0. );  */
    Sjp[0] = 0. ;
   
    while( (unsigned)tre.size() > 1){ 
    /*
       for( auto& x: tre){ std::cout << ' ' << x; }
        std::cout << '\n';
        */
   /* sample time and $x$ */
     
     timi = itime( (double)tre.size(), a, M, b, x, Sjp);

  /* given sampled time update spectrum */
     if (theta > 0. ){
       rpois.reset();
       rpois = std::poisson_distribution<unsigned>(theta * timi);}
     for( unsigned& y: tre){
       z = (theta > 0. ? (double)rpois(rng) : timi) ;
       g = (inRange(y, Kmin, Kmax) ? 2 : 3) ;
       ebi[g] = ebi[g] + z;
       g = (inRange(y, Imin, Imax) ? 0 : 3) ;
       ebi[g] = ebi[g] + z;
       g = (inRange(y, Jmin, Jmax) ? 1 : 3) ;
       ebi[g] = ebi[g] + z;}
      /*
      ebi[(y < Kmin ? 3 :(y > Kmax ? 3 : 0))] = ebi[(y < Kmin ? 3 :(y > Kmax ? 3 : 0))] + z ;
      ebi[ (y < Imin ? 3 : (y > Imax ? 3 : 1))] = ebi[ (y < Imin ? 3 : (y > Imax ? 3 : 1))] + z ;
      ebi[ (y < Jmin ? 3 : (y > Jmax ? 3 : 2))] = ebi[ (y < Jmin ? 3 : (y > Jmax ? 3 : 2))] + z ; }
      /*
   if( (unsigned)tre.size() > 2){
    /* first sample  blocks to merge */
    /* v.resize( (unsigned)tre.size(), 0); */
     kblocks = samplek( (double)tre.size(),  x[0]);

     v.clear();
     v.assign( tre.begin(), tre.end()); 
     std::shuffle( v.begin(), v.end(), rng );
     tre.assign( v.begin(), v.end());
    
     size = 0 ; 
     for( it = tre.begin(); it != std::next( tre.begin(), kblocks); ++it){
       assert( *it > 0 );
       size = size + (*it);
    /* label the block pointed to by $it$ to be removed */
       *it = 0; }
     tre.push_back(size);
     tre.remove(0); }
   /* generated one tree */
   tre.clear();
   std::list<unsigned>().swap(tre);
   v.clear();
   std::vector<unsigned>().swap(v);
   free( Sjp);
   free( x);
   ebi[0] = (ebi[2] > 0.0 ? ebi[0]/ebi[2] : 0.0);
   ebi[1] = (ebi[2] > 0.0 ? ebi[1]/ebi[2] : 0.0);
   /*
   else{
     tre.resize(1); }
    generated one tree */
   /*
   std::ofstream f;
   f.open(resskra, std::ios::app );
   f <<  (ebi[2] > 0. ? ebi[1]/ebi[2] : 0.) << ' ' << (ebi[2] > 0. ? ebi[1]/ebi[2] : 0.)  << ' ' << ebi[2] << std::endl ;
   f.close();
   */ 
/* print out average of spectrum */
       /* for(  double& x: ebi){ std::cout << x/ebi[0] << '\n' ; } */
       /*  std::cout << ebi[1]/ebi[0] << ' ' << ebi[2]/ebi[0] <<  '\n' ; */
  /*  
   FILE * f;
  f = fopen(resskra, "a");
  fprintf(f,  "%g %g\n", ebi[1]/ebi[0], ebi[2]/ebi[0]  );
  fclose(f); 
  std::ofstream f;
  f.open(resskra, std::ios::app );
  f <<  (ebi[0] > 0. ? ebi[1]/ebi[0] : 0.) << ' ' << (ebi[0] > 0. ? ebi[2]/ebi[0] : 0.)  << ' ' << ebi[0] << '\n' ;
  f.close();
  */
  /*  if working with matrix `res` : 
res(r-1, 0) =  ebi[1]/ebi[0] ;
  res(r-1, 1) =  ebi[2]/ebi[0];
  res(r-1, 2) =  ebi[0];
  */ 
  /* close while loop over $trials$ */
}
