
double  **alloc_2d_array(  int n)
{
  //
  double ** m = NULL ;
  m = (double **)calloc( n + 1, sizeof( double *)) ;
  int i ;
  for( i = 0 ; i <= n ; i++ ){
    m[i] = (double *)calloc( n+1, sizeof( double)) ; } 
    //    printf("%d  %d\n", i, sizeof(m[i]) ) ; }
  ////
  return( m ) ;
}

void free_2d_array( double **a, int n )
{
  //
  int i ;
  for( i = 0 ; i <= n ; i++ ){
    free( a[i] ); }
  free( a ) ; 
  // a = NULL ;
}




void prentagslfylki( int z,  gsl_matrix * m )
{

  int x, y ;
  for( x = 1 ; x <= z ; x++ ){
    for( y = 1 ; y <= z ; y++ ){
      printf("%g  ", gsl_matrix_get( m, x, y) ) ; }
    printf("\n") ; }
}



void syna( int k,  int y[] )
{

  int z ;
  for( z = 1 ; z <= k ; z++ ){
    printf( "%d ", y[z] ) ; }
  printf("\n") ; 

}


/*
double  **alloc_2d_array(  int n)
{
  //
  double ** m = NULL ;
  m = (double **)calloc( n + 1, sizeof( double *)) ;
  int i ;
  for( i = 0 ; i <= n ; i++ ){
    m[i] = (double *)calloc( n+1, sizeof( double)) ; }
  ////
  return( m ) ;
}



void free_2d_array( double **a, int n )
{
  //
  int i ;
  for( i = 0 ; i <= n ; i++ ){
    free( a[i] ); }
  free( a ) ; a = NULL ;
}
*/


static double ff( int a, int m )
{
  
  int svar = 1 ;
  int i ;
  if( m > 0 ){
  for( i = 0 ; i < m ; i++ ){
    svar = svar * ( a - i ) ; }  }
  else{
    svar = 1 ; }
  
  return( (double)svar ) ;

}





static double  Sconst( int n,  int x[] )
{
  // number of mergers  with same sizes as in x; simply reordering the mergers
  // allowing at most 4 sim mergers
  // x[] should be non-increasing  order
  
  double alpha  = 1. ;
  
  alpha =  ( x[4] > 1 ?    ( ((x[1] == x[2]) && (x[1] == x[3])) && (x[1] == x[4]) ? 24. : ( ((x[1] == x[2]) && (x[1] == x[3])) || ((x[2] == x[3]) && (x[2] == x[4])) ? 6. : (              (x[1] == x[2]) && (x[3] == x[4]) ? 4. : (((x[1] == x[2]) || (x[2] == x[3])) || (x[3] == x[4]) ? 2. : 1.)  ) ) ) :  (  x[3] > 1 ?                   ( ((x[1] == x[2]) && (x[1] == x[3])) ? 6. : ((x[1] == x[2]) || (x[2] == x[3]) ? 2. : 1. ) ) : ( x[2] > 1 ? ( x[1] == x[2] ? 2. : 1.) : 1. ) ) ) ;
  

  //  return( ff(n, ((x[1] > 1 ? x[1] : 0) + (x[2] > 1 ? x[2] : 0)  + (x[3] > 1 ? x[3] : 0)  +  (x[4] > 1 ? x[4] : 0))  ) / ( alpha * (x[1] > 1 ? gsl_sf_fact(x[1]) : 1.) * (x[2] > 1 ? gsl_sf_fact(x[2]) : 1.) * (x[3] > 1 ? gsl_sf_fact(x[3]) : 1.) * (x[4] > 1 ? gsl_sf_fact(x[4]) : 1.) ) )  ;
  // printf("alpha  %g\n", alpha ) ;
  return (   gsl_sf_choose( x[1] + x[2], x[2]) * gsl_sf_choose( x[1] + x[2] + x[3],  x[3]) * gsl_sf_choose( x[1] + x[2] + x[3] + x[4], x[4] ) * gsl_sf_choose( n, n -  x[1] - x[2] - x[3] - x[4]) / alpha ) ;
  //    return( gsl_sf_fact( n) /(  gsl_sf_fact( x[1]) * gsl_sf_fact(x[2]) * gsl_sf_fact( x[3]) * gsl_sf_fact(x[4]) * gsl_sf_fact( n -  (x[1] > 1 ? x[1] : 0) -  (x[2] > 1 ? x[2] : 0) -  (x[3] > 1 ? x[3] : 0) -  (x[4] > 1 ? x[4] : 0) ) * alpha ) ) ;

}





/* computes the rate lambda(Xi) of a Xi-coal on (psi/4, \psi/4, \psi/4, \psi/4, 0,0,...)*/
static double Cconst( int n,  double psi, int x[]  )
{

  
  // psi is \psi in the Dirac coalescent
  int r, l, s;
  double sum = 0. ;
  r = (x[1] > 1 ? 1 : 0) + (x[2] > 1 ? 1 : 0) + (x[3] > 1 ? 1 : 0) + (x[4] > 1 ? 1 : 0) ;
  
  s = n -  (x[1] > 1 ? x[1] : 0)  -  (x[2] > 1 ? x[2] : 0)  -  (x[3] > 1 ? x[3] : 0)  -  (x[4] > 1 ? x[4] : 0) ;


  for( l = 0 ; l <= (s < 4-r ? s : 4-r) ; l++){
    sum  =  sum  +  (gsl_sf_choose( s, l) * ff(4,  r+l)*gsl_sf_pow_int( 1. - psi, s-l) * gsl_sf_pow_int(psi, n-s + l)/gsl_sf_pow_int( 4, n-s+l)) ; }

  return( 4.* sum/(psi*psi)  ) ; 

  // the coalescence rate \lambda_{n;k} = Sconst * Cconst
}




/* compute the coal rate lambda_{b; k} where b are active blocks and k a vector of merger sizes; computes only the  rate for the specific ordering of the merger sizes, does not sum over all orderings */
static double  Diraccoalrate( int b,  double cstiki,  double psistiki,   int vki[] )
{

  /* if K > 0, then kingman part included in Dirac coalescent
  return(  cstiki  *  Cconst( b, psistiki, vki)  +  (K > 0. ? (vki[1] == 2 && (((vki[2] < 2) && (vki[3] < 2)) && (vki[4] < 2)) ? 1. : 0.) : 0.) ) ;
  */
  double prob =  (cstiki*psistiki*psistiki/4.)/(1. +   (cstiki*psistiki*psistiki/4.)); 
  return(  (prob*Cconst( b, psistiki, vki))  +    (vki[1] == 2 && (((vki[2] < 2) && (vki[3] < 2)) && (vki[4] < 2)) ? 1. : 0.)*(1. - prob) ); 
}



static double Betacoalrate( int b, double alpha, const double K,  const double Cconst,   int vki[] )
{

  // returns lambda_b;k for given ordered merger specified by vki 
  // see eq 3.3 in Cronj thesis, page 69
  /* when 0 < alpha_1 < 1 and alpha_2 \ge 2 then E[X_1] ->  2  +  pow(2, a)*pow(3, 1-a)/(a - 1) with a  = alpha_2 \ge 2 
     and at most one   individual gets  alpha_1 each time */  
  const double M = (K > 0.0 ? ( alpha < 1.0 ?  K/(K + 1.0 + (Cconst * alpha* pow( K, 1.0 - alpha))) : K/(1.0 + K) ) : 1.0) ;
  int l ;
  int r =  (vki[1] > 1 ? 1 : 0)  +  (vki[2] > 1 ? 1 : 0) + (vki[3] > 1 ? 1 : 0) + (vki[4] > 1 ? 1 : 0) ;
  int s =  b - (vki[1] > 1 ? vki[1] : 0)  -  (vki[2] > 1 ? vki[2] : 0) - (vki[3] > 1 ? vki[3] : 0) - (vki[4] > 1 ? vki[4] : 0) ;
  
  double svar = 0. ; 
  for( l = 0 ; l <= (s < 4-r ? s : 4-r) ; l++ ){

    
   svar +=  Cconst* gsl_sf_choose(s,l) * ff(4, r+l) * (K > 0.0 ? (gsl_sf_beta_inc(  ((double)(l+b-s)) - alpha,  (double)(s-l) + alpha, M) * gsl_sf_beta(  ((double)(l+b-s)) - alpha,  (double)(s-l) + alpha )/gsl_sf_beta_inc( 2.0 - alpha, alpha, M)) : gsl_sf_beta( ((double)(l+b-s)) - alpha, (double)(s-l) + alpha)/gsl_sf_beta(2. - alpha, alpha)) / gsl_sf_pow_int(  4. , l + b - s); }
     
    /*
    svar += gsl_sf_choose(s,l) * ff(4, r+l) * gsl_sf_beta( ((double)(l+b-s)) - alpha, (double)(s-l) + alpha)/( gsl_sf_beta(2. - alpha, alpha) * gsl_sf_pow_int(  4. , l + b - s ) ) ; }
    */
  return( svar ) ;

}

/* computes the lambda(Xi) rate for a Xi-measure on (psi/4, psi/4, psi/4, psi/4, 0,0,...) */
double XiDiracrate( int b,  double  psistiki,  int vki[] )
{

   int l ;
  int r =  (vki[1] > 1 ? 1 : 0)  +  (vki[2] > 1 ? 1 : 0) + (vki[3] > 1 ? 1 : 0) + (vki[4] > 1 ? 1 : 0) ;
  int s =  b - (vki[1] > 1 ? vki[1] : 0)  -  (vki[2] > 1 ? vki[2] : 0) - (vki[3] > 1 ? vki[3] : 0) - (vki[4] > 1 ? vki[4] : 0) ;
  
  double svar = 0. ; 
  for( l = 0 ; l <= (s < 4-r ? s : 4-r) ; l++ ){
    svar  =   svar  +  gsl_sf_choose(s,l) * ff(4, r+l) * ( gsl_sf_pow_int( psistiki, l + b - s - 2)*gsl_sf_pow_int(1. -  psistiki,  b - (l+b-s)))/gsl_sf_pow_int(  4. , l + b - s )  ; }

  return( svar ) ;

}





static double coalrate( int in,  double coalstikar[],  int Vki[] )
{
  //// if coalstikar[0] = 0 ( < 1 ) then Beta model
  //// if coalstikar[0] = 1 ( > 0 ) then Dirac model
  //// 
  ////  double Betacoalrate( int b, double alpha,  int vki[] )
  ////  double  Diraccoalrate( int b,   double cstiki,  double psistiki,  int vki[] )
  ////
  
  //return( coalstikar[0] < 1 ? Betacoalrate( in, coalstikar[1],  Vki ) :  Diraccoalrate( in, coalstikar[1], coalstikar[2],   Vki ) ) ;
  /*  return( coalstikar[0] < 1 ? Betacoalrate( in, coalstikar[1],  Vki ) :  XiDiracrate( in, coalstikar[1],  Vki ) ) ; */

  /* if Xi-Dirac then coalstikar[1] is c, coalstikar[2] is psi */
  /*  return( coalstikar[0] < 1 ? Betacoalrate( in, coalstikar[3],  Vki ) :  Diraccoalrate( in, coalstikar[1], coalstikar[2],  Vki ) ) ;  */
  /* sample size;  alpha; K; C; vk */
  return( Betacoalrate( in, coalstikar[1], coalstikar[2], coalstikar[3],  Vki ) );  
}




void partzs2( int n  )
{

  int x[n+1], i, m, h, r, j ;
  i = m = h = r = j =  0 ;
  
  for( i = 1 ; i <= n ; i++ ){
    x[i] = 1 ; }

  x[0] = -1 ; 
  x[1] = 2 ;
  h = 1 ;
  
  m =  n - 1 ;

  while( (x[1] != n) )  {
    if( m - h > 1 ){
      h = h + 1 ;
      x[h] = 2 ;
      m = m - 1 ; }
    else{
      j = m - 2 ;
      while( x[j] == x[m-1] ){
	x[j] = 1 ; 
	j  =  j - 1 ; }
      h = j + 1 ;
      x[h] = x[m-1] + 1 ; 
      r = x[m] + ( x[m-1]*(m-h-1) );
      x[m] = 1 ;
      x[m-1] = (m-h > 1 ? 1 : x[m-1] ) ;
      // original code
        m = h + r - 1 ; }
    //syna( m, x ) ;
    
  }
      
}



void mergers( int N, double psi )
{

  int k1, k2, k3, k4 ;

  
  for( k1 = 2 ; k1 <= N  ; k1++ ){
    //printf( "%d\n", k1 ) ; 
  }

  if( N > 3 ){
    for( k1 = 2; k1 <=  N - 2 ; k1++){
      for( k2 = 2 ; k2 <= (N-k1 < k1 ? N-k1 : k1) ; k2++ ){
	// printf( "%d %d\n", k1, k2 ) ; 
      } } }

  
  if( N > 5 ){
    for( k1 = 2 ; k1 <= N - 4 ; k1++ ){
      for( k2 = 2 ; k2 <= (N-k1-2 < k1 ? N - k1 - 2 : k1 ); k2++ ){
	for( k3 = 2 ; k3 <= (N-k1-k2 < k2 ? N-k1-k2 : k2); k3++ ){
	  //  printf( "%d %d %d\n", k1, k2, k3 ) ; 
	} } } }
  
  
  if( N > 7){
    for( k1 = 2 ; k1 <= N-6 ; k1++ ){
      for( k2 = 2 ; k2 <= (N-k1-4 < k1 ? N-k1-4 : k1); k2++ ){
	for( k3 = 2 ; k3 <= (N-k1-k2 - 2 < k2 ? N-k1-k2 - 2 : k2 ); k3++ ){
	  for( k4 = 2 ; k4 <= (N-k1-k2-k3 < k3 ? N-k1-k2-k3 : k3); k4++){
	    //  printf( "%d %d %d %d\n", k1, k2, k3, k4 ) ; 
	  } } } } }

}








void Mqp( int N, double icoals[],   gsl_matrix * Mq, gsl_matrix * Mp )
{

  // compute rate of going from N to m ;
  // double  Sconst( int n,  int x[] )
  // 
  ////  static double coalrate( int in,  double coalstikar[],  int Vki[] )
  ////  double coalrate( int in,  double coalstikar[],  int Vki[] )

  int vk[] = {0,0,0,0,0} ;
  double rate = 0. ;
  int n, k1, k2, k3, k4 ;
  for( n = 2 ; n <= N ; n++ ){

    for( k1 = 2 ; k1 <= n ; k1++ ){
      vk[1] = k1 ; vk[2] = vk[3] = vk[4] = 0 ;
      // rate = gsl_sf_choose( n, k1) * ( Cconst( n, psi, vk) + ( k1 == 2 ? 1. : 0. ) ) ;
      // icoals should be  {0, alpha, K, C}
      rate =    Sconst( n,  vk ) *  coalrate( n, icoals,  vk ) ;
      ////
      gsl_matrix_set( Mq, n,  n - k1 + 1,  gsl_matrix_get( Mq, n, n-k1 + 1) +  rate) ; 
      gsl_matrix_set( Mq, n, n,  gsl_matrix_get( Mq,n,n)  +  rate) ; 
      gsl_matrix_set( Mp, n, n - k1 + 1, gsl_matrix_get( Mp, n, n - k1 + 1) + rate) ; }

  if( (n > 3) ){
    for( k1 = 2; k1 <=  n - 2 ; k1++){
      for( k2 = 2 ; k2 <= (n-k1 < k1 ? n-k1 : k1) ; k2++ ){
	vk[1] = k1 ; vk[2] = k2 ; vk[3] = vk[4] = 0 ; 
	// rate =  Sconst( n, vk) * Cconst( n, psi, vk ) ;
	rate =    Sconst( n,  vk ) * coalrate( n, icoals,  vk ) ;
	gsl_matrix_set( Mq, n,  n - k1 - k2 + 2,  gsl_matrix_get( Mq, n, n - k1 - k2 + 2) + rate ) ;
	gsl_matrix_set( Mq, n, n, gsl_matrix_get( Mq, n,n)  +  rate ) ;
	gsl_matrix_set( Mp, n, n - k1 - k2 + 2,  gsl_matrix_get( Mp, n, n - k1 - k2 + 2) +  rate) ; 
      } } }

  
  if( n > 5 ){
    for( k1 = 2 ; k1 <= n - 4 ; k1++ ){
      for( k2 = 2 ; k2 <= (n-k1-2 < k1 ? n - k1 - 2 : k1 ); k2++ ){
	for( k3 = 2 ; k3 <= (n-k1-k2 < k2 ? n-k1-k2 : k2); k3++ ){
	  vk[1] = k1 ; vk[2] = k2 ; vk[3] = k3 ;  vk[4] = 0 ;
	  // rate =  Sconst( n, vk) * Cconst( n, psi, vk ) ;
	  rate =   Sconst( n,  vk ) * coalrate( n, icoals,  vk ) ;

	  gsl_matrix_set( Mq, n,  n - k1 - k2 - k3 + 3,  gsl_matrix_get( Mq, n, n - k1 - k2 - k3 + 3) + rate ) ;
	  gsl_matrix_set( Mq, n, n, gsl_matrix_get( Mq, n,n)  +  rate ) ;
	  gsl_matrix_set( Mp, n, n - k1 - k2 - k3 + 3,  gsl_matrix_get( Mp, n, n - k1 - k2 - k3 + 3) +  rate) ;
	} } } }
  
  
  if( n > 7){
    for( k1 = 2 ; k1 <= n-6 ; k1++ ){
      for( k2 = 2 ; k2 <= (n-k1-4 < k1 ? n-k1-4 : k1); k2++ ){
	for( k3 = 2 ; k3 <= (n-k1-k2 - 2 < k2 ? n-k1-k2 - 2 : k2 ); k3++ ){
	  for( k4 = 2 ; k4 <= (n-k1-k2-k3 < k3 ? n-k1-k2-k3 : k3); k4++){
	    vk[1] = k1 ; vk[2] = k2 ; vk[3] = k3 ; vk[4] = k4 ; 
	    //// rate =  Sconst( n, vk) * Cconst( n, psi, vk ) ;
	    rate  =   Sconst( n,  vk ) * coalrate( n, icoals,  vk ) ;

	  gsl_matrix_set( Mq, n,  n - k1 - k2 - k3 - k4 + 4,  gsl_matrix_get( Mq, n, n - k1 - k2 - k3 - k4 + 4) + rate ) ;
	  gsl_matrix_set( Mq, n, n, gsl_matrix_get( Mq, n,n)  +  rate ) ;
	  gsl_matrix_set( Mp, n, n - k1 - k2 - k3 - k4 + 4,  gsl_matrix_get( Mp, n, n - k1 - k2 - k3 - k4 + 4) +  rate) ;
	  } } } } } }
  
  // normalize the entries of Mp
   
  for( n = 2 ; n <= N ; n++ ){
    for( k1 = 1 ; k1 < n ; k1++){
      gsl_matrix_set( Mp, n, k1, gsl_matrix_get( Mp, n, k1)/gsl_matrix_get(Mq, n,n) ) ; 
      gsl_matrix_set( Mp, n,n, gsl_matrix_get( Mq, n, n) ) ; }}
  

}





int  draw_number_blocks(  int b,  double * M, gsl_rng * r )
{

  int svar = 0 ;

  gsl_ran_discrete_t * P =  gsl_ran_discrete_preproc( b,  M ) ; 
  svar =  gsl_ran_discrete( r,  P ) ; 
  gsl_ran_discrete_free( P ) ; 
  
  // svar is number of active blocks after jump
  return( b - svar + 1 ) ; 

}



static void totalrate( int N,  double cstikar[],  gsl_vector * vln )
{

  // double coalrate( int in,  double coalstikar[],  int Vki[] )
  //// if coalstikar[0] = 0 ( < 1 ) then Beta model
  //// if coalstikar[0] = 1 ( > 0 ) then Dirac model
  //// double  Sconst( int n,  int x[] )

  // return total coal rate vln[n] for n = 2 to N
  // Diraccoalrate( n, psi, cstiki,  iK,  vk)

  int vk[] = {0,0,0,0,0} ;
  int n, k1, k2, k3, k4 ;
  for( n = 2 ; n <= N ; n++ ){

    for( k1 = 2 ; k1 <= n ; k1++ ){
      vk[1] = k1 ; vk[2] = vk[3] = vk[4] = 0 ;
      gsl_vector_set( vln, n,  gsl_vector_get( vln, n)  +   (Sconst(n, vk) *  coalrate( n, cstikar, vk)) ) ; }

  if( (n > 3) ){
    for( k1 = 2; k1 <=  n - 2 ; k1++){
      for( k2 = 2 ; k2 <= (n-k1 < k1 ? n-k1 : k1) ; k2++ ){
	vk[1] = k1 ; vk[2] = k2 ; vk[3] = vk[4] = 0 ;
	// rate =  Sconst( n, vk) * cstiki * Cconst( n, psi, vk ) ;
	gsl_vector_set( vln, n, gsl_vector_get( vln,  n)  +   (Sconst(n, vk) *coalrate( n, cstikar, vk)) ) ; 
      } } }

  
  if( n > 5 ){
    for( k1 = 2 ; k1 <= n - 4 ; k1++ ){
      for( k2 = 2 ; k2 <= (n-k1-2 < k1 ? n - k1 - 2 : k1 ); k2++ ){
	for( k3 = 2 ; k3 <= (n-k1-k2 < k2 ? n-k1-k2 : k2); k3++ ){
	  vk[1] = k1 ; vk[2] = k2 ; vk[3] = k3 ; vk[4] = 0 ; 
	  //  rate =  Sconst( n, vk) * cstiki * Cconst( n, psi, vk ) ;
	  gsl_vector_set( vln, n,  gsl_vector_get( vln,  n)  +   (Sconst(n, vk) *coalrate( n, cstikar, vk)) ) ; 
	} } } }
  
  
  if( n > 7){
    for( k1 = 2 ; k1 <= n-6 ; k1++ ){
      for( k2 = 2 ; k2 <= (n-k1-4 < k1 ? n-k1-4 : k1); k2++ ){
	for( k3 = 2 ; k3 <= (n-k1-k2 - 2 < k2 ? n-k1-k2 - 2 : k2 ); k3++ ){
	  for( k4 = 2 ; k4 <= (n-k1-k2-k3 < k3 ? n-k1-k2-k3 : k3); k4++){
	    vk[1] = k1 ; vk[2] = k2 ; vk[3] = k3 ; vk[4] = k4 ; 
	    // rate =  Sconst( n, vk) * cstiki * Cconst( n, psi, vk ) ;
	    gsl_vector_set( vln, n,  gsl_vector_get( vln,  n)  +   (Sconst(n, vk) *coalrate( n, cstikar, vk)) ) ; 
	  } } } } } }

}


void Betatotalrate( int N,  const double  ialpha, const double iK, const double iCconst,     gsl_vector * vln )
{

  // return total rate associated with Beta coalescents
  /* 
 static double Betacoalrate( int b, double alpha, const double K,  const double Cconst,   int vki[] )
  */

  // double Betacoalrate( int b, double alpha,  int vki[] )
  //  Sconst( int n,  int x[] )

  int vk[] = {0, 0, 0, 0, 0} ;
  int n, k1, k2, k3, k4 ;
  for( n = 2 ; n <= N ; n++ ){

    for( k1 = 2 ; k1 <= n ; k1++ ){
      vk[1] = k1 ; vk[2] = vk[3] = vk[4] = 0 ;
      gsl_vector_set( vln, n,  gsl_vector_get( vln, n)  +    Sconst(n, vk) * Betacoalrate( n, ialpha, iK, iCconst,   vk) ); }

  if( (n > 3) ){
    for( k1 = 2; k1 <=  n - 2 ; k1++){
      for( k2 = 2 ; k2 <= (n-k1 < k1 ? n-k1 : k1) ; k2++ ){
	vk[1] = k1 ; vk[2] = k2 ; vk[3] = vk[4] = 0 ;
	// rate =  Sconst( n, vk) * cstiki * Cconst( n, psi, vk ) ;
	gsl_vector_set( vln, n, gsl_vector_get( vln,  n)  +    Sconst(n, vk) * Betacoalrate( n, ialpha, iK, iCconst,   vk)  ) ; 
      } } }

  
  if( n > 5 ){
    for( k1 = 2 ; k1 <= n - 4 ; k1++ ){
      for( k2 = 2 ; k2 <= (n-k1-2 < k1 ? n - k1 - 2 : k1 ); k2++ ){
	for( k3 = 2 ; k3 <= (n-k1-k2 < k2 ? n-k1-k2 : k2); k3++ ){
	  vk[1] = k1 ; vk[2] = k2 ; vk[3] = k3 ; vk[4] = 0 ; 
	  //  rate =  Sconst( n, vk) * cstiki * Cconst( n, psi, vk ) ;
	  gsl_vector_set( vln, n,  gsl_vector_get( vln,  n)  +    Sconst(n, vk) * Betacoalrate( n, ialpha, iK, iCconst,   vk)  ) ; 
	} } } }
  
  
  if( n > 7){
    for( k1 = 2 ; k1 <= n-6 ; k1++ ){
      for( k2 = 2 ; k2 <= (n-k1-4 < k1 ? n-k1-4 : k1); k2++ ){
	for( k3 = 2 ; k3 <= (n-k1-k2 - 2 < k2 ? n-k1-k2 - 2 : k2 ); k3++ ){
	  for( k4 = 2 ; k4 <= (n-k1-k2-k3 < k3 ? n-k1-k2-k3 : k3); k4++){
	    vk[1] = k1 ; vk[2] = k2 ; vk[3] = k3 ; vk[4] = k4 ; 
	    // rate =  Sconst( n, vk) * cstiki * Cconst( n, psi, vk ) ;
	    gsl_vector_set( vln, n,  gsl_vector_get( vln,  n)  +    Sconst(n, vk) * Betacoalrate( n, ialpha, iK, iCconst,   vk)  ) ; 
	  } } } } } }

}








void gnm( int N,   gsl_matrix * Mpp,   gsl_matrix * Mgnm )
{


  int n, m, k ; 
  double s = 0. ;
  for( n = 2 ; n <= N ; n++){
    gsl_matrix_set( Mgnm, n, n, 1. / gsl_matrix_get( Mpp, n, n) ) ;
    for( m = 2 ; m < n ; m++ ){
      s = 0. ;
      for( k = m ; k < n ; k++ ){
	s  =  s  +  (gsl_matrix_get( Mpp, n, k) * gsl_matrix_get( Mgnm, k, m )) ; }
      gsl_matrix_set( Mgnm, n, m, s ) ; } }

}



double tw_switch_number_zero(  int b, int m, double db, double dm,  int vki[],   gsl_matrix * mr  )
{
  /* k_1 > k_2 */
  double summa = 0. ;

  /*   (0, 0, b ) */

  
  summa = summa +   (b < m - 1 ? (dm - db)*(dm - db - 1.)*gsl_matrix_get( mr, m, b)/(dm*(dm-1.)) : 0.);

  
  /*  (1, 0, - ) */
  summa  =  summa  +       (b > vki[1] - 1 ?  (db - (double)vki[1]  +  1.)*(dm - db +  ((double)vki[1]) - 1.) * gsl_matrix_get( mr, m,  b - vki[1] + 1)/(dm*(dm - 1.)) : 0. );

  /*  (0, 1, - ) */

  summa = summa   +    (b > vki[2] - 1 ?  (db - (double)vki[2]  +  1.)*(dm - db +  ((double)vki[2]) - 1.) * gsl_matrix_get( mr, m,  b - vki[2] + 1)/(dm*(dm - 1.)) : 0.);


  /* (1, 1,  - )   */
  summa = summa   +    (b > vki[1] + vki[2] - 1 ?  (db -  ((double)(vki[1] + vki[2]))  + 2.)*(db -  ((double)(vki[1]  + vki[2]))  + 1.)*gsl_matrix_get( mr, m, b - vki[1] - vki[2] + 2 )/(dm*(dm-1.)) : 0.);

  return ( summa ) ;

}



double tw_switch_number_one(  int b, int m, double db, double dm,  int vki[],   gsl_matrix * mr  )
{
  /* k_1  ==  k_2 */
  double summa = 0. ;

  /*   (0, b ) */
  summa = summa +   (b < m - 1 ? (dm - db)*(dm - db - 1.)*gsl_matrix_get( mr, m, b)/(dm*(dm-1.)) : 0.);


  /* (1, - )   */
  summa = summa    +  ( b > vki[1] - 1 ?  2.*(db -  ((double)vki[1]) + 1.)*(dm - db  +   ((double)vki[1])  - 1.)*gsl_matrix_get(mr, m, b - vki[1] + 1)/( dm*(dm - 1.)) : 0.);


  /*  ( 2, - ) */
  summa = summa  +  (b >  2*vki[1] - 1 ? (db -  ((double)(vki[1] + vki[2])) + 2.)*(db -  ((double)(vki[1] + vki[2])) +  1.)*gsl_matrix_get( mr, m, b - 2*vki[1] + 2)/(dm*(dm - 1.)) :   0.);


  return ( summa ) ;

    }
  




int th_switch_number( int VK[] )
{

  /* compute switch number for 3 sim mergers */
  /* returns 3 if all mergers are distinct: k_1 > k_2  > k_3  */
  return(  (VK[1] ==  VK[2])  && (VK[2] == VK[3]) ? 0 : ( (VK[1] > VK[2]) && (VK[2] ==  VK[3]) ? 1 : ( (VK[1] == VK[2]) && (VK[2] >  VK[3]) ? 2 : 3 ) ) ) ;  

}


double th_switch_number_zero(  int b, int m, double db, double dm,  int vki[],   gsl_matrix * mr )
{

  /* k_1  == k_2 == k_3 */
  double summa = 0. ;

  /* (0, b ) */

  summa  =  summa  +   (b < m - 2 ? (1. -  db/dm)*(1. -  db/(dm - 1.))*(1. -  db/(dm - 2.))*gsl_matrix_get( mr,  m,  b) : 0.); 
  
  /* (1, - )  */
  summa = summa   +   (b > vki[1] - 1 ?  3.*(db -  (double)vki[1]   +  1.)*(dm - db +   (double)vki[1]  -  1.)*(dm - db +   (double)vki[1]   - 2.)*gsl_matrix_get(        mr, m, b - vki[1]    + 1)/(dm*(dm-1.)*(dm-2.)) : 0.) ;

  /* (2, - ) */
  summa = summa  +   (b > 2*vki[1] - 1 ?  3.*(db  -  ((double)(2*vki[1]))  + 2.)*(db  -  ((double)(2*vki[1]))  +  1.)*(dm  -  db  +  ((double)(2*vki[1]))   -   2.)*gsl_matrix_get( mr, m, b -  2*vki[1] + 2)/( dm*(dm-1.)*(dm-2.) ) : 0.); 

  /* (3, - ) */
  summa = summa  +   (b > 3*vki[1] - 1 ?   (db  -  ((double)(3*vki[1]))  + 3.)*(db  -  ((double)(3*vki[1]))  +  2.)*(db  -  ((double)(3*vki[1]))  +  1.)*gsl_matrix_get( mr, m,  b - 3*vki[1]  +  3)/( dm*(dm-1.)*(dm-2.) ) : 0.) ;


  return ( summa ) ; 
}



double th_switch_number_one(  int b, int m, double db, double dm,  int vki[],   gsl_matrix * mr )
{

  // return(  (VK[1] ==  VK[2])  && (VK[2] == VK[3]) ? 0 : ( (VK[1] > VK[2]) && (VK[2] ==  VK[3]) ? 1 : ( (VK[1] == VK[2]) && (VK[2] >  VK[3]) ? 2 : 3 ) ) ) ; 

  /* k_1 > k_2 = k_3 */

  double summa = 0. ;

  /* (0, 0, b ) */
  
  summa = summa  +   ( b < m-2 ? (1. -  db/dm)*(1. -  db/(dm - 1.))*(1. -  db/(dm - 2.))*gsl_matrix_get( mr,    m,  b) : 0.);

  /* (1, 0, - ) */
  summa = summa   +   ( b > vki[1] - 1 ?   (db -  ((double)vki[1]) + 1.)*(dm - db + ((double)vki[1]) - 1.)*(dm - db +  ((double)vki[1]) -  2.)*gsl_matrix_get(mr, m,    b - vki[1] + 1)/( dm*(dm-1.)*(dm-2.) ) : 0.);

  /*  (0, 1, - )   */
  summa = summa   +    ( b > vki[2] - 1 ? 2.*(db -  ((double)vki[2])   +  1.)*(dm  -  db  +  ((double)vki[2])    -  1.)*(dm  -  db  +  ((double)vki[2]) -  2.)*gsl_matrix_get( mr, m,  b - vki[2] + 1)/(dm*(dm-1.)*(dm-2.)) : 0.);


  /*  (0, 2, - )   */
  summa = summa   +   (b > 2*vki[2] - 1 ?  (db -  ((double)(2*vki[2]))  +  2.)*(db -  ((double)(2*vki[2]))  +  1.)*(dm  -  db  +  ((double)(2*vki[2]))  -  2.)*gsl_matrix_get( mr, m, b - 2*vki[2] + 2)/( dm*(dm-1.)*(dm-2.) ) : 0.);


  /* (1, 1, - )   */
  summa =  summa  +    (b > vki[1] + vki[2] - 1 ? 2.*(db  -  ((double)(vki[1] + vki[2])) + 2.)*(db  -  ((double)(vki[1] + vki[2])) + 1.)*(dm  -   db  +   ((double)(vki[1] + vki[2]))  -  2.)*gsl_matrix_get( mr, m,  b - vki[1] - vki[2] + 2)/( dm*(dm-1.)*(dm-2.) ) : 0.);



  /*  (1, 2, - ) */
  summa = summa  +  (b > vki[1]  +  2*vki[2]  - 1 ? (db -  ((double)(vki[1] +  2*vki[2]))  + 3.)*(db -  ((double)(vki[1] +  2*vki[2]))  + 2.)*(db -  ((double)(vki[1] +  2*vki[2]))  +  1.)*gsl_matrix_get( mr, m, b - vki[1] -  2*vki[2]   + 3)/( dm*(dm-1.)*(dm-2.) ) : 0.);


  return ( summa ) ;
}



double th_switch_number_two(  int b, int m, double db, double dm,  int vki[],   gsl_matrix * mr  )
{
 
 // return(  (VK[1] ==  VK[2])  && (VK[2] == VK[3]) ? 0 : ( (VK[1] > VK[2]) && (VK[2] ==  VK[3]) ? 1 : ( (VK[1] == VK[2]) && (VK[2] >  VK[3]) ? 2 : 3 ) ) ) ;  

  double summa = 0. ;
   
  /* (0, 0, b )  */
  summa = summa  +    (b < m-2 ? (1. -  db/dm)*(1. -  db/(dm - 1.))*(1. -  db/(dm - 2.))*gsl_matrix_get( mr, m, b) : 0.);

  
  /* (1, 0, - )   */
  summa = summa  +    ( b > vki[1] - 1 ?   2. * (db -  ((double)vki[1]) + 1.)*(dm - db + ((double)vki[1]) - 1.)*(dm - db + ((double)vki[1]) - 2.)*gsl_matrix_get(  mr, m,  b - vki[1] + 1)/( dm*(dm-1.)*(dm-2.) ) : 0.);

  /* (0, 1, - )   */
  summa = summa  +   ( b > vki[3] - 1 ?   (db -  ((double)vki[3]) + 1.)*(dm - db + ((double)vki[3]) - 1.)*(dm - db + ((double)vki[3]) - 2.)*gsl_matrix_get( mr, m,    b - vki[3] + 1)/( dm*(dm-1.)*(dm-2.) ) : 0.);

  /* (1, 1, - )   */
  summa =  summa   +  ( b > vki[1] + vki[3] - 1 ?  2.*( db -  ((double)(vki[1] + vki[3])) + 2.)*(db -  ((double)(vki[1] + vki[3])) + 1.)*(  dm - db  +  ((double)(vki[1] + vki[3])) - 2.) * gsl_matrix_get( mr, m, b - vki[1] - vki[3] + 2) /(dm*(dm-1.)*(dm-2.)) : 0.);


  /* (2, 1, - )   */
  summa = summa  +    (b >     2*vki[1]  + vki[3] - 1 ?  (db -  ((double)(2*vki[1] +  vki[3]))  + 3.)*(db -  ((double)(2*vki[1] +  vki[3]))  + 2.)*(db -  ((double)(2*vki[1] +  vki[3]))  + 1.)*gsl_matrix_get( mr, m, b - 2*vki[1]  -  vki[3] + 3)/( dm*(dm-1.)*(dm-2.) ) : 0.); 


  /* (2, 0, - )   */
  summa = summa  +   ( b > 2*vki[1] - 1 ?   (db  -  ((double)(2*vki[1])) +  2.)*(db  -  ((double)(2*vki[1]))  +  1.)*(dm  -  db  +  ((double)(2*vki[1]))  - 2.)*gsl_matrix_get( mr, m, b - 2*vki[1]  + 2)/( dm*(dm-1.)*(dm-2.) ) : 0.);
  

  return ( summa ) ; 

}



double th_switch_number_three(  int b, int m, double db, double dm,  int vki[],   gsl_matrix * mr  )
{
 
 // return(  (VK[1] ==  VK[2])  && (VK[2] == VK[3]) ? 0 : ( (VK[1] > VK[2]) && (VK[2] ==  VK[3]) ? 1 : ( (VK[1] == VK[2]) && (VK[2] >  VK[3]) ? 2 : 3 ) ) ) ;  

  /* k_1 > k_2 > k_3 */
  double summa = 0. ;
   
  /* (0, 0, 0, b )  */
  summa = summa  +    (b < m-2 ? (1. -  db/dm)*(1. -  db/(dm - 1.))*(1. -  db/(dm - 2.))*gsl_matrix_get( mr, m, b) : 0.);

  /* (1, 0, 0, - )  */
  summa = summa  +    (b > vki[1] - 1 ? ((db -  ((double)vki[1]) + 1.)*(dm - db +  ((double)vki[1])  - 1.)*(dm - db +  ((double)vki[1])  - 2.)/(dm*(dm-1.)*(dm-2.)))*gsl_matrix_get( mr, m, b - vki[1] + 1) : 0.);


   /* (0, 1, 0, - )  */
  summa = summa   +    (b > vki[2] - 1 ? ((db -  ((double)vki[2]) + 1.)*(dm - db +  ((double)vki[2])  - 1.)*(dm - db +  ((double)vki[2])  - 2.)/(dm*(dm-1.)*(dm-2.)))*gsl_matrix_get( mr, m, b - vki[2] + 1) : 0.);


  /* (0, 0, 1, - )   */
  
  summa = summa   +    (b > vki[3] - 1 ? ((db -  ((double)vki[3]) + 1.)*(dm - db +  ((double)vki[3])  - 1.)*(dm - db +  ((double)vki[3])  - 2.)/(dm*(dm-1.)*(dm-2.)))*gsl_matrix_get( mr, m, b - vki[3] + 1) : 0.);


  /* (1, 1, 0, - )    */
  summa = summa    +   (b > vki[1] + vki[2] - 1 ?  (db -   ((double)(vki[1] + vki[2])) + 2.)*    (db -  ((double)(vki[1] + vki[2]))  +1.)*(dm  -  db  + ((double)(vki[1] + vki[2])) - 2.) *   gsl_matrix_get( mr, m, b - vki[1] - vki[2] + 2)/( dm*(dm-1.)*(dm-2.)) :   0.);


  /*  (1, 0, 1, - ) */
  summa = summa    +   (b > vki[1] + vki[3] - 1 ?  (db -   ((double)(vki[1] + vki[3])) + 2.)*    (db -  ((double)(vki[1] + vki[3]))  +1.)*(dm  -  db  + ((double)(vki[1] + vki[3])) - 2.) *   gsl_matrix_get( mr, m, b - vki[1] - vki[3] + 2)/( dm*(dm-1.)*(dm-2.)) :   0.);



  /*  (0, 1, 1, - )   */
  summa = summa    +    (b > vki[3] + vki[2] - 1 ?  (db -   ((double)(vki[3] + vki[2])) + 2.)*    (db -  ((double)(vki[3] + vki[2]))  +1.)*(dm  -  db  + ((double)(vki[3] + vki[2])) - 2.) *   gsl_matrix_get( mr, m, b - vki[3] - vki[2] + 2)/( dm*(dm-1.)*(dm-2.)) :   0.);



  /* (1, 1, 1, - ) */
  summa = summa    +   (b > vki[1] + vki[2] + vki[3] - 1 ?  ( db -  ((double)(vki[1] + vki[2] + vki[3])) + 3.)*(db -  ((double)(vki[1] + vki[2] + vki[3]))  + 2.)*(db  -  ((double)(vki[1] + vki[2] + vki[3]))  + 1.)*gsl_matrix_get( mr, m,  b -  (vki[1] + vki[2] + vki[3]) + 3)/( dm*(dm-1.)*(dm-2.)) : 0.);

  return ( summa ) ;

}






int fm_switch_number( int VK[] )
{

  // returns switch number for 4 sim mergers
  // last number 7 is for k_1 = k_2 = k_3 = k_4; 
  // ie all 4 mergers of same size
  
  return( ((VK[1] > VK[2]) && (VK[2] > VK[3])) && (VK[3] > VK[4]) ? 0 : ( (VK[1] == VK[2]) && ((VK[2] > VK[3]) && (VK[3] > VK[4])) ? 1 : ( (VK[1] > VK[2]) &&      ((VK[2]  ==  VK[3]) && (VK[3] > VK[4])) ? 2 : ( (VK[1] > VK[2]) &&  ((VK[2] > VK[3]) && (VK[3]  ==  VK[4])) ? 3 : ( (VK[1]  ==  VK[2]) && ((VK[2] > VK[3])       && (VK[3]  ==  VK[4])) ? 4 : ( (VK[1] ==  VK[2]) && ((VK[2]  ==  VK[3]) && (VK[3] > VK[4])) ? 5 : ( (VK[1] > VK[2]) && ((VK[2]  ==  VK[3])  && (VK[3]  ==  VK[4])) ? 6 : 7   ) ) ) ) ) ) ) ;
}


double fm_switch_number_zero( int b, int m, double db, double dm,  int vki[],   gsl_matrix * mr )
{

  // returns addition when all mergers distinct
  // ie when  k_1  > k_2 > k_3 > k_4

  double summa  = 0 ; 

  /* (0, 0, 0, 0, b ) */
  summa  =  summa  +   (b < m - 3 ? (dm - db)*(dm - db - 1.)*(dm - db - 2.)*(dm - db - 3.)*gsl_matrix_get( mr,  m, b)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;

  /*  (1, 0, 0, 0, - ) */
  summa  =  summa  +    (b > vki[1] - 1 ?  (db  -  ((double)vki[1]) + 1.)*(dm  -   db  +  ((double)vki[1]) -  1.)*(dm  -   db  +  ((double)vki[1]) -  2.)*(dm  -   db  +  ((double)vki[1]) -  3.)*gsl_matrix_get( mr, m, b - vki[1] + 1)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;


  /*  (0, 1, 0, 0, - ) */
  summa  =  summa  +         ( b > vki[2] - 1 ?  (db  -  ((double)vki[2]) + 1.)*(dm  -   db  +  ((double)vki[2]) -  1.)*(dm  -   db  +  ((double)vki[2]) -  2.)*(dm  -   db  +  ((double)vki[2]) -  3.)*gsl_matrix_get( mr, m, b - vki[2] + 1)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;


  /* (0, 0, 1, 0, - )  */
  summa = summa  +      ( b > vki[3] - 1 ?  (db  -  ((double)vki[3]) + 1.)*(dm  -   db  +  ((double)vki[3]) -  1.)*(dm  -   db  +  ((double)vki[3]) -  2.)*(dm  -   db  +  ((double)vki[3]) -  3.)*gsl_matrix_get( mr, m, b - vki[3]  +  1)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.)   ;

  /* (0, 0, 0, 1, - )  */

  summa  =  summa    +       ( b > vki[4]  -  1 ?  (db  -  ((double)vki[4]) + 1.)*(dm  -   db  +  ((double)vki[4]) -  1.)*(dm  -   db  +  ((double)vki[4]) -  2.)*(dm  -   db  +  ((double)vki[4]) -  3.)*gsl_matrix_get( mr, m, b - vki[4]  +  1)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.)  ;


  /*  (1, 1, 0, 0, - )  */
  summa = summa   +   (b > vki[1] + vki[2] - 1 ?   (db -  ((double)(vki[1] + vki[2]))     + 2.)*(db -  ((double)(vki[1] + vki[2]))     +  1.)*(dm  -  db  +  ((double)(vki[1] + vki[2])) - 2.)*(dm  -  db  +  ((double)(vki[1] + vki[2])) -  3.)*gsl_matrix_get( mr, m,  b - vki[1] - vki[2] + 2)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;

  /* (1, 0, 1, 0, - )  */
  
  summa  =  summa  +     (b > vki[1] + vki[3] - 1 ?   (db -  ((double)(vki[1] + vki[3]))     + 2.)*(db -  ((double)(vki[1] + vki[3]))     +  1.)*(dm  -  db  +  ((double)(vki[1] + vki[3])) - 2.)*(dm  -  db  +  ((double)(vki[1] + vki[3])) -  3.)*gsl_matrix_get(  mr, m,  b - vki[1] - vki[3] + 2)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.)  ;

  /* (1, 0, 0, 1, - )   */

  summa  =  summa  +    (b > vki[1] + vki[4] - 1 ?   (db -  ((double)(vki[1] + vki[4]))     + 2.)*(db -  ((double)(vki[1] + vki[4]))     +  1.)*(dm  -  db  +  ((double)(vki[1] + vki[4])) - 2.)*(dm  -  db  +  ((double)(vki[1] + vki[4])) -  3.)*gsl_matrix_get(mr, m,  b - vki[1] - vki[4] + 2)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) :    0. )  ;
  
  
  /* (0, 1, 1, 0, - )  */

  summa  =  summa   +     (b > vki[3] + vki[2] - 1 ?   (db -  ((double)(vki[3] + vki[2]))     + 2.)*(db -  ((double)(vki[3] + vki[2]))     +  1.)*(dm  -  db  +  ((double)(vki[3] + vki[2])) - 2.)*(dm  -  db  +  ((double)(vki[3] + vki[2])) -  3.)*gsl_matrix_get(  mr, m,  b - vki[3] - vki[2] + 2)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;


  /*  (0, 1, 0, 1, - )   */
  summa =  summa   +    (b > vki[4] + vki[2] - 1 ?   (db -  ((double)(vki[4] + vki[2]))     + 2.)*(db -  ((double)(vki[4]  + vki[2]))     +  1.)*(dm  -  db  +  ((double)(vki[4]  +  vki[2])) - 2.)*(dm  -  db  +  ((double)(vki[4]  +  vki[2])) -  3.)*gsl_matrix_get(  mr, m,  b - vki[4] - vki[2] + 2)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;

  /*  (0, 0, 1, 1, - )  */
  
  summa  =  summa   +    (b > vki[3] + vki[4] - 1 ?   (db -  ((double)(vki[3] + vki[4]))     + 2.)*(db -  ((double)(vki[3] + vki[4]))     +  1.)*(dm  -  db  +  ((double)(vki[3] + vki[4])) - 2.)*(dm  -  db  +  ((double)(vki[3] + vki[4])) -  3.)*gsl_matrix_get(  mr, m,  b - vki[3] - vki[4] + 2)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.)  ;

  /* (0,  1, 1, 1, - )   */
  summa  =  summa   +   (b > vki[2] + vki[3] + vki[4] - 1 ?   (db  -  ((double)(vki[2] +  vki[3] + vki[4]))  + 3.)*(db  -  ((double)(vki[2] +  vki[3] + vki[4]))  + 2.)*(db  -  ((double)(vki[2] +  vki[3] + vki[4]))  + 1.)*(dm  -  db  +  ((double)(vki[2] +  vki[3] + vki[4]))  - 3.)*gsl_matrix_get( mr,  m,  b - vki[2] - vki[3] - vki[4] + 3)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;


  /* (1, 0, 1, 1, - )  */
  summa =  summa   +     (b > vki[1] + vki[3] + vki[4] - 1 ?   (db  -  ((double)(vki[1] +  vki[3] + vki[4]))  + 3.)*(db  -  ((double)(vki[1] +  vki[3] + vki[4]))  + 2.)*(db  -  ((double)(vki[1] +  vki[3] + vki[4]))  + 1.)*(dm  -  db  +  ((double)(vki[1] +  vki[3] + vki[4]))  - 3.)*gsl_matrix_get( mr, m, b - vki[1] - vki[3] - vki[4] + 3)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;
  
  /* (1, 1, 0, 1, - )  */
  summa  =  summa  +    (b > vki[2] + vki[1] + vki[4] - 1 ?   (db  -  ((double)(vki[2] +  vki[1] + vki[4]))  + 3.)*(db  -  ((double)(vki[2] +  vki[1] + vki[4]))  + 2.)*(db  -  ((double)(vki[2] +  vki[1] + vki[4]))  + 1.)*(dm  -  db  +  ((double)(vki[2] +  vki[1] + vki[4]))  - 3.)*gsl_matrix_get( mr, m, b - vki[2] - vki[1] - vki[4] + 3)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;

  
  /* (1, 1, 1, 0, - )  */
  summa =  summa         +    (b > vki[2] + vki[3] + vki[1] - 1 ?   (db  -  ((double)(vki[2] +  vki[3] + vki[1]))  + 3.)*(db  -  ((double)(vki[2] +  vki[3] + vki[1]))  + 2.)*(db  -  ((double)(vki[2] +  vki[3] + vki[1]))  + 1.)*(dm  -  db  +  ((double)(vki[2] +  vki[3] + vki[1]))  - 3.)*gsl_matrix_get( mr, m, b - vki[2] - vki[3] - vki[1] + 3)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;

  /* (1, 1, 1, 1, - ) */
  summa  =  summa  +    (b > vki[1] + vki[2] + vki[3] + vki[4] - 1 ?  (db -  ((double)(vki[1] + vki[2] + vki[3] + vki[4])) + 4.)*(db -  ((double)(vki[1] + vki[2] + vki[3] + vki[4])) +  3.)*(db -  ((double)(vki[1] + vki[2] + vki[3] + vki[4])) + 2.)*(db -  ((double)(vki[1] + vki[2] + vki[3] + vki[4]))  +  1.)*gsl_matrix_get( mr, m, b - (vki[1] + vki[2] + vki[3] + vki[4]) + 4)/( dm*(dm- 1.)*(dm - 2.)*(dm - 3.) ) : 0.); 


  return( summa ) ;
}


double fm_switch_number_one(  int b, int m, double db, double dm,  int  vki[],   gsl_matrix * mr )
{

  //// return 4 sim switch number 1
  //// when  k_1 = k_2 and k_2 > k_3  > k_4

  double summa = 0. ;
  
  /* (0, 0, 0, - )  */
  summa =  summa +    (b < m - 3 ? (dm - db)*(dm - db - 1.)*(dm - db - 2.)*(dm - db - 3.)*gsl_matrix_get( mr, m, b)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;

  /* (1, 0, 0, - )  */
  summa  =  summa  +   (b > vki[1] - 1 ?  2.*(db -  ((double)vki[1])  +  1.)*(dm  -  db  +  ((double)vki[1])  -  1.)*(dm  -  db  +  ((double)vki[1])  -  2.)*(dm  -  db  +  ((double)vki[1])  -  3.)*gsl_matrix_get( mr,  m,  b - vki[1] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;

  /* (0, 1, 0, - )  */
  summa =  summa  +   (b > vki[3] - 1 ?  (db -  ((double)vki[3])  +  1.)*(dm  -  db  +  ((double)vki[3])  -  1.)*(dm  -  db  +  ((double)vki[3])  -  2.)*(dm  -  db  +  ((double)vki[3])  -  3.)*gsl_matrix_get( mr, m,  b - vki[3] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ; 


  /* ( 0, 0, 1, - )  */
  summa =  summa  +       (b > vki[4] - 1 ?   (db -  ((double)vki[4])  +  1.)*(dm  -  db  +  ((double)vki[4])  -  1.)*(dm  -  db  +  ((double)vki[4])  -  2.)*(dm  -  db  +  ((double)vki[4])  -  3.)*gsl_matrix_get( mr, m,  b  - vki[4] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ; 

  /* (2, 0, 0, - ) */

  summa  =  summa    +    (b >  2*vki[1] - 1 ?  (db -   ((double)(vki[1] + vki[2])) + 2.)*(db -   ((double)(vki[1] + vki[2])) + 1.)*(dm  -  db  +   ((double)(vki[1] + vki[2])) -  2.)*(dm  -  db  +  ((double)(vki[1] + vki[2])) - 3.)*gsl_matrix_get(mr, m, b  -  2*vki[1]  + 2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;


  /* (1, 1, 0, - ) */

  summa  =  summa  +  (b > vki[1]  + vki[3] - 1 ?   2.*(db -  ((double)(vki[1] + vki[3]))  + 2.)*(db -  ((double)(vki[1] + vki[3]))  + 1.)*(dm  -  db   +   ((double)(vki[1] + vki[3]))      -  2.)*(dm  -  db  +  ((double)(vki[1] + vki[3]))  -  3.)*gsl_matrix_get( mr, m, b - vki[1] - vki[3] + 2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.)  ;


  /* (1, 0, 1, - ) */

  summa  =  summa  +  (b > vki[1]  + vki[4] - 1 ?  2.*(db -  ((double)(vki[1] + vki[4]))  + 2.)*(db -  ((double)(vki[1] + vki[4]))  + 1.)*(dm  -  db   +   ((double)(vki[1] + vki[4]))      -  2.)*(dm  -  db  +  ((double)(vki[1] + vki[4]))  -  3.)*gsl_matrix_get( mr, m, b - vki[1] - vki[4] + 2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.)  ;




  /* (0, 1, 1, - ) */


  summa = summa  +  (b > vki[4]  + vki[3] - 1 ?  (db -  ((double)(vki[4] + vki[3]))  + 2.)*(db -  ((double)(vki[4] + vki[3]))  + 1.)*(dm  -  db   +   ((double)(vki[4] + vki[3])) -  2.)*(dm  -  db  +  ((double)(vki[4] + vki[3]))  -  3.)*gsl_matrix_get( mr, m, b - vki[4] - vki[3] + 2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);

  /* (2, 1, 0, - ) */

  summa  =  summa  +  (b > 2*vki[1] + vki[3] - 1 ?   (db -  ((double)(2*vki[1] + vki[3]))  + 3.)*(db -  ((double)(2*vki[1] + vki[3]))  + 2.)*(db -  ((double)(2*vki[1] + vki[3]))  +  1.)*(dm  -  db +   ((double)(2*vki[1] + vki[3])) - 3.)*gsl_matrix_get( mr, m, b - 2*vki[1] -  vki[3] + 3)/(dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);


  /* (2, 0, 1, - ) */
 
  summa =  summa  +  (b > 2*vki[1] + vki[4] - 1 ?   (db -  ((double)(2*vki[1] + vki[4]))  + 3.)*(db -  ((double)(2*vki[1] + vki[4]))  + 2.)*(db -  ((double)(2*vki[1] + vki[4]))  +   1.)*(dm  -  db +   ((double)(2*vki[1] + vki[4]))  -  3.)*gsl_matrix_get( mr, m, b - 2*vki[1]  - vki[4] + 3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);  



  /* (1, 1, 1, - ) */

  summa =  summa   +    (b >  vki[1] + vki[3]  + vki[4] - 1 ? 2.*(db  -  ((double)(vki[1] + vki[3] + vki[4]))  +  3.)*(db  -  ((double)(vki[1] + vki[3] + vki[4]))  +  2.)*(db  -  ((double)(vki[1] + vki[3] + vki[4]))  +  1.)*(dm  -   db  +  ((double)(vki[1] + vki[3] + vki[4]))  -  3.)*gsl_matrix_get( mr, m,  b -  vki[1]  - vki[3]  - vki[4] +  3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.)  ;



   /* (2, 1, 1, - ) */

  summa = summa  +   (b > 2*vki[1]  + vki[3] + vki[4] - 1 ?  (db  -  ((double)(2*vki[1] + vki[3] + vki[4]))  + 4.)*(db  -  ((double)(2*vki[1] + vki[3] + vki[4]))  + 3.)*(db  -  ((double)(2*vki[1] + vki[3] + vki[4]))  + 2.)*(db  -  ((double)(2*vki[1] + vki[3] + vki[4]))  +  1.)*gsl_matrix_get(mr, m, b - 2*vki[1] - vki[3] - vki[4] + 4)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;


  return( summa ) ; 
}


double fm_switch_number_two(  int b, int m, double db, double dm,  int  vki[],   gsl_matrix * mr )
{

  // returns  summation when 
  // k_1 > k_2 = k_3 > k_4

  double summa = 0. ; 

  
  /* (0, 0, 0, b )  */
  summa = summa  +   (b < m - 3 ?   (dm - db)*(dm - db - 1.)*(dm - db - 2.)*(dm - db - 3.)*gsl_matrix_get( mr, m, b)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);

  
  /* (1, 0, 0, - )  */

  summa  =  summa   +   (b > vki[1] - 1 ?  (db  -  ((double)vki[1])  + 1.)*(dm  -  db  +  ((double)vki[1])  - 1.)*(dm  -  db  +  ((double)vki[1])  - 2.)*(dm  -  db  +  ((double)vki[1])  - 3.)*gsl_matrix_get( mr, m, b - vki[1] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.)  ;

  /* (0, 1, 0, - ) */
  
  summa  =  summa    +    (b > vki[2] - 1 ?   2.*(db  -  ((double)vki[2])  + 1.)*(dm  -  db  +  ((double)vki[2])  - 1.)*(dm  -  db  +  ((double)vki[2])  - 2.)*(dm  -  db  +  ((double)vki[2])  - 3.)*gsl_matrix_get( mr,  m,  b - vki[2] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;        


  /* (0, 0, 1, - )  */

  summa =  summa   +  (b > vki[4] - 1 ?  (db  -  ((double)vki[4])  + 1.)*(dm  -  db  +  ((double)vki[4])  - 1.)*(dm  -  db  +  ((double)vki[4])  - 2.)*(dm  -  db  +  ((double)vki[4])  - 3.)*gsl_matrix_get( mr,  m,   b - vki[4] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;


  /* (0, 2, 0, -  )  */

  summa = summa   +      (b >  2*vki[2]  -  1 ?  (db  -   ((double)(2*vki[2]))  + 2.)*(db  -   ((double)(2*vki[2]))  +  1.)*(dm  -  db   +   ((double)(2*vki[2])) -  2.)*(dm  -  db   +   ((double)(2*vki[2])) -  3.)*gsl_matrix_get( mr, m, b - 2*vki[2] + 2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;

  /* (1, 1, 0, - )  */

  summa  =  summa    +   (b > vki[1] + vki[2] - 1 ?  2.*(db  -  ((double)(vki[1] + vki[2])) + 2.)*(db  -  ((double)(vki[1] + vki[2])) + 1.)*(dm  -  db  +   ((double)(vki[1] + vki[2]))  -  2.)*(dm  -  db  +  ((double)(vki[1] + vki[2])) -  3.)*gsl_matrix_get( mr, m,  b - vki[1]  -  vki[2]  +  2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);

   /* (1, 0, 1, - )  */
  summa  =  summa  +   (b > vki[1] + vki[4] - 1 ?  (db  -  ((double)(vki[1] + vki[4])) + 2.)*(db  -  ((double)(vki[1] + vki[4])) + 1.)*(dm  -  db  +   ((double)(vki[1] +  vki[4]))  -  2.)*(dm  -  db  +  ((double)(vki[1]  +  vki[4])) -  3.)*gsl_matrix_get( mr,  m, b - vki[1] - vki[4]  + 2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);


  /*  (0, 1, 1, - )  */
  
  summa = summa  +  (b > vki[2]  +  vki[4] - 1 ?  2.*(db  -  ((double)(vki[2] + vki[4])) + 2.)*(db  -  ((double)(vki[2] + vki[4])) + 1.)*(dm  -  db  +   ((double)(vki[2] +  vki[4]))  -  2.)*(dm  -  db  +  ((double)(vki[2]  +  vki[4])) -  3.)*gsl_matrix_get( mr,  m, b - vki[2] - vki[4]  + 2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);

   /* (1, 1, 1, - )  */

  summa =  summa   +       (b > vki[1] + vki[2]  +  vki[4] - 1 ?  2.*(db  -  ((double)(vki[1] +  vki[2]  +  vki[4])) + 3.)*(db  -  ((double)(vki[1] + vki[2] +  vki[4])) + 2.)*(db  -  ((double)(vki[1] + vki[2] +  vki[4])) + 1.)*(dm  -  db  +   ((double)(vki[1] +  vki[2] +  vki[4]))  -  3.)*gsl_matrix_get( mr, m,  b - vki[1]  - vki[2] -   vki[4]  +  3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);

  /* (1, 2, 0, - )  */

  summa =  summa   +   (b > vki[1]  +  2*vki[2] - 1 ?   (db -   ((double)(vki[1] +   2*vki[2]))  + 3.)*(db -   ((double)(vki[1] +   2*vki[2]))  +  2.)*(db -   ((double)(vki[1] +   2*vki[2]))  +  1.)*(dm  -   db   +   ((double)(vki[1] +   2*vki[2]))  -  3.)*gsl_matrix_get( mr, m,  b - vki[1]  - 2*vki[2]  +  3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0. );

   /* (0, 2, 1, - )  */

  summa = summa   +   (b > vki[4]  +  2*vki[2] - 1 ?   (db -   ((double)(vki[4] +   2*vki[2]))  + 3.)*(db -   ((double)(vki[4] +   2*vki[2]))  +  2.)*(db -   ((double)(vki[4] +   2*vki[2]))  +  1.)*(dm  -   db   +   ((double)(vki[4] +   2*vki[2]))  -  3.)*gsl_matrix_get( mr, m,  b - vki[4]  -  2*vki[2]  +  3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);

  /* (1,2,1, - )  */

  summa = summa  +  (b >  vki[1] +  2*vki[2]  +  vki[4] - 1 ?   (db -   ((double)(vki[4] +   2*vki[2]   +   vki[1]))  +  4.)*(db -   ((double)(vki[4] +   2*vki[2]   +   vki[1]))  +  3.)*(db -   ((double)(vki[4] +   2*vki[2]   +   vki[1]))  +  2.)*(db -   ((double)(vki[4] +   2*vki[2]   +   vki[1]))  +  1.)*gsl_matrix_get(mr, m,    b - vki[1] - vki[2] - vki[3] - vki[4] + 4)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);

  return( summa ) ;

}



double fm_switch_number_three(  int b, int m, double db, double dm,  int  vki[],   gsl_matrix * mr )
{

  //  k_1 > k_2 >  k_3 == k_4
  double summa = 0. ;

  //gsl_matrix_set( Mnb, n, b, gsl_matrix_get( Mnb, n, b)  + 
			  /* (0, 0, 0, b ) */
  summa =  summa  +    (b < m - 3 ?  (dm - db)*(dm - db - 1.)*(dm - db - 2.)*(dm - db - 3.)*gsl_matrix_get(  mr, m, b)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;


			  /*  (1, 0, 0, - )   */
  summa =  summa  +  (b > vki[1] - 1 ?  (db  -  ((double)vki[1])  + 1.)*(dm  -  db  +  ((double)vki[1])  - 1.)*(dm  -  db  +  ((double)vki[1])  - 2.)*(dm  -  db  +  ((double)vki[1])  - 3.)*gsl_matrix_get( mr, m, b - vki[1] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ; 


			  /*  (0, 1, 0, - )  */
  summa  =  summa +   (b > vki[2] - 1 ?   (db  -  ((double)vki[2])  + 1.)*(dm  -  db  +  ((double)vki[2])  - 1.)*(dm  -  db  +  ((double)vki[2])  - 2.)*(dm  -  db  +  ((double)vki[2])  - 3.)*gsl_matrix_get( mr, m, b - vki[2] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;

  
  /* (0, 0, 1, - ) */
  summa  = summa  +   (b > vki[3] - 1 ?   2.*(db  -  ((double)vki[3])  + 1.)*(dm  -  db  +  ((double)vki[3])  - 1.)*(dm  -  db  +  ((double)vki[3])  - 2.)*(dm  -  db  +  ((double)vki[3])  - 3.)*gsl_matrix_get( mr, m, b - vki[3] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;



  /* (1, 1, 0, - ) */
  summa  =  summa  +  (b > vki[1] + vki[2] - 1 ?   (db  -  ((double)(vki[1] + vki[2])) + 2.)*(db  -  ((double)(vki[1] + vki[2])) + 1.)*(dm  -  db  +   ((double)(vki[1] + vki[2]))  -  2.)*(dm  -  db  +  ((double)(vki[1] + vki[2])) -  3.)*gsl_matrix_get( mr, m,  b - vki[1]  -  vki[2]  +  2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.)) : 0.);


  /* (1, 0, 1, - )  */
  summa =  summa  +  (b > vki[1] + vki[3] - 1 ?  2.*(db  -  ((double)(vki[1] + vki[3])) + 2.)*(db  -  ((double)(vki[1] + vki[3])) + 1.)*(dm  -  db  +  ((double)(vki[1] + vki[3]))  -  2.)*(dm  -  db  +  ((double)(vki[1] + vki[3])) -  3.)*gsl_matrix_get( mr, m,  b - vki[1]  -  vki[3]  +  2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.)) : 0.);




  /* (0, 1, 1, - ) */
  summa  = summa  +  (b > vki[3] + vki[2] - 1 ?  2.*(db  -  ((double)(vki[3] + vki[2])) + 2.)*(db  - ((double)(vki[3] + vki[2])) + 1.)*(dm  -  db  +   ((double)(vki[3] + vki[2]))  -  2.)*(dm  -  db  +  ((double)(vki[3] + vki[2])) -  3.)*gsl_matrix_get( mr, m,  b - vki[3]  -  vki[2]  +  2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.)) : 0.);  


  /* (0, 0, 2, - )  */
  summa =  summa  +    (b > vki[3] + vki[3] - 1 ?  (db  -  ((double)(vki[3] + vki[3])) + 2.)*(db  -  ((double)(vki[3] + vki[3])) + 1.)*(dm  -  db  +   ((double)(vki[3] + vki[3]))  -  2.)*(dm  -  db  +  ((double)(vki[3] + vki[3])) -  3.)*gsl_matrix_get( mr, m,  b - vki[3]  -  vki[3]  +  2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.)) : 0.);



  /* (1, 1, 1, - ) */
  summa =  summa  +  (b > vki[1] + vki[2] + vki[3]  -  1 ?  2.*(db  -  ((double)(vki[1] + vki[2]  +  vki[3])) +  3.)*(db  -  ((double)(vki[1] + vki[2] +   vki[3])) + 2.)*(db  -  ((double)(vki[1] + vki[2]  +  vki[3])) +  1.) *(dm  -  db  +   ((double)(vki[1] + vki[2] +  vki[3]))  -  3.)*gsl_matrix_get( mr, m,  b - vki[1] - vki[2]  -  vki[3]  +  3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ; 



  /* (1, 0, 2, - )  */
  summa = summa +   (b > vki[1] + vki[3] + vki[3]  -  1 ?  (db  -  ((double)(vki[1] + vki[3]  +  vki[3])) +  3.)*(db  -  ((double)(vki[1] + vki[3] +   vki[3])) + 2.)*(db  -  ((double)(vki[1] + vki[3]  +  vki[3])) +  1.) *(dm  -  db  +   ((double)(vki[1] + vki[3] +  vki[3]))  -  3.)*gsl_matrix_get( mr, m,  b - vki[1] - vki[3]  -  vki[3]  +  3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;

  
  /* (0, 1, 2, - )  */
  summa = summa  +   (b > vki[2] + vki[3] + vki[3]  -  1 ?  (db  -  ((double)(vki[2] + vki[3]  +  vki[3])) +  3.)*(db  -  ((double)(vki[2] + vki[3] +   vki[3])) + 2.)*(db  -  ((double)(vki[2] + vki[3]  +  vki[3])) +  1.) *(dm  -  db  +   ((double)(vki[2] + vki[3] +  vki[3]))  -  3.)*gsl_matrix_get( mr, m,  b - vki[2] - vki[3]  -  vki[3]  +  3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;


  /* (1, 1, 2, - )  */
  summa =  summa  +   (b > vki[1] + vki[2] +  2*vki[3]  -  1 ?  (db  -  ((double)(vki[1] + vki[2]  +   2*vki[3])) +  4.)*(db  -  ((double)(vki[1] +  vki[2] +   2*vki[3])) + 3.)*(db  -  ((double)(vki[1] + vki[2]  +  2*vki[3])) +  2.) *(db  -  ((double)(vki[1] +  vki[2] +   2*vki[3])) + 1.) *gsl_matrix_get( mr, m,  b - vki[1] - vki[2]  -   2*vki[3]  +  4)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ; 
  

  return( summa ) ;
}


double fm_switch_number_four( int b,  int m,   double db, double dm,    int VK[], gsl_matrix * mr  )
{

  // returns  addition when  k_1 = k_2 > k3 = k_4 
  // remains to multiply with  jump probability and  prob of hitting level k
  double summa = 0. ;
  
  /* (0, 0, b ) */
  summa =  summa +  ( b < m - 3 ? (dm - db)*(dm - db - 1.)*(dm - db - 2.)*(dm - db - 3.)*gsl_matrix_get(mr, m, b)/( dm*(dm - 1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;
  /* (1, 0, - ) */
  summa =  summa  +  (b > VK[1] - 1 ?  2.*(db  -  ((double)VK[1])  +  1.)*(dm  -  db  +  ((double)VK[1])  - 1.)*(dm  -  db  +  ((double)VK[1])  - 2.)*(dm  -  db  +  ((double)VK[1])  - 3.)*gsl_matrix_get( mr, m, b -  VK[1] +  1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;
  /* (0, 1, - ) */
  summa = summa  +   (b > VK[3] - 1 ?  2.*(db  -  ((double)VK[3])  +  1.)*(dm  -  db  +  ((double)VK[3])  - 1.)*(dm  -  db  +  ((double)VK[3])  - 2.)*(dm  -  db  +  ((double)VK[3])  - 3.)*gsl_matrix_get( mr, m, b -  VK[3] +  1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;
  /* (1, 1, - ) */
  summa  =  summa  +   (b > VK[1]  +  VK[3] - 1 ?  4.*(db  -  ((double)(VK[1]  +  VK[3]))   +   2.)*(db  -  ((double)(VK[1]  +  VK[3]))   +   1.)*(dm  -  db  +   ((double)(VK[1]  +  VK[3]))   - 2.)*(dm  -  db  +  ((double)VK[3] + VK[1] )  - 3.)*gsl_matrix_get( mr, m, b - VK[1] -  VK[3] +  2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;
  /*  (2, 0, - ) */
  summa  =  summa  +  (b > VK[1]  +  VK[1] - 1 ?  (db  -  ((double)(VK[1]  +  VK[1]))   +   2.)*(db  -  ((double)(VK[1]  +  VK[1]))   +   1.)*(dm  -  db  +   ((double)(VK[1]  +  VK[1]))   - 2.)*(dm  -  db  +  ((double)VK[1] + VK[1] )  - 3.)*gsl_matrix_get( mr, m, b - VK[1] -  VK[1] +  2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;
  /*  (0, 2, - ) */
  summa  =  summa  +  (b > VK[3]  +  VK[3] - 1 ?  (db  -  ((double)(VK[3]  +  VK[3]))   +   2.)*(db  -  ((double)(VK[3]  +  VK[3]))   +   1.)*(dm  -  db  +   ((double)(VK[3]  +  VK[3]))   - 2.)*(dm  -  db  +  ((double)VK[3] + VK[3] )  - 3.)*gsl_matrix_get( mr, m, b - VK[3] -  VK[3] +  2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;
  /* (2, 1, - ) */
  summa  =  summa  +   (b >  2*VK[1]  +   VK[3] - 1 ?   2.*(db  -  ((double)(2*VK[1]  +  VK[3]))   +   3.)*(db  -  ((double)(2*VK[1]  +  VK[3]))   +   2.)*(db  -  ((double)(2*VK[1]  +  VK[3]))   +   1.)*(dm  -  db  +   ((double)(2*VK[1]  +  VK[3]))   -  3.)*gsl_matrix_get( mr, m, b -  2*VK[1] -  VK[3] +  3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;
  /*  (1, 2, - )  */
  summa  =   summa  +   (b >  VK[1]  +    2*VK[3] - 1 ?   2.*(db  -  ((double)(VK[1]  +   2*VK[3]))   +   3.)*(db  -  ((double)(VK[1]  +   2*VK[3]))   +   2.)*(db  -  ((double)(VK[1]  +   2*VK[3]))   +   1.)*(dm  -  db  +   ((double)(VK[1]  +   2*VK[3]))   -  3.)*gsl_matrix_get( mr, m, b -  VK[1] -   2*VK[3] +   3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;
  /*  (2, 2, - ) */
  summa  =   summa  +   (b >  2*VK[1]  +    2*VK[3] - 1 ?   (db  -  ((double)(2*VK[1]  +   2*VK[3]))   +   4.)*(db  -  ((double)(2*VK[1]   +   2*VK[3]))   +   3.)*(db     -   ((double)(2*VK[1]  +   2*VK[3]))   +  2.)*(db  -  ((double)(2*VK[1]  +  2*VK[3]))  +  1.)*gsl_matrix_get( mr, m, b -  2*VK[1] -   2*VK[3] +  4)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;

  return( summa ) ; 
}



double fm_switch_number_five(  int b, int m, double db, double dm,  int  vki[],   gsl_matrix * mr )
{

  // k_1 == k_2 == k_3 > k_4

  double summa = 0. ;

  /* (0, 0, 0, b ) */
  summa = summa +  (b < m - 3 ?  (dm - db)*(dm - db - 1.)*(dm - db - 2.)*(dm - db - 3.)*gsl_matrix_get(mr, m, b)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;


     /**  (1,0, - )   */

  summa  =  summa  +   (b > vki[1] - 1 ?  3.*(db  -  ((double)vki[1])  + 1.)*(dm  -  db  +  ((double)vki[1])  - 1.)*(dm  -  db  +  ((double)vki[1])  - 2.)*(dm  -  db  +  ((double)vki[1])  - 3.)*gsl_matrix_get( mr, m, b - vki[1] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;

  
  /* (0, 1, - ) */
  summa = summa  +   (b > vki[4] - 1 ?  (db  -  ((double)vki[4])  + 1.)*(dm  -  db  +  ((double)vki[4])  -  1.)*(dm  -  db  +  ((double)vki[4])  - 2.)*(dm  -  db  +  ((double)vki[4])  - 3.)*gsl_matrix_get( mr, m, b - vki[4] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);


  /* (2, 0, - )   */
  summa = summa  +  (b > vki[1] + vki[2] - 1 ?   3.*(db  -  ((double)(vki[1] + vki[2])) + 2.)*(db  -  ((double)(vki[1] + vki[2])) + 1.)*(dm  -  db  +  ((double)(vki[1] + vki[2])) -  2.)*(dm  -  db  +  ((double)(vki[1] + vki[2])) -  3.)*gsl_matrix_get( mr, m,  b - vki[1]  -  vki[2]  +  2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.)) : 0.);



  /*  (1, 1, - )  */
  summa = summa +  (b > vki[1] + vki[4] - 1 ?   3.*(db  -  ((double)(vki[1] + vki[4])) + 2.)*(db  -  ((double)(vki[1] + vki[4])) + 1.)*(dm  -  db  +   ((double)(vki[1] + vki[4]))  -  2.)*(dm  -  db  +  ((double)(vki[1] + vki[4])) -  3.)*gsl_matrix_get( mr, m,  b - vki[1]  -  vki[4]  +  2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);


    
  /*  (3, 0, - )   */
  summa = summa  +   (b > 3*vki[1] - 1 ?  (db -  ((double)(3*vki[1])) + 3.)*(db -  ((double)(3*vki[1])) +  2.)*(db -  ((double)(3*vki[1])) + 1.)*(dm - db  +   ((double)(3*vki[1]))  - 3.)*gsl_matrix_get( mr, m, b -  3*vki[1]  + 3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.); 



   /*  (2, 1, - ) */
  summa = summa  +  (b > 2*vki[1]  + vki[4]  - 1 ?  3.*(db -  ((double)(2*vki[1]  +   vki[4])) + 3.)*(db -  ((double)(2*vki[1]  +   vki[4])) +  2.)*(db -  ((double)(2*vki[1]  +    vki[4])) + 1.)*(dm - db  +   ((double)(2*vki[1]  +   vki[4]))  - 3.)*gsl_matrix_get( mr, m, b -  2*vki[1]   -   vki[4]   + 3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);


   /*  (3, 1, - )  */
  summa = summa  +    ( b >  3*vki[1]  +  vki[4] - 1 ?   (db -  ((double)(3*vki[1]  +  vki[4]))  +  4. )*(db -  ((double)(3*vki[1]  +   vki[4])) +  3.)*(db -  ((double)(3*vki[1]  +   vki[4])) + 2.)*(db -  ((double)(3*vki[1]  +   vki[4])) +   1.)*gsl_matrix_get( mr, m,  b - (3*vki[1] +   vki[4])  +  4)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);


  return( summa ) ;
}


double fm_switch_number_six(  int b, int m, double db, double dm,  int  vki[],   gsl_matrix * mr )
{

  // k_1 > k_2 = k_3 = k_4
  
  double summa = 0. ;

  /* (0, 0, b )   */
  summa = summa  +  (b < m - 3 ?  (dm - db)*(dm - db - 1.)*(dm - db - 2.)*(dm - db - 3.)*gsl_matrix_get(  mr, m, b)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);

  /*  (1, 0, - ) */ 
  summa = summa   +   ( b > vki[1] - 1 ?  (db  -  ((double)vki[1])  + 1.)*(dm  -  db  +  ((double)vki[1])  - 1.)*(dm  -  db  +  ((double)vki[1])  - 2.)*(dm  -  db  +  ((double)vki[1])  - 3.)*gsl_matrix_get( mr, m, b - vki[1] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);


  /* (0, 1, - )  */
  summa = summa   +  ( b > vki[2] - 1 ?  3.*(db  -  ((double)vki[2])  + 1.)*(dm  -  db  +  ((double)vki[2])  - 1.)*(dm  -  db  +  ((double)vki[2])  - 2.)*(dm  -  db  +  ((double)vki[2])  - 3.)*gsl_matrix_get( mr, m, b - vki[2] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.); 


   /* (0, 2, - )  */

  summa = summa  +  (b >  vki[3] + vki[2] - 1 ?  3.*(db  -  ((double)(vki[3] + vki[2])) + 2.)*(db  -  ((double)(vki[3] + vki[2])) + 1.)*(dm  -  db  +   ((double)(vki[3] + vki[2]))  -  2.)*(dm  -  db  +  ((double)(vki[3] + vki[2])) -  3.)*gsl_matrix_get( mr, m,  b - vki[3] - vki[2]  +  2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);


  /*  (1, 1, - )  */

  summa = summa +  (b > vki[1] + vki[4] - 1 ?   3.*(db  -  ((double)(vki[1] + vki[4])) + 2.)*(db - ((double)(vki[1] + vki[4])) + 1.)*(dm  -  db  +   ((double)(vki[1] + vki[4]))  -  2.)*(dm  -  db  +  ((double)(vki[1] + vki[4])) -  3.)*gsl_matrix_get( mr, m,  b - vki[1] - vki[4]  +  2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);

  /*  (0, 3, - )  */
  summa = summa   +  (b > 3*vki[2] - 1 ?  (db -  ((double)(3*vki[2])) + 3.)*(db -  ((double)(3*vki[2])) +  2.)*(db -  ((double)(3*vki[2])) + 1.)*(dm - db  +   ((double)(3*vki[2]))  - 3.)*gsl_matrix_get( mr, m, b -  3*vki[2]  + 3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);  


  /*  (1, 2, - )  */
  summa = summa  +  (b > 2*vki[2]  + vki[1]  - 1 ?  3.*(db -  ((double)(2*vki[2]  +   vki[1])) + 3.)*(db -  ((double)(2*vki[2]  +   vki[1])) +  2.)*(db -  ((double)(2*vki[2]  +  vki[1])) + 1.)*(dm - db  +   ((double)(2*vki[2]  +   vki[1]))  - 3.)*gsl_matrix_get( mr, m, b -  2*vki[2]   -   vki[1]   + 3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);


  /*  (1, 3, - )   */
  summa = summa  +  ( b >  3*vki[2]  +  vki[1] - 1 ?   (db -  ((double)(3*vki[2]  +  vki[1]))  +  4. )*(db -  ((double)(3*vki[2]  +   vki[1])) +  3.)*(db -  ((double)(3*vki[2]  +   vki[1])) + 2.)*(db -  ((double)(3*vki[2]  +   vki[1])) +   1.)*gsl_matrix_get( mr, m,  b - (3*vki[2] +   vki[1])  +  4)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);

  return ( summa ) ;

}



double fm_switch_number_seven(  int b, int m, double db, double dm,  int  vki[],   gsl_matrix * mr )
{

  /*  k_1 = k_2  = k_3  =  k_4  */
  double summa = 0. ;
  
  /* (0, b)  */
  summa = summa  +    (b < m - 3 ?  (dm - db)*(dm - db - 1.)*(dm - db - 2.)*(dm - db - 3.)*gsl_matrix_get( mr, m, b)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);

  /*  (1, - ) */
  summa = summa  +    ( b > vki[1] - 1 ?  4.*(db  -  ((double)vki[1])  + 1.)*(dm  -  db  +  ((double)vki[1])  - 1.)*(dm  -  db  +  ((double)vki[1])  - 2.)*(dm  -  db  +  ((double)vki[1])  - 3.)*gsl_matrix_get( mr, m, b - vki[1] + 1)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.) ;

  /* (2, - )  */

  summa = summa + (b >  2*vki[1] - 1 ?  6.*(db -  ((double)(2*vki[1]))  +  2.)*(db -  ((double)(2*vki[1]))  +  1.)*(dm - db  +   ((double)(2*vki[1]))   -  2.)*(dm - db  +   ((double)(2*vki[1]))   -  3.) * gsl_matrix_get( mr, m, b -  (2*vki[1]) +  2)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);
  
  /* (3, - )  */

  summa = summa  +  (b >  3*vki[1]  -  1 ?  4.*(db -  ((double)(3*vki[1]))  +  3.)*(db -  ((double)(3*vki[1]))  +  2.)*(db -  ((double)(3*vki[1]))  +  1.)*(dm - db  +   ((double)(3*vki[1]))   -  3.) *gsl_matrix_get(mr,  m,     b - (3*vki[1]) + 3)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);

  /* (4, -)  */
  summa = summa +  (b >   4*vki[1]  -  1 ?  (db -  ((double)(4*vki[1]))  +  4.)*(db -  ((double)(4*vki[1]))  +  3.)*(db -  ((double)(4*vki[1]))  +  2.)*(db -   ((double)(4*vki[1]))    +  1.) * gsl_matrix_get( mr, m,  b - (4*vki[1]) + 4)/( dm*(dm-1.)*(dm - 2.)*(dm - 3.) ) : 0.);


  return( summa ) ;
}


static void pnkb( int N,   double cstikar[],     gsl_vector * lambdan,    gsl_matrix * Mgnm,   gsl_matrix * MpNkb  )
{

  //// if coalstikar[0] = 0 ( < 1 ) then Beta model
  //// if coalstikar[0] = 1 ( > 0 ) then Dirac model
  //// double coalrate( int in,  double coalstikar[],  int Vki[] ) 
    // lambdan is  vector of total rates of coalescence for each n 

  int m, k, n, b,  k2, k3, k4 ;
  int vki[5] = {0, 0, 0, 0, 0 } ;
  double  db, dm ;

  // double Betacoalrate( int b, double alpha,  int vki[] )
  
  for( k = 2 ; k < N ; k++ ){
    // temp matrix for fixed k
    gsl_matrix * Mnb  =  gsl_matrix_calloc(  N + 1, N + 1 ) ; 
    for( n = k ; n <= N ; n++ ){
      for( b = 1 ; b <= n - k + 1 ; b++ ){
	db = (double)b ; 
	gsl_matrix_set( Mnb, n, b,  (n == k ? (b == 1 ? 1. : 0. ) : 0. ) ) ;
	for( m = k ; m < n ; m++ ){
	  dm = (double)m ;
	  // add terms due to one merger
	  vki[1] = n - m + 1 ; 
	  vki[2] = vki[3] = vki[4] = 0 ;
	  //// 
	  ////
	  gsl_matrix_set( Mnb, n, b,  gsl_matrix_get( Mnb,n,b)  +   ( (b < m ? (dm - db) * gsl_matrix_get(Mnb, m, b)/dm : 0.)  +   (b > vki[1] - 1 ? (db - (double)vki[1]  +  1.)*gsl_matrix_get( Mnb, m, b - vki[1] + 1)/dm : 0.) ) * Sconst(n, vki)*coalrate( n, cstikar, vki)*gsl_matrix_get( Mgnm, m, k)/( gsl_matrix_get( Mgnm, n, k) * gsl_vector_get( lambdan, n) ) )  ; 
	  // check if can jump from n to m in 2 sim mergers
	  // need  r <=  m  <=  n - r
	  if( ((n - m > 1) && (m > 1)) && ((m < n-1) && (n > 3))  ){
	    for( k2 = 2 ; k2 <= (n - m + 2)/2 ; k2++ ){
	      // k1 =  n - m + 2 - k2
	      vki[1] = n - m + 2 - k2 ; 
	      vki[2] = k2 ; vki[3] = vki[4] = 0 ; 
	      // add terms due to 2 sim merger
	      // if merger \lambda_1 involves b lines, 
	      // \lambda_1 > \lambda_2 > 1
	      //printf("%d %d %d %d %d %d\n", n, m, k, b,  vki[1], vki[2] ) ;
	      switch( (vki[1] > vki[2] ? 0 : 1) ){
	      case 0 :
		{
		  gsl_matrix_set( Mnb, n, b,  gsl_matrix_get( Mnb, n, b)  +    tw_switch_number_zero(  b, m, db, dm,  vki, Mnb  ) * Sconst(n, vki)*coalrate( n, cstikar, vki)*gsl_matrix_get( Mgnm, m, k)/( gsl_matrix_get( Mgnm, n, k) * gsl_vector_get( lambdan, n) )) ; 
		}
		break ;
	      case 1 :
		{
		// vki[1] = vki[2] > 1
		  gsl_matrix_set( Mnb, n, b,  gsl_matrix_get( Mnb, n, b)  +   tw_switch_number_zero(  b, m, db, dm,  vki, Mnb  ) * Sconst(n,vki)*coalrate( n, cstikar, vki) *  gsl_matrix_get( Mgnm, m, k) / ( gsl_matrix_get( Mgnm, n, k)* gsl_vector_get( lambdan, n)) ) ; 
		}
		break ;
	      default :
		break ; 
	      }
	    } }
	  //////// done computations for 2 sim mergers  
	  ////
	  //// check if can jump from n to m in 3 sim merger
	  //// // need  r <=  m  <=  n - r
	  if( ((n - m + 3 > 5) && (2 < m)) && ((m < n-2) && (n > 5)) ){
	    // add terms due to 3 sim mergers
	    ////
	    for( k3 = 2 ; k3 <= (n-m+3)/3 ; k3++){
	      for( k2 = k3 ; k2 <= (n - m - k3 + 3)/2 ; k2++ ){
		vki[1]  =  n - m + 3 - k2 - k3 ; 
		vki[2]  =  k2 ; vki[3] = k3 ; vki[4] = 0 ; 
		// add terms due to 3 sim mergers
		// int th_switch_number( int VK[] )
		switch( th_switch_number( vki ) ) {
		case 0 :
		  {
		    // all mergers  same size
		    //  th_switch_number_zero(  int b, int m, double db, double dm,  int vki[],   gsl_matrix * mr )
		    gsl_matrix_set( Mnb, n, b,  gsl_matrix_get(Mnb, n, b) +   th_switch_number_zero( b, m, db, dm, vki, Mnb)   * Sconst(n,vki)*coalrate( n, cstikar, vki) *  gsl_matrix_get( Mgnm, m, k) / ( gsl_matrix_get( Mgnm, n, k)* gsl_vector_get( lambdan, n) ) ) ; 
		  }
		  break ;
		case 1 :
		  {
		    gsl_matrix_set( Mnb, n, b,  gsl_matrix_get( Mnb, n, b)  +    th_switch_number_one( b, m, db, dm, vki, Mnb) * Sconst(n,vki)*coalrate( n, cstikar, vki) *  gsl_matrix_get( Mgnm, m, k) / ( gsl_matrix_get( Mgnm, n, k)* gsl_vector_get( lambdan, n) ) ) ;  
		  }
		  break ;
		case 2 :
		  {
		    // \lambda_1  >  \lambda_2  =  \lambda_3
		    gsl_matrix_set( Mnb, n, b,  gsl_matrix_get( Mnb, n, b)  +    th_switch_number_two( b, m, db, dm, vki, Mnb) * Sconst(n,vki)*coalrate( n, cstikar, vki) *  gsl_matrix_get( Mgnm, m, k) / ( gsl_matrix_get( Mgnm, n, k)* gsl_vector_get( lambdan, n) ) ) ; 
		  }
		  break ;
		case 3 :
		  {
		      // \lambda_1  >  \lambda_2  >  \lambda_3
		    gsl_matrix_set( Mnb, n, b,  gsl_matrix_get( Mnb, n,b)  +    th_switch_number_three( b, m, db, dm, vki, Mnb)  * Sconst(n,vki)*coalrate( n, cstikar, vki) *  gsl_matrix_get( Mgnm, m, k) / ( gsl_matrix_get( Mgnm, n, k)* gsl_vector_get( lambdan, n) ) ) ; 
		  }
		  break ;
		default :
		  break ;
		}
	      } } }
	  
	  ////  done with computations for  3 sim mergers	      
	  ////
	  //// check if can jump from n to m in 4 sim mergers
	  // need  r <=  m  <=  n - r
	  if( ((n - m + 4 > 7)  &&  (3 < m))  &&  ((m < n - 3)  && (7 < n)) ){
	    for( k4 = 2 ; k4 <= (n-m+4)/4 ; k4++){
	      for( k3 = k4 ; k3 <= (n-m+4 - k4)/3 ; k3++){
		for( k2 = k3 ; k2 <= (n-m+4 -k4 - k3)/2 ; k2++ ){
		  vki[1] = n - m + 4 - k2 - k3 - k4 ; 
		  vki[2] = k2 ; vki[3] = k3 ; vki[4] = k4 ; 
		  ////
		  // printf("%d %d %d %d :  %d %d %d %d\n", n,  m, k, b,   vki[1], vki[2], vki[3], vki[4]);
		  //// add terms due to 4 sim mergers
		  /*     
		     return( ((VK[1] > VK[2]) && (VK[2] > VK[3])) && (VK[3] > VK[4]) ? 0 : ( (VK[1] == VK[2]) && ((VK[2] > VK[3]) && (VK[3] > VK[4])) ? 1 : ( (VK[1] > VK[2]) &&      ((VK[2]  ==  VK[3]) && (VK[3] > VK[4])) ? 2 : ( (VK[1] > VK[2]) &&  ((VK[2] > VK[3]) && (VK[3]  ==  VK[4])) ? 3 : ( (VK[1]  ==  VK[2]) && ((VK[2] > VK[3])       && (VK[3]  ==  VK[4])) ? 4 : ( (VK[1] ==  VK[2]) && ((VK[2]  ==  VK[3]) && (VK[3] > VK[4])) ? 5 : ( (VK[1] > VK[2]) && ((VK[2]  ==  VK[3])  && (VK[3]  ==  VK[4])) ? 6 : 7   ) ) ) ) ) ) ) ; */
		  // int fm_switch_number( int VK[] )
		  switch( fm_switch_number( vki ) ){
		  case 0 :
		    {
		      /* double fm_switch_number_zero( int b, int m, double db, double dm,  int vki[],   gsl_matrix * mr ) */
		      /* if( ((vki[1] > vki[2]) && (vki[2] > vki[3])) && (vki[3] > vki[4]) ) */
		      // \lambda_1 > \lambda_2 > \lambda_3 > \lambda_4 
		      gsl_matrix_set( Mnb, n, b, gsl_matrix_get( Mnb, n, b)  +  fm_switch_number_zero( b, m, db, dm, vki, Mnb ) * (  Sconst(n,vki)*coalrate( n, cstikar, vki) *  gsl_matrix_get( Mgnm, m, k) / ( gsl_matrix_get( Mgnm, n, k)* gsl_vector_get( lambdan, n))) ) ;
		    }
		    break ; 
		  case 1 :
		    {
		      //// \lambda_1 = \lambda_2  > \lambda_3 > \lambda_4
		      gsl_matrix_set( Mnb, n, b, gsl_matrix_get(Mnb, n, b)  +    fm_switch_number_one( b, m, db, dm, vki, Mnb ) *  Sconst(n,vki)*coalrate( n, cstikar, vki) *  gsl_matrix_get( Mgnm, m, k) / ( gsl_matrix_get( Mgnm, n, k)* gsl_vector_get( lambdan, n) )) ; 
		    }
		    break ;
		  case 2 :
		    {
		      //// \lambda_1 > \lambda_2  =  \lambda_3 > \lambda_4
		      gsl_matrix_set( Mnb, n, b,  gsl_matrix_get( Mnb, n, b)   +    fm_switch_number_two( b, m, db, dm, vki, Mnb ) *  Sconst(n,vki)*coalrate( n, cstikar, vki) *  gsl_matrix_get( Mgnm, m, k) / ( gsl_matrix_get( Mgnm, n, k)* gsl_vector_get( lambdan, n)) ) ; 
		    }
		    break ;
		  case 3 : 
		    {
		      //// \lambda_1 > \lambda_2  >  \lambda_3  ==  \lambda_4
			////
		      gsl_matrix_set( Mnb, n, b, gsl_matrix_get( Mnb, n,b)  +   fm_switch_number_three( b, m, db, dm, vki, Mnb ) * Sconst(n,vki)*coalrate( n, cstikar, vki) *  gsl_matrix_get( Mgnm, m, k) / ( gsl_matrix_get( Mgnm, n, k)* gsl_vector_get( lambdan, n) )) ; 
		    }
		    break ;
		  case 4 :
		    {
		      gsl_matrix_set( Mnb, n, b, gsl_matrix_get( Mnb, n, b)  +   fm_switch_number_four( b, m, db, dm, vki, Mnb ) *  Sconst(n,vki)*coalrate( n, cstikar, vki) *  gsl_matrix_get( Mgnm, m, k) / ( gsl_matrix_get( Mgnm, n, k)* gsl_vector_get( lambdan, n))) ;     
		    } 
		    break ;
		  case 5 :
		    {
		      gsl_matrix_set( Mnb, n, b, gsl_matrix_get( Mnb, n, b)  +   fm_switch_number_five( b, m, db, dm, vki, Mnb ) *  Sconst(n,vki)*coalrate( n, cstikar, vki) *  gsl_matrix_get( Mgnm, m, k) / ( gsl_matrix_get( Mgnm, n, k)* gsl_vector_get( lambdan, n))) ; 
		    }
		    break ;
		  case 6 :
		    {
		       gsl_matrix_set( Mnb, n, b, gsl_matrix_get( Mnb, n, b)  +   fm_switch_number_six( b, m, db, dm, vki, Mnb ) *  Sconst(n,vki)*coalrate( n, cstikar, vki) *  gsl_matrix_get( Mgnm, m, k) / ( gsl_matrix_get( Mgnm, n, k)* gsl_vector_get( lambdan, n))) ; 
		    }
		    break ;
		  case 7 :
		    {
		      gsl_matrix_set( Mnb, n, b, gsl_matrix_get( Mnb, n, b)  +   fm_switch_number_seven( b, m, db, dm, vki, Mnb ) *  Sconst(n,vki)*coalrate( n, cstikar, vki) *  gsl_matrix_get( Mgnm, m, k) / ( gsl_matrix_get( Mgnm, n, k)* gsl_vector_get( lambdan, n))) ; 
		    }
		    break ;
		  default :
		    break ;
		  }
	  // case of 4 mergers  computed
		} } } } } } } 
    //// update Nkb matrix for p^{( N )}[k, b] values 
    //  b = 1 ; b <= N - k + 1 ; b++
    for( b = 1 ; b <=  N  - k +  1 ; b++ ){
      gsl_matrix_set( MpNkb,  k, b,  gsl_matrix_get( Mnb, N, b) ) ; }
    gsl_matrix_free( Mnb ) ; }

  gsl_matrix_set( MpNkb,  N, 1, 1. ) ; 
}







void EBone( int N, double a )
{
  // expected length of external branches  by recursion

  gsl_matrix * E =  gsl_matrix_calloc( N+1, N + 1 ) ;
  // double Betacoalrate( int b, double alpha,  int vki[] )
  //  Sconst( int n,  int x[] )

  // vln is vector of sum of rates for each n
  int k1, k2, k3, k4, i, j, n ;
  int vki[5] = {0, 0, 0, 0, 0 } ;
  vki[1] = 2 ;
  gsl_matrix_set( E, 2, 1, 1./(Betacoalrate( 2, a, 1.0, 1.0,  vki)*Sconst( 2, vki )) ) ;
  gsl_matrix_set( E, 2, 2, 2./(Betacoalrate( 2, a, 1.0, 1.0,  vki)*Sconst( 2, vki )) ) ;

  gsl_vector * vln = gsl_vector_calloc( N + 2 ) ;

  // void Betatotalrate( int N,  double  ialpha,  gsl_vector * vln )
  // Betatotalrate( N, a,    vln ) ;
  
  
  for( n = 3 ; n <= N ; n++ ){
    for( i = 1 ; i <= n ; i++ ){
      gsl_matrix_set( E, n, i,  ((double)i)/gsl_vector_get( vln, n) ) ;
      for( k1 = 2 ; k1 <= n  ; k1++ ){
	//printf( "%d\n", k1 ) ;
	vki[1] = k1 ; vki[2] = vki[3] = vki[4] = 0 ;
	for( j = 0 ; j <= (k1 < i ? k1 : i) ; j++ ){
	  gsl_matrix_set( E, n, i, gsl_matrix_get( E, n, i)  +      (i - j > n - k1 ? 0. :  gsl_sf_choose( k1, j) *gsl_sf_choose( n-k1, i - j) * (Betacoalrate( n, a, 1.0, 1.0,  vki)*Sconst(n, vki)) * gsl_matrix_get(E, n - vki[1] + 1, i - j)/( gsl_sf_choose(n, i) * gsl_vector_get( vln, n)) ) )  ; } }


      if( n > 3 ){
	for( k1 = 2; k1 <=  n - 2 ; k1++){
	  for( k2 = 2 ; k2 <= (n-k1 < k1 ? n-k1 : k1) ; k2++ ){
	    // printf( "%d %d\n", k1, k2 ) ;
	    vki[1] = k1 ; vki[2] = k2 ; vki[3] = vki[4] = 0 ;
	    for( j = 0 ; j <= (k1 + k2 < i ? k1 + k2 : i) ; j++ ){
	      gsl_matrix_set( E, n, i, gsl_matrix_get( E, n, i)  +   (i - j > n - k1 - k2 ? 0. : gsl_sf_choose( k1 + k2, j)*gsl_sf_choose(n- k1 - k2, i - j) * (Betacoalrate( n, a, 1.0, 1.0,  vki)*Sconst(n, vki)) *  gsl_matrix_get(E, n - vki[1] - vki[2] + 2, i - j)/(gsl_sf_choose(n,i)*gsl_vector_get( vln, n) ) ) ) ; }}}}

      
      if( n > 5 ){
	for( k1 = 2 ; k1 <= n - 4 ; k1++ ){
	  for( k2 = 2 ; k2 <= (n-k1-2 < k1 ? n - k1 - 2 : k1 ); k2++ ){
	    for( k3 = 2 ; k3 <= (n-k1-k2 < k2 ? n-k1-k2 : k2); k3++ ){
	      vki[1] = k1 ; vki[2] = k2 ; vki[3] = k3 ; vki[4] = 0 ;
	      //  printf( "%d %d %d\n", k1, k2, k3 ) ; 
	      for( j = 0 ; j <= ( k1 + k2 + k3 < i ? k1 + k2 + k3 : i ) ;  j++ ){
		gsl_matrix_set( E, n, i, gsl_matrix_get( E, n, i)  +  (i - j > n - k1 - k2 - k3 ? 0. : gsl_sf_choose( k1 + k2+k3, j)*gsl_sf_choose(n- k1 - k2 - k3, i - j) * (Betacoalrate( n, a, 1.0, 1.0,  vki)*Sconst(n, vki)) * gsl_matrix_get(E, n - vki[1] - vki[2] - vki[3] + 3, i - j)/(gsl_sf_choose(n,i)*gsl_vector_get( vln, n) ) ) ) ; }}}}}
	   
  
  
      if( n > 7){
	for( k1 = 2 ; k1 <= n-6 ; k1++ ){
	  for( k2 = 2 ; k2 <= (n-k1-4 < k1 ? n-k1-4 : k1); k2++ ){
	    for( k3 = 2 ; k3 <= (n-k1-k2 - 2 < k2 ? n-k1-k2 - 2 : k2 ); k3++ ){
	      for( k4 = 2 ; k4 <= (n-k1-k2-k3 < k3 ? n-k1-k2-k3 : k3); k4++){
		//  printf( "%d %d %d %d\n", k1, k2, k3, k4 ) ;
		vki[1] = k1 ; vki[2] = k2 ; vki[3] = k3 ; vki[4] = k4 ;
		for( j = 0 ; j <= ( k1 + k2 + k3 + k4 < i ? k1 + k2 + k3 + k4 : i) ; j++){
		  gsl_matrix_set( E, n, i, gsl_matrix_get( E, n, i)  +  (i - j > n - k1 - k2 - k3 - k4 ? 0. : gsl_sf_choose( k1 + k2+k3 + k4, j)*gsl_sf_choose(n- k1 - k2 - k3 - k4, i - j) * (Betacoalrate( n, a, 1.0, 1.0,  vki)*Sconst(n, vki)) * gsl_matrix_get(E, n - vki[1] - vki[2] - vki[3] - vki[4] + 4, i - j)/(gsl_sf_choose(n,i)*gsl_vector_get( vln, n) ) ) ) ; }}}}} } } }

  printf("%g\n" ,  gsl_matrix_get( E, N, N) ) ; 
  
 
  gsl_matrix_free( E ) ; 
  gsl_vector_free( vln ) ;

}



static void profun( int Ni, double cs_0, double cs_1,  double cs_2,   double cs_3,  gsl_vector * vphi )
{

  /*  if cs_0  < 1 (= 0) then Beta model
      if cs_[0] > 1   then Dirac model 
  ------------
   if cs_0 > 1  then Dirac model, so then    
        cs_1 = c, cs_2 = psi 
   --------------------
  if cs_0 < 1 then Beta-model, so 
  cs_3 = alpha ( can put cs_1 = cs_2 = 0)
 */
  // int y[] = {0, 2, 0, 0, 0 } ;
  //double  Sconst( int n,  int x[] )
  //printf("%g\n",   Sconst( 4, y ) ) ;
  // double Cconst( int n,  double psi, int x[]  )
  //printf("%g\n",   Cconst( 4, 0.5,   y ) ) ;
  // void totalrate( int N, double psi,  double cstiki,   gsl_vector * vln )
  gsl_vector * vlambdan  =  gsl_vector_calloc( Ni + 2 ) ; 
  ////  void pnkb( int N,   gsl_matrix * Mgnm,   gsl_matrix * Mpnkb,  double cstikar[],   gsl_vector * lambdan )
  ////  void Mqp( int N,  double psi, gsl_matrix * Mq, gsl_matrix * Mp )
  ////  void totalrate( int N,  double cstikar[],  gsl_vector * vln )

  /* 0, alpha, K, C */
  double cst[4] = { 0.0, cs_1, cs_2, cs_3 } ;

  totalrate( Ni, cst,  vlambdan ); 

   // void Mqp( int N, double icoals[],   gsl_matrix * Mq, gsl_matrix * Mp )
   // void gnm( int N,   gsl_matrix * Mpp,   gsl_matrix * Mgnm )
  

  gsl_matrix * MQ  =  gsl_matrix_calloc( Ni + 1,  Ni + 1 ) ;
  gsl_matrix * MP  =  gsl_matrix_calloc( Ni + 1,  Ni + 1 ) ;
  gsl_matrix * MGNM  =  gsl_matrix_calloc( Ni + 1, Ni + 1 ) ;
  gsl_matrix * Mr  =  gsl_matrix_calloc( Ni + 1, Ni + 1 ) ;

  Mqp( Ni, cst,  MQ, MP ) ; 
  gnm( Ni, MP,  MGNM ) ;  
  
  //  prentagslfylki( int z,  gsl_matrix * m )
  // prentagslfylki( Ni,   MGNM )  ; 

  ////  void pnkb( int N,   double cstikar[],   gsl_vector * lambdan,    gsl_matrix * Mgnm,   gsl_matrix * Mpnkb  )
  pnkb( Ni, cst,  vlambdan,  MGNM, Mr ) ;  

  int i, k ;
  double s = 0. ;
  
   
  for( i = 1 ; i < Ni ; i++ ){
    gsl_vector_set( vphi, i, 0. ) ;
    for ( k = 2 ; k <= Ni - i + 1 ; k++){
      gsl_vector_set( vphi, i,  gsl_vector_get(vphi, i)  +     gsl_matrix_get( Mr, k, i) * ((double)k) * gsl_matrix_get( MGNM,  Ni, k ) ) ; }
    s = s  +  gsl_vector_get( vphi, i ) ; }

  // normalise with expected total length s = E[B]; or not 
  for( i = 1 ; i < Ni ; i++ ){
    gsl_vector_set( vphi, i,  gsl_vector_get( vphi, i) / s ) ; }
    
  

  
  /*
     int u, v ; 
     double  rowsum = 0. ;
     for( u = 2 ; u <= Ni ; u++ ){
       printf("%d : ", u ) ;
       rowsum = 0. ;
       for( v = 1 ; v < Ni ; v++ ){
	 rowsum = rowsum  +   gsl_matrix_get( Mr, u, v) ; 
	 printf("%g ", gsl_matrix_get( Mr, u, v) ) ; }
       printf(" : %g\n",  rowsum ) ; }
  */   

  gsl_vector_free( vlambdan ) ; 
  gsl_matrix_free( MQ ) ;
  gsl_matrix_free( MP ) ;
  gsl_matrix_free( MGNM ) ;
  gsl_matrix_free( Mr ) ;
}



// static void printvarphi( int lauf,  double process, double cparameter, double psiparameter, double alphaparameter )
static void printvarphi( int lauf,  const double alphaparameter, const double Kparam, const double Cparam )
{
  // print out results from profun for cwebfile
  // void profun( int Ni, double cs_0, double cs_1,  double cs_2,   double cs_3,  gsl_vector * vphi )
  //// if process  < 1  then Beta model
  //// otherwise    Dirac model
  

  gsl_vector * v = gsl_vector_calloc( lauf + 2 ) ;


  /* static void profun( int Ni, double cs_0, double cs_1,  double cs_2,   double cs_3,  gsl_vector * vphi ) */
  // profun( lauf, process, cparameter, psiparameter, alphaparameter, v ) ;
  profun( lauf, 0,  alphaparameter, Kparam,  Cparam,   v ) ; 

  int j ;
  for( j = 1 ; j < lauf ; j++ ){
    printf("%f\n", gsl_vector_get( v, j)) ; }

  gsl_vector_free( v ) ; 

}



void xivarphifylki( void )
{
  // void profun( int Ni, double cs_0, double cs_1,  double cs_2,   double cs_3,  gsl_vector * vphi )
  //  //// if coalstikar[0] = 0 ( < 1 ) then Beta model
  //// if coalstikar[0] = 1 ( > 0 ) then Dirac model
  // coalstikar[1] = cs_1 then is either \alpha or \psi

  gsl_vector * v  =  gsl_vector_calloc( 52 ) ;

  double a = 0.0 ;
  int i ;
  while( a < 0.99999 ){
    gsl_vector_set_zero( v ) ;
    profun( 10, 1, a, 0, 0, v ) ; 
    for( i = 1 ; i < 10 ; i++ ){
      printf("%g ", gsl_vector_get(v,i)) ; }
    printf("\n") ;
    a = a + 0.025 ; }

  gsl_vector_free( v) ;

}



/* modules from simxi.c  */



gsl_rng * rngtype ;
void setup_rng( unsigned long int seed )
{
  //
  // *gsl_rng_alloc( const gsl_rng_type *T )
  // const gsl_rng_type *T ;
  //gsl_rng_env_setup();
  // T = gsl_rng_default ;
  // rngtype = gsl_rng_alloc( T ) ;
  //
  gsl_rng_default_seed = seed ;
  rngtype = gsl_rng_alloc(gsl_rng_mt19937);
  // gsl_rng_default_seed = rngseed ;
}




void blocks_merge( int N,  gsl_vector_int * ni, gsl_vector_int * coallabel,   gsl_ran_discrete_t * teningur, const gsl_rng * r )
{

  // sample branches to merge
  // prob(X=0) = 1 - psi
  // prob(X=i) = psi/4 for 0 < i < 5.
  // N is number of leaves

  int i, k ;
  // vki records the number of lines associated with each of groups 1 to 4
  int vki[] = {0, 0, 0, 0, 0} ;
  int size[] = {0, 0, 0, 0, 0 } ;

  for( i = 1 ; i <= N ; i++ ){
    // coallabel[i] is label of the coalescing  group which block i becomes associated with
    gsl_vector_int_set(coallabel, i, ( gsl_vector_int_get( ni, i) > 0 ?  gsl_ran_discrete( r, teningur ) : 0 ) )  ; 
    vki[ gsl_vector_int_get( coallabel,i) ]  =  vki[ gsl_vector_int_get( coallabel,i) ] + 1 ; 
    size[ gsl_vector_int_get( coallabel,i) ] =  size[ gsl_vector_int_get( coallabel,i) ]  +  ( gsl_vector_int_get( coallabel,i) > 0 ? gsl_vector_int_get( ni, i) : 0 ) ; }
  
  if( ((vki[1] > 1) || (vki[2] > 1)) || ((vki[3] > 1) || (vki[4] > 1)) ){
    for( i = 1 ; i <= N ; i++){
      gsl_vector_int_set( ni, i,  ( gsl_vector_int_get( coallabel, i) > 0 ? (vki[ gsl_vector_int_get(coallabel, i) ] > 1 ? 0 : gsl_vector_int_get( ni, i)) :      gsl_vector_int_get( ni, i)) ) ; }

    
    i = 1;
    k = ( vki[1] > 1 ? 1 : (vki[2] > 1 ? 2 : (vki[3] > 1 ? 3 : 4 ))) ;
    while( (i <=  N)  &&  (k < 5) ){
      if( gsl_vector_int_get( ni,i) == 0 ){
	gsl_vector_int_set( ni, i,  size[k] ) ; 
	k = ( k == 1 ? (vki[2] > 1 ? 2 : (vki[3] > 1 ? 3 : (vki[4] > 1 ? 4 : 5))) : (k == 2 ? (vki[3] > 1 ? 3 : (vki[4] > 1 ? 4 : 5)) : ( k == 3 ? ( vki[4] > 1 ? 4   : 5) : 5 ) ) ) ; }
      i = i + 1 ; } 
    
    // update total number of active blocks
    gsl_vector_int_set( ni, 0,  gsl_vector_int_get( ni, 0)  -  (vki[1] > 1 ? vki[1] - 1 : 0) -   (vki[2] > 1 ? vki[2] - 1 : 0)  -  (vki[3] > 1 ? vki[3] - 1 : 0)  -  (vki[4] > 1 ? vki[4] - 1 : 0) ) ;    }
}


void tvo_blocks_merge( int N, gsl_vector_int * ni,   int * tni , int * ctveir, const gsl_rng * r )
{

  // small reproduction event, so forcing 2 blocks to merge
  // value of ni[i] is size of block number i

  int i, j ;
  j = 0 ;
  for( i = 1 ; i <= N ; i++){
    if( gsl_vector_int_get( ni,i) > 0 ){
      tni[j] =  i ; 
      j = j + 1 ; }}
 
  // draw 2 blocks to coalesce from; 
  // put the blocks in ctveir
  gsl_ran_choose( r, ctveir,  2,  tni,  j,  sizeof( int ) ) ;

  gsl_vector_int_set( ni, ctveir[0],  gsl_vector_int_get( ni, ctveir[0])  +  gsl_vector_int_get( ni, ctveir[1]) ) ; 

  gsl_vector_int_set( ni, ctveir[1], 0 ) ; 

  // update current number of active blocks
  gsl_vector_int_set( ni, 0,   gsl_vector_int_get( ni, 0)  - 1 ) ; 
}




void Diracsimx( int N, double coalparameters[],   int B,   const gsl_rng * r )
{

  // initially Dirac coalescent
  int biters = 0 ;
  int n ; 
  double p = coalparameters[2]/( 1.  +  coalparameters[2] ) ;
  double t = 0. ;

  gsl_vector * vbi = gsl_vector_calloc( N + 2 ) ; 
  gsl_vector * svbi = gsl_vector_calloc( N + 2 ) ; 
  gsl_vector_int * ni = gsl_vector_int_calloc( N + 2 ) ;
  gsl_vector * lambdan = gsl_vector_calloc( N + 2 ) ; 
  int * tempni  = (int *)calloc( N + 2,  sizeof( int ) ) ;  
  int * btveir  = (int *)calloc( 2,  sizeof( int ) ) ;  
  gsl_vector_int * coalmidi = gsl_vector_int_calloc( N + 2 ) ; 
  double * P = (double *)calloc( 5, sizeof(double) ) ; 
  P[1] = P[2] = P[3] = P[4]  =  coalparameters[1] / 4. ;
  P[0] =  1. - coalparameters[1] ; 
  gsl_ran_discrete_t * dice =  gsl_ran_discrete_preproc( 5,  P ) ; 
  
  
  // compute total coal rates
  // void totalrate( int N,  double cstikar[],  gsl_vector * vln )
  totalrate( N, coalparameters, lambdan ) ;
  for( n = 2 ; n <= N ; n++ ){
     printf("%d  %g\n", n,  gsl_vector_get( lambdan, n) ) ; }

  
  while( biters < B ){

    // iterate a genealogy
    // collect branch lengths

    // the value of ni[i] is size of block i
    for( n = 1 ; n <= N ; n++ ){
    gsl_vector_int_set( ni, n, 1 ) ; 
    gsl_vector_set( vbi, n, 0. ) ; }
    // value of ni[0] is current number of active blocks
    gsl_vector_int_set( ni, 0, N ) ; 

    while( gsl_vector_int_get( ni, 0)  >  1 ){
      // ni[0] is current number of blocks
       
      //  for( n = 0 ; n <= N ; n++ ){
      // printf("%d ", gsl_vector_int_get( ni, n) ) ; }
      //  printf("\n") ; 
      
      // draw time
      t = gsl_ran_exponential( r,  1. / gsl_vector_get( lambdan,  gsl_vector_int_get( ni, 0)  ) ) ;
      //printf("t %g\n", t ) ;
      
      // update branch lengths
      for( n = 1 ; n <= N ; n++ ){
	gsl_vector_set( vbi, gsl_vector_int_get(ni, n),  gsl_vector_get( vbi, gsl_vector_int_get(ni, n))  +  ( gsl_vector_int_get( ni, n) > 0 ?  t : 0. ) ) ; }
      
      // draw between small and large mergers
	// 0 means small event; 1  means  large  event
      if( (p > 0 ? gsl_ran_bernoulli( r,  p ) : 0) == 0 ){
	// small, or kingman, merger
	// draw 2 branches to merge 
	//  tvo_blocks_merge( int N, gsl_vector * ni,   int * tni , int * ctveir, const gsl_rng * r )
	tvo_blocks_merge( N, ni,  tempni,  btveir,  r ); }
      else{
	// large reproduction event
	//   blocks_merge( int N,  gsl_vector_int * ni, gsl_vector * coallabel,   gsl_ran_discrete_t * teningur, const gsl_rng * r )
	blocks_merge( N, ni,  coalmidi,  dice,  r ) ; }
     
      
 
    }
    // update the estimate of branch lengths
    for( n = 1 ; n <= N ; n++ ){
      gsl_vector_set( svbi, n,  gsl_vector_get( svbi, n)  +   gsl_vector_get( vbi, n) ) ; }


    // add to count of iterations
    biters =  biters  +  1 ; }


  for( n = 1 ; n < N ; n++ ){
    printf("%g\n", gsl_vector_get( svbi, n)/( (double)B ) ) ; }


  // clean up
  gsl_vector_free( lambdan ) ; 
  gsl_vector_free( vbi ) ;
  gsl_vector_free( svbi ) ; 
  gsl_vector_int_free( ni ) ; 
  gsl_vector_int_free( coalmidi ) ; 
  free( tempni ) ; 
  free( btveir ) ;
  free( P ) ;
  gsl_ran_discrete_free( dice ) ; 
}




void drawcoalvector(int n , double ialpha, int * sampledvki,   gsl_rng * r  )
{

  // draw a vector of k_i : number of blocks to coalesce in each group
  // n  is current  active number of blocks

  // first count  number of mergers
  int m = 0 ; 
  int k1, k2, k3, k4 ;
   for( k1 = 2 ; k1 <= n  ; k1++ ){
    //printf( "%d\n", k1 ) ;
     m = m + 1 ;
  }

  if( n > 3 ){
    for( k1 = 2; k1 <=  n - 2 ; k1++){
      for( k2 = 2 ; k2 <= (n-k1 < k1 ? n-k1 : k1) ; k2++ ){
	// printf( "%d %d\n", k1, k2 ) ; 
	m = m + 1 ;
      } } }

  
  if( n > 5 ){
    for( k1 = 2 ; k1 <= n - 4 ; k1++ ){
      for( k2 = 2 ; k2 <= (n-k1-2 < k1 ? n - k1 - 2 : k1 ); k2++ ){
	for( k3 = 2 ; k3 <= (n-k1-k2 < k2 ? n-k1-k2 : k2); k3++ ){
	  //  printf( "%d %d %d\n", k1, k2, k3 ) ; 
	  m = m + 1 ;
	} } } }
  
  
  if( n > 7){
    for( k1 = 2 ; k1 <= n-6 ; k1++ ){
      for( k2 = 2 ; k2 <= (n-k1-4 < k1 ? n-k1-4 : k1); k2++ ){
	for( k3 = 2 ; k3 <= (n-k1-k2-2 < k2 ? n-k1-k2 -2 : k2 ); k3++ ){
	  for( k4 = 2 ; k4 <= (n-k1-k2-k3 < k3 ? n-k1-k2-k3 : k3); k4++){
	    //  printf( "%d %d %d %d\n", k1, k2, k3, k4 ) ; 
	    m = m + 1 ;
	  } } } } }
  
  // m now holds number of mergers possible from  n  blocks 

  // printf("m  %d\n", m ) ;

  double * P =  (double *)calloc( m, sizeof(double) ) ;
  // add the merger rates
  //  Sconst(n, vk) * Betacoalrate( n, ialpha,  vk) 
  
  
  gsl_matrix_int * K =  gsl_matrix_int_calloc( m, 5 ) ;

  int vki[5] = {0, 0, 0, 0, 0} ;

  m = 0 ;
  for( k1 = 2 ; k1 <= n  ; k1++ ){
    //printf( "%d\n", k1 ) ;
    vki[1] = k1; vki[2] = vki[3] = vki[4] = 0 ;
    P[m]  =    Sconst(n, vki) * Betacoalrate( n, ialpha, 1.0, 1.0,   vki) ; 
    gsl_matrix_int_set( K, m,  1,  vki[1] ) ;
    m = m + 1 ; }
  

  if( n > 3 ){
    for( k1 = 2; k1 <=  n - 2 ; k1++){
      for( k2 = 2 ; k2 <= (n-k1 < k1 ? n-k1 : k1) ; k2++ ){
	// printf( "%d %d\n", k1, k2 ) ; 
	vki[1] = k1; vki[2] =  k2 ; 
	vki[3] = vki[4] = 0 ;
	P[m]  =    Sconst(n, vki) * Betacoalrate( n, ialpha, 1.0, 1.0,   vki) ; 
	gsl_matrix_int_set( K, m,  1,  vki[1] ) ;
	gsl_matrix_int_set( K, m,  2,  vki[2] ) ;
	m = m + 1 ;
      } } }

  
  if( n > 5 ){
    for( k1 = 2 ; k1 <= n - 4 ; k1++ ){
      for( k2 = 2 ; k2 <= (n-k1-2 < k1 ? n - k1 - 2 : k1 ); k2++ ){
	for( k3 = 2 ; k3 <= (n-k1-k2 < k2 ? n-k1-k2 : k2); k3++ ){
	  //  printf( "%d %d %d\n", k1, k2, k3 ) ; 
	  vki[1] = k1; vki[2] =  k2 ;  vki[3] = k3 ;
	  vki[4] = 0 ;
	  P[m]  =    Sconst(n, vki) * Betacoalrate( n, ialpha, 1.0, 1.0,  vki) ; 
	  gsl_matrix_int_set( K, m,  1,  vki[1] ) ;
	  gsl_matrix_int_set( K, m,  2,  vki[2] ) ;
	  gsl_matrix_int_set( K, m,  3,  vki[3] ) ;
	  m = m + 1 ;
	} } } }
  
  
  if( n > 7){
    for( k1 = 2 ; k1 <= n-6 ; k1++ ){
      for( k2 = 2 ; k2 <= (n-k1-4 < k1 ? n-k1-4 : k1); k2++ ){
	for( k3 = 2 ; k3 <= (n-k1-k2 < k2 ? n-k1-k2 : k2 ); k3++ ){
	  for( k4 = 2 ; k4 <= (n-k1-k2-k3 < k3 ? n-k1-k2-k3 : k3); k4++){
	    //  printf( "%d %d %d %d\n", k1, k2, k3, k4 ) ; 
	    vki[1] = k1; vki[2] =  k2 ;  vki[3] = k3 ; vki[4] = k4 ;
	    P[m]  =    Sconst(n, vki) * Betacoalrate( n, ialpha, 1.0, 1.0,   vki) ; 
	    gsl_matrix_int_set( K, m,  1,  vki[1] ) ;
	    gsl_matrix_int_set( K, m,  2,  vki[2] ) ;
	    gsl_matrix_int_set( K, m,  3,  vki[3] ) ;
	    gsl_matrix_int_set( K, m,  4,  vki[4] ) ;
	    m = m + 1 ;
	  } } } } }
  

  gsl_ran_discrete_t * g =  gsl_ran_discrete_preproc( m,  P ) ;
  m  =  (int)gsl_ran_discrete( r,  g ) ; 

  sampledvki[1] =  gsl_matrix_int_get( K, m, 1 ) ;
  sampledvki[2] =  gsl_matrix_int_get( K, m, 2 ) ;
  sampledvki[3] =  gsl_matrix_int_get( K, m, 3 ) ;
  sampledvki[4] =  gsl_matrix_int_get( K, m, 4 ) ;


  

  free( P ) ;
  gsl_ran_discrete_free( g ) ;
  gsl_matrix_int_free( K ) ;

}



void  Beta_blocks_merge(  int N,    int b,  gsl_vector_int * Ni,   gsl_ran_discrete_t * fjorir,  int * umrodun,    gsl_rng * r )
{

  // N is leaves
  // b is number of blocks  drawn from Beta-coal
  // to possibly merge in up to 4 groups
  // Ni[i] is size of block i 
  // fjorir is  P( X = i ) = 1/4 for 0 < i < 5

  
  int *  vki = (int *)calloc( 5, sizeof( int ) ) ;
  int size[] = {0, 0, 0, 0, 0} ;
  int groups[] = {0, 0, 0, 0, 0} ;

  int i, g ; 

  //  ( (vki[1] < 2) && (vki[2] < 2) ) && ((vki[3] < 2) && (vki[4] < 2)) 
  //   ( (vki[1] < 2) && (vki[2] < 2) ) && ((vki[3] < 2) && (vki[4] < 2))
  // 
  //while(  (vki[1] < 2 ? 0 : vki[1])  +  (vki[2] < 2 ? 0 : vki[2]) + (vki[3] < 2 ? 0 : vki[3])  + (vki[4] < 2 ? 0 : vki[4]) < 1   ){
    // retry  until at least 2 blocks merge
    // vki[1] = vki[2] = vki[3] = vki[4] = 0 ; 
  ////
  for(  i = 1 ; i <= b ; i++ ){
    // g is group index 
    g =  gsl_ran_discrete( r,  fjorir ) ; 
    vki[g] =  vki[g] + 1 ; }
  // void drawcoalvector(int n , double ialpha, int * sampledvki,   gsl_rng * r  )
  // drawcoalvector( b,  ialpha,  vki,  r ) ;
  
     
    if( ( (vki[1] > 1) || (vki[2] > 1) ) || ((vki[3] > 1) || (vki[4] > 1)) ){
      //// at least one merger
      // randomly order the  index of the blocks
      gsl_ran_shuffle( r, umrodun, N,  sizeof( int ) ) ; 
      i = 0 ;

      g = (vki[1] > 1 ? 1 : (vki[2] > 1 ? 2 : ( vki[3] > 1 ? 3 : 4) ) ) ;
      
      while( g < 5 ){
	if( gsl_vector_int_get( Ni, umrodun[i]) > 0 ){
	  // active block with label umrodun[i]
	  // add the size of the block to the total size of the group
	  size[g]  =  size[g]  +  gsl_vector_int_get( Ni, umrodun[i] ) ;
	  gsl_vector_int_set( Ni, umrodun[i], 0 ) ;
	  groups[g]  =  groups[g] + 1 ; 
	} 
	g =  (groups[g] < vki[g] ? g : (g == 1 ? (vki[2] > 1 ? 2 : (vki[3] > 1 ? 3 :( vki[4] > 1 ? 4 : 5))) : (g == 2 ? (vki[3] > 1 ? 3 : (vki[4] > 1 ? 4 : 5)) : (g == 3 ? (vki[4] > 1 ? 4 : 5) : 5 ) ) ) ) ;
	i = i + 1 ; }
      //// 
      i = 0 ;
      g = (vki[1] > 1 ? 1 : (vki[2] > 1 ? 2 : ( vki[3] > 1 ? 3 : 4) ) ) ;


      // randomly order the  index of the blocks
      gsl_ran_shuffle( r, umrodun, N,  sizeof( int ) ) ; 

      while ( g < 5 ){
	if( gsl_vector_int_get( Ni, umrodun[i] ) == 0 ){
	  gsl_vector_int_set( Ni, umrodun[i],  size[g] ) ; 
	  g = (g == 1 ? (vki[2] > 1 ? 2 : (vki[3] > 1 ? 3 :( vki[4] > 1 ? 4 : 5))) : (g == 2 ? (vki[3] > 1 ? 3 : (vki[4] > 1 ? 4 : 5)) : (g == 3 ? (vki[4] > 1 ? 4 : 5) : 5 ) ) ) ; }
	i = i + 1 ; } 

      //// update the number of  active blocks  
      gsl_vector_int_set( Ni, 0, gsl_vector_int_get( Ni, 0)  -  (vki[1] > 1 ? vki[1] - 1 : 0)  -  (vki[2] > 1 ? vki[2] - 1 : 0)  -  (vki[3] > 1 ? vki[3] - 1 : 0)  -  (vki[4] > 1 ? vki[4] - 1 : 0) ) ;  } 

    free( vki ) ;
}



void Betasimx( int N, double ialpha,    int B,    gsl_rng * r )
{

  // Beta-coalescent as Lambda-coalescent
  int biters = 0 ;
  int n,   fjoldi_blokka, currentn  ; 

  double t = 0. ;

  gsl_vector * vbi = gsl_vector_calloc( N + 2 ) ; 
  gsl_vector * svbi = gsl_vector_calloc( N + 2 ) ; 
  gsl_vector_int * ni = gsl_vector_int_calloc( N + 2 ) ;
  int * Umrodun  =  (int *)calloc( N, sizeof( int ) ) ;
  gsl_vector * lambdan = gsl_vector_calloc( N + 2 ) ; 
  int * tempni  = (int *)calloc( N + 2,  sizeof( int ) ) ;  
  int * btveir  = (int *)calloc( 2,  sizeof( int ) ) ;  
  
  double * P = (double *)calloc( 5, sizeof(double) ) ; 
  P[1] = P[2] = P[3] = P[4]  =  0.25  ;
  P[0] = 0. ;
  gsl_ran_discrete_t * dice =  gsl_ran_discrete_preproc( 5,  P ) ; 
  
  /* ************** Betasimx  */
  // double  **alloc_2d_array(  int n)
  double ** MBetalambda  =  alloc_2d_array( N ) ; 
  Umrodun[0]  = 1  ; 

  for( n = 2 ; n <= N ; n++ ){
    for( biters = 1 ; biters < n ; biters++){
      // biters is number of active blocks after jump from n blocks
      MBetalambda[n][biters]  =   gsl_sf_choose(n,n - biters + 1) * gsl_sf_beta( (double)(n-biters+1) - ialpha,  (double)(biters - 1)  +  ialpha)/gsl_sf_beta(2. -  ialpha,  ialpha) ;  }
    Umrodun[ n-1] =   n ; }

  /* ****************** Betasimx  */
  // compute total coal rates
  // void Betatotalrate( int N,  double  ialpha,  gsl_vector * vl )
  // Betatotalrate( N, ialpha, lambdan ) ;


  /*
  for( n = 2 ; n <= N ; n++ ){
       printf("%d  %g\n", n,  gsl_vector_get( lambdan, n) ) ; }
  */

  biters = 0 ; 
  while( biters < B ){

    // iterate a genealogy
    // collect branch lengths

    // the value of ni[i] is size of block i
    for( n = 1 ; n <= N ; n++ ){
    gsl_vector_int_set( ni, n, 1 ) ; 
    gsl_vector_set( vbi, n, 0. ) ; }
    // value of ni[0] is current number of active blocks
    gsl_vector_int_set( ni, 0, N ) ; 

    while( gsl_vector_int_get( ni, 0)  >  1 ){
      // ni[0] is current number of blocks
      

      /********
        printf("%d : ", gsl_vector_int_get( ni, 0) ) ;
       for( n = 1 ; n <= N ; n++ ){
	printf("%d ", gsl_vector_int_get( ni, n) ) ; }
      printf("\n") ; 
      **********/
 
      // draw time
      t = gsl_ran_exponential( r,  1. / gsl_vector_get( lambdan,  gsl_vector_int_get( ni, 0)  ) ) ;
      // printf("t %g\n", t ) ;
      
      // update branch lengths
      for( n = 1 ; n <= N ; n++ ){
	gsl_vector_set( vbi, gsl_vector_int_get(ni, n),  gsl_vector_get( vbi, gsl_vector_int_get(ni, n))  +  ( gsl_vector_int_get( ni, n) > 0 ?  t : 0. ) ) ; }
      
      /********************** Betasim ******/

      // draw total number of blocks  to  possibly merge in up to 4 groups
      // int  draw_number_blocks(  int b,  double * M, gsl_rng * r )
      currentn =  gsl_vector_int_get( ni, 0 ) ; 
      while( currentn == gsl_vector_int_get( ni, 0 ) ){
	fjoldi_blokka =  draw_number_blocks( gsl_vector_int_get( ni, 0), MBetalambda[ gsl_vector_int_get( ni, 0) ], r ) ; 
      
      //  printf("%d :  fjb %d \n", gsl_vector_int_get( ni, 0),    fjoldi_blokka ) ;
      
      // void  Beta_blocks_merge(  int N,     int b,  gsl_vector_int * Ni,   gsl_ran_discrete_t * fjorir,  int * umrodun,    gsl_rng * r )
	Beta_blocks_merge( N,    fjoldi_blokka,   ni,   dice,   Umrodun,   r ) ; } }
    // one realisation of branch lengths obtained
  
    // update the estimate of branch lengths
    for( n = 1 ; n <= N ; n++ ){
      gsl_vector_set( svbi, n,  gsl_vector_get( svbi, n)  +   gsl_vector_get( vbi, n) ) ; }

    // add to count of iterations
    biters =  biters  +  1 ; }


  for( n = 1 ; n < N ; n++ ){
    printf("%g\n", gsl_vector_get( svbi, n)/( (double)B ) ) ; }


  // clean up
  gsl_vector_free( lambdan ) ; 
  gsl_vector_free( vbi ) ;
  gsl_vector_free( svbi ) ; 
  gsl_vector_int_free( ni ) ; 
  free( Umrodun ) ; 
  
  free( tempni ) ; 
  free( btveir ) ;
  free( P ) ;
  gsl_ran_discrete_free( dice ) ; 
  // void free_2d_array( double **a, int n )
  free_2d_array( MBetalambda, N ) ; 
}



void XifourDiracsim(  int N, double ipsi,   int B,    gsl_rng * r  )
{
   // Dirac-coalescent as Lambda-coalescent
  int biters = 0 ;
  int n,   fjoldi_blokka, currentn  ; 

  double t = 0. ;

  gsl_vector * vbi = gsl_vector_calloc( N + 2 ) ; 
  gsl_vector * svbi = gsl_vector_calloc( N + 2 ) ; 
  gsl_vector_int * ni = gsl_vector_int_calloc( N + 2 ) ;
  int * Umrodun  =  (int *)calloc( N, sizeof( int ) ) ;
  gsl_vector * lambdan = gsl_vector_calloc( N + 2 ) ; 
  int * tempni  = (int *)calloc( N + 2,  sizeof( int ) ) ;  
  int * btveir  = (int *)calloc( 2,  sizeof( int ) ) ;  
  
  double * P = (double *)calloc( 5, sizeof(double) ) ; 
  P[1] = P[2] = P[3] = P[4]  =  0.25  ;
  P[0] = 0. ;
  gsl_ran_discrete_t * dice =  gsl_ran_discrete_preproc( 5,  P ) ; 
  
  /* **************  XifourDiracsim   */
  // double  **alloc_2d_array(  int n)
  double ** MBetalambda  =  alloc_2d_array( N ) ; 
  Umrodun[0]  = 1  ; 

  for( n = 2 ; n <= N ; n++ ){
    for( biters = 1 ; biters < n ; biters++){
      // biters is number of active blocks after jump from n blocks
      MBetalambda[n][biters]  =   gsl_sf_choose(n,n - biters + 1) * gsl_sf_pow_int( ipsi, n - biters - 1)*gsl_sf_pow_int(1. - ipsi, biters - 1) ;  }
    Umrodun[ n-1] =   n ; }

  /* ****************** XifourDiracsimx  */
  // compute total coal rates
  // void totalrate( int N,  double cstikar[],  gsl_vector * vln )
  double stikar[2] =  { 2. ,  ipsi } ;
  totalrate( N, stikar, lambdan ) ;


  /*
  for( n = 2 ; n <= N ; n++ ){
       printf("%d  %g\n", n,  gsl_vector_get( lambdan, n) ) ; }
  */

  biters = 0 ; 
  while( biters < B ){

    // iterate a genealogy
    // collect branch lengths

    // the value of ni[i] is size of block i
    for( n = 1 ; n <= N ; n++ ){
    gsl_vector_int_set( ni, n, 1 ) ; 
    gsl_vector_set( vbi, n, 0. ) ; }
    // value of ni[0] is current number of active blocks
    gsl_vector_int_set( ni, 0, N ) ; 

    while( gsl_vector_int_get( ni, 0)  >  1 ){
      // ni[0] is current number of blocks
      

      /********
        printf("%d : ", gsl_vector_int_get( ni, 0) ) ;
       for( n = 1 ; n <= N ; n++ ){
	printf("%d ", gsl_vector_int_get( ni, n) ) ; }
      printf("\n") ; 
      **********/
 
      // draw time
      t = gsl_ran_exponential( r,  1. / gsl_vector_get( lambdan,  gsl_vector_int_get( ni, 0)  ) ) ;
      // printf("t %g\n", t ) ;
      
      // update branch lengths
      for( n = 1 ; n <= N ; n++ ){
	gsl_vector_set( vbi, gsl_vector_int_get(ni, n),  gsl_vector_get( vbi, gsl_vector_int_get(ni, n))  +  ( gsl_vector_int_get( ni, n) > 0 ?  t : 0. ) ) ; }
      
       /* ****************** XifourDiracsimx  */

      // draw total number of blocks  to  possibly merge in up to 4 groups
      // int  draw_number_blocks(  int b,  double * M, gsl_rng * r )
      currentn =  gsl_vector_int_get( ni, 0 ) ; 
      while( currentn == gsl_vector_int_get( ni, 0 ) ){
	fjoldi_blokka =  draw_number_blocks( gsl_vector_int_get( ni, 0), MBetalambda[ gsl_vector_int_get( ni, 0) ], r ) ; 
      
      //  printf("%d :  fjb %d \n", gsl_vector_int_get( ni, 0),    fjoldi_blokka ) ;
      
      // void  Beta_blocks_merge(  int N,  double ialpha,   int b,  gsl_vector_int * Ni,   gsl_ran_discrete_t * fjorir,  int * umrodun,    gsl_rng * r )
	Beta_blocks_merge( N, fjoldi_blokka,  ni,  dice,  Umrodun, r) ; } }
    // one realisation of branch lengths obtained
  
    // update the estimate of branch lengths
    for( n = 1 ; n <= N ; n++ ){
      gsl_vector_set( svbi, n,  gsl_vector_get( svbi, n)  +   gsl_vector_get( vbi, n) ) ; }

    // add to count of iterations
    biters =  biters  +  1 ; }


  for( n = 1 ; n < N ; n++ ){
    printf("%g\n", gsl_vector_get( svbi, n)/( (double)B ) ) ; }


  // clean up
  gsl_vector_free( lambdan ) ; 
  gsl_vector_free( vbi ) ;
  gsl_vector_free( svbi ) ; 
  gsl_vector_int_free( ni ) ; 
  free( Umrodun ) ; 
  
  free( tempni ) ; 
  free( btveir ) ;
  free( P ) ;
  gsl_ran_discrete_free( dice ) ; 
  // void free_2d_array( double **a, int n )
  free_2d_array( MBetalambda, N ) ; 


}

