/***************************************************************************
 *   Copyright (C) 2011 by António Múrias dos Santos                       *
 *   amsantos@fc.up.pt                                                     *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 *                                                                         *
 *   IMPORTANT NOTICE!                                                     *
 *   The code here presented was derived from Richard Saucier's Random.h   *
 *   The original file can be found in http://ftp.arl.mil/random/Random.h  *
 *                                                                         *
 *   Here is its copyright notice                                          *
 *   Random.h: Definition and Implementation of Random Number Distribution *
 *   Class. Ref: Richard Saucier, "Computer Generation of Statistical      *
 *   Distributions," ARL-TR-2168, US Army Research Laboratory,             *
 *   Aberdeen Proving Ground, MD, 21005-5068, March 2000.                  * 
 *                                                                         *
 ***************************************************************************/
 
#include "random.h"

static const long   _M    = 0x7fffffff; // 2147483647 (Mersenne prime 2^31-1)
static const long   A_    = 0x10ff5;	 // 69621
static const long   Q_    = 0x787d;	 // 30845
static const long   R_    = 0x5d5e;	 // 23902
static const double _F    = 1. / _M;
static const short  _NTAB = 32;	 // arbitrary length of shuffle table
static const long   _DIV  = 1+(_M-1)/_NTAB;
long	      _table[ _NTAB ];  	// shuffle table of seeds
long	      _next;			// seed to be used as index into table
long	      _seed;			// current random number seed


void randreset(){   // reset the seed 
 _seed = (long) time(NULL);
 _seedTable();
}

void randreset(long s){   // reset the seed 
 _seed = s;
 _seedTable();
}

void _seedTable(){                   // seeds the shuffle table
  for ( int i = _NTAB + 7; i >= 0; i-- ) {     // first perform 8 warm-ups

  long k = _seed / Q_;                       // seed = ( A * seed ) % M
  _seed = A_ * ( _seed - k * Q_ ) - k * R_;  // without overflow by
  if ( _seed < 0 ) _seed += _M;              // Schrage's method

  if ( i < _NTAB ) _table[ i ] = _seed;      // store seeds into table
 }
 _next = _table[ 0 ];                          // used as index next time
}

double _u(){                             // uniform rng
 long k = _seed / Q_;                          // seed = ( A*seed ) % M
 _seed = A_ * ( _seed - k * Q_ ) - k * R_;     // without overflow by
 if ( _seed < 0 ) _seed += _M;                 // Schrage's method.

 int index = _next / _DIV;                     // Bays-Durham shuffle
 _next = _table[ index ];                      // seed used for next time
 _table[ index ] = _seed;                      // replace with new seed. 

 return _next * _F;                            // scale value within [0,1)
}
   

bool bernoulli( double p = 0.5 ){   // Bernoulli Trial
 assert( 0. <= p && p <= 1. );  
 return _u() < p;
}

double uniform( double xMin = 0., double xMax = 1. ){   // Uniform on [xMin,xMax)
 //assert( xMin < xMax );  
 return xMin + ( xMax - xMin ) * _u();
}

int binomial( int n, double p ){   // Binomial
 //assert( n >= 1 && 0. <= p && p <= 1. );   
 int sum = 0;
 for ( int i = 0; i < n; i++ ) sum += bernoulli( p );
 return sum;
}

unsigned int poisson( double mu )   // Poisson
{
 double a = exp( -mu );
 double b = 1.;
 int i;
 for ( i = 0; b >= a; i++ ) b *= _u();   
 return i - 1;
}

double normal( double mu = 0., double sigma = 1. ){   // Normal 
 static bool f = true;
 static double p, p1, p2;
 double q;  
 if ( f ) {
  do {
   p1 = uniform( -1., 1. );
   p2 = uniform( -1., 1. );
   p = p1 * p1 + p2 * p2;
  }while( p >= 1. );
  q = p1;
 }
 else q = p2;
 f = !f;
 return mu + sigma * q * sqrt( -2. * log( p ) / p );
}

/*
// Multinomial trials n, probability vector p, success vector count, number of disjoint events m
unsigned int multinomial(unsigned int n, double *p, unsigned int *count, int m ){
 unsigned int sum = 0;
 if(m>1){
  memset(count,0,m*sizeof(unsigned int));  // initialize
  // generate n uniform variates in the interval [0,1) and bin the results
  sum = 0;
  for (unsigned int i = 0; i < n; i++ ){
   double lower = 0., upper = 0., u = _u();
   for ( int bin = 0; bin < m; bin++ ) {
    // locate subinterval, which is of length p[ bin ],
    // that contains the variate and increment the corresponding counter
    lower = upper;
    upper += p[ bin ];
    if(bin==(m-1)) { count[ bin ]++; sum++; break; }
    if(lower <= u && u < upper) { count[ bin ]++; sum++; break; }
   }
  }
 }
 return sum;
}
*/

unsigned long multinomial(unsigned long  n, double *p, unsigned long *c, unsigned long m ){
 unsigned long i,bin,sum;
 double lower,upper,u;
 // c (where counts of classes are stored) must be set to zero outside this function
 // if not, each bin will accumulate counts along several uses of multinomial. This
 // behaviour is intended in Mute() and Recombine() functions
 sum =0;
 if(m>1){ 
  for (i = 0; i < n; i++ ){
   lower = 0.0; upper = 0.0; u = _u();
   for (bin = 0; bin < m; bin++ ) {
    // locate subinterval, which is of length p[ bin ],
    // that contains the variate and increment the corresponding counter
    lower = upper;
    upper += p[ bin ];
    if ( lower <= u && u < upper ) { c[ bin ]++; sum++; break; }
   }
  }
 }
 return sum;
}

unsigned long multinomial_one(double *p, unsigned long m ){
 if(m>1){
  double lower = 0., upper = 0.0, u = _u();
  for (unsigned long bin = 0; bin < m; bin++ ) {
   lower = upper;
   upper += p[ bin ];
   if ( lower <= u && u < upper ) return bin;
  }
 }
 return 0;
}

double erlang( double b, int c ){   // Erlang (b > 0. and c >= 1)
 if( b > 0. && c >= 1 ){  
  double prod = 1.;
  for ( int i = 0; i < c; i++ ) prod *= _u();      
  return -b * log( prod );
 }
 else return 0;    
}
