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

#ifndef RANDOM_H
#define RANDOM_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <cassert>

#define THRESHOLD  0.0000000001


void   _seedTable();
double _u();
bool   bernoulli(double);
double uniform( double , double );
void   randreset();
void   randreset(long);
int    binomial( int , double );
double normal(double, double );
unsigned int poisson( double  );
unsigned long multinomial(unsigned long, double *, unsigned long *, unsigned long  );
unsigned long multinomial_one(double *, unsigned long);
double erlang( double , int );

#endif
