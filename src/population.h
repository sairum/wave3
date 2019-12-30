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
 ***************************************************************************/
 
#ifndef POPULATION_H
#define POPULATION_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "config.h"

using std::string;
using std::ofstream;

class Population{
  long 		*alleles;	// Pointer to a matrix that hold haplotypes (line=left, column=right)
  double        *freqs;		// Pointer to frequencies of haplotypes
  long 		*clear;		// Pointer to an area of zeros that has half the number of total haplotypes 
  unsigned char	*right;		// Pointer to array with right haplotype attributes
  unsigned char	*left;		// Pointer to array with left haplotype attributes
  char		*htype;		// Pointer to array with haplotype type: A or S 
  char          nright;
  char          nleft;
  long 		size;		// Number of haplotypes in population 
  double	rsize;		// (and its double version for growing pops)
  long		psize;		// Same as above but this doesn't change and thus mantain original size if reset is needed
  long		maxpsize;	// Carrying capacity
  double	growthrate;	// as the name says...
  int           cleft;		// Number of possible different allele combinations at left side (0-256)
  int           cright;         // Number of possible different allele combinations at left side (0-256)
  int           sstart;
  long		tot_haplo;	// Total number of haplotypes in population
  double	rmr[MAXLOCI];	// Right mutation rate per locus
  double	lmr[MAXLOCI];	// Left mutation rate per locus
  double	rrf[MAXLOCI];	// Right recombination fraction per locus
  double	lrf[MAXLOCI];	// Left recombination fraction per locus
  double	aafit;
  double 	asfit;
  double	ssfit;
  double  	convrate;
  
public:
  Population();
  Population(long n, char nr, char nl);  
  ~Population();
  
  void 		Reset();
  void 		AddSHaplotype(long);
  long		GetSize();
  void          SetSize(long p);
  void          SetGrowthRate(double,long);
  string        SetMuteRate(double *rmr, double *lmr, int gamma);
  void          SetRecFraction(double *rrf, double *lrf);
  void		SetFitness(double aa, double as, double ss);
  void          SetConvertionRate(double cr);
  bool 		FixedS();
  bool 		FixedA();
  void          Mutate();
  void		Recombine();
  void          Convert();
  long		Mate();
  double	ComputeP();
  long		TotalS();
  double	ComputeSFreq();
  string	StatisticsHeader();
  string	SummaryStatistics();
  string	DetailedStatistics();
};

#endif
