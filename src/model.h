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
 
#ifndef MODEL_H
#define MODEL_H


#include "config.h"
#include "population.h"
#include "random.h"
#include <fstream>

using std::string;

class Model{
  Population    *p;
  
  char      nright;
  char      nleft;
  long      popsize;
  long      maxpopsize;
  double    growthrate;
  double    selcoef;
 
  long      simulations;
  double    aafitness;
  double    asfitness;
  double    ssfitness;
  
  double    rmutrate[MAXLOCI];
  double    lmutrate[MAXLOCI];
  
  double    rrecfrac[MAXLOCI];
  double    lrecfrac[MAXLOCI];

  
  long      sinit;
  int       *source1,*source2;
  long      *list;
  
  long      tot_haplo;
  long      *p1,*p2;
  double    pmin,pmax,smin,selectmin,selectmax;
  int       selectmodel,ns;
  long      maxgenerations;
  
  std::ifstream f;
  std::ofstream of,ef;
  
  char  token[100];
  int   current_parameter;
  bool  left_loci_defined, right_loci_defined, verbose;
  bool  logerrors;
  int   gamma;
  
  bool GetLeftLoci();
  bool GetRightLoci();
  bool GetSimulations();
  bool GetPmin();
  bool GetPmax();
  bool GetSmin();
  bool GetNs();
  bool GetSInit();
  bool GetSelectmin();
  bool GetSelectmax();
  bool GetSelectmodel();
  bool GetMaxGenerations();
  bool GetGametes();
  bool GetMaxPopSize();
  bool GetGrowthRate();
  bool GetSeed();
  bool GetPopSize();
  bool GetRightMutRate();
  bool GetLeftMutRate();
  bool GetRightRecRate();
  bool GetLeftRecRate();
  bool GetGamma();
  bool SetLogErrors();
  bool SetVerbose();
    
  bool GetValues(int);
  int  GetToken(const char *);
  void OutData();
    
public:
  Model(); 
  ~Model();
  

  bool ReadConf(string filename);
  
  void Run(string);
};

#endif
