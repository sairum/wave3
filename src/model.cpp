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

#include <cstring>
#include <fstream>
#include "model.h"



using namespace std;


Model::Model(){
 p=NULL;

 popsize=POPSIZE;
 maxpopsize=MAXPOPSIZE;
 aafitness=AAFITNESS;
 asfitness=ASFITNESS;
 ssfitness=SSFITNESS;
 nleft=NLEFT;
 nright=NRIGHT;
 simulations=SIMULATIONS;
 tot_haplo=0;
 pmin=PMIN;
 pmax=PMAX;
 smin=SMIN;
 ns=NS;
 sinit=SINIT;
 selectmodel=SELECTMODEL;
 selectmin=SELECTMIN;
 selectmax=SELECTMAX;
 selcoef=SELCOEF;
 gamma=GAMMA;
 logerrors=false;
 
 verbose=false; 
 maxgenerations=MAXGENERATIONS;
 memset(rmutrate,0,sizeof(rmutrate));
 memset(lmutrate,0,sizeof(lmutrate));
 memset(rrecfrac,0,sizeof(rrecfrac));
 memset(lrecfrac,0,sizeof(lrecfrac));
}

Model::~Model(){
 if(p) delete p;
 if(p1) delete [] p1;
 if(p2) delete [] p2;
 if(list) delete [] list;
}


void Model::Run(string output){
 long     i,g,s;
 double   tsize;
 bool	  stop,found; 
 double   pdiff,sthres,po,sf;
 string   mutrates;
 
 of.open(output.c_str());
 if(logerrors){
  output+=".err";
  ef.open(output.c_str());
 }
 
 if(!of){
  cerr << "Problems creating output file! See if there is space available!" << endl;
  exit(1);
 }
 
 if((logerrors)&&(!ef)){
  cerr << "Problems creating error output file! See if there is space available!" << endl;
  exit(1);
 }
 
 pdiff=pmax-pmin; 
 if((smin==0)||(pdiff==0)){
  cerr << "SMIN not defined (must be greater than 0) or PMIN equal to PMAX (they must be different)" << endl;
  exit(1);
 }
 
 
 tot_haplo=2*(int)pow(2.0,(double)nleft)*(int)pow(2.0,(double)nright);

 p= new Population(popsize,nright,nleft); 
 p->SetGrowthRate(growthrate,maxpopsize);
 p->SetRecFraction(rrecfrac,lrecfrac);
 
 of << p->StatisticsHeader() << endl;
  
 if(logerrors){
  ef << "simulation;";
  for(int j=0;j<nleft;j++) ef <<  "lmr" << j+1 << ";";
  for(int j=0;j<nright;j++) ef << "rmr" << j+1 << ";";
  ef << "sel_coef;s_min;s_freq;g;pop_size" << endl;
 }

 i=0;
 s=1;
 do{
  
  if(logerrors) ef << s << ";";
  
  p->Reset();
  p->AddSHaplotype(sinit);
  mutrates=p->SetMuteRate(rmutrate,lmutrate,gamma);
  
  
  
  if(logerrors) ef << mutrates;
  //Get selection coefficient from a flat prior
  
  selcoef=uniform(selectmin,selectmax);
  
  if(ns>0){ //Sample smin from binomial, otherwise leave it as it is
   sthres=(double)binomial(ns,smin)/ns;
  }
  else sthres=smin;
  
  //Compute fitness;
  
  if(selectmodel==ADDITIVE){ //Additive model
   aafitness=1.0-selcoef;
   asfitness=1.0-selcoef/2.0;
   ssfitness=1.0;
  }
  else{	//Dominant model
   aafitness=1.0-selcoef;
   asfitness=1.0;
   ssfitness=1.0;  
  }
  p->SetFitness(aafitness,asfitness,ssfitness);
   
  g=0;
  stop=false;
  found=false;
  
  po=sf=0;
  
  tsize=p->GetSize();
  
  
  if(logerrors) ef << selcoef << ";" << sthres << ";";
  do{
   g++; 
   tsize=(double)p->Mate();  
   //cerr << sthres << endl;
   if(!p->FixedS()){
    p->Reset();
    p->AddSHaplotype(sinit);
    g=0;
   }
   else{    
    sf=p->ComputeSFreq();
    if(sf>sthres){
     po=p->ComputeP();
     if((po>=pmin)&&(po<=pmax)) found=true;
     else found=false;
     stop=true;
    } 
   }  
  }while((g<maxgenerations)&&(tsize>0.0)&&(!stop));
  
  
  if(found){
   i++;
   of << p->SummaryStatistics() << ";" << selcoef << ";" << g << ";" << sf << endl;
   if(logerrors) ef << sf << ";" << g << ";" << tsize << ";" << i << " simulation(s) successful" << endl; 
  }
  else if(logerrors) ef << sf << ";" << g << ";" << tsize << ";" << endl; 
  s++;
 }while(i<simulations); 
 of.close();
 if(logerrors) ef.close();
}


