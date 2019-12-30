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

#include <string> 
#include <cstring>
#include <cctype>
#include <fstream>
#include "model.h"

using namespace std;

bool is_token(const char *str){
 if(isalpha(str[0])) return true;
 else return false; 
}

bool is_number(const char *str){
 if(isdigit(str[0])) return true;
 else return false; 
}

int  Model::GetToken(const char *str){
 int result=T_NEXT;
 if(strcmp(str,"SIMULATIONS")==0) result=T_SIMULATIONS;
 if(strcmp(str,"PMIN")==0) result=T_PMIN;
 if(strcmp(str,"PMAX")==0) result=T_PMAX;
 if(strcmp(str,"MAXGENERATIONS")==0) result=T_MAXGENERATIONS;
 if(strcmp(str,"LEFTLOCI")==0) result=T_LEFTLOCI;
 if(strcmp(str,"RIGHTLOCI")==0) result=T_RIGHTLOCI;
 if(strcmp(str,"POPSIZE")==0) result=T_POPSIZE;
 if(strcmp(str,"MAXPOPSIZE")==0) result=T_MAXPOPSIZE;
 if(strcmp(str,"GROWTHRATE")==0) result=T_GROWTHRATE; 
 if(strcmp(str,"SEED")==0) result=T_SEED;
 if(strcmp(str,"RIGHTMUTRATE")==0) result=T_RIGHTMUTRATE;
 if(strcmp(str,"LEFTMUTRATE")==0) result=T_LEFTMUTRATE;   
 if(strcmp(str,"RIGHTRECRATE")==0) result=T_RIGHTRECRATE; 
 if(strcmp(str,"LEFTRECRATE")==0) result=T_LEFTRECRATE; 
 if(strcmp(str,"SMIN")==0) result=T_SMIN;
 if(strcmp(str,"NS")==0) result=T_NS;
 if(strcmp(str,"SINIT")==0) result=T_SINIT;
 if(strcmp(str,"VERBOSE")==0) result=T_VERBOSE;
 if(strcmp(str,"SELECTMIN")==0) result=T_SELECTMIN;
 if(strcmp(str,"SELECTMAX")==0) result=T_SELECTMAX;
 if(strcmp(str,"SELECTMODEL")==0) result=T_SELECTMODEL;
 if(strcmp(str,"GAMMA")==0) result=T_GAMMA;
 if(strcmp(str,"LOGERRORS")==0) result=T_LOGERRORS;
 return result;
}

void  Model::OutData(){
 cerr << "LEFTLOCI      : " << (int)nleft << endl; 
 cerr << "RIGHTLOCI     : " << (int)nright << endl; 
 cerr << "SIMULATIONS   : " << simulations << endl; 
 cerr << "PMIN:         : " << pmin  << endl;
 cerr << "PMAX          : " << pmax << endl;
 cerr << "SMIN          : " << smin << endl;
 cerr << "NS           : " << ns << endl;
 cerr << "SINIT         : " << sinit << endl;
 cerr << "SELECTMIN     : " << selectmin << endl;
 cerr << "SELECTMAX     : " << selectmax << endl;
 cerr << "SELECTMODEL   : " << (selectmodel?"Additive":"Dominant") << endl;
 cerr << "MAXGENERATIONS: " << maxgenerations << endl;
 cerr << "MAXPOPSIZE    : " << maxpopsize << endl;
 cerr << "GROWTHRATE    : " << growthrate << endl;
 cerr << "POPSIZE       : " << popsize << endl;
 cerr << "GAMMA         : " << gamma << endl;
 cerr << "RIGHTMUTRATE  : " << endl;
 for(int j=0;j<nright;j++) cerr << rmutrate[j] << " "; 
 cerr << endl;
 
 cerr << "LEFTMUTRATE   : " << endl;
 for(int j=0;j<nleft;j++) cerr << lmutrate[j] << " "; 
 cerr << endl;

 cerr << "RIGHTRECRATE  : " << endl;
 for(int j=0;j<nright;j++) cerr << rrecfrac[j] << " "; 
 cerr << endl;

 cerr << "LEFTRECRATE   : " << endl;
 for(int j=0;j<nleft;j++) cerr << lrecfrac[j] << " "; 
 cerr << endl;
}
 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                Single inpute parameters                                  //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool Model::GetLeftLoci(){
 f >> token;
 if(is_number(token)){
  nleft=(char)atoi(token);
  #ifdef DEBUG_MODEL
  cerr << "LEFTLOCI: " << (int) nleft << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;
 return true;
}

bool Model::GetRightLoci(){
 f >> token;
 if(is_number(token)){
  nright=(char)atoi(token);
  #ifdef DEBUG_MODEL
  cerr << "RIGHTLOCI: " << (int) nright << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;
 return true;
}


bool Model::GetSimulations(){
 f >> token;
 if(is_number(token)){
  simulations=(long)atof(token);
  #ifdef DEBUG_MODEL
  cerr << "SIMULATIONS: " << simulations << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;   
 return true;
}
   
bool Model::GetPmin(){
 f >> token;
 if(is_number(token)){
  pmin=(double)atof(token);
  #ifdef DEBUG_MODEL
  cerr << "PMIN: " << pmin << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;   
 return true;
}

bool Model::GetPmax(){
 f >> token;
 if(is_number(token)){
  pmax=(double)atof(token);
  #ifdef DEBUG_MODEL
  cerr << "PMAX: " << pmax << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;   
 return true;
}

bool Model::GetSmin(){
 f >> token;
 if(is_number(token)){
  smin=(double)atof(token);
  #ifdef DEBUG_MODEL
  cerr << "SMIN: " << smin << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;   
 return true;
}

bool Model::GetNs(){
 f >> token;
 if(is_number(token)){
  ns=atoi(token);
  #ifdef DEBUG_MODEL
  cerr << "NS: " << ns << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;   
 return true;
}

bool Model::GetSInit(){
 f >> token;
 if(is_number(token)){
  sinit=(long)atoi(token);
  #ifdef DEBUG_MODEL
  cerr << "SINIT: " << sinit << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;   
 return true;
}

bool Model::GetSelectmin(){
 f >> token;
 if(is_number(token)){
  selectmin=(double)atof(token);
  #ifdef DEBUG_MODEL
  cerr << "SELECTMIN: " << selectmin << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;   
 return true;
}

bool Model::GetSelectmax(){
 f >> token;
 if(is_number(token)){
  selectmax=(double)atof(token);
  #ifdef DEBUG_MODEL
  cerr << "SELECTMAX: " << selectmax << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;   
 return true;
}

bool Model::GetSelectmodel(){
 f >> token;
 if(is_number(token)){
  selectmodel=(int)atof(token);
  #ifdef DEBUG_MODEL
  cerr << "SELECTMODEL: " << (selectmodel?"Additive":"Dominant") << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;   
 return true;
}

bool Model::GetMaxGenerations(){
 f >> token;
 if(is_number(token)){
  maxgenerations=(int)atof(token);
  #ifdef DEBUG_MODEL
  cerr << "MAXGENERATIONS: " << maxgenerations << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;   
 return true;
}



bool Model::GetSeed(){
 f >> token;
 if(is_number(token)){
  randreset((long)atof(token));
  #ifdef DEBUG_MODEL
  cerr << "SEED: " << (long)atof(token) << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;   
 return true;
}

bool Model::SetVerbose(){
 verbose=true; 
 current_parameter=T_NEXT;
 return verbose;
}

bool Model::GetGamma(){
 f >> token;
 if(is_number(token)){
  gamma=(long)atoi(token);
  #ifdef DEBUG_MODEL
  cerr << "GAMMA: " << gamma << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false;   
 return true;
}

bool Model::SetLogErrors(){
 logerrors=true;
 current_parameter=T_NEXT;
 return logerrors;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                Multiple inpute parameters                                //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool Model::GetMaxPopSize(){
 f >> token;
 if(is_number(token)){
  maxpopsize=(long)atof(token);
  #ifdef DEBUG_MODEL
  cerr << "MAXPOPSIZE: " << maxpopsize << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false; 
 return true;
}

bool Model::GetGrowthRate(){
 f >> token;
 if(is_number(token)){
  growthrate=(double)atof(token);
  #ifdef DEBUG_MODEL
  cerr << "GROWTHRATE: " << growthrate << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false; 
 return true;
}

bool Model::GetPopSize(){
 f >> token;
 if(is_number(token)){
  popsize=(long)atof(token);
  #ifdef DEBUG_MODEL
  cerr << "POPSIZE: " << popsize << endl;
  #endif
  current_parameter=T_NEXT;
 }
 else return false; 
 return true;
}

bool Model::GetRightMutRate(){
 int i=0;
 do{
  f >> token;
  if(is_number(token)) rmutrate[i]=(double)atof(token);
  else break;
  i++;
 }while((!f.eof())&&(i<nright));
 if(i<nright){
  cerr << "Problems in RIGHTMUTRATE! Missing values or non-numeric characters in vector" << endl;
  cerr << i << " values read, but should have been " << nright << endl;
  return false;
 }
 #ifdef DEBUG_MODEL
 cerr << i << " values read for RIGHTMUTRATE" << endl;
 #endif
 current_parameter=T_NEXT;
 return true;
}

bool Model::GetLeftMutRate(){
 int i=0;
 do{
  f >> token;
  if(is_number(token)) lmutrate[i]=(double)atof(token);
  else break;
  i++;
 }while((!f.eof())&&(i<nleft));
 if(i<nleft){
  cerr << "Problems in LEFTMUTRATE! Missing values or non-numeric characters in vector" << endl;
  cerr << i << " values read, but should have been " << nleft << endl;
  return false;
 }
 #ifdef DEBUG_MODEL
 cerr << i << " values read for LEFTMUTRATE" << endl;
 #endif
 current_parameter=T_NEXT;
 return true;
}


bool Model::GetRightRecRate(){
 int i=0;
 do{
  f >> token;
  if(is_number(token)) rrecfrac[i]=(double)atof(token);
  else break;
  i++;
 }while((!f.eof())&&(i<nright));
 if(i<nright){
  cerr << "Problems in RIGHTRECRATE! Missing values or non-numeric characters in vector" << endl;
  cerr << i << " values read, but should have been " << nright << endl;
  return false;
 }
 #ifdef DEBUG_MODEL
 cerr << i << " values read for RIGHTRECRATE" << endl;
 #endif
 current_parameter=T_NEXT;
 return true;
}

bool Model::GetLeftRecRate(){
 int i=0;
 do{
  f >> token;
  if(is_number(token)) lrecfrac[i]=(double)atof(token);
  else break;
  i++;
 }while((!f.eof())&&(i<nleft));
 if(i<nleft){
  cerr << "Problems in LEFTRECRATE! Missing values or non-numeric characters in vector" << endl;
  cerr << i << " values read, but should have been " << nleft << endl;
  return false;
 }
 #ifdef DEBUG_MODEL
 cerr << i << " values read for LEFTRECRATE" << endl;
 #endif
 current_parameter=T_NEXT;
 return true;
}
 
bool Model::GetValues(int param){
 bool result=true;
 switch(param){
  case  1: result=GetLeftLoci(); break;
  case  2: result=GetRightLoci(); break;
  case  3: result=GetSelectmodel(); break; 
  case  4: result=GetSimulations(); break;
  case  5: result=GetPmin(); break;
  case  6: result=GetPmax(); break;
  case  7: result=GetMaxGenerations(); break;  
  case  8: result=SetLogErrors(); break;
  case  9: result=GetGamma(); break; 
  case 10: result=GetMaxPopSize(); break;
  case 11: result=GetGrowthRate(); break;
  case 12: result=GetSeed(); break;
  case 13: result=GetPopSize(); break;
  case 14: result=GetRightMutRate(); break;
  case 15: result=GetLeftMutRate(); break;
  case 16: result=GetRightRecRate(); break;
  case 17: result=GetLeftRecRate(); break;
  case 18: result=GetSmin(); break; 
  case 19: result=GetNs(); break; 
  case 20: result=GetSInit(); break; 
  case 21: result=SetVerbose(); break; 
  case 22: result=GetSelectmin(); break; 
  case 23: result=GetSelectmax(); break; 
  
 } 
 return result;
}



bool Model::ReadConf(string filename){   
 bool result=true;
 left_loci_defined=right_loci_defined=false;
  
 f.open(filename.c_str());
 if(f){
  do{
   f >> token;
   if(current_parameter==T_NEXT){ //A parameter and its values have been read sucessfully. This should be a new parameter name or garbage
    if(is_token(token)){
     current_parameter=GetToken(token);
     result=GetValues(current_parameter);
    } 
    //#ifdef DEBUG_MODEL
    else cerr << "Unidentified garbage: " << token << " (but proceeding anyway...)" << endl;
    //#endif
   }
   //#ifdef DEBUG_MODEL
   else cerr << "Unidentified garbage: " << token << " (but proceeding anyway...)" << endl;
   //#endif
  }while((!f.eof())&&(result==true));
  f.close();
 }
 else cerr << "Cannot open parameter file " << filename.c_str() << endl; 
 //cerr << "Now outputing data" << endl;
 if(verbose) OutData();
 return result;
} 

