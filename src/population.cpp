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
 
#include <sstream> 
#include <cmath>
#include <cstring>
#include "random.h"
#include "population.h"
#include <iomanip>

using namespace std;

Population::Population(){
 alleles=NULL;
 freqs=NULL;
 right=NULL;
 left=NULL;
 clear=NULL;
 size=psize=maxpsize=0;
 rsize=(double)size;
 growthrate=0;
 cleft=0;
 cright=0;
 nright=0;
 nleft=0;
 sstart=0;
 tot_haplo=0;
 aafit=1;
 ssfit=1;
 asfit=1;
 convrate=0;
 memset(rmr,0,sizeof(rmr));
 memset(lmr,0,sizeof(lmr));
 memset(rrf,0,sizeof(rrf));
 memset(lrf,0,sizeof(lrf));
}

Population::Population(long n, char nr, char nl){
 long i,j;
 
 nright=nr;
 nleft=nl;
 cright=(int)pow(2.0,(double)nr);
 cleft=(int)pow(2.0,(double)nl);
 size=psize=maxpsize=n;
 rsize=(double)n;
 growthrate=0;
 sstart=cright*cleft;
 tot_haplo=2*cleft*cright;
 aafit=1;
 ssfit=1;
 asfit=1; 
 memset(rmr,0,sizeof(rmr));
 memset(lmr,0,sizeof(lmr));
 memset(rrf,0,sizeof(rrf));
 memset(lrf,0,sizeof(lrf));
 alleles  = new long[tot_haplo];
 if(alleles) memset(alleles,0,sizeof(long)*tot_haplo);
 freqs  = new double[tot_haplo];
 if(freqs) memset(freqs,0,sizeof(double)*tot_haplo);
 clear  = new long[sstart];
 if(clear) memset(clear,0,sizeof(long)*sstart);
 
 right = new unsigned char[tot_haplo];
 if(right){
  for(i=0;i<cleft;i++){
   for(j=0;j<cright;j++){
    right[i*cright+j]=(unsigned char)j;
    right[sstart+i*cright+j]=(unsigned char)j;
   } 
  } 
 }
 
 left = new unsigned char[tot_haplo];
 if(left){
  for(i=0;i<cleft;i++){
   for(j=0;j<cright;j++){
    left[i*cright+j]=(unsigned char)i;
    left[sstart+i*cright+j]=(unsigned char)i;
   } 
  } 
 }
 
 htype = new char[tot_haplo];
 if(htype){
  for(i=0;i<cleft;i++){
   for(j=0;j<cright;j++){
    htype[i*cright+j]='A';
    htype[sstart+i*cright+j]='S';
   } 
  } 
 }
}


Population::~Population(){
 if(alleles) delete [] alleles;
 if(freqs) delete [] freqs;
 if(right) delete [] right;
 if(left) delete [] left;
 if(htype) delete [] htype;
 if(clear) delete [] clear;
}

void Population::Reset(){
 long i;
 unsigned char l,r;
 if(alleles&&freqs){
  size=psize;
  rsize=(double)size;
  memset(alleles,0,sizeof(long)*tot_haplo);
  memset(freqs,0,sizeof(double)*tot_haplo);
  for(i=0;i<size;i++){
   l  = (unsigned char)uniform(0,cleft);  
   r  = (unsigned char)uniform(0,cright);
   alleles[l*cright+r]++;
  }
  if(size>0){
   for(i=0;i<tot_haplo;i++)
    freqs[i]=(double)alleles[i]/size;  
  } 
 } 
}

void Population::AddSHaplotype(long sn){
 long i,j;
 if(size>0){
  for(j=0;j<sn;j++){
   i=multinomial_one(freqs,tot_haplo);
   //cerr << i << "\t";
   alleles[tot_haplo-1]++;
   alleles[i]--;
   freqs[i]=(double)alleles[i]/(size-(j+1));
   //cerr << size-(j+1) << endl;
  }
  for(j=0;j<tot_haplo;j++) freqs[j]=(double)alleles[j]/size;
 }  
}

long Population::GetSize(){
 return size;
}

void Population::SetSize(long p){
 size=p;
 rsize=(double)size;
}

string Population::SetMuteRate(double *rmrate, double *lmrate, int gamma){
 ostringstream s (ostringstream::out);
 if(gamma==0){
  for(int i=0;i<nright;i++) rmr[i]=rmrate[i];
  for(int i=0;i<nleft;i++)  lmr[i]=lmrate[i];
 }
 else{
  for(int i=0;i<nright;i++) rmr[i]=erlang(rmrate[i]/gamma,gamma);
  for(int i=0;i<nleft;i++)  lmr[i]=erlang(lmrate[i]/gamma,gamma); 
 }
 //cerr << "lmutrate: ";
 for(int j=0;j<nleft;j++) s << setiosflags(ios::right|ios::fixed)  << lmr[j] <<  ";";
 //cerr << "rmutrate: ";
 for(int j=0;j<nright;j++) s << setiosflags(ios::right|ios::fixed) << rmr[j]  << ";";
 //cerr << endl;
 return s.str();
}

void Population::SetRecFraction(double *rrfrac, double *lrfrac){
 int i;
 for(i=0;i<nright;i++) rrf[i]=rrfrac[i];
 for(i=0;i<nleft;i++)  lrf[i]=lrfrac[i];
}

void Population::SetFitness(double aa, double as, double ss){
 aafit=aa;
 asfit=as;
 ssfit=ss;
}

void Population::SetConvertionRate(double cr){
 convrate=cr;
}

bool Population::FixedS(){
 long i;
 for(i=sstart;i<tot_haplo;i++) if(alleles[i]>0) return true;
 return false;
}

bool Population::FixedA(){
 long i;
 for(i=0;i<sstart;i++) if(alleles[i]>0) return true;
 return false;
}

void Population::SetGrowthRate(double gr, long mps){
 growthrate=gr;
 maxpsize=mps;
 if(maxpsize<size) maxpsize=size;
}
  
void Population::Mutate(){
 long i,n,h1,p;
 unsigned char j,k,l,r,m;
 char	       t;
 
 for(k=0;k<nleft;k++){
  n=(long)poisson(lmr[k]*size);
  for(i=0;i<n;i++){
   h1=multinomial_one(freqs,tot_haplo);
   if(alleles[h1]>0){ 
    m=0x01;
    for(j=0;j<k;j++) m=m<<1;
    l=left[h1];
    if(l&m){
     r=right[h1];
     t=htype[h1];
     l-=m;
     if(t=='A'){
      p=l*cright+r;
      alleles[p]++;  
      freqs[p]=(double)alleles[p]/size; 
     } 
     else{
      p=sstart+l*cright+r;
      alleles[p]++;     
      freqs[p]=(double)alleles[p]/size;
     } 
     alleles[h1]--;
     freqs[h1]=(double)alleles[h1]/size;
    }
   }     
  }
 }
    
 for(k=0;k<nright;k++){
  n=(long)poisson(rmr[k]*size);
  for(i=0;i<n;i++){
   h1=multinomial_one(freqs,tot_haplo);
   if(alleles[h1]>0){ 
    m=0x01;
    for(j=0;j<k;j++) m=m<<1;
    r=right[h1];
    if(r&m){ 
     l=left[h1];
     t=htype[h1];
     r-=m;
     if(t=='A'){
      p=l*cright+r;
      alleles[p]++;   
      freqs[p]=(double)alleles[p]/size;
     } 
     else{
      p=sstart+l*cright+r;
      alleles[p]++;
      freqs[p]=(double)alleles[p]/size;     
     } 
     alleles[h1]--;
     freqs[h1]=(double)alleles[h1]/size;
    }
   }      
  }
 }   
}

void Population::Recombine(){
 long n,q,h1,h2;
 unsigned char j,k,l1,r1,l2,r2,m,part;
 char	       t1,t2;

 for(k=0;k<nleft;k++){
  j=0xFF;
  for(m=0;m<k;m++) j=j<<1;
  m=255&j;
  n=(long)poisson(lrf[k]*size/2.0);
  q=0;
  if(n>0){
   do{
    h1=multinomial_one(freqs,tot_haplo);
    h2=multinomial_one(freqs,tot_haplo);   
    if((alleles[h1]>0)&&(alleles[h2]>0)){	   
     r1=right[h1];
     l1=left[h1];
     t1=htype[h1];
     r2=right[h2];
     l2=left[h2];
     t2=htype[h2];	 
     part=l1&m;
     alleles[h1]--;
     alleles[h2]--; 
     l1=(l1&(~j))|(l2&m);
     l2=(l2&(~j))|part;   
     if(t1=='A') alleles[l1*cright+r1]++;   
     else alleles[sstart+l1*cright+r1]++;     
     if(t2=='A') alleles[l2*cright+r2]++;   
     else alleles[sstart+l2*cright+r2]++; 
     freqs[h1]=(double)alleles[h1]/size;
     freqs[h2]=(double)alleles[h2]/size;
     if(t1=='A') freqs[l1*cright+r1]=(double)alleles[l1*cright+r1]/size;   
     else freqs[sstart+l1*cright+r1]=(double)alleles[sstart+l1*cright+r1]/size;     
     if(t2=='A') freqs[l2*cright+r2]=(double)alleles[l2*cright+r2]/size;   
     else freqs[sstart+l2*cright+r2]=(double)alleles[sstart+l2*cright+r2]/size;     
     q++;
    } 
   }while(q<n);
  } 
 }
 
 for(k=0;k<nright;k++){
  j=0xFF;
  for(m=0;m<k;m++) j=j<<1;
  m=255&j;
  n=(long)poisson(rrf[k]*size/2.0);
  q=0;
  if(n>0){
   do{
    h1=multinomial_one(freqs,tot_haplo);
    h2=multinomial_one(freqs,tot_haplo);
    if((alleles[h1]>0)&&(alleles[h2]>0)){   
     r1=right[h1];
     l1=left[h1];
     t1=htype[h1];
     r2=right[h2];
     l2=left[h2];
     t2=htype[h2];	  
     part=r1&m;
     alleles[h1]--;
     alleles[h2]--; 
     r1=(r1&(~j))|(r2&m);
     r2=(r2&(~j))|part;    
     if(t1=='A') alleles[l1*cright+r1]++;   
     else alleles[sstart+l1*cright+r1]++;    
     if(t2=='A') alleles[l2*cright+r2]++;   
     else alleles[sstart+l2*cright+r2]++;
     freqs[h1]=(double)alleles[h1]/size;
     freqs[h2]=(double)alleles[h2]/size;
     if(t1=='A') freqs[l1*cright+r1]=(double)alleles[l1*cright+r1]/size;   
     else freqs[sstart+l1*cright+r1]=(double)alleles[sstart+l1*cright+r1]/size;     
     if(t2=='A') freqs[l2*cright+r2]=(double)alleles[l2*cright+r2]/size;   
     else freqs[sstart+l2*cright+r2]=(double)alleles[sstart+l2*cright+r2]/size;
     q++; 
    }
   }while(q<n);
  } 
 } 
}

void Population::Convert(){
 long i,n,h1,h2;
 unsigned char l,r;
 if((convrate>0)&&(size>0)){
  if(FixedS()){
   n=(long)poisson(convrate*size/2);
   for(i=0; i<n;i++){
    h1=multinomial_one(freqs,tot_haplo);
    h2=multinomial_one(freqs,tot_haplo);
    if((alleles[h1]>0)&&(alleles[h2]>0)&&(htype[h1]!=htype[h1])){
     if(htype[h1]=='A'){ 
      alleles[h1]--;
      freqs[h1]=(double)alleles[h1]/size;
      r=right[h1];
      l=left[h1];
     }
     else{
      alleles[h2]--;
      freqs[h2]=(double)alleles[h2]/size;
      r=right[h2];
      l=left[h2];
     }
     alleles[sstart+l*cright+r]++;  
     freqs[sstart+l*cright+r]=(double)alleles[sstart+l*cright+r]/size; 
    }	   
   }
  } 
 }
}

long Population::Mate(){
 unsigned long h1,h2;
 long i,s;
 double af;
 bool survive;
 
 if(size>0){

  Mutate();
  
  Recombine();
   
  Convert();
   
  //cerr << rsize*(1.0+(growthrate*(maxpsize-rsize))/(double)maxpsize) << "\t";
  
  //memset(alleles,0,sizeof(long)*tot_haplo);
  for(i=0;i<tot_haplo;i++) alleles[i]=0;

  s=0;
  
  if(growthrate>0){
   do{
    h1=multinomial_one(freqs,tot_haplo);
    h2=multinomial_one(freqs,tot_haplo);
    survive=false;  
    af=uniform(0,1);
    if((htype[h1]=='A')&&(htype[h2]=='A')){ // This is an AA  
     if(af<aafit) survive=true;
     //else cerr << "*";
    }
    else{
     if(((htype[h1]=='A')&&(htype[h2]=='S'))||((htype[h1]=='S')&&(htype[h2]=='A'))){ // This is an AS
      if(af<asfit) survive=true;
     }
     else if(af<ssfit) survive=true;   // Only SS remains....
    }
    if(survive){
     alleles[h1]++; 
     alleles[h2]++; 
     s+=2;
    }
   }while(s<size);    
  }
  else{
   do{
    h1=multinomial_one(freqs,tot_haplo);
    h2=multinomial_one(freqs,tot_haplo);
    survive=false;  
    af=uniform(0,1);
    if((htype[h1]=='A')&&(htype[h2]=='A')){ // This is an AA  
     if(af<aafit) survive=true;
    }
    else{
     if(((htype[h1]=='A')&&(htype[h2]=='S'))||((htype[h1]=='S')&&(htype[h2]=='A'))){ // This is an AS
      if(af<asfit) survive=true;
     }
     else if(af<ssfit) survive=true;   // Only SS remains....
    }
    if(survive){
     alleles[h1]++; 
     alleles[h2]++; 
     s+=2;
    }	
   }while(s<size);  
  } 
  
  if(growthrate>0){
   rsize=rsize*(1.0+(growthrate*((maxpsize-rsize)/(double)maxpsize)));
   size=(long)rsize;
  }
  else size=psize;	   
 
  // recompute allele frequencies  
  for(i=0;i<tot_haplo;i++) freqs[i]=(double)alleles[i]/size; 

  //cerr << size << "\t" << rsize << "\t" << growthrate << "\t" << maxpsize <<  "\t" << ComputeSFreq() << "\t";
 }
 return size;
}



double Population::ComputeP(){
 long totals=0;
 long original;
 int i;
 for(i=sstart;i<tot_haplo;i++){
  if(alleles[i]>0) totals+=alleles[i];
 }
 original=alleles[tot_haplo-1];
 if(totals>0) return (double) original/totals;
 else return 0;
}

long Population::TotalS(){
 long totals=0;
 int i;
 for(i=sstart;i<tot_haplo;i++) if(alleles[i]>0) totals+=alleles[i];
 return totals;
}

double Population::ComputeSFreq(){
 long totals=0;
 int i;
 for(i=sstart;i<tot_haplo;i++) if(alleles[i]>0) totals+=alleles[i];
 if(size>0) return (double)totals/size;
 return 0;
}

string Population::StatisticsHeader(){
 ostringstream s (ostringstream::out);
 s << "s_number;different_s;original_s;heter;p;size;sel_coef;g;s_freq";
 return s.str(); 
}


  
string Population::SummaryStatistics(){
 long i,totals,diffs;
 unsigned char j;
 double het,p;
 
 ostringstream s (ostringstream::out);
 
 s << "";
 
 totals=0;
 
 //memset(ls,0,sizeof(ls));
 //memset(rs,0,sizeof(rs));


 for(i=sstart;i<tot_haplo;i++){
  if(alleles[i]>0){
   totals+=alleles[i];
  }
 }
 //cout << (int)nleft << " " << (int) nright << endl;
 
 
 diffs=0;
 if(totals>0){ 	//Compute number of different S haplotypes  
  for(i=0;i<cleft;i++)
   for(j=0;j<cright;j++) if(alleles[sstart+i*cright+j]>0) diffs++;
 }
 
 // Compute S expected heterozigousity
 het=0;
 if(totals>0){
  for(i=0;i<cleft;i++)
   for(j=0;j<cright;j++) if(alleles[sstart+i*cright+j]>0) het+=pow((double)alleles[sstart+i*cright+j]/totals,2);
  het=1-het;
 }
 
 if(totals>0) p=(double)alleles[tot_haplo-1]/totals; 
 else p=0;
 s << totals << ";" << diffs << ";"  << alleles[tot_haplo-1] << ";" << het << ";";
 s << p << ";" << size;

 return s.str();  
}


string Population::DetailedStatistics(){
 long total,totals,diffs;
 unsigned char k,l1,l2;
 int i,j;
 ostringstream s (ostringstream::out);
 
 s << "";
 // Compute frequencies of different A haplotypes

 totals=0; 
 diffs=0;
 for(i=sstart;i<tot_haplo;i++){
  totals+=alleles[i];
  if(alleles[i]>0) diffs++;
 } 
 
 if(totals>0){
  s << "\t" << diffs << endl; 
  // Compute frequencies of different S haplotypes 
  for(i=0;i<cleft;i++){
   for(j=0;j<cright;j++){
    total=alleles[sstart+i*cright+j];
    
    if(total>0){
     
     l1=(unsigned char)i;
     l2=(unsigned char)j;   

     s << i << "\t" << j << "\t";

     for(k=0x80;k;k=k>>1) (k&l1)?s<<"1":s<<"0";
     s << "S";
     // Print right side...
     for(k=0x01;k;k=k<<1) (k&l2)?s<<"1":s<<"0";
     
     s << "\t" << total << "\t" << (double)total/totals << endl;
     
     //s << i << "S" << j << " :" << total << " (" << (double)total/size << ")" << " " << freqs[sstart+i*cright+j] << endl;
  
    }
   }
  }
 }
 else s << "\t" << 0 << endl; 
 return s.str();
}
