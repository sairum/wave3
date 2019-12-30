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

#include "model.h"
#include "random.h"
#include <string>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {
    
 Model *m;
 
 string filename, output;
 
 randreset();
 
 m = new Model();
 if(m){
  if(argc>2){
   filename=argv[1];
   output=argv[2];
   if(m->ReadConf(filename)){
    m->Run(output);
   }   
  }
  else{
   cerr << "wave3 version 06112011" << endl;
   cerr << "Run as: wave3 <parameter_file> <result_file>" << endl;  
  } 
  delete m;
 } 

 return 0;
}


