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
 
#ifndef CONFIG_H
#define CONFIG_H

//#define DEBUG_MODEL

#define ADDITIVE        1
#define DOMINANT        0

#define MAXLOCI         8
#define MAXPOPSIZE      1000
#define	POPSIZE         1000
#define NLEFT           8
#define NRIGHT          8
#define SIMULATIONS     1
#define	PMIN            0
#define PMAX            0
#define SMIN            0
#define NS              0
#define MAXGENERATIONS  1000
#define	SINIT           1
#define SELECTMIN       0
#define SELECTMAX       1
#define SELECTMODEL     0
#define AAFITNESS       1
#define SSFITNESS       1
#define ASFITNESS       1
#define SELCOEF         0
#define GAMMA           0
    
#define T_NEXT              0
#define T_LEFTLOCI          1
#define T_RIGHTLOCI         2
#define T_SELECTMODEL       3
#define T_SIMULATIONS       4
#define T_PMIN              5
#define T_PMAX              6
#define T_MAXGENERATIONS    7
#define T_LOGERRORS         8
#define T_GAMMA             9
#define T_MAXPOPSIZE        10
#define T_GROWTHRATE        11 
#define T_SEED              12
#define T_POPSIZE           13
#define T_RIGHTMUTRATE      14
#define T_LEFTMUTRATE       15   
#define T_RIGHTRECRATE      16 
#define T_LEFTRECRATE       17
#define T_SMIN              18
#define T_NS                19
#define T_SINIT             20
#define T_VERBOSE           21
#define T_SELECTMIN         22
#define T_SELECTMAX         23


#endif

