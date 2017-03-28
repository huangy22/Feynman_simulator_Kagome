//
//  lattice.cpp
//  Feynman_Simulator
//
//  Created by yuan on 4/10/15.
//  Copyright (c) 2015 Kun Chen. All rights reserved.
//

#include "lattice.h"
#include "utility/utility.h"
using namespace lattice;

Lattice::Lattice(const int dim, const std::vector<int>& size, int NSublat)
{
    Dimension = dim;
    
    SublatVol = NSublat;
    L = size;
    NeighNum = 0;
    
    Vol = 1;
    for (int i = 0; i < Dimension; i++) {
     Vol *= L[i];
    }
    Vol *= SublatVol;

    for( int i=0; i<MxVol; i++){
	for(int j=0; j<MxDim; j++){
	    SiteVectors.at(i, j)= 0;
	}
    }

    for( int i=0; i<MxSub; i++)
	for(int j=0; j<MxSub; j++)
	    for(int x=0; x<MxL; x++)
		for(int y=0; y<MxL; y++)
		    Coupling.at(i,j,x,y) = 0.0;

}
