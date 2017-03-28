//
//  lattice.h
//  Feynman_Simulator
//
//  Created by yuan on 4/10/15.
//  Copyright (c) 2015 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__lattice__
#define __Feynman_Simulator__lattice__

#include <stdio.h>
#include "utility/vector.h"
#include "utility/Array2D.h"
#include "utility/Array4D.h"

namespace lattice {

const int MxL = 100;
const int MxVol = 3200;
const int MxNeighNum = 3200;
const int MxSub = 4;
const int MxDim = 3;
    
class Lattice {
public:
    int Dimension;
    int Vol;       ///Lx*Ly*Lz*SublatVol
    int SublatVol;
    int NeighNum;
    std::vector<int> L;

    Array2D<int, MxVol, MxDim> SiteVectors;
    int SiteSublat[MxVol]={0};
    Array4D<real, MxSub, MxSub, MxL, MxL> Coupling;

    Lattice(const int Dimension=3, const std::vector<int>& size={4,4,4}, int NSublat = 2);
};
}

#endif /* defined(__Feynman_Simulator__lattice__) */
