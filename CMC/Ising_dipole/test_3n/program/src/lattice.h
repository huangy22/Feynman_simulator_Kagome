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

namespace lattice {

const int MxVol=1000;
const int MxNeighNum = 100;
    
class Lattice {
public:
    int Dimension;
    int Vol;       ///Lx*Ly*Lz*SublatVol
    int SublatVol;
    int NeighNum;
    std::vector<int> L;

    int Neighbor[MxVol][MxNeighNum];
    real Coupling[MxVol][MxNeighNum];

    int TetrahedraNum;
    int Tetrahedra[MxVol][MxNeighNum];

    Lattice(const int Dimension=3, const std::vector<int>& size={4,4,4}, int NSublat = 2);
};
}

#endif /* defined(__Feynman_Simulator__lattice__) */
