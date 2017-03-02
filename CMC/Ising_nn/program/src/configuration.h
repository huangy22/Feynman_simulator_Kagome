//
//  configuration.h
//  Feynman_Simulator
//
//  Created by yuan on 4/14/15.
//  Copyright (c) 2015 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__configuration__
#define __Feynman_Simulator__configuration__

#include <stdio.h>
#include "utility/vector.h"
#include "lattice.h"

using namespace lattice;

namespace configuration{
    class Config {
    public:
        Vec<real> Spin[MxVol]={0.0};
        Config(const Lattice& lattice);
    };
}

#endif /* defined(__Feynman_Simulator__configuration__) */
