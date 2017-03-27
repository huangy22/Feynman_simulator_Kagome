//
//  configuration.cpp
//  Feynman_Simulator
//
//  Created by yuan on 4/14/15.
//  Copyright (c) 2015 Kun Chen. All rights reserved.
//

#include "configuration.h"
#include "lattice.h"
#include "utility/convention.h"
using namespace lattice;
using namespace configuration;

Config::Config(const Lattice& lattice)
{
    for(int i=0; i<lattice.Vol; i++){
        Spin[i][0]=1.0;
        for(int j=1; j<D; j++)
            Spin[i][j]=0.0;
    }
}
