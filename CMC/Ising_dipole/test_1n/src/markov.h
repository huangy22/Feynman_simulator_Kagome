//
//  markov.h
//  Feynman_Simulator
//
//  Created by yuan on 4/15/15.
//  Copyright (c) 2015 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__markov__
#define __Feynman_Simulator__markov__

#include <stdio.h>
#include "lattice.h"
#include "configuration.h"
#include "job/job.h"

using namespace para;
using namespace lattice;
using namespace configuration;

namespace markov{
    class Markov{
    public:
        long long* Counter;
        real Beta;
        real ExternalField;
        Lattice* Lat;
        Config* Conf;
        RandomFactory* RNG;
        
        Markov(Job&, Lattice&, Config&);
        void Hop(int);
        void PrintDetailBalanceInfo();
        
    private:
        real _Accepted;
        real _Proposed;
        Vec<real> _RandomPickSpin();
        int _RandomPickSite();
        Vec<real> _SumSpin(int);
	real _GetSign(int);
    };
}

#endif /* defined(__Feynman_Simulator__markov__) */
