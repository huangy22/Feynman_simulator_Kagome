//
//  measure.h
//  Feynman_Simulator
//
//  Created by yuan on 4/15/15.
//  Copyright (c) 2015 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__measure__
#define __Feynman_Simulator__measure__

#include <stdio.h>
#include "lattice.h"
#include "configuration.h"
#include "estimator/estimator.h"
#include "utility/dictionary.h"
#include "job/job.h"

using namespace para;
using namespace lattice;
using namespace configuration;

namespace measure{
    class Observables{
    public:
        Lattice* Lat;
        Config* Conf;
        real Beta;
        std::string LatName;
        int PID;
        
        EstimatorBundle<real> Quantities;
        EstimatorBundle<real> Correlation;
        
        Observables(Job&, Lattice&, Config&);
        
        void Measure();
        void AddStatistics();
        void Save();
    private:
        int _GetSiteDiff(int, int);
	real _GetSign(int);
    };
}

#endif /* defined(__Feynman_Simulator__measure__) */
