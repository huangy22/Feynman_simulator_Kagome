//
//  job.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/1/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__job__
#define __Feynman_Simulator__job__

#include "utility/convention.h"
#include "utility/vector.h"
#include "utility/rng.h"
#include <string>
#include <set>

namespace para {
class Job {
public:
    Job(std::string inputfile);
    Job(bool, bool, int);

    int Sample;
    int PID;
    
    real Beta;
    real ExternalField;
    int Dim;
    std::vector<int> L;
    int NSublat;
    std::string LatticeName;
    
    long long Counter;
    int Toss;
    int Sweep;
    int Seed;
    RandomFactory RNG;
    
    int PrinterTimer;
    int DiskWriterTimer;
    
    std::string ParaFile;
    std::string StatisticsFile;
    std::string LogFile;
    std::string InputFile;
};
}
#endif /* defined(__Feynman_Simulator__job__) */
