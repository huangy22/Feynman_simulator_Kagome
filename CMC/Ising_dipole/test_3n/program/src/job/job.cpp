//
//  job.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/1/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "job.h"
#include "utility/dictionary.h"
#include <iostream>

using namespace std;

para::Job::Job(string inputfile)
{

    Dictionary _Para;
    _Para.Load(inputfile);
    _Para = _Para.Get<Dictionary>("Job");

    GET(_Para, PID);
    GET(_Para, Sample);
    GET(_Para, Beta);
    GET(_Para, ExternalField);
    
    GET(_Para, LatticeName);
    GET(_Para, Dim);
    GET(_Para, L);
    auto Lnew = _Para.Get<std::vector<int> >("L");
    ASSERT_ALLWAYS(Dim == Lnew.size(), "Lattice dimension is " << Dim << ", not " << Lnew.size());
    for(int i=0; i<Dim; i++)
        L[i] = Lnew.data()[i];
    
    GET(_Para, NSublat);
    GET(_Para, Toss);
    GET(_Para, Sweep);
    GET_WITH_DEFAULT(_Para, Counter, 0);
    GET_WITH_DEFAULT(_Para, Seed, 0);
    if (_Para.HasKey("RNG"))
        GET(_Para, RNG);
    else
        RNG.Reset(Seed);
    auto _timer = _Para.Get<Dictionary>("Timer");
    GET(_timer, PrinterTimer);
    GET(_timer, DiskWriterTimer);
    
    string Prefix = ToString(PID)+"_MC";
    ParaFile = Prefix + "_para";
    StatisticsFile = Prefix + "_statis";
    LogFile = Prefix + ".log";
    InputFile = inputfile;
}
