//
//  main.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

/********************** include files *****************************************/
#include <iostream>
#include <unistd.h>
#include "test.h"
#include "lattice.h"
#include "configuration.h"
#include "markov.h"
#include "measure.h"
#include "utility/pyglue/pywrapper.h"
#include "utility/pyglue/pyarraywrapper.h"
#include "job/job.h"
#include "utility/timer.h"
#include "utility/dictionary.h"
#include "utility/vector.h"

using namespace std;
using namespace para;
using namespace lattice;
using namespace configuration;
using namespace markov;
using namespace measure;

const string HelpStr = "Usage:"
                       "-p N / --PID N   use N to construct input file path."
                       "or -f / --file PATH   use PATH as the input file path.";

void MonteCarlo(Job&, Lattice&, Config&);
Lattice LoadLattice(const Job&);

int main(int argc, const char* argv[])
{
    Python::Initialize();
    Python::ArrayInitialize();
    RunTest();
    string InputFile;

    if (strcmp(argv[1], "-p") == 0 || strcmp(argv[1], "--PID") == 0)
        InputFile = string("infile/_in_MC_") + argv[2];
    else if (strcmp(argv[1], "-f") == 0 || strcmp(argv[1], "--file") == 0)
        InputFile = argv[2];
    else
        ABORT("Unable to parse arguments!\n" + HelpStr);

    para::Job Job(InputFile);

    LOGGER_CONF(Job.LogFile, "MC", Logger::file_on | Logger::screen_on, INFO, INFO);

    Lattice Lat;
    Lat = LoadLattice(Job);
    Config Conf(Lat);
    
    MonteCarlo(Job, Lat, Conf);
    
    Python::Finalize();
    return 0;
}

void MonteCarlo(Job& Job, Lattice& Lat, Config& Conf)
{
    InterruptHandler Interrupt;

    LOG_INFO("Markov is started!");
    timer PrinterTimer, DiskWriterTimer;
    PrinterTimer.start();
    DiskWriterTimer.start();
    
    Markov Markov(Job, Lat, Conf);
    Observables Observables(Job, Lat, Conf);

    for (uint Step = 0; Step < Job.Toss; Step++) {
        Markov.Hop(Job.Sweep);
    }
    LOG_INFO("Markov Thermalization is ended!");

    uint Step = 0;
    Job.Counter = 0;
    while (true) {
        //Don't use Para.Counter as counter
        Step++;
        Markov.Hop(Job.Sweep);
        Observables.Measure();

        if (Step % 1000 == 0) {
            Observables.AddStatistics();
            if (PrinterTimer.check(Job.PrinterTimer)) {
                Markov.PrintDetailBalanceInfo();
            }

            if (DiskWriterTimer.check(Job.DiskWriterTimer)) {
                Interrupt.Delay();
                Observables.Save();
                Interrupt.Resume();
            }
        }
    }
    LOG_INFO("Markov is ended!");
}

Lattice LoadLattice(const Job& job){
    LOG_INFO("Loading lattice is started!");
    Dictionary dict;
    dict.Load(job.LatticeName);
    Lattice lat(job.Dim, job.L, job.NSublat);
    
    if(lat.Vol!=dict.Get<int>("Vol"))
        ABORT("LoadLattice: Vol error!");
    if(lat.SublatVol!=dict.Get<int>("NSublat"))
        ABORT("LoadLattice: NSublat error!");
    if(job.LatticeName!=dict.Get<string>("Name"))
        ABORT("LoadLattice: Lattice Name error!");
    
    lat.NeighNum = dict.Get<int>("NeighNum");
    
    auto siteSublat = dict.Get<Python::ArrayObject>("SiteSublat");
    real* sublatPointer = siteSublat.Data<real>();

    auto siteVectors = dict.Get<Python::ArrayObject>("SiteVectors");
    real* vectorPointer = siteVectors.Data<real>();

    auto couplingArray = dict.Get<Python::ArrayObject>("Coupling");
    real* couplingPointer = couplingArray.Data<real>();
    
    auto arrayShape = couplingArray.Shape();

    if(arrayShape[0]!=lat.SublatVol)
	ABORT("LoadLattice: lattice array shape error!");
    if(arrayShape[1]!=lat.SublatVol)
	ABORT("LoadLattice: lattice array shape error!");
    if(arrayShape[2]!=job.L[0])
	ABORT("LoadLattice: lattice array shape error!");
    if(arrayShape[3]!=job.L[1])
	ABORT("LoadLattice: lattice array shape error!");
    
    for(int i=0; i<arrayShape[0]; i++){
        for(int j=0; j<arrayShape[1]; j++){
	    for(int x=0;  x<arrayShape[2]; x++){
                for(int y=0; y<arrayShape[3]; y++){
		    lat.Coupling.at(i,j,x,y) = (*couplingPointer);
		    couplingPointer++;
		}
	    }
        }
    }

    auto siteShape = siteVectors.Shape();

    for(int i=0; i<siteShape[0]; i++){
	lat.SiteSublat[i] = int(*sublatPointer);
	sublatPointer++;
	for(int j=0; j<siteShape[1]; j++){
	    lat.SiteVectors.at(i,j) = int(*vectorPointer);
	    vectorPointer++;
	}
    }

    //if(job.LatticeName=="Pyrochlore"){
        //lat.TetrahedraNum = dict.Get<int>("TetrahedraNum");
        //auto TetrahedraArray = dict.Get<Python::ArrayObject>("TetrahedraList");
        //real* TetrahedraPointer = TetrahedraArray.Data<real>();

        //auto arrayShape = TetrahedraArray.Shape();
        //if(arrayShape[0]!=lat.TetrahedraNum)
            //ABORT("LoadLattice: tetrahedra array shape error!");

        //for(int site=0; site<arrayShape[0]; site++){
            //for(int nb=0; nb<arrayShape[1]; nb++){
                //lat.Tetrahedra[site][nb] = int(*TetrahedraPointer);
                //TetrahedraPointer++;
            //}
        //}
    //}
//
    LOG_INFO("Loading lattice is finished!");
    return lat;
}

