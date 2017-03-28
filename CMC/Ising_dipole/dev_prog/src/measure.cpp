//
//  measure.cpp
//  Feynman_Simulator
//
//  Created by yuan on 4/15/15.
//  Copyright (c) 2015 Kun Chen. All rights reserved.
//

#include "measure.h"
#include "lattice.h"
#include "configuration.h"
#include "job/job.h"
#include "utility/pyglue/pyarraywrapper.h"

using namespace measure;
using namespace lattice;
using namespace configuration;
using namespace para;
using namespace Python;

Observables::Observables(Job &job, Lattice &lattice, Config &config)
{
    Lat = &lattice;
    Conf = &config;
    Beta = job.Beta;
    LatName = job.LatticeName;
    PID = job.PID;

    Quantities.AddEstimator("Magnet^2");
    Quantities.AddEstimator("Magnet^4");
    Quantities.AddEstimator("StagMagnet^2");
    Quantities.AddEstimator("StagMagnet^4");
    Quantities.AddEstimator("BinderRatio");
    Quantities.AddEstimator("Energy");
    Quantities.AddEstimator("Energy^2");
    Quantities.AddEstimator("SpecificHeat");
    Quantities.AddEstimator("S0S1");
    Quantities.AddEstimator("S0SL/2");
    //if(D==1 && LatName=="Pyrochlore") Quantities.AddEstimator("Defect");

    for (int i = 0; i < Lat->SublatVol * Lat->Vol; i++)
        Correlation.AddEstimator(ToString(i));
}

void Observables::Measure()
{
    Vec<real> TotalSpin = Vec<real>(0.0);
    Vec<real> TotalStagSpin = Vec<real>(0.0);
    real Energy = 0.0;
    for (int site = 0; site < Lat->Vol; site++) {
        for (int sub = 0; sub < Lat->SublatVol; sub++) {
            Correlation[sub * Lat->Vol + site].Measure(Conf->Spin[sub] * Conf->Spin[site]);
        }
        TotalSpin += Conf->Spin[site] / Lat->Vol;
        TotalStagSpin += Conf->Spin[site] *_GetSign(site)/ Lat->Vol;
        for (int n = 0; n < Lat->NeighNum; n++) {
	    if(n!=site)
		Energy += Conf->Spin[site] * Conf->Spin[n];
        }
    }
    Quantities["S0S1"].Measure(Conf->Spin[0]*Conf->Spin[1]);
    int L = Lat->L[0];
    int site = (L/2+L/2*L+L/2*L*L)*Lat->SublatVol;
    Quantities["S0SL/2"].Measure(Conf->Spin[0]*Conf->Spin[site]);

    Quantities["Energy"].Measure(Energy / real(Lat->NeighNum) / real(Lat->Vol));
    Quantities["Energy^2"].Measure(pow(Energy / real(Lat->NeighNum) / real(Lat->Vol), 2.0));
    real Cv = Beta*Beta*(Quantities["Energy^2"].Estimate().Mean-pow(Quantities["Energy"].Estimate().Mean, 2.0));
    Quantities["SpecificHeat"].Measure(Cv);

    Quantities["Magnet^2"].Measure(TotalSpin * TotalSpin);
    Quantities["Magnet^4"].Measure(pow(TotalSpin * TotalSpin, 2.0));
    Quantities["StagMagnet^2"].Measure(TotalStagSpin * TotalStagSpin);
    Quantities["StagMagnet^4"].Measure(pow(TotalStagSpin * TotalStagSpin, 2.0));
    real bin = 1.0 - Quantities["Magnet^4"].Estimate().Mean / 3.0 / pow(Quantities["Magnet^2"].Estimate().Mean, 2.0);
    Quantities["BinderRatio"].Measure(bin);

    //if(D==1 && LatName=="Pyrochlore"){
        //real defect = 0.0;
        //for(int site=0; site<Lat->TetrahedraNum; site++) {
            //Vec<real> spin = Vec<real>(0.0);
            //for(int i =0; i < Lat->SublatVol; i++)
                //spin += Conf->Spin[Lat->Tetrahedra[site][i]];
            //if(fabs(spin[0])>1e-1)   defect+=1.0;
        //}
        //Quantities["Defect"].Measure(defect/real(Lat->TetrahedraNum));
    //}
}

real Observables::_GetSign(int site)
{
    real sign = 1.0;
    int x, y, z;
    if(Lat->Dimension==3 && Lat->SublatVol==1){
	z = site/Lat->L[0]/Lat->L[1];
	y = (site-z*Lat->L[0]*Lat->L[1])/Lat->L[0];
	x = site-z*Lat->L[0]*Lat->L[1]-y*Lat->L[0];
	if((x+y+z)%2==0) {
	    sign = 1.0;
	}else{
	    sign = -1.0;
	}
    }
    return sign;
}
void Observables::AddStatistics()
{
    Quantities.AddStatistics();
    Correlation.AddStatistics();
}

void Observables::Save()
{
    LOG_INFO("Start saving data...");

    Dictionary Quan;
    Quan = Quantities.ToDict();

    real Corr[Lat->SublatVol * Lat->Vol];
    for (int sub = 0; sub < Lat->SublatVol; sub++) {
        for (int i = 0; i < Lat->Vol; i++) {
            Corr[sub * Lat->Vol + i] = Correlation[sub * Lat->Vol + i].Estimate().Mean;
        }
    }
    uint Shape[2];
    Shape[0] = static_cast<unsigned int>(Lat->SublatVol);
    Shape[1] = static_cast<unsigned int>(Lat->Vol);
    ArrayObject corr = ArrayObject(Corr, Shape, 2);
    Quan["Correlation"] = corr;
    Quan.Save(ToString(Beta) + "_Quantities_" + ToString(PID), "w");
    
    LOG_INFO("Saving data is finished.");

    if(D==1){
        LOG_INFO("Snap Shot of spins saving...");
        Dictionary DictSpin;
        real Spin[Lat->Vol];
        for (int i = 0; i < Lat->Vol; i++) {
            Spin[i] = Conf->Spin[i][0];
        }

        uint SpinShape[1];
        SpinShape[0] = static_cast<unsigned int>(Lat->Vol);
        ArrayObject spin = ArrayObject(Spin, SpinShape, 1);
        DictSpin["Spin"] = spin;
        DictSpin.Save(ToString(Beta) + "_Snap_Shot_" + ToString(PID), "w");
        LOG_INFO("Saving Snap Shot is finished.");
    }
}
