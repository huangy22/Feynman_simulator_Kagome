//
//  markov.cpp
//  Feynman_Simulator
//
//  Created by yuan on 4/15/15.
//  Copyright (c) 2015 Kun Chen. All rights reserved.
//

#include <stdio.h>
#include "markov.h"
#include "lattice.h"
#include "configuration.h"
#include "job/job.h"
#include "utility/utility.h"
#include "utility/vector.h"
#include "utility/rng.h"
#include "utility/logger.h"
#include "math.h"

using namespace std;
using namespace lattice;
using namespace configuration;
using namespace markov;

Markov::Markov(Job &job, Lattice &lattice, Config &conf)
{
    Beta = job.Beta;
    ExternalField = job.ExternalField;
    Counter = &job.Counter;
    RNG = &job.RNG;
    Lat = &lattice;
    Conf = &conf;

    _Accepted = 0.0;
    _Proposed = 0.0;
}

void Markov::Hop(int sweep)
{
    for (int i = 0; i < sweep; i++) {
        int site = _RandomPickSite();
	real h = Beta*ExternalField;
        Vec<real> newspin = _RandomPickSpin();
        Vec<real> neighborspin = _SumSpin(site);
        real AccRatio = exp(-Beta * ((newspin - Conf->Spin[site]) * neighborspin));
	h = h* _GetSign(site);
	AccRatio = AccRatio* exp(h*(newspin[0]-Conf->Spin[site][0]));
        _Proposed += 1.0;
        if (AccRatio >= 1.0 || RNG->urn() < AccRatio) {
            _Accepted += 1.0;
            Conf->Spin[site] = newspin;
            (*Counter)++;
        }
    }
    return;
}

void Markov::PrintDetailBalanceInfo()
{
    string Output = "";
    Output = string(60, '=') + "\n";
    Output += "DiagCounter: " + ToString(*Counter) + "\n";
    Output += "AccepatnceRatio: " + ToString(_Accepted / _Proposed) + "\n";
    LOG_INFO(Output);
}

Vec<real> Markov::_SumSpin(int site)
{
    Vec<real> sum = Vec<real>(0.0);
    int neigh;
    for (int i = 0; i < Lat->NeighNum; i++) {
        neigh = Lat->Neighbor[site][i];
        sum += Conf->Spin[neigh] * (Lat->Coupling[site][i]);
    }
    return sum;
}

real Markov::_GetSign(int site)
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

int Markov::_RandomPickSite()
{
    return RNG->irn(0, Lat->Vol);
}

Vec<real> Markov::_RandomPickSpin()
{
    Vec<real> spin = Vec<real>(0.0);
    if (D == 3) {
        real z = RNG->urn()*2.0 -1.0;
        real phi = RNG->urn() * 2.0 * PI;
        spin[0] = sqrt(1.0-z*z) * cos(phi);
        spin[1] = sqrt(1.0-z*z) * sin(phi);
        spin[2] = z;
    }
    else if (D == 2) {
        real phi = RNG->urn() * 2.0 * PI;
        spin[0] = cos(phi);
        spin[1] = sin(phi);
    }
    else if (D == 1) {
        spin[0] = 2 * RNG->irn(0, 1) - 1;
    }
    else {
        LOG_ERROR("D is not 1, 2, or 3! Not implemented yet!");
    }
    return spin;
}
