#!/usr/bin/env python
from logger import *
import Model as lattice
import numpy as np
import os, matplotlib
import argparse

if "DISPLAY" not in os.environ:
    print "no DISPLAY detected, switch to Agg backend!"
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def PlotArray(array, Beta, Name, DoesSave=True):
    x=np.linspace(0, Beta, len(array))
    plt.figure()
    plt.plot(x,array,'-')
    if DoesSave:
        plt.savefig("{0}.jpg".format(Name))
    else:
        plt.show()
    plt.close()

def PlotSpatial(array, lattice, DoesSave=True):
    x=[]
    y=[]
    z=[]
    points, _=lattice.GetSitesList()
    for vec, coord, sub in points:
        if lattice.Dim==2 or (lattice.Dim==3 and coord[2]==0):
            x.append(vec[0])
            y.append(vec[1])
            if coord[0]==coord[1]==0:
                z.append(0.0)
            else:
                z.append()
    log.info("Max:{0}, Min: {1}".format(max(z), min(z)))
    plt.figure()
    plt.scatter(x,y,c=z, s=10, edgecolor="black", linewidth=0)
    c = plt.colorbar(orientation='horizontal', shrink=0.8, ticks=np.linspace(min(z),max(z),4))
    c.set_label("magnitude")
    plt.axis('equal')
    if DoesSave:
        plt.savefig("chi_spatial.pdf")
    else:
        plt.show()
    plt.close()

def PlotWeightvsR(Name, array, lattice, beta=None, DoesSave=True):
    plt.figure()
    x=[]
    y=[]

    SubOriginal = 0
    points,_=lattice.GetSitesList(HasOffset=False)
    for vec, coord, sub in points:
        i = lattice.GetSiteNumber(coord, sub)
        #if not (all(v==0 for v in coord) and sub==0) and all(v<l/2 for v,l in 
                #zip(coord, lattice.L)):
        if all(v<l/2 for v,l in zip(coord, lattice.L)):
            x.append(np.linalg.norm(vec))
            #y.append(array[SubOriginal, i])
            y.append(abs(array[SubOriginal, i]))

    #sort x,y according to the distance in x
    x,y = (list(x) for x in zip(*sorted(zip(x, y), key=lambda pair: pair[0])))
    #fitting
    #fitParams, fitCovariances = curve_fit(Exp, x, y)
    plt.plot(x,y, "o")
    #plt.plot(x, Exp(x, fitParams[0], fitParams[1]), '-',
            #label="fit with ${0}exp(-R/{1}a)$".format(fitParams[0], 1.0/fitParams[1]))
    plt.yscale("log")
    plt.xlabel("$R/a$")
    plt.ylim(1e-3, 1)
    plt.ylabel("$|{0}_s|$".format(Name))
    #plt.legend()
    if DoesSave:
        if beta is not None:
            plt.savefig("chivsR_Beta_{0}.pdf".format(beta))
        else:
            plt.savefig("chivsR.pdf")
    else:
        plt.show()
    plt.close()

#def PlotChiAlongPath(Chi, lat, DoesSave=True):
    #map=Chi.Map
    ##fig=plt.figure(figsize=(20, 10))
    #fig=plt.figure()
    #for BZcenter in lat.IndependtBZCenter:
        #x=[]
        #KList=[]
        #offset=0
        #ticks=[0]
        #for i in range(0, len(lat.Path)-1):
            #start, end=np.array(lat.Path[i]),np.array(lat.Path[i+1])
            #for k in lat.GetKVecAlongPath(start, end, BZcenter):
                #pos=offset+np.linalg.norm(k-np.array(BZcenter)-start)
                #x.append(pos)
                #KList.append(k)
            #offset+=np.linalg.norm(end-start)
            #ticks.append(offset)
        #_, y=lat.FourierTransformation(Chi.Data[0,:,0,:,:,omega]*map.Beta/map.MaxTauBin, KList, "Real")
        #y=[e.real for e in y]
        #BZstr=["{:.3f}".format(e) for e in BZcenter]
        ##x obtained previously may from big to small, so we have to reorder x here
        #x,y=zip(*sorted(zip(x, y)))
        #plt.plot(x,y,'o-', label="BZ:({0})".format(",".join(BZstr)))
    #plt.legend(loc='best', fancybox=True, framealpha=0.5)
    #plt.xticks(ticks, lat.PathName)
    #if DoesSave:
        #plt.savefig("chi_1D.pdf")
    #else:
        #plt.show()
    #plt.close()

def PlotChi_2D(Chi, lat, beta=None, DoesSave=True):
    omega=0

    KList=[]
    for i in range(-2*lat.L[0], 2*lat.L[0]+1):
        for j in range(-2*lat.L[1], 2*lat.L[1]+1):
            KList.append((i,j))
    k, ChiK=lat.FourierTransformation(Chi,
            KList, "Integer", bound=[[-20,20], [-20,20]])
    ChiK=[e.real for e in ChiK]
    k=np.array(k)
    plt.figure()
    plt.scatter(k[:, 0],k[:, 1],c=ChiK, s=6, edgecolor="black", linewidth=0)
    c = plt.colorbar(orientation='horizontal')
    c.set_label("magnitude")
    plt.axis('equal')

    if DoesSave:
        if beta is not None:
            plt.savefig("chiK_{0}_Beta_{1}.pdf".format(lat.Name, beta))
        else:
            plt.savefig("chiK_{0}.pdf".format(lat.Name))
    else:
        plt.show()
    plt.close()
    log.info("Plotting done!")


if __name__=="__main__":
    import Model as lattice
    import IO

    parser = argparse.ArgumentParser()
    parser.add_argument("beta", help="use beta value to find the data file")
    args = parser.parse_args()
    Beta=args.beta

    print Beta

    #l=lattice.Lattice("Kagome",[8,8])
    l=lattice.Lattice("Kagome",[16,16])

    Vol = l.NSublat
    for i in l.L:
        Vol =Vol*i

    workspace="./"

    FileList = [f for f in os.listdir(workspace) if os.path.isfile(os.path.join(workspace,f))]
    StatisFileList=[os.path.join(workspace, f) for f in FileList if f.find("{0}_Quantities_".format(Beta)) is not -1]

    print StatisFileList

    filenum=0
    Dicts=[]
    for f in StatisFileList:
        Dicts.append(IO.LoadDict(f))
        filenum +=1

    Chi=np.zeros((l.NSublat, Vol))
    for i in range(filenum):
        Chi += Dicts[i]["Correlation"]
    Chi = Chi/filenum

    Quantities = Dicts[0].keys()
    Dict={}
    Dict['Chi'] = Chi
    for e in Quantities:
        if e is not "Correlation":
            Dict[e]={}
            Data = 0.0
            Err = 0.0
            for i in range(filenum):
                Data += Dicts[i][e]['Estimation']['Mean']
                Err += Dicts[i][e]['Estimation']['Error']
            Data = Data/filenum
            Err = Err/filenum/np.sqrt(filenum)
            Dict[e]['Value']=Data
            Dict[e]['Error']=Err

    IO.SaveDict("{0}_Quantities".format(Beta), "w", Dict)

    #PlotChiAlongPath(Chi, l)
    PlotChi_2D(Chi, l, Beta)
    PlotWeightvsR("\chi", Chi,l, Beta)
    PlotChi_2D(Chi, l, DoesSave=False)
    PlotWeightvsR("\chi", Chi,l, DoesSave=False)

