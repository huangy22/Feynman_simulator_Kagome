#!/usr/bin/env python
from logger import *
import Model as lattice
import numpy as np
import os, matplotlib

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

def PlotWeightvsR(Name, array, lattice, DoesSave=True):
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
        plt.savefig("ChivsR.pdf")
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

def PlotChi_2D(Chi, lat, DoesSave=True):
    omega=0

    if lat.Name=="Pyrochlore":
        #####Pyrochlore
        KList_hhl=[]
        KList_hl0=[]
        for i in range(-lat.L[0]*4, lat.L[0]*4+1):
            for j in range(-lat.L[1]*4, lat.L[1]*4+1):
                for k in range(-lat.L[2]*4, lat.L[2]*4+1):
                    kpoint = i*lat.ReciprocalLatVec[0]+j*lat.ReciprocalLatVec[1]+ \
                            k*lat.ReciprocalLatVec[2]
                    if np.abs(kpoint[0]-kpoint[1])<1e-5:
                        KList_hhl.append((i,j,k))
                    if np.abs(kpoint[2])<1e-5:
                        KList_hl0.append((i,j,k))

        bound=[[-40,40],[-40,40]]
        ######hhl
        k_hhl, ChiK_hhl=lat.FourierTransformation(Chi, KList_hhl, "Integer", bound=bound)
        ChiK_hhl=[e.real for e in ChiK_hhl]

        x_hhl=[]
        y_hhl=[]
        for e in k_hhl:
            x_hhl.append(np.sqrt(2.0)*e[0])
            y_hhl.append(e[2])

        ######hl0
        k_hl0, ChiK_hl0=lat.FourierTransformation(Chi, KList_hl0, "Integer", bound=bound)
        ChiK_hl0=[e.real for e in ChiK_hl0]
        x_hl0=[]
        y_hl0=[]

        for e in k_hl0:
            x_hl0.append(e[0])
            y_hl0.append(e[1])

        plt.figure(1)
        ax1=plt.subplot(121,aspect='equal')
        plt.scatter(x_hhl,y_hhl,c=ChiK_hhl, s=29, edgecolor="black", linewidth=0)
        plt.xlabel("Direction [hh0]")
        plt.ylabel("Direction [00l]")
        plt.xlim(-30, 30)
        plt.ylim(-30, 30)
        label=np.linspace(min(ChiK_hhl),max(ChiK_hhl), 4)

        PI2=2*np.pi
        sqrt2 = np.sqrt(2.0)
        xlist = PI2*sqrt2*np.array([-0.75,-0.25, 0.25, 0.75, 0.25,-0.25,-0.75])
        ylist = PI2*np.array([          0,    1,    1,    0,   -1,   -1,    0])
        lc="white"
        plt.plot(xlist, ylist, color=lc)
        plt.plot(xlist, ylist+2*PI2, color=lc)
        plt.plot(xlist, ylist-2*PI2, color=lc)
        plt.plot(xlist+sqrt2*PI2, ylist+1*PI2, color=lc)
        plt.plot(xlist-sqrt2*PI2, ylist+1*PI2, color=lc)
        plt.plot(xlist+sqrt2*PI2, ylist-1*PI2, color=lc)
        plt.plot(xlist-sqrt2*PI2, ylist-1*PI2, color=lc)
        c = plt.colorbar(orientation='horizontal', shrink=0.8, ticks=label)
        c.set_label("magnitude")

        ax2=plt.subplot(122,aspect='equal')
        plt.scatter(x_hl0,y_hl0,c=ChiK_hl0, s=18, edgecolor="black", linewidth=0)
        plt.xlabel("Direction [h00]")
        plt.ylabel("Direction [0l0]")
        label=np.linspace(min(ChiK_hl0),max(ChiK_hl0), 4)
        plt.xlim(-40, 40)
        plt.ylim(-40, 40)
        plt.tight_layout()

        xlist = PI2*np.array([-1.0,-0.5, 0.5, 1.0, 1.0, 0.5,-0.5,-1.0,-1.0])
        ylist = PI2*np.array([ 0.5, 1.0, 1.0, 0.5,-0.5,-1.0,-1.0,-0.5, 0.5])
        plt.plot(xlist, ylist, color=lc)
        plt.plot(xlist+2*PI2, ylist, color=lc)
        plt.plot(xlist, ylist+2*PI2, color=lc)
        plt.plot(xlist+2*PI2, ylist+2*PI2, color=lc)
        c = plt.colorbar(orientation='horizontal', shrink=0.8, ticks=label)
        c.set_label("magnitude")

    elif lat.Name in ["3DCheckerboard", "Cubic"]:
        ####3D Checkerboard
        KList_hl0=[]
        KList_hhl=[]

        for i in range(-2*lat.L[0]+1, 2*lat.L[0]):
            for j in range(-2*lat.L[1]+1, 2*lat.L[1]):
                KList_hl0.append((i,j,0))

        k_hl0, ChiK_hl0=lat.FourierTransformation(Chi,
                KList_hl0, "Integer")
        ChiK_hl0=[e.real for e in ChiK_hl0]
        x_hl0=[]
        y_hl0=[]
        for e in k_hl0:
            x_hl0.append(e[0])
            y_hl0.append(e[1])

        plt.figure(1)
        plt.scatter(x_hl0,y_hl0,c=ChiK_hl0, s=10, edgecolor="black", linewidth=0)
        c = plt.colorbar(orientation='horizontal')
        c.set_label("magnitude")
        plt.axis('equal')

    elif lat.Dim==2:
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
        #Ktemp=[(-3,0),(3,0),(0,3),(0,-3)]
        #k, ChiK=lat.FourierTransformation(Chi,
                #Ktemp, "Integer")
    else:
        print "Lattice PlotChi_2D not implemented yet!"

    if DoesSave:
        plt.savefig("chiK_{0}.pdf".format(lat.Name))
    else:
        plt.show()
    plt.close()
    log.info("Plotting done!")


if __name__=="__main__":
    import Model as lattice
    import IO

    #Beta = 0.008
    Beta = 0.1
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
    PlotChi_2D(Chi, l, False)
    PlotWeightvsR("\chi", Chi,l, False)
    PlotChi_2D(Chi, l)
    PlotWeightvsR("\chi", Chi,l)

