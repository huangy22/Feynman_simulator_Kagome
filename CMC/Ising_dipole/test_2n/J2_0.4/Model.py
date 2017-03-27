#!/usr/bin/env python
import numpy as np
import math
import random
import itertools
import unittest
from logger import *
PI=np.pi

class Lattice:
    def __init__(self, Name, L):
        self.L=L
        self.Name=Name
        if Name=="Honeycomb":
            self.__Honeycomb()
        elif Name=="Triangular":
            self.__Triangular()
        elif Name=="Kagome":
            self.__Kagome()
        elif Name=="Square":
            self.__Square()
        elif Name=="Cubic":
            self.__Cubic()
        elif Name=="Pyrochlore":
            self.__Pyrochlore()
        else:
            Assert(False, "Lattice {0} has not been implemented!".format(self.Name))

    def __Kagome(self):
        self.Dim=2
        self.__AssertDim()
        self.NSublat=3 #####3 sublattices
        self.J1=1

        Lx, Ly = self.L[0], self.L[1]
        A, B, C = 0, 1, 2


        root3=math.sqrt(3)
        self.LatVec=np.array([[1, 0],
                             [0.5,root3/2]])
        #self.SubLatVec = np.array([[0.5, 0],
                                  #[1/4.0, root3/4],
                                  #[3.0/4.0, root3/4]])
        self.SubLatVec = np.array([[-1/4.0, -root3/12.0],
                                  [1/4.0, -root3/12.0],
                                  [0.0, root3/6.0]])
        self.ReciprocalLatVec=np.array([[2*PI, -2*PI/root3],[0, 4*PI/root3]])

        N = []
        for i in range(3):
            N.append([])
            for j in range(3):
                for e in range(Lx*Ly):
                    coord = self.__GetSiteCoordi(e, j)
                    dist,vec = self.GetRealDistance(coord, i, j)
                    if not(i==j and e==0):
                        if dist<=0.5*root3:
                            if dist ==0.5:
                                N[i].append((np.array(coord), j, self.J1))
                            else:
                                N[i].append((np.array(coord), j, self.J1*0.4));
        self.Neighbor=np.array([N[A], N[B], N[C]])


    #2D lattice
    def __Square(self):
        self.Dim=2
        self.__AssertDim()
        self.NSublat=1 #####1 sublattices
        self.J1, self.J2= 1, 0

        Lx, Ly =(self.L[0], self.L[1])

        ########Nearest Neighbor########################
        NA=[]
        NA.append((np.array([1, 0]), 0, self.J1))
        NA.append((np.array([0, 1]), 0, self.J1))
        NA.append((np.array([Lx-1, 0]), 0, self.J1))
        NA.append((np.array([0, Ly-1]), 0, self.J1))

        ########Next Nearest Neighbor########################
        NA.append((np.array([1, 1]), 0, self.J2))
        NA.append((np.array([1, Ly-1]), 0, self.J2))
        NA.append((np.array([Lx-1, 1]), 0, self.J2))
        NA.append((np.array([Lx-1, Ly-1]), 0, self.J2))

        self.Neighbor=np.array([NA])

        self.LatVec=np.array([[1,0],
                              [0,1]])
        self.SubLatVec=np.array([[0,0]])

    def __Triangular(self):
        self.Dim=2
        self.__AssertDim()
        self.NSublat=1 #####1 sublattices
        self.J1, self.J2= 1, 0

        Lx, Ly =(self.L[0], self.L[1])

        ########Nearest Neighbor########################
        NA=[]
        NA.append((np.array([1, 0]), 0, self.J1))
        NA.append((np.array([0, 1]), 0, self.J1))
        NA.append((np.array([Lx-1, 1]), 0, self.J1))
        NA.append((np.array([Lx-1, 0]), 0, self.J1))
        NA.append((np.array([0, Ly-1]), 0, self.J1))
        NA.append((np.array([1, Ly-1]), 0, self.J1))

        ########Next Nearest Neighbor########################
        #NA.append((np.array([1, 1]), 0, self.J2))
        #NA.append((np.array([1, Ly-1]), 0, self.J2))
        #NA.append((np.array([Lx-1, 1]), 0, self.J2))
        #NA.append((np.array([Lx-1, Ly-1]), 0, self.J2))

        self.Neighbor=np.array([NA])

        root3=math.sqrt(3)
        self.LatVec=np.array([[1,0],
                              [0.5,root3/2.0]])
        self.SubLatVec=np.array([[0,0]])
        self.ReciprocalLatVec=np.array([[2*PI, -2*PI/root3],[0, 4*PI/root3]])


    def __Honeycomb(self):
        self.Dim=2
        self.__AssertDim()
        self.NSublat=2 #####2 sublattices
        self.J1, self.J2=1, 0

        Lx, Ly = self.L[0], self.L[1]
        A, B = 0, 1

        NA=[]
        ########Nearest Neighbor########################
        NA.append((np.array([0, 0]), B, self.J1))
        NA.append((np.array([Lx-1, 0]), B, self.J1))
        NA.append((np.array([Lx-1, Ly-1]), B, self.J1))

        ########Next Nearest Neighbor########################
        NA.append((np.array([1, 0]), A, self.J2))
        NA.append((np.array([1, 1]), A, self.J2))
        NA.append((np.array([0, 1]), A, self.J2))
        NA.append((np.array([Lx-1, 0]), A, self.J2))
        NA.append((np.array([Lx-1, Ly-1]), A, self.J2))
        NA.append((np.array([0, Ly-1]), A, self.J2))

        NB=[]
        ########Nearest Neighbor########################
        NB.append((np.array([0, 0]), A, self.J1))
        NB.append((np.array([1, 0]), A, self.J1))
        NB.append((np.array([1, 1]), A, self.J1))

        ########Next Nearest Neighbor########################
        NB.append((np.array([1, 0]), B, self.J2))
        NB.append((np.array([1, 1]), B, self.J2))
        NB.append((np.array([0, 1]), B, self.J2))
        NB.append((np.array([Lx-1, 0]), B, self.J2))
        NB.append((np.array([Lx-1, Ly-1]), B, self.J2))
        NB.append((np.array([0, Ly-1]), B, self.J2))

        self.Neighbor=np.array([NA, NB])

        root3=math.sqrt(3)
        self.LatVec=np.array([[0,1],
                             [root3/2,-0.5]])
        self.SubLatVec = np.array([[0, 0],
                                  [1/2.0/root3, 0.5]])

    #3D lattice
    def __Cubic(self):
        self.Dim=3
        self.__AssertDim()
        self.NSublat=1 #####1 sublattices
        self.J1, self.J2=(1, 0)

        Lx, Ly, Lz = self.L[0], self.L[1], self.L[2]
        NA=[]
        ########Nearest Neighbor########################
        NA.append((np.array([1, 0, 0]), 0, self.J1))
        NA.append((np.array([0, 1, 0]), 0, self.J1))
        NA.append((np.array([0, 0, 1]), 0, self.J1))
        NA.append((np.array([Lx-1, 0, 0]), 0, self.J1))
        NA.append((np.array([0, Ly-1, 0]), 0, self.J1))
        NA.append((np.array([0, 0, Lz-1]), 0, self.J1))

        ########Next Nearest Neighbor########################
        #NA.append((np.array([ 1,    1, 0]), 0, self.J2))
        #NA.append((np.array([Lx-1,  1, 0]), 0, self.J2))
        #NA.append((np.array([ 1, Ly-1, 0]), 0, self.J2))
        #NA.append((np.array([Lx-1, Ly-1, 0]), 0, self.J2))
        #NA.append((np.array([ 1,  0,  1]), 0, self.J2))
        #NA.append((np.array([Lx-1,  0,  1]), 0, self.J2))
        #NA.append((np.array([ 1,  0, Lz-1]), 0, self.J2))
        #NA.append((np.array([Lx-1,  0, Lz-1]), 0, self.J2))
        #NA.append((np.array([ 0,  1,  1]), 0, self.J2))
        #NA.append((np.array([ 0, Ly-1,  1]), 0, self.J2))
        #NA.append((np.array([ 0,  1, Lz-1]), 0, self.J2))
        #NA.append((np.array([ 0, Ly-1, Lz-1]), 0, self.J2))

        self.Neighbor=np.array([NA])

        self.LatVec=np.array([[1,0,0],
                              [0,1,0],
                              [0,0,1]])
        self.SubLatVec=np.array([[0,0,0]])

    def __Pyrochlore(self):
        self.Dim=3
        self.__AssertDim()
        self.NSublat=4 #####4 sublattices
        self.J1, self.J2= 1.0, 0.0

        Lx, Ly, Lz = self.L[0], self.L[1], self.L[2]
        A, B, C, D = 0, 1, 2, 3


        NA=[]
        ########Nearest Neighbor########################
        NA.append((np.array([0, 0, 0]), B, self.J1))
        NA.append((np.array([0, 0, 0]), C, self.J1))
        NA.append((np.array([0, 0, 0]), D, self.J1))
        NA.append((np.array([Lx-1, 0, 0]), B, self.J1))
        NA.append((np.array([0, Ly-1, 0]), C, self.J1))
        NA.append((np.array([0, 0, Lz-1]), D, self.J1))

        NA.append((np.array([0, Ly-1, 0]), B, self.J2))
        NA.append((np.array([0, 0, Lz-1]), B, self.J2))
        NA.append((np.array([Lx-1, 1, 0]), B, self.J2))
        NA.append((np.array([Lx-1, 0, 1]), B, self.J2))
        NA.append((np.array([Lx-1, 0, 0]), C, self.J2))
        NA.append((np.array([0, 0, Lz-1]), C, self.J2))
        NA.append((np.array([1, Ly-1, 0]), C, self.J2))
        NA.append((np.array([0, Ly-1, 1]), C, self.J2))
        NA.append((np.array([Lx-1, 0, 0]), D, self.J2))
        NA.append((np.array([0, Ly-1, 0]), D, self.J2))
        NA.append((np.array([1, 0, Lz-1]), D, self.J2))
        NA.append((np.array([0, 1, Lz-1]), D, self.J2))

        NB=[]
        ########Nearest Neighbor########################
        NB.append((np.array([0, 0, 0]), A, self.J1))
        NB.append((np.array([0, 0, 0]), C, self.J1))
        NB.append((np.array([0, 0, 0]), D, self.J1))
        NB.append((np.array([1, 0, 0]), A, self.J1))
        NB.append((np.array([1, Ly-1, 0]), C, self.J1))
        NB.append((np.array([1, 0, Lz-1]), D, self.J1))

        NB.append((np.array([1, Ly-1, 0]), A, self.J2))
        NB.append((np.array([1, 0, Lz-1]), A, self.J2))
        NB.append((np.array([0, 0, 1]), A, self.J2))
        NB.append((np.array([0, 1, 0]), A, self.J2))
        NB.append((np.array([1, 0, 0]), C, self.J2))
        NB.append((np.array([1, 0, Lz-1]), C, self.J2))
        NB.append((np.array([0, Ly-1, 0]), C, self.J2))
        NB.append((np.array([0, Ly-1, 1]), C, self.J2))
        NB.append((np.array([1, Ly-1, 0]), D, self.J2))
        NB.append((np.array([1, 0, 0]), D, self.J2))
        NB.append((np.array([0, 0, Lz-1]), D, self.J2))
        NB.append((np.array([0, 1, Lz-1]), D, self.J2))

        NC=[]
        ########Nearest Neighbor########################
        NC.append((np.array([0, 0, 0]), A, self.J1))
        NC.append((np.array([0, 0, 0]), B, self.J1))
        NC.append((np.array([0, 0, 0]), D, self.J1))
        NC.append((np.array([0, 1, 0]), A, self.J1))
        NC.append((np.array([Lx-1, 1, 0]), B, self.J1))
        NC.append((np.array([0, 1, Lz-1]), D, self.J1))

        NC.append((np.array([Lx-1, 1, 0]), A, self.J2))
        NC.append((np.array([0, 1, Lz-1]), A, self.J2))
        NC.append((np.array([1, 0, 0]), A, self.J2))
        NC.append((np.array([0, 0, 1]), A, self.J2))
        NC.append((np.array([0, 1, 0]), B, self.J2))
        NC.append((np.array([0, 1, Lz-1]), B, self.J2))
        NC.append((np.array([Lx-1, 0, 0]), B, self.J2))
        NC.append((np.array([Lx-1, 0, 1]), B, self.J2))
        NC.append((np.array([0, 1, 0]), D, self.J2))
        NC.append((np.array([Lx-1, 1, 0]), D, self.J2))
        NC.append((np.array([0, 0, Lz-1]), D, self.J2))
        NC.append((np.array([1, 0, Lz-1]), D, self.J2))

        ND=[]
        ########Nearest Neighbor########################
        ND.append((np.array([0, 0, 0]), A, self.J1))
        ND.append((np.array([0, 0, 0]), B, self.J1))
        ND.append((np.array([0, 0, 0]), C, self.J1))
        ND.append((np.array([0, 0, 1]), A, self.J1))
        ND.append((np.array([Lx-1, 0, 1]), B, self.J1))
        ND.append((np.array([0, Ly-1, 1]), C, self.J1))

        ND.append((np.array([Lx-1, 0, 1]), A, self.J2))
        ND.append((np.array([0, Ly-1, 1]), A, self.J2))
        ND.append((np.array([1, 0, 0]), A, self.J2))
        ND.append((np.array([0, 1, 0]), A, self.J2))
        ND.append((np.array([0, 0, 1]), B, self.J2))
        ND.append((np.array([0, Ly-1, 1]), B, self.J2))
        ND.append((np.array([Lx-1, 0, 0]), B, self.J2))
        ND.append((np.array([Lx-1, 1, 0]), B, self.J2))
        ND.append((np.array([0, 0, 1]), C, self.J2))
        ND.append((np.array([Lx-1, 0, 1]), C, self.J2))
        ND.append((np.array([0, Ly-1, 0]), C, self.J2))
        ND.append((np.array([1, Ly-1, 0]), C, self.J2))

        self.Neighbor=np.array([NA, NB, NC, ND])

        self.LatVec=np.array([[0.0,0.5,0.5],
                              [0.5,0.0,0.5],
                              [0.5,0.5,0.0]])
        self.SubLatVec=np.array([[0,0,0],
                                 [0.0,0.25,0.25],
                                 [0.25,0.0,0.25],
                                 [0.25,0.25,0.0]])

        self.ReciprocalLatVec=np.array([[-2*PI, 2*PI, 2*PI],
                                        [2*PI, -2*PI, 2*PI],
                                        [2*PI, 2*PI, -2*PI]])

        P={"G": (0,0,0), "X":(0,2*PI,0),  "W":(PI,2*PI,0), \
           "K":(1.5*PI,1.5*PI,0),"L": (PI,PI,PI), "U": (PI/2,2*PI,PI/2)}
        L={"G":"$\Gamma$\n$(0,0,0)$", "X":"$X$\n$(0,2\pi,0)$", "W": "$W$\n$(\pi,2\pi,0)$", \
           "K": "$K$\n$(3\pi/2,3\pi/2,0)$", "L": "$L$\n$(\pi,\pi,\pi)$", "U":"$U$\n$(\pi/2,2\pi,0)$"}
        self.Path=[P["G"], P["X"], P["W"], P["K"],
                P["G"], P["L"], P["U"], P["W"], P["L"], P["K"], P["U"], P["X"]]
        self.PathName=[L["G"], L["X"], L["W"], L["K"],
                L["G"], L["L"], L["U"], L["W"], L["L"], L["K"], L["U"], L["X"]]
        self.IndependtBZCenter=[(0,0,0),(2*PI,2*PI,-2*PI),(2*PI,2*PI,2*PI),(4*PI,0,0)]

    def __AssertDim(self):
        Assert(len(self.L)==self.Dim, "Dimension {0} is expected for {1} Lattice, not {2}" \
                .format(self.Dim, self.Name, len(self.L)))

    def __Shift(self,Coordi):
        v=list(Coordi) #make a copy of Vec
        for i in range(len(Coordi)):
            if v[i]<0:
                v[i]+=self.L[i]
            if v[i]>=self.L[i]:
                v[i]-=self.L[i]
        return np.array(v)

    def __Distance(self,vec):
        return np.sum(np.abs(vec)**2,axis=-1)**(1./2)

    def GetRealDistance(self, Coordi, SubLatIn, SubLatOut):
        offset = np.array(self.L)/2-1
        v=self.__Shift(np.array(Coordi)+offset) - self.__Shift(offset)
        distvec = np.einsum("ij,i->j",self.LatVec,v)+self.SubLatVec[SubLatOut]-self.SubLatVec[SubLatIn]
        return self.__Distance(distvec), distvec

    def GetSiteNumber(self, Coordi, Sublat):
        Num=0
        for i in range(0, len(Coordi)):
            tmp = Coordi[i]
            for j in range(0, i):
                tmp *= self.L[j]
            Num+=tmp
        return int(self.NSublat*Num + Sublat)

    def GetSiteCoordi(self, Number):
        Num = int(Number/self.NSublat)
        Sublat = int(Number%self.NSublat)
        Coordi = np.zeros(self.Dim)
        for i in range(self.Dim-1,-1,-1):
            tmp = 1
            for j in range(0, i):
                tmp *= self.L[j]
            Coordi[i]=int(Num/tmp)
            Num = int(Num%tmp)
        return Coordi, Sublat

    def __GetSiteCoordi(self, Num, Sublat):
        """Only used in Kagome dipolar subroutine"""
        Coordi = np.zeros(self.Dim)
        for i in range(self.Dim-1,-1,-1):
            tmp = 1
            for j in range(0, i):
                tmp *= self.L[j]
            Coordi[i]=int(Num/tmp)
            Num = int(Num%tmp)
        return Coordi

    def GetSiteNeighbor(self, Number):
        Coord, Sublat = self.GetSiteCoordi(Number)
        neighlist = []
        neighcouplist = []
        neighnum = 0
        for i in self.Neighbor[Sublat]:
            neighnum += 1
            neiCoord = self.__Shift(i[0]+Coord)
            neiSublat = i[1]
            neiCoup = i[2]
            neighlist.append(self.GetSiteNumber(neiCoord, neiSublat))
            neighcouplist.append(neiCoup)
        return neighnum, neighlist, neighcouplist

    def GetTetrahedras(self):
        if self.Name is not "Pyrochlore" :
            Assert(False, "Lattice {0} doesn't have tetrahedras!".format(self.Name))

        Vol=1
        for i in range(self.Dim):
            Vol*=self.L[i]

        TetrahedraNum = Vol*2
        TetrahedraList = np.zeros((TetrahedraNum, self.NSublat))

        for x,y,z in itertools.product(range(self.L[0]), repeat=3):
            t = x+y*self.L[0]+z*self.L[0]*self.L[1]
            Coord = np.array([x, y, z])
            site = self.GetSiteNumber(Coord, 0)
            TetrahedraList[t][0] = site
            TetrahedraList[t+Vol][0] = site
            num, neighlist, neiCoup = self.GetSiteNeighbor(site)
            for i in range(num):
                neiCoord, neiSub = self.GetSiteCoordi(neighlist[i])
                if all(int(i)==0 for i in neiCoord-Coord) :
                    TetrahedraList[t][neiSub] = neighlist[i]
                else:
                    TetrahedraList[t+Vol][neiSub] = neighlist[i]
        return TetrahedraNum, TetrahedraList

    def GetNeighborArray(self):
        Vol=1
        for i in range(self.Dim):
            Vol*=self.L[i]
        Vol*=self.NSublat

        NeighborNum=np.zeros(Vol, dtype=int)
        for i in range(Vol):
            NeighborNum[i]=self.GetSiteNeighbor(i)[0]
        MaxNum = np.max(NeighborNum)

        Neighbor=np.zeros((Vol, MaxNum))
        NeighborCoupling=np.zeros((Vol, MaxNum))
        for i in range(Vol):
            num=NeighborNum[i]
            Neighbor[i][0:num]=np.array(self.GetSiteNeighbor(i)[1])
            Neighbor[i][num:MaxNum]=-1
            NeighborCoupling[i][0:num]=np.array(self.GetSiteNeighbor(i)[2])
            NeighborCoupling[i][num:MaxNum]=0
        return Vol, MaxNum, Neighbor, NeighborCoupling 

    def GetRealVec(self, Coordi, SubLatIn, SubLatOut, offset):
        '''
           Coordi: D-dimensional vector of coordinates
           SubLat: only the OUT sublattice is needed, IN sublattice is assumed to be 0
        '''
        v=self.__Shift(Coordi+offset)
        #v=self.__Shift(Coordi+offset) - self.__Shift(offset)
        return tuple(np.einsum("ij,i->j",self.LatVec,v)+self.SubLatVec[SubLatOut]-self.SubLatVec[SubLatIn])


    def GetSitesList(self, HasOffset=False, SubLatIn=0):
        """
        return: list of all sites, with format 
                [tuple of real space vectors of sites, tuple of integer vectors of coordinates, SubLat] 
        """
        if HasOffset:
            offset=self.L/2-1
        else:
            offset=np.array([0 for e in self.L])

        Vol=1
        for i in range(self.Dim):
            Vol*=self.L[i]
        Vol*=self.NSublat

        Points=[None,]*(Vol)
        LinesInUnitCell=[]

        for i in range(Vol):
            Coord, Sub = self.GetSiteCoordi(i)
            Points[i]=[self.GetRealVec(Coord, SubLatIn, Sub, offset), Coord,Sub]
            for subN in range(Sub, self.NSublat):
                LinesInUnitCell.append([(i, self.GetSiteNumber(Coord, subN)), Sub])
        
        return Points, LinesInUnitCell

    def GetInteraction(self):
        Vol=1
        for i in range(self.Dim):
            Vol*=self.L[i]
        Vol*=self.NSublat

        #offset=self.L/2-1
        offset=np.array([0 for e in self.L])

        Interaction=[]
        for origin in range(Vol):
            Origin, sub = self.GetSiteCoordi(origin)
            nnum, neigh, _ = self.GetSiteNeighbor(origin)
            for n in neigh:
                # vec is (real vector of in-site, real vector of out-site) 
                Coord, sub2 = self.GetSiteCoordi(n)
                flag = True
                for i in range(self.Dim):
                    if (Origin[i]==0 and Coord[i]==self.L[i]-1) or (Origin[i]==self.L[i]-1
                            and Coord[i]==0):
                        flag = False

                if flag:
                    vec=(self.GetRealVec(Origin, 0, sub, offset), \
                                    self.GetRealVec(Coord, 0, sub2, offset))
                    coord=(origin, n)
                    Interaction.append([vec, coord, sub])
        return Interaction
    
    def FourierTransformation(self, Data, KCoordi, KType="Integer", bound=None):
        """Fourier Transformation in real lattice vector space:
           Data: numpy array with shape [NSublat,VOL]
           KCoordi: list of momentums needed to do fourier transformation, with type KType
           Definition:
               F(K)=\sum_{i;j} G(i;j)*exp(-i*(r_{i}-r_{j})*K)
               where i,j are IN/OUT lattice indexies, a,b are IN/OUT sublattice indexies
            """
        DataK=[]
        K=[]

        Vol=1
        for i in range(self.Dim):
            Vol*=self.L[i]
        Vol*=self.NSublat

        SiteNum=Vol
        vec=np.zeros((SiteNum*self.NSublat, self.Dim))
        data=np.zeros(SiteNum*self.NSublat)
        index=0

        for sub in range(self.NSublat):
            points,_=self.GetSitesList(HasOffset=False, SubLatIn=sub)
            for i in range(SiteNum):
                RealVec, Coord, SubLatOut=points[i]
                vec[index,:]=RealVec
                data[index] = Data[sub,i]
                index+=1

        Assert(index==SiteNum*self.NSublat, "{0} is expected for the total number of sites, not {1}!".format(SiteNum, index))
                
        for p in KCoordi:
            if KType=="Integer":
                KVec=self.ReciprocalLatVec[0,:]*p[0]/self.L[0]
                for i in range(1,self.Dim):
                    KVec+=self.ReciprocalLatVec[i,:]*p[i]/self.L[i]
            elif KType=="Real":
                KVec=np.array(p)

            flag=True
            if bound is not None:
                bound=np.array(bound)
                for d in range(0, bound.shape[0]):
                    if KVec[d]<bound[d][0] or KVec[d]>bound[d][1]:
                        flag=False
                        break
            if flag:
                f=0
                f+=np.dot(data, np.exp(-1j*np.dot(vec[:,:], KVec)))
                K.append(KVec)
                DataK.append(f/self.NSublat)
        return K, DataK

class TestLattice(unittest.TestCase):
    def setUp(self):
        ####### 2D Lattices##############
        self.L=[8,8]
        self.Name="Triangular"
        #self.Name="Square"
        #self.Name="Honeycomb"
        #self.Name="Kagome"

        ####### 3D Lattices##############
        #self.Name="Cubic"
        #self.L=np.array([8,8,8])
        #self.Name="Pyrochlore"

    def test_site_map(self):
        lattice=Model(self.Name, self.L)
        if lattice.Name=="Square":
            print "Testing Square Lattice..."
            x,y = 6,7
            sub = 0
            num = lattice.NSublat*(y*self.L[0]+x)+sub

            #Test Number and (Coordinate, Sublattice) Mapping
            self.assertEqual(lattice.GetSiteNumber((x,y), sub), num) 
            Coordi, Sublat = lattice.GetSiteCoordi(num)
            self.assertTrue(np.array_equal(Coordi, np.array([x,y])))
            self.assertEqual(Sublat, sub)

            #Test Neighbor List
            num, neilist, couplist = lattice.GetSiteNeighbor(num)
            self.assertEqual(8, num)
            self.assertEqual(len(neilist), num)
            self.assertEqual(len(couplist), num)

            for i in range(num):
                neiCoord, neiSublat = lattice.GetSiteCoordi(neilist[i])
                print couplist[i], lattice.GetRealVec(neiCoord-Coordi, sub, neiSublat, [0])

        elif lattice.Name=="Triangular":
            print "Testing Triangular Lattice..."
            x,y = 6,7
            sub = 0
            num = lattice.NSublat*(y*self.L[0]+x)+sub

            #Test Number and (Coordinate, Sublattice) Mapping
            self.assertEqual(lattice.GetSiteNumber((x,y), sub), num) 
            Coordi, Sublat = lattice.GetSiteCoordi(num)
            self.assertTrue(np.array_equal(Coordi, np.array([x,y])))
            self.assertEqual(Sublat, sub)

            #Test Neighbor List
            num, neilist, couplist = lattice.GetSiteNeighbor(num)
            self.assertEqual(9, num)
            self.assertEqual(len(neilist), num)
            self.assertEqual(len(couplist), num)

            for i in range(num):
                neiCoord, neiSublat = lattice.GetSiteCoordi(neilist[i])
                print couplist[i], lattice.GetRealVec(neiCoord-Coordi, sub, neiSublat, [0])

        elif lattice.Name=="Honeycomb":
            print "Testing Honeycomb Lattice..."
            x,y = 6,7
            sub = 0
            num = lattice.NSublat*(y*self.L[0]+x)+sub

            #Test Number and (Coordinate, Sublattice) Mapping
            self.assertEqual(lattice.GetSiteNumber((x,y), sub), num) 
            Coordi, Sublat = lattice.GetSiteCoordi(num)
            self.assertTrue(np.array_equal(Coordi, np.array([x,y])))
            self.assertEqual(Sublat, sub)

            #Test Neighbor List
            num, neilist, couplist = lattice.GetSiteNeighbor(num)
            self.assertEqual(9, num)
            self.assertEqual(len(neilist), num)
            self.assertEqual(len(couplist), num)

            for i in range(num):
                neiCoord, neiSublat = lattice.GetSiteCoordi(neilist[i])
                print couplist[i], lattice.GetRealVec(neiCoord-Coordi, sub, neiSublat, [0])

        elif lattice.Name=="Kagome":
            print "Testing Kagome Lattice..."
            x,y = 3,7
            sub = 1
            num = lattice.NSublat*(y*self.L[0]+x)+sub

            #Test Number and (Coordinate, Sublattice) Mapping
            self.assertEqual(lattice.GetSiteNumber((x,y), sub), num) 
            Coordi, Sublat = lattice.GetSiteCoordi(num)
            self.assertTrue(np.array_equal(Coordi, np.array([x,y])))
            self.assertEqual(Sublat, sub)

            #Test Neighbor List
            num, neilist, couplist = lattice.GetSiteNeighbor(num)
            self.assertEqual(4, num)
            self.assertEqual(len(neilist), num)
            self.assertEqual(len(couplist), num)

            for i in range(num):
                neiCoord, neiSublat = lattice.GetSiteCoordi(neilist[i])
                print couplist[i], lattice.GetRealVec(neiCoord-Coordi, sub, neiSublat, [0])

        elif lattice.Name=="Cubic":
            print "Testing Cubic Lattice..."
            x,y,z = 6,7,5
            sub = 0
            num = lattice.NSublat*(z*self.L[0]*self.L[1]+y*self.L[0]+x)+sub

            #Test Number and (Coordinate, Sublattice) Mapping
            self.assertEqual(lattice.GetSiteNumber((x,y,z), sub), num) 
            Coordi, Sublat = lattice.GetSiteCoordi(num)
            self.assertTrue(np.array_equal(Coordi, np.array([x,y,z])))
            self.assertEqual(Sublat, sub)

            #Test Neighbor List
            num, neilist, couplist = lattice.GetSiteNeighbor(num)
            self.assertEqual(6, num)
            self.assertEqual(len(neilist), num)
            self.assertEqual(len(couplist), num)

            for i in range(num):
                neiCoord, neiSublat = lattice.GetSiteCoordi(neilist[i])
                print couplist[i], lattice.GetRealVec(neiCoord-Coordi, sub, neiSublat, [0,0])

        elif lattice.Name=="Pyrochlore":
            print "Testing Pyrochlore Lattice..."
            x,y,z = 3,7,6
            sub = 0
            num = lattice.NSublat*(z*self.L[0]*self.L[1]+y*self.L[0]+x)+sub

            #Test Number and (Coordinate, Sublattice) Mapping
            self.assertEqual(lattice.GetSiteNumber((x,y,z), sub), num) 
            Coordi, Sublat = lattice.GetSiteCoordi(num)
            self.assertTrue(np.array_equal(Coordi, np.array([x,y,z])))
            self.assertEqual(Sublat, sub)

            #Test Neighbor List
            num, neilist, couplist = lattice.GetSiteNeighbor(num)
            self.assertEqual(6, num)
            self.assertEqual(len(neilist), num)
            self.assertEqual(len(couplist), num)

            for i in range(num):
                neiCoord, neiSublat = lattice.GetSiteCoordi(neilist[i])
                print couplist[i], lattice.GetRealVec(neiCoord-Coordi, sub, neiSublat, [0,0])

if __name__=='__main__':
    #unittest.main()

    #L=np.array([8,8])
    #Name="Square"

    L=np.array([16,16])
    Name="Kagome"

    #L=np.array([4,4,4])
    #Name="Pyrochlore"

    #L=np.array([8,8,8])
    #Name="Cubic"

    lattice=Lattice(Name, L)
    nsublat = lattice.NSublat
    vol, maxnum, neighbor, coupling = lattice.GetNeighborArray()

    if Name=="Pyrochlore":
        tetranum, tetralist = lattice.GetTetrahedras()
        dict={"Vol": vol, "NSublat":nsublat, "Name": Name, "NeighNum": int(maxnum), "Neighbor": neighbor, "Coupling": coupling, "TetrahedraNum":tetranum, "TetrahedraList":tetralist}
        IO.SaveDict(Name, "w", dict)
    else:
        dict={"Vol": vol, "NSublat":nsublat, "Name": Name, "NeighNum": int(maxnum), "Neighbor": neighbor, "Coupling": coupling}
        IO.SaveDict(Name, "w", dict)

    points, lines=lattice.GetSitesList(SubLatIn=0)
    interaction = lattice.GetInteraction()

    dict={"Points": points, "Lines": lines, "Interaction": interaction}
    IO.SaveDict("Coordinate", "w", dict)

