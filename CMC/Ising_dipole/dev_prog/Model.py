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
                        N[i].append((np.array(coord), j, self.J1/(2.0*dist)**3.0));
        self.Neighbor=np.array([N[A], N[B], N[C]])

        self.Coupling = np.zeros((3,3,Lx, Ly))
        for i in range(3):
            for j in range(3):
                e = 0
                for y in range(Ly):
                    for x in range(Lx):
                        dist,vec = self.GetRealDistance(np.array([x, y]), i, j)
                        if not(i==j and e==0):
                            self.Coupling[i,j,x,y] = self.J1/(2.0*dist)**3.0
                        e = e + 1


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
            neighSite = self.GetSiteNumber(neiCoord, neiSublat)
            neighlist.append(neighSite)
            neighcouplist.append(neiCoup)
        return neighnum, neighlist, neighcouplist

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

    def GetNeighborArrays_New(self):
        Vol=1
        for i in range(self.Dim):
            Vol*=self.L[i]
        Vol*=self.NSublat

        SiteVectors=np.zeros((Vol,2))
        for i in range(Vol):
            SiteVectors[i]=self.GetSiteCoordi(i)[0]

        SiteSublat=np.zeros(Vol)
        for i in range(Vol):
            SiteSublat[i]=self.GetSiteCoordi(i)[1]

        Coupling = np.zeros((self.NSublat, self.NSublat, self.L[0], self.L[1]))
        Coupling = self.Coupling
        return Vol, SiteVectors, SiteSublat, Coupling 

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
        self.L=[4,4]
        #self.Name="Triangular"
        #self.Name="Square"
        #self.Name="Honeycomb"
        self.Name="Kagome"

        ####### 3D Lattices##############
        #self.Name="Cubic"
        #self.L=np.array([8,8,8])
        #self.Name="Pyrochlore"

    def test_site_map(self):

        lattice=Lattice(self.Name, self.L)

        if lattice.Name=="Kagome":
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
            vol, vectors, sublat, coupling_dr = lattice.GetNeighborArrays_New()
            vol, maxnum, neighbor, coupling = lattice.GetNeighborArray()

            for num1 in range(vol):
                for num2 in range(vol):
                    if num1!=num2:
                        vec1 = vectors[num1]
                        vec2 = vectors[num2]
                        sub1 = sublat[num1]
                        sub2 = sublat[num2]
                        dx = vec2[0]-vec1[0]
                        if dx<0:
                            dx = dx + self.L[0]
                        dy = vec2[1]-vec1[1]
                        if dy<0:
                            dy = dy + self.L[1]

                        coup = coupling_dr[sub1, sub2, dx, dy]
                        i = np.where(neighbor[num1] == num2)
                        if coup != coupling[num1, i]:
                            print vec1, sub1, vec2, sub2
                            print "correct: ", neighbor[num1, i], coupling[num1, i]
                            print "wrong: ", coup
                        self.assertEqual(coup,coupling[num1, i])
                    else:
                        vec1 = vectors[num1]
                        vec2 = vectors[num2]
                        sub1 = sublat[num1]
                        sub2 = sublat[num2]
                        coup = coupling_dr[sub1, sub2, 0, 0]
                        self.assertEqual(coup,0.0)



if __name__=='__main__':
    #unittest.main()

    L=np.array([16,16])
    Name="Kagome"

    lattice=Lattice(Name, L)
    nsublat = lattice.NSublat

    #vol, maxnum, neighbor, coupling = lattice.GetNeighborArray()

    #dict={"Vol": vol, "NSublat":nsublat, "Name": Name, "NeighNum": int(maxnum), "Neighbor": neighbor, "Coupling": coupling}
    #IO.SaveDict(Name, "w", dict)

    vol, vectors, sublat, coupling_dr = lattice.GetNeighborArrays_New()
    dict={"Vol": vol, "NSublat":nsublat, "Name": Name, "NeighNum": vol, "SiteVectors": vectors, "SiteSublat": sublat, "Coupling": coupling_dr}
    IO.SaveDict(Name, "w", dict)


    points, lines=lattice.GetSitesList(SubLatIn=0)
    interaction = lattice.GetInteraction()

    dict={"Points": points, "Lines": lines, "Interaction": interaction}
    IO.SaveDict("Coordinate", "w", dict)

    #2D lattice
    #def __Square(self):
        #self.Dim=2
        #self.__AssertDim()
        #self.NSublat=1 #####1 sublattices
        #self.J1, self.J2= 1, 0

        #Lx, Ly =(self.L[0], self.L[1])

        #########Nearest Neighbor########################
        #NA=[]
        #NA.append((np.array([1, 0]), 0, self.J1))
        #NA.append((np.array([0, 1]), 0, self.J1))
        #NA.append((np.array([Lx-1, 0]), 0, self.J1))
        #NA.append((np.array([0, Ly-1]), 0, self.J1))

        #########Next Nearest Neighbor########################
        #NA.append((np.array([1, 1]), 0, self.J2))
        #NA.append((np.array([1, Ly-1]), 0, self.J2))
        #NA.append((np.array([Lx-1, 1]), 0, self.J2))
        #NA.append((np.array([Lx-1, Ly-1]), 0, self.J2))

        #self.Neighbor=np.array([NA])

        #self.LatVec=np.array([[1,0],
                              #[0,1]])
        #self.SubLatVec=np.array([[0,0]])

    #def __Triangular(self):
        #self.Dim=2
        #self.__AssertDim()
        #self.NSublat=1 #####1 sublattices
        #self.J1, self.J2= 1, 0

        #Lx, Ly =(self.L[0], self.L[1])

        #########Nearest Neighbor########################
        #NA=[]
        #NA.append((np.array([1, 0]), 0, self.J1))
        #NA.append((np.array([0, 1]), 0, self.J1))
        #NA.append((np.array([Lx-1, 1]), 0, self.J1))
        #NA.append((np.array([Lx-1, 0]), 0, self.J1))
        #NA.append((np.array([0, Ly-1]), 0, self.J1))
        #NA.append((np.array([1, Ly-1]), 0, self.J1))

        #########Next Nearest Neighbor########################
        ##NA.append((np.array([1, 1]), 0, self.J2))
        ##NA.append((np.array([1, Ly-1]), 0, self.J2))
        ##NA.append((np.array([Lx-1, 1]), 0, self.J2))
        ##NA.append((np.array([Lx-1, Ly-1]), 0, self.J2))

        #self.Neighbor=np.array([NA])

        #root3=math.sqrt(3)
        #self.LatVec=np.array([[1,0],
                              #[0.5,root3/2.0]])
        #self.SubLatVec=np.array([[0,0]])
        #self.ReciprocalLatVec=np.array([[2*PI, -2*PI/root3],[0, 4*PI/root3]])


    #def __Honeycomb(self):
        #self.Dim=2
        #self.__AssertDim()
        #self.NSublat=2 #####2 sublattices
        #self.J1, self.J2=1, 0

        #Lx, Ly = self.L[0], self.L[1]
        #A, B = 0, 1

        #NA=[]
        #########Nearest Neighbor########################
        #NA.append((np.array([0, 0]), B, self.J1))
        #NA.append((np.array([Lx-1, 0]), B, self.J1))
        #NA.append((np.array([Lx-1, Ly-1]), B, self.J1))

        #########Next Nearest Neighbor########################
        #NA.append((np.array([1, 0]), A, self.J2))
        #NA.append((np.array([1, 1]), A, self.J2))
        #NA.append((np.array([0, 1]), A, self.J2))
        #NA.append((np.array([Lx-1, 0]), A, self.J2))
        #NA.append((np.array([Lx-1, Ly-1]), A, self.J2))
        #NA.append((np.array([0, Ly-1]), A, self.J2))

        #NB=[]
        #########Nearest Neighbor########################
        #NB.append((np.array([0, 0]), A, self.J1))
        #NB.append((np.array([1, 0]), A, self.J1))
        #NB.append((np.array([1, 1]), A, self.J1))

        #########Next Nearest Neighbor########################
        #NB.append((np.array([1, 0]), B, self.J2))
        #NB.append((np.array([1, 1]), B, self.J2))
        #NB.append((np.array([0, 1]), B, self.J2))
        #NB.append((np.array([Lx-1, 0]), B, self.J2))
        #NB.append((np.array([Lx-1, Ly-1]), B, self.J2))
        #NB.append((np.array([0, Ly-1]), B, self.J2))

        #self.Neighbor=np.array([NA, NB])

        #root3=math.sqrt(3)
        #self.LatVec=np.array([[0,1],
                             #[root3/2,-0.5]])
        #self.SubLatVec = np.array([[0, 0],
                                  #[1/2.0/root3, 0.5]])

    ##3D lattice
    #def __Cubic(self):
        #self.Dim=3
        #self.__AssertDim()
        #self.NSublat=1 #####1 sublattices
        #self.J1, self.J2=(1, 0)

        #Lx, Ly, Lz = self.L[0], self.L[1], self.L[2]
        #NA=[]
        #########Nearest Neighbor########################
        #NA.append((np.array([1, 0, 0]), 0, self.J1))
        #NA.append((np.array([0, 1, 0]), 0, self.J1))
        #NA.append((np.array([0, 0, 1]), 0, self.J1))
        #NA.append((np.array([Lx-1, 0, 0]), 0, self.J1))
        #NA.append((np.array([0, Ly-1, 0]), 0, self.J1))
        #NA.append((np.array([0, 0, Lz-1]), 0, self.J1))

        #########Next Nearest Neighbor########################
        ##NA.append((np.array([ 1,    1, 0]), 0, self.J2))
        ##NA.append((np.array([Lx-1,  1, 0]), 0, self.J2))
        ##NA.append((np.array([ 1, Ly-1, 0]), 0, self.J2))
        ##NA.append((np.array([Lx-1, Ly-1, 0]), 0, self.J2))
        ##NA.append((np.array([ 1,  0,  1]), 0, self.J2))
        ##NA.append((np.array([Lx-1,  0,  1]), 0, self.J2))
        ##NA.append((np.array([ 1,  0, Lz-1]), 0, self.J2))
        ##NA.append((np.array([Lx-1,  0, Lz-1]), 0, self.J2))
        ##NA.append((np.array([ 0,  1,  1]), 0, self.J2))
        ##NA.append((np.array([ 0, Ly-1,  1]), 0, self.J2))
        ##NA.append((np.array([ 0,  1, Lz-1]), 0, self.J2))
        ##NA.append((np.array([ 0, Ly-1, Lz-1]), 0, self.J2))

        #self.Neighbor=np.array([NA])

        #self.LatVec=np.array([[1,0,0],
                              #[0,1,0],
                              #[0,0,1]])
        #self.SubLatVec=np.array([[0,0,0]])

    #def __Pyrochlore(self):
        #self.Dim=3
        #self.__AssertDim()
        #self.NSublat=4 #####4 sublattices
        #self.J1, self.J2= 1.0, 0.0

        #Lx, Ly, Lz = self.L[0], self.L[1], self.L[2]
        #A, B, C, D = 0, 1, 2, 3


        #NA=[]
        #########Nearest Neighbor########################
        #NA.append((np.array([0, 0, 0]), B, self.J1))
        #NA.append((np.array([0, 0, 0]), C, self.J1))
        #NA.append((np.array([0, 0, 0]), D, self.J1))
        #NA.append((np.array([Lx-1, 0, 0]), B, self.J1))
        #NA.append((np.array([0, Ly-1, 0]), C, self.J1))
        #NA.append((np.array([0, 0, Lz-1]), D, self.J1))

        #NA.append((np.array([0, Ly-1, 0]), B, self.J2))
        #NA.append((np.array([0, 0, Lz-1]), B, self.J2))
        #NA.append((np.array([Lx-1, 1, 0]), B, self.J2))
        #NA.append((np.array([Lx-1, 0, 1]), B, self.J2))
        #NA.append((np.array([Lx-1, 0, 0]), C, self.J2))
        #NA.append((np.array([0, 0, Lz-1]), C, self.J2))
        #NA.append((np.array([1, Ly-1, 0]), C, self.J2))
        #NA.append((np.array([0, Ly-1, 1]), C, self.J2))
        #NA.append((np.array([Lx-1, 0, 0]), D, self.J2))
        #NA.append((np.array([0, Ly-1, 0]), D, self.J2))
        #NA.append((np.array([1, 0, Lz-1]), D, self.J2))
        #NA.append((np.array([0, 1, Lz-1]), D, self.J2))

        #NB=[]
        #########Nearest Neighbor########################
        #NB.append((np.array([0, 0, 0]), A, self.J1))
        #NB.append((np.array([0, 0, 0]), C, self.J1))
        #NB.append((np.array([0, 0, 0]), D, self.J1))
        #NB.append((np.array([1, 0, 0]), A, self.J1))
        #NB.append((np.array([1, Ly-1, 0]), C, self.J1))
        #NB.append((np.array([1, 0, Lz-1]), D, self.J1))

        #NB.append((np.array([1, Ly-1, 0]), A, self.J2))
        #NB.append((np.array([1, 0, Lz-1]), A, self.J2))
        #NB.append((np.array([0, 0, 1]), A, self.J2))
        #NB.append((np.array([0, 1, 0]), A, self.J2))
        #NB.append((np.array([1, 0, 0]), C, self.J2))
        #NB.append((np.array([1, 0, Lz-1]), C, self.J2))
        #NB.append((np.array([0, Ly-1, 0]), C, self.J2))
        #NB.append((np.array([0, Ly-1, 1]), C, self.J2))
        #NB.append((np.array([1, Ly-1, 0]), D, self.J2))
        #NB.append((np.array([1, 0, 0]), D, self.J2))
        #NB.append((np.array([0, 0, Lz-1]), D, self.J2))
        #NB.append((np.array([0, 1, Lz-1]), D, self.J2))

        #NC=[]
        #########Nearest Neighbor########################
        #NC.append((np.array([0, 0, 0]), A, self.J1))
        #NC.append((np.array([0, 0, 0]), B, self.J1))
        #NC.append((np.array([0, 0, 0]), D, self.J1))
        #NC.append((np.array([0, 1, 0]), A, self.J1))
        #NC.append((np.array([Lx-1, 1, 0]), B, self.J1))
        #NC.append((np.array([0, 1, Lz-1]), D, self.J1))

        #NC.append((np.array([Lx-1, 1, 0]), A, self.J2))
        #NC.append((np.array([0, 1, Lz-1]), A, self.J2))
        #NC.append((np.array([1, 0, 0]), A, self.J2))
        #NC.append((np.array([0, 0, 1]), A, self.J2))
        #NC.append((np.array([0, 1, 0]), B, self.J2))
        #NC.append((np.array([0, 1, Lz-1]), B, self.J2))
        #NC.append((np.array([Lx-1, 0, 0]), B, self.J2))
        #NC.append((np.array([Lx-1, 0, 1]), B, self.J2))
        #NC.append((np.array([0, 1, 0]), D, self.J2))
        #NC.append((np.array([Lx-1, 1, 0]), D, self.J2))
        #NC.append((np.array([0, 0, Lz-1]), D, self.J2))
        #NC.append((np.array([1, 0, Lz-1]), D, self.J2))

        #ND=[]
        #########Nearest Neighbor########################
        #ND.append((np.array([0, 0, 0]), A, self.J1))
        #ND.append((np.array([0, 0, 0]), B, self.J1))
        #ND.append((np.array([0, 0, 0]), C, self.J1))
        #ND.append((np.array([0, 0, 1]), A, self.J1))
        #ND.append((np.array([Lx-1, 0, 1]), B, self.J1))
        #ND.append((np.array([0, Ly-1, 1]), C, self.J1))

        #ND.append((np.array([Lx-1, 0, 1]), A, self.J2))
        #ND.append((np.array([0, Ly-1, 1]), A, self.J2))
        #ND.append((np.array([1, 0, 0]), A, self.J2))
        #ND.append((np.array([0, 1, 0]), A, self.J2))
        #ND.append((np.array([0, 0, 1]), B, self.J2))
        #ND.append((np.array([0, Ly-1, 1]), B, self.J2))
        #ND.append((np.array([Lx-1, 0, 0]), B, self.J2))
        #ND.append((np.array([Lx-1, 1, 0]), B, self.J2))
        #ND.append((np.array([0, 0, 1]), C, self.J2))
        #ND.append((np.array([Lx-1, 0, 1]), C, self.J2))
        #ND.append((np.array([0, Ly-1, 0]), C, self.J2))
        #ND.append((np.array([1, Ly-1, 0]), C, self.J2))

        #self.Neighbor=np.array([NA, NB, NC, ND])

        #self.LatVec=np.array([[0.0,0.5,0.5],
                              #[0.5,0.0,0.5],
                              #[0.5,0.5,0.0]])
        #self.SubLatVec=np.array([[0,0,0],
                                 #[0.0,0.25,0.25],
                                 #[0.25,0.0,0.25],
                                 #[0.25,0.25,0.0]])

        #self.ReciprocalLatVec=np.array([[-2*PI, 2*PI, 2*PI],
                                        #[2*PI, -2*PI, 2*PI],
                                        #[2*PI, 2*PI, -2*PI]])

        #P={"G": (0,0,0), "X":(0,2*PI,0),  "W":(PI,2*PI,0), \
           #"K":(1.5*PI,1.5*PI,0),"L": (PI,PI,PI), "U": (PI/2,2*PI,PI/2)}
        #L={"G":"$\Gamma$\n$(0,0,0)$", "X":"$X$\n$(0,2\pi,0)$", "W": "$W$\n$(\pi,2\pi,0)$", \
           #"K": "$K$\n$(3\pi/2,3\pi/2,0)$", "L": "$L$\n$(\pi,\pi,\pi)$", "U":"$U$\n$(\pi/2,2\pi,0)$"}
        #self.Path=[P["G"], P["X"], P["W"], P["K"],
                #P["G"], P["L"], P["U"], P["W"], P["L"], P["K"], P["U"], P["X"]]
        #self.PathName=[L["G"], L["X"], L["W"], L["K"],
                #L["G"], L["L"], L["U"], L["W"], L["L"], L["K"], L["U"], L["X"]]
        #self.IndependtBZCenter=[(0,0,0),(2*PI,2*PI,-2*PI),(2*PI,2*PI,2*PI),(4*PI,0,0)]

