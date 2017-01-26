import numpy as np
import random
import copy
from fitCurves import fitCurve
import playground3 as pg3
import sys
import traceback
import math
from intervaltree import IntervalTree as it
from intervaltree import Interval as iv



check = pg3.Check()
count = 0
pairs = {}
tracededges=[]
class fp(object):
    maxcount = 1
    def __init__(self,letter):
        #self.edges,self.its,self.cnts,self.vertices,self.convex_vert,self.concave_vert=
        self.modEdges(letter)
        self.untracededges=[]
        self.tracededges = []
        for k,_ in enumerate(self.edges):
            self.untracededges.append(k)



    def ccw(self,A,B,C):
            return (C[1]-A[1])*(B[0]-A[0]) > (B[1]-A[1])*(C[0]-A[0])

    def p(self,quad,t):
        if len(quad) ==3:
            return quad[0]*(1-t)**2 + quad[1]*2*t*(1-t) + quad[2]*t**2
        else:
            return quad[0]
    def pprime(self,quad,t):
        return 2*((quad[1]-quad[0])*(1-t) + (quad[2]-quad[1])*t)
    def f(self,quad1,t,quad2,r):
        p1 = self.p(quad1,t)
        p1prime = self.pprime(quad1,t)
        p2 = self.p(quad2,r)
        p2prime = self.pprime(quad2,r)
        return np.dot(p1prime,p1prime)*(np.dot(p1-p2,p2prime))**2 - np.dot(p2prime,p2prime)*(np.dot(p1-p2,p1prime))**2
    def fprime(self,quad1,t,quad2,r):
        p1 = self.p(quad1,t)
        p1prime = self.pprime(quad1,t)
        p2 = self.p(quad2,r)
        p2prime = self.pprime(quad2,r)
        p2dprime = 2*(quad2[2]+quad2[0]-2*quad2[1])
        term1 = np.dot(p1prime,p1prime)
        term2 = np.dot(p2prime,p2prime)
        term3 = np.dot(p1-p2,p1prime)
        term4 = np.dot(p1-p2,p2prime)
        term5 = np.dot(p2prime,p1prime)
        term6 = np.dot(p2prime,p2dprime)
        term7 = np.dot(p1-p2,p2dprime)

        term8 = 2*term1*term4*(term7-term2)
        term9 = 2*term6 * term3**2
        term10 = 2*term2*term3*(-1)*term5
        return term8 - term9 - term10
    def solve(self,quad1,t,quad2,r,err=10.0):
        p1 = self.p(quad1,t)
        p1prime = self.pprime(quad1,t)
        p2 = self.p(quad2,r)
        p2prime = self.pprime(quad2,r)
        normal1 = np.array([p1prime[1],-1*p1prime[0]])
        normal2 = np.array([p2prime[1],-1*p2prime[0]])
        denominator = np.dot(normal1,p2prime)
        beta = np.dot(p2-p1,p1prime)
        alpha = np.dot(p2-p1,p2prime)
        if abs(denominator) < 0.01: # parallel normals
            if abs(alpha) <err and abs(beta) < err : # overlapping normals
                #print "overlapping normal"
                return r,True
            else:
                #print "parallel normals- get out of this",val1,val2
                return r,False
                #return random.uniform(0,1),False #a random number between 0 and 1
            #print r
        else: # not parallel normal
            fp = self.fprime(quad1,t,quad2,r)
            if fp != 0:
                val1 = self.f(quad1,t,quad2,r)/fp
            else:
                #print "fp is zero"
                val1 = self.f(quad1,t,quad2,r)
            if  abs(val1) < 0.001:
                #print p1+alpha/denominator*normal1-p2-beta/denominator*normal2
                return r,True
            else:
            #    print "next iteration"
                #print r
                #print val1
                r = r - val1
                #print r
                if r>1 or r<0:
                    return random.uniform(0,1),False
                else:
                    return r,False
    def solve1(self,quad1,t,quad2,r,err=10.0):
        p1 = self.p(quad1,t)
        p1prime = self.pprime(quad1,t)
        p2 = self.p(quad2,r)
        p2prime = self.pprime(quad2,r)
        normal1 = np.array([p1prime[1],-1*p1prime[0]])
        normal2 = np.array([p2prime[1],-1*p2prime[0]])
        denominator = np.dot(normal1,p2prime)
        beta = np.dot(p2-p1,p1prime)
        alpha = np.dot(p2-p1,p2prime)
        if abs(denominator) < 0.01: # parallel normals
            if abs(alpha) <err and abs(beta) < err : # overlapping normals
                #print "overlapping normal"
                return r,True
            else:
                #print "parallel normals- get out of this",val1,val2
                return random.uniform(0,1),False #a random number between 0 and 1
            #print r
        else: # not parallel normal
            fp = self.fprime(quad1,t,quad2,r)
            if fp != 0:
                val1 = self.f(quad1,t,quad2,r)/fp
            else:
                val1 = abs(self.f(quad1,t,quad2,r))
            if  abs(val1) < 0.001:
                #print p1+alpha/denominator*normal1-p2-beta/denominator*normal2
                return r,True
            else:
            #    print "next iteration"
                #print r
                #print val1
                r = r - val1
                #print r
                if r>1 or r<0:
                    return random.uniform(0,1),False
                else:
                    return r,False
    def footprint(self,quad1,t,quad2,initialr=-1):
        p1 = self.p(quad1,t)
        p1prime = self.pprime(quad1,t)
        n1 = np.array([p1prime[1],-1*p1prime[0]])

        dist0,p2i,r1 =self.distance(quad2,p1) #randomly ..no partcular thought
        edges = []
        edges.append([p1,p2i])
        #print dist0,r1
        if dist0 ==float('inf'):
            #print "here"
            p2i = self.p(quad2,0.5);dist0 = np.linalg.norm(p2i-p1)

        for i in xrange(20):
            costheta1 = np.dot(p2i-p1,n1)/(np.linalg.norm(n1)*dist0)
            if abs(costheta1)>0.01:
                dist1 =  dist0/(2*costheta1)
            else:
                dist1 = dist0/2.0
            fpi = p1 +dist1*n1/np.linalg.norm(n1)
            dist2,p2i,r1 =self.distance(quad2,fpi) #randomly ..no partcular thought
            #print dist2,p2i,r1
            if dist2 ==float('inf'):
                break
            if abs(dist2-dist1)<1.0:
                break
            else:
                dist0 = np.linalg.norm(p2i-p1)
        edges.append([p1,fpi]);edges.append([p2i,fpi]);edges.append([p1]);edges+=[[p2i],[p1,p2i],[fpi]]
        print np.cross(fpi-p1,p1prime),abs(dist2-dist1)<1.0
        if np.cross(fpi-p1,p1prime)>=0 and abs(dist2-dist1)<1.0:
          return edges,r1,fpi
        else:
          return edges,-1,fpi



    def footprintDiscarded(self,quad1,t,quad2,initialr=-1,diagnostics=False,err=15.0,err2=10.0):
        edges = []

        def f1(quad1,t,quad2,r):
                p1 = self.p(quad1,t)
                p1prime = self.pprime(quad1,t)
                normal1 = np.array([p1prime[1],-1*p1prime[0]])
                p2 = self.p(quad2,r)
                p2prime = self.pprime(quad2,r)
                normal2 = np.array([p2prime[1],-1*p2prime[0]])
                denominator = np.dot(normal1,p2prime)
                if abs(denominator) >1.0 :
                    alpha = np.dot(p2-p1,p2prime)/denominator
                    beta = np.dot(p2-p1,p1prime)/denominator
                else:
                    alpha = beta = np.linalg.norm(p1-p2)/2.0
                    alpha = alpha/np.linalg.norm(p1prime)
                    beta = beta/np.linalg.norm(p2prime)
                #print alpha,beta

                fp1 = p1 + alpha * normal1
                fp2 = p2 + beta * normal2
                #dv = fp1-fp2
                #a1 = alpha*normal1;b1=beta*normal2

                diff = np.linalg.norm(fp1-fp2)
                #diff1 = abs(np.linalg.norm(alpha*normal1)-np.linalg.norm(beta*normal2))
                #return [abs(self.f(quad1,t,quad2,r))+5*diff,fp1,fp2,p1,p2,normal1,normal2,alpha,beta]
                return [abs(np.linalg.norm(alpha*normal1)-np.linalg.norm(beta*normal2)),fp1,fp2,p1,p2,normal1,normal2,alpha,beta]

                #return [abs(self.f(quad1,t,quad2,r))+(abs(alpha)-alpha+abs(beta)-beta)*+diff,fp1,fp2,p1,p2,normal1,normal2,alpha,beta]
        def f2(quad1,t,quad2,r):
             return np.dot(self.p(quad1,t)-self.p(quad2,r),self.pprime(quad2,r))
        def f3(quad1,t,quad2,r):
             return np.dot(self.p(quad1,t)-self.p(quad2,r),self.pprime(quad1,r))
        def f4(quad1,t,quad2,r):
            pdash = self.pprime(quad2,r)
            p1prime = self.pprime(quad1,t)
            n = np.array([pdash[1],-1*pdash[0]])
            return np.dot(p1prime,n)





        rs = []
        rs.append(0.0)
        rs.append(1.0)
        i = 2
        #fxi = [2.0,3.0]
        fxi = f1(quad1,t,quad2,rs[-1])
        fximinus1 =f1(quad1,t,quad2,rs[-2])
        while i<100:
            r = rs[-1]
            if abs(fxi[0]-fximinus1[0])<0.01 or (abs(fxi[0])<err   and r>=0 and r<=1 and fxi[-1]>0 and fxi[-2]>0):
                #print fxi[0]-fximinus1[0]
                break
            rn = rs[-1]-(fxi[0]*(rs[-1]-rs[-2]))/(fxi[0]-fximinus1[0])
            rs.append(rn)
            fximinus1 = fxi
            fxi = f1(quad1,t,quad2,rs[-1])
            i += 1
        rs1 = [];rs1.append(0.0);rs1.append(1.0)
        p1 = self.p(quad1,t)
        p1prime = self.pprime(quad1,t)
        n1 = np.array([p1prime[1],-1*p1prime[0]])

        dist1,p2i,r1 =self.distance(quad2,p1) #randomly ..no partcular thought
        costheta1 = np.dot(p2i-p1,n1)/(np.linalg.norm(n1)*dist)
        fpi = p1 + dist1/(2*costheta1)*n1/np.linalg.norm(n1)
        dist1,p2i,r2 =self.distance(quad2,fpi) #randomly ..no partcular thought


        gxi = self.f(quad1,t,quad2,rs1[-1])
        gximinus1 =self.f(quad1,t,quad2,rs1[-2])
        i =0
        while i<100:
            r1 = rs1[-1]
            print i,r1
            if abs(f4(quad1,t,quad2,r1))<0.001:
                print "first condition",abs(f2(quad1,t,quad2,r1)),abs(f3(quad1,t,quad2,r1))
                if abs(f2(quad1,t,quad2,r1))<0.01 and abs(f3(quad1,t,quad2,r1))<0.01:
                    print "second condition"
                    break
                else:
                    # how to generate next rn
                    print "parallel normals"
                    _,_,r1n =self.distance(quad2,self.p(quad1,t)) #randomly ..no partcular thought
                    rs1.append(r1n)

            else:
                fxi = f1(quad1,t,quad2,r1)
                if abs(gxi-gximinus1)<0.01:
                    if (abs(fxi[0])<0.01 and fxi[-1]>0 and fxi[-2]>0 and abs(gxi)<0.01   and r1>=0 and r1<=1):
                    #print fxi[0]-fximinus1[0]
                        break
                    else:
                        # get some  suggestion
                        r1n = random.uniform(0,1)
                        rs1.append(r1n)


                else:
                    r1n = rs1[-1]-(gxi*(rs1[-1]-rs1[-2]))/(gxi-gximinus1)
                    rs1.append(r1n)
                    gximinus1 = gxi
                    #gxi = self.f(quad1,t,quad2,rs1[-1])
                    gxi = self.f(quad1,t,quad2,rs1[-1])

            i += 1

        print fxi,r,gxi,r1
        #edges += [[fxi[1]],[fxi[2]],[fxi[3],fxi[4]],[fxi[3],fxi[3]+fxi[7]*fxi[5]],[fxi[4],fxi[4]+fxi[8]*fxi[6]],[fxi[1],fxi[2]]]
        fxi = f1(quad1,t,quad2,r1)
        edges += [[fxi[1]],[fxi[2]],[fxi[3],fxi[4]],[fxi[3],fxi[3]+fxi[7]*fxi[5]],[fxi[4],fxi[4]+fxi[8]*fxi[6]],[fxi[1],fxi[2]]]

        #initialr=-1,diagnostics=False
        #print "footprint,",fxi[0]
        if abs(fxi[0])<err and abs(gxi)<err and r>=0 and r<=1 and fxi[-1]>0 and fxi[-2] and abs(r-r1)<0.001>0:
            #return r,edges
            return edges,r,fxi[1]
        else:
            #print fxi[0]
            return edges,-1,fxi[1]

    def distance(self,quad,p0):#distance between a point and p
       if len(quad)==1: #single point
          return np.linalg.norm(p0-quad[0])
       t = 0.5
       pdprime = 2*(quad[2]+quad[0]-2*quad[1])

       for j in xrange(10): # maximum 100 iterations
            val = np.dot(self.pprime(quad,t),self.p(quad,t)-p0)
            if val>=0 and abs(val) <= 0.01:
                break
            else:
                fprime = np.linalg.norm(self.pprime(quad,t))**2 + np.dot(pdprime,self.p(quad,t)-p0)
                t = t - val/fprime
                if t<0 or t>1:
                    t = random.uniform(0,1)
       p1 = self.p(quad,t)
       p1prime = self.pprime(quad,t)
       #print "val is",val
       if np.cross(p0-p1,p1prime) >=0 and abs(val)<0.01: # point is on right side of the side
            return np.linalg.norm(p1-p0),p1,t
       else:
           return float('inf'),p1,t
    def distanceIt(self,ind,p0):#distance between a point and p
           quad = self.edges[ind]
           dist,p1,t = self.distance(quad,p0)
           if self.its[ind].search(t)!=set():
               return dist,p1,t
           else:
               return float("inf"),[],-1
    def distanceItn(self,ind,p0):
        quad = self.edges[ind]
        dist,p1,t = self.distance(quad,p0)
        if self.its[ind].search(t)!=set():
            return dist,p1,ind,t
        else:
            return float("inf"),[],-1,-1




    def friends(self,e,cnts,convex_vert):
        return [e]
        frds = []
        init = self.nextedge(e,cnts)
        frds.append(e)
        while not init in convex_vert and init !=e:
            frds.append(init)
            init = self.nextedge(init,cnts)
        init = self.prevedge(e,cnts)
        while not init in convex_vert and init!=e:
            frds.append(init)
            init = self.prevedge(init,cnts)
        return frds

    def manage_interval(self,mplast,mpcurr):
        # TODO:
        first ={}
        first[mplast[1]] = mplast[2]
        first[mplast[3]] = mplast[4]
        second = {}
        second[mpcurr[1]] = mpcurr[2]
        second[mpcurr[3]] = mpcurr[4]
        e11 = mplast[1];e12=mpcurr[1]
        e21= mplast[3];e22 =mpcurr[3]
        e1 = e11
        while (e1!=e12):
            if e1 in first:
                self.its[e1].chop(first[e1],self.its[e1].end())
            else:
                self.its[e1].chop(self.its[e1].begin(),self.its[e1].end())

            e1 = self.nextedge(e1)
        if e12 in first:
            inite1 = first[e12]
        else:
            inite1 = self.its[e12].begin()
        self.its[e12].chop(inite1,second[e12])
        e2 = e21
        while e2!=e21:
            if e2 in first:
                self.its[e1].chop(self.its[e2].begin(),first[e2])
            else:
                self.its[e1].chop(self.its[e2].begin(),self.its[e2].end())
            e2 = self.prevedge(e2)
        if e22 in first:
            ende2 = first[e22]
        else:
            ende2 = self.its[e22].end()

        self.its[e22].chop(second[e22],ende2)



    def insertReflex(self):#
        concave_vert = copy.copy(self.concave_vert)
        concave_vert.sort()
        #print cnts
        cnts1 = copy.copy(self.cnts)
        while len(concave_vert) >0:
            #print concave_vert
            i = concave_vert.pop(-1)
            #print i
            e2 = i
            e1 = self.prevedge(i)
            e1prime = self.p(self.edges[e1],0.7)
            e2prime = self.p(self.edges[e2],0.3)
            quad = []
            quad.append(e1prime)
            #print e1prime
            quad.append(self.vertices[i])
            quad.append(e2prime)
            #print e1,e2
            #lastts[e1][1] = 0.999
            #edges[e1][2] = e1prime
            #edges[e2][0] = e2prime
            self.vertices.insert(i,e1prime)
            self.vertices[i+1] = e2prime
            #lastts[i+1] = [0.001,1.0]

            self.edges.insert(i,quad)
            ind1 = 0
            ind3 = 0
            #print "before adding",i
            #print its

            for ind2,cnt in enumerate(self.cnts):
                for ind,k in enumerate(cnt):
                    if k > i:
                        cnt[ind] = k+1
                    elif k ==i:
                        ind1 = ind
                        ind3 = ind2
            for k in  xrange(len(self.edges)-1,i,-1):
                if k >i:

                    self.its[k+1] = self.its[k]
            self.its[i]=iv(0.0,1.001)
            lowerbound = self.its[e1].begin
            self.its[e1] = iv(lowerbound,0.7)
            upperbound = self.its[i+1].end
            self.its[i+1] =iv(0.3,upperbound)

            self.cnts[ind3].insert(ind1+1,i+1)
            #print "after adding",i
            #print its

            #for j,k in enumerate(concave_vert):
            #    if k>i:
            #        concave_vert[j] = k+1


        #print cnts
        #print its
        return self.its
    def nextt(self,e1,t,e2,stepsize=0.05):
        if self.its[e1].search(t+stepsize) != set():
            return e1,t+stepsize
        else:
            preve1 = e1
            e1 = self.nextedge(e1)
            #tracededges.append(preve1)
            if e1!=e2:
                #self.remove_interval(its,preve1,t,its[preve1].end())
                return e1,self.its[e1].begin()
            else:
                return -1,-1
    def prevr(self,e2,r,e1,stepsize=0.05):
        print "in prevr",e2,r,self.its[e2],r-stepsize
        if self.its[e2].search(r-stepsize) != set():
            return e2,r-stepsize
        else:
            preve2 = e2
            e2 = self.prevedge(e2)
            #tracededges.append(preve1)
            if e2!=e1:
                #self.remove_interval(its,preve2,its[preve2].begin(),r)

                #print "returning e2,",e2,self.its[e2].end()-0.001
                return e2,self.its[e2].end()-0.001
            else:
                return -1,-1
    def updateT(self,e1prev,tprev,e2prev,rprev,stepsize=0.05):
        # try using nextt,prevr if r fails, try using prevr,nextt, if both fails
        # then stop and return
        e1next,tnext = self.nextt(e1prev,tprev,e2prev,stepsize)
        if e1next == -1:
            return -1,-1,-1,-1,[]
        #print e1next,tnext

        _,rnext,mp,x = self.footprintIt(e1next,tnext,e2prev)
        if rnext == -1:
            #print "failed ",e1next,tnext,e2prev
            e2next,rnext= self.prevr(e2prev,rprev,e1prev,stepsize)
            #print "now trying ",e2next,rnext,e1prev
            if e2next ==-1:
            #    print "this also failed as hit end of contour"
                return -1,-1,-1,-1,[]

            _,tnextn,mp,x = self.footprintIt(e2next,rnext,e1prev)
            #print e2next,rnext,e1prev,tnextn
            if tnextn==-1:
            #    print "this also failed"
                if e1next==e1prev:
                    e1next = self.nextedge(e1next)
            #    print "now let us try ",e2next,rnext,e1next
                #try
                #e1next,tnext = self.nextt(e1next,tnext,its,cnts,e2next,stepsize)
                #print "try nextt",e1next,tnext," with e2 as ",e2next
                _,tnextn,mp,x = self.footprintIt(e2next,rnext,e1next)
                if tnextn ==-1:
            #        print "this also failed, so let us try ",e1next,tnext,e2next
                    _,rnext,mp,x = self.footprintIt(e1next,tnext,e2next)
                    if rnext == -1:
            #            print "this also failed so let us give up"
                        return -1,-1,-1,-1,[]
                    else:
                        return e1next,tnext,e2next,rnext,mp
                else:
                    return e1next,tnextn,e2next,rnext,mp


            else:
                return e1prev,tnextn,e2next,rnext,mp
        else:
            return e1next,tnext,e2prev,rnext,mp
    def iterate(self,e1,t,e2,r,mplast,points,stepsize=0.05,counts=10):
        if e1==e2:
            return
        print "beginning",e1,t,e2,r,len(points),counts
        e1in,tin,e2in,rin = e1,t,e2,r
        count = -1
        noofsteps = 0
        #mplastn = self.footprintIt(edges[e1],t,edges[e2],its[e1],its[e2],r)
        #mplast = [mplastn,e1,t,e2,r]
        while count < counts:
            count += 1
            noofsteps += 1
            e1prev,tprev,e2prev,rprev = e1,t,e2,r
            e1,t,e2,r,mp = self.updateT(e1,t,e2,r,stepsize)
            print count,e1,t,e2,r,mp
            if e1 == -1 or t==-1 or r==-1:
                break
            dist = np.linalg.norm(mp-self.p(self.edges[e1],t))
            #points.append([mp])
            points.append([mp,e1,t,e2,r,dist])
            '''
            R1 = self.radiusCurvature(edges[e1],t)
            R2 = self.radiusCurvature(edges[e2],r)

            if dist > R1:
                print "curvature failure at ",e1," dist=",dist,' R=',R1
            if dist > R2:
                print "curvature failure at ",e2," dist=",dist,' R=',R2
            '''
            mpcurr = [mp,e1,t,e2,r,dist]
            if noofsteps>2 :
                print "now doing distance check"
                pe,mindist,mink,dist = self.distCheck(e1,t,e2,r,mp)
                noofsteps = 0

                #print pe
                try:
                    print mplast
                    e1prev,tprev,e2prev,rprev = mplast[1:5]


                    if len(pe) >0 :
                        #mpcurr = [mp,e1,t,e2,r,dist]
                        print "dist failure with ",pe
                        #while noofsteps>1:
                        #for popcount in xrange(3):
                        #    points.pop(-1)
                        negdist = False
                        for (_,dist0,r) in pe:
                            if dist0 <0:
                                negdist = True
                                break
                        if negdist:
                            #mplast = mpcurr
                            #print mplast
                            continue

                        #print pes
                        #print its[mpcurr[1]],its[mpcurr[3]]
                        if mplast[1]!= mpcurr[1]:#e1 has changed
                            k = -2
                            curre = mpcurr[1]
                            while points[k][1]!=mplast[1]:
                                print k
                                if points[k][1]!=curre:
                                    curre = points[k][1]
                                    mpcurr = points[k]

                                mpcurr1 = points[k]
                                pe1,mindist1,mink1,dist1 = self.distCheck(mpcurr1[1],mpcurr1[2],mpcurr1[3],mpcurr1[4],mpcurr1[0])
                                if len(pe1)==0:
                                    mplast = mpcurr1
                                    break
                                k = k-1
                        mps = self.retract(mpcurr,mplast,mink,mindist)
                        mpfn = mps.pop(-1)
                        print "retracted..first",mps
                        e1,t,e2,r = mpfn[1:5]
                        #points.append([mpfn[0]])
                        points.append(mpfn)
                        for k,mp1 in enumerate(mps):
                            #if mp[1] in convex_vert or e3!=self.nextedge(mpfn[1],cnts):
                            print mp1
                            self.iterate(mp1[1],mp1[2],mp1[3],mp1[4],mp1,points,counts=counts)
                        if len(mps)>0:
                            break
                    else:
                        self.manage_interval(mplast,mpcurr)
                        '''
                        e1,t,e2,r = mpcurr[1:5]

                        if e1 == e1prev:
                            self.remove_interval(its,e1,tprev,t)
                        else:
                            self.remove_interval(its,e1prev,tprev,its[e1prev].end())
                            self.remove_interval(its,e1,its[e1].begin(),t)

                        if e2 == e2prev:
                            self.remove_interval(its,e2,r,rprev)
                        else:
                            self.remove_interval(its,e2prev,its[e2prev].begin(),rprev)
                            self.remove_interval(its,e2,r,its[e2].end())
                        '''

                except Exception:
                        traceback.print_exc()
                        return

                        #get next point and then retract
                        dist3 = float("Inf")
                        for e3 in pe:
                            dist3 = min(dist3,e3[1])
                        e1n,tn,e2n,rn,mpn = self.updateT(e1,t,e2,r,its,cnts,edges,stepsize)
                        #mpcurr = mplast
                        distn = np.linalg.norm(mpn-self.p(edges[e1n],tn))
                        mplast = [mpn,e1n,tn,e2n,rn,distn]
                        count = 0

                        while dist3 < distn and count<100 and distn <dist:
                            e1n,tn,e2n,rn,mpn = self.updateT(e1n,tn,e2n,rn,its,cnts,edges,stepsize=0.05)
                            mpcurr = mplast
                            print e1n,tn,e2n,rn,mpn
                            distn = np.linalg.norm(mpn-self.p(edges[e1n],tn))
                            dist3 = float('inf')
                            for e3 in pe:
                                distv,_,_ = self.distanceIt(e3[0],mpn,its,edges)
                                dist3 = min(dist3,distv[0])


                            mplast = [mpn,e1n,tn,e2n,rn,distn]
                            print distn,dist3
                            count += 1


                        print "now retracting..second"
                        mpfn,pes = self.retract(mpcurr,mplast,edges,its,pe,cnts,convex_vert)
                        print mpfn,pes
                        e1,t,e2,r = mpfn[1:5]

                        points.append([mpfn[0]])
                        for e3 in pes:
                            self.iterate(mpfn[1],mpfn[2],e3,pes[e3][2],its,cnts,edges,points,convex_vert,stepsize,counts=counts)
                            self.iterate(e3,pes[e3][2],mpfn[3],mpfn[4],its,cnts,edges,points,convex_vert,stepsize,counts=counts)
                        break



                mplast = mpcurr




        return points
    def nextfp(self,mp,stepsize):
        e1 = mp[1];t0=mp[2]
        e2 = mp[3];r0=mp[4]
        # first increment e1 by stepsize
        t1 = min(self.its[e1].end(),t0+stepsize)
        fp = []
        dist = float('inf')
        _,r1,fp= self.footprint(self.edges[e1],t1,self.edges[e2])
        if r1 == -1:
            r1 = max(r0-stepsize,self.its[e2].begin())
            _,t1,fp= self.footprint(self.edges[e2],r1,self.edges[e1])
            if t1!=-1:
                dist = np.linalg.norm(self.p(self.edges[e2],r1)-fp)
        else:
            dist = np.linalg.norm(self.p(self.edges[e2],r1)-fp)
        return [fp,e1,t1,e2,r1,dist]
    def distCheckn(self,points):
        if len(points) <2:
            return []
        # returns a set of footprints, which have equal distance
        dist0 = points[-1][-1]
        mp = points[-1][0]
        pe = []
        mindist = float('inf')
        mps = []
        e1 = points[-1][1];e2 = points[-1][3]
        for k,_ in enumerate(self.edges):
            if k==e1 or k==e2:
                continue
            dist,p1,ind,t = self.distanceItn(k,mp)
            if dist < dist0:
                pe.append(k)
                if dist < mindist:
                    mink = k
                    mindist = dist
        if len(pe) >0:
            print "dist failure with ",pe
            mpl = points.pop(-1)
            for k in xrange(self.maxcount-1):
                points.pop(-1)
            points.append(mpl)

            mps=self.retractn(points,mindist,mink,pe)
        return mps
    def retractn(self,points,mindist,mink,pe):
        if len(points)<2:
            return []
        mps = []
        e1 = points[-2][1]
        tf = points[-2][2];rf = points[-2][4]
        tl = points[-1][2];rl = points[-1][4]
        e2 = points[-2][3]
        dist = points[-1][-1]
        assert e2 == points[-1][3] and e1 == points[-1][1]
        count = 0
        mindistn = mindist;distn = dist;fp = points[-1]
        tn = tl
        #print dist,mindist
        while abs(mindistn-distn)>2.0  and count<20:
            #print tf,tn,tl,distn,mindistn
            count += 1
            tn = (tf + tl)/2
            _,rn,fp= self.footprint(self.edges[e1],tn,self.edges[e2])
            if rn == -1:
                rn = (rf + rl)/2
                _,tn,fp= self.footprint(self.edges[e2],rn,self.edges[e1])
            assert tn!=-1 and rn!=-1
            mindistn,p1,e3,te3 = self.distanceItn(mink,fp)
            distn = np.linalg.norm(self.p(self.edges[e1],tn)-fp)
            if mindistn < distn:
                tl = tn;rl = rn
            else:
                tf = tn;rf = rn
        mpn = [fp,e1,tn,e2,rn,distn]
        #print distn,mindistn,mink
        if abs(mindistn-distn)<2.0:
            points.pop(-1)
            mp1 = [fp,e1,tn,mink,te3,distn]
            mp2 = [fp,e2,rn,mink,te3,distn]
            points.append(mpn)
            mps.append(mp1);mps.append(mp2)
            #check if other sides in pe also have same distn
            for k in pe:
                if k!= mink:
                    dist0,p1,k,newt = self.distanceItn(k,fp)
                    if abs(dist0-mindistn)<2.0:
                        mps.append([fp,e1,tn,k,newt,distn])
                        mps.append([fp,e2,rn,k,newt,distn])
            return mps
        elif mindistn <distn:
            points.pop(-1)
            return [points[-1]]
        else:
            points.pop(-1)
            points.append(mpn)
            return [mpn]
    def trace(self,mplast,stepsize=0.05):
        # mplast shall have two footprints
        # get next mp
        points = [mplast]
        breakpoint = False
        endofedge = False
        count = 0
        e1 = mplast[1];e2 = mplast[3];mpn = mplast
        print "tracing ",e1,e2
        mps = []
        if self.its[e1].end()-self.its[e1].begin() <= 0.001:
            breakpoint = True
            #self.tracededges.append(e1);self.untracededges.remove(e1)
        if self.its[e2].end()-self.its[e2].begin() <= 0.001:
            breakpoint = True
            #self.tracededges.append(e2);self.untracededges.remove(e2)
        # distance check
        # if fail-stop
        # if convex_vert stop
        # if curvature failure stop
        # if edges change stop
        while not breakpoint:
            mpn = self.nextfp(mpn,stepsize)
            if mpn[2]!=-1 and mpn[4]!=-1:
            #print mpn
                points.append(mpn)
            else:
                break
            if mpn[2]==self.its[e1].end() or mpn[4] == self.its[e2].begin():
                breakpoint = True
            print len(points),e1,e2,count,self.maxcount
            count += 1
            if count == self.maxcount:
                count = 0
                mps = self.distCheckn(points)
                #print mps
                if len(mps)>0:
                    breakpoint = True
                if len(mps) == 1:
                    mps = []
                    breakpoint = False
        mpn=points[-1]
        r0 = mpn[4];t0=mpn[2]
        if self.its[e1].end()-t0 <= 0.001:
            endofedge = True
            #self.tracededges.append(e1);self.untracededges.remove(e1)
            e1 = self.nextedge(e1,self.cnts);t0 = 0.0
        if r0-self.its[e2].begin() <= 0.001:
            breakpoint=True
            endofedge=True
            e2 = self.prevedge(e2,self.cnts);r0=1.0
            #self.tracededges.append(e2);self.untracededges.remove(e2)
        '''
        if endofedge == True:
            _,r0,fp,_ = self.footprintItn(e1,t0,e2)
            if r0==-1:
                _,t0,fp,_ = self.footprintItn(e2,r0,e1)
            if r0!=-1 and t0!=-1:
                mps.append([fp,e1,t0,e2,r0,np.linalg.norm(fp-self.p(self.edges[e1],t0))])
        '''
        if len(points) >0:
            rl = mplast[4];tf = mplast[2]
            rf = points[-1][4];tl = points[-1][2]
            self.its[e1].chop(tf,tl);self.its[e2].chop(rf,rl)
        return mps,points
    def test(self):
        mps = []
        points = []
        for cv in self.convex_vert:
            mps.append([self.vertices[cv],cv,0.0,self.prevedge(cv,self.cnts),1.0,0.0])
        while len(mps) > 0:
            mp = mps.pop(-1)
            mpsn,pts = self.trace(mp)
            points = points + pts
            mps = mps + mpsn
        eds = []
        for point in points:
            eds.append([point[0]])
        return eds,self.edges
    def test1(self,points):
        mps = []
        for cv in self.convex_vert:
            mps.append([self.vertices[cv],cv,0.0,self.prevedge(cv),1.0,0.0])
        while len(mps) > 0:
            mp = mps.pop(-1)
            pts = self.iterate(mp[1],mp[2],mp[3],mp[4],mp,points,0.001,500)
            print len(points)

        return points,self.edges
    def distCheck(self,e1,t,e2,r,mp):

        #r is not used at all
        dist = np.linalg.norm(self.p(self.edges[e1],t)-mp)
        pe = []
        mindist3 = float('inf')
        mink = -1
        for k,_ in enumerate(self.edges):
            #if k in self.friends(e1,cnts,convex_vert) or k in self.friends(e2,cnts,convex_vert):
            if k==e1 or k==e2:
            #or k==self.nextedge(e1,cnts) or k==self.prevedge(e2,cnts) or k==self.prevedge(e1,cnts) or k==self.prevedge(e2,cnts):
            # or k in tracededges:
                continue
            dist3,p,t = self.distanceIt(k,mp)
            #dist3,p,t = self.distance(self.edges[k],mp)


            #print k,dist3,dist,t
            if dist3 <dist:
                pe.append([k,dist3,t])
                if dist3<mindist3:
                    mink = k
                    mindist3 = dist3
        return pe,mindist3,mink,dist
    def retract(self,mpcurr,mplast,mink,mindist):
           print "retract begin mplast,mpcurr ",mplast,mpcurr
           count = -1
           pes ={} # to store footprints
           e1 = mpcurr[1]
           tl = mpcurr[2]
           e2 = mpcurr[3]
           rl = mpcurr[4]
           e11 = mplast[1];e21=mplast[3]
           dist= mpcurr[-1]
           tn = tl
           rn = rl
           distk,_,_ = self.distance(self.edges[mink],mpcurr[0])
           tf = mplast[2]
           mpn = mplast[0]
           while abs(distk-dist)>2.0  and count < 20:
                count += 1
                tn = (tf + tl)/2.0
                _,rn,mpn,x = self.footprintIt(e1,tn,e2)

                #print "rn =",rn
                if rn == -1:
                    #check if there are two values of e2
                    if e2!=e21:
                        _,rn,mpn,x = self.footprintIt(e1,tn,e21)
                    else:# it should not happen
                        print "possible bug in retract or footprint ",e1,tn,e2,mplast,mpcurr
                        #return tn,rn,-1
                        break
                else:
                    dist = np.linalg.norm(mpn-self.p(self.edges[e1],tn))
                    distk,_,_ = self.distance(self.edges[mink],mpn)
                    if distk < dist:
                        tl = tn
                    else:
                        tf = tn

           if rn!=-1:
                mpfn = [mpn,e1,tn,e2,rn,dist]
           else:
               mpfn = mpcurr
           if abs(distk-dist)>2.0:
               print "possible bug in retraction, new dist difference is ",abs(distk-dist),mpcurr,mplast
           mps = []
           #mps.append([mpfn[0],mink,])
           for i,_ in enumerate(self.edges):
               if i==mpfn[1] or i==mpfn[3]:
                   continue
               disti,fpi,ti = self.distanceIt(i,mpn)
               if abs(disti-distk)<2.0:
                   mpfn1 = [mpn,i,ti,mpfn[1],mpfn[2],disti]
                   mpfn2 = [mpn,i,ti,mpfn[3],mpfn[4],disti]
                   mps.append(mpfn1);mps.append(mpfn2)



               #return tn,rn,-1
           mps.append(mpfn)
           self.manage_interval(mplast,mpfn)
           return mps
    def modEdges(self,letter):
        self.edges,self.vertices,self.vert_edges,self.cnts = check.fontOutline(letter)
        self.vertAnalyse()
        self.its = {} # store as interval tree
        for k in xrange(len(self.edges)+len(self.concave_vert)):
            self.its[k]=iv(0.0,1.001)
        #its =
        self.insertReflex()
        for k in xrange(len(self.edges)):
            self.its[k] = it([self.its[k]])
        self.vertAnalyse()
        return self.edges,self.its,self.cnts,self.vertices,self.convex_vert,self.concave_vert
    def remove_interval(self,its,e1,t0,tn):
        its[e1].chop(t0,tn)
    #def manage_interval(self,e11,tf,)
    def footprintn(self,quad1,t,quad2,initialr=-1,diagnostics=False):
        if initialr == -1:
            r = 0.5
        else:
            r = initialr
        edges = []
        p1 = self.p(quad1,t)
        p1prime = self.pprime(quad1,t)
        normal1 = np.array([p1prime[1],-1*p1prime[0]])

        edges.append([p1])
        edges.append([p1,p1+p1prime*50/np.linalg.norm(p1prime)])
        #print np.linalg.norm(p1prime)


        for count in xrange(100):
                p2 = self.p(quad2,r)
                p2prime = self.pprime(quad2,r)
                normal2 = np.array([p2prime[1],-1*p2prime[0]])
                edges.append([p2])
                edges.append([p2,p2+normal2*50/np.linalg.norm(normal2)])

                denominator = np.dot(normal1,p2prime)
                if abs(denominator) > 0.01 :
                    alpha = np.dot(p2-p1,p2prime)/denominator
                    beta = np.dot(p2-p1,p1prime)/denominator
                else: #parallel found
                    alpha = beta = np.linalg.norm(p1-p2)/2.0
                    alpha = alpha/np.linalg.norm(p1prime)
                    beta = beta/np.linalg.norm(p2prime)
                r,result = self.solve1(quad1,t,quad2,r,1.0)
                #print r
                if result:
                    break
        if result and alpha>=0 and beta>=0 :
            edges.append([p2,p2+normal2*50/np.linalg.norm(normal2)])
            edges.append([p1,p1 + alpha*normal1])
            edges.append([p2,p2 + beta*normal2])



            return edges,r,p1+alpha*normal1
        else:

            return edges,-1,None



    def footprinto(self,quad1,t,quad2,initialr=-1,diagnostics=False):
        edges = []
        p1 = self.p(quad1,t)
        p1prime = self.pprime(quad1,t)
        normal1 = np.array([p1prime[1],-1*p1prime[0]])

        #edges.append([p1])
        #edges.append([p1,p1+p1prime*50/np.linalg.norm(p1prime)])
        #print np.linalg.norm(p1prime)


        if initialr <0:
            dist,_,r = self.distance(quad2,p1)
            p2 = self.p(quad2,r)

            p2prime = self.pprime(quad2,r)
            denominator = np.dot(normal1,p2prime)
            print "initial denominator is ",denominator
            if abs(denominator)>0.01:# not horizontal
               #fp = p1 + normal1/np.linalg.norm(normal1)*dist
               #distn,_,rn = self.distance(quad2,fp)
               r = 1-r



        else:
            r = initialr
        print "initial values",r,dist
        if r==-1:
            r = 0.5
        result = False
        for count in xrange(100):
            if  r != -1 and r>=0.0 and r<=1.0:
                p2 = self.p(quad2,r)
                p2prime = self.pprime(quad2,r)
                normal2 = np.array([p2prime[1],-1*p2prime[0]])
                #edges.append([p2])
                #edges.append([p2,p2+normal2*50/np.linalg.norm(normal2)])

                denominator = np.dot(normal1,p2prime)
                print "denominator is ",denominator
                if abs(denominator) > 0.01 :
                    alpha = np.dot(p2-p1,p2prime)/denominator
                    beta = np.dot(p2-p1,p1prime)/denominator

                    #print alpha>0,beta>0



                else: #parallel found
                    print "parallels"
                    alpha = beta = np.linalg.norm(p1-p2)/2.0
                    alpha = alpha/np.linalg.norm(p1prime)
                    beta = beta/np.linalg.norm(p2prime)
                    #print alpha,beta




                r,result = self.solve(quad1,t,quad2,r)
                print r
                if result:
                    break
        if result and alpha>=0 and beta>=0 :
            print alpha,beta
            #edges.append([p2,p2+normal2*50/np.linalg.norm(normal2)])
            edges.append([p1 + alpha*normal1])
            edges.append([p2 + beta*normal2])

            edges.append([p1,p1 + alpha*normal1])
            edges.append([p2,p2 + beta*normal2])



            if diagnostics:
                return edges,r,p2 + beta*normal2,t,r
            else:
                return edges,r,p2+beta*normal2
        else:
            if diagnostics:
                return edges,-1,None,t,r
            else:
                return edges,-1,None
    def footprintItn(self,e1,t,e2,initialr=-1):#using interval tree
        if( self.its[e1].search(t) == set()):
            return [],-1,[],None
        e,r,mp = self.footprint(self.edges[e1],t,self.edges[e2],initialr)
        x = it2.search(r)
        if x!=set():
            assert len(x) == 1 # if not then we have not managed intervals better
            return e,r,mp,x
        else:
            return e,-1,[],None

    def footprintIt(self,e1,t,e2,initialr=-1):#using interval tree
        quad1=self.edges[e1];quad2=self.edges[e2];it1=self.its[e1];it2=self.its[e2]
        if( it1.search(t) == set()):
            return [],-1,[],None

        e,r,mp = self.footprint(quad1,t,quad2,initialr)
        #print "line 1057:r is",r
        x = it2.search(r)
        #y = it2.search(r)
        if x!=set():
            assert len(x) == 1 # if not then we have not managed intervals better
            return e,r,mp,x
        else:
            return e,-1,[],None


    def bezIntersect(self,bez1,bez2):
        return
    def bisector(self,quad1,quad2,counts=100):
        finedges = []
        points = []
        for q in xrange(counts+1):
            t = q/(counts*1.0)
            edges,r,point = self.footprint(quad1,t,quad2)
            #r,edges = self.solve2(quad1,t,quad2)
            #point = edges[0][0]
            #print "t=,r=",t,r


            if r>=0 and r<=1.0:
                print "t=,r=",t,r


                finedges +=  edges
                points.append(point)

        #print points
        if len(points)>2:
            beziers = fitCurve(points,10.0)
            print beziers
            return  finedges,beziers
        else:
            return finedges,[[points]]
    def findNewMP(self,edges,firstedge,secondedge,thirdedge,lastt,currentt,dist):
       for j in xrange(20): # max 20 iterations
           newt = (lastt + currentt)/2.0
           _,r,mp = self.footprint(edges[firstedge],newt,edges[secondedge])
           if r == -1:
               lasst = newt
               continue
           newdist,_,_ = self.distance(edges[thirdedge],mp)
           if abs(newdist-dist)<0.1:
               break
           lastt = newt



       return mp,newt
    def nextedge(self,e): #gives the next edge given an index of edge
       for cnt in self.cnts:
           for i,pt in enumerate(cnt):
               if pt == e:
                   if i <len(cnt)-1:
                       return cnt[i + 1]
                   else:
                       return cnt[0]
    def prevedge(self,e): #gives the next edge given an index of edge
       for cnt in self.cnts:
           for i,c in enumerate(cnt):
               if c == e:
                   return cnt[i-1]

    def medaxn():
        '''
          steps:
          1. initialize all edges as interval tree with initial interval 0.0 1.0
          2. start tracing e1,e2
          3. At breakpoint, delete traced intervals from e1 and e2
          4. Trace next pair
          Trace Pair e1,e2
            repeat  till success
              try first point of first interval of e11, e1 - footprint e2x
                 if success record
                     try last point of the segment
                     if failure then overlapping parts are,  e2 (e2x, last point of the interval) e1 (e11,footprint of lastpoint of interval)
                     else
                       overlapping portions are e1(first interval)
                 if failure
                      try last part of segment
                      if sucess e12 then overlapping parts are e2(beginning of this interval(x), fp(e2)) and e1(fp(x),end of interval)
            if success:
                then start for three intervals and then do distance check and then retract
                store these points and convert as bezier curve and store
        '''
        return
    def printMa(self,ma):
        points = []
        for k,ed in enumerate(ma):
            for point in ed:
                if point != None:
                    if len(point)>0:
                        points.append([point[0]])
        return points
    def vertAnalyse(self):
        self.concave_vert = []
        self.convex_vert = []
        print len(self.vertices),len(self.edges)
        for i,_ in enumerate(self.vertices):
            j = self.prevedge(i)
            cross = np.cross(self.pprime(self.edges[i],0.0),self.pprime(self.edges[j],1.0))
            #print i,j,cross

            if cross>0.01:
                #print "convex"
                self.convex_vert.append(i)
            elif cross<-0.01:
                #print "concave"
                self.concave_vert.append(i)
        return self.convex_vert,self.concave_vert


    def radiusCurvature(self,quad,t):
         #R=((x^('2)+y^('2))^(3/2))/(|x^'y^('')-y^'x^('')|)
         # R = norm(pprime)^3/norm(dot(pprime,normal of pdoubleprime)
         pprime = self.pprime(quad,t)
         pdprime = 2*(quad[2]+quad[0]-2*quad[1])
         pdprimen = np.array([pdprime[1],-1*pdprime[0]])
         if np.linalg.norm(pdprime ==0):
             R = float('inf')
         else:
             R = math.pow(np.linalg.norm(pprime),3)/np.linalg.norm(np.dot(pprime,pdprimen))
         return R
    def strokes(self,edges,vertices,vert_edges):
        # edges are identfied by vertex index, vertex index+1 ...0 is also len(edges)
        # there is another in which vertices are mapped to edges they contain
        # take first convex vertex
        # find edge e1 = 0-1 and e2=25-0 then e2 is traced 0-25, e1 must be straight and e2 must be opposite
        # find footpoint c1 on e1 and correspoind point c2 on e2 and mp
        # radius = norm(mp- c2)
        # radius > R at c2 #curvature criterion is satisfied
               # handle curvature failure --create a branch point and edge pair is AD1 and D1B of same
               # edge AB first MAT is centre of curvature and it meets the last point traced
        # distance check with other edges which are not present edges and which have not been traced
        # if dist > radius :
              # handle branch point (edge1,edge2,edge3)
              # trace with edge1,edge2 and edge1 and edge3

        #----trace u1 till one gets a convex vertex, reflex vertex or a branch point




        return 1
