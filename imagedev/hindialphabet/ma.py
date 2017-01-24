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
                return random.uniform(0,1),False #a random number between 0 and 1
            #print r
        else: # not parallel normal
            fp = self.fprime(quad1,t,quad2,r)
            if fp != 0:
                val1 = self.f(quad1,t,quad2,r)/fp
            else:
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
    def distance(self,quad,p0):#distance between a point and p
       if len(quad)==1: #single point
          return np.linalg.norm(p0-quad[0])
       t = 0.5
       pdprime = 2*(quad[2]+quad[0]-2*quad[1])

       for j in xrange(100): # maximum 100 iterations
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
       if np.cross(p0-p1,p1prime) >=0: # point is on right side of the side
            return np.linalg.norm(p1-p0),p1,t
       else:
           return float('inf'),p1,t
    def distanceIt(self,ind,p0,its,edges):#distance between a point and p
           quad = edges[ind]
           dist,p1,t = self.distance(quad,p0)
           if its[ind].search(t)!=set():
               return dist,p1,t
           else:
               return float("inf"),[],-1


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


    def insertReflex(self,edges,vertices,cnts,concave_vert,its):#
        concave_vert = copy.copy(concave_vert)
        concave_vert.sort()
        #print cnts
        cnts1 = copy.copy(cnts)
        while len(concave_vert) >0:
            #print concave_vert
            i = concave_vert.pop(-1)
            #print i
            e2 = i
            e1 = self.prevedge(i,cnts)
            e1prime = self.p(edges[e1],0.9)
            e2prime = self.p(edges[e2],0.1)
            quad = []
            quad.append(e1prime)
            #print e1prime
            quad.append(vertices[i])
            quad.append(e2prime)
            #print e1,e2
            #lastts[e1][1] = 0.999
            #edges[e1][2] = e1prime
            #edges[e2][0] = e2prime
            vertices.insert(i,e1prime)
            vertices[i+1] = e2prime
            #lastts[i+1] = [0.001,1.0]

            edges.insert(i,quad)
            ind1 = 0
            ind3 = 0
            #print "before adding",i
            #print its

            for ind2,cnt in enumerate(cnts):
                for ind,k in enumerate(cnt):
                    if k > i:
                        cnt[ind] = k+1
                    elif k ==i:
                        ind1 = ind
                        ind3 = ind2
            for k in  xrange(len(edges)-1,i,-1):
                if k >i:

                    its[k+1] = its[k]
            its[i]=iv(0.01,1.001)
            lowerbound = its[e1].begin
            its[e1] = iv(lowerbound,0.9)
            upperbound = its[i+1].end
            its[i+1] =iv(0.1,upperbound)

            cnts[ind3].insert(ind1+1,i+1)
            #print "after adding",i
            #print its

            for j,k in enumerate(concave_vert):
                if k>i:
                    concave_vert[j] = k+1


        #print cnts
        #print its
        return its
    def nextt(self,e1,t,its,cnts,e2,stepsize=0.05):
        if its[e1].search(t+stepsize) != set():
            return e1,t+stepsize
        else:
            preve1 = e1
            e1 = self.nextedge(e1,cnts)
            #tracededges.append(preve1)
            if e1!=e2:
                self.remove_interval(its,preve1,t,its[preve1].end())
                return e1,its[e1].begin()
            else:
                return -1,-1
    def prevr(self,e2,r,its,cnts,e1,stepsize=0.05):
        print "in prevr",e2,r,its[e2],r-stepsize
        if its[e2].search(r-stepsize) != set():
            return e2,r-stepsize
        else:
            preve2 = e2
            e2 = self.prevedge(e2,cnts)
            #tracededges.append(preve1)
            if e2!=e1:
                self.remove_interval(its,preve2,its[preve2].begin(),r)

                print "returning e2,",e2,its[e2].end()-0.001
                return e2,its[e2].end()-0.001
            else:
                return -1,-1
    def updateT(self,e1prev,tprev,e2prev,rprev,its,cnts,edges,stepsize=0.05):
        # try using nextt,prevr if r fails, try using prevr,nextt, if both fails
        # then stop and return
        e1next,tnext = self.nextt(e1prev,tprev,its,cnts,e2prev,stepsize)
        if e1next == -1:
            return -1,-1,-1,-1,[]
        #print e1next,tnext

        _,rnext,mp,x = self.footprintIt(edges[e1next],tnext,edges[e2prev],its[e1next],its[e2prev],-1)
        if rnext == -1:
            print "failed ",e1next,tnext,e2prev
            e2next,rnext= self.prevr(e2prev,rprev,its,cnts,e1prev,stepsize)
            print "now trying ",e2next,rnext,e1prev
            if e2next ==-1:
                print "this also failed as hit end of contour"
                return -1,-1,-1,-1,[]

            _,tnextn,mp,x,tfin,rfin = self.footprintIt(edges[e2next],rnext,edges[e1prev],its[e2next],its[e1prev],-1,True)
            print e2next,rnext,e1prev,tnextn,tfin,rfin
            if tnextn==-1:
                print "this also failed"
                if e1next==e1prev:
                    e1next = self.nextedge(e1next,cnts)
                print "now let us try ",e2next,rnext,e1next
                #try
                #e1next,tnext = self.nextt(e1next,tnext,its,cnts,e2next,stepsize)
                #print "try nextt",e1next,tnext," with e2 as ",e2next
                _,tnextn,mp,x,tfin,rfin = self.footprintIt(edges[e2next],rnext,edges[e1next],its[e2next],its[e1next],-1,True)
                if tnextn ==-1:
                    print "this also failed, so let us try ",e1next,tnext,e2next
                    _,rnext,mp,x = self.footprintIt(edges[e1next],tnext,edges[e2next],its[e1next],its[e2next],-1)
                    if rnext == -1:
                        print "this also failed so let us give up"
                        return -1,-1,-1,-1,[]
                    else:
                        return e1next,tnext,e2next,rnext,mp
                else:
                    return e1next,tnextn,e2next,rnext,mp


            else:
                return e1prev,tnextn,e2next,rnext,mp
        else:
            return e1next,tnext,e2prev,rnext,mp
    def iterate(self,e1,t,e2,r,mplast,its,cnts,edges,points,convex_vert,stepsize=0.05,counts=10):
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
            e1,t,e2,r,mp = self.updateT(e1,t,e2,r,its,cnts,edges,stepsize)
            print count,e1,t,e2,r,mp


            if e1 == -1:
                break



            points.append([mp])
            dist = np.linalg.norm(mp-self.p(edges[e1],t))
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
                pe,mindist,mink,dist = self.distCheck(e1,t,e2,r,mp,edges,its,cnts,convex_vert)
                noofsteps = 0

                #print pe
                try:
                    e1prev,tprev,e2prev,rprev = mplast[1:5]


                    if len(pe) >0 :
                        #mpcurr = [mp,e1,t,e2,r,dist]
                        print "dist failure with ",pe
                        #while noofsteps>1:
                        for popcount in xrange(3):
                            points.pop(-1)
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
                        print its[mpcurr[1]],its[mpcurr[3]]
                        mpfn,pes = self.retract(mpcurr,mplast,edges,its,pe,cnts,convex_vert)
                        print "retracted..first"

                        e1,t,e2,r = mpfn[1:5]
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



                        points.append([mpfn[0]])
                        for (e3,_,r) in pe:
                            if e3 in convex_vert or e3!=self.nextedge(mpfn[1],cnts):
                                self.iterate(mpfn[1],mpfn[2],e3,r,mpfn,its,cnts,edges,points,convex_vert,counts=counts)
                            if e3==17:
                                print mpfn[3]
                            self.iterate(e3,r,mpfn[3],mpfn[4],mpfn,its,cnts,edges,points,convex_vert,counts=counts)
                        if len(pe)>0:
                            break
                    else:
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
                                distv = self.distanceIt(e3[0],mpn,its,edges)
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
    def distCheck(self,e1,t,e2,r,mp,edges,its,cnts,convex_vert):

        #r is not used at all
        dist = np.linalg.norm(self.p(edges[e1],t)-mp)
        pe = []
        mindist3 = float('inf')
        mink = -1
        for k,_ in enumerate(edges):
            if k in self.friends(e1,cnts,convex_vert) or k in self.friends(e2,cnts,convex_vert):
            #or k==self.nextedge(e1,cnts) or k==self.prevedge(e2,cnts) or k==self.prevedge(e1,cnts) or k==self.prevedge(e2,cnts):
            # or k in tracededges:
                continue
            dist3,p,t = self.distanceIt(k,mp,its,edges)
            #print k,dist3,dist,t
            if dist3 <dist-1.0:
                pe.append([k,dist3,t])
                if dist3<mindist3:
                    mink = k
                    mindist3 = dist3
        return pe,mindist3,mink,dist
    def retract(self,mpcurr,mplast,edges,its,pe,cnts,convex_vert):
           e1 = mpcurr[1]
           tl = mpcurr[2]
           e2 = mpcurr[3]
           rl = mpcurr[4]
           dist= mpcurr[-1]
           tn = tl
           rn = rl
           # TODO: if mpcurr and mplast belongs to different edges
           if mpcurr[3]!=mplast[3]:#both e2s are not same
              e21 = mplast[3]
              e22 = mpcurr[3]
              e12 = mpcurr[1]
              e11 = mplast[1]
              _,t1,mp1,_ = self.footprintIt(edges[e22],its[e22].end()-0.001,edges[e12],its[e22],its[e12])
              if t1==-1:
                   _,t1,mp1,_ = self.footprintIt(edges[e21],its[e21].begin(),edges[e11],its[e22],its[e11])
                   if t1==-1:
                        _,t1,mp1,_ = self.footprintIt(edges[e21],its[e21].begin(),edges[e12],its[e22],its[e12])
              # distancecheck here
              pe1,dist30,mink,dist = self.distCheck(e12,t1,e22,its[e22].end()-0.001,mp1,edges,its,cnts,convex_vert)
              if len(pe1)>0:
                  #failed
                  #last e2 is not feasible with e1 so maybe return the

                  pes = {}
                  dist3 = float('inf')
                  _,t1,mp1,_ = self.footprintIt(edges[e21],its[e21].begin(),edges[e12],its[e21],its[e12])
                  e1 = e12
                  if t1==-1:
                      _,t1,mp1,_ = self.footprintIt(edges[e21],its[e21].begin(),edges[e11],its[e21],its[e11])
                      e1= e11
                  if t1==-1:
                      print "bug here at line#465"
                      return
                  dist = np.linalg.norm(mp1-self.p(edges[e1],t1))
                  mpfn = [mp1,e1,t1,e21,its[e21].beign(),dist]


                  for (e3,_,_) in pe:
                      dist1,pind,tind = self.distanceIt(e3,mpfn[0],its,edges)
                      #print dist1
                      pes = {}
                      if dist1<mpfn[-1]:
                          #print e3
                          pes[e3] = [dist1,pind,tind]
                      dist3 = min(dist3,dist1)
                  return mpfn,pes

                  '''
                  _,t1,mp1,_ = self.footprintIt(edges[e21],its[e21].begin(),edges[e11],its[e21],its[e11])
                  if t1!=-1:
                      dist = np.linalg.norm(mp1-self.p(edges[e11],t1))
                      mpcurr = [mp1,e11,t1,e21,0.0,dist]
                  else:
                      print "error here check footprint of ",e21,its[e21].begin(),e11
                      return
                  #dist3 = self.distance(edge[mink],mp1)
                  '''
              else:
                  mplast = [mp1,e12,t1,e22,1.0,dist]
           elif mpcurr[1] == mplast[1]:#both e1s same and both e2s same
               tf = mplast[2]
           else:
               #both e2s are same
                   #first try if beginning of current edge distance check is violated
                   #if not change tf to beginning of current edge
                   # if distance check is violated then make tl end of previous edge
                   _,r0,mp0,_ = self.footprintIt(edges[mpcurr[1]],0.0,edges[mpcurr[3]],its[mpcurr[1]],its[mpcurr[3]])
                   if r0 !=-1:
                       pe1,dist30,mink,dist = self.distCheck(mpcurr[1],0.0,mpcurr[3],r0,mp0,edges,its,cnts,convex_vert)
                       if len(pe1)>0:
                           '''
                           e1 = mplast[1]
                           e2 = mplast[3]
                           _,r1,mp1,_ = self.footprintIt(edges[e1],1.0,edges[e2],its[e1],its[e2])
                           if r1!=-1:
                               tl = 1.0
                               dist = np.linalg.norm(self.p(edges[e1],1.0)-mp1)
                               mpcurr = [mp1,e1,1.0,e2,r1,dist]
                               #dist3 = self.distance(edges[mink],mp1)
                           else:
                               print "error in retract",mplast,mpcurr
                          '''
                       else:
                           #dist3 = dist30
                           mplast = [mp0,e1,0.0,e2,r0,dist]
                           #tf = 0.0





                   else:
                       print "error in retracting",mpcurr,mplast
                       tf = its[e1].begin()
           count = -1
           pes ={} # to store footprints
           e1 = mpcurr[1]
           tl = mpcurr[2]
           e2 = mpcurr[3]
           rl = mpcurr[4]
           dist= mpcurr[-1]
           tn = tl
           rn = rl
           dist3 = float('Inf')
           tf = mplast[2]


           while abs(dist3-dist)>1.0 and abs(tf-tl)>0.001 and count < 20:
                count += 1
                tn = (tf + tl)/2.0
                _,rn,mpn,x = self.footprintIt(edges[e1],tn,edges[e2],its[e1],its[e2],-1)

                #print "rn =",rn
                if rn == -1: # it should not happen
                    print "possible bug in retract or footprint ",e1,e2,mplast,mpcurr
                    #return tn,rn,-1
                    break
                else:
                    dist = np.linalg.norm(mpn-self.p(edges[e1],tn))
                    dist3 = float('Inf')

                    for (e3,_,_) in pe:
                        dist1,pind,tind = self.distanceIt(e3,mpn,its,edges)
                        #print dist1
                        pes = {}
                        if dist1<dist:
                            #print e3
                            pes[e3] = [dist1,pind,tind]
                        dist3 = min(dist3,dist1)
                    print tf,tn,tl,rn,dist3,dist,mpn

                    if dist3 < dist:
                        tl = tn
                    else:
                        tf = tn

           if rn!=-1:
                mpfn = [mpn,e1,tn,e2,rn,dist]
           else:
               mpfn = mpcurr
           if dist3<dist and abs(dist3-dist)>1.0:
               print "possible bug in retraction, new dist difference is ",abs(dist3-dist),mpcurr,mplast
               #return tn,rn,-1
           return mpfn,pes
    def modEdges(self,letter):
        edges,vertices,vert_edges,cnts = check.fontOutline(letter)
        convex_vert,concave_vert = self.vertAnalyse(vertices,edges,cnts)
        its = {} # store as interval tree
        for k in xrange(len(edges)+len(concave_vert)):
            its[k]=iv(0.0,1.001)
        its = self.insertReflex(edges,vertices,cnts,concave_vert,its)
        for k in xrange(len(edges)):
            its[k] = it([its[k]])
        convex_vert,concave_vert = self.vertAnalyse(vertices,edges,cnts)
        return edges,its,cnts,vertices,convex_vert,concave_vert
    def remove_interval(self,its,e1,t0,tn):
        its[e1].chop(t0,tn)

    def footprint(self,quad1,t,quad2,initialr=-1,diagnostics=False):
        edges = []
        p1 = self.p(quad1,t)
        p1prime = self.pprime(quad1,t)
        normal1 = np.array([p1prime[1],-1*p1prime[0]])

        edges.append([p1])
        edges.append([p1,p1+p1prime*50/np.linalg.norm(p1prime)])
        #print np.linalg.norm(p1prime)


        if initialr <0:
            _,_,r = self.distance(quad2,p1)
        else:
            r = initialr
        #print r
        if r==-1:
            r = 0.5
        result = False
        for count in xrange(100):
            if r != -1 and r>=0.0 and r<=1.0:
                p2 = self.p(quad2,r)
                p2prime = self.pprime(quad2,r)
                normal2 = np.array([p2prime[1],-1*p2prime[0]])
                #edges.append([p2])
                edges.append([p2,p2+normal2*50/np.linalg.norm(normal2)])

                denominator = np.dot(normal1,p2prime)
                if abs(denominator) > 0.01 :
                    alpha = np.dot(p2-p1,p2prime)/denominator
                    beta = np.dot(p2-p1,p1prime)/denominator
                    #print alpha>0,beta>0



                else: #parallel found
                    alpha = beta = np.linalg.norm(p1-p2)/2.0
                    alpha = alpha/np.linalg.norm(p1prime)
                    beta = beta/np.linalg.norm(p2prime)




                r,result = self.solve(quad1,t,quad2,r)
            #    print r
                if result:
                    break
        if result and alpha>=0 and beta>=0 :
            #edges.append([p2,p2+normal2*50/np.linalg.norm(normal2)])
            edges.append([p1,p1 + alpha*normal1])
            edges.append([p2,p2 + beta*normal2])



            if diagnostics:
                return edges,r,p1 + alpha*normal1,t,r
            else:
                return edges,r,p1+alpha*normal1
        else:
            if diagnostics:
                return edges,-1,None,t,r
            else:
                return edges,-1,None
    def footprintIt(self,quad1,t,quad2,it1,it2,initialr=-1,diagnostics=False):#using interval tree
        if( it1.search(t) == set()):
            return [],-1,[],None
        if diagnostics:
            e,r,mp,tfin,rfin = self.footprint(quad1,t,quad2,initialr,diagnostics)
        else:
            e,r,mp = self.footprint(quad1,t,quad2,initialr,False)



        #print "line 263:r is",r
        x = it2.search(r)
        #y = it2.search(r)

        if x!=set():
            assert len(x) == 1 # if not then we have not managed intervals better
            if diagnostics:
                return e,r,mp,x,tfin,rfin
            else:
                return e,r,mp,x
        else:
            print x,r,"in footprint it",it2

            if diagnostics:
                return e,-1,[],None,tfin,rfin
            else:
                return e,-1,[],None


    def bezIntersect(self,bez1,bez2):
        return
    def bisector(self,quad1,quad2):
        finedges = []
        points = []
        for q in xrange(101):
            t = q/100.0
            edges,r,point = self.footprint(quad1,t,quad2)
            if r!=-1:
                print t,r

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
    def nextedge(self,e,cnts): #gives the next edge given an index of edge
       for cnt in cnts:
           for i,pt in enumerate(cnt):
               if pt == e:
                   if i <len(cnt)-1:
                       return cnt[i + 1]
                   else:
                       return cnt[0]
    def prevedge(self,e,cnts): #gives the next edge given an index of edge
       for cnt in cnts:
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
    def medialaxis(self,edges1,vertices1,vert_edges1,cnts1,stepsize=0.05):

        def distance(ind,p0,its,edges):#distance between a point and p
           quad = edges[ind]
           pdprime = 2*(quad[2]+quad[0]-2*quad[1])
           t=t0 = its[ind].begin()
           tf = its[ind].end()
           result = False
           for j in xrange(100): # maximum 100 iterations
                #print t
                val = np.dot(self.pprime(quad,t),self.p(quad,t)-p0)
                if val>=0 and abs(val) <= 0.01:
                    result = True
                    break
                else:
                    fprime = np.linalg.norm(self.pprime(quad,t))**2 + np.dot(pdprime,self.p(quad,t)-p0)
                    t = t - val/fprime
                    if t<t0 or t>tf:
                        t = random.uniform(t0,tf)
           if result:
               p1 = self.p(quad,t)
               p1prime = self.pprime(quad,t)
               if np.cross(p0-p1,p1prime) >=0: # point is on right side of the side
                    return np.linalg.norm(p1-p0),p1,t


           return float("inf"),[],-1
        def reflexVertexOld(i):
            quad = []
            linked_edges = vert_edges[i]
            pt = vertices[i]

            if len(linked_edges) ==2:
                linked_edges.sort()
                quad.append(self.p(edges[linked_edges[0]],0.99))
                quad.append(pt)

                quad.append(self.p(edges[linked_edges[1]],0.01))
                return quad
            else:
                quad.append(self.p(edges[-1],0.99))
                quad.append(pt)

                quad.append(self.p(edges[linked_edges[0]],0.01))
            return quad
        def reflexVertex(i): # store as normals
            e1 = vert_edges[i][0]
            if len(vert_edges[i]) ==1:
                e2 = len(edges) -1
            else:
                e1 = vert_edges[i][1]
            p1prime = self.pprime(edges[e1],1.0)
            n1 = np.array([p1prime[1],-1*p1prime[0]])
            p2prime = self.pprime(edges[e2],1.0)
            n2 = np.array([p2prime[1],-1*p2prime[0]])
            return np.array([n1,n2])
        def remove_interval(its,e1,t0,tn):
            its[e1].slice(t0)
            its[e1].slice(tn)
            its[e1].remove_overlap(t0,tn)
        def retract(simpleedge,mpcurrent,mpdistl,e3,dist3,edges,its):
               e1 = mpdistl[1]
               e2 = mpdistl[3]
               rf = mpdistl[4]
               tf = mpdistl[2]
               rl = mpcurrent[4]
               tl = mpcurrent[2]
               assert e1 == mpcurrent[1]
               assert e2 == mpcurrent[3]
               dist = mpcurrent[-1] # current dist which needs to be reduced
               #print mpdistl
               #print mpcurrent
               #print "before retract distance from edges e1,e2,e3",dist,dist3,tf,tl
               #print simpleedge
               simpleedge.pop(-1)
               count = 0
               rn = rl
               tn = tl
               mpn = mpcurrent[0]


               while abs(dist3-dist)>1.0 and count < 10:
                    count += 1
                    tn = (tf + tl)/2.0
                    #print tf,tn,tl,rl,dist3,dist
                    _,rn,mpn,x = self.footprintIt(edges[e1],tn,edges[e2],its[e1],its[e2])
                    #print "rn =",rn
                    if rn == -1: # it should not happen
                        print "possible bug in retract or footprint ",e1,e2
                        return tn,rn,-1
                        break
                    else:
                        dist3,pind,tind = distance(e3,mpn,its,edges)
                        #dist3o,_,_ = self.distance(edges[e3],mpn)
                        #print dist3o,dist3
                        #assert abs(dist3o-dist3) < 1.0

                        dist = np.linalg.norm(mpn-self.p(edges[e1],tn))
                        #diste2 = np.linalg.norm(mpn-self.p(edges[e2],rn))
                        #print "diste1,diste2",diste1,diste2
                        #dist = max(diste1,diste2)

                        if dist3 < dist:
                            tl = tn
                        else:
                            tf = tn

               if rn!=-1:
                    mpfn = np.array([mpn,e1,tn,e2,rn,dist])
                    simpleedge.append(mpfn)
                    #t0 = simpleedge[0][2]
                    #r2 = simpleedge[0][4]
                    #remove_interval(its,e1,t0,tn)
                    #remove_interval(its,e2,rn,r2)
                    #print "after retract distance from edges e1,e2,e3",dist,dist3
                    #print simpleedge
               if count ==10 and abs(dist3-dist)>1.0:
                   print "possible bug in retraction, new dist difference is ",abs(dist3-dist)
                   return tn,rn,-1
               else:
                   return tn,rn,0
        def hp(e1,e2,simpleedge,maxis,mpl,its):
            if e1 == e2:
                return
            #get overlap intervals
            # start incrementing e1
            # after three increments,
            #distance check
            # if distance fail, retract and proceed to new pairs
            global count
            print "handling ",e1,e2
            print its[e1],its[e2]
            breakpoint = False
            success = False
            #print mpl
            for interval_obj in its[e1]:
                if interval_obj.end-interval_obj.begin <=.01:
                    continue
                print "interval is ",interval_obj
                t1 = interval_obj.begin
                #print t1
                _,r2,mp1,x=self.footprintIt(edges[e1],t1,edges[e2],its[e1],its[e2],initialr=-1)
                #print "r2 is ",r2
                success = False
                if r2==-1: # fail, now we shall try e2
                    # let us try t2 as end of  interval
                     t2 = interval_obj.end - 0.0001
                     _,r1,mp2,x1 =self.footprintIt(edges[e1],t2,edges[e1],its[e2],its[e1],initialr=-1)
                     if r1 == -1:
                         # we have to rule out if any e2 interval is totally inside the interval or not
                         print "may overlap "
                     else: #found the interval of e2
                         r2 = [y.end for y in x1]
                         r2 = r2[0]
                         success = True

                else: # we have found t1
                   #let us try t2 as last of current e1 interval
                   success = True
                   t2 = interval_obj.end - 0.0001
                 #  print "line 523:t2 is ",t2
                   _,r1,mp2,x1 =self.footprintIt(edges[e1],t2,edges[e2],its[e1],its[e2],initialr=-1)
                  # print "line 525: r1 is",r1
                   if r1 == -1:
                       # e2 interval is smaller
                       # r1 = x begin
                       r1 = [y.begin for y in x]
                       r1 = r1[0]
                    #   print "line 535:r1 is",r1
                       #we need to find t2
                       _,t2,mp2,x2 =self.footprintIt(edges[e2],r1,edges[e1],its[e2],its[e1],initialr=-1)





                if success:
                    break
            if not success:
                print "checking from e2 side"
                for interval_obj in its[e2]:
                    if interval_obj.end-interval_obj.begin<0.01:
                        continue
                    r2 = interval_obj.end-0.001
                    _,t1,mp1,x=self.footprintIt(edges[e2],r2,edges[e1],its[e2],its[e1],initialr=-1)
                    if t1!=-1:
                        r1 = interval_obj.begin
                        _,t2,mp2,x=self.footprintIt(edges[e2],r1,edges[e1],its[e2],its[e1],initialr=-1)
                        success = True
                        break




                return #nothing to be done
            if not success:
                #hp(e1,self.prevedge(e2,cnts),simpleedge,maxis,[],its)

                return
            else:
                # we have found first footprint
                print "t1,t2,r1,r2",t1,t2,r1,r2
                mpcurrent = np.array([mp1,e1,t1,e2,r2,np.linalg.norm(self.p(edges[e1],t1))])
                simpleedge.append(mpcurrent)
                print count,e1,e2
                count +=1
                distFailure = False
                t = t1
                noofsteps = 0
                dist = mpcurrent[-1]
                mpdistl = mpcurrent
                nextpairs = []
                while t<t2 and not distFailure:
                    t = min(t2,t + deltat)
                    _,r,mpcurrent1,_=self.footprintIt(edges[e1],t,edges[e2],its[e1],its[e2],initialr=-1)
                    noofsteps += 1
                    #count +=1
                    if r==-1:
                        print "let's investigate",e1,t,e2
                        break
                    mpcurrent = np.array([mpcurrent1,e1,t,e2,r,np.linalg.norm(mpcurrent1-self.p(edges[e1],t))])
                    simpleedge.append(mpcurrent)

                    #print mpcurrent
                    dist = mpcurrent[-1]
                    #print dist,noofsteps,t
                    if True or noofsteps >= maxnoofsteps:
                        possibleedges = []
                        noofsteps=0
                        for ind,_ in enumerate(edges):
                            if ind == e1 or ind == e2:
                                continue

                            dist1,pind,tind = distance(ind,mpcurrent[0],its,edges)
                            #print dist,dist1,e1,e2,ind

                            if dist1 <= dist: # found a third edge
                                 #TODO retract mp and delete additional points
                                 #retract(simpleedge,mp,)
                                 #print "line 584:",dist1,dist,e1,ind,mpcurrent,mpdistl
                                 #return
                                 print "distance failure of e1,e2,ind",e1,e2,ind,dist1,dist
                                 tn,rn,retractval=retract(simpleedge,mpcurrent,mpdistl,ind,dist1,edges,its)
                                 if retractval==-1:
                                     print "line 588:",dist1,dist,e1,ind,mpcurrent,mpdistl
                                     #return



                                 #return
                                 possibleedges.append(ind)
                                 #print dist,dist1,ind,tind
                                 break
                        if len(possibleedges)>0:

                            distFailure = True
                            maxis.append(simpleedge)
                            simpleedge = []
                            for ind in possibleedges:
                                nextpairs.append([ind,e1,copy.copy(mpl)])
                                nextpairs.append([ind,e2,copy.copy(mpl)])
                        mpdistl = copy.copy(mpcurrent)
                #print t1,t2,r1,r2
                if not distFailure:
                    if len(simpleedge)>0:
                        maxis.append(simpleedge)
                    simpleedge=[]
                    remove_interval(its,e1,t1,t2)
                    remove_interval(its,e2,r1,r2)
                    if its[e1].is_empty():
                        if not its[e2].is_empty():
                            nextpairs.append([self.nextedge(e1,cnts),e2,mpcurrent[0]])
                        else:
                            nextpairs.append([self.nextedge(e1,cnts),self.prevedge(e2,cnts),mpcurrent[0]])
                    else:
                        nextpairs.append([e1,self.prevedge(e2,cnts),mpcurrent[0]])



                    #nextpairs.append([e1,self.prevedge(e2,cnts),mpcurrent[0]])
                    #nextpairs.append([self.nextedge(e1,cnts),self.prevedge(e2,cnts),mpcurrent[0]])
                    #

                else:
                        print tn,rn
                        remove_interval(its,e1,t1,tn)
                        remove_interval(its,e2,rn,r2)
                print its[e1],its[e2]




                for pair in nextpairs:
                    hp(pair[0],pair[1],simpleedge,maxis,pair[2],its)
                    #print "we shall do it later",pair[0],pair[1]
                return

        maxnoofsteps = 3
        cnts = copy.copy(cnts1)
        edges = copy.copy(edges1)
        vertices = copy.copy(vertices1)
        vert_edges = copy.copy(vert_edges1)
        convex_vert,concave_vert = self.vertAnalyse(vertices,edges,cnts)
        len(edges)
        print convex_vert,concave_vert
        its = {} # store as interval tree
        for k in xrange(len(edges)+len(concave_vert)):
            its[k]=iv(0.0,1.0)


        its = self.insertReflex(edges,vertices,cnts,concave_vert,its)
        for k in xrange(len(edges)):
            its[k] = it([its[k]])
        convex_vert,concave_vert = self.vertAnalyse(vertices,edges,cnts)
        print convex_vert,concave_vert,its
        len(edges)
        #return





        maxis = []
        simpleedge =[]
        count = 0
        edge_count = len(edges) -1
        deltat = stepsize
        alledges = []
        tracededges = []
        print cnts
        #handlePair(0,17,simpleedge,maxis,vertices[0])

        for cvert in convex_vert:
            firstvertex = cvert

            secondedge = self.prevedge(firstvertex,cnts)
            firstedge = firstvertex
            #lastts[firstedge][0] = 0.0

            mp = [vertices[firstvertex],firstedge,0.0,secondedge,1.0]

            #print firstedge,secondedge
            nextpairs = hp(firstedge,secondedge,simpleedge,maxis,mp,its)
            #break

            #print nextpairs


        return maxis,edges
    def printMa(self,ma):
        points = []
        for k,ed in enumerate(ma):
            for point in ed:
                if point != None:
                    if len(point)>0:
                        points.append([point[0]])
        return points
    def fillPointsAtCorner(self,edge1,edge2,dist):
        slope1 = self.pprime(edge1,1)
        slope2 = self.pprime(edge2,0)
        normal1 = np.array([slope1[1],-1*slope1[0]])
        normal2 = np.array([slope2[1],-1*slope2[0]])
        pts = []
        for q in xrange(11):
            t = q/10.0
            normal = normal1*(1-t) + normal2 * t
            pts.append(edge1[-1] +normal*dist/np.linalg.norm(normal1) )
        return pts
    def vertAnalyse(self,vertices,edges,cnts):
        concave_vert = []
        convex_vert = []
        print len(vertices),len(edges)
        for i,_ in enumerate(vertices):
            j = self.prevedge(i,cnts)
            cross = np.cross(self.pprime(edges[i],0.0),self.pprime(edges[j],1.0))
            #print i,j,cross

            if cross>0.01:
                #print "convex"
                convex_vert.append(i)
            elif cross<-0.01:
                #print "concave"
                concave_vert.append(i)
        return convex_vert,concave_vert


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
    def test(self,letter,stepsize=0.05):
        edges,vertices,vert_edges,cnts = check.fontOutline(letter)
        maa,mod = self.medialaxis(edges,vertices,vert_edges,cnts,stepsize)
        eda = self.printMa(maa)
        check.printEdges(mod)
        check.printEdges(mod+eda,0,False)
        return maa,mod,eda
    def test1(self,quad1,quad2):
        edges =[]
        for q in xrange(11):
            t = q/10.0
            p1 = self.p(quad1,t)
            p2 = self.p(quad2,t)
            p1prime = self.pprime(quad1,t)
            p2prime = self.pprime(quad2,t)
            n1 = np.array([p1prime[1],-1*p1prime[0]])
            n2 = np.array([p2prime[1],-1*p2prime[0]])
            n = n2 - np.dot(p2,n2)/np.dot(p1,n1)*n1
            edges.append([p1+n])
            #edges.append([p2+n])
        c= pg3.Check()
        c.printEdges(edges+[quad1,quad2],0,False)
