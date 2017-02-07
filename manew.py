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
import debug
import scipy
import bezier
def sgolay2d ( z, window_size, order, derivative=None):
    """
    """
    # number of terms in the polynomial expression
    n_terms = ( order + 1 ) * ( order + 2)  / 2.0

    if  window_size % 2 == 0:
        raise ValueError('window_size must be odd')

    if window_size**2 < n_terms:
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2

    # exponents of the polynomial.
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ...
    # this line gives a list of two item tuple. Each tuple contains
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [ (k-n, n) for k in range(order+1) for n in range(k+1) ]

    # coordinates of points
    ind = np.arange(-half_size, half_size+1, dtype=np.float64)
    dx = np.repeat( ind, window_size )
    dy = np.tile( ind, [window_size, 1]).reshape(window_size**2, )

    # build matrix of system of equation
    A = np.empty( (window_size**2, len(exps)) )
    for i, exp in enumerate( exps ):
        A[:,i] = (dx**exp[0]) * (dy**exp[1])

    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2*half_size, z.shape[1] + 2*half_size
    Z = np.zeros( (new_shape) )
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] =  band -  np.abs( np.flipud( z[1:half_size+1, :] ) - band )
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band  + np.abs( np.flipud( z[-half_size-1:-1, :] )  -band )
    # left band
    band = np.tile( z[:,0].reshape(-1,1), [1,half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs( np.fliplr( z[:, 1:half_size+1] ) - band )
    # right band
    band = np.tile( z[:,-1].reshape(-1,1), [1,half_size] )
    Z[half_size:-half_size, -half_size:] =  band + np.abs( np.fliplr( z[:, -half_size-1:-1] ) - band )
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z

    # top left corner
    band = z[0,0]
    Z[:half_size,:half_size] = band - np.abs( np.flipud(np.fliplr(z[1:half_size+1,1:half_size+1]) ) - band )
    # bottom right corner
    band = z[-1,-1]
    Z[-half_size:,-half_size:] = band + np.abs( np.flipud(np.fliplr(z[-half_size-1:-1,-half_size-1:-1]) ) - band )

    # top right corner
    band = Z[half_size,-half_size:]
    Z[:half_size,-half_size:] = band - np.abs( np.flipud(Z[half_size+1:2*half_size+1,-half_size:]) - band )
    # bottom left corner
    band = Z[-half_size:,half_size].reshape(-1,1)
    Z[-half_size:,:half_size] = band - np.abs( np.fliplr(Z[-half_size:, half_size+1:2*half_size+1]) - band )

    # solve system and convolve
    if derivative == None:
        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, m, mode='valid')
    elif derivative == 'col':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -c, mode='valid')
    elif derivative == 'row':
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid')
    elif derivative == 'both':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid'), scipy.signal.fftconvolve(Z, -c, mode='valid')
def stitch(bezier1,bezier2,ini1=1,ini2=0,distthreshold=10.0,deltat=0.01):
    bz1 = bezier.Curve(np.array(bezier1),3);bz2=bezier.Curve(np.array(bezier2),3)
    currdist = np.linalg.norm(bezier1[0-ini1]-bezier2[0-ini2])
    if currdist>distthreshold:
        return bezier1,bezier2
    lastt = 1.0;firstt = 0.0;count =0
    dist = currdist - 0.001
    while dist>0.0001 and dist< currdist and count<10:
        currdist = dist
        lastt += deltat;firstt -= deltat
        if ini1==1:
            bz1n = bz1.specialize(0.0,lastt)
        else:
            bz1n = bz1.specialize(firstt,1.0)
        if ini2==1:
            bz2n = bz2.specialize(0.0,lastt)
        else:
            bz2n = bz2.specialize(firstt,1.0)
        dist = np.linalg.norm(bz1n.nodes[0-ini1]-bz2n.nodes[0-ini2])
        count += 1
        print dist
    return bz1n.nodes,bz2n.nodes



check = pg3.Check()
count = 0
pairs = {}
tracededges=[]
class Node(object):
    def __init__(self, value, left=None, right=None):
        self.value = value
        self.left = left
        self.right = right
class fp(object):
    maxcount = 1
    def __init__(self,letter,fontfile="DevanagariSangamMN.ttf"):
        #self.edges,self.its,self.cnts,self.vertices,self.convex_vert,self.concave_vert=
        self.debug = True;self.debug2=False
        self.reflex_edges=[]
        self.stroke_branches={}
        self.stroke_count=0
        self.maxis=[]

        self.modEdges(letter,fontfile)
        self.untracededges=[]
        self.tracededges = []
        self.reflex_objects = {}
        for k,_ in enumerate(self.edges):
            self.untracededges.append(k)

        for k in self.reflex_edges:
            r = Reflex()
            r.update(self,k)
            self.reflex_objects[k]=r
        #for k in self.reflex_edges:

            #r = self.reflex_objects[k]
            #r.updateRad(self,k)

    def secant(self,f,rs):
        fxi = f(rs[-1])
        fximinus1 = f(rs[-2])
        #print fxi,fximinus1
        rn = rs[-1]
        #return
        for i in xrange(50):
            if (abs(rs[-1]-rs[-2])<0.000000001 or abs(fxi-fximinus1)<0.0000000001 or abs(fxi)<0.000000001):

                break
            ri = rs[-1]
            rn = ri -fxi *(rs[-1]-rs[-2])/(fxi-fximinus1)
            #print rn
            rs.append(rn)
            fximinus1 = fxi
            fxi = f(rn)
        return rn
    def ccw(self,A,B,C):
            return (C[1]-A[1])*(B[0]-A[0]) > (B[1]-A[1])*(C[0]-A[0])

    def p(self,quad,t):
        if len(quad) ==3:
            return quad[0]*(1-t)**2 + quad[1]*2*t*(1-t) + quad[2]*t**2
        else:
            return quad[0]
    def pprime(self,quad,t):
        if len(quad) == 1:
            return []
        return 2*((quad[1]-quad[0])*(1-t) + (quad[2]-quad[1])*t)
    def normal(self,quad,t,unit=False):
        if len(quad)==1:
            return []
        pp = self.pprime(quad,t)
        n0 = np.array([pp[1],-1*pp[0]])
        if unit:
            return n0/np.linalg.norm(pp)
        else:
            return n0
    def footprint(self,e1,t,e2,initialr=-1):
        if e1 in self.reflex_edges:
            e1r = self.reflex_objects[e1];return e1r.footprintn(e2,t,-1)
        if e2 in self.reflex_edges:
            e2r = self.reflex_objects[e2];return e2r.footprintn(e1,-1,t)


        quad1 = self.edges[e1];quad2 = self.edges[e2]
        p1 = self.p(quad1,t)
        p1prime = self.pprime(quad1,t)
        n1 = np.array([p1prime[1],-1*p1prime[0]])
        normn1 = np.linalg.norm(n1)
        dist0,p2i,r1 =self.distance(e2,p1) #randomly ..no partcular thought
        edges = []
        edges.append([p1,p2i])
        #print dist0,r1
        if r1 <0.0:
            #print "here"
            p2i = self.p(quad2,0.0);dist0 = np.linalg.norm(p2i-p1)
        elif r1>1.0:
            p2i = self.p(quad2,1.0);dist0 = np.linalg.norm(p2i-p1)


        for i in xrange(20):
            costheta1 = np.dot(p2i-p1,n1)/(np.linalg.norm(n1)*dist0)
            if abs(costheta1)>0.01:
                dist1 =  dist0/(2*costheta1)
            else:
                dist1 = dist0/2.0
            fpi = p1 +dist1*n1/normn1
            dist2,p2i,r1 =self.distance(e2,fpi) #randomly ..no partcular thought
            '''
            if r1>1:
                fpi = p1 + dist1/2.0*n1/normn1
                dist2,p2i,r1 =self.distance(e2,fpi) #randomly ..no partcular thought

            elif r1<0:
                fpi = p1 + dist1*1.5*n1/normn1
                dist2,p2i,r1 =self.distance(e2,fpi) #randomly ..no partcular thought
            '''

            #print dist2,p2i,r1
            #if r1<0.0 or r1>1.0:
            #    break
            if abs(dist2-dist1)<0.00001:
                break
            else:
                dist0 = np.linalg.norm(p2i-p1)
        edges.append([p1,fpi]);edges.append([p2i,fpi]);edges.append([p1]);edges+=[[p2i],[p1,p2i],[fpi]]

        #print r1,np.cross(fpi-p1,p1prime),abs(dist2-dist1)
        if np.cross(fpi-p1,p1prime)>=0 and abs(dist2-dist1)<1.0:
          return edges,r1,fpi
        else:
          return edges,-1,fpi
    def distance(self,e1,p0,debug=False):
        if e1 in self.reflex_edges:
            e1r = self.reflex_objects[e1]
            return e1r.distance(p0)
        quad = self.edges[e1]
        def f(t):
                #if t>=1:
                #    return np.dot(self.pprime(quad,1.0),self.p(quad,1.0)-p0)
                #elif t<=0:
                #    return np.dot(self.pprime(quad,0.0),self.p(quad,0.0)-p0)
                #else:
                    return np.dot(self.pprime(quad,t),self.p(quad,t)-p0)
        rs = [];rs.append(0.4);rs.append(0.6)
        rn = self.secant(f,rs)
        p1 = self.p(quad,rn);p1prime=self.pprime(quad,rn)
        if debug:
            ep=[];ep.append([p0,p1]);ep.append([p1,p1+p1prime]);check.printEdges(ep+self.edges,0,False)
        #print abs(f(rn))
        if np.cross(p0-p1,p1prime) >=0 and abs(f(rn))<1.0: # point is on right side of the side
             return np.linalg.norm(p1-p0),p1,rn
        elif np.cross(p0-p1,p1prime) <0 and abs(f(rn))<1.0:
            return -1*np.linalg.norm(p1-p0),p1,rn
        else:
            return float('inf'),p1,rn
    def distanceIt(self,ind,p0,stepsize=0.01):#distance between a point and p
           #print "distanceIt,ind,p0",ind,p0
          # if ind==37:
            #   print self.its[37]
           if self.its[ind].end()-self.its[ind].begin()<0.001+stepsize:
               return float("inf"),[],-1
           if ind in self.reflex_edges:
               indr = self.reflex_objects[ind]
               return indr.distance(p0)
           dist,p1,t = self.distance(ind,p0)
           if len(self.its[ind].search(t-0.0001))!=0 or len(self.its[ind].search(t+0.0001))!=0:
               return dist,p1,t
           else:
               return float("inf"),[],-1
    def manage_interval(self,mplast,mpcurr):
        assert mplast[1]==mpcurr[1] and mplast[3]==mpcurr[3]
        if self.debug: print "now managing interval of",mplast,mpcurr
        self.its[mplast[1]].chop(mplast[2],mpcurr[2])
        self.its[mplast[3]].chop(mpcurr[4],mplast[4])

    def nextt(self,e1,t,e2,stepsize=0.05):
        #print "now in nextt"
        iv1 = sorted(self.its[e1].search(t))
        if len(iv1) == 0:
            if self.debug: print "not found t",t," in",self.its[e1]
            return -1,-1

        elif t+stepsize <= iv1[-1].end-0.001:
            return e1,t+stepsize
        #elif t<iv1[-1].end-0.001:
        #    return e1,iv1[-1].end-0.001
        else:
            #iterate with all intervals
            q1 = sorted(self.its[e1])
            for k,iv0 in enumerate(q1):
                if t>=iv0.begin and t<=iv0.end:
                    j=k
                    while j!=len(q1)-1:
                        if self.debug: print "this happened in nextt"
                        if q1[j+1].end-q1[j+1].begin>0.001:
                            return e1,q1[j+1].begin
                        j = j+1
            preve1 = e1
            e1 = self.nextedge(e1)
            if e1!=e2:# and not e1 in self.concave_vert:
                #self.remove_interval(its,preve1,t,its[preve1].end())
                if self.debug: print "from nextt change of edge",e1,self.its[e1].begin()
                return e1,self.its[e1].begin()
            #elif e1 in self.concave_vert:
            #    return -e1,0.0
            else:
                return -1,-1
    def prevr(self,e2,r,e1,stepsize=0.05):
        iv1 = sorted(self.its[e2].search(r))
        if len(iv1) == 0:
            return -1,-1

        elif r-stepsize >= iv1[0].begin:
            return e2,r-stepsize
        else:
            #iterate with all intervals
            q1 = sorted(self.its[e2])
            for k in xrange(len(q1)-1,-1,1):
                if r>=q1[k].begin and r<=q1[k].end-0.001:
                    if k!=0:
                        return e2,q1[k-1].end-0.001
            preve2 = e2
            e2 = self.prevedge(e2)
            if e2!=e1:# and not e2 in self.concave_vert:
                #self.remove_interval(its,preve1,t,its[preve1].end())
                return e2,self.its[e2].end()-0.001
            #elif e2 in self.concave_vert:
            #    return -e2,0.0
            else:
                return -1,-1
    def updateTn(self,e1prev,tprev,e2prev,rprev,stepsize=0.05):

        e1next,tnext = self.nextt(e1prev,tprev,e2prev,stepsize)
        #print "in update ",e1prev,tprev,e1next,tnext
        if tnext == -1 and e1next==-1:
            return -1,-1,-1,-1,[]
        _,rnext,mp,x = self.footprintIt(e1next,tnext,e2prev)
        if x!=None:
            return e1next,tnext,e2prev,rnext,mp
        else:
            if self.debug: print "this has failed..",e1next,tnext,e2prev,rnext,mp
            if rnext>self.its[e2prev].end():
                _,t1next,mp,x = self.footprintIt(e2prev,self.its[e2prev].end(),e1next)
                if x!=None:
                    return e1next,tnext,e2prev,self.its[e2prev].end(),mp
            e2next = self.prevedge(e2prev)
            if self.debug: print "so now trying",e1next,tnext,e2next
            _,rnext,mp,x = self.footprintIt(e1next,tnext,e2next)
            if rnext>=0.0 and rnext<=1.0 and x!=None:
                return e1next,tnext,e2next,rnext,mp
            else:
                return -1,-1,-1,-1,[]
    def negdistcheck(self,p0):
        negdist = True
        for cnt in self.cnts:
            negdist = False
            for i in cnt:
                dist,_,_ = self.distanceIt(i,p0)
                #print "in negdist",dist
                negdist = negdist or dist>0 or dist==float('inf')
            if not negdist:
                print "found neg distance"
                break
        return not negdist

    def handleReflex(self,e,pts,stepsize=0.001):
        #pts = []
        er = self.reflex_objects[e]
        e1 =er.rad[-4];r1 = er.rad[-2];e2=er.rad[-3];r2=er.rad[-1]
        e1next = e1;rnext = r1;rad1 = er.rad[0];rad2=er.rad[1]
        print e1,e2,r1,r2
        t=1.0;noofstep = 0;finishededges = [e]
        while t>=0:
            print "in handlereflex",e1next,rnext
            if not e1next in finishededges:
                finishededges.append(e1next)
            q=er.footprintn(e1next,-1,rnext,rad2)
            t = q[1]
            mp = q[-1];
            noofstep += 1
            fpexists = False
            if not fpexists and noofstep>2:# and not e1next in self.reflex_edges:
                noofstep=0
                for enew,_ in enumerate(self.edges):
                    if enew ==e1next or enew in finishededges: # or enew in self.reflex_edges
                        continue
                    fpexists,t0=self.checkFpPoint(e1next,enew,mp)
                    if fpexists==True:
                        break
            if fpexists:
                if self.debug: print "dist failure with enew",enew
                finishededges.append(enew)

                e1next= enew;rnext= t0;continue



            e1next,rnext = self.nextt(e1next,rnext,e2,stepsize)

            print t,e1next,rnext
            if e1next==-1 and rnext==-1:
                break

            pts.append([mp])
        tfin=t;e2prev = e2;rprev = r2;t=0.0;noofstep=0;#finishededges = []
        while t<=tfin:
            print "in handlereflex",e2prev,rprev
            if not e2prev in finishededges:
                finishededges.append(e2prev)


            q=er.footprintn(e2prev,-1,rprev,rad1)
            t = q[1]
            mp = q[-1];

            fpexists = False;noofstep += 1
            '''
            for e0 in self.reflex_edges:
                if e0==e:
                    continue

                fpexists,t0=self.checkFpPoint(e2prev,e0,mp)
                if fpexists==True:
                    break

            if fpexists:
                print "dist failure with e0",e0
                e2prev = e0;rprev=t0;continue
            '''
            if not fpexists and noofstep>2 and not e2prev in self.reflex_edges :
                for enew,_ in enumerate(self.edges):
                    if enew ==e2prev  or enew in finishededges:#or enew in self.reflex_edges
                        continue
                    fpexists,t0=self.checkFpPoint(e2prev,enew,mp)
                    if fpexists==True:
                        break
            if fpexists:
                print "dist failure with enew",enew
                finishededges.append(enew)

                e2prev= enew;rprev= t0;continue
            pts.append([mp])

            e2prev,rprev = self.prevr(e2prev,rprev,e1,stepsize)
            print t,e2prev,rprev

            if e2prev==-1 and rprev==-1:
                break

        return pts
    def checkFpPoint(self,e1,e2,mp):
        fpexists=False;t0=-1
        #print "checking",e1,e2,mp
        dist,_,_ = self.distanceIt(e1,mp)
        dist3,_,t0 = self.distanceIt(e2,mp)
        if dist3>0 and abs(dist3)<=dist:
            fpexists=True
        return fpexists,t0

        if e2 in self.reflex_edges:
            e2r = self.reflex_objects[e2]
            dist3 = min(e2r.rad[0],e2r.rad[1])
            if dist3<=dist:
                vect = mp-self.vertices[e2]
                theta1 =np.arctan2(vect[1],vect[0])
                if theta1<0:
                    theta1 += 2*np.pi + theta1
                t0 = (theta1-e2r.theta1)/(e2r.theta2-e2r.theta1)
                if t0>=0.0 and t0<=1.0:
                    fpexists=True
        else:
            dist3,_,t0 = self.distanceIt(e2,mp)
            if abs(dist3)<=dist:
                fpexists=True
        return fpexists,t0
    def iterate(self,e1,t,e2,r,mplast,points,stepsize=0.05,counts=10,MinDistBeforeDistCheck=20.0):
        if e1==e2:
            return
        if self.debug: print "beginning",e1,t,e2,r,counts
        e1in,tin,e2in,rin = e1,t,e2,r
        count = -1
        noofsteps = 0
        countlast = -1
        #points.append(mplast)
        #maxis.append(points)
        savedmplast = mplast
        #mplastn = self.footprintIt(edges[e1],t,edges[e2],its[e1],its[e2],r)
        #mplast = [mplastn,e1,t,e2,r]
        while count < counts:
            count += 1
            noofsteps += 1
            e1prev,tprev,e2prev,rprev = e1,t,e2,r
            edgechange = False
            e1,t,e2,r,mp = self.updateTn(e1,t,e2,r,stepsize)
            if (e1!=e1prev or e2!=e2prev) and e1!=-1 and e2!=-1:
                edgechange = True
                tcheck = t;rcheck=r
                if e1!=e1prev:
                    self.its[e1prev].chop(mplast[2],self.its[e1prev].end())
                    if t>self.its[e1].begin():
                        self.its[e1].chop(self.its[e1].begin(),t)
                    #tcheck = self.its[e1].begin()
                    #_,rcheck,mp0 = self.footprint(e1,tcheck,e2)
                if e2!=e2prev:
                    self.its[e2prev].chop(self.its[e2prev].begin(),mplast[4])
                    if r<self.its[e2].end():
                        self.its[e2].chop(r,self.its[e2].end())
                dist = np.linalg.norm(self.p(self.edges[e1],t)-mp)
                mplast = [mp,e1,t,e2,r,dist]
                countlast = count
            if self.debug: print 'points: #{0} {1},{2} {3},{4},{5}'.format(count,e1,e2,t,r,mp)
            if e1 == -1:# or t<0.0 or t>1.0 or r<0.0  or r>1.0:
                if self.debug: print "breaking ",e1,r,t
                break
            dist = np.linalg.norm(mp-self.p(self.edges[e1],t))
            mpcurr = [mp,e1,t,e2,r,dist]
            points.append(mpcurr)
            if dist>MinDistBeforeDistCheck and ( noofsteps>5 or edgechange):# or noofsteps>2 :
                if self.debug: print "now doing distance check"
                negdist = False
                while self.negdistcheck(mp):
                    mpcurr = points.pop(-1)
                    mp = mpcurr[0]
                    negdist = True
                if negdist:
                    return

                pe,mindist,mink,dist = self.distCheck(e1,t,e2,r,mp,stepsize)
                noofsteps = 0
                e1prev,tprev,e2prev,rprev = mplast[1:5]
                if len(pe) >0 :
                    if self.debug: print "dist failure with ",pe
                    del points[countlast-count:0]
                    mps = self.retract(mpcurr,mplast,mink,mindist)
                    mpfn = mps.pop(-1)
                    #self.manage_interval(savedmplast,mpfn)
                    print "retracted..first",mps
                    e1,t,e2,r = mpfn[1:5]
                    points.append(mpfn)
                    kpoints = {}
                    stroke_count = len(self.maxis)-1
                    self.stroke_branches[stroke_count]=[]
                    for k,mp1 in enumerate(mps):
                        kpoints[k]=[mpfn];self.maxis.append(kpoints[k]);self.stroke_branches[stroke_count].append(len(self.maxis)-1)
                    for k,mp1 in enumerate(mps):
                        mp1[0] = mpfn[0]
                        self.iterate(mp1[1],mp1[2],mp1[3],mp1[4],mp1,kpoints[k],stepsize=stepsize,counts=counts)
                    #if len(mps)>0:
                    #    maxis.append(points)
                        #points = []
                    break
                else:
                    self.manage_interval(mplast,mpcurr)
                    pass
                mplast = mpcurr
                countlast = count
        #self.maxis.append(points)
        #self.stroke_count += 1




        #return points\
    def iteraten(self,mplast,parent=-1,stepsize=0.05,counts=10,MinDistBeforeDistCheck=20.0):
        e1,t,e2,r = mplast[1:5]
        if e1==e2:
            return
        if self.debug: print "beginning",e1,t,e2,r,counts
        e1in,tin,e2in,rin = e1,t,e2,r
        count = -1
        noofsteps = 0
        countlast = -1
        points = [mplast]
        #points.append(mplast)
        self.maxis.append(points)
        currentindex = len(self.maxis)-1
        if parent !=-1:
            if parent in self.stroke_branches:
                self.stroke_branches[parent].append(currentindex)
            else:
                self.stroke_branches[parent] = [currentindex]
        savedmplast = mplast
        #mplastn = self.footprintIt(edges[e1],t,edges[e2],its[e1],its[e2],r)
        #mplast = [mplastn,e1,t,e2,r]
        while count < counts:
            count += 1
            noofsteps += 1
            e1prev,tprev,e2prev,rprev = e1,t,e2,r
            edgechange = False
            e1,t,e2,r,mp = self.updateTn(e1,t,e2,r,stepsize)
            if (e1!=e1prev or e2!=e2prev) and e1!=-1 and e2!=-1:
                edgechange = True
                tcheck = t;rcheck=r
                if e1!=e1prev:
                    self.its[e1prev].chop(mplast[2],self.its[e1prev].end())
                    if t>self.its[e1].begin():
                        self.its[e1].chop(self.its[e1].begin(),t)
                    #tcheck = self.its[e1].begin()
                    #_,rcheck,mp0 = self.footprint(e1,tcheck,e2)
                if e2!=e2prev:
                    self.its[e2prev].chop(self.its[e2prev].begin(),mplast[4])
                    if r<self.its[e2].end():
                        self.its[e2].chop(r,self.its[e2].end())
                dist = np.linalg.norm(self.p(self.edges[e1],t)-mp)
                mplast = [mp,e1,t,e2,r,dist]
                countlast = count
            if self.debug: print 'points: #{0} {1},{2} {3},{4},{5}'.format(count,e1,e2,t,r,mp)
            if e1 == -1:# or t<0.0 or t>1.0 or r<0.0  or r>1.0:
                if self.debug: print "breaking ",e1,r,t
                break
            dist = np.linalg.norm(mp-self.p(self.edges[e1],t))
            mpcurr = [mp,e1,t,e2,r,dist]
            points.append(mpcurr)
            if dist>MinDistBeforeDistCheck and ( noofsteps>5 or edgechange):# or noofsteps>2 :
                if self.debug: print "now doing distance check"
                negdist = False
                while self.negdistcheck(mp):
                    mpcurr = points.pop(-1)
                    mp = mpcurr[0]
                    negdist = True
                if negdist:
                    return

                pe,mindist,mink,dist = self.distCheck(e1,t,e2,r,mp,stepsize)
                noofsteps = 0
                e1prev,tprev,e2prev,rprev = mplast[1:5]
                if len(pe) >0 :
                    if self.debug: print "dist failure with ",pe
                    del points[countlast-count:0]
                    mps = self.retract(mpcurr,mplast,mink,mindist)
                    mpfn = mps.pop(-1)
                    #self.manage_interval(savedmplast,mpfn)
                    print "retracted..first",mps
                    e1,t,e2,r = mpfn[1:5]
                    points.append(mpfn)
                    #self.maxis.append(points)
                    #kpoints = {}
                    #stroke_count = len(self.maxis)-1
                    #self.stroke_branches[stroke_count]=[]
                    #for k,mp1 in enumerate(mps):
                    #    kpoints[k]=[mpfn];self.maxis.append(kpoints[k]);self.stroke_branches[stroke_count].append(len(self.maxis)-1)
                    for k,mp1 in enumerate(mps):
                        mp1[0] = mpfn[0]
                        #self.iterate(mp1[1],mp1[2],mp1[3],mp1[4],mp1,kpoints[k],stepsize=stepsize,counts=counts)
                        self.branches.append([currentindex,mp1])
                    #if len(mps)>0:
                    #    maxis.append(points)
                        #points = []
                    break
                else:
                    self.manage_interval(mplast,mpcurr)
                    pass
                mplast = mpcurr
                countlast = count
        #self.maxis.append(points)
        #self.stroke_count += 1




        #return points\

    def analyse(self):
        preve1='';preve2=''
        for ed in self.maxis:
            for pt in ed:
                print pt[1],pt[3]," ",pt[2],pt[4]
                #if pt[1]!=preve1:
                #    print pt[1],pt[3];preve1=pt[1]
                #if pt[3]!=preve2:
                #    print pt[1],pt[3];preve2=pt[3]

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
        e1 = mplast[1];e2 = mplast[3];mpn = mplast[0]
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
    def test2(self,e1,t,e2,r,pts):
        if e2 in self.reflex_edges:
            e2r = self.reflex_objects[e2]
            q1 = e2r.footprint(e1,-1,t)
            mp = q1[2]
        else:
            q1 = self.footprint(e1,t,e2)
            mp = q1[2]
        dist = np.linalg.norm(mp-self.p(self.edges[e1],t))
        mpn = [mp,e1,t,e2,r,dist]

        self.iterate(mpn[1],mpn[2],mpn[3],mpn[4],mpn,pts,0.001,500)

    def test1(self,stepsize=0.01,noofpoints=10000):
        mps = []
        self.maxis=[]
        for cv in self.convex_vert:
        #cv=self.convex_vert[0]
            mps.append([self.vertices[cv],cv,0.0,self.prevedge(cv),1.0,0.0])
        while len(mps) > 0:
            mp = mps.pop(-1)
            points = [mp];self.maxis.append(points)
            self.iterate(mp[1],mp[2],mp[3],mp[4],mp,points,stepsize,noofpoints)
            #print len(points)
        '''
        for i,edge1 in enumerate(self.edges):
            if self.its[i].end()-self.its[i].begin()<=0.001:
                continue
            for j in xrange(i,len(self.edges),1):
                if self.its[j].end()-self.its[j].begin()<=0.001:
                    continue
                e1 = i;e2 =j;t=self.its[e1].begin();r=self.its[e2].end()
                q1 = self.footprint(self.edges[e1],t,self.edges[e2])
                if q1[1]!=-1:
                    mp = q1[2]
                    dist = np.linalg.norm(self.p(self.edges[e1],t)-mp)
                    mpn = [mp,e1,t,e2,r,dist]

                    pts1=self.iterate(mpn[1],mpn[2],mpn[3],mpn[4],mpn,pts,0.01,500)
        return points,self.edges
        '''
    def test1n(self,stepsize=0.01,noofpoints=10000):
        mps = []
        self.maxis=[]
        self.branches=[]
        for cv in self.convex_vert:
        #cv=self.convex_vert[0]
            self.branches.append([-1,[self.vertices[cv],cv,0.0,self.prevedge(cv),1.0,0.0]])
            count = 0
            while len(self.branches) > 0:
                mp = self.branches.pop(0)
            #points = [mp];self.maxis.append(points)
                self.iteraten(mp[1],mp[0],stepsize,noofpoints)
                count += 1
            #print len(points)
        '''
        for i,edge1 in enumerate(self.edges):
            if self.its[i].end()-self.its[i].begin()<=0.001:
                continue
            for j in xrange(i,len(self.edges),1):
                if self.its[j].end()-self.its[j].begin()<=0.001:
                    continue
                e1 = i;e2 =j;t=self.its[e1].begin();r=self.its[e2].end()
                q1 = self.footprint(self.edges[e1],t,self.edges[e2])
                if q1[1]!=-1:
                    mp = q1[2]
                    dist = np.linalg.norm(self.p(self.edges[e1],t)-mp)
                    mpn = [mp,e1,t,e2,r,dist]

                    pts1=self.iterate(mpn[1],mpn[2],mpn[3],mpn[4],mpn,pts,0.01,500)
        return points,self.edges
        '''

    def testcv(self,cv,points,stepsize=0.01,noofpoints=10000):
        assert cv in self.convex_vert
        mp=self.vertices[cv],cv,0.0,self.prevedge(cv),1.0,0.0
        pts = self.iterate(mp[1],mp[2],mp[3],mp[4],mp,points,stepsize,noofpoints)
        return pts
    def ppold(self,maxis): #post processing
        maxis1 = copy.copy(maxis)
        for ed in maxis1:
            if ed[0][1] in self.reflex_edges or not ed[0][3] in self.reflex_edges:
                continue
            rad = ed[0][-1]
            for pt in ed:
                if pt[3] in self.reflex_edges:
                    pp = self.pprime(self.edges[pt[1]],pt[2])
                    if len(pp)==0:
                        continue
                    n = np.array([pp[1],-1*pp[0]])
                    n = n/np.linalg.norm(pp)
                    p = self.p(self.edges[pt[1]],pt[2])
                    mpnew = p + rad*n
                    pt[0] = mpnew
        return maxis1





    def test3(self,letter,pts):
        fl = fp(letter)
        try:
            pts1=fl.test1(pts)
        except:
            q1=fl.printMa(pts);check.printEdges(q1[0]+fl.edges,0,False)
            check.printEdges(fl.edges)
            check.printEdges(q1[1]+fl.edges,0,False)


    def distCheck(self,e1,t,e2,r,mp,stepsize=0.01):

        #r is not used at all
        #dist = max(np.linalg.norm(self.p(self.edges[e1],t)-mp),np.linalg.norm(self.p(self.edges[e2],r)-mp))
        diste1,_,_= self.distanceIt(e1,mp,stepsize)
        diste2,_,_= self.distanceIt(e2,mp,stepsize)
        dist = min(diste1,diste2)
        if dist == float('inf'):
            print "infinite distance ",e1,e2,t,r,diste1,diste2,mp
            return [],float('inf'),-1,float('inf')
        pe = []
        mindist3 = float('inf')
        mink = -1
        for k,_ in enumerate(self.edges):
            #if k in self.friends(e1,cnts,convex_vert) or k in self.friends(e2,cnts,convex_vert):
            if k==e1 or k==e2:
            #or k==self.nextedge(e1,cnts) or k==self.prevedge(e2,cnts) or k==self.prevedge(e1,cnts) or k==self.prevedge(e2,cnts):
            # or k in tracededges:
                continue
            dist3,p,t = self.distanceIt(k,mp,stepsize)
            if dist3>0 and abs(dist3)< dist:
                #print k,dist3,dist

                pe.append([k,dist3,t])
                if abs(dist3)<mindist3:
                    mink = k
                    mindist3 = dist3
        return pe,mindist3,mink,dist
    def retract(self,mpcurr,mplast,mink,mindist):
          if self.debug: print "retract begin mplast,mpcurr ",mplast,mpcurr
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
          distk,p1,r0 = self.distance(mink,mpcurr[0])
          tf = mplast[2];rf=mplast[4]
          mpn = mplast[0]
          while abs(distk-dist)>0.00001  and count < 20 and tf!=tl:
                count += 1
                tn = (tf + tl)/2.0
                rn = (rf+rl)/2.0
                _,rn0,mpn,x = self.footprintIt(e1,tn,e2)
                if rn0>=0.0 and rn0<=1.0:
                    rn = rn0
                    distk,p1,r0 = self.distanceIt(mink,mpn)
                    dist = np.linalg.norm(mpn-self.p(self.edges[e2],rn))
                else:
                    _,tn,mpn,x = self.footprintIt(e2,rn,e1)
                    if self.debug2: print "e2,rn,e1,tn",e2,rn,e1,tn
                    r0=-1
                    if x!=None:
                        distk,p1,r0 = self.distanceIt(mink,mpn)
                        dist = np.linalg.norm(mpn-self.p(self.edges[e1],tn))

                    if r0 == -1:
                        print mplast,mpcurr
                        break
                        #raise Exception("do not how to proceed from here")

                if distk < dist:
                        tl = tn
                else:
                        tf = tn
          if rn!=-1:
                mpfn = [mpn,e1,tn,e2,rn,dist]
          else:
               mpfn = mpcurr
          if abs(distk-dist)>0.0001:
               print "possible bug in retraction, new dist difference is ",abs(distk-dist),mpcurr,mplast
          mps = []
           #mps.append([mpfn[0],mink,])
          #distk,fpk,tk = self.distanceIt(mink,mpn)
          #_,tk,mpk = self.footprint(mpfn[1],mpfn[2],mink)
          #_,rtk,mpk1=self.footprint(mink,tk,mpfn[3],mink)
          #mpfn[4] = rtk
          #try: assert np.linalg.norm(mpk-mpn)<2.0 and np.linalg.norm(mpk1-mpn)<2.0
          #except:   print "in retract:",abs(tk-tk1),np.linalg.norm(mpk-mpn),np.linalg.norm(mpk1-mpn);#raise Exception('possible error in retract')

          mpfn1 = [mpn,mpfn[1],mpfn[2],mink,r0,distk]
          mpfn2 = [mpn,mink,r0,mpfn[3],mpfn[4],distk]
          mps.append(mpfn1);mps.append(mpfn2)


          '''
          for i,_ in enumerate(self.edges):
               if i==mpfn[1] or i==mpfn[3]:
                   continue
               disti,fpi,ti = self.distanceIt(i,mpn)
               if abs(disti-distk)<2.0:
                   mpfn1 = [mpn,mpfn[1],mpfn[2],i,ti,disti]
                   mpfn2 = [mpn,i,ti,mpfn[3],mpfn[4],disti]
                   mps.append(mpfn1);mps.append(mpfn2)
                      #return tn,rn,-1
          '''
          mps.append(mpfn)
          self.manage_interval(mplast,mpfn)
          return mps
    def modEdges(self,letter,fontfile="DevanagariSangamMN.ttf"):
        self.edges,self.vertices,self.vert_edges,self.cnts = check.fontOutline(letter,fontfile)
        self.vertAnalysen()
        self.its = {} # store as interval tree
        for k in xrange(len(self.edges)+len(self.concave_vert)):
            self.its[k]=iv(0.0,1.001)
        #its =
        self.insertReflexn()
        for k in xrange(len(self.edges)):
            self.its[k] = it([self.its[k]])
        #self.vertAnalysen()
        return self.edges,self.its,self.cnts,self.vertices,self.convex_vert,self.concave_vert
    def insertReflexn(self):#
        concave_vert = copy.copy(self.concave_vert)
        concave_vert.sort()
        self.reflex_edges = []
        #print cnts
        cnts1 = copy.copy(self.cnts)
        while len(concave_vert) >0:
            #print concave_vert
            i = concave_vert.pop(-1)
            #print "inserting",i
            e2 = i
            e1 = self.prevedge(i)
            e1prime = self.p(self.edges[e1],0.6)
            e2prime = self.p(self.edges[e2],0.4)
            quad = []
            #quad.append(e1prime)
            #print e1prime
            quad.append(self.vertices[i])
            #quad.append(e2prime)
            #print e1,e2
            #lastts[e1][1] = 0.999
            #edges[e1][2] = e1prime
            #edges[e2][0] = e2prime
            self.vertices.insert(i,self.vertices[i])
            #self.vertices[i+1] = e2prime
            #lastts[i+1] = [0.001,1.0]

            self.edges.insert(i,quad)
            #self.reflex_edges.append(i)
            ind1 = 0
            ind3 = 0
            #print "before adding",i
            #print self.its

            for ind2,cnt in enumerate(self.cnts):
                for ind,k in enumerate(cnt):
                    if k > i:
                        cnt[ind] = k+1
                    elif k ==i:
                        ind1 = ind
                        ind3 = ind2
            for k in  xrange(len(self.edges)-1,i-1,-1):
                if k >=i:
                    self.its[k+1] = self.its[k]
            self.its[i]=iv(0.0,1.001)
            lowerbound = self.its[e1].begin
            #print e1,lowerbound,"lb"
            self.its[e1] = iv(lowerbound,1.001)
            upperbound = self.its[i+1].end
            #print i+1,upperbound,"ub"
            self.its[i+1] =iv(0.0,upperbound)

            self.cnts[ind3].insert(ind1+1,i+1)
            #print "after adding",i
            #print self.its

            for j,k in enumerate(self.convex_vert):
                if k>i:
                    self.convex_vert[j] = k+1
        for k,edge in enumerate(self.edges):
            if len(edge)==1:
                self.reflex_edges.append(k)


        #print cnts
        #print its
        return self.its
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
    def pp(self,q):
         def offset(e1,t,offset):
             #print e1,t,offset
             #n = np.array(fka.normal(fka.edges[e1],t,True))
             pprime = self.pprime(self.edges[e1],t)
             norm = np.linalg.norm(pprime)
             normal = np.array([pprime[1],-pprime[0]])
             n = normal/norm
             p = self.p(self.edges[e1],t)
             return p + offset*n
         rad=60.0
         rads={}
         #c.printEdges(fka.edges)
         for j,ed in enumerate(q):
           preve2='';preve1='';rad=61.0
           for i,pt in enumerate(ed):
                if pt[3] in self.reflex_edges and pt[1] not in self.reflex_edges:
                   if pt[3]!=preve2 and pt[1]==preve1:
                       preve2=pt[3];
                       if i>0:
                           print "rad here"
                           rad = ed[i-1][-1];
                       else:
                           rad = pt[-1]
                       preve1=pt[1]
                       print "rad for e1,e2=",pt[1],pt[3],rad,j
                   else:
                       pt[0] = offset(pt[1],pt[2],rad);preve2 = pt[3];preve1=pt[1]
                else:
                       preve2 = pt[3];preve1=pt[1]
           for i in xrange(len(ed)-1,-1,-1):
               pt = ed[i]
               if pt[1] in self.reflex_edges and pt[3] not in self.reflex_edges:
                  if pt[1]!=preve1 and pt[3]==preve2:
                      preve2=pt[3];preve1=pt[1]
                      if i<len(ed)-1:
                          print "rad here"
                          rad = ed[i+1][-1]
                      else:
                          rad = pt[-1];
                      print "rad for e1,e2=",pt[1],pt[3],rad,j
                  else:
                      pt[0] = offset(pt[3],pt[4],rad);preve2 = pt[3];preve1=pt[1]
               else:
                      preve2 = pt[3];preve1=pt[1]
           print "-----"
         pointindices=[]
         for i,ed in enumerate(q):
             if len(ed)==1: pointindices.append(i)
             if i<len(q)-1:
                 if np.linalg.norm(ed[-1][0]-q[i+1][0][0])<30.0:
                     ed[-1][0]=q[i+1][0][0]
         q1 = [q[i] for i,ed in enumerate(q) if i not in pointindices]

         #q1 =fka.printMa(q1);
         return q1
         c.printEdges(q1[1]+fka.edges,0,False)
         c.printEdges(q1[1])
         c.printEdges(q2[1]+fka.edges,0,False)




    def footprintIt(self,e1,t,e2,initialr=-1):#using interval tree
        if e1<0:
            assert initialr!=-1
            return footprintneg(e2,initialr,e1)
        elif e2<0:
            return footprintneg(e1,t,e2)
        if( self.its[e1].search(t) == set()):
            return [],-1,[],None

        if e1 in self.reflex_edges: # it is reflex object
            ref = self.reflex_objects[e1]
            e,r,mp=ref.footprintn(e2,t,-1)
        elif e2 in self.reflex_edges:
            ref = self.reflex_objects[e2]
            e,r,mp=ref.footprintn(e1,-1,t)
        else:
            e,r,mp = self.footprint(e1,t,e2,initialr)
        x = self.its[e2].search(r)
        if x!=set():
            try: assert len(x) == 1 # if not then we have not managed intervals better
            except: print e2,self.its[e2],r,x;raise Exception('failure in managing interval')
            return e,r,mp,x
        else:
            return e,r,mp,None


    def bezIntersect(self,bez1,bez2):
        return
    def bisectorn(self,e1,e2,counts=100):
        finedges = []
        points = []
        for q in xrange(counts+1):
            t = q/(counts*1.0)
            if e1 in self.reflex_edges:
                e1r = self.reflex_objects[e1]
                edges,r,point= e1r.footprintn(e2,t,-1)
            elif e2 in self.reflex_edges:
                e2r = self.reflex_objects[e2]
                edges,r,point= e2r.footprintn(e1,-1,t)
            else:
                quad1 = self.edges[e1];quad2=self.edges[e2]
                edges,r,point = self.footprintn(quad1,t,quad2)



            #r,edges = self.solve2(quad1,t,quad2)
            #point = edges[0][0]
            #print "t=,r=",t,r


            if r>=-0.0001 and r<=1.0001:
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



    def bisector(self,quad1,quad2,counts=100):
        finedges = []
        points = []
        for q in xrange(counts+1):
            t = q/(counts*1.0)
            edges,r,point = self.footprint(quad1,t,quad2)
            #r,edges = self.solve2(quad1,t,quad2)
            #point = edges[0][0]
            #print "t=,r=",t,r


            if r>=-0.0001 and r<=1.0001:
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
    def printMa(self,ma,threshold=10):
        def color(k):
            return (k*25%255)/255.0,(k*40%255)/255.0,(k*10%255)/255.0
        points = []
        bzs = []
        man = self.pp(ma)
        indexedbzs=[]
        for k,ed in enumerate(man):
            crvs = [];offsetvect = np.array([0.0,0.0])

            for m,point in enumerate(ed):
                if point != None:
                    if (point[1] in self.convex_vert and point[3] == self.prevedge(point[1])):
                        continue
                    if (point[3] in self.convex_vert and point[1] == self.prevedge(point[3])):
                        continue

                    if len(point)>0:
                        if m>=2:
                            vect1 = point[0]-ed[m-1][0];vect2=ed[m-1][0]-ed[m-2][0]
                            theta1=np.arctan2(vect1[1],vect1[0]);theta2 = np.arctan2(vect2[1],vect2[0])
                            if np.tan(theta2-theta1) >threshold:
                                print np.tan(theta2-theta1),threshold,theta1,theta2,m
                                offsetvect = -1*vect2

                        pt = point[0]-offsetvect
                        points.append(pt)
                        crvs.append(pt)
            if len(crvs)>2:
                #pass
                crvs = self.smoothPoints(crvs)
                #print "began",crvs
                #crvs = sgolay2d(np.array(crvs), window_size=11, order=2)
                #crvs=map(np.array,zip(*crvs))
                #print "finished",crvs
                if len(crvs)>2:
                    bz = fitCurve(crvs,10.0);indexedbzs.append(bz)
                    bz0=[];offset=0
                    for l,bzc in enumerate(bz):
                        if l>0: offset = np.linalg.norm(bz[l][0]-bz[l-1][-1]);offsetvec1 = bz[l-1][-1]-bz[l][0]
                        if offset>0:
                            print "doing offset";
                            for bzci,_ in enumerate(bzc):
                                bzc[i] = bzc[bzci] + offsetvec1
                        bz0 += [[bzc,color(k),0,k]]
                        bzs += bz0
            elif len(crvs)==2:
                bzc = [[np.array([crvs[0],(crvs[0]+crvs[1])/2.0,crvs[1]])]];
                indexedbzs.append(bzc);
                bzc=[bzc,color(k),0,k]
                bzs += bzc


        return points,bzs,indexedbzs
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
    def vertAnalysen(self):
        print len(self.vertices),len(self.edges)
        if len(self.vertices)==len(self.edges):
            self.concave_vert = []
            self.convex_vert = []
        else:
            return self.convex_vert,self.concave_vert



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
    def pts(self,pts):
        pts1 = []
        for pt in pts:
            pts1.append([pt[0]])
        check.printEdges(self.edges+pts1,0,False)
    def verifyFootprint(self):
        for i,edge1 in enumerate(self.edges):
            #p1 = self.p(edge,1.0)
            rs = []
            for k,edge2 in enumerate(self.edges):
                if k==i:
                    continue
                _,r,fp = self.footprint(edge1,1.0,edge2)
                if r!=-1:
                    rs.append([r,i,k])
            if len(rs)>1:
                print "error?:",rs
    def bisectorr(self,e1,e2):
        assert e1 in self.reflex_edges and e2 in self.reflex_edges
        e1r = self.reflex_objects[e1];e2r = self.reflex_objects[e2];
        pts = []
        ep=[]
        for q in xrange(21):
            t = q/40.0
            q1 = e1r.footprintn(e2,t,-1)
            ep=ep+q1[0]
            if q1[1]!=-1:
                pts.append([q1[-1]])
        check.printEdges(pts+self.edges+ep,0,False)


    '''
    Based on paper by Klimenko on Handwritten fonts modelling
    '''
    def smoothPoints(self,pts):
        pts1 = []
        k1 = 1/3.0;k2=1/3.0;k3=1/3.0
        #print len(pts)
        for i in xrange(1,len(pts)-2,1):
            #print pts[i-1][0]
            pt = np.array(k1*pts[i-1] + pts[i]*k2 + pts[i+1]*k3)
            #print pt

            pts1.append(pt)
        return np.array(pts1)
    '''
    n = n1*(1-t) + n2*t or t =norm(n-n1)/norm(n2-n1)
    Proposed if edge as e1 and reflext object as e2:
       radius of e1 and next edge of e2 at 0.0 till midway that is t =0.5 or n1/2 and n2/2
    so as t is increased we go on finding n and hence t if t becomes less than 1/2 we stop
    similar method to the  previous is applied if it is edge
    '''
class Reflex(object):
        def p(self,t):
            return self.p
        def pprime(self,t):
            theta = self.theta1*(1-t)+self.theta2*t
            return np.array([np.cos(theta+np.pi/2),np.sin(theta+np.pi/2)])
        def n(self,t):
            theta = self.theta1*(1-t)+self.theta2*t
            return np.array([np.cos(theta),np.sin(theta)])

        def updateo(self,fp,vert_ind):
            p2prime = fp.pprime(fp.edges[fp.nextedge(vert_ind)],0.0)
            p1prime = fp.pprime(fp.edges[fp.prevedge(vert_ind)],1.0)
            n1 = np.array([p1prime[1],-1*p1prime[0]])
            n2 = np.array([p2prime[1],-1*p2prime[0]])

            self.theta2 = np.arctan2(n2[1],n2[0])
            if self.theta2<0:
                self.theta2 += 2*np.pi
            self.n1 = n1;self.n2=n2
            self.theta1 = np.arctan2(n1[1],n1[0])
            if self.theta1<0:
                self.theta1 += 2*np.pi

            self.p = fp.vertices[vert_ind]
            self.fp = fp
            self.i = vert_ind
            self.rads={}
            e2 = fp.nextedge(vert_ind)
            e1 = fp.prevedge(vert_ind)
            #print fp.reflex_edges
            savedrad1 = float('inf');savedrad2=float('inf');savede1=-1;savede2=-2
            savedmp1=[];savedmp2=[];savedr1=-1;savedr2=-1

            for i,edge in enumerate(fp.edges):
                if i in fp.reflex_edges or i ==fp.nextedge(vert_ind) or i==fp.prevedge(vert_ind):
                    continue
                #print "i,e1,e2",i,e2,e1
                rad1=rad2=float('inf');mp1=mp2=[];r1=r2=-1

                _,r,mp = fp.footprint(e2,0.0,i)
                #print "with e2=0.0",e2,r
                if r!=-1:
                    rad2 = np.linalg.norm(mp-fp.p(fp.edges[e2],0.0))
                    if rad2<savedrad2:
                        savedrad2 = rad2
                        savedmp2=mp2
                        savedr2 =r
                        savede2 =i


                    mp2 = mp
                    r2=r;

                _,r,mp = fp.footprint(e1,1.0,i)
                #print "with e1=1.0",e2,r

                if r!=-1:
                    rad1 = np.linalg.norm(mp-fp.p(fp.edges[e1],1.0))
                    if rad1<savedrad1:
                        savedmp1=mp
                        savede1 = i
                        savedrad1=rad1
                        savedr1=r



                    mp1 = mp
                    r1=r;

                self.rads[i]=[r1,r2,rad1,rad2,mp1,mp2]
            self.rad=[savedrad1,savedrad2,savedmp1,savedmp2,savede1,savede2,savedr1,savedr2]
        def update(self,fp,vert_ind):
            p2prime = fp.pprime(fp.edges[fp.nextedge(vert_ind)],0.0)
            p1prime = fp.pprime(fp.edges[fp.prevedge(vert_ind)],1.0)
            n1 = np.array([p1prime[1],-1*p1prime[0]])
            n2 = np.array([p2prime[1],-1*p2prime[0]])

            theta2 = np.arctan2(n2[1],n2[0])
            #if theta2<0:
            #    theta2 += 2*np.pi
            self.n1 = n1;self.n2=n2
            theta1 = np.arctan2(n1[1],n1[0])
            self.thetae=[theta1,theta2]
            self.thetaflag=False #if true add 2*np.pi if angle <0
            #if theta1<0 and theta2>0 and theta2-theta1>np.pi:
            #    theta1 += 2*np.pi;self.thetaflag=True
            if theta2<0  and theta1>0:
                theta2 += 2*np.pi;self.thetaflag=True
            if theta2<0 and theta1<0 and abs(theta2)>abs(theta1):
                theta2 += 2*np.pi;self.thetaflag=True
            self.theta1 = theta1
            self.theta2 = theta2
            self.degree1 = theta1*180/np.pi;self.degree2 = theta2*180/np.pi
            self.p = fp.vertices[vert_ind]
            self.fp = fp
            self.i = vert_ind
            self.rads={}
        def updateRad(self,fp,vert_ind):
            e2 = fp.nextedge(vert_ind)
            e1 = fp.prevedge(vert_ind)
            #print fp.reflex_edges
            savedrad1 = float('inf');savedrad2=float('inf');savede1=-1;savede2=-2
            savedmp1=[];savedmp2=[];savedr1=-1;savedr2=-1

            for i,edge in enumerate(fp.edges):
                if i ==fp.nextedge(vert_ind) or i==fp.prevedge(vert_ind):
                    continue
                #print "i,e1,e2",i,e2,e1
                rad1=rad2=float('inf');mp1=mp2=[];r1=r2=-1

                _,r,mp = self.footprintn(i,1.0,-1)
                #print "with e2=0.0",e2,r
                if r>=0.0 and r<=1.0:
                    rad2 = np.linalg.norm(mp-self.p)
                    if rad2<savedrad2:
                        savedrad2 = rad2
                        savedmp2=mp2
                        savedr2 =r
                        savede2 =i


                    mp2 = mp
                    r2=r;

                _,r,mp = self.footprintn(i,0.0,-1)
                #print "with e1=1.0",e2,r

                if r>=0.0 and r<=1.0:
                    rad1 = np.linalg.norm(mp-self.p)
                    if rad1<savedrad1:
                        savedmp1=mp
                        savede1 = i
                        savedrad1=rad1
                        savedr1=r



                    mp1 = mp
                    r1=r;

                self.rads[i]=[r1,r2,rad1,rad2,mp1,mp2]
            self.rad=[savedrad1,savedrad2,savedmp1,savedmp2,savede1,savede2,savedr1,savedr2]

        def solve(self,t,e0):
            def f(r):
                p = ft.p(e[e0],r)
                pprime = ft.pprime(e[e0],r)
                norm1 =np.linalg.norm(pprime)
                n = np.array([pprime[1],-1*pprime[0]])

                #print p,pprime,n
                mp = p + rad*n*1/norm1
                ep.append([v[self.i],mp])
                d1 = mp-v[self.i]
                return np.cross(d1,vectr)
            rs=[]
            ep=[]
            ft = self.fp
            e = self.fp.edges
            v = ft.vertices
            #r0=self.rad[e0][1];rad=self.rad[e0][0]
            vectr = self.n1*(1-t)+self.n2*t
            e0prev = self.fp.prevedge(e0);e0next = self.fp.nextedge(e0)
            rad = -1
            if self.rads[e0][0]!=-1:
                rad = self.rads[e0][2]
            elif self.rads[e0][1]!=-1:
                rad = self.rads[e0][3]
            elif self.rads[e0prev][0]!=-1:
                rad = self.rads[e0prev][2]
            elif self.rads[e0prev][1]!=-1:
                rad = self.rads[e0prev][3]
            elif self.rads[e0next][0]!=-1:
                rad = self.rads[e0next][2]
            elif self.rads[e0next][1]!=-1:
                rad = self.rads[e0next][3]
            if rad==-1:
                return ep,-1,[]




            if t>=0.5:
                rad=self.rad[1]
                rs.append(0.0)
                rs.append(1.0)
            else:
                rad = self.rad[0]
                rs.append(1.0)
                rs.append(0.0)


            fxi = f(rs[-1])
            fximinus1 = f(rs[-2])
            print fxi,fximinus1
            rn = rs[-1]
            #return
            for i in xrange(20):
                if abs(rs[-1]-rs[-2])<0.001 or abs(fxi-fximinus1)<0 or abs(fxi)<0.001:

                    break
                ri = rs[-1]
                rn = ri -fxi *(rs[-1]-rs[-2])/(fxi-fximinus1)
                #print rn
                rs.append(rn)
                fximinus1 = fxi
                fxi = f(rn)
            if rn>=0 and rn<=1:
                p = ft.p(e[e0],rn)
                n = ft.normal(e[e0],rn,True)
                mp = p + rad*n
                ep.append([[v[self.i],mp],'r'])
                ep.append([[p,mp],'g'])
                return ep,rn,mp
            else:
               return ep,-1,[]

        def footprint(self,other,t,r,rad=67.0):#t is reflex edge and r is of nonedge
            ep=[]
            if other in self.fp.reflex_edges:
                refother = self.fp.reflex_objects[other]
                if t!=-1:
                    p1p2 = refother.p-self.p
                    d = np.linalg.norm(p1p2)/2.0
                    n = self.n(t)
                    costheta = np.dot(n,p1p2)/np.linalg.norm(p1p2) #n is unit length
                    d1 = d/costheta
                    mp = self.p + d1*n
                    vect = mp-refother.p
                    thetan = np.arctan2(vect[1],vect[0])
                    if thetan<0 and refother.flag:
                        thetan = 2*np.pi+thetan

                    ep.append([self.p,mp])
                    ep.append([mp])
                    ep.append([refother.p,mp])
                    r = (thetan-refother.theta1)/(refother.theta2-refother.theta1)
                    return ep,r,mp
                    #pass
                else:
                    p1p2 = self.p-refother.p
                    d = np.linalg.norm(p1p2)/2.0
                    n = refother.n(r)
                    costheta = np.dot(n,p1p2)/np.linalg.norm(p1p2) #n is unit length
                    d1 = d/costheta
                    mp = refother.p + d1*n
                    vect = mp-self.p
                    thetan = np.arctan2(vect[1],vect[0])
                    #if thetan<0:
                    #    thetan = 2*np.pi+thetan
                    #print thetan,self.theta1,self.theta2
                    t = (thetan-self.theta1)/(self.theta2-self.theta1)
                    ep.append([self.p,mp])
                    ep.append([mp])
                    ep.append([refother.p,mp])
                    return ep,t,mp

                    #pass

            elif t!=-1:
                return self.solve(t,other)
            else:
                p = self.fp.p(self.fp.edges[other],r)
                ep.append([p])

                pprime = self.fp.pprime(self.fp.edges[other],r)
                normal = np.array([pprime[1],-1*pprime[0]])
                norm = np.linalg.norm(normal)
                q0=other;found=False
                mp = p +1/norm*normal*rad
                ep.append([p,mp])
                ep.append([self.fp.vertices[self.i],mp])
                vectr = mp - self.fp.vertices[self.i]
                thetan = np.arctan2(vectr[1],vectr[0])
                #if thetan<0:
                #    thetan = thetan + 2*np.pi
                t = (thetan-self.theta1)/(self.theta2-self.theta1)
                return ep,t,mp
                #get unit normal vector
            return
        def footprintn(self,other,t,r,rad=0):
            def f(r):
                p1 = self.fp.p(self.fp.edges[other],r)
                dist0 = np.linalg.norm(p1-self.p)/2.0
                mp = self.p+n1*dist0
                dist1,_,_ = self.fp.distance(other,mp)
                return abs(dist1-dist)

            ep=[]
            if other==self.i:
                return [],-1,[]
            if other in self.fp.reflex_edges:
                refother = self.fp.reflex_objects[other]
                if t!=-1:
                    p1p2 = refother.p-self.p
                    d = np.linalg.norm(p1p2)/2.0
                    n = self.n(t)
                    costheta = np.dot(n,p1p2)/np.linalg.norm(p1p2) #n is unit length
                    if abs(costheta)<0.00001:
                        return [],-1,[]
                    d1 = d/costheta
                    mp = self.p + d1*n
                    vect = mp-refother.p
                    thetan = np.arctan2(vect[1],vect[0])
                    if thetan<0 and refother.thetaflag:
                        thetan = 2*np.pi+thetan

                    ep.append([self.p,mp])
                    ep.append([mp])
                    ep.append([refother.p,mp])
                    r = (thetan-refother.theta1)/(refother.theta2-refother.theta1)
                    return ep,r,mp
                    #pass
                else:
                    p1p2 = self.p-refother.p
                    d = np.linalg.norm(p1p2)/2.0
                    n = refother.n(r)
                    costheta = np.dot(n,p1p2)/np.linalg.norm(p1p2) #n is unit length
                    d1 = d/costheta
                    if abs(costheta)<0.00001:
                        return [],-1,[]
                    mp = refother.p + d1*n
                    vect = mp-self.p
                    thetan = np.arctan2(vect[1],vect[0])
                    if thetan<0 and self.thetaflag:
                        thetan = 2*np.pi+thetan
                    #print thetan,self.theta1,self.theta2
                    t = (thetan-self.theta1)/(self.theta2-self.theta1)
                    ep.append([self.p,mp])
                    ep.append([mp])
                    ep.append([refother.p,mp])
                    return ep,t,mp

                    #pass

            elif t!=-1:
                n1 = self.n(t)
                dist,fp1,r0 =self.fp.distance(other,self.p);dist0=float('inf')
                r = r0;count=0
                while abs(dist0-dist)>0.0001 and count<20:
                    count += 1
                    p1 = self.fp.p(self.fp.edges[other],r)
                    distp1p2 = np.linalg.norm(p1-self.p)
                    dist0 = np.linalg.norm(p1-self.p)/2.0
                    costheta = np.dot(n1,p1-self.p)/distp1p2
                    dist0 = distp1p2/(2.0*costheta)
                    mp = self.p +n1*dist0
                    dist,fp,r=self.fp.distance(other,mp)
                    if dist==float('inf') or dist<0:
                        break
                if dist0!=float('inf') and abs(dist0-dist)<0.0001:
                    ep.append([self.p,mp]);ep.append([p1,mp]);ep.append([mp])
                    return ep,r,mp
                else:
                    return [],-1,[]
            else:
                p = self.fp.p(self.fp.edges[other],r)
                pprime = self.fp.pprime(self.fp.edges[other],r)
                normal = np.array([pprime[1],-1*pprime[0]])
                assert np.dot(pprime,normal)==0
                norm = np.linalg.norm(normal)

                distp1p2 = np.linalg.norm(p-self.p)
                costheta = np.dot(normal,self.p-p)/(distp1p2*norm)
                if abs(costheta)<0.0000001:
                    return [],-1,[]

                d0 = distp1p2/(2.0*costheta)
                mp = p +normal*d0/norm
                ep.append([mp])
                ep.append([p,mp])
                ep.append([self.p,mp])
                ep.append([self.p,self.p+20*self.n1])
                ep.append([self.p,self.p+20*self.n2])

                vectr = mp - self.p
                try: assert abs(np.linalg.norm(vectr)-d0)<0.001
                except: return [],-1,[]
                costheta1 = np.dot(vectr,p-self.p)/(np.linalg.norm(vectr)*distp1p2)
                #print costheta1,costheta
                assert abs(costheta1-costheta)<0.0001
                thetan = np.arctan2(vectr[1],vectr[0])
                if thetan<0 and self.thetaflag:
                    thetan = thetan + 2*np.pi
                #if self.debug: print "in footprint", thetan,self.theta1,self.theta2
                t = (thetan-self.theta1)/(self.theta2-self.theta1)
                return ep,t,mp
                #get unit normal vector
            return
        def printNormals(self):
            ep = []
            for q in xrange(11):
                t = q/10.0
                ep.append([self.p+20*self.n(t),self.p])
                #theta =self.theta1*(1-t) + self.theta2*t
            check.printEdges(ep,0)







        def distance(self,pi):

            #dist,t,pi = self.distance(edge,self.p)
            ep = []
            ep.append([self.p,self.p+self.n1])
            ep.append([self.p,self.p+self.n2])
            ep.append([self.p,pi])
            vect =pi-self.p
            theta = np.arctan2(vect[1],vect[0])
            #if theta<0:
            #    theta += np.pi*2
            t = (theta-self.theta1)/(self.theta2-self.theta1)
            n1 = self.n(0.0);n2=self.n(1.0)
            if np.cross(n1,pi-self.p)>0 and np.cross(pi-self.p,n2)>0:
                return np.linalg.norm(pi-self.p),self.p,t
            else:
                return -1*np.linalg.norm(pi-self.p),self.p,t
