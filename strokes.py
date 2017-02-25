import numpy as np
import operator as op
import math
import bezier
import ttfquery,ttfquery.describe,ttfquery.glyphquery,ttfquery.glyph
from intervaltree import IntervalTree as it
from intervaltree import Interval as iv

from matplotlib import pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
class bz(object):
    @staticmethod
    def ncr(n, r):
        r = min(r, n - r)
        if r == 0: return 1
        numer = reduce(op.mul, xrange(n, n - r, -1))
        denom = reduce(op.mul, xrange(1, r + 1))
        return numer // denom
    @staticmethod
    def B(n,i,u):
        return bz.ncr(n,i)*math.pow(u,i)*math.pow(1-u,n-i)
    @staticmethod
    def p(nodes,u):
        n = len(nodes) - 1
        sum = 0.0
        for i in xrange(n+1):
            sum += bz.B(n,i,u)*nodes[i]
        return sum
    @staticmethod
    def pprime(nodes,u):
        n = len(nodes) - 1
        sum = 0.0
        for i in xrange(n):
            sum += bz.B(n-1,i,u)*n*(nodes[i+1]-nodes[i])
        return sum
    @staticmethod
    def normal(nodes,u,unit=False):
        pp = bz.pprime(nodes,u)
        normpp = np.linalg.norm(pp)
        if unit and normpp!=0:
            return 1/normpp*np.array([pp[1],-1*pp[0]])
        return np.array([pp[1],-1*pp[0]])
    @staticmethod
    def secant(f,rs):
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
    @staticmethod
    def distance(quad,p0,debug=False):
        def f(t):
            return np.dot(bz.pprime(quad,t),p0-bz.p(quad,t))
        rs = [];rs.append(0.4);rs.append(0.6)
        rn = bz.secant(f,rs)
        p1 = bz.p(quad,rn);p1prime=bz.pprime(quad,rn)
        if True or debug:
            ep=[];ep.append([p0]);ep.append([p0,p1]);ep.append([p1]);ep.append([p1,p1+p1prime])
        #print abs(f(rn)),np.cross(p0-p1,p1prime)
        if np.cross(p0-p1,p1prime) >=0 and abs(f(rn))<1e-2: # point is on right side of the side
             return np.linalg.norm(p1-p0),p1,rn,ep
        #elif np.cross(p0-p1,p1prime) <0 and abs(f(rn))<1.0:
            #return -1*np.linalg.norm(p1-p0),p1,rn
        else:
            return float('inf'),p1,rn,ep
    def distancen(nodes,pi):
        def f(t):
            return np.dot(pi-bz.p(nodes,t),bz.pprime(nodes,t))
        def secant(f,x0,xminus):
            fximinus = f(xminus)
            fxi = f(x0)
            count = 0
            xn = x0
            while abs(fxi)>1e-6 and abs(fxi-fximinus)>1e-6 and count<50:
                xn = x0 - fxi * (x0-xminus)*1/(fxi-fximinus)
                fximinus = fxi
                fxi = f(xn)
                count += 1
                #print xn
            if abs(fxi)<1e-6:
                return xn
            else:
                return float('inf')
        t0 = secant(f,0.6,0.4)
        if t0 != float('inf'):
            dist = np.linalg.norm(pi-bz.p(nodes,t0))
        else:
            dist = float('inf')
        return dist,t0
    @staticmethod
    def footprint(quad1, t, quad2):
        p1 = bz.p(quad1, t)
        p1prime = bz.pprime(quad1, t)
        n1 = bz.normal(quad1,t,True)
        dist0, p2i, r1,_ = bz.distance(quad2, p1,True)  # randomly ..no partcular thought
        edges = []
        edges.append([p1]);edges.append([p1, p2i])
        #print dist0,r1,p1
        '''
        if r1 < 0.0:
            # print "here"
            p2i = bz.p(quad2, 0.0);
            dist0 = np.linalg.norm(p2i - p1)
        elif r1 > 1.0:
            p2i = bz.p(quad2, 1.0);
            dist0 = np.linalg.norm(p2i - p1)
        '''


        for i in xrange(20):
            costheta1 = np.dot(p2i - p1, n1) / dist0
            #print costheta1
            if abs(costheta1) > 0.01:
                dist1 = dist0 / (2 * abs(costheta1))
            else:
                dist1 = dist0 / 2.0
            fpi = p1 + dist1 * n1
            dist2, p2i, r1,_ = bz.distance(quad2, fpi,True)  # randomly ..no partcular thought
            #print "dist2,",dist2,fpi
            if abs(dist2 - dist1) < 0.00001:
                break
            else:
                dist0 = np.linalg.norm(p2i - p1)
        edges.append([p1, fpi]);
        edges.append([p2i, fpi]);
        edges.append([p1]);
        edges += [[p2i], [p1, p2i], [fpi]]

        #print r1,np.cross(fpi-p1,p1prime),abs(dist2-dist1),dist1,dist2
        #bz.pe([[edges, 'g', '']])

        if np.cross(fpi - p1, p1prime) >= 0 and abs(dist2 - dist1) < 1e-2:
            return edges, r1, fpi
        else:
            return edges, float('inf'), fpi
    #return r,mp
    @staticmethod
    def pe(edges,n=False,l=False):
        def bbox(pts):
            x = zip(*pts)[0];y=zip(*pts)[1]
            return min(x),min(y),max(x),max(y)

        def vc(eds,ax,ec):
          verts = [];codes=[]
          for ind,ed in enumerate(eds):
              if len(ed[0])==0: continue
              verts.append(ed[0]);codes.append(Path.MOVETO)
              if n and len(ed[0])>0:
                  #print ed
                  ax.text(ed[0][0], ed[0][1], str(ind), style='italic',
                          bbox={'facecolor': ec, 'alpha': 0.5, 'pad': 10})
              if len(ed)==1:
                  #print ed
                  ax.plot(ed[0][0], ed[0][1], 'x-', color=ec)
              elif len(ed)==2:
                  verts.append(ed[1]);codes.append(Path.LINETO)
              elif len(ed)==3:
                  verts.append(ed[1]);codes.append(Path.CURVE3);verts.append(ed[2]);codes.append(Path.CURVE3)
              elif len(ed)==4:
                  verts.append(ed[1]);codes.append(Path.CURVE4);verts.append(ed[2]);codes.append(Path.CURVE4)
                  verts.append(ed[3]);codes.append(Path.CURVE4)
          return verts,codes

        # list of edges
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ec = 'r';label='';pts=[];
        xmin=0.0;ymin=0.0;xmax = 0.0;ymax=0.0
        for eg in edges:
            eg0=eg[0]
            if len(eg)==3:
                ec = eg[1];label = eg[2]
            verts,codes = vc(eg0,ax,ec)
            #print verts,codes
            xmin1,ymin1,xmax1,ymax1 = bbox(verts)
            xmin = min(xmin,xmin1);xmax=max(xmax,xmax1);ymin=min(ymin,ymin1);ymax=max(ymax,ymax1)
            path = Path(verts, codes)
            patch = patches.PathPatch(path, facecolor='none',ec=ec, lw=2)
            ax.add_patch(patch)
        ax.set_xlim(xmin - 30, xmax + 30)
        ax.set_ylim(ymin - 30, ymax + 30)
        plt.show(block=False)
class Font(object):
    def __init__(self,gl="A",fontfile="DevanagariSangamMN.ttf",returnobj=False):
        font= ttfquery.describe.openFont(fontfile)
        glo = ttfquery.glyph.Glyph(gl)
        glo.compile(font,steps=3)
        outlines = glo.outlines
        contours = glo.contours
        self.edges,self.vertices,self.vert_edges,self.cnts = self.edgesFromOutline(outlines, contours, False)
        self.vertAnalyse()
    def edgesFromOutline(self,outlines,contours,returnobj=False):

        codes = []
        edges = []
        offcurve = 0

        points = []
        vert_edges={}
        edge_index = 0
        pt_index = 0
        cnts = []
        for i,contour in enumerate(contours):

           start = True
           prevzero = False
           prevone = False
           edge = []
           cnt = []
           for j,point in enumerate(contour):

               if point[1] == 1: # on-curve
                   edge.append(np.array(point[0]).astype(float))
                   prevzero = False
                   if j!=len(contour)-1:
                       points.append(np.array(point[0]))
                       cnt.append(pt_index)
                       pt_index += 1

                   if prevone:
                      # if len(edge) ==2:
                       #    edge.insert(1,(edge[1]+edge[0])/2.0)

                       #edge = np.array(edge)
                       edges.append(edge)

                       edge_index += 1
                       edge = []
                       edge.append(np.array(point[0]).astype(float))

                       if not pt_index in vert_edges:
                            vert_edges[pt_index] = []
                       vert_edges[pt_index].append(edge_index)


                   start = False
                   prevone = True
               elif point[1] == 0:
                   if prevzero:
                       missing_point = (np.array(point[0]).astype(float) + np.array(prevpoint[0]).astype(float))/2
                       points.append(missing_point)
                       cnt.append(pt_index)

                       edge.append(missing_point)
                       #edge = np.array(edge)
                       edges.append(edge)
                       if not pt_index in vert_edges:
                           vert_edges[pt_index] = []
                       vert_edges[pt_index].append(edge_index)
                       edge_index += 1
                       edge = []
                       vert_edges[pt_index].append(edge_index)

                       edge.append(missing_point)
                       pt_index += 1

                   edge.append(np.array(point[0]).astype(float))
                   prevpoint = point
                   prevzero = True
           cnts.append(cnt)
        if returnobj:
            return edges,outlines,contours,points
        else:
            return edges,points,vert_edges,cnts
    def p(self,e,t):
        return bz.p(self.edges[e],t)
    def pprime(self,e,t):
        return bz.pprime(self.edges[e],t)
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
    def vertAnalyse(self):
        self.concave_vert = []
        self.convex_vert = []
        self.its={}
        #print len(self.vertices),len(self.edges)
        for e,_ in enumerate(self.edges):
            self.its[e]=it([iv(0.0,1.0001)])
        for i,_ in enumerate(self.vertices):
            j = self.prevedge(i)
            cross = np.cross(self.pprime(i,0.0),self.pprime(j,1.0))
            #print i,j,cross

            if cross>0.01:
                #print "convex"
                self.convex_vert.append(i)
            elif cross<-0.01:
                #print "concave"
                self.concave_vert.append(i)
                #self.its[i]=it([iv(-1.0,1.0)])
                #self.its[self.prevedge(i)]=it([iv(0.0,2.0)])
        return self.convex_vert,self.concave_vert
    def footprint(self,e1,t,e2):
        nodes1 = self.edges[e1];nodes2 = self.edges[e2]
        q=bz.footprint(nodes1,t,nodes2) #r,mp,d

        es=[];
        es.append([self.p(e1,t),q[2]]);
        es.append([self.p(e2,q[1]),q[2]])
        #bz.pe([[q[0],'g',''],[self.edges,'r','']])
        #bz.pe([[[self.edges[e1],self.edges[e2]],'r',''],[es,'g',''],[self.edges,'b','']],n=True)
        return q

    def bisector(self,e1,e2):
        pts=[];es=[]
        for q in xrange(11):
            t = q/10.0
            edges,r,mp = self.footprint(e1,t,e2)
            pts.append([mp])
            es += edges
        bz.pe([[pts,'r',''],[[self.edges[e1],self.edges[e2]],'g','']])
        return pts

    def distance(self,e,p,norangecheck=False):
        edge = self.edges[e]
        dist,pi,rn,ep= bz.distance(edge,p)
        if norangecheck:
            return dist,pi,rn
        if len(self.its[e].search(rn))>0:
            #bz.pe([[ep,'g',''],[self.edges,'r','']])
            return dist,pi,rn
        else:
            return float('inf'),pi,rn
    def diffMp(self,e1,e2,e3):
        def f(t):
            _,r,mp = self.footprint(e1,t,e2)
            _,r1,mp1=self.footprint(e1,t,e3)
            return np.linalg.norm(mp-mp1)
        rs=[0.0,0.1]
        t0=self.secant(f,rs)
        _,r, mp = self.footprint(e1, t0, e2)
        _,r1, mp1 = self.footprint(e1, t0, e3)
        print r,r1

        if np.linalg.norm(mp-mp1)<0.01:
            e = [];
            e.append([self.p(e1, t0), mp]);
            e.append([self.p(e2, r), mp]);
            e.append([self.p(e3, r1), mp1]);
            e.append([mp]);
            e.append([mp1])
            #bz.pe([[self.edges, 'r', ''], [e, 'g', '']])

            return t0,r,r1
        else:
            return float('inf'),float('inf'),float('inf')
    def secant(self,f,rs):
        fxi = f(rs[-1])
        fximinus1 = f(rs[-2])
        #print fxi,fximinus1
        rn = rs[-1]
        #return
        for i in xrange(5):
            if (abs(rs[-1]-rs[-2])<0.000000001 or abs(fxi-fximinus1)<0.0000000001 or abs(fxi)<0.000000001):

                break
            ri = rs[-1]
            rn = ri -fxi *(rs[-1]-rs[-2])/(fxi-fximinus1)
            #print rn
            rs.append(rn)
            fximinus1 = fxi
            fxi = f(rn)
        return rn
    def distCheck(self,e1,e2,t0=0.0):
        es=[]
        for e3,_ in enumerate(self.edges):
            if e3 ==e1 or e3==e2:
                continue
            q=self.diffMp(e1,e2,e3)
            if len(sorted(self.its[e1].search(q[0])))>0 and len(sorted(self.its[e2].search(q[1])))>0 and len(sorted(self.its[e3].search(q[2])))>0:
                es.append([e3,q])
        mint=float('inf');mine=[]
        for e in es:
            if e[1][0]<mint and e[1][0]>t0:
                mint = e[1][0];mine=e
        return es,mint,mine
    def nearestEdge(self,p,tracededges):
        mindist = float('inf');mine = -1;minr=float('inf')
        for e,_ in enumerate(self.edges):
            if e in self.tracededges or e in tracededges: continue
            #print "calculating distance for",e
            dist,pi,r = self.distance(e,p)
            if dist<mindist:
                mine=e;mindist=dist;minr=r
        return mine,mindist,minr
    def retract(self,mplast,mpcurr,currdist,e3):
        print "retract:",mplast,mpcurr
        tf = mplast[2];tl=mpcurr[2];e1=mplast[1];e2=mplast[3];radius=mpcurr[-1]
        assert mplast[1]==mpcurr[1] and mplast[3]==mpcurr[3]
        count = 0;mpn=mpcurr;prevcurrdist=float('inf')
        while abs(currdist-radius)>1e-2 and count<50:
            tn = (tf+tl)/2.0
            _,rn,mp = self.footprint(e1,tn,e2)
            currdist,p3,r3 = self.distance(e3,mp,True)
            if abs(currdist-prevcurrdist)<1e-4:
                break
            prevcurrdist = currdist
            #print "currdist",currdist
            p1 = self.p(e1,tn)
            radius = np.linalg.norm(mp-p1)
            mpn = [mp,e1,tn,e2,rn,radius]
            if currdist <radius:
                tf = tn
            else:
                tl=tn
            count += 1
        print abs(currdist-radius)
        return mpn
    def removeEndpoints(self,i):
        #mps = []
        cv=self.convex_vert[i]
        e1 = cv;
        e2 = self.prevedge(cv)
        print e1,e2
        #self.tracededges.append(e1);
        #self.tracededges.append(e2)
        deltat = 0.01;
        t = 0.01
        edges = []
        mplast = [self.vertices[e1], e1, self.its[e1].begin(), e2, self.its[e2].end(), 0.0]
        mps = mplast
        while t < self.its[e1].end():
            _, r, mp = self.footprint(e1, t, e2)
            radius = np.linalg.norm(mp - self.p(e1, t))
            mpcurr = [mp, e1, t, e2, r, radius]
            mine, mindist,_ = self.nearestEdge(mp, [e1,e2])
            mini, distconv = self.nearestConcaveVert(mp)

            if mindist < radius:
                mpcurr = self.retract(mplast, mpcurr, mindist, mine)
                self.its[e1].chop(self.its[e1].begin(), mpcurr[2])
                self.its[e2].chop(mpcurr[4], self.its[e2].end())
                _, r3, mp3 = self.footprint(e1, mpcurr[2], mine)
                radius = np.linalg.norm(mp3 - self.p(mine, r3))
                print np.linalg.norm(mpcurr[0] - mp3)
                mps = [mp3, mine, r3, e2, mpcurr[4], radius]
                break
            if distconv < radius:
                print "dist failure with concavevert", mini
                mpcurr = self.retract(mplast, mpcurr, distconv, self.prevedge(mini))
                _, r3, mp3 = self.footprint(mpcurr[1], mpcurr[2], self.prevedge(mini))
                #_, t0, mp0 = self.footprint(self.prevedge(mini), r3, self.nextedge(mpcurr[1]))
                radius = np.linalg.norm(mp3 - self.p(mpcurr[3], mpcurr[4]))
                mps=[mp3, mpcurr[1],mpcurr[2],self.prevedge(mini), r3, radius]
                print mps
                break

            mplast = mpcurr
            t += deltat

        print mps
        return mps
    def distConcaveVert1(self,i,mp):
        e1 = i;e2=self.prevedge(e1)
        dist1,p1,r1 = self.distance(e1,mp,True)
        dist2, p2, r2 = self.distance(e2, mp, True)
        return max(dist1,dist2)
        #return dist1
    def  distConcaveVert(self,i,mp):
        return np.linalg.norm(self.vertices[i]-mp)
    def nearestConcaveVert(self,mp,e1,e2,e1concaveedgefound,e2concaveedgefound):
        dist0=float('inf');mini = -1
        for i in self.concave_vert:
            if i in self.handledconcaveedges:
                continue
            if e1concaveedgefound or e2concaveedgefound:
                dist = self.distConcaveVert(i,mp)
            else:
                dist = self.distConcaveVert1(i,mp)
            if dist<dist0:
                dist0=dist;mini=i
        return mini,dist0
    def checke1ore2(self,e1,t1,e2): # check if two verts can become e1 or e2
        _,r1,_ = self.footprint(e1,t1,e2)
        if r1==float('inf'):
            return False
        else:
            _,r2,_ = self.footprint(e1,t1+1e-4,e2)
            if r2<r1:
                return True
            else:
                return False
    def proceed(self,mps,edges):
         e1 = mps[1];
         t = mps[2];
         e2 = mps[3];
         deltat = 0.01
         e1concaveedgefound = False;e2concaveedgefound = False;

         mplast = mps
         edgechange = False
         distFailure = False
         currconcaveedge = -1
         convexedgefound=False
         if e1 in self.tracededges or e2 in self.tracededges:
             return

         terminateloop = False
         print e1,e2
         noofsteps=0

         while not convexedgefound and not terminateloop:
             # while t < self.its[e1].end():
             noofsteps += 1

             _,r,mp = self.footprint(e1, t, e2)
             if r == float('inf'):
                 break
             radius = np.linalg.norm(mp - self.p(e1, t))
             mpcurr = [mp, e1, t, e2, r, radius]
             #print e1,e2,radius

             if r < 0.0 and not e2concaveedgefound:
                 #print "r negative??"
                 if e2 in self.convex_vert:
                     convexedgefound=True
                     break
                 if e2 in self.concave_vert:
                     e2concaveedgefound = True
                     currconcaveedge = e2
                     self.handledconcaveedges.append(e2)
                     #continue
                 else:
                     self.tracededges.append(e2)
                     e2p = self.prevedge(e2);
                     e2=e2p
                     print e1,e2
                     #if self.prevedge(e2) in self.concave_vert:
                      #   self.handledconcaveedges.append(self.prevedge(e2)) # to prevent distance check with next edge
                     if e2 in self.tracededges:
                         terminateloop=True
                     #self.tracededges.append(e2)

                     continue
             #print "e1,t,e2", e1, t, e2, r

             if edgechange:
                 mplast=mpcurr;edgechange=False
             noofsteps += 1

             if  noofsteps>5 or edgechange or e1concaveedgefound or e2concaveedgefound :
                 noofsteps=0
                 mine, mindist,rmine = self.nearestEdge(mp, [e1,e2])
                 mini, distconv = self.nearestConcaveVert(mp,e1,e2,e1concaveedgefound,e2concaveedgefound)
                 concavedist=False
                 if mindist < radius:
                     print "distance failure with ", mine, mindist,rmine
                     if self.nextedge(mine) in self.concave_vert:
                         concavedist = True
                     #print "retracting from",mpcurr," between ",mplast
                     if not edgechange:
                        mpcurr = self.retract(mplast, mpcurr, mindist, mine)
                     #print "retracted to ",mpcurr
                     _, r3, mp3 = self.footprint(mpcurr[3], mpcurr[4], mine)

                     if r3<0:
                         print "breaking ",mpcurr[1],mpcurr[2],mine
                     distFailure = True
                     if e1concaveedgefound:
                        _,t,mp4 = self.footprint(mine,r3,self.nextedge(mpcurr[1]))
                        #e1 = self.nextedge(mpcurr[1]);e2=mine;e1concaveedgefound=False
                        #print "continuing with current branch", e1,r3, e2
                        if concavedist:
                            if not e2concaveedgefound:
                                e1 = mine;t=r3
                                print "continue with ",e1,r3,e2
                            dist0,pi,t0 = self.distance(self.nextedge(mpcurr[1]),mpcurr[0],True)

                            #_,t0, mp5 = self.footprint(mpcurr[1], mpcurr[2], self.nextedge(mine))
                            _, t1, mp6 = self.footprint(self.nextedge(mpcurr[1]), t0, self.prevedge(mine))

                            #radius = np.linalg.norm(mp5 - self.p(mpcurr[1], mpcurr[2]))
                            radius = dist0
                            self.mps.append([mp3, self.nextedge(mpcurr[1]), t0,self.prevedge(mine) , t1, radius])
                            e1concaveedgefound = False
                            print "branching to",self.nextedge(mpcurr[1]),  t0,self.prevedge(mine)
                            continue
                        else:
                            _, t0, mp5 = self.footprint(mine, r3, self.nextedge(mpcurr[1]))
                            radius = np.linalg.norm(mp5 - self.p(mine, r3))
                            self.mps.append([mp3, self.nextedge(mpcurr[1]), t0,mine, r3, radius])
                            e1concaveedgefound = False
                            print "from e1 side branching to", self.nextedge(mpcurr[1]),mine,e2concaveedgefound

                            if not e2concaveedgefound: break

                     if e2concaveedgefound:
                        if concavedist:
                            if True or not e1concaveedgefound:
                                e2 = self.prevedge(mine);
                                print "continue with ",e1,mpcurr[2],e2
                            dist0,pi,t0 = self.distance(self.nextedge(mine),mpcurr[0],True)

                            #_,t0, mp5 = self.footprint(mpcurr[1], mpcurr[2], self.nextedge(mine))
                            _, t1, mp6 = self.footprint(self.nextedge(mine), t0, self.prevedge(mpcurr[3]))

                            #radius = np.linalg.norm(mp5 - self.p(mpcurr[1], mpcurr[2]))
                            radius = dist0
                            self.mps.append([mp3, self.nextedge(mine), t0, self.prevedge(mpcurr[3]), t1, radius])
                            e2concaveedgefound = False
                            print "branching to", self.nextedge(mine), t0,self.prevedge(mpcurr[3])
                            continue
                        else:
                            _, t0, mp5 = self.footprint(mine, r3, self.prevedge(mpcurr[3]))
                            radius = np.linalg.norm(mp5 - self.p(mine, r3))
                            self.mps.append([mp3, mine, r3, self.prevedge(mpcurr[3]), t0, radius])
                            e2concaveedgefound = False
                            print "from e2 side branching to", mine, self.prevedge(mpcurr[3])

                            break


                     continue
                     # break
                 if distconv < radius:
                     print "dist failure with concavevert", mini,distconv,radius,e1concaveedgefound,e2concaveedgefound
                     if not edgechange:
                        while distconv>mplast[-1] and len(edges)>1:
                            mplast = edges.pop(-1)
                            distconv = self.distance(mini,mplast[0],True)
                        mpcurr = self.retract(mplast, mpcurr, distconv, self.prevedge(mini))
                     if e1concaveedgefound and not e2concaveedgefound:
                        _, r3, mp3 = self.footprint(mpcurr[1], mpcurr[2], self.prevedge(mini))

                        _,t0,mp0 = self.footprint(self.prevedge(mini),r3,self.nextedge(mpcurr[1]))
                        if t0==float('inf'): #our assumption about it wrong..
                            #self.mps.append([mp0,mine,r3,mpcurr[3],mpcurr[4],mpcurr[-1]])
                            e2=mine;
                            continue
                        else:
                            radius = np.linalg.norm(mp0-self.p(self.nextedge(mpcurr[1]),t0))
                            self.mps.append([mp3,self.nextedge(mpcurr[1]),t0,self.prevedge(mini),r3,radius])
                            e1 = mini;_,t,mp=self.footprint(mpcurr[3],mpcurr[4],e1)
                            e2=mpcurr[3]
                            print "current branch with",e1,t,e2," and another branch with",self.nextedge(mpcurr[1]),t0,self.prevedge(mini)
                     elif e2concaveedgefound:
                         _, r3, mp3 = self.footprint(mpcurr[1], mpcurr[2], mini)
                         _, t0, mp0 = self.footprint(mini, r3, self.prevedge(mpcurr[3]))
                         radius = np.linalg.norm(mp0 - self.p(self.prevedge(mpcurr[3]), t0))
                         if not e1concaveedgefound:
                            e2=self.prevedge(mini)
                            print "continue current branch with",e1,e2

                         self.mps.append([mp3,mini,r3,self.prevedge(mpcurr[3]),t0,radius])
                         print "new branch with ",mini,self.prevedge(mpcurr[3])
                     else:
                         e2=self.prevedge(mini)
                         _,r3,mp3 = self.footprint(mpcurr[3],mpcurr[4],mini)
                         radius = np.linalg.norm(mp3-self.p(mini,r3))
                         self.mps.append([mp3,mini,r3,mpcurr[3],mpcurr[4],radius])
                         continue


                     e2concaveedgefound=False

                     #print mpcurr

                     distFailure = True
                     # break
                 mplast = mpcurr
             #if not (mpcurr[1] in self.convex_vert and mpcurr[3] == self.prevedge(mpcurr[1])):
             edges.append(mpcurr)

             #mplast = mpcurr
             t += deltat
             edgechange = False
             #print "reached here", e1, t, self.its[e1].end()
             if t >= self.its[e1].end():
                 self.tracededges.append(e1)

                 e1p = self.nextedge(e1);

                 if e1p in self.convex_vert:
                     convexedgefound=True
                     break
                 if e1p in self.concave_vert:
                     e1concaveedgefound = True
                     currconcaveedge = e1p
                     self.handledconcaveedges.append(e1p)
                 else:

                     t = self.its[e1].begin()
                     edgechange = True
                     e1 = e1p
                     print e1,e2
                     if e1 in self.tracededges:
                         terminateloop=True
                     #self.tracededges.append(e1)


    def iterate(self,edges):
        #take a convex verex- take its previous edge and this edge
        #do distcheck find e1,e3
        #increment e1 and find footprint
        #till you reach a concave_vertex ( either e1 or e3)
        # when you reach a concave vertex
        # create a branch and then proceed

        self.tracededges=[]
        self.handledconcaveedges=[]
        convc = self.concave_vert[0]
        convexvert=[]
        for convc in self.concave_vert:
            e0 = self.nextedge(convc)
            while not e0 in self.convex_vert and e0!=convc:
                e0 = self.nextedge(e0)
            try:
                ind = self.convex_vert.index(e0)
                if ind>=0:
                    convexvert.append(ind)
            except: pass
        for ind in self.convex_vert:
            # mps= self.removeEndpoints(ind)
             if ind in self.tracededges or self.prevedge(ind) in self.tracededges: continue
             mps = [self.vertices[ind],ind,self.its[ind].begin(),self.prevedge(ind),1.0,0.0]

             self.mps=[mps]
             print mps
             while len(self.mps)>0:
                    mps = self.mps.pop(0)
                    edge=[];edges.append(edge)
                    self.proceed(mps,edge)
            #bz.pe([[edges, 'r', ''], [self.edges, 'g', '']])

        return





        qs= self.distCheck(e1,e2,0.0)
        #mint=qs[0][1][0];minq=qs[0]
        #if len(qs)>1:
         #   for i,_ in enumerate(qs):
          ##         mint=qs[i][1][0];minq=qs[i]
        print "would iterate ",qs[2][0],e2
        tin = qs[2][1][2];r=float('inf')
        while r>0.0 and tin<1.0:
            _,r,mp = self.footprint(qs[2][0],tin,e2)
            edges.append([mp])
            tin += deltat
        if tin>1.0: # t ended
            qss = self.distCheck(qs[2][0],e2,1.0)

            while tin<qss[1]:
                _, r, mp = self.footprint(qs[2][0], tin, e2)
                edges.append([mp])
                tin += deltat
            preve1 = e1;preve2 = e2;prevr = r
            e1 = self.nextedge(qs[2][0]);
            e2 = qss[2][0]
            print "would iterate ", e1, e2

            r = qss[2][1][2]
            _,tin,mp = self.footprint(e2,r,e1)
            edges.append([mp]);tin += deltat
            while tin<min(self.its[e1].end(),1.0):
                _,r,mp = self.footprint(e1,tin,e2)
                edges.append([mp])
                tin += deltat
            if self.nextedge(e2) in self.concave_vert:
                e1=self.nextedge(e2);
                print "would now iterate",e1,preve2
                _,tin,mp = self.footprint(preve2,prevr,e1)
                tin = tin +deltat
                while tin<min(1.0,self.its[e1].end()):
                    _,r,mp = self.footprint(e1,tin,preve2)
                    edges.append([mp])
                    tin += deltat




        #qss=self.distCheck(minq[0],e2)
        #print edges






        return
    def pe1(self,edges):
        x = []
        for ed in edges:
            print ed
            if len(ed)==0 :
                continue
            ed = zip(*ed)[0]
            eds=[]
            for pt in ed:
                eds.append([pt])
            x.append([eds, 'r', ''])
        x.append([self.edges, 'g', ''])
        bz.pe(x)


class strokes(object):
    def f(self,e1,t,e2):
        r,mp,d=bz.footprint()
    def iterate(self):
        #take pair of convex edges
        # determine a point at which their distance with third edge becomes constant
        return


def test():
    f = Font()

    #bz.pe([f.edges,'r',''])
    return






