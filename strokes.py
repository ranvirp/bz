import numpy as np
import operator as op
import math
import bezier
import ttfquery,ttfquery.describe,ttfquery.glyphquery,ttfquery.glyph
from intervaltree import IntervalTree as it
from intervaltree import Interval as iv

from matplotlib import pyplot as plt
from matplotlib.path import Path
from fitCurves import fitCurve
import matplotlib.patches as patches
from matplotlib.collections import PathCollection
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
                #print abs(fxi-fximinus1),abs(fxi)

                break
            ri = rs[-1]
            rn = ri -fxi *(rs[-1]-rs[-2])/(fxi-fximinus1)
            #print "rn=",rn
            rs.append(rn)
            fximinus1 = fxi
            fxi = f(rn)
        return rn
    @staticmethod
    def distance(quad,p0,rin=0.4,debug=False):
        def f(t):
            #print t,bz.pprime(quad,t)
            vect1 = bz.pprime(quad,t)
            vect2 = p0-bz.p(quad,t)
            #norm1 = np.linalg.norm(vect1)
            #norm2=np.linalg.norm(vect2)
            #return np.dot(vect1,vect2)/(norm1*norm2)
            return np.dot(vect1,vect2)
            #return np.cross(bz.normal)
        rs = [];rs.append(rin);rs.append(rin+0.2)
        rn = bz.secant(f,rs)
        #print rn,f(rn)
        p1 = bz.p(quad,rn);p1prime=bz.pprime(quad,rn)
        if True or debug:
            ep=[];ep.append([p0]);ep.append([p0,p1]);ep.append([p1,p1+p1prime])
        #print "f(rn)",abs(f(rn)),abs(np.cross(p0-p1,p1prime)),abs(np.cross(p0-p1,p1prime)) <0.0001, np.cross(p0-p1,p1prime)>0
        if (abs(np.cross(p0-p1,p1prime)) <0.0001 or np.cross(p0-p1,p1prime)>0 )and abs(f(rn))<0.01: # point is on right side of the side
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
        dist0, p2i, r1,_ = bz.distance(quad2, p1)  # randomly ..no partcular thought
        rprev=r1
        if dist0 ==float('inf'):
            #print "here"
            p2i = bz.p(quad2,0.5);dist0 = np.linalg.norm(p2i-p1);rprev=0.5

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
            dist2, p2i, r1,_ = bz.distance(quad2, fpi,rprev)  # randomly ..no partcular thought
            rprev=r1
            '''
            if dist2==float('inf'):
                count2=0
                while dist2==float('inf') and count2<30:
                    dist1=dist1/2.0;count2+=1
                    fpi=p1+dist1*n1
                    dist2, p2i, r1, _ = bz.distance(quad2, fpi, rprev)  # randomly ..no partcular thought
                    print "dist2,dist1",dist2,dist1
                    rprev=r1
            '''

            #print "dist2,",dist2,fpi,dist1,abs(dist2-dist1)
            if abs(dist2 - dist1) < 0.00001:
                break
            else:
                dist0 = np.linalg.norm(p2i - p1)
        edges.append([p1, fpi]);
        edges.append([p2i, fpi]);
        edges.append([p1]);
        edges += [[p2i], [p1, p2i], [fpi]]

        #print r1,np.cross(fpi-p1,p1prime),abs(dist2-dist1),dist1,dist2
        #bz.pe([[edges, 'g', ''],[[quad1,quad2],'b','']])

        if np.cross(fpi - p1, p1prime) >= 0 and abs(dist2 - dist1) < 1e-2:
            #print "return here"
            return edges, r1, fpi
        else:
            return edges, float('inf'), fpi
    #return r,mp
    @staticmethod
    def pe(edges,n=False,l=False):
        def onpick(event):
            thispath = event.artist
            paths =thispath.get_paths()
            #xdata,ydata=event.get_data()
            for path in paths:
                if path.contains_point(event.mouseevent.xdata,event.mouseevent.ydata):
                    print path
            #path=thispath.get_path()
            #print path
            #for vertices,code in path.iter_segments():
             #   print len(vertices)
            #print('onpick points:', points)

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
        patches1=[]
        for eg in edges:
            ec = 'r';
            label = '';
            eg0=eg[0]
            if len(eg)==3:
                ec = eg[1];label = eg[2]
            verts,codes = vc(eg0,ax,ec)
            #print verts,codes
            xmin1,ymin1,xmax1,ymax1 = bbox(verts)
            xmin = min(xmin,xmin1);xmax=max(xmax,xmax1);ymin=min(ymin,ymin1);ymax=max(ymax,ymax1)
            path = Path(verts, codes)
            patches1=[path]
            colors = 100 * np.random.rand(len(patches1))
            p = PathCollection(patches1,facecolors=["none"], edgecolors=[ec])
            p.set_picker(True)
            #p.set_array(np.array(colors))
            ax.add_collection(p)

            #patch = patches.PathPatch(path, facecolor='none',ec=ec, lw=2)
            #patch.picker=True
            #path.set_picker(True)
            #patches1 += path

            #q=ax.add_patch(patch)

            #q.picker=True
        ax.set_xlim(xmin - 30, xmax + 30)
        ax.set_ylim(ymin - 30, ymax + 30)
        cid=fig.canvas.mpl_connect('pick_event', onpick)
        #print cid
        plt.show(block=False)
        #plt.show()
    def penew(edges,n=False,l=False):
        def onpick(event):
            thispath = event.artist
            print thispath.get_path().vertices
            #print('onpick points:', points)

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
            #patch.picker=True
            #patch.set_picker(True)

            q=ax.add_patch(patch)

            #q.picker=True
        ax.set_xlim(xmin - 30, xmax + 30)
        ax.set_ylim(ymin - 30, ymax + 30)
        cid=fig.canvas.mpl_connect('pick_event', onpick)
        #print cid
        plt.show(block=False)
        #plt.show()
    def merge(self,eds):
        #if end of one edge and end of other edge ..e1 is same
        def mergeeds(ed1,ed2):
            if ed1[-1][1]==ed2[-1][1]:
                ed1=ed1+ed2;ed2=[]
            elif ed1[0][1]==ed2[-1][1]:
                ed1=reversed(ed1)+ed2;ed2=[]
            elif ed2[-1][1]==ed1[-1][1]:
                ed2=ed1+ed2;ed1=[]
            elif ed2[0][1]==ed1[-1][1]:
                ed2=reversed(ed2)+ed1;ed1=[]
            return ed1,ed2

        for ind1,ed1 in enumerate(eds):
            for ind2,ed2 in enumerate(eds):
                if ind2<=ind1: continue




        return
class Font(object):
    def __init__(self,gl="A",fontfile="DevanagariSangamMN.ttf",returnobj=False):
        font= ttfquery.describe.openFont(fontfile)
        glo = ttfquery.glyph.Glyph(gl)
        glo.compile(font,steps=3)
        outlines = glo.outlines
        contours = glo.contours
        self.edges,self.vertices,self.vert_edges,self.cnts = self.edgesFromOutline(outlines, contours, False)
        self.vertAnalyse()
        self.debug3=False;self.debug=False
        self.pairs={}
        self.tracededges=[]
        self.cvs={}
        self.maxradius = 0.0
        for cv in self.concave_vert:
            self.cvs[cv]={}
            self.updateo(cv)
    def show(self):
        bz.pe([[self.edges]],n=True)

    def apairs(self,e1,e2):
        if not str(e1)+"-"+str(e2) in self.pairs:
            self.pairs[str(e1)+"-"+str(e2)] = True
            self.pairs[str(e2) + "-" + str(e1)] = True

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
    def normal(self,e,t,unit=False):
        return bz.normal(self.edges[e],t,unit)
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
            #cross = np.cross(self.pprime(i,0.0),self.pprime(j,1.0))
            cross=np.cross(self.normal(i,0.0,True),self.normal(j,1.0,True))
            #print i,j,cross

            if cross>0.2:
                #print "convex"
                self.convex_vert.append(i)
            elif cross<-0.2:
                #print "concave"
                self.concave_vert.append(i)
                #self.its[i]=it([iv(-1.0,1.0)])
                #self.its[self.prevedge(i)]=it([iv(0.0,2.0)])
            '''
            Ri = self.radiusCurvature(i,0.0)
            Rj = self.radiusCurvature(j,1.0)
            print "curvaturedifference",i,j,Ri,Rj,Ri-Rj
            '''
        return self.convex_vert,self.concave_vert
    def footprint(self,e1,t,e2):
        nodes1 = self.edges[e1];nodes2 = self.edges[e2]
        q=bz.footprint(nodes1,t,nodes2) #r,mp,d

        es=[];
        #es.append([self.p(e1,t),q[2]]);
        #es.append([self.p(e2,q[1]),q[2]])
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
        bz.pe([[pts,'r',''],[[self.edges[e1],self.edges[e2]],'g',''],[self.edges,'b','']])
        return pts

    def distance(self,e,p,norangecheck=False,rin=0.4):
        edge = self.edges[e]
        dist,pi,rn,ep= bz.distance(edge,p,rin)
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
            #print e,dist,pi,r
            if dist<mindist:
                mine=e;mindist=dist;minr=r
        return mine,mindist,minr

    def nearestEdge1(self,e1, t, tracededges):
        mindist = float('inf');
        mine = -1;
        minr = float('inf')
        for e, _ in enumerate(self.edges):
            if e in self.tracededges or e in tracededges: continue
            # print "calculating distance for",e
            _,r,mp = self.footprint(e1, t,e)
            if r!=float('inf'):
                dist = np.linalg.norm(mp-self.p(e1,t))
            # print e,dist,pi,r
                if dist < mindist:
                    mine = e;
                    mindist = dist;
                    minr = r
        return mine, mindist, minr
    def resultMp(self,cv,e2,radius1,e1=True):
        e10 = cv
        if e1:
            e10 = self.nextedge(cv)
            # else:
            #   e10 = self.prevedge(cv)

        def n(t):
            theta = self.cvs[e10]['theta1'] * (1 - t) + self.cvs[e10]['theta2'] * t
            return np.array([np.cos(theta), np.sin(theta)])

        theta0 = (self.cvs[e10]['theta1'] - self.cvs[e10]['theta2']) / 2.0
        radius = radius1 / abs(np.cos(theta0))
        x=[]
        res= self.vertices[e10]+ n(0.5)*radius
        x.append([res])
        #bz.pe([[x],[self.edges,'g','']])
        return res

    def printNormals(self,cv):
        ep = []
        def n(t):
            theta = self.cvs[cv]['theta1'] * (1 - t) + self.cvs[cv]['theta2'] * t
            return np.array([np.cos(theta), np.sin(theta)])

        for q in xrange(11):
            t = q / 10.0
            ep.append([self.vertices[cv] + 20 * n(t), self.vertices[cv]])
            # theta =self.theta1*(1-t) + self.theta2*t
        bz.pe([[ep],[self.edges,'g','']])

    def findCounterPart1(self,cv,e2,radius1,e1=True):
        e10=cv
        if e1:
            e10 = self.nextedge(cv)
        #else:
         #   e10 = self.prevedge(cv)

        def n(t):
            theta = self.cvs[e10]['theta1']*(1-t)+self.cvs[e10]['theta2']*t
            return np.array([np.cos(theta),np.sin(theta)])
        theta0 = (self.cvs[e10]['theta1'] - self.cvs[e10]['theta2']) / 2.0
        #radius = radius1 / abs(np.cos(theta0))
        radius=radius1
        breakloop = False
        mindist = float('inf');
        mindc = float('inf')
        rmine=-1
        minc=float('inf')
        mine=-1
        edgedist=[];cvdist=[]

        for q in xrange(1):
            if breakloop: break
            if e1:
                t = 1.0
            else:
                t=0.0
            mp=self.vertices[e10] + n(t)*radius
            for e,_ in enumerate(self.edges):
                if e in self.tracededges:
                    continue
                if e==e2 or e==e10  or e==self.nextedge(e10) or e==self.prevedge(e10): continue
                dist,fp,r = self.distance(e,mp)
                if dist==float('inf'): continue
                if dist!=float('inf'):
                    #print "edge",e,dist
                    pass
                edgedist.append([dist,e,r])
            for i in self.concave_vert:
                if i==e10 or i==e2 : continue
                if i in self.handledconcaveedges: continue

                dist1 = np.linalg.norm(mp-self.vertices[i])
                #print "concave vert",i,dist1
                cvdist.append([dist1,i])

                if dist1<mindc:
                    mindc=dist1;minc = i
        #print mine,mindist,minc,mindc
        #return mine,mindist,rmine,minc,mindc
        edgefinal=[];cvfinal=[];cvfound=False
        edgedist.sort(key=lambda x:x[0])
        for ed in edgedist:
            if e1:
                if self.checkPair(e10,0.0,ed[1]):
                    edgefinal = ed
                    break
            else:
                if self.checkPair(ed[1],ed[2],self.prevedge(e10)):
                    edgefinal = ed
                    break
        cvdist.sort(key=lambda x:x[0])
        for cved in cvdist:
            if e1:
                if self.checkPair(e10,0.0,self.prevedge(ed[1])):
                    cvfinal = cved
                    break
            else:
                if self.checkPair(cved[1],0.0,self.prevedge(e10)):
                    cvfinal = cved
                    break
        #print "cp,",cvfinal,edgefinal,self.maxradius
        if len(cvfinal)>0 and len(edgefinal)>0:
            if abs(cvfinal[0]-edgefinal[0])<30.0 and cvfinal[0]<self.maxradius*2.0:
                cvfound=True
                return cvfound,cvfinal
            elif edgefinal[0]<cvfinal[0] and edgefinal[0]<self.maxradius:
                cvfound=False
                return cvfound,edgefinal
            elif cvfinal[0]<edgefinal[0] and cvfinal[0]<self.maxradius*2.0:
                cvfound = True
                return cvfound,cvfinal
            else:
                return True,[]
        elif len(edgefinal)>0 and len(cvfinal)==0 and edgefinal[0]<self.maxradius:
            return False,edgefinal
        elif len(cvfinal)>0 and len(edgefinal)==0 and cvfinal[0]<self.maxradius*2.0:
            return True,cvfinal
        else:
            return False,[]


        #return edgedist,cvdist

    def findCounterPart(self,cv,e2,radius1,e1=True):
        e10=cv
        if e1:
            e10 = self.nextedge(cv)
        #else:
         #   e10 = self.prevedge(cv)

        def n(t):
            theta = self.cvs[e10]['theta1']*(1-t)+self.cvs[e10]['theta2']*t
            return np.array([np.cos(theta),np.sin(theta)])
        theta0 = (self.cvs[e10]['theta1'] - self.cvs[e10]['theta2']) / 2.0
        radius = radius1 / abs(np.cos(theta0))
        #radius=radius1
        breakloop = False
        mindist = float('inf');
        mindc = float('inf')
        rmine=-1
        minc=float('inf')
        mine=-1

        for q in xrange(1):
            if breakloop: break
            if e1:
                t = 1.0
            else:
                t=0.0
            mp=self.vertices[e10] + n(t)*radius
            for e,_ in enumerate(self.edges):
                if e in self.tracededges:
                    continue
                if e==e2 or e==e10  or e==self.nextedge(e10) or e==self.prevedge(e10): continue
                dist,fp,r = self.distance(e,mp)
                if dist!=float('inf'):
                    #print "edge",e,dist
                    pass
                if dist<mindist:
                    mine=e;mindist=dist;rmine=r
            for i in self.concave_vert:
                if i==e10 or i==e2 : continue
                if i in self.handledconcaveedges: continue

                dist1 = np.linalg.norm(mp-self.vertices[i])
                #print "concave vert",i,dist1
                if dist1<mindc:
                    mindc=dist1;minc = i
        #print mine,mindist,minc,mindc
        return mine,mindist,rmine,minc,mindc

        return
        res=[]
        t=0.0;mindist=float('inf');cp=-1;minr=float('inf')
        if not e1: e1=self.prevedge(cv);e2=cv;t=1.0
        else: e1=self.nextedge(cv);e2=cv
        n=self.normal(e1,t,True)
        mp = self.p(e1,t) + n*radius
        for e,_ in enumerate(self.edges):
            if e==e1 or e==e2: continue
            dist,fp,r = self.distance(e,mp)
            if dist<mindist and dist<radius:
                cp=e;mindist=dist;minr=r

        print "counterpart of",cv,"is",cp
        return cp,minr,mindist


    def radiusCurvature(self, e, t):
        # R=((x^('2)+y^('2))^(3/2))/(|x^'y^('')-y^'x^('')|)
        # R = norm(pprime)^3/norm(dot(pprime,normal of pdoubleprime)
        quad=self.edges[e]
        if len(quad)==2: return float('inf')
        pprime = self.pprime(e, t)
        pdprime = 2 * (quad[2] + quad[0] - 2 * quad[1])
        pdprimen = np.array([pdprime[1], -1 * pdprime[0]])
        if np.linalg.norm(pdprime == 0):
            R = float('inf')
        else:
            R = math.pow(np.linalg.norm(pprime), 3) / np.linalg.norm(np.dot(pprime, pdprimen))
        return R
    def checkPair(self,e1,t,e2):
        _,r,mp = self.footprint(e1,t,e2)
        _,r1,mp1= self.footprint(e1,t+0.001,e2)
        if r1<r and r!=float('inf') and r1!=float('inf'):
            return True
        else:
            #print "check failed",e1,e2
            return False

    def retract(self,mplast,mpcurr,currdist,e3):
        if self.debug3: print "retract:",mplast,mpcurr
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
        #print abs(currdist-radius)
        return mpn

    def retractConcave(self, mplast, mpcurr, currdist, e3):
        if self.debug3: print "retract concave:", mplast, mpcurr
        tf = mplast[2];
        tl = mpcurr[2];
        e1 = mplast[1];
        e2 = mplast[3];
        radius = mpcurr[-1]
        assert mplast[1] == mpcurr[1] and mplast[3] == mpcurr[3]
        count = 0;
        mpn = mpcurr;
        prevcurrdist = float('inf')
        while abs(currdist - radius) > 1e-2 and count < 50:
            tn = (tf + tl) / 2.0
            _, rn, mp = self.footprint(e1, tn, e2)
            currdist = self.distConcaveVert(e3, mp)
            if abs(currdist - prevcurrdist) < 1e-4:
                break
            prevcurrdist = currdist
            # print "currdist",currdist
            p1 = self.p(e1, tn)
            radius = np.linalg.norm(mp - p1)
            mpn = [mp, e1, tn, e2, rn, radius]
            if currdist < radius:
                tf = tn
            else:
                tl = tn
            count += 1
        #print abs(currdist - radius)
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
    def distConcaveVert(self,i,mp):
        e1 = i;e2=self.prevedge(e1)
        dist1,p1,r1 = self.distance(e1,mp,True)
        dist2, p2, r2 = self.distance(e2, mp, True)
        return max(dist1,dist2)
        #return dist1
    def  distConcaveVert1(self,i,mp):
        return np.linalg.norm(self.vertices[i]-mp)
    def nearestConcaveVert(self,mp,e1,e2,e1concaveedgefound,e2concaveedgefound,skipconcaves):
        dist0=float('inf');mini = -1
        for i in self.concave_vert:
            if i in self.handledconcaveedges:
                continue
            if i in skipconcaves:
                continue
            if i==self.nextedge(e1) or i ==e2:
                continue
            if e1concaveedgefound or e2concaveedgefound:
                dist = self.distConcaveVert1(i,mp)
            else:
                dist = self.distConcaveVert(i,mp)
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

    def updateo(self,vert_ind):
        p2prime = self.pprime(vert_ind, 0.0)
        p1prime = self.pprime(self.prevedge(vert_ind), 1.0)
        n1 = np.array([p1prime[1], -1 * p1prime[0]])
        n2 = np.array([p2prime[1], -1 * p2prime[0]])
        theta2  = np.arctan2(n2[1], n2[0])
        theta1 = np.arctan2(n1[1], n1[0])
        thetaflag = False
        if theta2 < 0 and theta1 > 0:
            theta2 += 2 * np.pi;
            thetaflag = True
        if theta2 < 0 and theta1 < 0 and abs(theta2) > abs(theta1):
            theta2 += 2 * np.pi;
            thetaflag = True

        self.cvs[vert_ind]['theta2'] = theta2
        self.cvs[vert_ind]['n1'] = n1;
        self.cvs[vert_ind]['n2'] = n2
        self.cvs[vert_ind]['theta1'] = theta1
        self.cvs[vert_ind]['thetaflag'] = thetaflag

    def extendCv(self,e1prev,e2prev,t0,mpfin,radius,edges):
        mpcurr = edges[-1]
        _, _, tfin = self.distance(e1prev, mpfin, True)
        if tfin==float('inf'): return edges[-1]
        deltat1 = 0.01;
        count1 = 0
        while t0 <= tfin:
            _, r1, mp1 = self.footprint(e1prev, t0, e2prev)
            #print e1prev, t0, e2prev, r1
            mpcurr = [mp1, e1prev, t0, e2prev, r1, radius]
            t0 += deltat1;
            count1 += 1
            edges.append(mpcurr)
        return mpcurr

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
         if str(e1)+"-"+str(e2) in self.pairs:
             return
         if e1 in self.tracededges:
             print "sorry",e1," traced";return
         if e2 in self.tracededges:
             print "sorry",e2,"traced";return

         terminateloop = False
         e1start = e1;e2start=e2;tin=t;
         rin=-float('inf')
         print e1,e2
         self.apairs(e1,e2)
         noofsteps=0

         while not convexedgefound and not terminateloop:
             # while t < self.its[e1].end():
             e1plus = self.nextedge(e1);
             e2minus = self.prevedge(e2)
             if e1plus in self.concave_vert:
                 e1concaveedgefound = True
                 #self.handledconcaveedges.append(e1plus)
             else:
                    e1concaveedgefound = False
             #if e2 in self.concave_vert:
              #   e2concaveedgefound = True

             noofsteps += 1

             _,r,mp = self.footprint(e1, t, e2)
             #print e1,t,e2,r
             if t>3.0:
                 return
             if r==float('inf'):
                 return
                 #_,t,mp = self.footprint(e2,1.0,e1)
                 #print "inf"
                 t+=deltat;
                 #e2=self.prevedge(e2)
                 continue
             if t>tin:
                self.its[e1].chop(tin,t)
             if rin>r:
                 self.its[e2].chop(r,rin)
             if rin==-float('inf'): rin = r
             radius = max(np.linalg.norm(mp - self.p(e1, t)),np.linalg.norm(mp-self.p(e2,r)))
             mpcurr = [mp, e1, t, e2, r, radius]
             #print e1,t,e2,radius

             if r < self.its[e2].begin():# and not e2concaveedgefound:
                 if e2concaveedgefound and r<-1.5:
                     e2=self.prevedge(e2)
                     e2concaveedgefound=False
                     continue

                 #print "r negative??"
                 if e2 in self.convex_vert:
                     convexedgefound=True
                     break
                 if e2 in self.concave_vert:
                     e2concaveedgefound = True
                     currconcaveedge = e2
                     #self.handledconcaveedges.append(e2)
                     #continue
                 else:
                     if self.its[e2].end() - self.its[e2].begin() < 0.1:
                         self.tracededges.append(e2)

                     e2p = self.prevedge(e2);
                     e2=e2p
                     print e1,e2
                     self.apairs(e1, e2)

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
             #else:
              #   e2concaveedgefound = False

             if  self.prevedge(e1)==e2 or e1concaveedgefound or e2concaveedgefound:# or noofsteps>5 or edgechange:
                 noofsteps=0
                 skipedges=[e1,e2]
                 if not e1plus  in self.convex_vert: skipedges.append(e1plus)
                 if not e2minus in self.convex_vert: skipedges.append(e2minus)
                 mine, mindist,rmine = self.nearestEdge(mp, skipedges)
                 print mine,midist,rmine,mp

                 skipconcaves=[]
                 # we are not going to check distances with concaveedges of next branch
                 if not e2concaveedgefound:
                     e0 = e2minus
                     while not e0 in self.convex_vert and e0!=e2:
                         if e0 in self.concave_vert:
                             skipconcaves.append(e0)
                         e0=self.prevedge(e0)
                 if not e1concaveedgefound:
                     e0 = e1plus
                     while not e0 in self.convex_vert and e0 != e1:
                         if e0 in self.concave_vert:
                             skipconcaves.append(e0)
                         e0= self.nextedge(e0)

                 mini, distconv = self.nearestConcaveVert(mp,e1,e2,e1concaveedgefound,e2concaveedgefound,skipconcaves)
                 concavedist=False;
                 #print "e1,t,mindist,distconv,radius",e1,t,mindist,distconv,radius
                 if e1concaveedgefound:
                     if distconv<np.linalg.norm(mpcurr[0]-self.vertices[e1plus]): concavedist=True
                 elif e2concaveedgefound:
                     if distconv < np.linalg.norm(mpcurr[0] - self.vertices[e2]): concavedist = True
                 else:
                     if  distconv<radius: concavedist = True


                 #print e1,t,concavedist,e1concaveedgefound,e2concaveedgefound


                 if mindist < radius and not concavedist:
                     print "distance failure with ", mine,mindist,radius
                     if self.debug3:
                         bz.pe([[[[self.p(mine,rmine)]],'g',''],[[[mpcurr[0]]],'r',''],[[[self.p(e1,t)]],'b',''],[self.edges]])
                     '''
                     if mine in self.convex_vert:
                         if mine ==e1plus:
                            self.tracededges.append(e1)
                         elif mine==e2minus:
                             self.tracededges.append(e2)
                     '''
                     if not edgechange:
                        mpcurr = self.retract(mplast, mpcurr, mindist, mine)
                     if e1concaveedgefound and e2concaveedgefound: #case #2A
                         _,r1,mp1 = self.footprint(mine,rmine,e1plus)
                         _,r2,mp2 = self.footprint(mine,rmine,e2minus)
                         self.mps.append([mp2,mine,rmine,e2minus,r2,mpcurr[-1]])
                         #self.mps.append([mp2, mine, rmine, e2, mpcurr[4], mpcurr[-1]])
                         #print "branch to",mine,rmine,e2
                         print "branch to",mine,rmine,e2minus
                         e1=e1plus;t=r1;e2=mine;tin=t;rin=-float('inf')
                         mplast[1] = e1;
                         mplast[2] = t;
                         mplast[3] = e2;

                         print "continue with ",e1,t,e2
                         #self.mps.append([mp1, e1plus, r1, mine, rmine, mpcurr[-1]])

                         e1concaveedgefound =False;e2concaveedgefound = False
                         continue
                     elif e2concaveedgefound and not e1concaveedgefound:
                         e1=mine;t=rmine;e2=e2minus;tin=t;rin=-float('inf')
                         mplast[1] = e1;
                         mplast[2] = t;
                         mplast[3] = e2

                         print "continue with ",e1,t,e2
                         e2concaveedgefound = False;#edgechange=True
                         continue
                     elif e1concaveedgefound and not e2concaveedgefound:
                         _,t,mp4 = self.footprint(mine,rmine,e1plus)
                         e1=e1plus;e2=mine;tin=t;rin=-float('inf')
                         mplast[1] = e1;
                         mplast[2] = t;
                         mplast[3] = e2

                         print "continue with ",e1,t,e2
                         e1concaveedgefound = False;#edgechange=True
                         continue

                     else:#case #2
                         self.mps.append([mpcurr[0],mine,rmine,e2,r,mpcurr[-1]])
                         e2 = mine;rin=-float('inf')
                         mplast[1] = e1;
                         mplast[2] = t;
                         mplast[3] = e2

                         #edgechange=True
                         continue

                 elif concavedist:
                     print "dist failure with concavevert", mini,distconv,radius,e1concaveedgefound,e2concaveedgefound
                     if not edgechange:
                        mpcurr = self.retractConcave(mplast, mpcurr, distconv, mini)
                     if e1concaveedgefound and not e2concaveedgefound:#case # 3
                        _, r3, mp3 = self.footprint(mpcurr[3], mpcurr[4],mini)
                        _, p0,t0 = self.distance(e1plus,mpcurr[0],True)
                        self.mps.append([mpcurr[0],mini,r3,mpcurr[3],mpcurr[4],mpcurr[-1]])
                        print "branching to ",mini,r3,mpcurr[3]
                        e1 = e1plus;t=t0;e2=self.prevedge(mini);tin=t;rin=-float('inf')
                        mplast[1]=e1;mplast[2]=t;mplast[3]=e2
                        print "proceeding with",e1,t,e2
                        e1concaveedgefound=False;e2concaveedgefound=False;self.handledconcaveedges.append(mini)
                        #edgechange=True
                        continue

                     elif e2concaveedgefound and not e1concaveedgefound:#case #4
                         _, r3, mp3 = self.footprint(mpcurr[1], mpcurr[2], self.prevedge(mini))
                         radius, p0,t0 = self.distance(mini, mpcurr[0],True)
                         _,r4,mp4=self.footprint(mini,t0,e2minus)
                         self.mps.append([mpcurr[0],mini,t0,e2minus,r4,radius])
                         print "branching to",mini,t0,e2minus
                         e2=self.prevedge(mini)
                         mplast[3]=e2
                         rin = -float('inf')
                         e2concaveedgefound=False;self.handledconcaveedges.append(mini)
                         print "continue with",e1,t,e1concaveedgefound,e2concaveedgefound
                         #edgechange=True
                         continue
                     elif not e1concaveedgefound and not e2concaveedgefound: #case#1
                         _,r3,mp3=self.footprint(e2,r,mini)
                         self.mps.append([mpcurr[0],mini,r3,e2,r,mpcurr[-1]])
                         print "branching to ",mini,r3,e2
                         e2=self.prevedge(mini)
                         mplast[3]=e2
                         rin = -float('inf')
                         print "continue with ",e1,t,e2
                         self.handledconcaveedges.append(mini)
                         #edgechange = True
                         continue

                     elif e1concaveedgefound and e2concaveedgefound:
                         #check e1,e3
                         self.handledconcaveedges.append(mini)
                         dot = np.dot(self.normal(e1,t,True),self.normal(mini,0.0,True))
                         print dot
                         if dot<0.01:
                             #_,r3,mp3=self.footprint(e1,t,self.prevedge(mini))
                             radius,fp,t0=self.distance(mini,mpcurr[0],True)
                             _,r4,mp4 = self.footprint(mini,t0,e2minus)
                             self.mps.append([mpcurr[0], mini, t0, e2minus, r4, radius])
                             print "branching to ",mini,t0,e2minus
                             e2=self.prevedge(mini);
                             rin = -float('inf')
                             e2concaveedgefound = False
                             mplast[3] = e2

                             #edgechange=True
                             print "continuing with ",e1,t,e2
                             continue
                         else:
                             radius,pi,t0=self.distance(e1plus,mpcurr[0],True,0.0)
                             _,r0,mp4 = self.footprint(e1plus,t0,self.prevedge(mini))
                             self.mps.append([mpcurr[0],e1plus,t0,self.prevedge(mini),r0,radius])
                             print "branching to ", e1plus, t0, self.prevedge(mini)
                             radius, pi, t = self.distance(mini, mpcurr[0], True)

                             e1 = mini;tin=t;
                             #edgechange=True;
                             mplast[1]=e1;mplast[2]=t;rin=-float('inf')
                             e1concaveedgefound=False
                             print "continue with",e1,t,e2
                             continue


                     #e2concaveedgefound=False

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
             if t >= self.its[e1].end() and not e1concaveedgefound:

                 e1p = self.nextedge(e1);

                 if e1p in self.convex_vert:
                     convexedgefound=True
                     break
                 if e1p in self.concave_vert:
                     e1concaveedgefound = True
                     currconcaveedge = e1p
                     self.handledconcaveedges.append(e1p)
                 else:
                     if self.its[e1].end()-self.its[e1].begin()<0.1:
                         self.tracededges.append(e1)
                     #edgechange = True
                     e1 = e1p
                     t = self.its[e1].begin()
                     mplast[1]=e1;mplast[2]=t
                     tin=t
                     print e1,e2
                     self.apairs(e1, e2)

                     if e1 in self.tracededges:
                         terminateloop=True
                     #self.tracededges.append(e1)
    def proceed3(self, mps, edges1,eds):
        edges=edges1
        e1 = mps[1];
        t = mps[2];
        e2 = mps[3];
        deltat = 0.01
        e1concaveedgefound = False;
        e2concaveedgefound = False;

        mplast = mps
        edgechange = False
        distFailure = False
        currconcaveedge = -1
        convexedgefound = False
        if str(e1) + "-" + str(e2) in self.pairs:
            print "sorry this pair traced", e1, e2
            return

        if e1 in self.tracededges:
            if self.debug: print "sorry", e1, " traced";
            return
        if e2 in self.tracededges:
            if self.debug: print "sorry", e2, "traced";
            return
        convexedge = False
        if e2 == self.prevedge(e1): convexedge=True

        terminateloop = False
        e1start = e1;
        e2start = e2;
        tin = t;
        rin = -float('inf')
        print e1, e2
        noofsteps = 0

        while not convexedgefound and not terminateloop:
            # while t < self.its[e1].end():
            self.apairs(e1, e2)

            if e1 in self.tracededges or e2 in self.tracededges:
                break
            e1plus = self.nextedge(e1);
            e2minus = self.prevedge(e2)
            noofsteps += 1
            _, r, mp = self.footprint(e1, t, e2)
            #print e1,t,e2,r
            # Here, we take care of more than one intervals and search if t exists within intervals or not
            if t>=0.0 and t<=1.0 and len(self.its[e1].search(t))==0:
                print e1,t,self.its[e1]
                if self.its[e1].end()-self.its[e1].begin()>0.01:
                    sortedtree = sorted(self.its[e1])
                    i=0;tfound=False
                    for i,iv in enumerate(sortedtree):
                        if iv.begin>t:
                            t=iv.begin
                            tfound=True
                            print "t=",t
                            break
                    if tfound: continue
                    else: print "t not found";break
                else:
                    break



            if t > 3.0:
                print "t>3"
                return
            if r == float('inf'):
                print "r inf",e1,t,e2
                return
            if t > tin:
                #print "e1",e1,"chopped"
                self.its[e1].chop(tin, t)
            if rin > r:
                #print "e2", e2, "chopped"

                self.its[e2].chop(r, rin)
            if rin == -float('inf'): rin = r
            radius = max(np.linalg.norm(mp - self.p(e1, t)), np.linalg.norm(mp - self.p(e2, r)))
            if radius>self.maxradius*1.5 and self.maxradius>1.0: print "radius high";break
            mpcurr = [mp, e1, t, e2, r, radius]
            if r <= 0.0 and not e2concaveedgefound :
                '''
                if e2concaveedgefound and r < -1.5:
                    e2 = self.prevedge(e2)
                    e2concaveedgefound = False
                    continue
                '''

                # print "r negative??"
                if e2 in self.convex_vert:
                    convexedgefound = True
                    break
                if e2 in self.concave_vert :
                    if e2  in self.handledconcaveedges and t>1.0: break
                    if not(e2 in self.handledconcaveedges and t<=1.0):
                        e2concaveedgefound = True
                        currconcaveedge = e2
                        self.handledconcaveedges.append(e2)
                        #self.findCounterPart(e2,radius,False)
                        # continue
                else:
                    if self.its[e2].end() - self.its[e2].begin() < 0.1:
                        self.tracededges.append(e2)
                    e2p = self.prevedge(e2);
                    e2 = e2p
                    print e1, e2
                    self.apairs(e1, e2)
                    if e2 in self.tracededges:
                        terminateloop = True
                    continue

            #print "e1,t,e2", e1, t, e2, r

            if edgechange:
                mplast = mpcurr;
                edgechange = False
            noofsteps += 1
            #print e1,e2,t,r,radius
            edges.append(mpcurr)
            t += deltat
            mplast = mpcurr
            e1prev=e1;e2prev=e2;t0=t
            noofsteps += 1


            if  e1concaveedgefound or e2concaveedgefound:  # or noofsteps>5 or edgechange:
                noofsteps = 0
                if e1concaveedgefound:
                    resultmp1 = self.resultMp(e1,e2,radius,True)
                    mpcurr = self.extendCv(e1,e2,t,resultmp1,radius,edges)

                elif e2concaveedgefound:
                    resultmp1 = self.resultMp(e2,e1,radius,False)
                    mpcurr = self.extendCv(e1, e2, t, resultmp1, radius, edges)

                if e2concaveedgefound and not e1concaveedgefound:
                    e2concaveedgefound=False

                    #print "checking"
                    cp=self.findCounterPart1(e2,e1, radius,False)
                    print e1,e2,"counterpart new",cp

                    #mine2, mindist2,rmine2, minc2, mindistc2 = self.findCounterPart(e2,e1, radius,False)

                    resultMp=self.resultMp(e2,e1, radius,False)

                        #print "adding here"

                    if cp[0] and len(cp[1])>0 : # A counterpart concave edge has been found
                        minc2=cp[1][1]
                        #self.mps.append([mpcurr[0],e1, t, self.prevedge(minc2), 1.0,  mpcurr[-1]])

                        concavedist = True
                        print "e2:will branch to concavec", minc2,self.prevedge(e2)
                        self.handledconcaveedges.append(minc2)
                        _,_,t=self.distance(minc2,resultMp,True)
                        e1=minc2;e2=self.prevedge(e2)
                        _,r1,mpfin = self.footprint(e1,t,e2)
                        mpcurr = self.extendCv(e1prev, e2prev, t0, mpfin, radius, edges)
                        if self.checkPair(e1,t,e2):
                            self.mps.append([mpcurr, e1, t, e2, r1, radius])
                            print "e2:shall have a branch", e1, e2, t


                        e1n=e1prev;t1=mpcurr[2];e2n=self.prevedge(minc2)
                        if  self.checkPair(e1n,t1,e2n):
                            _,rn,mpn=self.footprint(e1n,t1,e2n)
                            e1=e1n;e2==e2n;t=t1
                            tin = t;
                            rin = rn
                            edges = [];
                            eds.append(edges)
                            print "e2:continue with ",e1n,e2n,t1
                            self.apairs(e1,e2)
                            continue
                        else: break
                    elif  not cp[0] and len(cp[1])>0:# An edge has been found which would take it rightwards
                        mine2=cp[1][1];
                        _,_,rmine2=self.distance(mine2,resultMp)
                        #rmine2=cp[1][2]
                        e2=self.prevedge(e2)
                        print "e2:wil branch to edge",mine2,e2,rmine2

                        e1=mine2;t=rmine2;
                        tin = t;
                        rin = -float('inf')
                        count1 = 0;

                        if self.checkPair(e1,t,e2):
                            #self.apairs(e1, e2)
                            _, r1, mpfin = self.footprint(e1, t, e2)
                            mpcurr = self.extendCv(e1prev, e2prev, t0, mpfin, radius, edges)
                            self.mps.append([mpcurr, e1, t, e2, r1, radius])
                        else:
                            #print "could not verify",e1,e2
                            self.mps.append([mpcurr, e1, 0.0, e2, -float('inf'), radius])
                        if not mine2 in self.convex_vert:
                            e1=e1prev;t=mpcurr[2];e2=self.prevedge(mine2)
                            if self.checkPair(e1,t,e2):
                                continue
                            else:
                                break
                        else: break
                    else: break
                elif e1concaveedgefound and not e2concaveedgefound:
                    #print "here"
                    e1concaveedgefound=False
                    cp = self.findCounterPart1(e1, e2, radius, True)
                    print e1, e2, "counterpart new", cp

                    #mine1, mindist1,rmine1, minc1, mindistc1 = self.findCounterPart(e1,e2, radius,True)
                    #print mine1, mindist1, rmine1, minc1, mindistc1,radius
                    resultMp = self.resultMp(e1, e2, radius, True)
                    #self.extendCv(e1prev, e2prev, t0, resultMp, radius, edges)

                    #print "e1:abs difference",e1,minc1,mine1,abs(mindistc1-mindist1)
                    if cp[0] and len(cp[1])>0: # a counterpart concave edge has been found
                        minc1=cp[1][1]
                        self.handledconcaveedges.append(minc1)

                        e1n=self.nextedge(e1);e2n=self.prevedge(minc1)
                        _, _, t1 = self.distance(e1n, resultMp, True)
                        _, r1, mpfin = self.footprint(e1n, t1, e2n)
                        mpcurr = self.extendCv(e1prev, e2prev, t0, mpfin, radius, edges)
                        self.mps.append([mpcurr[0], e1n, t1, e2n, r1, mpcurr[-1]])

                        concavedist = True
                        print "e1:will branch to concave", e1n, e2n
                        #TODO
                        e1n= minc1
                        _,_,t3=self.distance(e1n,mpfin,True)
                        if self.checkPair(e1n,t3,e2prev):
                            _,rn,_ = self.footprint(e1n,t3,e2prev)
                            #self.mps.append([mpfin,e1n,t3,e2prev,rn,radius])
                            e1=e1n;t=t3;tin=t;e2=e2prev;rin=rn
                            print "e1:will have a branch of ",e1n,e2prev,t3
                            self.apairs(e1, e2)

                            continue
                        else:
                            print "failed ",e1n,e2prev
                            break
                    elif  not cp[0] and len(cp[1])>0:

                        mine1=cp[1][1];
                        _,_,rmine1=self.distance(mine1,resultMp)
                        #rmine1=cp[1][2]
                        e1=self.nextedge(e1);e2=mine1
                        print "e1:branch to edge", e1, mine1

                        _, _, t = self.distance(e1, resultMp, True)
                        tin = t;
                        rin = -float('inf')
                        _, rfin, mpfin = self.footprint(e1, t, e2)
                        mpcurr = self.extendCv(e1prev, e2prev, t0, mpfin, radius, edges)
                        self.mps.append([mpcurr, e1, t, e2, rfin, radius])

                        if not e2 in self.convex_vert:
                            if rmine1>1.0:
                                e1n=self.nextedge(mine1)
                                _,_,tn=self.distance(e1n,resultMp,True)
                            else:
                                e1n=mine1;tn=rmine1
                            print "e1:try continuing from here",e1n,e2prev,tn
                            e1=e1n;e2=e2prev;t=tn;
                            if self.checkPair(e1,t,e2):
                                self.apairs(e1, e2)
                                continue
                            else:
                                break
                        else:
                            break
                    else: break
                elif e1concaveedgefound and e2concaveedgefound:
                    e2concaveedgefound = False
                    e1concaveedgefound = False
                    e1prev=e1;e2prev=e2

                    # print "checking"

                    #mine2, mindist2, rmine2, minc2, mindistc2 = self.findCounterPart(e2, e1, radius, False)
                    #mine1, mindist1, rmine1, minc1, mindistc1 = self.findCounterPart(e1, e2, radius, True)
                    cp1 = self.findCounterPart1(e1, e2, radius, True)
                    print e1, e2, "counterpart new", cp1
                    cp2 = self.findCounterPart1(e2, e1, radius, False)
                    print e1, e2, "counterpart new", cp2

                    resultMp = self.resultMp(e2, e1, radius, False)


                    if cp2[0] and len(cp2[1])>0: #counterpart concave edge of cp2 found
                        minc2=cp2[1][1]
                        concavedist = True
                        self.handledconcaveedges.append(minc2)

                        print "e1e2:will branch to", minc2, self.prevedge(e2)
                        _, _, t = self.distance(minc2, resultMp, True)
                        e1 = minc2;
                        e2 = self.prevedge(e2)
                        _, rfin, mpfin = self.footprint(e1, t, e2)
                        self.extendCv(e1prev,e2prev,t0,mpfin,radius,edges)
                        self.mps.append([mpcurr[0],e1, t, e2, 1.0,  mpcurr[-1]])

                    elif  not cp2[0] and len(cp2[1])>0:
                        mine2=cp2[1][1];
                        _,_,rmine2=self.distance(mine2,resultMp)
                        #rmine2=cp2[1][2]
                        print "e1e2:wil branch to", mine2, self.prevedge(e2)
                        e1 = mine2;
                        t = rmine2;
                        e2=self.prevedge(e2)
                        tin = t;
                        rin = -float('inf')
                        _, rfin, mpfin = self.footprint(e1, t, e2)
                        self.extendCv(e1prev,e2prev,t0,mpfin,radius,edges)
                        self.mps.append([mpcurr[0], mine2, rmine2, e2, 0.9, mpcurr[-1]])
                        #edges = [];
                        #eds.append(edges)

                        #continue

                    if cp1[0] and len(cp1[1])>0:
                        minc1=cp1[1][1]
                        concavedist = True
                        self.handledconcaveedges.append(minc1)

                        print "e1e2:will branch to concave", self.nextedge(e1prev), self.prevedge(minc1)
                        #TODO

                        #self.mps.append([mpcurr[0],self.nextedge(e1),0.0,minc1,1.0,mpcurr[-1]])
                        e2 = self.prevedge(minc1);

                        e1 = self.nextedge(e1prev)
                        _, _, t = self.distance(e1, resultMp, True)
                        _, rfin, mpfin = self.footprint(e1, t, e2)
                        mpcurr = self.extendCv(e1prev, e2prev, t0, mpfin, radius, edges)
                        self.mps.append([mpcurr[0], e1, t, e2, 1.0, mpcurr[-1]])

                        break


                    elif not cp1[0] and len(cp1[1])>0:
                        mine1=cp1[1][1];rmine1=cp1[1][2]
                        print "e1e2:branch to edge",self.nextedge(e1prev),mine1
                        e1=self.nextedge(e1prev);e2=mine1
                        _, _, t = self.distance(e1, resultMp, True)
                        tin = t;
                        _, rfin, mpfin = self.footprint(e1, t, e2)
                        mpcurr=self.extendCv(e1prev,e2prev,t0,mpfin,radius,edges)
                        self.mps.append([mpcurr[0], e1, t, e2, 1.0, mpcurr[-1]])

                        break
                    continue
            elif  convexedge:# or  e2==self.prevedge(e1):
                #print "we check distance"

                noofsteps = 0
                skipedges = [e1, e2]
                if not e1plus in self.convex_vert: skipedges.append(e1plus)
                #if not e2minus in self.convex_vert: skipedges.append(e2minus)
                mine, mindist, rmine = self.nearestEdge(mp, skipedges)
                #print mine,mindist,rmine,mp
                skipconcaves = []
                # we are not going to check distances with concaveedges of next branch
                if not e2concaveedgefound:
                    e0 = e2minus
                    while not e0 in self.convex_vert and e0 != e2:
                        if e0 in self.concave_vert:
                            skipconcaves.append(e0)
                        e0 = self.prevedge(e0)
                if not e1concaveedgefound:
                    e0 = e1plus
                    while not e0 in self.convex_vert and e0 != e1:
                        if e0 in self.concave_vert:
                            skipconcaves.append(e0)
                        e0 = self.nextedge(e0)

                mini, distconv = self.nearestConcaveVert(mp, e1, e2, e1concaveedgefound, e2concaveedgefound,
                                                         skipconcaves)
                concavedist = False;
                 #print "e1,t,e2,mindist,distconv,radius",e1,t,mindist,distconv,radius
                #print "radius=",radius
                if e1concaveedgefound:
                    if distconv < np.linalg.norm(mpcurr[0] - self.vertices[e1plus]): concavedist = True
                elif e2concaveedgefound:
                    if distconv < np.linalg.norm(mpcurr[0] - self.vertices[e2]): concavedist = True
                else:
                    if distconv < radius: concavedist = True

                #print e1,t,e2,concavedist,e1concaveedgefound,e2concaveedgefound



                #print "radius=",radius

                if  mindist < radius and not concavedist:
                    convexedge=False
                    print "distance failure with ", mine,radius,mindist
                    if self.maxradius<1.0 or (self.maxradius<radius*1.2 and radius<self.maxradius*1.5):
                        print self.maxradius,radius
                        self.maxradius=radius*1.2
                    if self.debug3:
                        bz.pe(
                            [[[[self.p(mine, rmine)]], 'g', ''], [[[mpcurr[0]]], 'r', ''], [[[self.p(e1, t)]], 'b', ''],
                             [self.edges]])
                    '''
                    if mine in self.convex_vert:
                        if mine ==e1plus:
                           self.tracededges.append(e1)
                        elif mine==e2minus:
                            self.tracededges.append(e2)
                    '''
                    if not edgechange:
                        mpcurr = self.retract(mplast, mpcurr, mindist, mine)

                    if True:  # case #2
                        self.mps.append([mpcurr[0], mine, rmine, e2, r, mpcurr[-1]])
                        e2 = mine;
                        rin = -float('inf')
                        mplast[1] = e1;
                        mplast[2] = t;
                        mplast[3] = e2

                        # edgechange=True
                        self.apairs(e1, e2)

                        continue

                elif  concavedist:
                    convexedge=False

                    print "dist failure with concavevert", mini, distconv, radius, e1concaveedgefound, e2concaveedgefound
                    if self.maxradius < 1.0 or (self.maxradius < radius * 1.2 and radius < self.maxradius * 1.5):
                        print self.maxradius, radius
                        self.maxradius = radius * 1.2

                    if not edgechange:
                        mpcurr = self.retractConcave(mplast, mpcurr, distconv, mini)
                    if True:  # case#1
                        _, r3, mp3 = self.footprint(e2, r, mini)
                        self.mps.append([mpcurr[0], mini, r3, e2, r, mpcurr[-1]])
                        print "branching to ", mini, r3, e2
                        e2 = self.prevedge(mini)
                        mplast[3] = e2
                        rin = -float('inf')
                        print "continue with ", e1, t, e2
                        self.handledconcaveedges.append(mini)
                        self.apairs(e1, e2)

                        continue


                    mplast = mpcurr
            #print "adding",mpcurr[0],mpcurr[1],mpcurr[3],mpcurr[2]
            if t >= 1.0 and not e1concaveedgefound :
                e1p = self.nextedge(e1);

                if e1p in self.convex_vert:
                    convexedgefound = True
                    break
                if e1p in self.concave_vert :
                    if e1p  in self.handledconcaveedges and r<0.0:
                        break
                    if not (e1p in self.handledconcaveedges and r >= 0.0):

                        e1concaveedgefound = True
                        currconcaveedge = e1p
                        self.handledconcaveedges.append(e1p)
                        # self.findCounterPart(e1,radius)
                else:
                    if self.its[e1].end() - self.its[e1].begin() < 0.1:
                        self.tracededges.append(e1)
                    # edgechange = True
                    e1 = e1p
                    t = self.its[e1].begin()
                    mplast[1] = e1;
                    mplast[2] = t
                    tin = t;rin=-float('inf')
                    print e1, e2
                    self.apairs(e1, e2)

                    if e1 in self.tracededges:
                        terminateloop = True
                        # self.tracededges.append(e1)
        return

    def proceed2(self, mps, edges1,eds):
        edges=edges1
        e1 = mps[1];
        t = mps[2];
        e2 = mps[3];
        deltat = 0.01
        e1concaveedgefound = False;
        e2concaveedgefound = False;

        mplast = mps
        edgechange = False
        distFailure = False
        currconcaveedge = -1
        convexedgefound = False
        if str(e1) + "-" + str(e2) in self.pairs:
            print "sorry this pair traced",e1,e2
            return
        if e1 in self.tracededges:
            if self.debug: print "sorry", e1, " traced";
            return
        if e2 in self.tracededges:
            if self.debug: print "sorry", e2, "traced";
            return

        terminateloop = False
        e1start = e1;
        e2start = e2;
        tin = t;
        rin = -float('inf')
        print e1, e2
        self.apairs(e1, e2)
        noofsteps = 0

        while not convexedgefound and not terminateloop:
            # while t < self.its[e1].end():
            if e1 in self.tracededges or e2 in self.tracededges:
                break
            e1plus = self.nextedge(e1);
            e2minus = self.prevedge(e2)
            noofsteps += 1
            _, r, mp = self.footprint(e1, t, e2)
            #print e1,t,e2,r
            if t>=0.0 and t<=1.0 and self.its[e1].search(t)==set(): break
            if t > 3.0:
                print "t>3"
                return
            if r == float('inf'):
                print "r inf",e1,t,e2
                return
            if t > tin:
                self.its[e1].chop(tin, t)
            if rin > r:
                self.its[e2].chop(r, rin)
            if rin == -float('inf'): rin = r
            radius = max(np.linalg.norm(mp - self.p(e1, t)), np.linalg.norm(mp - self.p(e2, r)))
            if radius>self.maxradius and self.maxradius>5.0: break
            mpcurr = [mp, e1, t, e2, r, radius]
            if r <= 0.0 and not e2concaveedgefound:
                if e2concaveedgefound and r < -1.5:
                    e2 = self.prevedge(e2)
                    e2concaveedgefound = False
                    continue

                # print "r negative??"
                if e2 in self.convex_vert:
                    convexedgefound = True
                    break
                if e2 in self.concave_vert:
                    e2concaveedgefound = True
                    currconcaveedge = e2
                    self.handledconcaveedges.append(e2)
                    #self.findCounterPart(e2,radius,False)
                    # continue
                else:
                    if self.its[e2].end() - self.its[e2].begin() < 0.1:
                        self.tracededges.append(e2)
                    e2p = self.prevedge(e2);
                    e2 = e2p
                    print e1, e2
                    self.apairs(e1, e2)
                    if e2 in self.tracededges:
                        terminateloop = True
                    continue

            #print "e1,t,e2", e1, t, e2, r

            if edgechange:
                mplast = mpcurr;
                edgechange = False
            noofsteps += 1
            #print e1,e2,t,r,radius
            edges.append(mpcurr)
            t += deltat
            mplast = mpcurr
            e1prev=e1;e2prev=e2;t0=t
            if t >= 1.0 and not e1concaveedgefound:
                e1p = self.nextedge(e1);

                if e1p in self.convex_vert:
                    convexedgefound = True
                    break
                if e1p in self.concave_vert :
                    e1concaveedgefound = True
                    currconcaveedge = e1p
                    self.handledconcaveedges.append(e1p)
                    #self.findCounterPart(e1,radius)
                else:
                    if self.its[e1].end() - self.its[e1].begin() < 0.1:
                        self.tracededges.append(e1)
                    # edgechange = True
                    e1 = e1p
                    t = self.its[e1].begin()
                    mplast[1] = e1;
                    mplast[2] = t
                    tin = t
                    print e1, e2
                    self.apairs(e1, e2)

                    if e1 in self.tracededges:
                        terminateloop = True
                        # self.tracededges.append(e1)
            noofsteps += 1


            if  e1concaveedgefound or e2concaveedgefound:  # or noofsteps>5 or edgechange:
                noofsteps = 0
                if e2concaveedgefound and not e1concaveedgefound:
                    e2concaveedgefound=False

                    #print "checking"
                    print e1,e2,"counterpart new",self.findCounterPart1(e2,e1, radius,False)

                    mine2, mindist2,rmine2, minc2, mindistc2 = self.findCounterPart(e2,e1, radius,False)

                    resultMp=self.resultMp(e2,e1, radius,False)

                        #print "adding here"

                    if mindistc2 <= mindist2 or abs(mindistc2-mindist2)<30.0:
                        #self.mps.append([mpcurr[0],e1, t, self.prevedge(minc2), 1.0,  mpcurr[-1]])

                        concavedist = True
                        print "e2:will branch to concavec", minc2,self.prevedge(e2)
                        self.handledconcaveedges.append(minc2)
                        _,_,t=self.distance(minc2,resultMp,True)
                        e1=minc2;e2=self.prevedge(e2)
                        _,r1,mpfin = self.footprint(e1,t,e2)
                        tin = t;
                        rin = -float('inf')
                        count1 = 0;
                        #t0 += deltat
                        print t0
                        _,_,tfin = self.distance(e1prev,mpfin,True)
                        #tempedges=[]
                        #tempedges.append([mpfin])
                        #bz.pe([[tempedges],[self.edges,'g','']])
                        deltat1=0.01;count1=0
                        while t0<tfin and count1<50:
                            _, r1, mp1 = self.footprint(e1prev, t0, e2prev)
                            print e1prev,t0,e2prev,r1
                            mpcurr = [mp1, e1prev, t0, e2prev, r1, radius]
                            t0 += deltat1;
                            count1 += 1
                            edges.append(mpcurr)
                        edges=[];eds.append(edges)
                        e1n=e1prev;t1=t0;e2n=self.prevedge(minc2)
                        if self.checkPair(e1n,t1,e2n):
                            _,rn,mpn=self.footprint(e1n,t1,e2n)
                            self.mps.append([mpfin,e1n,t1,e2n,rn,radius])
                            print "e2:shall have a branch",e1n,e2n
                        continue
                    elif  mindist2<mindistc2:
                        print "e2:wil branch to edge",mine2,e2
                        e1=mine2;t=rmine2;
                        tin = t;
                        rin = -float('inf')
                        count1 = 0;
                        _, r1, mpfin = self.footprint(e1, t, e2)

                        t0 += deltat
                        _,_,tfin = self.distance(e1,mpfin,True)
                        while t0<tfin and count1 < 50:
                            _, r1, mp1 = self.footprint(e1prev, t0, e2prev)
                            mpcurr = [mp1, e1prev, t0, e2prev, r1, radius]
                            t0 += deltat;
                            count1 += 1
                            edges.append(mpcurr)

                        #self.mps.append([mpcurr[0], mine2, rmine2, self.prevedge(e2), 1.0, mpcurr[-1]])
                        edges = [];
                        eds.append(edges)
                        if self.checkPair(e1,t,e2):

                            continue
                        else:
                            break
                elif e1concaveedgefound and not e2concaveedgefound:
                    #print "here"
                    e1concaveedgefound=False

                    #print "checking"
                    print e1, e2, "counterpart new", self.findCounterPart1(e1, e2, radius, True)

                    mine1, mindist1,rmine1, minc1, mindistc1 = self.findCounterPart(e1,e2, radius,True)
                    #print mine1, mindist1, rmine1, minc1, mindistc1,radius
                    resultMp = self.resultMp(e1, e2, radius, True)
                    print "e1:abs difference",e1,minc1,mine1,abs(mindistc1-mindist1)
                    if mindistc1 <= mindist1 or abs(mindistc1-mindist1)<30.0:
                        self.handledconcaveedges.append(minc1)

                        e1n=self.nextedge(e1);e2n=self.prevedge(minc1)
                        _, _, t1 = self.distance(e1n, resultMp, True)
                        if not self.checkPair(e1n,t1,e2n):
                            break

                        #self.mps.append([mpcurr[0], minc1, t1, e2, 1.0, mpcurr[-1]])

                        concavedist = True
                        print "e1:will branch to concave", e1n, e2n
                        #TODO

                        #self.mps.append([mpcurr[0],self.nextedge(e1),0.0,minc1,1.0,mpcurr[-1]])
                        e2 = e2n;

                        e1 = e1n;t=t1
                        _,rfin,mpfin = self.footprint(e1,t,e2)
                        count1 = 0;
                        #t0 += deltat
                        deltat1=0.01
                        _,_,tfin=self.distance(e1prev,mpfin,True)
                        while t0<tfin and count1 < 100:
                            _, r1, mp1 = self.footprint(e1prev, t0, e2prev)
                            mpcurr = [mp1, e1prev, t0, e2prev, r1, radius]
                            t0 += deltat1;
                            count1 += 1
                            edges.append(mpcurr)
                        tin=t;rin=-float('inf')
                        edges = [];
                        eds.append(edges)
                        e1n= minc1
                        _,_,t0=self.distance(e1n,mpfin,True)
                        if self.checkPair(e1n,t0,e2prev):
                            _,_,t1=self.distance(e1n,mpfin,True)
                            _,rn,_ = self.footprint(e1n,t0,e2)
                            self.mps.append([mpfin,e1n,t0,e2prev,rn,radius])
                            print "will have a branch of ",e1n,e2prev
                        else:
                            print "failed ",e1n,e2prev

                        continue


                    elif  mindist1<mindistc1:
                        print "e1:branch to edge",self.nextedge(e1),mine1
                        e1=self.nextedge(e1);e2=mine1
                        _, _, t = self.distance(e1, resultMp, True)
                        tin = t;
                        _, rfin, mpfin = self.footprint(e1, t, e2)
                        count1 = 0;
                        #t0 += deltat
                        deltat1=0.01
                        _,_,tfin=self.distance(e1prev,mpfin,True)
                        while t0<tfin and count1 < 100:
                            _, r1, mp1 = self.footprint(e1prev, t0, e2prev)
                            mpcurr = [mp1, e1prev, t0, e2prev, r1, radius]
                            t0 += deltat1;
                            count1 += 1
                            edges.append(mpcurr)
                        rin = -float('inf')
                        edges = [];
                        eds.append(edges)

                        continue
                #if e1concaveedgefound and e2concaveedgefound:
                 #   e1concaveedgefound=e2concaveedgefound=False
                  #  break
                elif e1concaveedgefound and e2concaveedgefound:
                    e2concaveedgefound = False
                    e1concaveedgefound = False
                    e1prev=e1;e2prev=e2

                    # print "checking"

                    mine2, mindist2, rmine2, minc2, mindistc2 = self.findCounterPart(e2, e1, radius, False)
                    mine1, mindist1, rmine1, minc1, mindistc1 = self.findCounterPart(e1, e2, radius, True)

                    resultMp = self.resultMp(e2, e1, radius, False)


                    if mindistc2 <= mindist2 or abs(mindistc2-mindist2)<30.0:
                        concavedist = True
                        self.handledconcaveedges.append(minc2)

                        print "e1e2:will branch to", minc2, self.prevedge(e2)
                        _, _, t = self.distance(minc2, resultMp, True)
                        e1 = minc2;
                        e2 = self.prevedge(e2)
                        _, rfin, mpfin = self.footprint(e1, t, e2)
                        count1 = 0;
                        #t0 += deltat
                        deltat1=0.01
                        _,_,tfin=self.distance(e1prev,mpfin,True)
                        while t0<tfin and count1 < 100:
                            _, r1, mp1 = self.footprint(e1prev, t0, e2prev)
                            mpcurr = [mp1, e1prev, t0, e2prev, r1, radius]
                            t0 += deltat1;
                            count1 += 1
                            edges.append(mpcurr)

                        #tin = t;
                        #rin = -float('inf')

                        self.mps.append([mpcurr[0],e1, t, e2, 1.0,  mpcurr[-1]])
                        #edges = [];
                        #eds.append(edges)

                        #continue
                    elif  mindist2<mindistc2:
                        print "e1e2:wil branch to", mine2, self.prevedge(e2)
                        e1 = mine2;
                        t = rmine2;
                        e2=self.prevedge(e2)
                        tin = t;
                        rin = -float('inf')
                        _, rfin, mpfin = self.footprint(e1, t, e2)
                        count1 = 0;
                        #t0 += deltat
                        deltat1=0.01
                        _,_,tfin=self.distance(e1prev,mpfin,True)
                        while t0<tfin and count1 < 100:
                            _, r1, mp1 = self.footprint(e1prev, t0, e2prev)
                            mpcurr = [mp1, e1prev, t0, e2prev, r1, radius]
                            t0 += deltat1;
                            count1 += 1
                            edges.append(mpcurr)

                        self.mps.append([mpcurr[0], mine2, rmine2, e2, 1.0, mpcurr[-1]])
                        #edges = [];
                        #eds.append(edges)

                        #continue

                    if mindistc1 <= mindist1 or abs(mindistc1-mindist1)<30.0:
                        concavedist = True
                        self.handledconcaveedges.append(minc1)

                        print "e1e2:will branch to concave", self.nextedge(e1prev), self.prevedge(minc1)
                        #TODO

                        #self.mps.append([mpcurr[0],self.nextedge(e1),0.0,minc1,1.0,mpcurr[-1]])
                        e2 = self.prevedge(minc1);

                        e1 = self.nextedge(e1prev)
                        _, _, t = self.distance(e1, resultMp, True)
                        _, rfin, mpfin = self.footprint(e1, t, e2)
                        count1 = 0;
                        #t0 += deltat
                        deltat1=0.01
                        _,_,tfin=self.distance(e1prev,mpfin,True)
                        while t0<tfin and count1 < 100:
                            _, r1, mp1 = self.footprint(e1prev, t0, e2prev)
                            mpcurr = [mp1, e1prev, t0, e2prev, r1, radius]
                            t0 += deltat1;
                            count1 += 1
                            edges.append(mpcurr)

                        tin=t;rin=-float('inf')
                        edges = [];
                        eds.append(edges)

                        continue


                    elif mindist1<mindistc1:
                        print "e1e2:branch to edge",self.nextedge(e1prev),mine1
                        e1=self.nextedge(e1prev);e2=mine1
                        _, _, t = self.distance(e1, resultMp, True)
                        tin = t;
                        _, rfin, mpfin = self.footprint(e1, t, e2)
                        count1 = 0;
                        #t0 += deltat
                        deltat1=0.01
                        tfin=self.distance(e1prev,mpfin,True)
                        while t0<tfin and count1 < 100:
                            _, r1, mp1 = self.footprint(e1prev, t0, e2prev)
                            mpcurr = [mp1, e1prev, t0, e2prev, r1, radius]
                            t0 += deltat1;
                            count1 += 1
                            edges.append(mpcurr)

                        rin = -float('inf')
                        edges = [];
                        eds.append(edges)

                        continue
                    continue
            elif   e2==self.prevedge(e1):
                noofsteps = 0
                skipedges = [e1, e2]
                if not e1plus in self.convex_vert: skipedges.append(e1plus)
                if not e2minus in self.convex_vert: skipedges.append(e2minus)
                mine, mindist, rmine = self.nearestEdge(mp, skipedges)
                #print mine,mindist,rmine,mp
                skipconcaves = []
                # we are not going to check distances with concaveedges of next branch
                if not e2concaveedgefound:
                    e0 = e2minus
                    while not e0 in self.convex_vert and e0 != e2:
                        if e0 in self.concave_vert:
                            skipconcaves.append(e0)
                        e0 = self.prevedge(e0)
                if not e1concaveedgefound:
                    e0 = e1plus
                    while not e0 in self.convex_vert and e0 != e1:
                        if e0 in self.concave_vert:
                            skipconcaves.append(e0)
                        e0 = self.nextedge(e0)

                mini, distconv = self.nearestConcaveVert(mp, e1, e2, e1concaveedgefound, e2concaveedgefound,
                                                         skipconcaves)
                concavedist = False;
                 #print "e1,t,e2,mindist,distconv,radius",e1,t,mindist,distconv,radius
                #print "radius=",radius
                if e1concaveedgefound:
                    if distconv < np.linalg.norm(mpcurr[0] - self.vertices[e1plus]): concavedist = True
                elif e2concaveedgefound:
                    if distconv < np.linalg.norm(mpcurr[0] - self.vertices[e2]): concavedist = True
                else:
                    if distconv < radius: concavedist = True

                #print e1,t,e2,concavedist,e1concaveedgefound,e2concaveedgefound



                #print "radius=",radius

                if  mindist < radius and not concavedist:
                    print "distance failure with ", mine,radius,mindist
                    if self.debug3:
                        bz.pe(
                            [[[[self.p(mine, rmine)]], 'g', ''], [[[mpcurr[0]]], 'r', ''], [[[self.p(e1, t)]], 'b', ''],
                             [self.edges]])
                    '''
                    if mine in self.convex_vert:
                        if mine ==e1plus:
                           self.tracededges.append(e1)
                        elif mine==e2minus:
                            self.tracededges.append(e2)
                    '''
                    if not edgechange:
                        mpcurr = self.retract(mplast, mpcurr, mindist, mine)
                    if e1concaveedgefound and e2concaveedgefound:  # case #2A
                        _, r1, mp1 = self.footprint(mine, rmine, e1plus)
                        _, r2, mp2 = self.footprint(mine, rmine, e2minus)
                        self.mps.append([mp2, mine, rmine, e2minus, r2, mpcurr[-1]])
                        # self.mps.append([mp2, mine, rmine, e2, mpcurr[4], mpcurr[-1]])
                        # print "branch to",mine,rmine,e2
                        print "branch to", mine, rmine, e2minus
                        e1 = e1plus;
                        t = r1;
                        e2 = mine;
                        tin = t;
                        rin = -float('inf')
                        mplast[1] = e1;
                        mplast[2] = t;
                        mplast[3] = e2;

                        print "continue with ", e1, t, e2
                        # self.mps.append([mp1, e1plus, r1, mine, rmine, mpcurr[-1]])

                        e1concaveedgefound = False;
                        e2concaveedgefound = False
                        continue
                    elif e2concaveedgefound and not e1concaveedgefound:
                        e1 = mine;
                        t = rmine;
                        e2 = e2minus;
                        tin = t;
                        rin = -float('inf')
                        mplast[1] = e1;
                        mplast[2] = t;
                        mplast[3] = e2

                        print "continue with ", e1, t, e2
                        e2concaveedgefound = False;  # edgechange=True
                        continue
                    elif e1concaveedgefound and not e2concaveedgefound:
                        _, t, mp4 = self.footprint(mine, rmine, e1plus)
                        e1 = e1plus;
                        e2 = mine;
                        tin = t;
                        rin = -float('inf')
                        mplast[1] = e1;
                        mplast[2] = t;
                        mplast[3] = e2

                        print "continue with ", e1, t, e2
                        e1concaveedgefound = False;  # edgechange=True
                        continue

                    else:  # case #2
                        self.mps.append([mpcurr[0], mine, rmine, e2, r, mpcurr[-1]])
                        e2 = mine;
                        rin = -float('inf')
                        mplast[1] = e1;
                        mplast[2] = t;
                        mplast[3] = e2

                        # edgechange=True
                        continue

                elif  concavedist:
                    print "dist failure with concavevert", mini, distconv, radius, e1concaveedgefound, e2concaveedgefound
                    if not edgechange:
                        mpcurr = self.retractConcave(mplast, mpcurr, distconv, mini)
                    if e1concaveedgefound and not e2concaveedgefound:  # case # 3
                        _, r3, mp3 = self.footprint(mpcurr[3], mpcurr[4], mini)
                        _, p0, t0 = self.distance(e1plus, mpcurr[0], True)
                        if mini not in self.concave_vert:
                            self.mps.append([mpcurr[0], mini, r3, mpcurr[3], mpcurr[4], mpcurr[-1]])
                        else:
                            self.mps.append([mpcurr[0], mpcurr[1], mpcurr[2], mpcurr[3], mpcurr[4], mpcurr[-1]])

                        print "branching to ", mini, r3, mpcurr[3]
                        e1 = e1plus;
                        t = t0;
                        e2 = self.prevedge(mini);
                        tin = t;
                        rin = -float('inf')
                        mplast[1] = e1;
                        mplast[2] = t;
                        mplast[3] = e2
                        print "proceeding with", e1, t, e2
                        e1concaveedgefound = False;
                        self.handledconcaveedges.append(mini)
                        continue

                    elif e2concaveedgefound and not e1concaveedgefound:  # case #4
                        _, r3, mp3 = self.footprint(mpcurr[1], mpcurr[2], self.prevedge(mini))
                        radius, p0, t0 = self.distance(mini, mpcurr[0], True)
                        _, r4, mp4 = self.footprint(mini, t0, e2minus)
                        self.mps.append([mpcurr[0], mini, t0, e2minus, r4, radius])
                        print "branching to", mini, t0, e2minus
                        e2 = self.prevedge(mini)
                        mplast[3] = e2
                        rin = -float('inf')
                        e2concaveedgefound = False;
                        self.handledconcaveedges.append(mini)
                        print "continue with", e1, t, e1concaveedgefound, e2concaveedgefound
                        continue
                    elif not e1concaveedgefound and not e2concaveedgefound:  # case#1
                        _, r3, mp3 = self.footprint(e2, r, mini)
                        self.mps.append([mpcurr[0], mini, r3, e2, r, mpcurr[-1]])
                        print "branching to ", mini, r3, e2
                        e2 = self.prevedge(mini)
                        mplast[3] = e2
                        rin = -float('inf')
                        print "continue with ", e1, t, e2
                        self.handledconcaveedges.append(mini)
                        continue

                    elif e1concaveedgefound and e2concaveedgefound:
                        # check e1,e3
                        self.handledconcaveedges.append(mini)
                        dot = np.dot(self.normal(e1, t, True), self.normal(mini, 0.0, True))
                        print dot
                        if dot < 0.01:
                            # _,r3,mp3=self.footprint(e1,t,self.prevedge(mini))
                            radius, fp, t0 = self.distance(mini, mpcurr[0], True)
                            _, r4, mp4 = self.footprint(mini, t0, e2minus)
                            self.mps.append([mpcurr[0], mini, t0, e2minus, r4, radius])
                            print "branching to ", mini, t0, e2minus
                            e2 = self.prevedge(mini);
                            rin = -float('inf')
                            e2concaveedgefound = False
                            mplast[3] = e2

                            # edgechange=True
                            print "continuing with ", e1, t, e2
                            continue
                        else:
                            radius, pi, t0 = self.distance(e1plus, mpcurr[0], True, 0.0)
                            _, r0, mp4 = self.footprint(e1plus, t0, self.prevedge(mini))
                            self.mps.append([mpcurr[0], e1plus, t0, self.prevedge(mini), r0, radius])
                            print "branching to ", e1plus, t0, self.prevedge(mini)
                            radius, pi, t = self.distance(mini, mpcurr[0], True)

                            e1 = mini;
                            tin = t;
                            # edgechange=True;
                            mplast[1] = e1;
                            mplast[2] = t;
                            rin = -float('inf')
                            e1concaveedgefound = False
                            print "continue with", e1, t, e2
                            continue

                    mplast = mpcurr
            #print "adding",mpcurr[0],mpcurr[1],mpcurr[3],mpcurr[2]

    def proceed1(self,ed):
        self.cvcp={};bzs=[];testedges=[]

        for cv1 in self.concave_vert:
            edges = [];
            ed.append(edges)
            mine2, mindist2, rmine2, minc2, mindistc2 = self.findCounterPart(cv1, -1, 59.0, False)
            resultMp = self.resultMp(cv1, -1, 59.0, False)
            if mindistc2 <= mindist2:
                concavedist = True
                print "e1:will branch to concavec", minc2, self.prevedge(cv1)
                _, _, t = self.distance(minc2, resultMp, True)
                e1 = minc2;
                e2 = self.prevedge(cv1)
                e10=e1;deltat=0.01
                while not(e1 != e10 and e1  in self.concave_vert):
                    _,r,mp = self.footprint(e1,t,e2)
                    edges.append(mp)
                    if r<0.0:
                        e2=self.prevedge(e2)
                    t += deltat
                    if t>1.0:
                        e1=self.nextedge(e1)
                if e1 in self.concave_vert:
                    resultMp1 = self.resultMp(e1, -1, 59.0, False)
                    e1 = self.prevedge(e1)
                    _, _, t1 = self.distance(e1, resultMp1, True)
                    #print t,t1
                    while t1!=float('inf') and t<=t1:
                        _,r,mp = self.footprint(e1,t,e2)
                        print e1,t,e2
                        t+=deltat
                        edges.append(mp)

                if len(edges)>2:
                    bzs += fitCurve(edges,10.0)
        bz.pe([[bzs],[self.edges,'g','']])

    def test1(self,edges,ind):
        self.tracededges = []
        self.handledconcaveedges = []
        if ind in self.tracededges or self.prevedge(ind) in self.tracededges: pass
        else:
            mps = [self.vertices[ind], ind, self.its[ind].begin(), self.prevedge(ind), 1.0, 0.0]

            self.mps = [mps]
            print mps
            while len(self.mps) > 0:
                mps = self.mps.pop(-1)
                edge = [];
                edges.append(edge)
                print mps
                self.proceed(mps, edge)

    def iterate(self,edges):
        #take a convex verex- take its previous edge and this edge
        #do distcheck find e1,e3
        #increment e1 and find footprint
        #till you reach a concave_vertex ( either e1 or e3)
        # when you reach a concave vertex
        # create a branch and then proceed

        self.tracededges=[]
        self.handledconcaveedges=[]
        if len(self.concave_vert)>0:
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
             e1nodes = self.edges[ind];e2nodes=self.edges[self.prevedge(ind)]
             lenratio=np.linalg.norm(e1nodes[-1]-e1nodes[0])/np.linalg.norm(e2nodes[-1]-e2nodes[0])
             print lenratio
             if lenratio>5.0 :
                 continue
             if ind in self.tracededges or self.prevedge(ind) in self.tracededges: continue
             l1 = np.linalg.norm(self.edges[ind][-1]-self.edges[ind][0])
             l2 = np.linalg.norm(self.edges[self.prevedge(ind)][-1]-self.edges[self.prevedge(ind)][0])
             if l1/l2 >5.0: continue

             mps = [self.vertices[ind],ind,0.0,self.prevedge(ind),1.0,0.0]

             self.mps=[mps]
             #print mps
             while len(self.mps)>0:
                    mps = self.mps.pop(0)
                    edge=[];edges.append(edge)
                    #if self.its[mps[1]].search(mps[2])==set(): continue
                    print "proceed on",mps[1],mps[2],mps[3]
                    self.proceed3(mps,edge,edges)
            #bz.pe([[edges, 'r', ''], [self.edges, 'g', '']])

        return
    def test(self,pairs):
        deltat=0.01;bzs=[]
        for pair in pairs:
            e1=pair[0];t=0.0;e2=pair[1];pts=[]
            e10=pair[0];e20=pair[1]
            start=True
            while start or e1==e10 or e2==e20 or (not e1 in self.convex_vert and not e2 in self.convex_vert and not e1 in self.concave_vert and not e2 in self.concave_vert):
                start = False
                _,r,mp = self.footprint(e1,t,e2)
                print e1,t,e2,r
                if r<0.0: e2=self.prevedge(e2)
                t=t+ deltat
                if t>1.0:
                    e1=self.nextedge(e1);t=0.0
                if r==float('inf'):e2=self.prevedge(e2);continue; break
                pts.append(mp)
            if len(pts)>2:
                bz1=fitCurve(pts,10.0)
                bzs += bz1
        #print bzs
        bz.pe([[bzs],[self.edges,'g','']])
        return



        pts=[]
        for q in xrange(40):
            t=1.0+q/10.0
            _,r,mp = self.footprint(pair1[0],t,pair1[1])
            if r==float('inf'): break
            pts.append(mp)
        bzs1 = fitCurve(pts,10.0)
        pts=[]
        for q in xrange(20):
            t = 1.0 + q / 10.0
            _, r, mp = self.footprint(pair2[0], t, pair2[1])
            if r == float('inf'): break

            pts.append(mp)
        bzs2 = fitCurve(pts, 10.0)
        curve1 = bezier.Curve(np.array(bzs1[-1]),3)
        curve1=curve1.specialize(0.0,2.0)
        curve2=bezier.Curve(np.array(bzs2[-1]),3)
        bz.pe([[[curve1.nodes]+bzs2,'g',''],[self.edges,'r','']])
        print curve1.intersect(curve2)

        return
    def pe1(self,edges):

        x = [];bzs=[];eds=[]
        for ed in edges:
            #print ed
            if len(ed)==0 :
                continue
            #ed = zip(*ed)[0]
            ed1=[]
            for point in ed:
                #print point[1],point[3]
                if (point[1] in self.convex_vert and point[3] == self.prevedge(point[1])):
                    continue
                if (point[3] in self.convex_vert and point[1] == self.prevedge(point[3])):
                    continue
                ed1.append(point)
            if len(ed1)==0: continue
            ed =zip(*ed1)[0]

            if len(ed)>2:
                bzs += fitCurve(ed,10.0)
            for pt in ed:
                eds.append([pt])
            #x.append([eds, 'r', ''])
        x.append([self.edges, 'g', ''])
        x.append([bzs,'r',''])
        bz.pe(x)
        return bzs,eds
    def pfile(self,ed):
        fout = open("ed.out", "w")
        fout.write(repr(ed))
        fout.close()


class strokes(object):
    def f(self,e1,t,e2):
        r,mp,d=bz.footprint()
    def iterate(self):
        #take pair of convex edges
        # determine a point at which their distance with third edge becomes constant
        return








