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
        dist0, p2i, r1,_ = bz.distance(quad2, p1)  # randomly ..no partcular thought
        edges = []
        edges.append([p1, p2i])
        #print dist0,r1
        '''
        if r1 < 0.0:
            # print "here"
            p2i = bz.p(quad2, 0.0);
            dist0 = np.linalg.norm(p2i - p1)
        elif r1 > 1.0:
            p2i = bz.p(quad2, 1.0);
            dist0 = np.linalg.norm(p2i - p1)
        '''


        for i in xrange(30):
            costheta1 = np.dot(p2i - p1, n1) / dist0
            #print costheta1
            if abs(costheta1) > 0.01:
                dist1 = dist0 / (2 * abs(costheta1))
            else:
                dist1 = dist0 / 2.0
            fpi = p1 + dist1 * n1
            dist2, p2i, r1,_ = bz.distance(quad2, fpi)  # randomly ..no partcular thought
            #print "dist2,",dist2
            if abs(dist2 - dist1) < 0.00001:
                break
            else:
                dist0 = np.linalg.norm(p2i - p1)
        edges.append([p1, fpi]);
        edges.append([p2i, fpi]);
        edges.append([p1]);
        edges += [[p2i], [p1, p2i], [fpi]]

        #print r1,np.cross(fpi-p1,p1prime),abs(dist2-dist1)
        #bz.pe([[edges, 'g', '']])

        if np.cross(fpi - p1, p1prime) >= 0 and abs(dist2 - dist1) < 1e-2:
            return edges, r1, fpi
        else:
            return edges, -1.0, fpi
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
                self.its[i]=it([iv(-1.0,1.0)])
                self.its[self.prevedge(i)]=it([iv(0.0,2.0)])
        return self.convex_vert,self.concave_vert
    def footprint(self,e1,t,e2):
        nodes1 = self.edges[e1];nodes2 = self.edges[e2]
        q=bz.footprint(nodes1,t,nodes2) #r,mp,d

        es=[];
        es.append([self.p(e1,t),q[2]]);
        es.append([self.p(e2,q[1]),q[2]])

       # bz.pe([[[self.edges[e1],self.edges[e2]],'r',''],[es,'g','']])
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

    def distance(self,e,p):
        edge = self.edges[e]
        dist,pi,rn,ep= bz.distance(edge,p)
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
            bz.pe([[self.edges, 'r', ''], [e, 'g', '']])

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
    def distCheck(self,e1,e2):
        es=[]
        for e3,_ in enumerate(self.edges):
            if e3 ==e1 or e3==e2:
                continue
            q=self.diffMp(e1,e2,e3)
            if q[2]!=float('inf'):
                es.append(e3)
        return es
    def nearestEdge(self,p,edges):
        mindist = float('inf');mine = -1
        for e in edges:
            print "calculating distance for",e
            dist,pi,r = self.distance(e,p)
            if dist<mindist:
                mine=e;mindist=dist
        return mine,mindist

    def iterate(self):
        #take a convex verex- take its previous edge and this edge
        #do distcheck find e1,e3
        #increment e1 and find footprint
        #till you reach a concave_vertex ( either e1 or e3)
        # when you reach a concave vertex
        # create a branch and then proceed
        cv = self.convex_vert[0]
        e1=cv;e2=self.prevedge(cv)
        p=self.p(cv,0.01)
        edges = [i for i,_ in enumerate(self.edges) if i!=e1 and i!=e2]
        e3=self.nearestEdge(p,edges)[0]
        print e1,e2,e3
        t0,r,r1=self.diffMp(e1,e2,e3)
        print t0,r,r1

        return
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






