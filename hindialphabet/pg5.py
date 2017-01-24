import numpy as np
import playground3 as pg3
from fitCurves import *
#import svgpathtools as spt
offset = 59
class EdgeTool(object):
    def __init__(self):
        self.c = pg3.Check()
    def ccw(self,A,B,C):
        return (C[1]-A[1])*(B[0]-A[0]) > (B[1]-A[1])*(C[0]-A[0])
    def intersect(self,A,B,C,D):
        return self.ccw(A,C,D) != self.ccw(B,C,D) and self.ccw(A,B,C) != self.ccw(A,B,D)

    def facingEdges(self,edges,offset):
        # take an edge -- generate rays from 1.5*offset  and for each edge
        # check for each edge if it intersects line between first and last points
        pairs = []
        testedges = []
        uncoveredt = []
        for i,edge1 in enumerate(edges):
            pair = []
            unct = []
            testedge = []
            for q in xrange(11):
                t = q/10.0
                if len(edge1) == 2:
                    p = edge1[0] * (1-t) + edge1[1] * t
                    pdash = edge1[1] -edge1[0]
                    pdashnormal = (pdash[1],-1*pdash[0])
                    costheta = pdashnormal[0]* 1/np.linalg.norm(pdashnormal)
                    sintheta = pdashnormal[1]* 1/np.linalg.norm(pdashnormal)
                elif len(edge1) == 3:
                    p1 = edge1[0];pc=edge1[1];p2=edge1[2]
                    #p = edge1[0] * (1-t)**2 + 2*edge1[1]*t*(1-t) + edge1[2] * t**2
                    p = (1-t)**2 * p1 + 2 *t*(1-t)*pc + t**2*p2
                    pdash = 2*(t*((p2-pc)-(pc-p1)) + (pc-p1))
                    pdashnormal = (pdash[1],-1*pdash[0])
                    costheta = pdashnormal[0]* 1/np.linalg.norm(pdashnormal)
                    sintheta = pdashnormal[1]* 1/np.linalg.norm(pdashnormal)
                elif len(edge1) == 4:
                    p1 = edge[0];pc1=edge[1];pc2=edge[2];p2=edge[3]
                    p = (1-t)**3 * p1 + 3*t*(1-t)**2 * pc1 + 3*t**2*(1-t)*pc2 + t**3 * p2
                    pdash = 3*(pc1-p1)*(1-t)**2 + 6*(pc2-pc1)*t*(1-t) + 3*(p2-pc2)*t**2
                    pdashnormal = (pdash[1],-1*pdash[0])
                    costheta = pdashnormal[0]* 1/np.linalg.norm(pdashnormal)
                    sintheta = pdashnormal[1]* 1/np.linalg.norm(pdashnormal)
                pnew = p + np.array([costheta,sintheta]) * 2.1 * offset
                print [p,pnew]
                testedge.append([p,pnew])
                #print pnew
                #prev = -1
                uncovered = True
                for j in xrange(len(edges)):
                    if j == i:
                        continue
                    edge2 = edges[j]

                    if self.intersect(p,pnew,edge2[0],edge2[-1]) and not np.array_equal(p,edge2[0]) and not np.array_equal(p,edge2[-1]):
                        testedges.append([pnew,p] + [edge2[0],edge2[-1]])
                        uncovered = False
                        if j not in pair:
                            pair.append(j)
                if uncovered:
                        unct.append(t)

            testedges += [testedge]
            pairs.append(pair)
            uncoveredt.append(unct)
        return pairs,testedges,uncoveredt
    def offset(self,edges,offset):
        offsetedges = []
        fe,_,uct = self.facingEdges(edges,offset) #uncoveredt
        indices = sorted(range(len(fe)), reverse=True,key=lambda k: len(fe[k]))
        done = []
        for _,i in enumerate(indices):
            if i in done:
                continue
            for j in fe[i]:
                p,_ = self.offsetNPEdges(edges[i],edges[j],uct[i],uct[j],offset)
                if len(p) >1:
                    offsetedges.append(p)
                done.append(j)
        oe = []
        for e in offsetedges:
            for p in e:
                oe.append(p)

        return oe







    def offsetNPEdges(self,edge1,edge2,uct1,uct2,offset): # offset non parallel edges
       # for 100 points on edge1 -
       #   find tangent and normal to that point
       #          draw a normal from a point of offset -5 to edge2 and see if it is equal to offset -2 and go upto offset +5
       #               to find that point
       points = []
       testpoints = []
       for q in xrange(31):
           t = q/30.0
           if len(edge1) == 2:
               p = edge1[0] * (1-t) + edge1[1] * t
               pdash = edge1[1] -edge1[0]
               pdashnormal = (pdash[1],-1*pdash[0])
               costheta = pdashnormal[0]* 1/np.linalg.norm(pdashnormal)
               sintheta = pdashnormal[1]* 1/np.linalg.norm(pdashnormal)
           elif len(edge1) == 3:
               p1 = edge1[0];pc=edge1[1];p2=edge1[2]
               #p = edge1[0] * (1-t)**2 + 2*edge1[1]*t*(1-t) + edge1[2] * t**2
               p = (1-t)**2 * p1 + 2 *t*(1-t)*pc + t**2*p2
               pdash = 2*(t*((p2-pc)-(pc-p1)) + (pc-p1))
               pdashnormal = (pdash[1],-1*pdash[0])
               costheta = pdashnormal[0]* 1/np.linalg.norm(pdashnormal)
               sintheta = pdashnormal[1]* 1/np.linalg.norm(pdashnormal)
           elif len(edge1) == 4:
               p1 = edge[0];pc1=edge[1];pc2=edge[2];p2=edge[3]
               p = (1-t)**3 * p1 + 3*t*(1-t)**2 * pc1 + 3*t**2*(1-t)*pc2 + t**3 * p2
               pdash = 3*(pc1-p1)*(1-t)**2 + 6*(pc2-pc1)*t*(1-t) + 3*(p2-pc2)*t**2
               pdashnormal = (pdash[1],-1*pdash[0])
               costheta = pdashnormal[0]* 1/np.linalg.norm(pdashnormal)
               sintheta = pdashnormal[1]* 1/np.linalg.norm(pdashnormal)
           diff = 50.0
           for i in xrange(-5,5):
                if diff < 5:
                   break
                pnew = p + np.array([costheta,sintheta]) * (offset+i)
                testpoints.append([p,pnew])
                dist,tp = self.normalFromPoint(edge2,pnew)
                testpoints += tp
                if dist == -1:
                    continue
                diff = abs(offset +i - dist)
                #print diff
           if diff < 50.0:
               points.append(pnew)
       #print uct1,uct2
       bezier1 = self.c.offsetEdge(edge1,offset,uct1)
       bezier2 = self.c.offsetEdge(edge1,offset,uct1)
       #print bezier1,bezier2,len(bezier1),len(bezier2)


       if len(points) >1:
           beziers = fitCurve(points,10.0)
           return beziers + bezier1 + bezier2,testpoints
       else:
           return bezier1 + bezier2 ,testpoints


    def normalFromPoint(self,edge,point):
        found = False
        testpoints = []
        for q in xrange(31):
            t = q/30.0
            if len(edge) == 3:
                p1 = edge[0];pc=edge[1];p2=edge[2]
                p = (1-t)**2 * p1 + 2 *t*(1-t)*pc + t**2*p2
                pdash = 2*(t*((p2-pc)-(pc-p1)) + (pc-p1))
            elif len(edge) == 2:
                p = edge[0] * (1-t) + edge[1] * t
                pdash = edge[1] -edge[0]
            elif len(edge) == 4:
                p1 = edge[0];pc1=edge[1];pc2=edge[2];p2=edge[3]
                p = (1-t)**3 * p1 + 3*t*(1-t)**2 * pc1 + 3*t**2*(1-t)*pc2 + t**3 * p2
                pdash = 3*(pc1-p1)*(1-t)**2 + 6*(pc2-pc1)*t*(1-t) + 3*(p2-pc2)*t**2
            normpdash = pdash/np.linalg.norm(pdash)
            normdiff = (p-point)/np.linalg.norm(p-point)
            testpoints.append([p,point])
            #print t,abs(np.dot(normdiff,normpdash))
            if abs(np.dot(normdiff,normpdash)) < 0.1:
                found = True
                return np.linalg.norm(point-p),testpoints
        if not found:
            return -1,testpoints
c= pg3.Check()
edges = c.fontOutline('dn_ra')
slopes = []
for edge in edges:
    slope = []
    p1 = edge[1]-edge[0]
    slope.append(180/np.pi* np.arctan2(p1[1],p1[0]))

    if len(edge) > 2:
        p2 = edge[2] - edge[1]
    else:
        p2 = p1
    slope.append(180/np.pi* np.arctan2(p2[1],p2[0]))
    slopes.append(slope)
