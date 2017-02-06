import numpy as np
import playground3 as pg3
import matplotlib.pyplot as plt
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
    def facingEdges(self,edges,offset,noofpoints=100):
        pairs = []
        testedges = []
        uncoveredt = []
        for i,edge1 in enumerate(edges):
          uncovered = []
          pair = []
          for q in xrange(noofpoints+1):
              t = q/(noofpoints * 1.0)
              if len(edge1) == 2:
                  p = edge1[0] * (1-t) + edge1[1] * t
              elif len(edge1) == 3:
                  p1 = edge1[0];pc=edge1[1];p2=edge1[2]
                  #p = edge1[0] * (1-t)**2 + 2*edge1[1]*t*(1-t) + edge1[2] * t**2
                  p = (1-t)**2 * p1 + 2 *t*(1-t)*pc + t**2*p2
              elif len(edge1) == 4:
                  p1 = edge[0];pc1=edge[1];pc2=edge[2];p2=edge[3]
                  p = (1-t)**3 * p1 + 3*t*(1-t)**2 * pc1 + 3*t**2*(1-t)*pc2 + t**3 * p2
              for j in xrange(len(edges)):
                    if j == i:
                        continue
                    edge2 = edges[j]

                    #if self.intersect(p,pnew,edge2[0],edge2[-1]) and not np.array_equal(p,edge2[0]) and not np.array_equal(p,edge2[-1]):
                    dist,minslope,savedp,_ = self.normalFromPoint(edge2,p)
                    if dist!=-1 and  dist <= 2.5 * offset and minslope <0.5:
                        if j not in pair:
                            pair.append(j)
                    else:
                        uncovered.append(t)
          uncoveredt.append(uncovered)
          pairs.append(pair)
        return pairs,uncoveredt
    def offsetEdges(self,edges):
        curves = []
        for k,edge1 in enumerate(edges):
            for p in range(len(edges)-1,0,-1):
                if k==p:
                    continue
                curve= self.bisectorQuad1(edge1,edges[p])
                if len(curve) >0 :
                    curves += curve
        return curves





    def facingEdgesOld(self,edges,offset):
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
                pnew = p + np.array([costheta,sintheta]) * 2.2 * offset
                #print [p,pnew]
                #testedges.append([p,pnew])
                #print pnew
                #prev = -1
                uncovered = True
                for j in xrange(len(edges)):
                    if j == i:
                        continue
                    edge2 = edges[j]

                    #if self.intersect(p,pnew,edge2[0],edge2[-1]) and not np.array_equal(p,edge2[0]) and not np.array_equal(p,edge2[-1]):
                    dist,minslope,savedp,_ = self.normalFromPoint(edge2,pnew)
                    testedges.append([[p,savedp],'blue'])
                    vect = (savedp - p)/np.linalg.norm(savedp - p)
                    vectnorm = np.array([vect[1],-1*vect[0]])
                    vect1 = (pnew-p)/np.linalg.norm(pnew-p)

                    #print dist,np.dot(vect1,vectnorm),np.linalg.norm(p-pnew),np.linalg.norm(p-savedp)
                    if dist>0 and  np.linalg.norm(p-pnew)> np.linalg.norm(p-savedp) and abs(np.dot(vect1,vectnorm)) <0.1:
                        #testedges.append([pnew,p])
                    #if dist <
                        uncovered = False
                        if j not in pair:
                            pair.append(j)
                if uncovered:
                        unct.append(t)

            pairs.append(pair)
            uncoveredt.append(unct)
        return pairs,testedges,uncoveredt
    def offset(self,edges,offset,fe=None,uct=None):
        offsetedges = []
        if fe == None:
           fe,uct = self.facingEdges(edges,offset) #uncoveredt
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
    def bisectorQuad1(self,edge1,edge2):
        if len(edge1)!=3:
            return [[]]
        p11,pdash11,pnormal11,costheta11,sintheta11 = self.curveValues(edge1,0)
        p12,pdash12,pnormal12,costheta12,sintheta12 = self.curveValues(edge1,1)
        p21,pdash21,pnormal21,costheta21,sintheta21 = self.curveValues(edge2,0)
        p22,pdash22,pnormal22,costheta22,sintheta22 = self.curveValues(edge2,1)


        points = []
        for q in xrange(101):
            t = q/100.0
            p1,pdash1,pnormal1,costheta1,sintheta1 = self.curveValues(edge1,t)
            #r = (np.dot(pdash11-pdash21,pdash22-pdash21) + t* np.dot(pdash12-pdash11,pdash22-pdash21))/(np.linalg.norm(pdash22-pdash21))**2
            r = np.dot(pnormal1,-1*pdash21)/np.dot(pnormal1,pdash22-pdash21)
            #print r
            if r>1 or r<0:
                continue

            p2,pdash2,pnormal2,costheta2,sintheta2 = self.curveValues(edge2,r)
            #print abs(np.dot(pdash2,p2-p1)/(np.linalg.norm(p2-p1)*np.linalg.norm(pdash2))), abs(np.dot(pdash2,p2-p1)/(np.linalg.norm(p2-p1)*np.linalg.norm(pdash1)))
            if not(abs(np.dot(pdash2,p2-p1)/(np.linalg.norm(p2-p1)*np.linalg.norm(pdash2) ))<0.2 and abs(np.dot(pdash2,p2-p1)/(np.linalg.norm(p2-p1)*np.linalg.norm(pdash1))) <0.2):
                continue
            offset = np.linalg.norm(p1-p2)/2
            print offset
            points.append(p1+np.array([costheta1,sintheta1])*offset)
            #print offset
        if len(points) > 2:
            return fitCurve(points,10.0)
        elif len(points) ==2:
            return [points]
        else:
            return []





    def bisectorQuad(self,edge1,edge2):
        #Method 1:
        #costheta = normal.p1-p2 /normal.p1-p2
        #norm(p1-p2)*norm(p1-p2).norm(p1)/2*(normal.p1-p2)
        #offset each point
        points = []
        for q in xrange(101):
            t = q/100.0
            p1,pdash1,pnormal1,costheta1,sintheta1 = self.curveValues(edge1,t)
            for m in xrange(101):
                r = 1 - m/100.0
                p2,pdash2,pnormal2,costheta2,sintheta2 = self.curveValues(edge2,r)
                cosangle1 = np.dot(p2-p1,pnormal1)/(np.linalg.norm(p1-p2)*np.linalg.norm(pnormal1))
                cosangle2 = np.dot(p1-p2,pnormal2)/(np.linalg.norm(p1-p2)*np.linalg.norm(pnormal2))
                print cosangle2,cosangle1
                if abs(cosangle2 - cosangle1) > 0.1:
                    continue
                offset = ((np.linalg.norm(p1-p2))**2 * np.linalg.norm(pnormal1))/(2*np.dot(pnormal1,p2-p1))
                print offset
                if offset < 0:
                    continue
                else:
                    points.append(p1 + np.array([costheta1,sintheta1]) * offset)
                    break
        if len(points)>2:
            beziers = fitCurve(points,10.0)
            return beziers
        else:
            return [points]
    @classmethod
    def curveValues(self,edge1,t):
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
        return p,pdash,pdashnormal,costheta,sintheta

    def offsetNPEdges(self,edge1,edge2,uct1,uct2,offset): # offset non parallel edges
       # for 100 points on edge1 -
       #   find tangent and normal to that point
       #          draw a normal from a point of offset -5 to edge2 and see if it is equal to offset -2 and go upto offset +5
       #               to find that point
       points = []
       testpoints = []
       for q in xrange(101):
           t = q/100.0
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
                dist,_,_,tp = self.normalFromPoint(edge2,pnew)
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


    def fillPoints(self,edges,noofpoints=100):
        outline = []
        for edge1 in edges:
            for q in xrange(noofpoints+1):
                t = q/(noofpoints * 1.0)
                if len(edge1) == 2:
                    p = edge1[0] * (1-t) + edge1[1] * t
                elif len(edge1) == 3:
                    p1 = edge1[0];pc=edge1[1];p2=edge1[2]
                    #p = edge1[0] * (1-t)**2 + 2*edge1[1]*t*(1-t) + edge1[2] * t**2
                    p = (1-t)**2 * p1 + 2 *t*(1-t)*pc + t**2*p2
                elif len(edge1) == 4:
                    p1 = edge[0];pc1=edge[1];pc2=edge[2];p2=edge[3]
                    p = (1-t)**3 * p1 + 3*t*(1-t)**2 * pc1 + 3*t**2*(1-t)*pc2 + t**3 * p2
                outline.append(p)
        return outline

    def normalFromPoint(self,edge,point):
        found = False
        testpoints = []
        minslope = 100.0
        savedp = None
        dist = -1
        testpoints.append([[point],'r'])
        for q in xrange(5):
            t = q/4.0
            if len(edge) == 3:
                p1 = edge[0];pc=edge[1];p2=edge[2]
                p = (1-t)**2 * p1 + 2 *t*(1-t)*pc + (t**2)*p2
                pdash = 2*(t*((p2-pc)-(pc-p1)) + (pc-p1))
            elif len(edge) == 2:
                p = edge[0] * (1-t) + edge[1] * t
                pdash = edge[1] -edge[0]
            elif len(edge) == 4:
                p1 = edge[0];pc1=edge[1];pc2=edge[2];p2=edge[3]
                p = ((1-t)**3) * p1 + 3*t*(1-t)**2 * pc1 + 3*t**2*(1-t)*pc2 + (t**3) * p2
                pdash = 3*(pc1-p1)*(1-t)**2 + 6*(pc2-pc1)*t*(1-t) + 3*(p2-pc2)*t**2
            normpdash = pdash/np.linalg.norm(pdash)
            normdiff = (p-point)/np.linalg.norm(p-point)
            testpoints.append([[p,point],'g'])
            #print p,point
            #print t,abs(np.dot(normdiff,normpdash)),np.linalg.norm(p-point),p
            if abs(np.dot(normdiff,normpdash)) < minslope:
                minslope = abs(np.dot(normdiff,normpdash))
            #    savedp = p
                dist= np.linalg.norm(point-p)
                #print p


            if abs(np.dot(normdiff,normpdash)) < 0.3:
                found = True
                dist= np.linalg.norm(point-p)
                savedp = p
                #print p

        return dist,minslope,savedp,testpoints
    def sfv(vor): # skeleton from voronoi
        for vpair in vor.ridge_vertices:
            if vpair[0] >= 0 and vpair[1] >= 0:
                v0 = vor.vertices[vpair[0]]
                v1 = vor.vertices[vpair[1]]
                # Draw a line from v0 to v1.
                plt.plot([v0[0], v1[0]], [v0[1], v1[1]], 'k', linewidth=2)
        plt.show(block=False)

plt.show()
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
