import numpy as np
import random
import copy
from fitCurves import fitCurve
import playground3 as pg3
import sys
check = pg3.Check()

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
           return float("inf"),[],-1


    def footprintReflex(reflex,quad):
        def f(p1,p1prime,p2,p2prime):
           return np.dot(p1prime,p1prime)*(np.dot(p1-p2,p2prime))**2 - np.dot(p2prime,p2prime)*(np.dot(p1-p2,p1prime))**2

        def fprime(p1,p1prime,p2,p2prime,p2dprime):
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
        def solve(p1,p1prime,normal1,p2,p2prime,normal2,p2dprime):
            denominator = np.dot(normal1,p2prime)
            val1 = np.dot(p1-p2,p1prime)
            val2 = np.dot(p1-p2,p2prime)
            if abs(denominator) < 0.01: # parallel normals
                if abs(val1) <0.01 and abs(val2) < 0.01: # overlapping normals
                    #print "overlapping normal"
                    return r,True
                else:
                    #print "parallel normals- get out of this",val1,val2
                    return random.uniform(0,1),False #a random number between 0 and 1
                #print r
            else: # not parallel normal
                val = f(p1,p1prime,p2,p2prime)
                #print val
                if  abs(val) < 1:
                    #print "reached destination"
                    return r,True
                else:
                    #print "next iteration"
                    r = r - val/fprime(p1,p1prime,normal1,p2,p2prime,normal2,p2dprime)
                    if r>1 or r<0:
                        return random.uniform(0,1),False
                    else:
                        return r,False
        p1 = reflex[0] #
    def insertReflex(self,edges,vertices,cnts,concave_vert,lastts):#
        concave_vert = copy.copy(concave_vert)
        concave_vert.sort()
        #print cnts
        cnts1 = copy.copy(cnts)
        for k,_ in enumerate(edges):
            if not k in lastts:
                lastts[k] = [0.0,1.0]

        while len(concave_vert) >0:
            print concave_vert
            i = concave_vert.pop(-1)
            print i
            e2 = i
            e1 = self.prevedge(i,cnts)
            e1prime = self.p(edges[e1],0.999)
            e2prime = self.p(edges[e2],0.001)
            quad = []
            quad.append(e1prime)
            print e1prime
            quad.append(vertices[i])
            quad.append(e2prime)
            print e1,e2
            lastts[e1][1] = 0.999
            edges[e1][2] = e1prime
            edges[e2][0] = e2prime
            vertices.insert(i,e1prime)
            vertices[i+1] = e2prime
            lastts[i+1] = [0.001,1.0]

            edges.insert(i,quad)
            ind1 = 0
            ind3 = 0

            for ind2,cnt in enumerate(cnts):
                for ind,k in enumerate(cnt):
                    if k > i:
                        cnt[ind] = k+1
                    elif k ==i:
                        ind1 = ind
                        ind3 = ind2

            cnts[ind3].insert(ind1+1,i+1)
            for k,cv in enumerate(concave_vert):
                if cv > i:
                    concave_vert[k] = cv + 1
            for k,_ in enumerate(edges):
                if not k in lastts:
                    lastts[k] = [0.0,1.0]

        print cnts











    def footprint(self,quad1,t,quad2,initialr=-1):
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
        print r
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



            return edges,r,p1 + alpha*normal1
        else:
            return edges,-1,None
    def bisector(self,quad1,quad2):
        finedges = []
        points = []
        for q in xrange(101):
            t = q/100.0
            edges,r,point = self.footprint(quad1,t,quad2)
            if r!=-1:
                #print t,r

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

    def medialaxis(self,edges1,vertices1,vert_edges1,cnts1,stepsize=0.05):

        cnts = copy.copy(cnts1)
        edges = copy.copy(edges1)
        vertices = copy.copy(vertices1)
        vert_edges = copy.copy(vert_edges1)
        convex_vert,concave_vert = self.vertAnalyse(vertices,edges,cnts)
        print convex_vert,concave_vert
        lastts = {} # store as [float,float] initial,final
        for k,_ in enumerate(edges):
            lastts[k] = [0.0,1.0]

        self.insertReflex(edges,vertices,cnts,concave_vert,lastts)
        for k,_ in enumerate(edges):
            if not k in lastts:
                lastts[k] = [0.0,1.0]

        #Sreturn edges
        convex_vert,concave_vert = self.vertAnalyse(vertices,edges,cnts)



        maxis = []
        simpleedge =[]
        count = 0
        edge_count = len(edges) -1
        deltat = stepsize
        alledges = []
        tracededges = []

        def distance(ind,p0,lastts,edges):#distance between a point and p
           quad = edges[ind]
           pdprime = 2*(quad[2]+quad[0]-2*quad[1])
           t = lastts[ind][0]
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
                    if t<lastts[ind][0] or t>lastts[ind][1]:
                        t = random.uniform(lastts[ind][0],lastts[ind][1])
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
        def retract(simpleedge,mpcurrent,mpdistl,e3,dist3,edges,lastts):
               e1 = mpdistl[1]
               e2 = mpdistl[3]
               rf = mpdistl[4]
               tf = mpdistl[2]
               rl = mpcurrent[4]
               tl = mpcurrent[2]
               dist = mpcurrent[-1]
               #print mpdistl
               #print mpcurrent
               #print "before retract distance from edges e1,e2,e3",dist,dist3,tf,tl
               #print simpleedge
               simpleedge.pop(-1)
               count = 0
               if tf ==tl:
                   tf = tl -0.02
               rn = rl
               tn = tl
               mpn = mpcurrent[0]


               while abs(dist3-dist)>1.0 and count < 10:
                    count += 1
                    tn = (tf + tl)/2.0
                    #print tf,tn,tl,rl
                    _,rn,mpn = self.footprint(edges[e1],tn,edges[e2])
                    #print "rn =",rn
                    if rn == -1: # it should not happen
                        print "possible bug in retract or footprint "
                        break
                    else:
                        dist3,pind,tind = distance(e3,mpn,lastts,edges)

                        dist = np.linalg.norm(mpn-self.p(edges[e1],tn))
                        #diste2 = np.linalg.norm(mpn-self.p(edges[e2],rn))
                        #print "diste1,diste2",diste1,diste2
                        #dist = max(diste1,diste2)

                        if dist3 < dist:
                            tl = tn
                        else:
                            tf = tn

               if rn!=-1:
                    mpfn = [mpn,e1,tn,e2,rn]
                    simpleedge.append(mpfn)
                    lastts[e1][0] = tn
                    lastts[e2][1] = rn
                    #print "after retract distance from edges e1,e2,e3",dist,dist3
                    #print simpleedge





        def hp(e1,e2,simpleedge,maxis,mpl):
            print "handling ",e1,e2
            breakpoint = False
            t = lastts[e1][0] # go forward
            print 'initial t',t
            u = lastts[e2][1] #go backward
            nextpairs = []
            mp0 = mpl[0]
            mpdistl = copy.copy(mpl)
            #print mpl
            if t>=lastts[e1][1] or u <=lastts[e2][0]:
                return
            while  breakpoint == False:
                e1false = True
                e2false = True

                if lastts[e1][1]<=lastts[e1][0]: # e1 exhausted
                    e4= self.nextedge(e1,cnts)
                    if e4==e1:
                        break
                    else:
                        e1 = e4
                else:
                    _,r,mp1 = self.footprint(edges[e1],t,edges[e2],lastts[e2][1])
                    print e1,t,e2,r,lastts[e1],lastts[e2]


                    if r!=-1 and r<lastts[e2][1]:
                        t = min(lastts[e1][1],t + deltat)

                        lastts[e1][0] = t
                        lastts[e2][1] = r

                        mp0 = mp1

                        simpleedge.append(mpl)
                        dist = np.linalg.norm(mp0-self.p(edges[e2],r))
                        mpl=[mp0,e1,t,e2,r,dist]

                        print '#',e1,e2
                        e2false = False
                        #count += 1
                if e2==7 and e1==8:
                    print lastts[e2],t,r
                if lastts[e2][1]<=lastts[e2][0]:
                    e4 = self.prevedge(e2,cnts)
                    if e2 == e4:
                        break
                    else:
                        e2 = e4
                else:
                    u = lastts[e2][1]
                    _,r,mp1 = self.footprint(edges[e2],u,edges[e1],lastts[e1][0])
                    print e1,r,e2,u,lastts[e1],lastts[e2]

                    if r!=-1 and r<lastts[e1][1]:
                        u = max(lastts[e2][0],u-deltat)
                        lastts[e1][0] = r
                        lastts[e2][1] = u
                        mp0 = mp1
                        dist = np.linalg.norm(mp0-self.p(edges[e1],r))
                        mpl=[mp0,e1,r,e2,u,dist]
                        simpleedge.append(mpl)


                        e1false = False
                        #count += 1

                if e1false and e2false:
                    break
                noofsteps = 0
                possibleedges = []
                dist = mpl[-1]
                for ind,_ in enumerate(edges):
                    if ind == e1 or ind == e2:
                        continue

                    dist1,pind,tind = distance(ind,mpl[0],lastts,edges)

                    if ind==-1:
                        print dist,dist1,ind,mpl[0],lastts[6]
                    if dist1 <= dist: # found a third edge
                         #TODO retract mp and delete additional points
                         #retract(simpleedge,mp,)
                         retract(simpleedge,mpl,mpdistl,ind,dist1,edges,lastts)
                         #return
                         possibleedges.append(ind)
                         #print dist,dist1,ind,tind
                         break
                mpdistl = copy.copy(mpl)
                if len(possibleedges)>0:

                    breakpoint = True
                    maxis.append(simpleedge)
                    simpleedge = []
                    for ind in possibleedges:
                        nextpairs.append([ind,e1,copy.copy(mpl)])
                        nextpairs.append([ind,e2,copy.copy(mpl)])
            for pair in nextpairs:
                hp(pair[0],pair[1],simpleedge,maxis,pair[2])
            return











        def handlePair(e1,e2,simpleedge,maxis,mp):
            print "handling edges ",e1,e2
            if e1 == e2:
                return
            nexte2 = self.prevedge(e2,cnts)
            nexte1 = self.nextedge(e1,cnts)
            t = lastts[e1][0]
            tfinal = lastts[e1][1]
            if t>=tfinal:
                return
            ri = lastts[e2][0]
            rf = lastts[e2][1]
            if ri>=rf:
                return
            print t
            count = 0
            noofsteps = 0
            maxnoofsteps = 3 # after every so many steps we would check distance check
            while t < lastts[e1][1] and t<=1.0 and t>=0.0 and lastts[e2][1]>lastts[e2][0]:
                _,r,mp1 = self.footprint(edges[e1],t,edges[e2])
                print e1,t,e2,r,lastts[e1],lastts[e2],nexte1,nexte2

                if r!=-1:
                    mp = mp1
                    simpleedge.append(mp)
                    print e1,e2,count
                    count += 1
                    noofsteps += 1
                    dist = np.linalg.norm(mp-self.p(edges[e1],t))
                    t =min (t+deltat,1.0)
                    lastts[e1][0] = t
                    lastts[e2][1] = r
                else:
                    # vertex of edge 2 has been reached
                    r = max(lastts[e2][0],lastts[e2][1]-deltat)
                    _,t,mp1 = self.footprint(edges[e2],r,edges[e1])
                    print e1,t,e2,r,lastts[e1],lastts[e2],nexte1,nexte2

                    if t!=-1:
                        lastts[e1][0] =t
                        simpleedge.append(mp)
                        mp =mp1
                        count += 1
                        noofsteps += 1
                        dist = np.linalg.norm(mp-self.p(edges[e1],t))
                        t = min(1.0,t+deltat)
                        lastts[e1][0] = t
                        lastts[e2][1] = r

                    if lastts[e1][0] >=lastts[e1][1]:
                        handlePair(nexte1,e2,simpleedge,maxis,mp)
                        return
                    elif lastts[e2][1]<=lastts[e2][0]:
                        handlePair(e1,nexte2,simpleedge,maxis,mp)
                        return
                    else:
                        handlePair(nexte1,nexte2,simpleedge,maxis,mp)
                        return



                distcheck = True
                if e2==e1+1 and e2 in convex_vert: # if both share a convex edge do not do distance check
                    distcheck = False
                if True:
                #noofsteps >= maxnoofsteps and distcheck: # check distance criterion
                    noofsteps = 0
                    possibleedges = []
                    for ind,_ in enumerate(edges):
                        if ind == e1 or ind == e2:
                            continue

                        dist1,_,_ = distance(ind,mp,lastts,edges)
                        #print dist,dist1,ind,mp
                        if dist1 < dist: # found a third edge
                             #TODO retract mp and delete additional points
                             possibleedges.append(ind)
                             print ind,e1," after distance"
                             print ind,e2," after distance"


                    maxis.append(simpleedge)
                    simpleedge = []
                    if len(possibleedges)>0:
                        for _,ind in enumerate(possibleedges):
                             handlePair(ind,e2,simpleedge,maxis,mp)
                             handlePair(e1,ind,simpleedge,maxis,mp)
                        return

            #e1 exhausted
            lastts[e1] = [1.1,1.0]
            maxis.append(simpleedge)
            simpleedge = []
            if lastts[nexte1][0]<lastts[nexte1][1]:
                handlePair(nexte1,e2,simpleedge,maxis,mp)
                return
                    #reflexedges = {}
        #for i in concave_vert:
            #reflexedges[i] = reflexVertex(i)

        print cnts
        #handlePair(0,17,simpleedge,maxis,vertices[0])

        for cvert in convex_vert:
            firstvertex = cvert

            secondedge = self.prevedge(firstvertex,cnts)
            firstedge = firstvertex
            #lastts[firstedge][0] = 0.0

            mp = [vertices[firstvertex],firstedge,0.0,secondedge,1.0]

            #print firstedge,secondedge
            nextpairs = hp(firstedge,secondedge,simpleedge,maxis,mp)
            print nextpairs


        '''
        for k,_ in enumerate(edges):
            for l,_ in enumerate(edges):
                if lastts[k][0]<lastts[k][1] and lastts[l][0]<lastts[l][1]:
                    handlePair(k,l,simpleedge,maxis,mp)
        '''

        return maxis,edges
    def printMa(self,ma):
        points = []
        for k,ed in enumerate(ma):
            for point in ed:
                if point != None:
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
         R = math.power(np.linalg.norm(pprime),3)/np.linalg.norm(np.dot(pprime,pdprimen))
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
        return mod
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
