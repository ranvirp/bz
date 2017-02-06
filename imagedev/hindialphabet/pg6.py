import numpy as np
from numpy import array
import playground3 as pg3
import manew
import copy
reload(manew)
c = pg3.Check()
f1 = open('output')
q=f1.read().rstrip('\n');q=eval(q)
fka =manew.fp('dn_ka')
#c.printEdges(fka.edges)
def pp():
 def offset(e1,t,offset):
     #print e1,t,offset
     #n = np.array(fka.normal(fka.edges[e1],t,True))
     pprime = fka.pprime(fka.edges[e1],t)
     norm = np.linalg.norm(pprime)
     normal = np.array([pprime[1],-pprime[0]])
     n = normal/norm
     p = fka.p(fka.edges[e1],t)
     return p + offset*n
 q2=fka.printMa(q)
 rad=60.0
 rads={}
 #c.printEdges(fka.edges)
 for j,ed in enumerate(q):
   preve2='';preve1='';rad=61.0
   for i,pt in enumerate(ed):
        if pt[3] in fka.reflex_edges and pt[1] not in fka.reflex_edges:
           if pt[3]!=preve2 and pt[1]==preve1:
               preve2=pt[3];
               rad = pt[-1];preve1=pt[1]
               print "rad for e1,e2=",pt[1],pt[3],rad,j
           else:
               pt[0] = offset(pt[1],pt[2],rad);preve2 = pt[3];preve1=pt[1]
        else:
               preve2 = pt[3];preve1=pt[1]
   for i in xrange(len(ed)-1,-1,-1):
       pt = ed[i]
       if pt[1] in fka.reflex_edges and pt[3] not in fka.reflex_edges:
          if pt[1]!=preve1 and pt[3]==preve2:
              preve2=pt[3];preve1=pt[1]
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
         if np.linalg.norm(ed[-1][0]-q[i+1][0][0])<3.0:
             ed[-1][0]=q[i+1][0][0]
 q1 = [q[i] for i,ed in enumerate(q) if i not in pointindices]

 q1 =fka.printMa(q1);
 c.printEdges(q1[1]+fka.edges,0,False)
 c.printEdges(q1[1])
 c.printEdges(q2[1]+fka.edges,0,False)


def test():
    ep = []
    for q in xrange(int(r2*100),101,1):
        t = q/100.0
        p = ft.p(e[7],t)
        pprime = ft.pprime(e[7],t)
        n = [pprime[1],-1*pprime[0]]
        mp  = p + n/np.linalg.norm(n)*rad2
        #ep.append([mp])
        n0 = mp-v[5]
        t1 = np.linalg.norm(n0-n01)/np.linalg.norm(n02-n01)
        print t1
        theta = np.arctan2(n0[1],n0[0])
        #print "theta=",theta*180/np.pi,(theta2+theta1)*180/np.pi*1/2.0


        if theta<(theta1+theta2)/2 or np.linalg.norm(mp-v[2])<=rad2:
            break
    for q in xrange(101,51,-1):
        t = q/100.0
        #t = 1-t
        p = v[5]+n01*(1-t)+n02*t
        ep.append([[p],'b'])
        #print "theta=",theta*180/np.pi,(theta2+theta1)*180/np.pi*1/2.0


    c.printEdges(ep+e+[[v[5]]],0,False)
def test1(e1,e2,t,v1,e0):
    _,r2,mp2 = ft.footprint(e[v1],0.0,e[e1])
    n02 = mp2-v[v1]
    rad2 = np.linalg.norm(n02)
    eprev = ft.prevedge(v1)
    print "eprev=",eprev

    _,r1,mp1 = ft.footprint(e[eprev],1.0,e[e2])
    print "r1 =",r1
    n01 = mp1-v[v1]
    rad1 = np.linalg.norm(n01)
    #theta2 = np.arctan2(n02[1],n02[0])
    #theta1 = np.arctan2(n01[1],n01[0])
    ep=[]
    #theta = theta1*(1-t)+theta2*t

    #vectr = np.array([np.cos(theta),np.sin(theta)])
    for q in xrange(101,51,-1):
        t0 = q/100.0
        #t = 1-t
        p = v[v1]+n01*(1-t0)+n02*t0
        ep.append([[p],'b'])

    vectr = n01*(1-t)+n02*t
    print "t=",t

    if t>=0.5:
        rad = rad2
        r0=r2
    else:
        rad = rad1
        r0=r1
        print "r0,r1",r1,r0
    ep.append([v[v1],v[v1]+vectr])
    print r0, rad,rad1,rad2

    def f(r):
        p = ft.p(e[e0],r)
        pprime = ft.pprime(e[e0],r)
        norm1 =np.linalg.norm(pprime)
        n = np.array([pprime[1]/norm1,-1*pprime[0]/norm1])

        #print p,pprime,n
        mp = p + rad*n
        ep.append([v[v1],mp])
        d1 = mp-v[v1]
        return np.cross(d1,vectr)

    rs=[]
    rs.append(r0)
    rs.append(1.0)
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
        print rn
        rs.append(rn)
        fximinus1 = fxi
        fxi = f(rn)
    if rn>=0 and rn<=1:
        p = ft.p(e[e0],rn)
        n = ft.normal(e[e0],rn,True)
        mp = p + rad*n
        ep.append([[v[v1],mp],'r'])
        ep.append([[p,mp],'g'])
    else:
       print "failed"
    c.printEdges(ep+ft.edges,0,False)
