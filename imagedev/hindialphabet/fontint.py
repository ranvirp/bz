import numpy as np
from matplotlib import pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

def quadCurve(p0,m0,pg1,m1,pg2,m2):
    A = [[2*m1,m1,2,1],[0,m2,0,1],[m0,0,1,0],[-1*m2,m2,-1,1]]
    try:
        invA = np.linalg.inv(A)

        B = [4*pg1[0]-p0[0]+4*m1*pg1[1]-m1*p0[1],pg2[0]+m2*pg2[1],m0*p0[1]+p0[0],0]
        C = np.dot(invA,B)
        return C[0],C[1],C[2],C[3]
    except:
        return -1,-1,-1,-1
def printCurve(curve,swt):
    lmt = (len(curve)-1)/2 # we take first point and then two points in pair to generate quad bezier curves
    #plt.gca().invert_yaxis()

    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.set_xlim(0,swt.shape[1])
    #ax.set_ylim(swt.shape[0],0)
    p2 = curve[0]
    m2 = np.tan(np.pi/180*swt[curve[0][1],curve[0][0]])

    codes1 = []
    #for point in curve:
    #    codes1.append(Path.LINETO)
    #codes1[0] = Path.MOVETO
    #path = Path(curve,codes1)
    #patch = patches.PathPatch(path,facecolor='red', lw=1)
    #ax.add_patch(patch)
    #plt.show()
    for i in xrange(2):
        p0 = p2
        m0 = m2
        pg1 = curve[2*i+1]
        m1 = np.tan(np.pi/180*swt[curve[2*i+1][1],curve[2*i+1][0]])
        pg2 = curve[2*i+2]
        m2 = np.tan(np.pi/180*swt[curve[2*i+2][1],curve[2*i+2][0]])
        y1,y2,x1,x2 = quadCurve(p0,m0,pg1,m1,pg2,m2)
        if y1<0:
            continue
        p2 = (x2,y2)
        p1 = (x1,y1)
        verts=[p0,p1,p2]
        print verts,m0,m1,m2,pg1,pg2
    #    codes = [Path.MOVETO,Path.CURVE3,Path.CURVE3]
    #    path = Path(verts, codes)
    #    patch = patches.PathPatch(path,facecolor='none', lw=1)
        #xs, ys = zip(*verts)
        #ax.plot(xs, ys, 'x--', lw=2, color='black', ms=10)

        #ax.add_patch(patch)
    #plt.show()
    #fig.savefig("quadbezier.png")
