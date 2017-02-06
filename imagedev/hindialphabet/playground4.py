import sys
import cv2
import numpy as np
from skimage.morphology import black_tophat, skeletonize, convex_hull_image
from playground3 import Check
from findstroke import FindStroke
from matplotlib import pyplot as plt
def diskSE(center,radius):
    size = (2*radius,2*radius)
    circkernel = np.zeros(size,dtype=np.uint8)
    circkernel = cv2.circle(circkernel,center,radius,1,-1)
    return circkernel

def removeNbrs(img):
    indices = np.where(img>0)
    for i,y in enumerate(indices[0]):
        x = indices[1][i]
        nbrs = [[x+1,y+1],[x+1,y],[x+1,y-1],[x+1,y],[x+1,y],[x+1,y],[x+1,y]]
fs = FindStroke("x",".",".")
c= Check()

img =cv2.imread("letters-small/916.png",0)
_,img = cv2.threshold(img,40,255,cv2.THRESH_BINARY_INV)


res = cv2.resize(img,None,fx=2, fy=2, interpolation = cv2.INTER_CUBIC) #scale by two to get even stroke width
edges1 = cv2.Canny(res,40,255)
sobelx64f = cv2.Sobel(res,cv2.CV_64F,1,0,ksize=3)
sobely64f = cv2.Sobel(res,cv2.CV_64F,0,1,ksize=3)
theta_normal = np.arctan2(sobely64f,sobelx64f)

#theta_degree = (theta_normal+np.pi) *255/(2*np.pi)
theta_degree = theta_normal * 180/np.pi
laplacian = cv2.Laplacian(edges1,cv2.CV_64F)


srekhaheightn = fs.findShirorekhaHeight(res,edges1)
swt,rays = c.swtNew(theta_degree,res,c.boundary(edges1),laplacian,maxlength=srekhaheightn/2)
#curves = fs._trace(swt,"write_new.swift",srekhaheightn,2)

cv2.imwrite("swt-improved.jpg",swt)
sys.exit()
#res[np.where((sobely64f==0) & (edges1>0))] = 0
#res[np.where((sobelx64f==0) & (edges1>0))] = 0
circSE1 =   diskSE ((srekhaheightn/2,srekhaheightn/2),srekhaheightn/2)
eroded1 = cv2.erode(res,circSE1)
eroded = eroded1.copy()
_,eroded1 = cv2.threshold(eroded1,40,255,cv2.THRESH_BINARY)
eroded1_copy = eroded1.copy()
#eroded1 = np.uint64(eroded1) # to prevent overflow
#eroded2 = np.uint64(eroded2) # to prevent overflow
for r in xrange(8):
    circSE2 =   diskSE ((srekhaheightn/2,srekhaheightn/2),srekhaheightn/2+r+2)
    circSE2[srekhaheightn/2+1,2]=1
    circSE2[2,srekhaheightn/2+1]=1
    circSE2[2,2]=1
    circSE2[srekhaheightn/2+1,srekhaheightn/2+1]=1



    eroded2 = cv2.erode(res,circSE2)
    _,eroded2 = cv2.threshold(eroded2,40,255,cv2.THRESH_BINARY)

    eroded1 = eroded1-eroded2
    _,eroded1 = cv2.threshold(eroded1,40,255,cv2.THRESH_BINARY)
eroded1_n = eroded1_copy.copy()
for r in xrange(8):
    circSE2 =   diskSE ((srekhaheightn/2,srekhaheightn/2),srekhaheightn/2+r+2)



    eroded2 = cv2.erode(res,circSE2)
    _,eroded2 = cv2.threshold(eroded2,40,255,cv2.THRESH_BINARY)

    eroded1_n = eroded1_n-eroded2
    _,eroded1_n = cv2.threshold(eroded1,40,255,cv2.THRESH_BINARY)
result = eroded1_n - eroded1_copy
_,result = cv2.threshold(result,0,255,cv2.THRESH_BINARY)

#result = eroded1_copy - result
#opening = cv2.morphologyEx(eroded_copy, cv2.MORPH_OPEN, circ1kernel)
#eroded_copy[np.where(eroded_copy==255)] =0
#bdry = c.boundary(eroded)
#bdry = c.boundary(bdry)
#eroded = eroded & eroded1
result1 = (eroded1_copy-result)

sk = c.skeleton(result1)
sk = c.boundary(sk)
sk = sk & eroded


raylength = 20
sobelx64f = cv2.Sobel(res,cv2.CV_64F,1,0,ksize=-1)
sobely64f = cv2.Sobel(res,cv2.CV_64F,0,1,ksize=-1)
theta_normal = np.arctan2(sobely64f,sobely64f)
theta_degree = (theta_normal+np.pi) *255/(2*np.pi)

img1 = c.updateVaues(res,edges1,theta_degree,raylength)
#ret,thresh = cv2.threshold(img1,40,1,cv2.THRESH_BINARY)
#img2 = skeletonize(thresh)
#img2 = c.skeleton(img)
#plt.imshow(img2)
#plt.show(block=False)
indices = np.where(sk>0)
for i,y in enumerate(indices[0]):
    x = indices[1][i]
    sk[y,x] = img1[y,x]
curves = fs._trace(sk,"write_new.swift",15,2)

'''
sobelx64f = cv2.Sobel(img,cv2.CV_64F,1,0,ksize=-1)
sobely64f = cv2.Sobel(img,cv2.CV_64F,0,1,ksize=-1)
sobelx_n = sobelx64f.copy()
sobelx_n[np.where(sobelx_n!=0)] = 255
sobely_n = sobely64f.copy()
sobely_n[np.where(sobely_n!=0)] = 255
horvertlines = sobely_n -sobelx_n
horvertlines[np.where(horvertlines<0)] = 255
horvertlines = np.uint8(horvertlines)
inv = np.zeros(img.shape)
inv[:] = 255
inv = np.uint8(inv)
horvertmask = inv-horvertlines # mask
img = img & horvertmask

theta_normal = np.arctan2(sobely64f,sobely64f)
theta_degree = (theta_normal+np.pi) *255/(2*np.pi)

cv2.imwrite("theta_degree.png",theta_degree)
cv2.imwrite("theta_normal.png",theta_normal)

'''
