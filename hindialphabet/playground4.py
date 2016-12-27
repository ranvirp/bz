import cv2
import numpy as np
from skimage.morphology import black_tophat, skeletonize, convex_hull_image
from playground3 import Check
from findstroke import FindStroke
from matplotlib import pyplot as plt

img =cv2.imread("letters-small/91d.png",0)
inv_img = np.zeros(img.shape)
inv_img[:] = 255
inv_img = np.uint8(inv_img)
img = inv_img - img
sobelx64f = cv2.Sobel(img,cv2.CV_64F,1,0,ksize=-1)
sobely64f = cv2.Sobel(img,cv2.CV_64F,0,1,ksize=-1)
theta_normal = np.arctan2(sobely64f,sobely64f)
theta_degree = (theta_normal+np.pi) *255/(2*np.pi)

cv2.imwrite("theta_degree.png",theta_degree)
cv2.imwrite("theta_normal.png",theta_normal)
c= Check()
edge = c.boundary(img)
raylength = 20
img1 = c.updateVaues(img,edge,theta_degree,raylength)
ret,thresh = cv2.threshold(img,40,1,cv2.THRESH_BINARY)
img2 = skeletonize(thresh)
#img2 = c.skeleton(img)
#plt.imshow(img2)
#plt.show(block=False)
indices = np.where(img2>0)
for i,y in enumerate(indices[0]):
    x = indices[1][i]
    img2[y,x] = img1[y,x]
fs = FindStroke("x",".",".")
fs._trace(img2,"write_new.swift",15,50)
