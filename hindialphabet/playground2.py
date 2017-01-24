
import sys
sys.path.append('/Users/mac/.virtualenvs/cv/lib/python2.7/site-packages')
import numpy as np
import cv2
import time
import math
import operator
from matplotlib import pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
filename = "letters/91b.png"
swtfile = "test0x091B/swt.jpg"
#swtfile = "swt-exp.jpg"
img1 = cv2.imread(filename)
swt = cv2.imread(swtfile)
img = swt
for x in xrange(img.shape[1]):
    for y in xrange(img.shape[0]):
        if np.array_equal(img[y,x],[255,255,255]):
            img[y,x] = [0,0,0]
kernel = np.array([[1,1,1],[1,-8,1],[1,1,1]])

imgLaplacian = cv2.filter2D(img,cv2.CV_32F,kernel)
imgResult = img - imgLaplacian
imgResult = np.uint8(imgResult)
imgLaplacian = np.uint8(imgLaplacian)
img = imgResult
bw = cv2.cvtColor(img,  cv2.COLOR_BGR2GRAY)
_,bw = cv2.threshold(bw, 40, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)
dist=    cv2.distanceTransform(bw,  cv2.DIST_L2, 3);
#cv2.imwrite("dist-cha.jpg",dist)
dist=cv2.normalize(dist, 0, 1., cv2.NORM_MINMAX);
#plt.imshow(dist,'gray')
#plt.show(block=False)

_,dist = cv2.threshold(dist,0.4,1.,cv2.THRESH_BINARY)
kernel1 = np.ones([3,3])
dist = cv2.dilate(dist,kernel1)
dist = np.uint8(dist)
contours,hier = cv2.findContours(dist,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
markers = np.zeros(dist.shape)
markers = np.uint8(markers)
markers = cv2.cvtColor()
#plt.imshow(dist)
#plt.show()
