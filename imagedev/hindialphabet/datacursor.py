import sys
sys.path.append('/Users/mac/.virtualenvs/cv/lib/python2.7/site-packages')
import cv2
import matplotlib.pyplot as plt
import numpy as np
from mpldatacursor import datacursor
img = cv2.imread("ka.jpg",0)
edges = cv2.Canny(img, 175, 320, apertureSize=3)
#edges = findSinglePixelEdge(img)
# Create gradient map using Sobel
sobelx64f = cv2.Sobel(img,cv2.CV_64F,1,0,ksize=3)
sobely64f = cv2.Sobel(img,cv2.CV_64F,0,1,ksize=3)

theta = np.arctan2(sobely64f, sobelx64f)
theta_visible = (theta + np.pi)*255/(2*np.pi)
#images.append(edges)
#images.append(np.absolute(sobelx64f))
#images.append(np.absolute(sobely64f))
#images.append(theta_visible)

#if diagnostics:
##     cv2.imwrite('sobelx64f.jpg', np.absolute(sobelx64f))
#     cv2.imwrite('sobely64f.jpg', np.absolute(sobely64f))
    # amplify theta for visual inspection
   #
data1 = cv2.imread("swt-exp.jpg",0)
#data1 = data1 * 255/np.pi
data = theta_visible
fig, axes = plt.subplots(ncols=1)
#axes[0].imshow(data, interpolation='nearest', origin='upper')
axes.imshow(data1, interpolation='nearest', origin='upper')#,
#                     extent=[200, 300, 400, 500])
datacursor(display='single')

fig.suptitle('Click anywhere on the image')

plt.show()
