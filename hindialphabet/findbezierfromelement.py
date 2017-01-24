import sys
sys.path.append('/Users/mac/.virtualenvs/cv/lib/python2.7/site-packages')
import numpy as np
import cv2
import time
import math
from matplotlib import pyplot as plt

class FindBezier(object):
    @classmethod
    def _create_derivative(cls, img,filenamePrefix="",diagnostics=False):
         #img = cv2.imread(filepath,0)
         edges = cv2.Canny(img, 175, 320, apertureSize=3)
         # Create gradient map using Sobel
         sobelx64f = cv2.Sobel(img,cv2.CV_64F,1,0,ksize=3)
         sobely64f = cv2.Sobel(img,cv2.CV_64F,0,1,ksize=3)
         theta = np.arctan2(sobely64f, sobelx64f)
         theta_visible = (theta + np.pi)*255/(2*np.pi)

         if diagnostics:
             cv2.imwrite(filenamePrefix+'edges.jpg',edges)
             cv2.imwrite(filenamePrefix+'sobelx64f.jpg', np.absolute(sobelx64f))
             cv2.imwrite(filenamePrefix+'sobely64f.jpg', np.absolute(sobely64f))
             # amplify theta for visual inspection

             cv2.imwrite(filenamePrefix+'theta.jpg', theta_visible)
         return (edges, sobelx64f, sobely64f, theta,theta_visible)

    @classmethod
    def findHorizontal(self,img):
        (edges, sobelx64f, sobely64f, theta,theta_visible) = self._create_derivative(img)
        lines = {} # yvalues of lines encountered
        img = np.zeros(edges.shape)
        points = []

        count = 0
        #abssobel = np.absolute(-1 * sobelx64f)
        #ret,thresh1 = cv2.threshold(sobelx64f,10,255,cv2.THRESH_BINARY_INV)
        theta_copy = theta_visible.copy()
        theta_copy = np.asarray(theta_copy)
        vertical1 = 3/2 *255
        vertical2 = 1/2 * 255
        delta = 5
        a_mask1 = (theta_copy >vertical1-delta) & (theta_copy < vertical1+delta)
        a_mask2 = (theta_copy >vertical2-delta) & (theta_copy < vertical2+delta)


        # | ((theta_copy <128) & (theta_copy > 120))
        theta_visible[a_mask1] =255
        theta_visible[a_mask2] =255
        theta_visible[~a_mask1] =0
        theta_visible[~a_mask2] =0

        cv2.imwrite("theta.jpg",theta_visible)


src ="ka.jpg"
img = cv2.imread(src,0)
FindBezier.findHorizontal(img)
