
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
from fitCurves import *

class Check(object):
    def correctSize(self,img):
        try:
            indices = np.where(img==0)
            leftx = min(indices[1])
            rightx = max(indices[1])
            topy = min(indices[0])
            bottomy = max(indices[0])
            if leftx>5:
                leftx = leftx -5
            if topy>5:
                topy = topy -5


            print leftx,rightx,topy,bottomy
            img = img[topy:bottomy+5,leftx:rightx+5]
            return img
        except:
            return img

    def drawRays(self,theta,img):
         indices = np.where(theta > 0)
         theta_color = cv2.cvtColor(theta,cv2.COLOR_GRAY2RGB)
         rays = []
         #cv2.imwrite("img.jpg",img)
         #plt.subplot(1,2,1),plt.imshow(img)
         #plt.subplot(1,2,2),plt.imshow(theta_degrees)
         #plt.show()
         for i,y in enumerate(indices[0]):
             x = indices[1][i]
             ray = []
             ray.append((x,y))
             theta_color[y,x] = (0,0,255)
             cur_y = y
             cur_x = x
             ray_length = 10
             angle = theta[y,x]
             radian = angle*np.pi/180


             while ray_length < 100:
                 #print img[cur_y,cur_x],ray_length

                 delta_x = ray_length * np.cos(radian)
                 delta_y =  ray_length * np.sin(radian)
                 #delta_x = ray_length * grad_x
                 #delta_y = ray_length * grad_y

                 cur_x = int(round(x + delta_x))
                 cur_y = int(round(y + delta_y))

                 try:
                     if img[cur_y,cur_x] ==0:
                         break
                     if abs(theta[cur_y,cur_x] - theta[y,x])< 15:
                         ray.append((cur_x,cur_y))
                         rays.append(ray)

                         break
                 except:
                     break
                 ray_length += 1

         cv2.imwrite('theta_color.jpg',theta_color)
         return rays





    def swt(self,theta1,img,edge,laplacian,maxlength = 20,errMax=5):

      theta = theta1.astype('int64')
      #plt.imshow(theta)
      #plt.show()
      indices = np.where(theta >0)
      #print indices.shape
      edge_color = cv2.cvtColor(edge,cv2.COLOR_GRAY2RGB)

      swt = np.zeros(img.shape)
      swt[:] = 0
      swt = np.uint8(swt)
      #swt = cv2.cvtColor(swt,cv2.COLOR_GRAY2RGB)
      #swt[:] = [255,255,255]

      print len(indices[0])
      print maxlength
      rays = []
      for i,y in enumerate(indices[0]):
          x = indices[1][i]
          ray_length = 5
          cur_x = x
          cur_y = y
          prev_x = x
          prev_y = y
          angle = float(theta[y,x])

          radian = (angle/180)*np.pi
          #if theta[y,x] ==90:
            #  print radian

          edge_color[y,x] = [0,0,255]
          #print abs(theta[y,x]-theta[cur_y,cur_x])
          #sys.exit()
          #print errMax, maxlength
          ray = []
          ray.append((x,y))
         # grad_x = grad_x_g[y,x]
          #grad_y = grad_y_g[y,x]
          while ray_length < maxlength:
              delta_x = ray_length * np.cos(radian)
              delta_y =  ray_length * np.sin(radian)

              #delta_x = ray_length * grad_x
              #delta_y = ray_length * grad_y

              cur_x = int(round(x + delta_x))
              cur_y = int(round(y + delta_y))

              #print x,y,cur_x,cur_y,delta_x,":",delta_y,':',ray_length,img[cur_y,cur_x]

              try:
                  if cur_x ==x and cur_y ==y:
                      break
                  if img[cur_y,cur_x] ==0:
                     break
                  #if edge[cur_y,cur_x] > 0 and ray_length > 10 :
                  #if ray_length ==maxlength or (abs(float(theta[cur_y,cur_x]) - float(theta[y,x]))< errMax or abs(abs(float(theta[cur_y,cur_x]) - float(theta[y,x]))-180)<errMax) :
                  if abs(theta[cur_y,cur_x]-theta[y,x]) <errMax:
                      median_x = (cur_x + x)/2
                      median_y = (cur_y + y)/2
                      if theta[median_y,median_x] >0:
                          break
                      median_theta = int(round((float(theta[cur_y,cur_x])+float(theta[y,x]))/2))
                      if median_theta ==0:
                          median_theta = 180
                      swt[median_y,median_x] = median_theta
                      ray.append((cur_x,cur_y))
                      rays.append(ray)

                      break
              except:
                  break
              ray_length += 1

      cv2.imwrite('edge-color.jpg',edge_color)
      self.group(swt)
      return swt,rays
    def updateVaues(self,img,edge,theta_normal,raylength):
        img1 = img.copy()
        indices = np.where(edge>0)
        for i,y in enumerate(indices[0]):
            x = indices[1][i]
            theta = theta_normal[y,x]
            img1[y,x] = theta
            for length in xrange(raylength):
                try:
                    cur_x = x + length*cos(theta)
                    cur_y = y + length*sin(theta)
                    #print cur_x
                    cur_x = int(round(cur_x))
                    cur_y = int(round(cur_y))
                    img1[cur_y,cur_x] = theta
                except:
                    break

        return img1


    def printRays(self,rays,edges,filename):
         k = 0
         edges_copy = cv2.cvtColor(edges,cv2.COLOR_GRAY2RGB)
         for ray in rays:
             k = k+ 1
             #cv2.polylines(edges,ray,(255,0,0),3)
             if k%1 == 0:
                 cv2.line(edges_copy,ray[0],ray[1],(0,0,255),1)
         cv2.imwrite("rays-"+filename,edges_copy)
    def printCurves(self,img,curves):
         img = img.astype('uint8')

         img_copy = cv2.cvtColor(img,cv2.COLOR_GRAY2RGB)
         for i,curve in enumerate(curves):
             for k,point in enumerate(curve):
                  img_copy[point[1],point[0]] = (0,255,0)
         cv2.imwrite("curves.jpg",img_copy)
    def skeleton(self,img):
        #skeleton = union of all skeletons
        # each skeleton = (A eroded k times by B) and then negative by (A eroded k times by B) opening by B
        # opening is dilation followed by erosion
        kernel = np.ones((3,3),np.uint8)
        #kernel[0,0]=0
        #kernel[0,2]=0
        #kernel[2,0]=0
        #kernel[2,2]=0

        result = np.zeros(img.shape)
        tophat = cv2.morphologyEx(img, cv2.MORPH_TOPHAT, kernel)
        result = result + tophat
        kerosion = img.copy()
        k=0
        while len(np.where(kerosion>0)[0])>0:
            kerosion = cv2.erode(kerosion,kernel,iterations = 1) # k = 4
            k += 1
        kerosion = img.copy()
        for i in xrange(k):
            kerosion = cv2.erode(kerosion,kernel,iterations = 1) # k = 4
            #tophat
            #opening = cv2.morphologyEx(kerosion, cv2.MORPH_OPEN, kernel)
            tophat = cv2.morphologyEx(kerosion, cv2.MORPH_TOPHAT, kernel)
            result = result + tophat

        #plt.imshow(result)
        #plt.show()
        result = np.uint8(result)
        cv2.imwrite('result.jpg',result)
        return result
        '''
        edge_result = cv2.Canny(result,100,255)
        im2, contours, hierarchy = cv2.findContours(edge_result,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
        img1 = img.copy()
        img1[:] = 0
        img1= cv2.cvtColor(img1,cv2.COLOR_GRAY2RGB)
        print len(contours)
        print len(hierarchy)
        #cv2.drawContours(img1,contours,10,255,1)
        print hierarchy[0]
        #print contours[0]
        for i,cnt in enumerate(contours):
            if hierarchy[0][i][3] != -1:
                continue
            clr1 = int(i%2*255)
            clr2 = int(i%3*255)
            font = cv2.FONT_HERSHEY_PLAIN
            if i%5 ==0:
                cv2.putText(img1,str(i),(cnt[0][0][0] +2 ,cnt[0][0][1]), font, 0.6,(255,255,255),1,cv2.LINE_AA,False)
            cv2.drawContours(img1,[cnt],0,(clr1,clr2,255),1)
        plt.imshow(img1)
        plt.show()
        '''
    def group(self,img):
        #produce contonuous chain of points - starting at source
        # first finds all straight lines
        #scan horizontally
        img = img.astype('int64')
        rows = img.shape[0]
        cols = img.shape[1]
        indices = np.argsort(img,axis=None)
        curves = []
        visited = np.zeros(img.shape)
        visited = np.uint8(visited)
        noofiterations = 1
        while len(indices)>0 and noofiterations <10:
            indices1 = indices.copy()
            curve = []
            prev_x = -1
            prev_y = -1
            for l,no in enumerate(indices):
                y = no/cols
                x = no%cols
                #print l,no,y,x
                #sys.exit()
                if img[y,x] ==0:
                    np.delete(indices1,l)
                    #prev_x = x
                    #prev_y = y
                    #visited[y,x] = 255
                    continue

                if visited[y,x] > 0:
                    np.delete(indices1,l)
                    #prev_x = x
                    #prev_y = y
                    continue
                #visited[y,x] = 255
                dist = 100
                if prev_x > 0:
                    if np.linalg.norm(np.array([x,y])-np.array([prev_x,prev_y])) > 10 or abs(img[y,x]-img[prev_y,prev_x])>5:
                        #prev_x = x
                        #prev_y = y
                        continue
                    else:
                        np.delete(indices1,l)
                        visited[y,x] = 255
                        curve.append([x,y])
                        prev_x = x
                        prev_y = y

                else:
                    curve.append([x,y])
                    visited[y,x] = 255
                    np.delete(indices1,l)
                    prev_x = x
                    prev_y = y

            if len(curve) == 0:
                #curves.append(indices)
                break
            else:
                curves.append(curve)

                indices = indices1.copy()
            noofiterations += 1


            #nbrs = [[y+1,x],[y-1,x],[y-1,x+1],[y+1,x+1],[y+1,x-1],[y-1,x-1]]
        #print curves
        #self._printBezier([curves[0]],img)
        self.printCurves(img,[curves[0]])
        cv2.imwrite('visited.jpg',visited)
        return curves
    def boundary(self,img):
        kernel = np.ones((3,3),np.uint8)
        erosion = cv2.erode(img,kernel,iterations=1)
        return img-erosion
    def hitormiss(self,img,b1,b2=None):
        erosion = cv2.erode(img,b1)
        if b2 != None:
          dilation= cv2.dilate(img,b2)
          return erosion & dilation
        else:
            return erosion

    def intersection(self,img1,img2):
        assert img1.shape == img2.shape
        inv_img = np.zeros(img2.shape)
        inv_img[:] = 255
        inv_img = inv_img - img2
        return img1 & img2

    def thinning (self,img):
        assert 1 ==1




    def _printBezier(self,curves,img):
        plt.gca().invert_yaxis()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim(0,img.shape[1])
        ax.set_ylim(img.shape[0],0)
        for curve in curves:
            #print curve
            curve=np.array(curve)
            #print curve
            #target.write("path.move(to:CGPoint(x:"+str(curve[0][0])+",y:"+str(curve[0][1])+"))\n")

            if len(curve)==2: # straight line
            #    target.write("path.addLine(to:CGPoint(x:"+str(curve[1][0])+",y:"+str(curve[1][1])+"))\n")
                verts=[curve[0],curve[1]]
                codes = [Path.MOVETO,Path.LINETO]
                path = Path(verts, codes)
                patch = patches.PathPatch(path,facecolor=None ,lw=2)
                ax.add_patch(patch)

            else:
                beziers = fitCurve(curve,10.0)
                #beziers = fitCurve(currentContour, 10)
                #print(len(beziers))
                for bezier in beziers:
                    #print bezier
                    #print bezier[0][0]
                    if not (math.isnan(bezier[1][0])) and not (math.isnan(bezier[2][0])):
                        #print("\n")
            #            target.write("path.addCurve(to:CGPoint(x:"+str(bezier[3][0])+",y:"+str(bezier[3][1])+"),controlPoint1:CGPoint(x:"+str(bezier[1][0])+",y:"+str(bezier[1][1])+"),controlPoint2:CGPoint(x:"+str(bezier[2][0])+",y:"+str(bezier[2][1])+"))\n")
                        verts=[bezier[0],bezier[1],bezier[2],bezier[3]]
                        codes = [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4]
                        path = Path(verts, codes)
                        patch = patches.PathPatch(path,facecolor='none', lw=1)
                        #xs, ys = zip(*verts)
                        #ax.plot(xs, ys, 'x--', lw=2, color='black', ms=10)

                        ax.add_patch(patch)



        plt.show()
        fig.savefig("bezier.png")
    def plot_comparison(original, filtered, filter_name):

        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(8, 4), sharex=True, sharey=True)
        ax1.imshow(original, cmap=plt.cm.gray)
        ax1.set_title('original')
        ax1.axis('off')
        ax1.set_adjustable('box-forced')
        ax2.imshow(filtered, cmap=plt.cm.gray)
        ax2.set_title(filter_name)
        ax2.axis('off')
        ax2.set_adjustable('box-forced')
        plt.show(block=False)



'''
img = cv2.imread("letters/91b.png",0)
#img = cv2.imread("result.jpg",0)

#plt.imshow(img)
#plt.show()

c = Check()
img = c.correctSize(img)

_,inv_img = cv2.threshold(img,40,255,cv2.THRESH_BINARY_INV)
cv2.imwrite("invimg.jpg",inv_img)
edges = cv2.Canny(inv_img,20,255)
bdry = c.boundary(inv_img)
edges = bdry
sobelx = cv2.Sobel(edges,cv2.CV_64F,1,0,ksize=-1)
sobely = cv2.Sobel(edges,cv2.CV_64F,0,1,ksize=-1)
laplacian = cv2.Laplacian(edges,cv2.CV_64F)

theta = np.arctan2(sobely,sobelx)
theta_degrees = np.uint8(theta * 180/np.pi)
theta_degrees[np.where(theta_degrees<0)] = 0
cv2.imwrite("thetadegrees.jpg",theta_degrees)

#theta_degrees[np.where(theta_degrees<0)] = 0
#c.skeleton(inv_img)
#cv2.imwrite('bdry.jpg',bdry)
swt,rays=c.swt(theta_degrees,inv_img,edges,laplacian,maxlength=50,errMax=5)
#c.printRays(rays,theta_degrees,"test.jpg")

#c.printRays(c.drawRays(theta_degrees,inv_img),theta_degrees,"test.jpg")
cv2.imwrite("swt-new.jpg",swt)
#plt.imshow(swt,'gray')
#plt.show()
'''
#
