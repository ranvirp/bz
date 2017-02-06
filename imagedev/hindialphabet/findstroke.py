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
import os
import errno

from fitCurves import *
from scipy.interpolate import interp1d



def compare (item1, item2):
    return cmp(item1[-1][0], item2[0][0])
def compare1 (item1,item2):
    return cmp(len(item2), len(item1))
def compare2 (item1, item2):
    return cmp(item2[1], item1[1])
def compare3 (item1, item2):
    return cmp(item1[0], item2[0])
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

class FindStroke(object):
    #unicodeno = 0x097F
    #workingdir = ""
    #lettersdir = ""
    #diagnostics = ""
    @classmethod
    def __init__(self,unicodeno,lettersdir="letters",workingdir=".",diagnostics=False):
        self.unicodeno = unicodeno
        self.lettersdir = lettersdir
        self.workingdir = workingdir
        self.diagnostics = diagnostics
        make_sure_path_exists(workingdir)

    @classmethod
    def process(self):
        filename = self.lettersdir + "/" + self.findUnicodeStringForCodeWord(self.unicodeno)+".png"
        print filename
        img = cv2.imread(filename,0)
        imgcolor = cv2.cvtColor(img,cv2.COLOR_GRAY2RGB)
        (imgbeginx,imgendx,imgbeginy,imgendy) = self.find_bounding_box(img)
        print (imgbeginx,imgendx,imgbeginy,imgendy)
        #sys.exit()
        #img = img[imgbeginy-1:imgendy+1,imgbeginx-1:imgendx+1]
        #imgcolor = imgcolor[imgbeginy-1:imgendy+1,imgbeginx-1:imgendx+1]
        (edges, sobelx64f, sobely64f, theta_normal,theta_along_edge,theta_visible) = self._create_derivative(img)
        inv_image = np.empty(img.shape)
        inv_image[:] = 255
        inv_img = inv_image - img
        cv2.imwrite('inv.jpg',inv_img)
        inv_img = cv2.imread('inv.jpg',0)
        ret,thresh1 = cv2.threshold(inv_img,10,255,cv2.THRESH_BINARY)
        edges = cv2.Canny(thresh1,10,255)
        kernel = np.ones((1,img.shape[1]/0.6),np.uint8)
        #kernel = np.ones((5,5),np.uint8)
        erosion = cv2.erode(edges,kernel,iterations = 1)
        #erosion = cv2.erode(edges,kernel,iterations = 1)
        dilation = cv2.dilate(erosion,kernel,iterations = 1)
        edge_without_horizontal = edges - dilation
        kernel = np.ones((edges.shape[0]/0.6,1),np.uint8)
        #kernel = np.ones((5,5),np.uint8)
        erosion = cv2.erode(edge_without_horizontal,kernel,iterations = 1)
        #erosion = cv2.erode(edges,kernel,iterations = 1)
        dilation = cv2.dilate(erosion,kernel,iterations = 1)
        edge_without_horizontal_vertical = edge_without_horizontal - dilation
        edges = edge_without_horizontal_vertical
        #kernel = np.ones((10,10),np.uint8)
        #kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(5,5))
        #img1 = cv2.imread("swt-new-with-theta.jpg",0)
        #img2=cv2.morphologyEx(img1, cv2.MORPH_CLOSE, kernel)
        #cv2.imwrite("morphed.jpg",img2)
        #cv2.imwrite("morphed-before.jpg",img1)
        if self.diagnostics:
            cv2.imwrite(self.workingdir+'/edges-after-extend.jpg',edges)
            cv2.imwrite(self.workingdir+'/firstprocessing.jpg',img)

        shirorekhaheight = self.findShirorekhaHeight(inv_image,edges)
        print shirorekhaheight
        thickness = shirorekhaheight * 1.5
        self.extendLines(inv_img,edges,theta_normal)
        #edges = self.extendEdge(inv_img,edges,theta_normal,thickness)

        (img,swt,rays,rays_others) = self.swt1(theta_normal, edges , sobelx64f, sobely64f,inv_img,initialthickness=thickness,slopeErr=1)
        swt = np.uint8(swt)
        #swt = cv2.distanceTransform(swt,cv2.DIST_L2,3)
        #swt = cv2.normalize(swt,0,255,cv2.NORM_MINMAX)
        #_,swt =cv2.threshold(swt,10,255,cv2.THRESH_BINARY_INV)
        #cv2.imwrite(self.workingdir+"/"+"swt-distance.jpg",swt)
        #self.remVerticalHor(swt)
        #swt = cv2.GaussianBlur(swt,(5,5),1.5);
        #swt = cv2.Canny(swt, 1, 320, apertureSize=3)
        if self.diagnostics:
          cv2.imwrite(self.workingdir+'/swt-exp.jpg',swt)
        #swt = cv2.imread('result.jpg',0)
        self._trace(swt,self.workingdir+"/write.swift",distMax=thickness/2,errMax=250/50)


    @classmethod
    def findUnicodeStringForCodeWord(self,no):
        return "{:02x}".format(no)
    @classmethod
    def find_bounding_box(self,img):
        imgbeginy = 0
        imgendy = 0
        imgbeginx = 0
        imgendx = 0

        for j in xrange(img.shape[0]):
            #vert.append(
            c = np.count_nonzero(img[j,:])
            if c > 10:
                imgendy = j
                if imgbeginy == 0:
                   imgbeginy = j


        for j in xrange(img.shape[1]):
            #vert.append(
            c = np.count_nonzero(img[:,j])
            if c > 10:
                imgendx = j
                if imgbeginx == 0:
                   imgbeginx = j
        return (imgbeginx,imgendx,imgbeginy,imgendy)
    @classmethod
    def _create_derivative(self, img):
         edges = cv2.Canny(img, 175, 320, apertureSize=3)
         sobelx64f = cv2.Sobel(img,cv2.CV_64F,1,0,ksize=-1)
         sobely64f = cv2.Sobel(img,cv2.CV_64F,0,1,ksize=-1)
         theta_along_edge = np.arctan2(sobelx64f,-1*sobely64f)
         theta_normal = np.arctan2(sobely64f,sobely64f)
         theta_visible = (theta_normal + np.pi)*255/(2*np.pi) # to ensure none of it is zero
         if self.diagnostics:
            cv2.imwrite(self.workingdir+'/sobelx64f.jpg', np.absolute(sobelx64f))
            cv2.imwrite(self.workingdir+'/sobely64f.jpg', np.absolute(sobely64f))
            cv2.imwrite(self.workingdir+'/edges.jpg', edges)
            cv2.imwrite(self.workingdir+'/img.jpg', img)


             #amplify theta for visual inspection
            cv2.imwrite(self.workingdir+'/theta.jpg', theta_visible)
         return (edges, sobelx64f, sobely64f, theta_normal,theta_along_edge,theta_visible)
    @classmethod
    def printRays(self,rays,edges,filename):
         k = 0
         edges_copy = edges.copy()
         for ray in rays:
             k = k+ 1
             #cv2.polylines(edges,ray,(255,0,0),3)
             if k%2 == 0:
                 cv2.line(edges_copy,ray[0],ray[1],(255,0,0),3)
         cv2.imwrite(self.workingdir+"/"+"rays-"+filename,edges_copy)
    @classmethod
    def findWidthOfStroke(self,point,theta,edges,maxWidth):
        i = 0
        x = point[0]
        y = point[1]
        prev_x = x
        prev_y = y
        while i < maxWidth:
            i += 1
            cur_x = int(math.floor(x + i * np.cos(theta) ))
            cur_y = int(math.floor(y + i * np.sin(theta) ))
            if edges[cur_y,cur_x] > 10:
                return [np.linalg.norm(np.array((cur_x,cur_y))-np.array((x,y))),[point,(cur_x,cur_y)]]
        return [] #used in dealWithOtherRays
    @classmethod
    def dealWithOtherRays(self,swt,rays_others,edges,theta_visible,maxWidth):

         # Try moving the ray 10 degrees plus and minus and see if the length of the ray increases or decreases,
         # if length increases both sides then keep the ray otherwise keep moving in the direction where it is decreasing
         # to find the lowest length
         # find angle of ray- get a direction +10 and -10 ..s
         # find length of stroke at angle theta from point
         for i,ray in enumerate(rays_others):
             length = np.linalg.norm(np.array(ray[0])-np.array(ray[1]))
             theta = np.arctan2(ray[0][1]-ray[1][1],ray[0][0]-ray[1][0])
             dist1 = self.findWidthOfStroke(ray[0],theta + 10/180*np.pi,edges,maxWidth)
             dist2 = self.findWidthOfStroke(ray[0],theta - 10/180*np.pi,edges,maxWidth)
             try:
                 if dist1[0] > length and dist2[0] > length :
                     # this ray should be included
                     median = (np.array(ray[1]) + np.array(ray[0]))/2
                     swt[median[1],median[0]] = max(theta_visible[ray[0][1],ray[0][0]],theta_visible[ray[1][1],ray[1][0]])
                     rays_others.pop(i)
                 elif dist1[0]>length and dist2[0] < length:
                      median = (np.array(dist2[1]) + np.array(ray[0]))/2
                      swt[median[1],median[0]] = max(theta_visible[ray[0][1],ray[0][0]],theta_visible[dist2[1][1][1],dist2[1][1][0]])
                      rays_others.pop(i)
             except IndexError:
                continue




         return (swt,rays_others)
    @classmethod
    def swt(self,theta, edges, sobelx64f, sobely64f,img,initialthickness=75,slopeErr=0.1):
                 theta_visible = (theta + np.pi)*255/(2*np.pi)
                 newedges = np.empty_like(edges)
                 newedges[:] = edges
                 # create empty image, initialized to infinity
                 swt = np.zeros(theta.shape)
                 #swt[:] = np.Infinity
                 rays = []
                 rays_others = []

                 #print time.clock() - t0

                 # now iterate over pixels in image, checking Canny to see if we're on an edge.
                 indices = np.where((edges>0)) # & (sobelx64f>=0))
                 for k,y in enumerate(indices[0]):
                    # for y in xrange(edges.shape[0]):
                             x = indices[1][k]
                    #     if edges[y, x] > 0 and sobelx64f[y,x] >= 0:
                             #step_x = step_x_g[y, x]
                             #step_y = step_y_g[y, x]
                             #mag = mag_g[y, x]
                             #grad_x = grad_x_g[y, x]
                             #grad_y = grad_y_g[y, x]
                             grad_x =  int(round(np.cos(np.pi/2+theta[y,x])))

                             grad_y =  int(round(np.sin(np.pi/2+ theta[y,x])))


                             ray = []
                             ray.append((x, y))
                             newedges[y,x] = 0
                             prev_x, prev_y, i = x, y, 0
                             while True:
                                 i += 1

                                 cur_x = int(math.floor(x + grad_x * i))
                                 cur_y = int(math.floor(y + grad_y * i))
                                 if cur_x != prev_x or cur_y != prev_y:
                                     # we have moved to the next pixel!
                                     try:
                                         if edges[int(cur_y), int(cur_x)] > 0  :
                                             # found edge,
                                             #print("found edge")


                                             ray.append((int(cur_x), int(cur_y)))
                                             theta_point = theta[y, x]
                                             alpha = theta[int(cur_y), int(cur_x)]
                                             thickness = math.sqrt( (cur_x - x) * (cur_x - x) + (cur_y - y) * (cur_y - y) )
                                             #if abs(thickness -initialthickness)/initialthickness
                                             #if abs(abs(alpha - theta_point) - np.pi) < slopeErr :

                                             if abs(alpha + theta_point) -np.pi  < slopeErr:

                                             #if True: or thickness < initialthickness:
                                                 if  (thickness < initialthickness):

                                                     midx = int ((x+cur_x)/2)
                                                     midy = (int (y+cur_y)/2)
                                                     #print theta_visible[y,x]
                                                     if sobelx64f[y,x] >=0 :
                                                        swt[midy,midx] = theta_visible[y,x]
                                                     else:
                                                        swt[midy,midx] = theta_visible[cur_y,cur_x]

                                                     #swt[midy,midx] = max(theta_visible[int(cur_y),int(cur_x)],theta_visible[y,x])
                                                     rays.append(ray)
                                                     newedges[int(cur_y),int(cur_x)] = 0


                                                 #for (rp_x, rp_y) in ray:
                                                   #  swt[rp_y, rp_x] = min(thickness, swt[rp_y, rp_x])
                                                 #rays.append(ray)
                                             else:
                                                 print "alpha =",alpha, " theta", theta_point," thickness=",thickness

                                                 rays_others.append(ray)

                                             break

                                             break
                                         # this is positioned at end to ensure we don't add a point beyond image boundary
                                         #ray.append((int(cur_x), int(cur_y)))
                                     except IndexError:
                                         # reached image boundary
                                         break
                                     except ValueError:
                                         break
                                     prev_x = cur_x
                                     prev_y = cur_y

                 (swt,rays_others) = dealWithOtherRays(swt,rays_others,edges,theta_visible,initialthickness*1.1)
                 print "now printing"
                 printRays(rays,edges,"edge-experimental.jpg")
                 printRays(rays_others,edges_copy,"others-edge-experimental.jpg")

                 return (img,swt,rays,rays_others)
    @classmethod
    def swt1(self,theta, edges, sobelx64f, sobely64f,img,initialthickness=75,slopeErr=0.1):
                 theta_visible = (theta + np.pi)*205/(2*np.pi) + 50
                 newedges = np.empty_like(edges)
                 newedges[:] = edges
                 # create empty image, initialized to infinity
                 swt = np.zeros(theta.shape)
                 #swt[:] = np.Infinity
                 rays = []
                 rays_others = []

                 #print time.clock() - t0

                 # now iterate over pixels in image, checking Canny to see if we're on an edge.
                 # if we are, follow a normal a ray to either the next edge or image border
                 # edgesSparse = scipy.sparse.coo_matrix(edges)
                 step_x_g = -1 * sobelx64f
                 step_y_g = -1 * sobely64f
                 mag_g = np.sqrt( step_x_g * step_x_g + step_y_g * step_y_g )
                 grad_x_g = step_x_g / mag_g
                 grad_y_g = step_y_g / mag_g
                 #normal_g = -1* sobely64f/sobelx64f
                 indices = np.where((edges > 0))# & ((sobelx64f >0)  | (sobely64f < 0)))
                 #print len(indices[0])
                 #print indices
                 #sys.exit()
                 #for x in xrange(edges.shape[1]):
                #     for y in xrange(edges.shape[0]):
                 for loopvariable,value in enumerate(indices[0]):
                             y = value
                             x = indices[1][loopvariable]
                            # if edges[y, x] > 0 and (sobelx64f[y,x] >0  or sobely64f[y,x] < 0):
                             #step_x = step_x_g[y, x]
                             #step_y = step_y_g[y, x]
                             #mag = mag_g[y, x]
                             grad_x = grad_x_g[y, x]
                             grad_y = grad_y_g[y, x]
                             if math.isnan(grad_x):
                                  continue


                             ray = []
                             ray.append((x, y))
                             newedges[y,x] = 0
                             prev_x, prev_y, i = x, y, 0
                             while True:
                                 i += 1

                                 cur_x = int(math.floor(x + grad_x * i))
                                 cur_y = int(math.floor(y + grad_y * i))
                                 if cur_x != prev_x or cur_y != prev_y:
                                     # we have moved to the next pixel!
                                     try:
                                         theta_point = theta[y, x]
                                         alpha = theta[int(cur_y), int(cur_x)]
                                         thickness = math.sqrt( (cur_x - x) * (cur_x - x) + (cur_y - y) * (cur_y - y) )

                                         if edges[int(cur_y), int(cur_x)] > 0  :
                                             # found edge,
                                             #print("found edge")


                                             ray.append((int(cur_x), int(cur_y)))
                                             #if abs(thickness -initialthickness)/initialthickness
                                             if abs(abs(alpha - theta_point) - np.pi) < slopeErr  or ( abs(alpha - theta_point) < slopeErr and thickness > 10):
                                             #if math.acos(grad_x * -grad_x_g[cur_y, cur_x] + grad_y * -grad_y_g[cur_y, cur_x]) < np.pi/2.0:


                                             #if abs(alpha + theta_point) -np.pi  < slopeErr and thickness > 5:

                                             #if True: or thickness < initialthickness:
                                                 if  thickness < initialthickness:

                                                     midx = int ((x+cur_x)/2)
                                                     midy = (int (y+cur_y)/2)
                                                     #print theta_visible[y,x]
                                                     if sobelx64f[y,x] >= 0 :
                                                        swt[midy,midx] = theta_visible[y,x]
                                                     else:
                                                        swt[midy,midx] = theta_visible[cur_y,cur_x]

                                                     #swt[midy,midx] = max(theta_visible[int(cur_y),int(cur_x)],theta_visible[y,x])
                                                     rays.append(ray)
                                                     newedges[int(cur_y),int(cur_x)] = 0
                                                     break


                                                 #for (rp_x, rp_y) in ray:
                                                   #  swt[rp_y, rp_x] = min(thickness, swt[rp_y, rp_x])
                                                 #rays.append(ray)
                                             else:
                                                 print "alpha =",alpha, " theta", theta_point," thickness=",thickness

                                                 rays_others.append(ray)

                                             #break

                                             #break
                                         # this is positioned at end to ensure we don't add a point beyond image boundary
                                         #ray.append((int(cur_x), int(cur_y)))
                                     except IndexError:
                                         # reached image boundary
                                         break
                                     except ValueError:
                                         break
                                     prev_x = cur_x
                                     prev_y = cur_y

                 (swt,rays_others) = self.dealWithOtherRays(swt,rays_others,edges,theta_visible,initialthickness*1.1)
                 if self.diagnostics:
                     self.printRays(rays,edges,"edge-experimental.jpg")
                     self.printRays(rays_others,edges,"others-edge-experimental.jpg")
                     cv2.imwrite(self.workingdir+"/"+"swt.jpg",swt)

                 return (img,swt,rays,rays_others)
    @classmethod

    def extendEdge(self,img,edges,theta,distMax): # theta of original image
        # find corners
        # draw line till edge is found or till length 3 * thickness
        q = self.corners(edges,theta,"tmp-corners.jpg")
        edges_color = cv2.cvtColor(edges,cv2.COLOR_GRAY2RGB)
        edges_new =np.copy(edges)
        print len(q)
        for cornerpoints in q:
            for p in cornerpoints:

                i =0
                if theta[p[1],p[0]] > 0:
                    grad_x = np.cos(theta[p[1],p[0]])
                    grad_y = -np.sin(theta[p[1],p[0]])
                else:
                    grad_x = -np.cos(theta[p[1],p[0]])
                    grad_y = np.sin(theta[p[1],p[0]])


                while i < 2000:
                    i += 1
                    x = math.floor(p[0] + i * grad_x)
                    y = math.floor(p[1] + i * grad_y)
                    #x = int(round(x))
                    #y = int(round(x))
                    try:
                        edges_color[y,x] = (0,0,255)
                        edges_new[y,x] = 255

                        dist = np.linalg.norm(np.array((x,y))-np.array((p[0],p[1])))

                        if edges[y,x] > 20 or dist > distMax:


                            break
                        else:
                           #edges[y,x] = 255
                           c=1
                    except IndexError:
                        break

        if self.diagnostics:
            cv2.imwrite(self.workingdir+"/tmp-edges-extended.jpg",edges_color)
            cv2.imwrite(self.workingdir+"/tmp-gray-edges-extended.jpg",edges_new)

        return edges_new
    @classmethod
    def extendLines(self,img,edges,theta):
        # if there are vertical lines then go up and down till you reach the limit
        indices = np.where((edges >0) & ((abs(theta) <=0.5) | (abs(theta - np.pi)<=0.5) | (abs(theta +np.pi)<0.5)))
        #sys.exit()
        for i,y in enumerate(indices[0]):
            k = y
            while k < img.shape[0] and img[k,indices[1][i]] >0  :
                edges[k,indices[1][i]] = 255
                theta[k,indices[1][i]] = theta[y,indices[1][i]]
                k = k + 1
            k = y
            while k>=0 and img[k,indices[1][i]] >0 :
                edges[k,indices[1][i]] = 255
                theta[k,indices[1][i]] = theta[y,indices[1][i]]

                k = k -1
        indices = np.where((edges >0) & ( (abs(theta - np.pi/2)<=0.5) | (abs(theta +np.pi/2)<0.5)))
        #sys.exit()
        for i,x in enumerate(indices[1]):
            k = x
            while k < img.shape[1] and img[indices[0][i],k] >0  :
                edges[indices[0][i],k] = 255
                theta[indices[0][i],k] = theta[indices[0][i],x]
                k = k + 1
            k = x
            while k>=0 and img[indices[0][i],k] >0 :
                edges[indices[0][i],k] = 255
                theta[indices[0][i],k] = theta[indices[0][i],x]

                k = k -1


    @classmethod
    def remVerticalHor(self,img):
        for x in xrange(img.shape[1]):
            for y in xrange(img.shape[0]):
                try:
                    if img[y,x] > 0:
                        img[y,x+1] = 0
                        img[y,x-1] = 0
                        img[y+1,x] = 0
                        img[y-1,x] = 0
                except IndexError:
                    continue
    @classmethod
    def findBestLabelForAppending(self,point,pointslope,labels,distMax,errMax):
        # get distance array for all labels
        point = np.array((point[0],point[1]))
        distarray =[]
        for i,label in enumerate(labels):
            dist1 = np.linalg.norm(point-np.array((label[0][1],label[0][2]))) # from beginning of curve
            dist2 = np.linalg.norm(point-np.array((label[-1][1],label[-1][2]))) # from end of curve
            distarray.append((dist2,i,-1))
            distarray.append((dist1,i,0))
        distarray.sort(compare3)
        connectionElements =[]
        for i,distelement in enumerate(distarray):
            if distelement[0] < distMax and abs(pointslope-labels[distelement[1]][distelement[2]][0])<errMax:
                connectionElements.append((abs(pointslope-labels[distelement[1]][distelement[2]][0]),distelement[1],distelement[2])) # collect all those labels which are candidate for adding
                break
        #connectionElements.sort(compare3)
        if len(connectionElements) > 0:
            return connectionElements[0]
        else:
            return []
    #def printBezierCurves(img,curves):
    @classmethod
    def _trace(self,img,filename,distMax=15,errMax = 5.0):

        labels = []
        noofpoints = 0
        img1 = img
        img = img.copy()
        target = open(filename,'w')
        indices = np.where(img>0)
        for l,y in enumerate(indices[0]):
                #if noofpoints > 25000:
                #    break
                x = indices[1][l]
                found = False
                k = 0
                if img[y,x] >0:
                    noofpoints += 1
                    #print img[y,x]

                    for label in labels:
                        k += 1
                        a = np.array((x,y))
                        b = np.array((label[-1][1],label[-1][2]))
                        c = np.array((label[0][1],label[0][2]))

                        if abs(img[y,x] - label[-1][0]) < errMax and np.linalg.norm(a-b) < distMax:
                            label.append([img[y,x],x,y])
                            found = True
                            #print "appended at the end of label no " + str(k)
                            break
                        elif abs(label[0][0]-img[y,x]) <errMax  and np.linalg.norm(a-c)<distMax:
                             label.insert(0,[img[y,x],x,y])
                             #print "inserted at the beginning of label no " + str(k)

                             found = True
                             break
                    '''
                    connectionElement = findBestLabelForAppending((y,x),img[y,x],labels[:],distMax,errMax);
                    if len(connectionElement) >0:
                        found = True
                        if connectionElement[2] == -1:
                            labels[connectionElement[1]].append([img[y,x],x,y])
                        else:
                            labels[connectionElement[1]].insert(0,[img[y,x],x,y])
                    '''
                    if not found:
                        labels.append([[img[y,x],x,y]])

                        #print "created new label#",len(labels)


        #sort labels
        labels.sort(compare)

        for label in labels:
            if abs((label[0][0]-label[-1][0])/255)<=0.001: # straight line
                print label[0][0],label[0][1],label[0][2]
                for j,point in enumerate(label):
                     if j!=0 and j<len(label)-1:
                        label.pop(j)
                #sys.exit(0)

        i=0

        j = 0
        clr1=0
        clr2=0
        curves = []
        labels.sort(compare1)
        self.stitchLabelsNew(labels,distMax,errMax)
        for i,label1 in enumerate(labels):
            if len(label1) < 10 and len(label1) > 2:
                labels.pop(i)

        for label in labels:
            curve = []
            for point in label:
                curve.append(np.array((point[1],point[2])))
            if len(curve) > 10:
                curves.append(curve)

        plt.gca().invert_yaxis()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim(0,img.shape[1])
        ax.set_ylim(img.shape[0],0)
        for curve in curves:
            #print curve
            curve=np.array(curve)
            #print curve
            target.write("path.move(to:CGPoint(x:"+str(curve[0][0])+",y:"+str(curve[0][1])+"))\n")

            if len(curve)==2: # straight line
                target.write("path.addLine(to:CGPoint(x:"+str(curve[1][0])+",y:"+str(curve[1][1])+"))\n")
                verts=[curve[0],curve[1]]
                codes = [Path.MOVETO,Path.LINETO]
                #path = Path(verts, codes)
                #patch = patches.PathPatch(path,facecolor=None ,lw=2)
                ax.add_patch(patch)

            else:
                beziers = fitCurve(curve,15.0)
                #beziers = fitCurve(currentContour, 10)
                #print(len(beziers))
                for bezier in beziers:
                    #print bezier
                    #print bezier[0][0]
                    if not (math.isnan(bezier[1][0])) and not (math.isnan(bezier[2][0])):
                        #print("\n")
                        target.write("path.addCurve(to:CGPoint(x:"+str(bezier[3][0])+",y:"+str(bezier[3][1])+"),controlPoint1:CGPoint(x:"+str(bezier[1][0])+",y:"+str(bezier[1][1])+"),controlPoint2:CGPoint(x:"+str(bezier[2][0])+",y:"+str(bezier[2][1])+"))\n")
                        verts=[bezier[0],bezier[1],bezier[2],bezier[3]]
                        codes = [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4]
                        path = Path(verts, codes)
                        patch = patches.PathPatch(path,facecolor='none', lw=1)
                        #xs, ys = zip(*verts)
                        #ax.plot(xs, ys, 'x--', lw=2, color='black', ms=10)

                        ax.add_patch(patch)



        plt.show()
        fig.savefig(self.workingdir+"/"+"bezier.png")
        j = 0
        img1 = img.copy()
        img1 = cv2.cvtColor(img1,cv2.COLOR_GRAY2RGB)
        img2 = img1.copy()
        make_sure_path_exists(self.workingdir+"/labels")
        for label in labels:
            img2[:] = (0,0,0)

            i = 0
            j += 1
            #img1 = np.copy(img3)

            #img1 = img3
            for points in label:
                if i == 0 or i == len(label)-1:
                    font = cv2.FONT_HERSHEY_PLAIN
                    cv2.putText(img1,str(j)+","+str(i),(points[1],points[2]), font, 0.6,(255,255,255),1,cv2.LINE_AA,False)

                if j%2 ==0:
                    clr1 =0
                    clr2 = 255
                else:
                    clr1 = 255
                    clr2 = 0
                img1[points[2],points[1]] = (clr1 ,0,clr2)
                img2[points[2],points[1]] = (clr1 ,0,clr2)

                i +=1

            cv2.imwrite(self.workingdir+"/labels/"+str(j)+".png",img1)
        cv2.imwrite(self.workingdir+"/labels.png",img1)
        return curves
    @classmethod

    def stitchLabelsNew(self,labels,distMax,errMax):
        #labels = np.uint64(labels)

        for i,label in enumerate(labels):
            for j,point in enumerate(label):
               label[j] = [int(k) for k in label[j]]

        k = 0
        while k <len(labels):
            i = k
            j = i+1
            while j <len(labels) and j >i:
                a = np.array((labels[i][0][1],labels[i][0][2]))
                b = np.array((labels[i][-1][1],labels[i][-1][2]))
                c = np.array((labels[j][0][1],labels[j][0][2]))
                d = np.array((labels[j][-1][1],labels[j][-1][2]))
                distances = {}
                distances['ac'] = np.linalg.norm(a-c)
                distances['ad'] = np.linalg.norm(a-d)
                distances['bc'] = np.linalg.norm(b-c)
                distances['bd'] = np.linalg.norm(b-d)
                slopediff={}
                slopediff['ac'] = abs(labels[i][0][0]-labels[j][0][0])
                slopediff['ad'] = abs(labels[i][0][0]-labels[j][-1][0])
                slopediff['bc'] = abs(labels[i][-1][0]-labels[j][0][0])
                slopediff['bd'] = abs(labels[i][-1][0]-labels[j][-1][0])
                sorted_slope = sorted(slopediff, key=slopediff.get)
                if True or slopediff[sorted_slope[0]]  < errMax:
                    if distances[sorted_slope[0]] < distMax:
                        if sorted_slope[0] == 'ac':
                            labels[i].reverse()
                            labels[i] = labels[i] + labels[j]
                            labels.pop(j)
                            k = -1
                            j = i

                        elif sorted_slope[0] == 'ad':
                            labels[i] = labels[j] + labels[i]
                            labels.pop(j)
                            k =-1
                            j = i
                        elif sorted_slope[0] == 'bc':
                            labels[i] = labels[i] + labels[j]
                            labels.pop(j)
                            k = -1
                            j = i
                        elif sorted_slope[0] == 'bd':
                            labels[j].reverse()
                            labels[i] = labels[i] + labels[j]
                            labels.pop(j)
                            k = -1
                            j = i
                j = j + 1

            k += 1

    @classmethod

    def stitchLabels(self,labels,distMax,errMax):
        #stitch labels
        k = 0
        while k <len(labels):
            i = k
            j = i+1
            while j <len(labels) and j >i:
                a = np.array((labels[i][0][1],labels[i][0][2]))
                b = np.array((labels[i][-1][1],labels[i][-1][2]))
                c = np.array((labels[j][0][1],labels[j][0][2]))
                d = np.array((labels[j][-1][1],labels[j][-1][2]))
                distances = {}
                distances['ac'] = np.linalg.norm(a-c)
                distances['ad'] = np.linalg.norm(a-d)
                distances['bc'] = np.linalg.norm(b-c)
                distances['bd'] = np.linalg.norm(b-d)
                slopediff={}
                slopediff['ac'] = abs(labels[i][0][0]-labels[j][0][0])
                slopediff['ad'] = abs(labels[i][0][0]-labels[j][-1][0])
                slopediff['bc'] = abs(labels[i][-1][0]-labels[j][0][0])
                slopediff['bd'] = abs(labels[i][-1][0]-labels[j][-1][0])
                sorted_slope = sorted(slopediff, key=slopediff.get)
                if slopediff[sorted_slope[0]]  < errMax:
                    if distances[sorted_slope[0]] < distMax:
                        if sorted_slope[0] == 'ac':
                            labels[i].reverse()
                            labels[i] = labels[i] + labels[j]
                            labels.pop(j)
                            k = -1
                            j = i

                        elif sorted_slope[0] == 'ad':
                            labels[i] = labels[j] + labels[i]
                            labels.pop(j)
                            k =-1
                            j = i
                        elif sorted_slope[0] == 'bc':
                            labels[i] = labels[i] + labels[j]
                            labels.pop(j)
                            k = -1
                            j = i
                        elif sorted_slope[0] == 'bd':
                            labels[j].reverse()
                            labels[i] = labels[i] + labels[j]
                            labels.pop(j)
                            k = -1
                            j = i
                j = j + 1

            k += 1

    @classmethod
    def moveAlongEdge(self,img,filename):
         swt = img.copy()

         img2 = np.zeros(img.shape)
         img1 = np.uint8(swt)
         img1 = cv2.cvtColor(img1,cv2.COLOR_GRAY2RGB) #get color image
         #groupOfPoints =[]
         label = 0
         curves = []

         #sys.exit()
         for y in xrange(img.shape[0]):
             for x in xrange(img.shape[1]):

                 if  img[y,x] > 0 and img2[y,x]<255 : #found an edge
                    points = [(x,y)]
                    img2[y,x] = 255
                    #if (label % 10 ==0):
                       #cv2.putText(img1,str(label),(x+4,y),cv2.FONT_HERSHEY_SIMPLEX,0.4,(255,0,0),1,cv2.LINE_AA)
                    i = 1 #delta_y
                    angle = img[y,x]/255 *2*np.pi -np.pi

                    grad_x = np.sin(angle)
                    grad_y = -np.cos(angle)
                    prev_y = y
                    prev_x = x
                    count = 0 # no of edge points traversed
                    located = True
                    #while cur_y>0 and cur_x >0  and i<edges.shape[0] and cur_y < edges.shape[0] and cur_x < edges.shape[1]:
                    while prev_y < img.shape[0] and prev_x < img.shape[1] and located and i < 500:

                        angle = (img[prev_y,prev_x]-10)/255 *2*np.pi -np.pi

                        grad_x = int(round(-1*np.cos(angle)))
                        grad_y = int(round(np.sin(angle)))
                        #print grad_x,':',grad_y
                        cur_x = prev_x +  grad_x
                        cur_y = prev_y + grad_y
                        if (cur_x != prev_x or cur_y != prev_y) and cur_y> 0 and cur_x>0 and cur_y < img.shape[0] and cur_x < img.shape[1] and img2[cur_y,cur_x] < 255:
                            # we have moved to the next pixel!

                            #print cur_x,'-',cur_y,' angle=',angle
                            img1[cur_y,cur_x] = (0,0,255)

                            try:


                                   located = False

                                   if cur_y>0 and cur_x >0 and cur_y < img.shape[0] and cur_x < img.shape[1]  and img[cur_y,cur_x] >0:
                                      img1[cur_y,cur_x] = (0,255,0)
                                      img2[cur_y,cur_x] = 255 #label visited
                                      located = True
                                      points.append((cur_x,cur_y))
                                      count += 1
                                      prev_x = cur_x
                                      prev_y = cur_y
                                   elif cur_y>0 and cur_x >0 and cur_y < img.shape[0] and cur_x < img.shape[1]:
                                       # try locating nearby edge
                                       for k in xrange(-10,10,1):
                                           for t in xrange(-10,10,1):

                                              check_x = cur_x + k
                                              check_y = cur_y + t
                                              if check_x > 0 and check_x < img.shape[1] and check_y > 0 and check_y < img.shape[1]  and img[check_y,check_x] > 0:
                                                  img1[check_y,check_x] = (0,255,0)
                                                  points.append((check_x,check_y))
                                                  #print "checked point at label #",label
                                                  img2 [check_y,check_x] = 255 #label visited
                                                  located = True
                                                  count += 1
                                                  prev_x = check_x
                                                  prev_y = check_y
                                                  break
                                       #break
                            except IndexError:
                                break
                            except ValueError:
                                continue
                            except:
                                continue
                        i = i + 1
                        if not located:
                                print "could not locate"
                                img1[cur_y,cur_x] = (0,0,255)

                    #img[y,x] = np.Infinity

                    if count > 0:
                            label += 1

                            curves.append(points)
                            #if label % 5 ==0:
                             #  cv2.putText(img1,str(label),(x+4,y),cv2.FONT_HERSHEY_SIMPLEX,0.4,(255,0,0),1,cv2.LINE_AA)


         cv2.imwrite(filename,img2)
         #img1 = cv2.imread("swt-exp.jpg")
         j=0
         print len(curves)
         for curve in curves:
             i = 0
             j += 1
             #img1 = np.copy(img3)

             #img1 = img3
             for points in curve:
                 if i == 0:
                     font = cv2.FONT_HERSHEY_PLAIN
                     cv2.putText(img1,str(j),(int(round(points[0])),int(round(points[1]))), font, 0.8,(255,255,255),1,cv2.LINE_AA,False)

                 if j%2 ==0:
                     clr1 =0
                     clr2 = 255
                 else:
                     clr1 = 255
                     clr2 = 0
                 img1[points[1],points[0]] = (clr1 ,0,clr2)
                 i +=1

             #cv2.imwrite("labels/"+str(j)+".jpg",img1)
         cv2.imwrite("labels.png",img1) # rival of _trace - would use in future
    @classmethod
    def findShirorekhaHeight(self,img,edges):
        #let us find horizontal projection
        vert = []
        imgbeginy = 0
        imgendy = 0
        imgbeginx = 0
        imgendx = 0
        lineheight=[]
        for j in xrange(img.shape[0]):
            #vert.append(
            c = np.count_nonzero(edges[j,:])
            lineheight.append((j,c))
            if c > 10:
                imgendy = j
                if imgbeginy == 0:
                   imgbeginy = j

        for j in xrange(img.shape[1]):
            #vert.append(
            c = np.count_nonzero(edges[:,j])
            if c > 10:
                imgendx = j
                if imgbeginx == 0:
                   imgbeginx = j
        lineheight.sort(compare2)
        shirorekhaheight = abs(lineheight[0][0]-lineheight[1][0])
        return shirorekhaheight
    @classmethod
    def corners(self,img,theta,filename):
         diffMin = 0.2
         #(img2, sobelx64f, sobely64f, theta) = FindStrokes._create_derivative(img)

         #img1 = np.zeros((img.shape[0],img.shape[1],1),np.uint8)
         img1 = cv2.cvtColor(img,cv2.COLOR_GRAY2RGB) #get color image
         p = []
         for x in xrange(img.shape[1]):
             for y in xrange(img.shape[0]):
                 neighbours = [[x+1,y],[x+1,y+1],[x+1,y-1],[x,y+1],[x,y-1],[x-1,y-1],[x-1,y],[x+1,y+1]]
                 cornerfound = False
                 if img[y,x] > 20 :
                     count = 0
                     q =[]
                     for neighbor in neighbours:
                         if abs(theta[neighbor[1],neighbor[0]] - theta[y,x]) > diffMin and img[neighbor[1],neighbor[0]]>20:
                             count  += 1
                             q.append((neighbor[0],neighbor[1]))

                     if count >=2  :
                          img1[y,x] = (0,0,255)
                          p.append(q)
         cv2.imwrite(self.workingdir+"/"+filename,img1)
         return p
