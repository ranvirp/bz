
import numpy as np
import cv2
#import time
#import math
#import operator
from matplotlib import pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.collections as collections
from ttfquery import describe
from ttfquery import glyph,glyphquery
from fitCurves import *

class Check(object):
    def edgeLink(self,points,slopes,nbrsize = 3): # Should consist of slopeValues
          pointsimage = np.zeros(slopes.shape)
          for point in points:
              pointsimage[point[1],point[0]] = 255
          count = len(points)
          pointscovered = 0
          joins = []
          for point in points:
              y = point[1]
              x = point[0]
              if pointsimage[y,x] ==0:
                  continue
              likely = []
              likely.append([x,y])
              pointscovered += 1
              found = True
              while found and pointscovered < count:
                  found = False
                  for i in xrange(-1*nbrsize/2,nbrsize/2):
                      if found:
                          break
                      for j in xrange(-1*nbrsize/2,nbrsize/2):
                          if  pointsimage[y+j,x+i] >0:
                              if i==0:
                                  slope = np.pi/2
                              else:
                                  slope = np.arctan2(j,i)
                              pointsimage[y+j,x+i] = 0
                              pointscovered += 1
                              print slope,slopes[y,x]
                              if abs(slope - slopes[y,x]) < 0.01:
                                   likely.append([x+i,y+j])
                                   y = y +j
                                   x = x+i
                                   found = True

                                   break
              joins.append(likely)
          return joins
    def offsetEdge(self,edge,offset,ts=None):
        if len(edge) ==2:
            return self.lineOffset(edge,offset)
        elif len(edge) == 3: #quadCurve
            return self.bezierQuadOffset(edge,offset,ts)
        elif len(edge) == 4:
            return self.bezierCubicOffset(edge,offset,ts)
    def bezierCubicOffset(self,cubic,offset,ts=None):
        p1 = cubic[0]
        pc1 = cubic[1]
        pc2 = cubic[2]
        p2 = cubic[3]
        # p = (1-t)**3 * p1 + 3*t*(1-t)**2 * pc1 + 3*t**2*(1-t)*pc2 + t**3 * p2
        # pdash = 3*(pc1-p1)*(1-t)**2 + 6*(pc2-pc1)*t*(1-t) + 3*(p2-pc2)*t**2
        points = []
        if ts == None:
            ts= [step/100.0 for step in xrange(101)]
        for t in ts: # 100 points on the curve
            p = (1-t)**3 * p1 + 3*t*(1-t)**2 * pc1 + 3*t**2*(1-t)*pc2 + t**3 * p2
            pdash = 3*(pc1-p1)*(1-t)**2 + 6*(pc2-pc1)*t*(1-t) + 3*(p2-pc2)*t**2
            pdashnormal = (pdash[1],-1*pdash[0])
            costheta = pdashnormal[0]* 1/np.linalg.norm(pdashnormal)
            sintheta = pdashnormal[1]* 1/np.linalg.norm(pdashnormal)
            pnew = p + np.array([costheta,sintheta]) * offset
            points.append(pnew)
        #print points
        if len(points) >1:
            beziers = fitCurve(points,10.0)
            return beziers
        else:
            return [points]

    def lineOffset(self,line,offset,ts=None):
        # p = p1*(1-t) + p2 *t
        # pdash = -p1 + p2
        p1 = line[0]
        p2 = line[1]
        if ts != None:
            p1 = p1*(1-ts[0]) + p2 * ts[0]
            p2 = p1*(1-ts[-1]) + p2 * ts[-1]
        pdash = -p1 + p2
        pdashnormal = (pdash[1],-1*pdash[0])
        costheta = pdashnormal[0]* 1/np.linalg.norm(pdashnormal)
        sintheta = pdashnormal[1]* 1/np.linalg.norm(pdashnormal)
        p1 = p1 + np.array([costheta,sintheta]) * offset
        p2 = p2 + np.array([costheta,sintheta]) * offset

        return [[p1,p2]]



    def bezierQuadOffset(self,quad,offset,ts=None):
        #print quad
        p1 = quad[0]
        pc = quad[1]
        p2 = quad[2]
        # p = (1-t)**2 * p1 + 2 *t*(1-t)*pc + t**2*p2
        #pdash = 2*(t*((p2-pc)-(pc-p1)) + (pc-p1))

        points = []
        if ts == None:
            ts= [step/100.0 for step in xrange(101)]
        for t in ts: # 100 points on the curve
            p = (1-t)**2 * p1 + 2 *t*(1-t)*pc + t**2*p2
            pdash = 2*(t*((p2-pc)-(pc-p1)) + (pc-p1))
            pdashnormal = (pdash[1],-1*pdash[0])
            costheta = pdashnormal[0]* 1/np.linalg.norm(pdashnormal)
            sintheta = pdashnormal[1]* 1/np.linalg.norm(pdashnormal)
            pnew = p + np.array([costheta,sintheta]) * offset
            points.append(pnew)
        if len(points) > 1:
            beziers = fitCurve(points,10.0)

            return beziers
        else:
            return [points]
    def offsetQuadBezier(self,quad,offsets):
        fig =plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_xlim(50,200)
        ax1.set_ylim(50,200)
        code1 = [Path.MOVETO,Path.CURVE3,Path.CURVE3]
        path1 = Path(quad,code1)
        patch1 = patches.PathPatch(path1,facecolor='none',lw=1)
        ax1.add_patch(patch1)
        for offset in offsets:
            points = self.bezierQuadOffset(quad,offset)
            beziers = fitCurve(points,10.0)
            print len(beziers)
            for bezier in beziers:
                verts=[bezier[0],bezier[1],bezier[2],bezier[3]]
                codes = [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4]
                path = Path(verts, codes)
                patch = patches.PathPatch(path,facecolor='none', lw=1)
                ax1.add_patch(patch)
        plt.show(block=False)
    def offsetCubicBezier(self,cubic,offsets):
        fig =plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_xlim(cubic[0][0],cubic[3][0])
        ax1.set_ylim(cubic[0][1],cubic[3][1])
        code1 = [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4]
        path1 = Path(cubic,code1)
        patch1 = patches.PathPatch(path1,facecolor='none',lw=1)
        ax1.add_patch(patch1)
        for offset in offsets:
            points = self.bezierCubicOffset(cubic,offset)
            beziers = fitCurve(points,10.0)
            print len(beziers)
            for bezier in beziers:
                verts=[bezier[0],bezier[1],bezier[2],bezier[3]]
                codes = [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4]
                path = Path(verts, codes)
                patch = patches.PathPatch(path,facecolor='none', lw=1)
                ax1.add_patch(patch)
        plt.show(block=False)

    def compQuadOffsets(self,quad,offsets):
        #fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(8, 4), sharex=True, sharey=True)
        fig =plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_xlim(50,200)
        ax1.set_ylim(50,200)
        code1 = [Path.MOVETO,Path.CURVE3,Path.CURVE3]
        path1 = Path(quad,code1)
        patch1 = patches.PathPatch(path1,facecolor='none',lw=1)
        ax1.add_patch(patch1)

        for offset in offsets:
            points = self.bezierQuadOffset(quad,offset)
            code2 = []
            for pt in points:
                code2.append(Path.LINETO)
            code2[0] = Path.MOVETO
            path2 = Path(points,code2)
            patch2 = patches.PathPatch(path2,facecolor='none',lw=1)
            ax1.add_patch(patch2)
        plt.show(block=False)


        return
    def glFromChar(self,unicode=None,fontfile="DevanagariSangamMN.ttf"):
        if unicode==None:
            return ''
        else:
            font = describe.openFont(fontfile)
            return glyphquery.glyphName(font, unichr(unicode))





    def fontOutline(self,gl,fontfile="DevanagariSangamMN.ttf",returnobj=False):
        font = describe.openFont(fontfile)
        gl1 = glyph.Glyph(gl)
        gl1.compile(font,steps=3)
        outlines = gl1.outlines
        contours = gl1.contours
        return self.edgesFromOutline(outlines,contours,returnobj)
    def edgesFromOutline(self,outlines,contours,returnobj=False):

        codes = []
        edges = []
        offcurve = 0

        points = []
        vert_edges={}
        edge_index = 0
        pt_index = 0
        cnts = []
        for i,contour in enumerate(contours):

           start = True
           prevzero = False
           prevone = False
           edge = []
           cnt = []
           for j,point in enumerate(contour):

               if point[1] == 1: # on-curve
                   edge.append(np.array(point[0]).astype(float))
                   prevzero = False
                   if j!=len(contour)-1:
                       points.append(np.array(point[0]))
                       cnt.append(pt_index)
                       pt_index += 1

                   if prevone:
                       if len(edge) ==2:
                           edge.insert(1,(edge[1]+edge[0])/2.0)

                       edges.append(edge)
                       edge_index += 1
                       edge = []
                       edge.append(np.array(point[0]).astype(float))

                       if not pt_index in vert_edges:
                            vert_edges[pt_index] = []
                       vert_edges[pt_index].append(edge_index)


                   start = False
                   prevone = True
               elif point[1] == 0:
                   if prevzero:
                       missing_point = (np.array(point[0]).astype(float) + np.array(prevpoint[0]).astype(float))/2
                       points.append(missing_point)
                       cnt.append(pt_index)

                       edge.append(missing_point)
                       edges.append(edge)
                       if not pt_index in vert_edges:
                           vert_edges[pt_index] = []
                       vert_edges[pt_index].append(edge_index)
                       edge_index += 1
                       edge = []
                       vert_edges[pt_index].append(edge_index)

                       edge.append(missing_point)
                       pt_index += 1

                   edge.append(np.array(point[0]).astype(float))
                   prevpoint = point
                   prevzero = True
           cnts.append(cnt)
        if returnobj:
            return edges,outlines,contours,points
        else:
            return edges,points,vert_edges,cnts
    def fontSkeleton(self,edgegroups):
        bzs=[]

        for k,edgegroup in enumerate(edgegroups):
            prevedge = -2
            edges = edgegroup[0]
            offset = edgegroup[1]
            for j,edge in enumerate(edges):
                if prevedge > -2:
                    bz1 = self.offsetEdge(edges[prevedge],offset)
                    #print bz1,offset
                    #print np.linalg.norm(bz1[0]-edges[prevedge][0])
                    offset = np.linalg.norm(bz1[0][-1]-edge[0])
                    #print offset
                prevedge = k
                if offset >0:
                    bz = self.offsetEdge(edge,offset)
                else:
                    bz = [edge]
                bzs += bz
        return bzs




    def printEdgesN(self,edges,offset=0,numbering=True,labeling=False):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        verts = []
        codes = []
        minx = maxx = miny = maxy = 0.0
        patchess = []
        prevedge = -2
        colours = ['r','b','g']
        mainoffset = offset
        for k,edge in enumerate(edges):
            #print len(edge)
            edgeinfo = [colours[k%3]];offset = mainoffset;label=''

            if isinstance(edge,(list,np.ndarray)) and isinstance(edge[0],(list,np.ndarray)) and  isinstance(edge[0][0],(list,np.ndarray)):
                if len(edge)>1:
                    edgeinfo = edge[1:]
                    #print edgeinfo

                edge = edge[0]
                if len(edgeinfo) >=2:
                    offset = edgeinfo[1]
                    label  = edgeinfo[2]
                    #print label
            #    print edgeinfo
            #print edge
            if len(edge) ==0:
                continue
            #print prevedge
            if prevedge > -2:
                bz1 = self.offsetEdge(edges[prevedge],offset)
                #print bz1,offset
                #print np.linalg.norm(bz1[0]-edges[prevedge][0])
                offset = np.linalg.norm(bz1[0][-1]-edge[0])
                #print offset
            prevedge = k
            if offset >0:
                bz = self.offsetEdge(edge,offset)
                try:
                    if len(bz) > 0:
                        edge = bz[0]
                except:
                        print edge,bz
            if numbering:
                    ax.text(edge[0][0]+0.01, edge[0][1]-0.01, str(k), style='italic',bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
            if labeling:
                    ax.text(edge[0][0]-0.01, edge[0][1]-0.01, label, style='italic',bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})


            for pt in edge:
                minx = min(minx,pt[0])
                miny = min(miny,pt[1])
                maxx = max(maxx,pt[0])
                maxy = max(maxy,pt[1])
                if len(edge) >1:
                    verts.append(pt)

            if len(edge) == 1:# a point
                ax.plot(edge[0][0],edge[0][1],'x-',edgecolor=edgeinfo[0])

            elif len(edge) == 2:
                #verts.append(edge[0])
                #verts.append(edge[1])
                codes.append(Path.MOVETO)
                codes.append(Path.LINETO)
            elif len(edge) == 3:
                #verts.append(edge[0])
                #verts.append(edge[1])
                #verts.append(edge[2])
                codes.append(Path.MOVETO)
                codes.append(Path.CURVE3)
                codes.append(Path.CURVE3)
            elif len(edge) == 4:
                #verts.append(edge[0])
                #verts.append(edge[1])
                #verts.append(edge[2])
                #verts.append(edge[3])
                codes.append(Path.MOVETO)
                codes.append(Path.CURVE4)
                codes.append(Path.CURVE4)
                codes.append(Path.CURVE4)
        #print len(verts)
        #print len(codes)
            if len(edge) >1:
                path = Path(verts,codes)
                #print edgeinfo[0]
                #edgeinfo[0]='k'
                patch = patches.PathPatch(path,facecolor='none',lw=1,ec=edgeinfo[0])
                #print edgeinfo[0]
                #patch.setEdgeColor(edgeinfo[0])
                patchess.append(patch)
        p = collections.PatchCollection(patchess) #, cmap=matplotlib.cm.jet)
        ax.add_collection(p)

        ax.set_xlim(minx-30,maxx+30)
        ax.set_ylim(miny-30,maxy+30)
        plt.show(block=False)







        return

    def printEdges(self,edges,offset=0,numbering=True,show=True,name='dn_ra',labeling=False):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        verts = []
        codes = []
        minx = maxx = miny = maxy = 0.0
        patchess = []
        colours = ['g','b','r']
        mainoffset = offset;label=''
        for k,edge in enumerate(edges):
            #print len(edge)
            verts = []
            codes = []

            if len(edge) == 0:
                continue
            edgeinfo = [colours[k%3],mainoffset,'']
            try:
                if isinstance(edge,(list,np.ndarray)) and isinstance(edge[0],(list,np.ndarray)) and  isinstance(edge[0][0],(list,np.ndarray)):
                    if len(edge)>1:
                        edgeinfo = edge[1:]
                    edge = edge[0]
            except:
                edge = edge[0]
            #    print edgeinfo
            #print edge
            if len(edge) ==0:
                continue

            if offset >0:
                bz = self.offsetEdge(edge,offset)
                #print np.linalg.norm(bz[0][0] - edge[0])
                try:
                    if len(bz) > 0:
                        edge = bz[0]
                except:
                        print edge,bz
            if numbering:
                    ax.text(edge[0][0], edge[0][1], str(k), style='italic',bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
            if labeling and edgeinfo[2]!=label:
                    label = edgeinfo[2]
                    ax.text(edge[0][0]-0.01, edge[0][1]-0.01, edgeinfo[2], style='italic',bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})

            for pt in edge:
                minx = min(minx,pt[0])
                miny = min(miny,pt[1])
                maxx = max(maxx,pt[0])
                maxy = max(maxy,pt[1])
                if len(edge) >1:
                    verts.append(pt)

            if len(edge) == 1:# a point
                ax.plot(edge[0][0],edge[0][1],'x-',color=edgeinfo[0])

            elif len(edge) == 2:
                #verts.append(edge[0])
                #verts.append(edge[1])
                codes.append(Path.MOVETO)
                codes.append(Path.LINETO)
            elif len(edge) == 3:
                #verts.append(edge[0])
                #verts.append(edge[1])
                #verts.append(edge[2])
                codes.append(Path.MOVETO)
                codes.append(Path.CURVE3)
                codes.append(Path.CURVE3)
            elif len(edge) == 4:
                #verts.append(edge[0])
                #verts.append(edge[1])
                #verts.append(edge[2])
                #verts.append(edge[3])
                codes.append(Path.MOVETO)
                codes.append(Path.CURVE4)
                codes.append(Path.CURVE4)
                codes.append(Path.CURVE4)
        #print len(verts)
        #print len(codes)

            if len(edge) >1:
                path = Path(verts,codes)
                #print edgeinfo[0]
                patch = patches.PathPatch(path,facecolor='none',fill=False,lw=3,ec=edgeinfo[0])
                #print edgeinfo[0]
                patch.set_ec(edgeinfo[0])
                #patchess.append(patch)
                ax.add_patch(patch)
        #p = collections.PatchCollection(patchess) #, cmap=matplotlib.cm.jet)
        #ax.add_collection(p)

        ax.set_xlim(minx-30,maxx+30)
        ax.set_ylim(miny-30,maxy+30)
        if show==True:
            plt.show(block=False)
        else:
            fig.save(name+'.png')







        return





    def showGlyph(self,gl,fontfile="DevanagariSangamMN.ttf"):
        font = describe.openFont(fontfile)
        gl1 = glyph.Glyph(gl)
        gl1.compile(font,steps=3)
        codes = []
        verts=[]
        offcurve = 0
        for contour in gl1.contours:
           for point in contour:
               if point[1] == 1: # on-curve

                   verts.append(point[0])
                   if offcurve == 0:
                      codes.append(Path.LINETO)
                   else:
                       codes.append(Path.CURVE4)
                       offcurve = 0

               else:
                   offcurve += 1
                   verts.append(point[0])
                   codes.append(Path.CURVE4)
        codes[0] = Path.MOVETO
        codes.pop(-1)
        verts.pop(-1)
        fig =plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim(-50,1200)
        ax.set_ylim(-50,1200)
        path = Path(verts,codes)
        patch = patches.PathPatch(path,facecolor='none',lw=1)
        ax.add_patch(patch)
        print verts
        plt.show()

    def orderedEdge(self,img):
        sobelx = cv2.Sobel(img,cv2.CV_64F,1,0,ksize=3)
        sobely = cv2.Sobel(img,cv2.CV_64F,0,1,ksize=3)
        theta = np.arctan2(sobelx,-1*sobely)
        printed = np.zeros(img.shape)
        printed[:] = 255
        printed = np.uint8(printed)
        printed = cv2.cvtColor(printed,cv2.COLOR_GRAY2RGB)
        visited = np.zeros(img.shape)
        ordered_edges = []
        indices = np.where(((sobelx>0) | (sobely>0)) )
        count = 1
        #print indices
        for i,y in enumerate(indices[0]):
            x = indices[1][i]
            if visited[y,x] >0:
                continue
            prev_x = x
            prev_y = y
            cur_x = -1
            cur_y = -1
            ordered_edge = []

            while x !=cur_x or y != cur_y:
                cur_x = prev_x + 1 * cos(theta[y,x])
                cur_y = prev_y + 1 * sin(theta[y,x])
                cur_x = int(round(cur_x))
                cur_y = int(round(cur_y))
                if sobely[cur_y,cur_x] > 0 or sobelx[cur_y,cur_x] > 0:
                    ordered_edge.append([cur_x,cur_y])
                    visited[cur_y,cur_x] = 255
                    printed[cur_y,cur_x] = (0,0,count)
                    #print x,y,cur_x,cur_y
                    #print len(ordered_edge)
                    prev_x = cur_x
                    prev_y = cur_y
                else:
                    ordered_edges.append(ordered_edge)
                    count += 1
                    break
        cv2.imwrite("ordered-edges.jpg",printed)
        cv2.imwrite("ordered-edges-visited.jpg",visited)


        return ordered_edges


    def correctSize(self,img):
        try:
            indices = np.where(img > 0)
            leftx = min(indices[1])
            rightx = max(indices[1])
            topy = min(indices[0])
            bottomy = max(indices[0])
            img1 = np.zeros((bottomy-topy+5,rightx-leftx+5))
            img1[:] = 0
            img1[:bottomy-topy,:rightx-leftx] = img[topy:bottomy,leftx:rightx]
            return img1
        except:
            return img

    def drawRays(self,theta,img):
         indices = np.where(theta > 0)
         theta_c = np.uint8(theta)
         theta_color = cv2.cvtColor(theta_c,cv2.COLOR_GRAY2RGB)
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
    def findStrokeWidth(self,img):
        sobelx = cv2.Sobel(img,cv2.CV_64F,1,0,ksize=3)
        sobely = cv2.Sobel(img,cv2.CV_64F,1,0,ksize=3)
        theta = np.arctan2(sobely,sobelx)
        #theta = theta + np.pi
        edge = cv2.Canny(img,40,255)
        edge_copy = cv2.cvtColor(edge,cv2.COLOR_GRAY2RGB)
        indices = np.where(edge > 0)
        raylengths = []
        for i,y in enumerate(indices[0]):
            x = indices[1][i]
            ray_length = 1
            found = False
            while ray_length < 100:
                try:
                    cur_x = x + ray_length * cos(theta[y,x])
                    cur_y = y + ray_length * sin(theta[y,x])
                    cur_x = int(round(cur_x))
                    cur_y = int(round(cur_y))
                    print x,y,cur_x,cur_y,theta[y,x]
                    if edge[cur_y,cur_x] > 0 and ray_length > 5:
                        found = True
                        print "found"
                        cv2.circle(edge_copy,(x,y),2,(0,0,255))
                        cv2.circle(edge_copy,(cur_x,cur_y),2,(255,0,0))

                        cv2.line(edge_copy,(x,y),(cur_x,cur_y),(0,255,0))

                        raylengths.append(ray_length)
                        break
                    else:
                        ray_length += 1
                except:
                    break
        cv2.imwrite("found-edge.jpg",edge_copy)

        if len(raylengths) ==0:
               return -1
        else:
            return np.median(raylengths)









    def swtNew(self,img,maxlength = 20,errMax=5):
      #edge = cv2.Canny(img,40,255)
      edge = self.boundary(img)
      sobelx64f = cv2.Sobel(img,cv2.CV_64F,1,0,ksize=3)
      sobely64f = cv2.Sobel(img,cv2.CV_64F,0,1,ksize=3)
      theta = np.arctan2(sobely64f,sobelx64f)
      visited = np.zeros(edge.shape)

      #theta = theta1.astype('int64')
      #plt.imshow(theta)
      #plt.show()
      indices = np.where(edge >0)
      #print indices.shape
      edge_color = cv2.cvtColor(edge,cv2.COLOR_GRAY2RGB)

      swt = np.zeros(img.shape)
      swt[:] = 0
      #swt = np.uint8(swt)
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
          if visited[y,x] > 0 :
              continue
          radian = float(theta[y,x])
          edge_color[y,x] = [0,0,255]
          ray = []
          #ray.append((x,y))
          delta_x = maxlength * np.cos(radian)
          delta_y =  maxlength * np.sin(radian)
          cur_x = int(round(x + delta_x))
          cur_y = int(round(y + delta_y))
          ray_lenth = maxlength
          cur_x_1 = cur_x
          cur_y_1 = cur_y
          try:
              while edge[cur_y_1,cur_x_1] ==0:
                  ray_length += 1

                  delta_x_1 = ray_length * np.cos(radian)
                  delta_y_1 =  ray_length * np.sin(radian)
                  cur_x_1 = int(round(x + delta_x_1))
                  cur_y_1 = int(round(y + delta_y_1))


             # if img[cur_y,cur_x] ==0:
            #      continue
              print ray_length
              if ray_length < 3 * maxlength:
                  if theta[y,x] == 0:
                         swt[cur_y,cur_x] = np.pi
                  else:
                       swt[cur_y,cur_x] = theta[y,x]
                  ray.append((x,y))
                  ray.append((cur_x,cur_y))
                  rays.append(ray)
                  visited[cur_y_1,cur_x_1] = 255


          except:
              continue

      self.printRays(rays,edge,"test2.jpg")
      cv2.imwrite('edge-color.jpg',edge_color)
      #self.group(swt)
      return swt,rays





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
      #self.group(swt)
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
         font = cv2.FONT_HERSHEY_PLAIN

         img_copy = cv2.cvtColor(img,cv2.COLOR_GRAY2RGB)
         img_copy[:] = (0,0,0)
         for i,curve in enumerate(curves):
             for k,point in enumerate(curve):
                  if k == 0:
                         cv2.putText(img_copy,str(i),(point[0] +2 ,point[1]), font, 0.6,(255,255,255),1,cv2.LINE_AA,False)
                  img_copy[point[1],point[0]] = (0,255,0)
         cv2.imwrite("curves.jpg",img_copy)
    def skeleton(self,img,kernel=None):
        #skeleton = union of all skeletons
        # each skeleton = (A eroded k times by B) and then negative by (A eroded k times by B) opening by B
        # opening is dilation followed by erosion
        if kernel == None:
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
        while len(indices)>0 and noofiterations <50:
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
                    if np.linalg.norm(np.array([x,y])-np.array([prev_x,prev_y])) > 10 :
                    #or abs(img[y,x]-img[prev_y,prev_x])>5:
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
        #self._printBezier(curves,img)
        #self.printCurves(img,[curves[0]])
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




    def _printBezier(self,curves):
        plt.gca().invert_yaxis()
        maxy = 0
        maxx = 0
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for curve in curves:
            #print curve
            curve=np.array(curve)
            #print curve
            #target.write("path.move(to:CGPoint(x:"+str(curve[0][0])+",y:"+str(curve[0][1])+"))\n")

            if len(curve)==2: # straight line
            #    target.write("path.addLine(to:CGPoint(x:"+str(curve[1][0])+",y:"+str(curve[1][1])+"))\n")
                verts=[curve[0],curve[1]]
                maxx = max(maxx,curve[0][0],curve[1][0])
                maxy = max(maxy,curve[0][1],curve[1][1])

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
                        maxx = max(maxx,bezier[0][0],bezier[3][0])
                        maxy = max(maxy,bezier[0][1],bezier[3][1])
                        codes = [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4]
                        path = Path(verts, codes)
                        patch = patches.PathPatch(path,facecolor='none', lw=1)
                        #xs, ys = zip(*verts)
                        #ax.plot(xs, ys, 'x--', lw=2, color='black', ms=10)

                        ax.add_patch(patch)


        ax.set_xlim(0,maxx)
        ax.set_ylim(maxy,0)
        print maxx
        print maxy

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
