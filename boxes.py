import sys
sys.path.append('/Users/mac/.virtualenvs/cv/lib/python2.7/site-packages')
import cv2

img = cv2.imread('/Users/mac/imagedev/hindi1.jpeg')
edges = cv2.Canny(img,100,200)
#gray=cv2.cvtColor(edges,cv2.COLOR_BGR2GRAY)
#bounding boxes
image,contours,hierarchy = cv2.findContours(edges,cv2.RETR_LIST,cv2.CHAIN_APPROX_SIMPLE)
idx =0
for cnt in contours:
    idx += 1
    x,y,w,h = cv2.boundingRect(cnt)
    print x,":",y,"-","w=",w,"h=",h
    roi=img[y:y+h,x:x+w]
    cv2.imwrite(str(idx) + '.jpg', roi)
#bounding boxes
    #cv2.rectangle(im,(x,y),(x+w,y+h),(200,0,0),2)
#cv2.imshow('img',img)
#cv2.waitKey(0)
