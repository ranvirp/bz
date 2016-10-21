#find loops in an image using opencv
# let us try probablistic Hough transform
# ro and theta for a loop image
#use on a binary image
# cv.FindContours(image, storage, mode=CV_RETR_LIST, method=CV_CHAIN_APPROX_SIMPLE, offset=(0, 0))  contours
#cv.DrawContours(img, contour, external_color, hole_color, max_level, thickness=1, lineType=8, offset=(0, 0))  None
import sys
sys.path.append('/Users/mac/.virtualenvs/cv/lib/python2.7/site-packages')
print "importing cv2"
import cv2
print "import finished"
img = cv2.imread('/Users/mac/imagedev/hindi1.jpeg')
edges = cv2.Canny(img,100,200)
(im1,contours,hierarchy) = cv2.findContours(edges,cv2.RETR_LIST,cv2.CHAIN_APPROX_SIMPLE)
cv2.drawContours(img,contours,-1,(0,255,0),3)
cv2.imwrite('/Users/mac/imagedev/hindi1-contours.jpeg',img)
cv2.namedWindow( "Contours", CV_WINDOW_AUTOSIZE );
cv2.imshow('Contours',img)
cv2.waitKey(0)# this does not close after pressing 0
