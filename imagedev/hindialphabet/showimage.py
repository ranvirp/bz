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
filename = sys.argv[1]
img = cv2.imread(filename)
if len(sys.argv) ==2:
    plt.imshow(img,sys.argv[2])
else:
      plt.imshow(img)
plt.show()
