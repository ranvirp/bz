Thinning
========

Thinning is the operation that takes a binary image and contracts the foreground until only single-pixel wide lines remain. It is also known as skeletonization.

This package implements the thinning algorithm by Guo and Hall[1] for Numpy arrays. It is thus compatible with OpenCV. The algorithm is implemented in C and fairly fast.

[1] http://dx.doi.org/10.1145/62065.62074