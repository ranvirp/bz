/*
1. 10 files of hindi text images- by tearing pages from newspapers and scanning using iphone
2. Perform various operations - convolutions, hough transform etc
3. also evaluate scantailor https://github.com/scantailor/scantailor
4. also evaluate https://github.com/tmbdev/ocropy
3. use tensorflow to train this data ---predict using hough transform by applying a model

Learning
-- installed opencv using http://www.pyimagesearch.com/2015/06/15/install-opencv-3-0-and-python-2-7-on-osx/
--cmake -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_INSTALL_PREFIX=/usr/local \
	-D PYTHON2_PACKAGES_PATH=~/.virtualenvs/cv/lib/python2.7/site-packages \
	-D PYTHON2_LIBRARY=/usr/local/Cellar/python/2.7.12_2/Frameworks/Python.framework/Versions/2.7/bin \
	-D PYTHON2_INCLUDE_DIR=/usr/local/Frameworks/Python.framework/Headers \
	-D INSTALL_C_EXAMPLES=ON -D INSTALL_PYTHON_EXAMPLES=ON \
	-D BUILD_EXAMPLES=ON \
	-D OPENCV_EXTRA_MODULES_PATH=~/opencv_contrib/modules ..


*/
