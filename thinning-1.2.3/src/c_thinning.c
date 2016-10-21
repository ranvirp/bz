#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include "arrayobject.h"
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <limits.h>

static PyObject *guo_hall_thinning(PyObject *self, PyObject *args);
int _guo_hall_thinning(unsigned char* binary_image, int width, int height);
void initthinning(void);

/* ==== Set up the methods table ====================== */
static PyMethodDef thinningMethods[] = {
	{"guo_hall_thinning",guo_hall_thinning, METH_VARARGS,
	"Takes a 2D numpy UBYTE array in C-order and thins it in place using the algorithm by Guo and Hall."
	"Images that come out of cv2.cvtColor(img, cv2.COLOR_BGR2GRAY) have the right format."
	"\n\n"
	"We assume that the dimensions of the image fit into an int on your platform. If your computer for some"
	"reason has a 2 byte int and lots of memory so that the image can become too large, bad things can happen."
	"\n\n"
	"interface:\n"
	"\tguo_hall_thinning(segmented_image)"
	"\tsegmented_image is a NumPy matrix,"
	"\treturns the same NumPy matrix (thinned)"},
	{NULL, NULL, 0, NULL}     /* Sentinel - marks the end of this structure */
};

/* ==== Initialize the C_test functions ====================== */
void initthinning()  {
	PyObject* module = Py_InitModule3("thinning",thinningMethods, "Thinning of segmented images. See https://bitbucket.org/adrian_n/thinning.");
	PyModule_AddStringConstant(module, "__author__", "Adrian Neumann <adrian_neumann@gmx.de>");
	PyModule_AddStringConstant(module, "__version__", "1.2.3");
	import_array();  // Must be present for NumPy.  Called first after above line.
}

/* ==== Guo Hall Thinning =========
	Takes a 2D numpy UBYTE array in C-order and thins it in place using the algorithm by Guo and Hall.
	Images that come out of cv2.cvtColor(img, cv2.COLOR_BGR2GRAY) have the right format.

	We assume that the dimensions of the image fit into an int on your platform. If your computer for some
	reason has a 2 byte int and lots of memory so that the image can become too large, bad things can happen.

	interface:  guo_hall_thinning(segmented_image)
				segmented_image is a NumPy matrix,
				returns the same NumPy matrix (thinned)
*/
static PyObject *guo_hall_thinning(PyObject *self, PyObject *args)
{
	PyArrayObject *segmented_image;

	/* Parse tuples separately since args will differ between C fcns */
	if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &segmented_image)) {
		return NULL;
	}
	if (NULL == segmented_image) {
		PyErr_SetString(PyExc_TypeError, "Parameter is not a valid image");
		return NULL;
	}
	if (PyArray_TYPE(segmented_image) != NPY_UBYTE || !PyArray_CHKFLAGS(segmented_image, NPY_ARRAY_CARRAY)) {
		PyErr_SetString(PyExc_TypeError, "Parameter is not a grayscale image");
		return NULL;
	}

	npy_intp* shape = PyArray_DIMS(segmented_image);

	int height = (int)shape[0];
	int width = (int)shape[1];

	unsigned char *in_data = PyArray_DATA(segmented_image);

	if (height>=3 && width>=3) {
		int ok = _guo_hall_thinning(in_data, width, height);
		if (ok<0) {
			return PyErr_NoMemory();
		}
	}
	Py_INCREF(segmented_image);
	return (PyObject*)segmented_image;
}

int nonzero_clever(const unsigned char* arr, unsigned int start, unsigned int len) {
	/* find the first nonzero element from arr[start] to arr[start+len-1] (inclusive)
	   look at a long long at a time to be faster on 64 bit cpus */
	const unsigned int step=sizeof(unsigned long long)/sizeof(unsigned char);
	unsigned int i=start;
	//unsigned types should throw exceptions on under/overflow...
	while(len>step && i<len-step) {
			if (*((unsigned long long*)(arr +i))==0) {
				i+=step;
			} else {
				int j=0;
				while(arr[i+j]==0) j++;
				return i+j;
			}
	}
	while(i<len) {
		if (arr[i]!=0) { return i;}
		i++;
	}
	return len;
}


int guo_hall_iteration(const unsigned char* binary_image, unsigned char* mask, const unsigned int width, const unsigned int height, const int iteration) {
		/* one iteration of the algorithm by guo and hall. see their paper for an explanation.
		   We only consider nonzero elemets of the image. We never reinitialize the mask, once a pixel is
		   black, it will never become white again anyway. */
		unsigned int changed = 0;
		for (unsigned int j = 1; j < height-1; j++) {
			const unsigned char* line = binary_image+j*width;
			unsigned int start=0;
			const int len = width-1;

			while(start+1<len) {
				start = nonzero_clever(line, start+1, len);
				if (start==len) break;

				const unsigned int i = start;
				assert(line[i]!=0);
				assert(binary_image[i + j*width]!=0);

				const bool p2 = binary_image[i-1 + width*j];
				const bool p6 = binary_image[i+1 + width*j];

				const bool p9 = binary_image[i-1 + width*(j-1)];
				const bool p8 = binary_image[i   + width*(j-1)];
				const bool p7 = binary_image[i+1 + width*(j-1)];

				const bool p3 = binary_image[i-1 + width*(j+1)];
				const bool p4 = binary_image[i   + width*(j+1)];
				const bool p5 = binary_image[i+1 + width*(j+1)];
				const unsigned int C = ((!p2 && (p3 || p4)) +
					(!p4 && (p5 || p6)) +
					(!p6 && (p7 || p8)) +
					(!p8 && (p9 || p2)));
				// printf("%d %d %d %d %d %d %d %d\n",p2,p3,p4,p5,p6,p7,p8,p9);
				if (C==1) {
					const unsigned int N1 = (p9 || p2) + (p3 || p4) + (p5 || p6) + (p7 || p8);
					const unsigned int N2 = (p2 || p3) + (p4 || p5) + (p6 || p7) + (p8 || p9);
					const unsigned int N = N1 < N2 ? N1 : N2;
					unsigned int m;

					if (iteration == 0)
						{m = (p8 && (p6 || p7 || !p9));}
					else
						{m = (p4 && (p2 || p3 || !p5));}

					if (2 <= N && N <= 3 && m == 0)   {
						mask[i + width*j] = 0;
						changed += 1;
					}
				}
			}

		}
		return changed;
}

void andImage(unsigned char* image, const unsigned char* mask, const int size) {
	/* calculate image &=mask.
	   to be faster on 64 bit cpus, we do this one long long at a time */
	const int step = sizeof(unsigned long long)/sizeof(unsigned char);
	unsigned long long* image_l = (unsigned long long*)image;
	const unsigned long long* mask_l = (unsigned long long*) mask;
	unsigned int i=0;
	for(; size/step>2 && i<size/step-2; i+=2) {
		image_l[i] = image_l[i] & mask_l[i];
		image_l[i+1] = image_l[i+1] & mask_l[i+1];
	}
	for(i=i*step; i<size; ++i) {
		image[i] = image[i] & mask[i];
	}
}

int _guo_hall_thinning(unsigned char* binary_image, int width, int height) {
	/* return -1 if we can't allocate the memory for the mask, else 0 */
	int changed;
	unsigned char* mask = (unsigned char*) malloc(width*height*sizeof(unsigned char));
	if (mask==NULL) {
		return -1;
	}

	memset(mask, UCHAR_MAX, width*height);
	do {
		changed = guo_hall_iteration(binary_image, mask, width, height, 0);
		andImage(binary_image, mask, width*height);

		changed += guo_hall_iteration(binary_image, mask, width, height, 1);
		andImage(binary_image, mask, width*height);
	} while (changed != 0);
	free(mask);

	return 0;
}

