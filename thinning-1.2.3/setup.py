import os.path
import setuptools
from setuptools import setup, Extension

#this is a collosal hack and I hate myself for having to use it.
setuptools.dist.Distribution(dict(setup_requires='numpy>=1.7'))
__builtins__.__NUMPY_SETUP__ = False
import numpy

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the relevant file
with open(os.path.join(here, 'README.txt')) as f:
	long_description = f.read()


setup(
	name='thinning',
	version='1.2.3',
	author = "Adrian Neumann",
	author_email = "adrian_neumann@gmx.de",
	url = "https://bitbucket.org/adrian_n/thinning",
	description = "C extension for thinning binary images.",
	ext_modules=[
		Extension(
			'thinning',
			['src/c_thinning.c'],
			include_dirs=[numpy.get_include(), os.path.join(numpy.get_include(),"numpy")],
			extra_compile_args = ["-std=c99"]
		)
	],
	classifiers=[
		"Programming Language :: Python",
		"Programming Language :: Python :: 3", #I think
		"Programming Language :: Python :: 2",
		"Programming Language :: C",
		"License :: OSI Approved :: BSD License",
		"Topic :: Scientific/Engineering :: Image Recognition",
		"Intended Audience :: Developers",
		"Development Status :: 4 - Beta"
	],
	keywords=["image processing", "thinning", "guo hall", "skeletonization"],
	long_description = long_description,
	install_requires=['numpy>=1.7'],
	license='BSD'
)