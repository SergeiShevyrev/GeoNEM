#how to do
# https://cython.readthedocs.io/en/latest/src/tutorial/cython_tutorial.html#the-basics-of-cython
#
#inclusion of libraries
#https://stackoverflow.com/questions/14657375/cython-fatal-error-numpy-arrayobject-h-no-such-file-or-directory

from setuptools import setup
from Cython.Build import cythonize
import numpy

#1 main_b
setup(
    ext_modules = cythonize("main_b.pyx"),
    include_dirs=[numpy.get_include()]
)
#2 extension
setup(
   ext_modules = cythonize("uilib.pyx"),
   include_dirs=[numpy.get_include()]
)


#command line of your system
#firstly build
#python setup.py build
#
#next install
#
#python setup.py install

#build in local directory 
#python setup.py build_ext --inplace

#launch, firstly launch python , than 
#>>> import Intrusion_emplacement