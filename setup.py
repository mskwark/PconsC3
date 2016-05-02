from __future__ import division
from setuptools import setup, Extension

from Cython.Build import cythonize

#setup(name='predict_tree', ext_modules=cythonize('predict_tree.pyx'))
#setup(name='predict_tree2', ext_modules=cythonize('predict_tree2.pyx'))


extension = Extension('_predict_parallel', ['_predict_parallel.pyx'],
                      extra_compile_args='-O2 -march=native -pipe -g0 -fopenmp -std=c11'.split(),
                      extra_link_args='-O2 -march=native -pipe -g0 -fopenmp -std=c11'.split())
setup(name='_predict_parallel', ext_modules=cythonize(extension))

