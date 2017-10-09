from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension("adj_list", ["adj_list.pyx"])
]

setup(
    ext_modules = cythonize(extensions),
    include_dirs=[np.get_include()]
)