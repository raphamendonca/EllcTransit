from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy


ext_modules = [
    Extension(
        "ellc.lcOpenMp",
        ["lcOpenMp.pyx"],
        #libraries=["gomp"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
        include_dirs=[numpy.get_include()],
        language='python'
    )
]

setup(
    name = 'ellc.lcOpenMp',
    ext_modules = cythonize(ext_modules)
)
