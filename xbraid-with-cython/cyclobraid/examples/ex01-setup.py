from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os
# import numpy

os.environ["CC"] = "mpicc"
os.environ["LDSHARED"] = "mpicc -shared"

ex01_extension = Extension(
    name="ex01",
    sources=["ex01.pyx"],
    libraries=["braid"],
    library_dirs=["../braid"],
    include_dirs=["../braid"], ##, numpy.get_include()],
    extra_compile_args=["-Wno-incompatible-pointer-types"] ##, "-fPIC"]   
)
setup(
    name="ex01",
    ext_modules=cythonize([ex01_extension])
)


