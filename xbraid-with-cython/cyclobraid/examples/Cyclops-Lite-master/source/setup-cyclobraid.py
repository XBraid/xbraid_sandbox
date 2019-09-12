from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os
# import numpy

os.environ["CC"] = "mpicc"
os.environ["LDSHARED"] = "mpicc -shared"

cyclobraid_extension = Extension(
    name="cyclobraid",
    sources=["cyclobraid.pyx"],
    libraries=["braid"],
    library_dirs=["../../../braid"],
    include_dirs=["../../../braid"], ##, numpy.get_include()],
    extra_compile_args=["-Wno-incompatible-pointer-types"] ##, "-fPIC"]   
)
setup(
    name="cyclobraid",
    ext_modules=cythonize([cyclobraid_extension])
)


