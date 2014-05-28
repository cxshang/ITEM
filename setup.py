from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


ext_modules=[
    Extension("docitem",
        ["docitem.pyx"],libraries=["cityhash","stdc++"],
        library_dirs=[r"/bwfs/home/cxshang/cityhash/lib",r"/usr/lib64"],
        include_dirs=[r"/bwfs/home/cxshang/sparsehash/include",r"/bwfs/home/cxshang/python-2.7.5/lib/python2.7/site-packages/numpy/core/include",r"/bwfs/home/cxshang/cityhash/include",r"/usr/include"],
        language="c++",
	extra_compile_args=['-O3', '-funsafe-math-optimizations', '-ffast-math',  '-fopenmp', '-mtune=generic', '-pipe','-msse',  '-msse2',  '-msse3',  '-mfpmath=sse'],
	extra_link_args=['-fopenmp'],
    )
]

setup(
    name = "docitem",
    cmdclass = {"build_ext": build_ext},
    ext_modules = ext_modules
)

