ITEM
====

an overlapping community detecting algorithm using information theory and expectation and maximization process

1. denpending softwares <br>
In my experiments, ITEM running Red Hat Enterprise Linux Server release 5.3 (Tikanga), <br>
python2.7.5 <br>
numpy1.7.1  <br>
cython0.19.1 <br>
cityhash, which codes can be found at https://code.google.com/p/cityhash <br>
igraph0.6.5 and python-igraph 0.6.5 <br>
openmp <br>
C++ STL <br>
<br>
ITEM also uses the log function in fmath.hpp, which can be found at https://github.com/herumi/fmath <br>
You should download fmath.hpp and save it in this directory. <br>
<br>
Above softwares (may be higher versions) must be install before running ITEM. <br>
<br>
if you want use the google's sparsehash to replace corresponding ones of C++ STL,  <br>
then sparesehash must be installed, which can be found at https://code.google.com/p/sparsehash/ , <br>
and must modify the code at 50th line of docitem.pyx from <br>
DEF USEING_SPARSE_MAPANDSET = False <br>
to <br>
DEF USEING_SPARSE_MAPANDSET = True <br>
<br>
meanwhile, please modify the code at 3th line of stdalgorithm.pxd from <br>
DEF USEING_SPARSE_MAPANDSET = False <br>
to <br>
DEF USEING_SPARSE_MAPANDSET = True <br>

<br>
2. compile and link  <br>
To compile and link, please first modify setup.py an compile.sh to make the linking and including options <br>
correct according your software configures.  Then run ./compile on console. <br>

#++++++++++++++++++++++++++++++++++++
3. run item
For the usages of item, just type ./item on console to see help.
