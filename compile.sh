#!/bin/bash

python setup.py build_ext --inplace
gcc -o item item.c ./docitem.so  -I/bwfs/home/cxshang/python-2.7.5/include/python2.7 -L /bwfs/home/cxshang/python-2.7.5/lib -lpython2.7

