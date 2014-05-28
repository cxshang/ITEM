
from __future__ import division

cdef extern from *:
	ctypedef char* const_char_ptr "const char*"
	
cdef extern from "stdint.h" :
	ctypedef long int 	int32_t
	ctypedef long long int	int64_t
	ctypedef unsigned long int 	uint32_t
	ctypedef unsigned long long int	uint64_t

cdef extern from "city.h" :
	ctypedef uint32_t uint32
	ctypedef uint64_t uint64
	uint64 CityHash64(const_char_ptr buf, int len) nogil
	uint32 CityHash32(const_char_ptr buf, int len) nogil

cdef  inline uint32  cityhash32(int64_t  keyvalue):
	pystr= str(keyvalue)
	cdef char* keystr=pystr
	cdef int xlen=len(keystr)
	return  CityHash32(keystr,xlen)
	
cdef  inline uint64 cityhash64(int64_t  keyvalue) nogil:
	'''
	cdef char* maps=['0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
	cdef char* keystr=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	cdef int xlen=0
	cdef int i,j,temp
	if keyvalue != 0:
		while keyvalue:
			keystr[xlen]=maps[keyvalue%10]
			xlen+=1
			keyvalue=keyvalue//10
		for i in range(xlen//2):
			temp=keystr[i]
			keystr[i]=keystr[xlen-i-1]
			keystr[xlen-i-1]=temp
	else:
		keystr[0]=maps[0]
		xlen=1
	return  CityHash64(keystr,xlen)
	'''

	cdef char* maps=['0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
	cdef char* keystr=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	cdef int xlen=0
	cdef int i,j,temp
	if keyvalue != 0:
		while keyvalue:
			keystr[xlen]=maps[keyvalue & 15]
			xlen+=1
			keyvalue>>=4
		for i in range(xlen//2):
			temp=keystr[i]
			keystr[i]=keystr[xlen-i-1]
			keystr[xlen-i-1]=temp
	else:
		keystr[0]=maps[0]
		xlen=1

	return  CityHash64(keystr,xlen)


