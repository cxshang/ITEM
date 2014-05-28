
cdef extern from "stdint.h" :
	ctypedef unsigned long long int	uint64_t


cdef extern from "popcount.hpp" :
	ctypedef uint64_t uint64
	int popcount_1(uint64 x) nogil
	int popcount_2(uint64 x) nogil
	int popcount_3(uint64 x) nogil
	int hammdist(uint64 x, uint64 y) nogil
