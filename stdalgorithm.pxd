
from libc.stdint cimport *
DEF USEING_SPARSE_MAPANDSET = False  #indicate whether using google spare map and set or cpp map and set

from libcpp.vector cimport  vector as cppvector
ctypedef cppvector[int32_t].iterator iteratorofcppvectorint32
ctypedef cppvector[int64_t].iterator iteratorofcppvectorint64
cdef extern from "algorithm" namespace "std":
	cdef void sortcppvectorint32 "std::sort"(iteratorofcppvectorint32,iteratorofcppvectorint32)
	cdef void sortcppvectorint64 "std::sort"(iteratorofcppvectorint64,iteratorofcppvectorint64)

from libcpp.set cimport  set as cppset
ctypedef cppset[int32_t].iterator iteratorofcppsetint32
ctypedef cppset[int64_t].iterator iteratorofcppsetint64
cdef extern from "algorithm" namespace "std":
	cdef void sortcppsetint32 "std::sort"(iteratorofcppsetint32,iteratorofcppsetint32)
	cdef void sortcppsetint64 "std::sort"(iteratorofcppsetint64,iteratorofcppsetint64)

IF USEING_SPARSE_MAPANDSET:
	from sset cimport sparse_hash_set as sset
	ctypedef sset[int32_t].iterator iteratorofssetint32
	ctypedef sset[int64_t].iterator iteratorofssetint64
	cdef extern from "algorithm" namespace "std":
		cdef void sortssetint32 "std::sort"(iteratorofssetint32,iteratorofssetint32)
		cdef void sortssetint64 "std::sort"(iteratorofssetint64,iteratorofssetint64)

cdef extern from "algorithm" namespace "std":
	#for std set
	cdef iteratorofcppsetint32 cppset_union_32 "std::set_union"\
		(iteratorofcppsetint32,iteratorofcppsetint32,iteratorofcppsetint32,iteratorofcppsetint32,iteratorofcppsetint32)
	cdef iteratorofcppsetint64 cppset_union_64 "std::set_union"\
		(iteratorofcppsetint64,iteratorofcppsetint64,iteratorofcppsetint64,iteratorofcppsetint64,iteratorofcppsetint64)
	
	cdef iteratorofcppsetint32 cppset_intersection_32 "std::set_intersection"\
		(iteratorofcppsetint32,iteratorofcppsetint32,iteratorofcppsetint32,iteratorofcppsetint32,iteratorofcppsetint32)
	cdef iteratorofcppsetint64 cppset_intersection_64 "std::set_intersection"\
		(iteratorofcppsetint64,iteratorofcppsetint64,iteratorofcppsetint64,iteratorofcppsetint64,iteratorofcppsetint64)
	
	cdef iteratorofcppsetint32 cppset_difference_32 "std::set_difference"\
		(iteratorofcppsetint32,iteratorofcppsetint32,iteratorofcppsetint32,iteratorofcppsetint32,iteratorofcppsetint32)
	cdef iteratorofcppsetint64 cppset_difference_64 "std::set_difference"\
		(iteratorofcppsetint64,iteratorofcppsetint64,iteratorofcppsetint64,iteratorofcppsetint64,iteratorofcppsetint64)
	

	#for sparse_hash_set 
	IF USEING_SPARSE_MAPANDSET:
		cdef iteratorofssetint32 sset_union_32 "std::set_union"\
			(iteratorofssetint32,iteratorofssetint32,iteratorofssetint32,iteratorofssetint32,iteratorofssetint32)
		cdef iteratorofssetint64 sset_union_64 "std::set_union"\
			(iteratorofssetint64,iteratorofssetint64,iteratorofssetint64,iteratorofssetint64,iteratorofssetint64)
		
		cdef iteratorofssetint32 sset_intersection_32 "std::set_intersection"\
			(iteratorofssetint32,iteratorofssetint32,iteratorofssetint32,iteratorofssetint32,iteratorofssetint32)
		cdef iteratorofssetint64 sset_intersection_64 "std::set_intersection"\
			(iteratorofssetint64,iteratorofssetint64,iteratorofssetint64,iteratorofssetint64,iteratorofssetint64)
			
		cdef iteratorofssetint32 sset_difference_32 "std::set_difference"\
			(iteratorofssetint32,iteratorofssetint32,iteratorofssetint32,iteratorofssetint32,iteratorofssetint32)
		cdef iteratorofssetint64 sset_difference_64 "std::set_difference"\
			(iteratorofssetint64,iteratorofssetint64,iteratorofssetint64,iteratorofssetint64,iteratorofssetint64)
			
