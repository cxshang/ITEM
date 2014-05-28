
cdef extern from "arrcontainer.h":
        cdef cppclass arrcontainer[T]:
                arrcontainer()
                arrcontainer(int)
                void init(int)
                T operator[](int) nogil
                void push_back(T*) nogil
                int size()
