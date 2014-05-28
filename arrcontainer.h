
#ifndef ARRCONTAINER_H
#define ARRCONTAINER_H
 
#include <vector>
 
template <typename T>
class arrcontainer
{
private:
	std::vector<T*> container;
public:
	arrcontainer() {};
	
	arrcontainer(int n):container() { container.reserve(n);
		for(int j=0; j<n; ++j) 	{ container.push_back(new T); } 
					}
	
	~arrcontainer(){
		for ( unsigned int j=0; j<container.size(); ++j)
                                        { delete container[j]; }  } 

	inline void init(int n){ for(int j=0; j<n; ++j)  { container.push_back(new T);} }
    
	inline void push_back( T *const val){ container.push_back(val); }
	
	inline T&  operator[](const int j)
                        { return (*(container[j]));  }

	inline int size(){ return int(container.size()); }

};
 
#endif

