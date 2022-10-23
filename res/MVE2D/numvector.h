#ifndef NUMVECTOR_H_
#define NUMVECTOR_H_

#include <vector>

using namespace std;

template<typename T, int n>
class numvector
{
private:
	T r[n];
public:
	T& operator[](int j) { return r[j]; }
	const T& operator[](int j) const { return r[j]; }

	int size() const { return n; } 

	numvector<T, n>& operator=(const numvector<T, n>& vec) {
		for (int i = 0; i < n; ++i)
			r[i] = vec.r[i];
		return *this;
	}

	numvector()	{};

	numvector(const numvector<T, n>& vec){
		for (int i = 0; i < n; ++i)
			r[i] = vec.r[i];
	}

	numvector(const initializer_list<T>& z){
		//if (z.size() > n) exit(0);
		//for (size_t i = 0; i < z.size(); ++i)
		//	r[i] = *(z.begin() + i); 
		//for (size_t i = z.size(); i < n; ++i)
		//	r[i] = 0;

		for (size_t i = 0; i < n; ++i)
			r[i] = *(z.begin() + i);
	}
};



#endif