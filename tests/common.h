#ifndef _TEST_COMMON_H
#define _TEST_COMMON_H

#include "../src/common.h"

template<class T>
int assert_almost_equal(T a, T b, int prec){
	std::cout<<"Assert not equal: a = "<<a<<" b = "<<b<<std::endl;
	T abs_error = fabs(a - b);
	assert(abs_error < pow(10.0, -prec));
	return 0;
}

#endif
