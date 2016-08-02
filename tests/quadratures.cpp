#include "../src/common.h"
#include "../src/shapefunctions.h"
#include "../src/quadratures.h"

#include "common.h"

void test_constant_functions(void){
	ShapeFunction1D shapefunctions = ShapeFunction1D();
	GaussQuadrature1D quadrature = GaussQuadrature1D();
	
	auto func_1 = [](double chi){return 0.0;};
	assert_almost_equal(quadrature.integral(func_1), 0.0, 13);
	
	auto func_2 = [](double chi){return 1.0;};
	assert_almost_equal(quadrature.integral(func_2), 2.0, 13);
	
	auto func_3 = [](double chi){return -1.0;};
	assert_almost_equal(quadrature.integral(func_3), -2.0, 13);
}

void test_linear_functions(void){
	ShapeFunction1D shapefunctions = ShapeFunction1D();
	GaussQuadrature1D quadrature = GaussQuadrature1D();
	
	auto func_1 = [](double chi){return chi;};
	assert_almost_equal(quadrature.integral(func_1), 0.0, 13);
	
	auto func_2 = [](double chi){return chi+1.0;};
	assert_almost_equal(quadrature.integral(func_2), 2.0, 13);
	
	auto func_3 = [](double chi){return chi-1.0;};
	assert_almost_equal(quadrature.integral(func_3), -2.0, 13);
}


void test_multiple_functions(void){
	ShapeFunction1D shapefunctions = ShapeFunction1D();
	GaussQuadrature1D quadrature = GaussQuadrature1D();
	
	auto func_1 = [](double chi){return 1.0;};
	auto func_2 = [](double chi){return chi;};
	assert_almost_equal(quadrature.integral(func_1, func_2), 0.0, 13);

	assert_almost_equal(quadrature.integral(func_2, func_2), 2.0/3.0, 13);
	
	auto func_3 = [](double chi){return chi+1.0;};
	assert_almost_equal(quadrature.integral(func_2, func_3), 2.0/3.0, 13);
	assert_almost_equal(quadrature.integral(func_3, func_2), 2.0/3.0, 13);

	
	auto func_4 = [](double chi){return chi-1;};
	assert_almost_equal(quadrature.integral(func_2, func_3, func_4), 0.0, 13);

	auto func_5 = [](double chi){return chi+1;};
	assert_almost_equal(quadrature.integral(func_2, func_3, func_5), 4.0/3.0, 13);
}


int main(void){
	test_constant_functions();
	test_linear_functions();
	test_multiple_functions();
	return 0;
}
