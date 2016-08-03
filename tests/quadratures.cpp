#include "../src/common.h"
#include "../src/shapefunctions.h"
#include "../src/quadratures.h"

#include "common.h"

void test_functions_1d(){
	double integral = 0.0;
	GaussQuadrature1D quadrature = GaussQuadrature1D();
	
	auto f0 = [](double *chi){return 1;}; 
	integral = 0.0;
	for(int i=0; i<quadrature.__neval_points; i++){
		integral += f0(quadrature.__eval_points[i])*quadrature.__weights[i];
	}
	assert_almost_equal(integral, 2.0, 1e-10);


	auto f1 = [](double *chi){return chi[0];}; 
	integral = 0.0;
	for(int i=0; i<quadrature.__neval_points; i++){
		integral += f1(quadrature.__eval_points[i])*quadrature.__weights[i];
	}
	assert_almost_equal(integral, 0.0, 1e-10);


	auto f2 = [](double *chi){return chi[0]*chi[0];}; 
	integral = 0.0;
	for(int i=0; i<quadrature.__neval_points; i++){
		integral += f2(quadrature.__eval_points[i])*quadrature.__weights[i];
	}
	assert_almost_equal(integral, 2.0/3.0, 1e-10);

	auto f3 = [](double *chi){return chi[0]*chi[0]*chi[0];}; 
	integral = 0.0;
	for(int i=0; i<quadrature.__neval_points; i++){
		integral += f3(quadrature.__eval_points[i])*quadrature.__weights[i];
	}
	assert_almost_equal(integral, 0.0, 1e-10);

}

void test_functions_2d(){
	double integral = 0.0;
	GaussQuadrature2D quadrature = GaussQuadrature2D();
	
	auto f0 = [](double *chi){return 1;}; 
	integral = 0.0;
	for(int i=0; i<quadrature.__neval_points; i++){
		integral += f0(quadrature.__eval_points[i])*quadrature.__weights[i];
	}
	assert_almost_equal(integral, 4.0, 1e-10);


	auto f1 = [](double *chi){return chi[0];}; 
	integral = 0.0;
	for(int i=0; i<quadrature.__neval_points; i++){
		integral += f1(quadrature.__eval_points[i])*quadrature.__weights[i];
	}
	assert_almost_equal(integral, 0.0, 1e-10);


	auto f2 = [](double *chi){return chi[1];}; 
	integral = 0.0;
	for(int i=0; i<quadrature.__neval_points; i++){
		integral += f2(quadrature.__eval_points[i])*quadrature.__weights[i];
	}
	assert_almost_equal(integral, 0.0, 1e-10);

	auto f3 = [](double *chi){return chi[0]+chi[1];}; 
	integral = 0.0;
	for(int i=0; i<quadrature.__neval_points; i++){
		integral += f3(quadrature.__eval_points[i])*quadrature.__weights[i];
	}
	assert_almost_equal(integral, 0.0, 1e-10);

	auto f4 = [](double *chi){return chi[0]*chi[1];}; 
	integral = 0.0;
	for(int i=0; i<quadrature.__neval_points; i++){
		integral += f4(quadrature.__eval_points[i])*quadrature.__weights[i];
	}
	assert_almost_equal(integral, 0.0, 1e-10);

}




int main(void){
	test_functions_1d();
	test_functions_2d();
	return 0;
}
