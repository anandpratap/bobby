#ifndef _QUADRATURES_H
#define _QUADRATURES_H

#include "common.h"

class GaussQuadrature1D{
 private:
	const double __eval_points[2] = {-sqrt(1.0/3.0), sqrt(1.0/3.0)};
	const double __weights[2] = {1.0, 1.0};
	const unsigned int __neval_points = 2;
 public:
	GaussQuadrature1D(){};
	double integral(double (*func)(double)){
		double integral_value = 0.0;
		for(int n=0; n<__neval_points; n++){
			integral_value += func(__eval_points[n])*__weights[n];
		}
		return integral_value;
	};

	double integral(double (*func_1)(double), double (*func_2)(double)){
		double integral_value = 0.0;
		for(int n=0; n<__neval_points; n++){
			double prod = func_1(__eval_points[n])*func_2(__eval_points[n]);
			integral_value += prod*__weights[n];
		}
		return integral_value;
	};


	double integral(double (*func_1)(double), double (*func_2)(double), double (*func_3)(double)){
		double integral_value = 0.0;
		for(int n=0; n<__neval_points; n++){
			double prod = func_1(__eval_points[n])*func_2(__eval_points[n])*func_3(__eval_points[n]);
			integral_value += prod*__weights[n];
		}
		return integral_value;
	};

	double integral(double (*func_1)(double), double (*func_2)(double), double (*func_3)(double), double (*func_4)(double)){
		double integral_value = 0.0;
		for(int n=0; n<__neval_points; n++){
			double prod = func_1(__eval_points[n])*func_2(__eval_points[n])*func_3(__eval_points[n])*func_4(__eval_points[n]);
			integral_value += prod*__weights[n];
		}
		return integral_value;
	};

};


#endif
