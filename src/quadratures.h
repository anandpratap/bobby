#ifndef _QUADRATURES_H
#define _QUADRATURES_H

#include "common.h"

class GaussQuadrature1D{
 public:
	const double __eval_points[2] = {-sqrt(1.0/3.0), sqrt(1.0/3.0)};
	const double __weights[2] = {1.0, 1.0};
	const unsigned int __neval_points = 2;
};


class GaussQuadrature2D{
 public:
	double **__eval_points;
	double *__weights;
	const unsigned int __neval_points = 9;

	GaussQuadrature2D(){
		__weights = new double[__neval_points]();
		__eval_points = allocate_2d_array<double>(__neval_points, 2);
		
		double w_0 = 5.0/9.0;
		double w_1 = 8.0/9.0;

		double x_0 = sqrt(3.0/5.0);

		__eval_points[0][0] = -x_0;
		__eval_points[0][1] = -x_0;
		__weights[0] = w_0*w_0;

		__eval_points[1][0] = 0.0;
		__eval_points[1][1] = -x_0;
		__weights[1] = w_1*w_0;

		__eval_points[2][0] = x_0;
		__eval_points[2][1] = -x_0;
		__weights[2] = w_0*w_0;

		
		__eval_points[3][0] = -x_0;
		__eval_points[3][1] = 0.0;
		__weights[3] = w_0*w_1;

		__eval_points[4][0] = 0.0;
		__eval_points[4][1] = 0.0;
		__weights[4] = w_1*w_1;

		__eval_points[5][0] = x_0;
		__eval_points[5][1] = 0.0;
		__weights[5] = w_0*w_1;
		
		__eval_points[6][0] = -x_0;
		__eval_points[6][1] = x_0;
		__weights[6] = w_0*w_0;

		__eval_points[7][0] = 0.0;
		__eval_points[7][1] = x_0;
		__weights[7] = w_1*w_0;

		__eval_points[8][0] = x_0;
		__eval_points[8][1] = x_0;
		__weights[8] = w_0*w_0;
	};

	~GaussQuadrature2D(){
		release_2d_array(__eval_points, __neval_points, 2);
		delete[] __weights;
	};
};


#endif
