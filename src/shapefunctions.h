#ifndef _SHAPEFUNCTIONS_H
#define _SHAPEFUNCTIONS_H

#include "common.h"


class ShapeFunction1D{
 private:
	double *__chi_min;
	double *__chi_max;
 public:
	ShapeFunction1D(){
		__chi_min = new double[1]();
		__chi_max = new double[1]();
		__chi_min[0] = -1.0;
		__chi_max[0] = 1.0;
	}

	~ShapeFunction1D(){
		delete[] __chi_min; delete[] __chi_max;
	}

	double value(double *chi, unsigned int node = 0){
		assert(node == 0 || node == 1);
		
		if(node == 0){
			return 0.5*(1.0 - chi[0]);
		}
		else if(node == 1){
			return 0.5*(1.0 + chi[0]);
		}
	};
	double derivative(double *chi, unsigned int node = 0, unsigned int direction = 0){
		assert(node == 0 || node == 1);
		assert(direction == 0);
		
		if(direction == 0){
			if(node == 0){
				return -0.5;
			}
			else if(node == 1){
				return 0.5;
			}
		}
	};
};

class ShapeFunction2D{
 private:
	double *__chi_min;
	double *__chi_max;
 public:
	ShapeFunction2D(){
		__chi_min = new double[2]();
		__chi_max = new double[2]();
		__chi_min[0] = -1.0;
		__chi_min[1] = -1.0;

		__chi_max[0] = 1.0;
		__chi_max[1] = 1.0;

	}

	~ShapeFunction2D(){
		delete[] __chi_min; delete[] __chi_max;
	}

	double value(double *chi, unsigned int node = 0){
		assert(node == 0 || node == 1 || node == 2 || node == 3);
		
		if(node == 0){
			return 0.25*(1.0 - chi[0])*(1.0 - chi[1]);
		}
		else if(node == 1){
			return 0.25*(1.0 + chi[0])*(1.0 - chi[1]);
		}
		else if(node == 2){
			return 0.25*(1.0 + chi[0])*(1.0 + chi[1]);
		}
		else if(node == 3){
			return 0.25*(1.0 - chi[0])*(1.0 + chi[1]);
		}
		
	};
	double derivative(double *chi, unsigned int node = 0, unsigned int direction = 0){
		assert(node == 0 || node == 1 || node == 2 || node == 3);
		assert(direction == 0 || direction == 1);

		if(direction == 0){
			if(node == 0){
				return -0.25*(1.0 - chi[1]);
			}
			else if(node == 1){
				return 0.25*(1.0 - chi[1]);
			}
			else if(node == 2){
				return 0.25*(1.0 + chi[1]);
			}
			else if(node == 3){
				return -0.25*(1.0 + chi[1]);
			}
	
		}
		else if (direction == 1){
			if(node == 0){
				return -0.25*(1.0 - chi[0]);
			}
			else if(node == 1){
				return -0.25*(1.0 + chi[0]);
			}
			else if(node == 2){
				return 0.25*(1.0 + chi[0]);
			}
			else if(node == 3){
				return 0.25*(1.0 - chi[0]);
			}
		}
	};
};


#endif
