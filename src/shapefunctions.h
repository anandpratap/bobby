#ifndef _SHAPEFUNCTIONS_H
#define _SHAPEFUNCTIONS_H

#include "common.h"

class ShapeFunction1D{
 private:
	double __chi_min;
	double __chi_max;
	double (*__function_array[2]) (double chi);
	double (*__dfunction_array[2]) (double chi);
	
 public:
	ShapeFunction1D(){
		
		auto f_0 = [](double chi){return 0.5*(1.0-chi);};
		__function_array[0] = f_0;
		auto f_1 = [](double chi){return 0.5*(1.0+chi);};
		__function_array[1] = f_1;

		auto df_0 = [](double chi){return -0.5;};
		__dfunction_array[0] = df_0;
		auto df_1 = [](double chi){return 0.5;};
		__dfunction_array[1] = df_1;

	}

	double value(double chi, unsigned int mode = 0){
		assert(mode == 0 || mode == 1);
		return (*__function_array[mode])(chi);
	};
	double derivative(double chi, unsigned int mode = 0){
		assert(mode == 0 || mode == 1);
		return (*__dfunction_array[mode])(chi);
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

	double value(double *chi, unsigned int mode = 0){
		assert(mode == 0 || mode == 1 || mode == 2 || mode == 3);
		
		if(mode == 0){
			return 0.25*(1.0 - chi[0])*(1.0 - chi[1]);
		}
		else if(mode == 1){
			return 0.25*(1.0 + chi[0])*(1.0 - chi[1]);
		}
		else if(mode == 2){
			return 0.25*(1.0 + chi[0])*(1.0 + chi[1]);
		}
		else if(mode == 3){
			return 0.25*(1.0 - chi[0])*(1.0 + chi[1]);
		}
		
	};
	double derivative(double *chi, unsigned int mode = 0, unsigned int direction = 0){
		assert(mode == 0 || mode == 1 || mode == 2 || mode == 3);
		assert(direction == 0 || direction == 1);

		if(direction == 0){
			if(mode == 0){
				return -0.25*(1.0 - chi[1]);
			}
			else if(mode == 1){
				return 0.25*(1.0 - chi[1]);
			}
			else if(mode == 2){
				return 0.25*(1.0 + chi[1]);
			}
			else if(mode == 3){
				return -0.25*(1.0 + chi[1]);
			}
	
		}
		else if (direction == 1){
			if(mode == 0){
				return -0.25*(1.0 - chi[0]);
			}
			else if(mode == 1){
				return -0.25*(1.0 + chi[0]);
			}
			else if(mode == 2){
				return 0.25*(1.0 + chi[0]);
			}
			else if(mode == 3){
				return 0.25*(1.0 - chi[0]);
			}
		}
	};
};


#endif
