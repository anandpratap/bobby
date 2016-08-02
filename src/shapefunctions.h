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


#endif
