#ifndef _EQUATIONS_H
#define _EQUATIONS_H

#include "common.h"

class Equation{
 private:
	void allocate(unsigned int nvar){
		flux.resize(nvar);
		wavespeed.resize(nvar);
		
		dflux.resize(nvar);
		for (int i = 0; i < nvar; ++i){
			dflux[i].resize(nvar);
		}
	};
 public:
	unsigned int nvar;
	std::vector<func_pointer> flux;
	std::vector<std::vector<func_pointer>> dflux;
	std::vector<func_pointer> wavespeed;
	
	Equation(unsigned int val_nvar = 1){
		nvar = val_nvar;
		allocate(nvar);
	};

	void setup(void){
		set_relations();
	}

	virtual void set_relations(void){
		auto f_1 = [](double *u){return u[0];};
		flux[0] = f_1;
		
		auto df_1 = [](double *u){return 1.0;};
		dflux[0][0] = df_1;

		auto wavespeed_1 = [](double *u){return 1.0;};
		wavespeed[0] = wavespeed_1;
	}
};



class BurgersEquation: public Equation{
 public:
	virtual void set_relations(void){
		auto f_1 = [](double *u){return u[0]*u[0]/2.0;};
		flux[0] = f_1;
		
		auto df_1 = [](double *u){return u[0];};
		dflux[0][0] = df_1;

		auto wavespeed_1 = [](double *u){return u[0];};
		wavespeed[0] = wavespeed_1;
	}
};



class SystemEquation: public Equation{
 public:
 SystemEquation()
	 : Equation(2){
	}

	virtual void set_relations(void){
		auto f_0 = [](double *u){return u[0];};
		flux[0] = f_0;

		auto f_1 = [](double *u){return u[1];};
		flux[1] = f_1;

		auto df_00 = [](double *u){return 1.0;};
		dflux[0][0] = df_00;

		auto df_01 = [](double *u){return 0.0;};
		dflux[0][1] = df_01;

		auto df_10 = [](double *u){return 0.0;};
		dflux[1][0] = df_10;
		
		auto df_11 = [](double *u){return 1.0;};
		dflux[1][1] = df_11;

		auto wavespeed_0 = [](double *u){return 1.0;};
		wavespeed[0] = wavespeed_0;

		auto wavespeed_1 = [](double *u){return 1.0;};
		wavespeed[1] = wavespeed_1;

	}
};


#endif
