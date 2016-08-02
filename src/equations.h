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

	virtual double get_tau(double *u, double dt, int ivar, double dchidx){
		double tau_dt_term = 4.0/dt/dt;
		double tau_advection_term = 4.0*pow(wavespeed[ivar](u), 2.0)*(dchidx*dchidx);
		double tau_chi =  1.0/sqrt(tau_dt_term + tau_advection_term);
		return tau_chi;
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

		auto f_1 = [](double *u){return 0.5*u[1];};
		flux[1] = f_1;

		auto df_00 = [](double *u){return 1.0;};
		dflux[0][0] = df_00;

		auto df_01 = [](double *u){return 0.0;};
		dflux[0][1] = df_01;

		auto df_10 = [](double *u){return 0.0;};
		dflux[1][0] = df_10;
		
		auto df_11 = [](double *u){return 0.5;};
		dflux[1][1] = df_11;

		auto wavespeed_0 = [](double *u){return 1.0;};
		wavespeed[0] = wavespeed_0;

		auto wavespeed_1 = [](double *u){return 0.5;};
		wavespeed[1] = wavespeed_1;

	}
};



class EulerEquation: public Equation{
 public:
 EulerEquation()
	 : Equation(3){
	}

	virtual void set_relations(void){
		auto f_0 = [](double *u){return u[1];};
		flux[0] = f_0;

		auto f_1 = [](double *u){return pow(u[1], 2)/u[0] + (GAMMA-1)*(u[2] - 0.5*pow(u[1], 2)/u[0]);};
		flux[1] = f_1;

		auto f_2 = [](double *u){return u[1]/u[0]*(u[2] + (GAMMA-1)*(u[2] - 0.5*pow(u[1], 2)/u[0]));};
		flux[2] = f_2;

		auto df_00 = [](double *u){return 0.0;};
		dflux[0][0] = df_00;

		auto df_01 = [](double *u){return 1.0;};
		dflux[0][1] = df_01;

		auto df_02 = [](double *u){return 0.0;};
		dflux[0][2] = df_02;

		
		auto df_10 = [](double *u){return -0.5*(3-GAMMA)*(pow(u[1],2)/pow(u[0],2));};
		dflux[1][0] = df_10;
		
		auto df_11 = [](double *u){return (3 - GAMMA)*u[1]/u[0];};
		dflux[1][1] = df_11;

		auto df_12 = [](double *u){return (GAMMA-1.0);};
		dflux[1][2] = df_12;


		auto df_20 = [](double *u){return -GAMMA*u[1]*u[2]/pow(u[0],2) + (GAMMA-1)*pow(u[1]/u[0], 3);};
		dflux[2][0] = df_20;
		
		auto df_21 = [](double *u){return GAMMA*u[2]/u[0] - 1.5*(GAMMA-1)*pow(u[1]/u[0], 2);};
		dflux[2][1] = df_21;

		auto df_22 = [](double *u){return GAMMA*u[1]/u[0];};
		dflux[2][2] = df_22;


		auto wavespeed_0 = [](double *u){return u[1]/u[0] + sqrt(GAMMA*((GAMMA-1)*(u[2] - 0.5*pow(u[1],2)/u[0]))/u[0]);};
		wavespeed[0] = wavespeed_0;

		auto wavespeed_1 = [](double *u){return u[1]/u[0] + sqrt(GAMMA*((GAMMA-1)*(u[2] - 0.5*pow(u[1],2)/u[0]))/u[0]);};
		wavespeed[1] = wavespeed_1;

		auto wavespeed_2 = [](double *u){return u[1]/u[0] + sqrt(GAMMA*((GAMMA-1)*(u[2] - 0.5*pow(u[1],2)/u[0]))/u[0]);};
		wavespeed[2] = wavespeed_2;



	}
};


#endif
