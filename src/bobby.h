#ifndef _BOBBY_H
#define _BOBBY_H

#include "common.h"
#include "quadratures.h"
#include "shapefunctions.h"
const double __eval_points[2] = {-sqrt(1.0/3.0), sqrt(1.0/3.0)};
const double __weights[2] = {1.0, 1.0};
const unsigned int __neval_points = 2;

class Bobby{
 public:
	int n, nvar, nt;
	
	double tf, cfl, dt;
	int periodic, step_count_max, implicit;

	double *x;
	double *q, *q_old, *dqdt, *dq;
	double *global_rhs, **global_lhs;
	double **global_mass, **global_linearization;
	
	double ***local_mass, ***local_linearization, **local_rhs;
	arma::mat lhs_arma;
	arma::mat rhs_arma, dq_arma;

	ShapeFunction1D *shapefunction;
	GaussQuadrature1D *quadrature;
	Equation *equation;

	void allocate(){
		x = new double[n]();
		q = new double[nt]();
		q_old = new double[nt]();
		dqdt = new double[nt]();
		dq = new double[nt]();
		global_rhs = new double[nt]();

		global_lhs = allocate_2d_array<double>(nt, nt);
		global_mass = allocate_2d_array<double>(nt, nt);
		global_linearization = allocate_2d_array<double>(nt, nt);

		
		local_rhs = allocate_2d_array<double>(nt-1,2*nvar);
		local_mass = allocate_3d_array<double>(nt-1,2*nvar,2*nvar);
		local_linearization = allocate_3d_array<double>(nt-1,2*nvar,2*nvar);
		
		lhs_arma = arma::mat(nt, nt);
		rhs_arma = arma::mat(nt, 1);
		dq_arma = arma::mat(nt, 1);
	   	}

	void dellocate(){
		delete[] x;
		delete[] q; delete[] q_old; delete[] dqdt;
		delete[] dq; delete[] global_rhs;
		release_2d_array(global_lhs, nt, nt);
		release_2d_array(global_mass, nt, nt);
		release_2d_array(global_linearization, nt, nt);
		release_2d_array(local_rhs, nt-1, 2*nvar);
		release_3d_array(local_mass, nt-1, 2*nvar, 2*nvar);
		release_3d_array(local_linearization, nt-1, 2*nvar, 2*nvar);

	}

	double get_dxdchi(unsigned int el){
		double h = x[el+1] - x[el];
		return h/2.0;
	}

	double get_dchidx(unsigned int el){
		return 1.0/get_dxdchi(el);
	}

	double get_volume(unsigned int el){
		double h = x[el+1] - x[el];
		return h/2.0;
	}
	
	double calc_dt(void){
		return 1e-3;
	}

	double get_local_variable(double *var, unsigned int el, unsigned int local_node, unsigned int ivar){
		return var[(el+local_node)*nvar+ivar];
	}

	
	void global_assembly(unsigned int el){
		for(int ivar=0; ivar<nvar; ivar++){
			for(int jvar=0; jvar<nvar; jvar++){
				global_mass[el*nvar+ivar][el*nvar+jvar] += local_mass[el][0+2*ivar][0+2*jvar];
				global_mass[el*nvar+ivar][(el+1)*nvar+jvar] += local_mass[el][0+2*ivar][1+2*jvar];
				global_mass[(el+1)*nvar+ivar][el*nvar+jvar] += local_mass[el][1+2*ivar][0+2*jvar];
				global_mass[(el+1)*nvar+ivar][(el+1)*nvar+jvar] += local_mass[el][1+2*ivar][1+2*jvar];
				/* global_linearization[el*nvar+ivar][el*nvar+jvar] += local_linearization[el][0+2*ivar][0+2*jvar]; */
				/* global_linearization[el*nvar+ivar][(el+1)*nvar+jvar] += local_linearization[el][0+2*ivar][1+2*jvar]; */
				/* global_linearization[(el+1)*nvar+ivar][el*nvar+jvar] += local_linearization[el][1+2*ivar][0+2*jvar]; */
                /* global_linearization[(el+1)*nvar+ivar][(el+1)*nvar+jvar] += local_linearization[el][1+2*ivar][1+2*jvar]; */
			}
			global_rhs[el*nvar + ivar] += local_rhs[el][0 + 2*ivar];
            global_rhs[(el+1)*nvar + ivar] += local_rhs[el][1 + 2*ivar];
		}
	}

	void step(){
		array_set_values(nt, nt, global_mass, 0.0);
		array_set_values(nt, nt, global_linearization, 0.0);
		array_set_values(nt, nt, global_lhs, 0.0);
		array_set_values(nt, global_rhs, 0.0);
		
		array_set_values(nt-1,2*nvar,2*nvar, local_mass, 0.0);
		array_set_values(nt-1,2*nvar, local_rhs, 0.0);

		for(int el=0; el<n-1; el++){
			elemental(el);
		}
		
		for(int el=0; el<n-1; el++){
			global_assembly(el);
		}


	}

	void explicit_step(){
		step();
		unsigned int start = nvar;
		unsigned int end = nt-1;
		for(int i=0; i<nt; i++){
			for(int j=0; j<nt; j++){
				lhs_arma(i,j) = global_mass[i][j];
			}
			rhs_arma(i,0) = global_rhs[i];
		}
		dq_arma = arma::solve(lhs_arma, rhs_arma);
		for(int i=0; i<nt; i++){
			q[i] += dq_arma(i)*dt;
		}
		std::cout<<arma::norm(dq_arma)<<std::endl;
	}

	void elemental(unsigned int el){
		for(unsigned int quad_idx=0; quad_idx<2; quad_idx++){
			elemental_quad(el, quad_idx);
		}
	}

	void elemental_quad(unsigned int el, unsigned int quad_idx){
		double weight =  __weights[quad_idx];
		double chi = __eval_points[quad_idx];
		double dchidx = get_dchidx(el);
		double dxdchi = get_dxdchi(el);
		
		double **q_local;
		q_local = allocate_2d_array<double>(2, nvar);

		double **dqdt_local;
		dqdt_local = allocate_2d_array<double>(2, nvar);

		for(int local_node=0; local_node<2; local_node++){
			for(int ivar=0; ivar<nvar; ivar++){
				q_local[local_node][ivar] = get_local_variable(q, el, local_node, ivar);
				dqdt_local[local_node][ivar] = get_local_variable(dqdt, el, local_node, ivar);
			}
		}
		
		auto N_0 = shapefunction->value(chi, 0);
		auto N_1 = shapefunction->value(chi, 1);

		auto dN_0 = shapefunction->derivative(chi, 0)*dchidx;
		auto dN_1 = shapefunction->derivative(chi, 1)*dchidx;

		
		for(int ivar=0; ivar<nvar; ivar++){
			local_mass[el][0+2*ivar][0+2*ivar] += N_0*N_0*dxdchi*weight; // integral
			local_mass[el][0+2*ivar][1+2*ivar] += N_0*N_1*dxdchi*weight; // integral
			local_mass[el][1+2*ivar][0+2*ivar] += N_1*N_0*dxdchi*weight; // integral
			local_mass[el][1+2*ivar][1+2*ivar] += N_1*N_1*dxdchi*weight; // integral
		}

		
		double *q_chi = new double[nvar]();
		double *fq_chi = new double[nvar]();
		double *residual_chi = new double[nvar]();
		
		for(int ivar = 0; ivar< nvar; ivar++){
			q_chi[ivar] = q_local[0][ivar]*N_0 + q_local[1][ivar]*N_1;
		}

		for(int ivar = 0; ivar< nvar; ivar++){
			fq_chi[ivar] = equation->flux[ivar](q_chi);
		}

		for(int ivar=0; ivar<nvar; ivar++){
			local_rhs[el][0+2*ivar] += dN_0*fq_chi[ivar]*dxdchi*weight; // integral
			local_rhs[el][1+2*ivar] += dN_1*fq_chi[ivar]*dxdchi*weight; // integral
		}

		for(int ivar = 0; ivar< nvar; ivar++){
			for(int jvar = 0; jvar< nvar; jvar++){
				residual_chi[ivar] += (q_local[0][jvar]*dN_0 + q_local[1][jvar]*dN_1)*equation->dflux[ivar][jvar](q_chi);
			}
		}

			//equation->get_tau(q_chi, dt, dchidx, kvar);
		

		for(int ivar = 0; ivar< nvar; ivar++){
			double tau_dt_term = 4.0/dt/dt;
			double tau_advection_term = 4.0*pow(equation->wavespeed[ivar](q_chi), 2.0)*(dchidx*dchidx);
			double tau_chi =  1.0/sqrt(tau_dt_term + tau_advection_term);
				
			for(int jvar = 0; jvar< nvar; jvar++){
				local_rhs[el][0+2*ivar] -= equation->dflux[ivar][jvar](q_chi)*dN_0*tau_chi*residual_chi[ivar]*dxdchi*weight; // integral
				local_rhs[el][1+2*ivar] -= equation->dflux[ivar][jvar](q_chi)*dN_1*tau_chi*residual_chi[ivar]*dxdchi*weight; // integral
			}
		}

		release_2d_array(q_local, 2, nvar);
		release_2d_array(dqdt_local, 2, nvar);
		delete[] q_chi;
		delete[] fq_chi;
		delete[] residual_chi;
	}

	void write(){
		FILE *fp;
		fp = fopen("out.dat", "w");
		for(int i=0; i<nt; i++){
			fprintf(fp, "%.14e %.14e %.14e\n", x[i], q[i], global_rhs[i]);
		}
		fclose(fp);
	}

	void run(){
		dt = 1e-3;
		array_linspace(n, x, 0.0, 1.0);
		for(int i=0; i<nt; i++){
			if(fabs(x[i]-0.5) < 0.3){
				q[i] = 1.0;
			}
		}
		BurgersEquation eq = BurgersEquation();
		equation = &eq;
		equation->setup();
		std::cout<<"asdasdaD"<<std::endl;
		ShapeFunction1D sh = ShapeFunction1D();
		shapefunction = &sh;
		for(int i=0; i<100; i++){
			std::cout<<i<<std::endl;
			explicit_step();
			write();
		}
	}

	Bobby(){
		nvar = 1;
		n = 100;
		nt = n*nvar;
		allocate();
	}

	~Bobby(){
		dellocate();
	}
};

#endif