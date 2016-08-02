#ifndef _BOBBY_H
#define _BOBBY_H

#include "common.h"
#include "quadratures.h"
#include "shapefunctions.h"
#include "equations.h"

const double __eval_points[2] = {-sqrt(1.0/3.0), sqrt(1.0/3.0)};
const double __weights[2] = {1.0, 1.0};
const unsigned int __neval_points = 2;
class Config{
 private:
	int periodic, step_count_max, implicit;
	double cfl, tf;
	
	Config(){
		periodic = 0;
		step_count_max = 100000;
		implicit = 0;
		cfl = 0.8;
		tf = 0.0;
	}
 public:
	void set_periodic(int val_periodic){
		periodic = val_periodic;
	}
};
class Bobby{
 public:
	int n, nvar, nt;
	int ndim;
	int nnode, nlocal_node;
	double tf, cfl, dt;
	int periodic, step_count_max, implicit;


	int **elm;
	int **bcm;

	
	double *x;
	double *q, *q_old, *dqdt, *dq;
	double *global_rhs, **global_lhs;
	double **global_mass, **global_linearization;
	
	double ***local_mass, ***local_linearization, **local_rhs;
	arma::sp_mat lhs_arma;
	arma::mat rhs_arma, dq_arma, q_arma;

	ShapeFunction1D *shapefunction;
	GaussQuadrature1D *quadrature;
	Equation *equation;

	void allocate(){
		x = new double[n*ndim]();
		q = new double[nt]();
		q_old = new double[nt]();
		dqdt = new double[nt]();
		dq = new double[nt]();
		global_rhs = new double[nt]();

		elm = allocate_2d_array<int>(n-1, nlocal_node);
		bcm = allocate_2d_array<int>(n, 3);

		global_lhs = allocate_2d_array<double>(nt, nt);
		global_mass = allocate_2d_array<double>(nt, nt);
		global_linearization = allocate_2d_array<double>(nt, nt);

		
		local_rhs = allocate_2d_array<double>(nt-1,2*nvar);
		local_mass = allocate_3d_array<double>(nt-1,2*nvar,2*nvar);
		local_linearization = allocate_3d_array<double>(nt-1,2*nvar,2*nvar);
		
		lhs_arma = arma::sp_mat(nt, nt);
		rhs_arma = arma::mat(nt, 1);
		dq_arma = arma::mat(nt, 1);
		q_arma = arma::mat(nt, 1);
	   	}

	void deallocate(){
		delete[] x;
		delete[] q; delete[] q_old; delete[] dqdt;
		delete[] dq; delete[] global_rhs;
		release_2d_array(global_lhs, nt, nt);
		release_2d_array(global_mass, nt, nt);
		release_2d_array(global_linearization, nt, nt);
		release_2d_array(local_rhs, nt-1, 2*nvar);
		release_3d_array(local_mass, nt-1, 2*nvar, 2*nvar);
		release_3d_array(local_linearization, nt-1, 2*nvar, 2*nvar);
		release_2d_array(elm, n-1, nlocal_node);
		release_2d_array(bcm, n, 3);

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
		dt = 1e0;
		double q__[nvar];
		double dxmin = 1e10;
		
		for(int i=0; i<n-1; i++){
			dxmin = fmin(dxmin, x[i+1] - x[i]);
		}
		
		for(int i=0; i<n; i++){
			for(int ivar=0; ivar<nvar; ivar++){q__[ivar] = q[i*nvar+ivar];}
			for(int ivar=0; ivar<nvar; ivar++){
				double w = equation->wavespeed[ivar](q__);
				dt = fmin(dt, dxmin/fabs(w)*cfl);
			}
		}
		std::cout<<"dt == "<<dt<<std::endl;
		return dt;
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
				
				global_linearization[el*nvar+ivar][el*nvar+jvar] += local_linearization[el][0+2*ivar][0+2*jvar];
				global_linearization[el*nvar+ivar][(el+1)*nvar+jvar] += local_linearization[el][0+2*ivar][1+2*jvar];
				global_linearization[(el+1)*nvar+ivar][el*nvar+jvar] += local_linearization[el][1+2*ivar][0+2*jvar];
                global_linearization[(el+1)*nvar+ivar][(el+1)*nvar+jvar] += local_linearization[el][1+2*ivar][1+2*jvar];
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
		
		
		array_set_values(nt-1,2*nvar,2*nvar, local_linearization, 0.0);
		array_set_values(nt-1,2*nvar,2*nvar, local_mass, 0.0);
		array_set_values(nt-1,2*nvar, local_rhs, 0.0);

		for(int el=0; el<n-1; el++){
			elemental(el);
		}
		
		for(int el=0; el<n-1; el++){
			global_assembly(el);
		}

		// bc
		double *qb = new double[nvar]();

		// rhs

		for(int ivar=0; ivar<nvar; ivar++){
			for(int j=0; j<nvar; j++){qb[j] = q[(n-1)*nvar+j];}
			if(ivar == 0 || ivar == 1){
				global_rhs[nt-ivar-1] -= equation->flux[ivar](qb);
				global_linearization[nt-ivar-1][nt-ivar-1] -=  equation->dflux[ivar][ivar](qb);
			}
			else{
				global_rhs[nt-ivar-1] = qb[ivar] - 0.1/(GAMMA-1.0) - 0.5*qb[1]*qb[1]/qb[0];
				global_linearization[nt-ivar-1][nt-0-1] =  -0.5*qb[1]*qb[1]/qb[0]/qb[0];
				global_linearization[nt-ivar-1][nt-1-1] =  -qb[1]/qb[0];
				global_linearization[nt-ivar-1][nt-2-1] =  1.0;
			}

		}


		// lhs

		for(int ivar=0; ivar<nvar; ivar++){
			for(int j=0; j<nvar; j++){qb[j] = q[j];}
			if(ivar == 0 || ivar == 1){
				global_rhs[ivar] += equation->flux[ivar](qb);
				global_linearization[ivar][ivar] +=  equation->dflux[ivar][ivar](qb);
			}
			else{
				global_rhs[ivar] = qb[ivar] - 1.0/(GAMMA-1.0) - 0.5*qb[1]*qb[1]/qb[0];
				global_linearization[ivar][0] =  -0.5*qb[1]*qb[1]/qb[0]/qb[0];
				global_linearization[ivar][1] =  -qb[1]/qb[0];
				global_linearization[ivar][2] =  1.0;
			}
		}

		// dirchlet
		

		delete[] qb;
}


	void newton_update(){
		int start = nvar;
		int end = nt-nvar;

		for(int i=0; i<nt; i++){
			rhs_arma(i,0) = global_rhs[i];
		}

		for(int i=0; i<nt; i++){
			for(int j=0; j<nt; j++){
				double tmp = -global_mass[i][j]/dt + global_linearization[i][j];
				if(fabs(tmp) > 1e-14)
					lhs_arma(i,j) = tmp;
				rhs_arma(i,0) = rhs_arma(i,0) - global_mass[i][j]*(q[j] - q_old[j])/dt;
			}
		}

		dq_arma = arma::spsolve(lhs_arma, rhs_arma);
		for(int i=start; i<end; i++){
			q[i] = q[i] - dq_arma(i);
			q_arma(i) = q[i];
		}
		std::cout<<arma::norm(dq_arma)/arma::norm(q_arma)<<std::endl;

	}	

	void implicit_step(){
		for(int i=0; i<nt; i++){
			q_old[i] = q[i];
		}


		for(int k=0; k<10; k++){
			step();
			newton_update();
			if(arma::norm(dq_arma)/arma::norm(q_arma) < 1e-2){
				break;
			}
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
		dq_arma = arma::spsolve(lhs_arma, rhs_arma);
		for(int i=start; i<nt; i++){
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
		double dchidx = get_dchidx(el);
		double dxdchi = get_dxdchi(el);
		double volume = get_volume(el);
		double weight =  __weights[quad_idx]*volume;
		double chi = __eval_points[quad_idx];

		double N[nlocal_node] = {shapefunction->value(chi, 0), shapefunction->value(chi, 1)};
		double dN[nlocal_node] = {shapefunction->derivative(chi, 0)*dchidx, shapefunction->derivative(chi, 1)*dchidx};

		
		double q_local[nlocal_node][nvar];
		double dqdt_local[nlocal_node][nvar];
		
		for(int local_node=0; local_node<nlocal_node; local_node++){
			for(int ivar=0; ivar<nvar; ivar++){
				q_local[local_node][ivar] = get_local_variable(q, el, local_node, ivar);
				dqdt_local[local_node][ivar] = get_local_variable(dqdt, el, local_node, ivar);
			}
		}
		

		
		for(int ivar=0; ivar<nvar; ivar++){
			for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
				for(int jlocal_node=0; jlocal_node<nlocal_node; jlocal_node++){
					local_mass[el][ilocal_node+nlocal_node*ivar][jlocal_node+nlocal_node*ivar] += N[ilocal_node]*N[jlocal_node]*weight; // integral
				}
			}
		}

		
		double q_chi[nvar] = {0.0};
		double fq_chi[nvar] = {0.0};
		double residual_chi[nvar] = {0.0};
		
		for(int ivar = 0; ivar< nvar; ivar++){
			q_chi[ivar] = q_local[0][ivar]*N[0] + q_local[1][ivar]*N[1];
		}

		for(int ivar = 0; ivar< nvar; ivar++){
			fq_chi[ivar] = equation->flux[ivar](q_chi);
		}

		for(int ivar=0; ivar<nvar; ivar++){
			for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
				local_rhs[el][ilocal_node+nlocal_node*ivar] += dN[ilocal_node]*fq_chi[ivar]*weight; // integral

			}//			local_rhs[el][1+nlocal_node*ivar] += dN[1]*fq_chi[ivar]*weight; // integral
		}

		for(int ivar = 0; ivar< nvar; ivar++){
			for(int jvar = 0; jvar< nvar; jvar++){
				for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
					residual_chi[ivar] += q_local[ilocal_node][jvar]*dN[ilocal_node]*equation->dflux[ivar][jvar](q_chi);
				}
			}
		}

			//equation->get_tau(q_chi, dt, dchidx, kvar);
		

		for(int ivar = 0; ivar< nvar; ivar++){
			double tau_chi =  equation->get_tau(q_chi, dt, ivar, dchidx);

			double A_ = 0.0;
			for(int kvar=0; kvar < nvar; kvar++){
				A_ += equation->dflux[ivar][kvar](q_chi);
			}

			for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
				local_rhs[el][ilocal_node+nlocal_node*ivar] -= A_*dN[ilocal_node]*tau_chi*residual_chi[ivar]*weight; // integral
			}
		}

		if(implicit){
			for(int ivar=0; ivar<nvar; ivar++){
				for(int jvar=0; jvar<nvar; jvar++){
					for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
						for(int jlocal_node=0; jlocal_node<nlocal_node; jlocal_node++){
							local_linearization[el][ilocal_node+nlocal_node*ivar][jlocal_node+nlocal_node*jvar] 
								+= dN[ilocal_node]*equation->dflux[ivar][jvar](q_chi)*N[jlocal_node]*weight; //integral

							double tau_chi =  equation->get_tau(q_chi, dt, ivar, dchidx);
							double A_ = 0.0;
							for(int kvar=0; kvar < nvar; kvar++){
								A_ += equation->dflux[ivar][kvar](q_chi);
							}
							
							local_linearization[el][ilocal_node+nlocal_node*ivar][jlocal_node+nlocal_node*jvar] 
								-= A_*dN[ilocal_node]*tau_chi*equation->dflux[ivar][jvar](q_chi)*dN[jlocal_node]*weight; //integral
						}
					}

				}
			}

		}
	}

	
	void set_solution(int val_ni, double *val_qi){
		assert(val_ni == nt);
		for(int i=0; i<nt; i++){
			q[i] = val_qi[i];
		}
	}
	
	void get_solution(int val_no, double *val_qo){
		assert(val_no == nt);
		for(int i=0; i<nt; i++){
			val_qo[i] = q[i];
		}
	}

	void get_mesh(int val_no, double *val_xo){
		assert(val_no == n);
		for(int i=0; i<n; i++){
			val_xo[i] = x[i];
		}
	}

	void read_restart(){
		int tmp_ndim, nsize, nelem;
		std::ifstream restart_file("restart.in");
		restart_file >> tmp_ndim;
		restart_file >> nsize;
		restart_file >> nelem;

		ndim = tmp_ndim;
		n = nsize;
		nt = n*nvar;
		allocate();

		int label;
		for(int i=0; i<nsize; i++){
			for(int d=0; d<ndim; d++){
				restart_file >> label;
				restart_file >> x[label+d*ndim];
			}
		}

		for(int i=0; i<nelem; i++){
			restart_file >> label;
			for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
				restart_file >> elm[label][ilocal_node];
			}
		}

		for(int i=0; i<nsize; i++){
			restart_file >> label;
			for(int ilocal=0; ilocal<2; ilocal++){
				restart_file >> bcm[label][ilocal];
			}
		}



		restart_file.close();

		
		
	}

	void write(){
		std::ofstream solution_file("solution.dat");
		solution_file.setf(std::ios::fixed);
		solution_file.precision(14);
		for(int i=0; i<n; i++){
			solution_file << x[i];
			
			for(int ivar =0; ivar<nvar; ivar++){
				solution_file << "\t" << q[i*nvar + ivar];
			}
			solution_file << std::endl;
		}
		solution_file.close();
	}

	void run(){
		cfl = 0.01;
		array_linspace(n, x, 0.0, 1.0);
		for(int i=0; i<n; i++){
			if(x[i] < 0.5){
				q[nvar*i] = 1.0;
				q[nvar*i+1] = 0.0;
				q[nvar*i+2] = 2.5;
			}
			else{
				q[nvar*i] = 0.125;
				q[nvar*i+1] = 0.0;
				q[nvar*i+2] = 0.25;
			}
		}
		EulerEquation eq = EulerEquation();
		equation = &eq;
		equation->setup();
		std::cout<<"asdasdaD"<<std::endl;
		ShapeFunction1D sh = ShapeFunction1D();
		shapefunction = &sh;
		for(int i=0; i<100000; i++){
			calc_dt();
			std::cout<<i<<std::endl;
			if(implicit){
				implicit_step();
			}
			else{
				explicit_step();
			}
			write();
		}
	}

	Bobby(){
		nvar = 3;
		nlocal_node = 2;
		read_restart();
		implicit = 1;

		//allocate();
	}

	~Bobby(){
		deallocate();
	}
};

#endif
