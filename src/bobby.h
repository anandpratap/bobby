#ifndef _BOBBY_H
#define _BOBBY_H

#include "common.h"
#include "quadratures.h"
#include "shapefunctions.h"
#include "equations.h"
#include "adolc/adolc.h"
#include "adolc/sparse/sparsedrivers.h"

typedef adouble atype;

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
	double tf, cfl, dt, t;
	int step_count;
	int periodic, step_count_max, implicit;

	int nelem;
	int **elm;
	int **bcm;
	int nbface;
	int **bface;
	int *face_tag;

	double *x;
	atype *q, *q_old;
	double *d_q, *d_global_rhs;
	double *dqdt, *dq;
	atype *global_rhs, **global_lhs;
	double **global_mass, **global_linearization;
	
	double ***local_mass, ***local_linearization;
	atype **local_rhs;

	double ***dxdchi;
	double ***dchidx;
	double *volume;
	double *dsmin;
	arma::sp_mat lhs_arma, global_lhs_arma, global_mass_arma, global_linearization_arma;
	arma::mat rhs_arma, dq_arma, q_arma;
	arma::mat dt_arma;
	ShapeFunction1D *shapefunction_boundary;
	ShapeFunction2D *shapefunction;
	GaussQuadrature1D *quadrature_boundary;
	GaussQuadrature2D *quadrature;
	EulerEquation2D *equation;
	adouble azero;

	int nnz;
	unsigned int *rind = nullptr;
	unsigned int *cind = nullptr;
	double *values = nullptr;
	int options[4] = {0,0,0,0};

	void allocate(){
		x = new double[n*ndim]();
		q = new adouble[nt]();
		d_q = new double[nt]();
		d_global_rhs = new double[nt]();

		q_old = new adouble[nt]();
		dqdt = new double[nt]();
		dq = new double[nt]();
		global_rhs = new adouble[nt]();

		elm = allocate_2d_array<int>(nelem, nlocal_node);
		//bcm = allocate_2d_array<int>(n, 3);
		bface = allocate_2d_array<int>(nbface, 2);
		face_tag = new int[nbface]();

		//global_lhs = allocate_2d_array<double>(nt, nt);
		//global_mass = allocate_2d_array<double>(nt, nt);
		//global_linearization = allocate_2d_array<double>(nt, nt);

		dxdchi = allocate_3d_array<double>(nelem,2,2);
		dchidx = allocate_3d_array<double>(nelem,2,2);
		volume = new double[nelem]();
		dsmin = new double[n]();
		array_set_values(n, dsmin, 1e10);
		local_rhs = allocate_2d_array<atype>(nelem,nlocal_node*nvar);
		local_mass = allocate_3d_array<double>(nelem,nlocal_node*nvar,nlocal_node*nvar);
		local_linearization = allocate_3d_array<double>(nelem,nlocal_node*nvar,nlocal_node*nvar);
		
		lhs_arma = arma::sp_mat(nt, nt);
		global_lhs_arma = arma::sp_mat(nt, nt);
		global_mass_arma = arma::sp_mat(nt, nt);
		global_linearization_arma = arma::sp_mat(nt, nt);

		rhs_arma = arma::mat(nt, 1);
		dq_arma = arma::mat(nt, 1);
		q_arma = arma::mat(nt, 1);
		dt_arma = arma::mat(nt, 1);
	   	}

	void deallocate(){
		delete[] x;
		delete[] q; delete[] q_old; delete[] dqdt;
		delete[] dq; delete[] global_rhs;
		//release_2d_array(global_lhs, nt, nt);
		//release_2d_array(global_mass, nt, nt);
		//release_2d_array(global_linearization, nt, nt);
		release_2d_array(local_rhs, nelem, nlocal_node*nvar);
		release_3d_array(local_mass, nelem, nlocal_node*nvar, nlocal_node*nvar);
		release_3d_array(local_linearization, nelem, nlocal_node*nvar, nlocal_node*nvar);
		release_2d_array(elm, nelem, nlocal_node);
		//		release_2d_array(bcm, n, 3);
		release_2d_array(bface, nbface, 2);

		release_3d_array(dxdchi, nelem, 2, 2);
		release_3d_array(dchidx, nelem, 2, 2);
		delete[] volume;
		delete[] face_tag;
		delete[] d_q;
		delete[] d_global_rhs;
	}

	double calc_metrics(){
		double chi[2];
		double x_tmp[4], y_tmp[4];
		double dxdchi_tmp, dxdeta_tmp, dydchi_tmp, dydeta_tmp;
		for(int el=0; el<nelem; el++){
			double ds_min_local = 1e10;
			dxdchi_tmp = 0.0; dxdeta_tmp=0.0;
			dydchi_tmp = 0.0; dydeta_tmp=0.0;
			
			for(int ilocal_node=0; ilocal_node < nlocal_node; ilocal_node++){
				//				std::cout<<ilocal_node<<std::endl;
				x_tmp[ilocal_node] = x[elm[el][ilocal_node]*ndim + 0];
				y_tmp[ilocal_node] = x[elm[el][ilocal_node]*ndim + 1];
			}


			for(int ilocal_node=0; ilocal_node < nlocal_node; ilocal_node++){
				for(int jlocal_node=0; jlocal_node < nlocal_node; jlocal_node++){
					if(ilocal_node != jlocal_node){
						double ds = pow(x_tmp[ilocal_node]-x_tmp[jlocal_node], 2.0) + 
							pow(y_tmp[ilocal_node]-y_tmp[jlocal_node], 2.0);
						ds_min_local = fmin(ds, ds_min_local);
					}
				}
			}

			for(int ilocal_node=0; ilocal_node < nlocal_node; ilocal_node++){
				//	std::cout<<ds_min_local<<std::endl;
				dsmin[elm[el][ilocal_node]] = fmin(ds_min_local, dsmin[elm[el][ilocal_node]]);
			}

			for(int i=0; i<4; i++){
				if(i == 0){chi[0] = -1.0; chi[1]=-1.0;}
				if(i == 1){chi[0] = 1.0; chi[1]=-1.0;}
				if(i == 2){chi[0] = 1.0; chi[1]=1.0;}
				if(i == 3){chi[0] = -1.0; chi[1]=1.0;}
				
				dxdchi_tmp += x_tmp[i]*shapefunction->derivative(chi, i, 0);
				dxdeta_tmp += x_tmp[i]*shapefunction->derivative(chi, i, 1);

				dydchi_tmp += y_tmp[i]*shapefunction->derivative(chi, i, 0);
				dydeta_tmp += y_tmp[i]*shapefunction->derivative(chi, i, 1);
			}
			/* std::cout<<dydeta_tmp<<std::endl; */
			/* std::cout<<dxdchi_tmp<<std::endl; */
			/* std::cout<<dydchi_tmp<<std::endl; */
			/* std::cout<<dxdeta_tmp<<std::endl; */

			volume[el] = dxdchi_tmp*dydeta_tmp - dydchi_tmp*dxdeta_tmp;
			//std::cout<<"vol "<<volume[el]<<std::endl;
			dxdchi[el][0][0] = dxdchi_tmp;
			dxdchi[el][0][1] = dxdeta_tmp;
			dxdchi[el][1][0] = dydchi_tmp;
			dxdchi[el][1][1] = dydeta_tmp;
			
			dchidx[el][0][0] = dydeta_tmp/volume[el];
			dchidx[el][0][1] = -dydchi_tmp/volume[el];
			dchidx[el][1][0] = -dxdeta_tmp/volume[el];
			dchidx[el][1][1] = dxdchi_tmp/volume[el];
		}

		
		/* for(int iglobal_node=0; iglobal_node < n; iglobal_node++){ */
		/* 	std::cout<<dsmin[iglobal_node]<<std::endl; */
		/* } */

	}

	double calc_dt(void){
		dt_arma.fill(1e10);
		atype q_tmp[nvar];
		double dxmin = 1e10;
		atype speed[2];
		for(int iglobal_node=0; iglobal_node<n; iglobal_node++){
			for(int ivar=0; ivar<nvar; ivar++){q_tmp[ivar] = q[iglobal_node*nvar+ivar];}
			for(int ivar=0; ivar<nvar; ivar++){
				equation->calc_wavespeed(q_tmp, speed);
				for(int idim=0; idim<ndim; idim++){
					dt_arma(iglobal_node*nvar+ivar,0) = fmin(dt_arma(iglobal_node*nvar+ivar,0), dsmin[iglobal_node]/fabs(speed[idim].value())*cfl);
				}
			}
			//			std::cout<<dsmin[iglobal_node]<<std::endl;
		}
		//		dt_arma.print();
		dt = dt_arma.min();
		std::cout<<"dt == "<<dt<<std::endl;
		return dt;
	}
template<class T>
	T get_local_variable(T *var, unsigned int el, unsigned int ilocal_node, unsigned int ivar){
		return var[elm[el][ilocal_node]*nvar+ivar];
	}

	
	void global_assembly(unsigned int el){
		for(int ivar=0; ivar<nvar; ivar++){
			for(int jvar=0; jvar<nvar; jvar++){
				for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
					for(int jlocal_node=0; jlocal_node<nlocal_node; jlocal_node++){
						//global_mass[elm[el][ilocal_node]*nvar+ivar][elm[el][jlocal_node]*nvar+jvar] 
						//	+= local_mass[el][ilocal_node+nlocal_node*ivar][jlocal_node+nlocal_node*jvar];
						global_mass_arma(elm[el][ilocal_node]*nvar+ivar,elm[el][jlocal_node]*nvar+jvar)
							+= local_mass[el][ilocal_node+nlocal_node*ivar][jlocal_node+nlocal_node*jvar];


						if(implicit){
							global_linearization_arma(elm[el][ilocal_node]*nvar+ivar,elm[el][jlocal_node]*nvar+jvar)
								+= local_linearization[el][ilocal_node+nlocal_node*ivar][jlocal_node+nlocal_node*jvar];
						}
						//std::cout<<el<<std::endl;
					}
				}

				/* global_mass[el*nvar+ivar][el*nvar+jvar] += local_mass[el][0+2*ivar][0+2*jvar]; */
				/* global_mass[el*nvar+ivar][(el+1)*nvar+jvar] += local_mass[el][0+2*ivar][1+2*jvar]; */
				/* global_mass[(el+1)*nvar+ivar][el*nvar+jvar] += local_mass[el][1+2*ivar][0+2*jvar]; */
				/* global_mass[(el+1)*nvar+ivar][(el+1)*nvar+jvar] += local_mass[el][1+2*ivar][1+2*jvar]; */

				/* global_linearization[el*nvar+ivar][el*nvar+jvar] += local_linearization[el][0+2*ivar][0+2*jvar]; */
				/* global_linearization[el*nvar+ivar][(el+1)*nvar+jvar] += local_linearization[el][0+2*ivar][1+2*jvar]; */
				/* global_linearization[(el+1)*nvar+ivar][el*nvar+jvar] += local_linearization[el][1+2*ivar][0+2*jvar]; */
                /* global_linearization[(el+1)*nvar+ivar][(el+1)*nvar+jvar] += local_linearization[el][1+2*ivar][1+2*jvar]; */
			}
			
			for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
				global_rhs[elm[el][ilocal_node]*nvar + ivar] += local_rhs[el][ilocal_node + nlocal_node*ivar];
			}
		}
	}

	void step(){
		//global_lhs_arma.fill(0.0); // set zero
		//global_mass_arma.fill(0.0); // ze
		//array_set_values(nt, nt, global_mass, 0.0);
		//array_set_values(nt, nt, global_linearization, 0.0);
		//array_set_values(nt, nt, global_lhs, 0.0);
		for(int i=0; i<nt; i++){
			d_q[i] = q[i].value();
		}

		trace_on(1);

		for(int i=0; i<nt; i++){
			q[i] <<= d_q[i];
		}

		global_mass_arma.zeros();
		global_linearization_arma.zeros();
		global_lhs_arma.zeros();
		array_set_values(nt, global_rhs, azero);
		
		
		array_set_values(nelem,nlocal_node*nvar,nlocal_node*nvar, local_linearization, 0.0);
		array_set_values(nelem,nlocal_node*nvar,nlocal_node*nvar, local_mass, 0.0);
		array_set_values(nelem,nlocal_node*nvar, local_rhs, azero);

		for(int el=0; el<nelem; el++){
			elemental(el);
		}

		
		for(int el=0; el<nelem; el++){
			global_assembly(el);
		}

		for(int face=0; face<nbface; face++){
			facets(face);
		}


		for(int i=0; i<nt; i++){
			global_rhs[i] >>= d_global_rhs[i];
		}

		trace_off();


		for(int i=0; i<nt; i++){
			d_global_rhs[i] = global_rhs[i].value();
		}

		// 

		/* for(int inode=0; inode<n; inode++){ */
		/* 	if(fabs(x[inode*ndim+0]+1) < 2e-2 || fabs(x[inode*ndim+0]-1) < 2e-2 || fabs(x[inode*ndim+1]-1) < 2e-2 || fabs(x[inode*ndim+1]+1) < 2e-2 ){ */
		/* 		for(int ivar=0; ivar<nvar;ivar++){ */
		/* 		for(int jnode=0; jnode<n; jnode++){ */
		/* 			global_mass[inode*nvar+ivar][jnode] = 0.0; */
		/* 		}} */

		/* 		global_rhs[inode*nvar+0] = q[inode*nvar+0]- 1.0; */
		/* 		global_mass[inode*nvar+0][inode*nvar+0] = 1.0; */

		/* 		global_rhs[inode*nvar+1] = q[inode*nvar+1] - 0.0; */
		/* 		global_mass[inode*nvar+1][inode*nvar+1] = 1.0; */
				
		/* 		global_rhs[inode*nvar+2] = q[inode*nvar+2] - 0.0; */
		/* 		global_mass[inode*nvar+2][inode*nvar+2] = 1.0; */

		/* 		global_rhs[inode*nvar+3] = q[inode*nvar+3] - 1.0/(GAMMA-1.0); */
		/* 		global_mass[inode*nvar+3][inode*nvar+3] = 1.0; */


		/* 	} */

		/* } */
	}


	void newton_update(){
		int start = 0;
		int end = nt;

		for(int i=0; i<nt; i++){
			rhs_arma(i,0) = global_rhs[i].value();
			global_linearization_arma(i,i) = -1/dt_arma(i,0) + global_linearization_arma(i,i);
		}

		//global_lhs_arma = -global_mass_arma/dt + global_linearization_arma;
		
		/* for(int i=0; i<nt; i++){ */
		/* 	//			for(int j=0; j<nt; j++){ */
		/* 		rhs_arma(i,0) = rhs_arma(i,0);// - global_mass_arma(i,j)*(q[j] - q_old[j])/dt; */
		/* 		//			} */
		/* } */
		
		//-global_mass_arma/dt + global_linearization_arma
		dq_arma = arma::spsolve(-global_linearization_arma, rhs_arma);
		for(int i=start; i<end; i++){
			q[i] = q[i] + dq_arma(i);
			q_arma(i) = q[i].value();
		}
		std::cout<<arma::norm(dq_arma)/arma::norm(q_arma)<<std::endl;
	}	

	void implicit_step(){
		for(int i=0; i<nt; i++){
			q_old[i] = q[i];
		}


		for(int k=0; k<1; k++){
			step();
			newton_update();
			if(arma::norm(dq_arma)/arma::norm(q_arma) < 1e-2){
				break;
			}
		}
	}

	void explicit_step(){
		step();
		unsigned int start = 0;
		unsigned int end = nt;
		for(int i=0; i<nt; i++){
			//for(int j=0; j<nt; j++){
			//lhs_arma(i,i) = global_mass[i][i];
				//}
			rhs_arma(i,0) = global_rhs[i].value();
		}
		if(step_count % 1 == 0){
			sparse_jac(1,nt,nt,0,d_q,&nnz,&rind,&cind,&values,options);
		for(int i=0; i<nnz; i++){
			lhs_arma(rind[i],cind[i]) = -values[i];
		}
		
		for(int i=0; i<nt; i++){
			lhs_arma(i,i) += 1.0/dt_arma(i,0);
			//lhs_arma(i,i) += 1.0/dt;
		}

		free(rind); rind=nullptr;
		free(cind); cind=nullptr;
		free(values); values=nullptr;
		}

		//		lhs_arma.print();
		//rhs_arma.print();
		dq_arma = arma::spsolve(lhs_arma, rhs_arma);
		for(int i=0; i<nt; i++){
			//if(fabs(x[i/nvar*ndim+0]+1) < 1e-1 || fabs(x[i/nvar*ndim+0]-1) < 1e-1 || fabs(x[i/nvar*ndim+1]-1) < 1e-1 || fabs(x[i/nvar*ndim+1]+1) < 1e-1 ){
			//}else{
			//				std::cout<<x[i/nvar*ndim+0]<<" "<<x[i/nvar*ndim+1]<<std::endl;
			q[i] += 0.5*dq_arma(i);//*dt_arma(i,0);
				//}
		}
		std::cout<<"norm(dq)"<<arma::norm(dq_arma)<<std::endl;
		std::cout<<"norm(rhs)"<<arma::norm(rhs_arma)<<std::endl;

	}

	// void calc_forces(){
	// 	double lift = 0.0;
		
	// 	for(int face=0; face<nbface; face++){
	// 		if(face_tag[face] == 0){
	// 			int node[2];
				
	// 			node[0] = bface[face][0];
	// 			node[1] = bface[face][1];

	// 			double x_tmp[2], y_tmp[2];
				
	// 			for(int i=0; i<2; i++){
	// 				x_tmp[i] = x[node[i]*ndim+0];
	// 				y_tmp[i] = x[node[i]*ndim+1];
	// 			}
	// 			double normal[2];
	// 			double dx = x_tmp[1] - x_tmp[0];
	// 			double dy = y_tmp[1] - y_tmp[0];
	// 			double ds = sqrt(dx*dx + dy*dy);
	// 			normal[0] = dy/ds;
	// 			normal[1] = -dx/ds;
				
	// 			atype q_local[2][nvar];
				
	// 			for(int ilocal_node=0; ilocal_node<2; ilocal_node++){
	// 				for(int ivar=0; ivar<nvar; ivar++){
	// 					q_local[ilocal_node][ivar] = q[bface[face][ilocal_node]*nvar + ivar];
	// 				}
	// 			}
				
	// 			double volume_face = ds;
	// 			static atype *q_chi = new atype[nvar];
	// 			static atype **fq_chi = allocate_2d_array<atype>(ndim, nvar);
	// 			static atype **local_face_rhs = allocate_2d_array<atype>(2, nvar);
	// 			array_set_values(2, nvar, local_face_rhs, azero);
				
	// 			for(unsigned int quad_idx=0; quad_idx<quadrature_boundary->__neval_points; quad_idx++){
	// 				double weight = quadrature_boundary->__weights[quad_idx]*volume_face;
	// 				double *chi = quadrature_boundary->__eval_points[quad_idx];
	// 				double N[2];
	// 				double dN[2];
	// 				for(int i=0; i<2; i++){
	// 					N[i] = shapefunction_boundary->value(chi, i);
	// 					dN[i] = shapefunction_boundary->derivative(chi, i, 0);
	// 				}
					
	// 				for(int ivar=0; ivar<nvar; ivar++){
	// 					q_chi[ivar] = q_local[0][ivar]*shapefunction_boundary->value(chi, 0) + 
	// 						q_local[1][ivar]*shapefunction_boundary->value(chi, 1);
	// 				}
	// 				atype qsq2 = q_chi[1]*q_chi[1] + q_chi[2]*q_chi[2];
	// 				atype p = (GAMMA-1.0)*(q_chi[3] - 0.5*qsq2/q_chi[0]);
	// 				lift += p*normal
	// 			}
	// 		}
	// 	}
	// }
	void elemental(unsigned int el){
		for(unsigned int quad_idx=0; quad_idx<quadrature->__neval_points; quad_idx++){
			elemental_quad(el, quad_idx);
		}
	}


	void facets(unsigned int face){
	
		int node[2];

		node[0] = bface[face][0];
		node[1] = bface[face][1];

		double x_tmp[2], y_tmp[2];

		for(int i=0; i<2; i++){
			
			x_tmp[i] = x[node[i]*ndim+0];
			y_tmp[i] = x[node[i]*ndim+1];
		}
		double normal[2];
		double dx = x_tmp[1] - x_tmp[0];
		double dy = y_tmp[1] - y_tmp[0];
		double ds = sqrt(dx*dx + dy*dy);
		normal[0] = dy/ds;
		normal[1] = -dx/ds;

		atype q_local[2][nvar];

		for(int ilocal_node=0; ilocal_node<2; ilocal_node++){
			for(int ivar=0; ivar<nvar; ivar++){
				q_local[ilocal_node][ivar] = q[bface[face][ilocal_node]*nvar + ivar];
			}
		}

		double volume_face = ds;
		static atype *q_chi = new atype[nvar];
		static atype **fq_chi = allocate_2d_array<atype>(ndim, nvar);
		static atype **local_face_rhs = allocate_2d_array<atype>(2, nvar);
		array_set_values(2, nvar, local_face_rhs, azero);
		
		for(unsigned int quad_idx=0; quad_idx<quadrature_boundary->__neval_points; quad_idx++){
			double weight = quadrature_boundary->__weights[quad_idx]*volume_face;
			double *chi = quadrature_boundary->__eval_points[quad_idx];
			double N[2];
			double dN[2];
			for(int i=0; i<2; i++){
				N[i] = shapefunction_boundary->value(chi, i);
				dN[i] = shapefunction_boundary->derivative(chi, i, 0);
			}

			for(int ivar=0; ivar<nvar; ivar++){
				q_chi[ivar] = q_local[0][ivar]*shapefunction_boundary->value(chi, 0) + 
					q_local[1][ivar]*shapefunction_boundary->value(chi, 1);
			}
			equation->calc_flux(q_chi, fq_chi);

			if(face_tag[face] == 0){
				array_set_values(ndim, nvar, fq_chi, azero);
				atype qsq2 = q_chi[1]*q_chi[1] + q_chi[2]*q_chi[2];
				atype p = (GAMMA-1.0)*(q_chi[3] - 0.5*qsq2/q_chi[0]);
				fq_chi[0][1] = p;
				fq_chi[1][2] = p;
				//std::cout<<"wall bc "<<face<<std::endl;
			}

				for(int ilocal_node=0; ilocal_node<2; ilocal_node++){
					for(int idim=0; idim<ndim; idim++){
						for(int ivar=0; ivar<nvar; ivar++){
							local_face_rhs[ilocal_node][ivar] += N[ilocal_node]*fq_chi[idim][ivar]*normal[idim]*weight;
						}
					}
				}
		}

		
		for(int ilocal_node=0; ilocal_node<2; ilocal_node++){
			for(int ivar=0; ivar<nvar; ivar++){
				//				std::cout<<"ivar "<<ivar<<" inode "<<ilocal_node<<"before "<<global_rhs[bface[face][ilocal_node]*nvar + ivar]<<" "<<local_face_rhs[ilocal_node][ivar]<<std::endl;
				global_rhs[bface[face][ilocal_node]*nvar + ivar] -= local_face_rhs[ilocal_node][ivar]; 
				//				std::cout<<"after "<<global_rhs[bface[face][ilocal_node]*nvar + ivar]<<std::endl;
				//for(int jlocal_node=0; jlocal_node<2; jlocal_node++){
					//global_linearization_arma(bface[face][ilocal_node]*nvar + ivar, bface[face][jlocal_node]*nvar + ivar) -= 
				//}
			}
		}


		if(face_tag[face] == 1){
			for(int ilocal_node=0; ilocal_node<2; ilocal_node++){
				int ivar = 0;
				global_rhs[bface[face][ilocal_node]*nvar + ivar] = q[bface[face][ilocal_node]*nvar + ivar] - 1.0;
				ivar = 1;
				global_rhs[bface[face][ilocal_node]*nvar + ivar] = q[bface[face][ilocal_node]*nvar + ivar] - 0.3;
				ivar = 2;
				global_rhs[bface[face][ilocal_node]*nvar + ivar] = q[bface[face][ilocal_node]*nvar + ivar] - 0.0;
				ivar = 3;
				global_rhs[bface[face][ilocal_node]*nvar + ivar] = q[bface[face][ilocal_node]*nvar + ivar] - 1.0/(GAMMA-1.0);
				
			}
		}
	}
	

	void elemental_quad(unsigned int el, unsigned int quad_idx){
		//		double dchidx = get_dchidx(el);
		//double dxdchi = get_dxdchi(el);
		//double volume = get_volume(el);
		
		double weight = quadrature->__weights[quad_idx]*volume[el];
		double *chi = quadrature->__eval_points[quad_idx];

		double N[nlocal_node] = {shapefunction->value(chi, 0), shapefunction->value(chi, 1), 
								 shapefunction->value(chi, 2),shapefunction->value(chi, 3)};

		double dN[ndim][nlocal_node];
		for(int dim=0; dim<ndim; dim++){
			for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
				//std::cout<<dim<<" "<<ilocal_node<<" "<<dchidx[el][dim][0]<<std::endl;
				dN[dim][ilocal_node] = shapefunction->derivative(chi, ilocal_node, 0)*dchidx[el][dim][0] 
					+ shapefunction->derivative(chi, ilocal_node, 1)*dchidx[el][dim][1];
			}
		}
		
		atype q_local[nlocal_node][nvar];
		atype dqdt_local[nlocal_node][nvar];
		
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
					//					std::cout<<local_mass[el][ilocal_node+nlocal_node*ivar][jlocal_node+nlocal_node*ivar]<<std::endl;
				}
			}
		}
		
		static atype *residual_chi = new atype[nvar];
		static atype *q_chi = new atype[nvar];
		static atype **fq_chi = allocate_2d_array<atype>(ndim, nvar);
		static atype ***dfq_chi = allocate_3d_array<atype>(ndim, nvar, nvar);
		
		for(int i=0; i<nvar; i++){
			residual_chi[i] = 0.0;
			q_chi[i] = 0.0;
			for(int dim=0;dim<ndim; dim++){
				fq_chi[dim][i] = 0.0;
				for(int j=0; j<nvar; j++){
					dfq_chi[dim][i][j] = 0.0;
				}
			}
		}

		
		for(int ivar = 0; ivar< nvar; ivar++){
			for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
				q_chi[ivar] += q_local[ilocal_node][ivar]*N[ilocal_node];
			}
		}

		
		equation->calc_flux(q_chi, fq_chi);
		equation->calc_dflux(q_chi, dfq_chi);


		for(int ivar=0; ivar<nvar; ivar++){
			for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
				for(int dim=0; dim<ndim; dim++){
					local_rhs[el][ilocal_node+nlocal_node*ivar] += dN[dim][ilocal_node]*fq_chi[dim][ivar]*weight; // integral
				}
			}
		}

		for(int ivar = 0; ivar< nvar; ivar++){
			for(int jvar = 0; jvar< nvar; jvar++){
				for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
					for(int dim=0; dim<ndim; dim++){
						residual_chi[ivar] += q_local[ilocal_node][jvar]*dN[dim][ilocal_node]*dfq_chi[dim][ivar][jvar];
					}
				}
			}
		}

			//equation->get_tau(q_chi, dt, dchidx, kvar);
		

		for(int ivar = 0; ivar< nvar; ivar++){
			atype tau_dt = 1/dt/dt;
			atype tau_u = 0.0;
			atype u_[2] = {q_chi[1]/q_chi[0], q_chi[2]/q_chi[0]};
			for(int idim=0; idim<ndim; idim++){
				for(int jdim=0; jdim<ndim; jdim++){
					for(int kdim=0; kdim<ndim; kdim++){
						atype ujdim = u_[jdim];
						atype ukdim = u_[kdim];
						tau_u += dchidx[el][idim][jdim]*ujdim*dchidx[el][idim][kdim]*ukdim;
					}
				}
			}

			
			atype tau_chi = 1.0/sqrt(tau_dt + tau_u);
			atype A_[ndim] = {0.0};
			for(int kvar=0; kvar < nvar; kvar++){
				for(int idim=0; idim<ndim; idim++){
					A_[idim] += dfq_chi[idim][ivar][kvar];
				}
			}

			for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
				for(int idim=0; idim<ndim; idim++){
					local_rhs[el][ilocal_node+nlocal_node*ivar] -= A_[idim]*dN[idim][ilocal_node]*tau_chi*residual_chi[ivar]*weight; // integral
				}
			}
		}

		if(implicit){
			for(int ivar=0; ivar<nvar; ivar++){
				for(int jvar=0; jvar<nvar; jvar++){
					for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
						for(int jlocal_node=0; jlocal_node<nlocal_node; jlocal_node++){
							for(int idim=0; idim<ndim; idim++){
								local_linearization[el][ilocal_node+nlocal_node*ivar][jlocal_node+nlocal_node*jvar]
									+= dN[idim][ilocal_node]*dfq_chi[idim][ivar][jvar].value()*N[jlocal_node]*weight; //integral
							}
							/* double tau_chi =  equation->get_tau(q_chi, dt, ivar, dchidx); */
							/* double A_ = 0.0; */
							/* for(int kvar=0; kvar < nvar; kvar++){ */
							/* 	A_ += equation->dflux[ivar][kvar](q_chi); */
							/* } */
							
							/* local_linearization[el][ilocal_node+nlocal_node*ivar][jlocal_node+nlocal_node*jvar] */
							/* 	-= A_*dN[ilocal_node]*tau_chi*equation->dflux[ivar][jvar](q_chi)*dN[jlocal_node]*weight; */ //integral
						}
					}

				}
			}

		}
	}
	
	template<class T>
	void set_solution(int val_ni, T *val_qi){
		assert(val_ni == nt);
		for(int i=0; i<nt; i++){
			q[i] = val_qi[i];
		}
	}
	
	template<class T>
	void get_solution(int val_no, T *val_qo){
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
		int tmp_ndim, nsize;
		std::ifstream restart_file("restart.in");
		restart_file >> tmp_ndim;
		restart_file >> nsize;
		restart_file >> nelem;
		restart_file >> nbface;

		ndim = tmp_ndim;
		std::cout<<tmp_ndim<<std::endl;
		std::cout<<nsize<<std::endl;
		std::cout<<nelem<<std::endl;
			
		n = nsize;
		nt = n*nvar;
		allocate();

		int label;
		for(int i=0; i<nsize; i++){
			restart_file >> label;
			for(int d=0; d<ndim; d++){
				restart_file >> x[label*ndim+d];
			}
		}

		for(int i=0; i<nelem; i++){
			restart_file >> label;
			for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
				//std::cout<<ilocal_node<<" "<<label<<std::endl;
				restart_file >> elm[label][ilocal_node];
			}
		}

			
		for(int i=0; i<nbface; i++){
			restart_file >> label;
			for(int ilocal=0; ilocal<2; ilocal++){
				restart_file >> bface[label][ilocal];
			}
			restart_file >> face_tag[label];
		}


		for(int i=0; i<n; i++){
			restart_file >> label;
			for(int ivar=0; ivar<nvar; ivar++){
				restart_file >> d_q[label*nvar+ivar];
			}
		}


		restart_file.close();
		calc_metrics();
	}

	void write_restart(){
		std::string s = std::to_string(step_count);
		std::string ss= "restart_";
		std::string dest = ss.append(std::string(10-s.length(), '0').append(s).append(".out"));
		std::ofstream restart_file(dest);
		restart_file.setf(std::ios::fixed);
		restart_file.precision(14);
		
		restart_file << ndim << std::endl;
		restart_file << n << std::endl;
		restart_file << nelem << std::endl;
		restart_file << nbface << std::endl;
		int label;
		for(int i=0; i<n; i++){
			label = i;
			restart_file << label;
			for(int d=0; d<ndim; d++){
				restart_file <<" "<<x[label*ndim+d];
			}
			restart_file << std::endl;
		}

		for(int i=0; i<nelem; i++){
			label = i;
			restart_file << label;
			for(int ilocal_node=0; ilocal_node<nlocal_node; ilocal_node++){
				//std::cout<<ilocal_node<<" "<<label<<std::endl;
				restart_file <<" "<< elm[label][ilocal_node];
			}
			restart_file << std::endl;
		}

			
		for(int i=0; i<nbface; i++){
			label = i;
			restart_file << label;
			for(int ilocal=0; ilocal<2; ilocal++){
				restart_file <<" "<< bface[label][ilocal];
			}
			restart_file <<" "<< face_tag[label];
			restart_file << std::endl;
		}

		for(int i=0; i<n; i++){
			label = i;
			restart_file << label;
			for(int ivar=0; ivar<nvar; ivar++){
				restart_file <<" "<<q[label*nvar+ivar].value();
			}
			restart_file << std::endl;
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
		cfl = 10.0;
	   
		//array_linspace(n, x, 0.0, 1.0);
		// for(int i=0; i<n; i++){
		// 	double r2 = pow(x[i*ndim+0], 2.0) + pow(x[i*ndim+1], 2.0);
		// 	q[nvar*i] = 1.0;
		// 	q[nvar*i+1] = 0.3;
		// 	double p = 1.0;
		// 	q[nvar*i+2] = 0.0;
		// 	q[nvar*i+3] = p/(GAMMA-1);
			
		// }

		for(int i=0; i<nt; i++){
			q[i] = d_q[i];
		}
		EulerEquation2D eq = EulerEquation2D();
		equation = &eq;
		std::cout<<"asdasdaD"<<std::endl;
		ShapeFunction2D sh = ShapeFunction2D();
		ShapeFunction1D sh1 = ShapeFunction1D();

		GaussQuadrature2D g = GaussQuadrature2D();
		GaussQuadrature1D g1 = GaussQuadrature1D();

		shapefunction = &sh;
		shapefunction_boundary = &sh1;

		quadrature = &g;
		quadrature_boundary = &g1;
		
		write_tecplot();
		write_restart();
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
			write_tecplot();
			write_restart();
			step_count += 1;
			t += dt;
			std::cout<<"Time: "<<t<<std::endl;
		}
	}

	Bobby(){
		nvar = 4;
		nlocal_node = 4;
		step_count = 0;
		read_restart();
		implicit = 0;
		t = 0.0;
		azero = 0.0;
		//allocate();
	}

	~Bobby(){
		deallocate();
	}


	void write_tecplot(){
		FILE *fp;
		std::string s = std::to_string(step_count);
		std::string ss= "tmp_";
		std::string dest = ss.append(std::string(10-s.length(), '0').append(s).append(".tec"));
		std::cout<<dest<<std::endl;
		fp = fopen(dest.c_str(), "w");
		fprintf(fp, "TITLE = \"%s\"\n", "Example file");
		fprintf(fp, "VARIABLES = \"X\", \"Y\",");
		for(int i=0; i<nvar; i++){
			fprintf(fp, "\"Var_%d\",", i);
		}
		fprintf(fp, "\n");
		fprintf(fp, "ZONE T=\"test\", DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FEQUADRILATERAL\n", n, nelem);
		for(int i=0; i<n; i++){
			for(int d=0; d<ndim; d++)
				fprintf(fp, " %.14f", x[i*ndim+d]);
			for(int j=0; j<nvar; j++){
				fprintf(fp, " %.14f", q[i*nvar + j].value());
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		for(int i=0; i<nelem; i++){
			fprintf(fp, "%d %d %d %d\n", elm[i][0]+1,  elm[i][1]+1, elm[i][2]+1, elm[i][3]+1);
		}
		
		fclose(fp);
	}
};

#endif
