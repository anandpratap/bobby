#include "../src/common.h"
#include "../src/equations.h"
#include "common.h"

void test_jacobian(Equation *equation){
	double dq = 1e-6;
	int nvar = equation->nvar;
	double *q = new double[nvar]();
	double *F = new double[nvar]();
	double *Fp = new double[nvar]();

	for(int ivar=0; ivar < nvar; ivar++){
		q[ivar] = 200.0;
	}
	
	double **dFdq_exact = allocate_2d_array<double>(nvar, nvar);
	double **dFdq_fd = allocate_2d_array<double>(nvar, nvar);
	
	for(int ivar=0; ivar < nvar; ivar++){
		for(int jvar=0; jvar < nvar; jvar++){
			std::cout<<ivar<<" "<<jvar<<std::endl;
			dFdq_exact[ivar][jvar] = equation->dflux[ivar][jvar](q);
		}
	}

	for(int ivar=0; ivar < nvar; ivar++){
		F[ivar] = equation->flux[ivar](q);
	}

	for(int ivar=0; ivar < nvar; ivar++){
		q[ivar] += dq;
		std::cout<<q[ivar]<<std::endl;

		for(int jvar=0; jvar < nvar; jvar++){
			Fp[jvar] = equation->flux[jvar](q);
			dFdq_fd[jvar][ivar] = (Fp[jvar] - F[jvar])/dq;
		}
		q[ivar] -= dq;
	}


	for(int ivar=0; ivar < nvar; ivar++){
		for(int jvar=0; jvar < nvar; jvar++){
			assert_almost_equal(dFdq_fd[ivar][jvar], dFdq_exact[ivar][jvar], 5);
		}
	}

	release_2d_array(dFdq_exact, nvar, nvar);
	release_2d_array(dFdq_fd, nvar, nvar);
	delete[] q;
	delete[] F;
	delete[] Fp;
}

int main(void){
	Equation equation = Equation();
	equation.setup();
	test_jacobian(&equation);

	BurgersEquation burgers_equation = BurgersEquation();
	burgers_equation.setup();
	test_jacobian(&burgers_equation);


	
	SystemEquation system_equation = SystemEquation();
	system_equation.setup();
	test_jacobian(&system_equation);

	return 0;
}
