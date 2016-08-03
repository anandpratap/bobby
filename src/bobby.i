%module bobby
%{
#define SWIG_FILE_WITH_INIT
#include "bobby.h"
%}

%include "numpy.i"
%init %{
	import_array();
%}

%apply (int DIM1, double* IN_ARRAY1) {(int val_ni, double *val_qi)}
%apply (int DIM1, double* ARGOUT_ARRAY1) {(int val_no, double *val_qo)}
%apply (int DIM1, double* ARGOUT_ARRAY1) {(int val_no, double *val_xo)}

%include "bobby.h"
