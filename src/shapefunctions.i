%module shapefunctions
%{
#define SWIG_FILE_WITH_INIT
#include "shapefunctions.h"
%}

%include "numpy.i"
%init %{
	import_array();
%}

%include "shapefunctions.h"
