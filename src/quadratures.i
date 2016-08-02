%module quadratures
%{
#define SWIG_FILE_WITH_INIT
#include "quadratures.h"
%}

%include "numpy.i"
%init %{
	import_array();
%}

%include "quadratures.h"
