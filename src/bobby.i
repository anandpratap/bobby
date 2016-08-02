%module bobby
%{
#define SWIG_FILE_WITH_INIT
#include "bobby.h"
%}

%include "numpy.i"
%init %{
	import_array();
%}

%include "bobby.h"
