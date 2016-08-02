CFLAGS=-std=c++11 -O0

CXX=g++

INC=-I${HOME}/local/anaconda2/include/python2.7/
INC+=-I${HOME}/local/anaconda2/lib/python2.7/site-packages/numpy/core/include/

all:
	swig -c++ -python -o shapefunctions_wrap.cpp shapefunctions.i
	${CXX} -fPIC -c shapefunctions_wrap.cpp -o shapefunctions_wrap.o  $(CFLAGS) $(INC)
	${CXX} -shared shapefunctions_wrap.o -o _shapefunctions.so $(CFLAGS) $(LIB)

	swig -c++ -python -o quadratures_wrap.cpp quadratures.i
	${CXX} -fPIC -c quadratures_wrap.cpp -o quadratures_wrap.o  $(CFLAGS) $(INC)
	${CXX} -shared quadratures_wrap.o -o _quadratures.so $(CFLAGS) $(LIB)

clean:
	rm *.o
