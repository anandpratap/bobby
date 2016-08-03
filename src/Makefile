CXX=g++
CFLAGS=-std=c++17 -O0 -g

LIB=-L${HOME}/local/armadillo/usr/lib -larmadillo
INC=-I${HOME}/local/armadillo/usr/include

INC+=-I${HOME}/local/anaconda2/include/python2.7/
INC+=-I${HOME}/local/anaconda2/lib/python2.7/site-packages/numpy/core/include/

python: pre-build wrapper

pre-build:
#	swig -c++ -python -o shapefunctions_wrap.cpp shapefunctions.i
#	swig -c++ -python -o quadratures_wrap.cpp quadratures.i
	swig -c++ -python -o bobby_wrap.cpp bobby.i

wrapper:
	# $(CXX) -fPIC -c quadratures_wrap.cpp -o quadratures_wrap.o  $(CFLAGS) $(INC)
	# ${CXX} -shared quadratures_wrap.o -o _quadratures.so $(CFLAGS) $(LIB)
	# $(CXX) -fPIC -c shapefunctions_wrap.cpp -o shapefunctions_wrap.o  $(CFLAGS) $(INC)
	# ${CXX} -shared shapefunctions_wrap.o -o _shapefunctions.so $(CFLAGS) $(LIB)
	$(CXX) -fPIC -c bobby_wrap.cpp -o bobby_wrap.o  $(CFLAGS) $(INC)
	${CXX} -shared bobby_wrap.o -o _bobby.so $(CFLAGS) $(LIB)

test: test_bobby test_quadratures test_equations
test_bobby:
	${CXX} ../tests/bobby.cpp $(CFLAGS) $(INC) $(LIB) -o test_bobby
test_quadratures:
	${CXX} ../tests/quadratures.cpp $(CFLAGS) $(INC) $(LIB) -o test_quadratures
test_equations:
	${CXX} ../tests/equations.cpp $(CFLAGS) $(INC) $(LIB) -o test_equations

clean:
	rm *.o
