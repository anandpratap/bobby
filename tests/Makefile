CFLAGS=-std=c++11 -O0 -g

CXX=g++

LIB=-L/home/asingh/local/armadillo/usr/lib -larmadillo
INC=-I/home/asingh/local/armadillo/usr/include

all:
	${CXX} bobby.cpp -o test_bobby $(CFLAGS) $(LIB) $(INC)
	${CXX} quadratures.cpp -o test_quadratures $(CFLAGS) $(LIB) $(INC)
	${CXX} equations.cpp -o test_equations $(CFLAGS) $(LIB) $(INC)

clean:
	rm *.o
