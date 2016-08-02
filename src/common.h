#ifndef _COMMON_H
#define _COMMON_H

#include <cstdarg>
#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>
#include <armadillo>
#include <stdio.h>

typedef double (*func_pointer)(double *);
typedef double (*func_double)(double);

template <class T>
T*** allocate_3d_array(int nx, int ny, int nz){
    T*** A = new T**[nx];
    for(int i(0); i < nx; ++i){
			A[i] = new T*[ny];
			for(int j(0); j < ny; ++j){
				A[i][j] = new T[nz];					
				for(int k(0); k < nz; ++k){
					A[i][j][k]= 0.;
				}
			}
	}
    return A;
}
template <class T>
void release_3d_array(T*** A, int nx, int ny, int nz){
    for (int i = 0; i < nx; ++i){
			for (int j = 0; j < ny; ++j){
				delete[] A[i][j];
			}
			delete[] A[i];
	}
    delete[] A;
}

template <class T>
T** allocate_2d_array(int nx, int ny){
    T** A = new T*[nx];
    for(int i(0); i < nx; ++i){
		A[i] = new T[ny];
	}
    return A;
}

template <class T>
void release_2d_array(T** A, int nx, int ny){
    for (int i = 0; i < nx; ++i){
		delete[] A[i];
	}
    delete[] A;
}


template<class T>
void array_set_values(int size, T *array, T val){
	for(int i=0; i<size; i++){
		array[i] = val;
	}
}

template<class T>
void array_set_values(int ni, int nj, T **array, T val){
	for(int i=0; i<ni; i++){
		for(int j=0; j<nj; j++){
			array[i][j] = val;
		}
	}
}

template<class T>
void array_linspace(int size, T *array, T start, T end){
	T darray = (end - start)/(size - 1);
	array[0] = start;
	for(int i=1; i<size; i++){
		array[i] = array[i-1] + darray;
	}
}


template<class T>
void array_set_values(int ni, int nj, int nk, T ***array, T val){
for(int i=0; i<ni; i++){
for(int j=0; j<nj; j++){
for(int k=0; k<nk; k++){
				array[i][j][k] = val;
			}
		}
	}
}

#endif
