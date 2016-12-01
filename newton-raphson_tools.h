/*
 Routines used in simulation.cpp to perform Newton-Raphson method for solving PDEs numerically 
 */
#include <iostream>
#include <math.h>

using namespace std;

void elimination(double **, double *, double *, int);    // Gaussan elimination of 5 diagonal matrix, gives q out of Aq=r

void calculate_matrix_A(double **, int, double *, double *, double *, double, double, double);

double calculate_volume(int, double *, double);

void check_volume(int, double, double *, double *, double, double);

int check_second_time_derivative(int, double *, double *, double *, double *, double, int *);

int check_newton_convergence(int, double *, double *, double *, double *, double, int *, int *, int, double *);

