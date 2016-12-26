/*
 Routines used in simulation.cpp to perform Newton-Raphson method for solving PDEs numerically 
 */
#include <iostream>
#include <math.h>

using namespace std;

void elimination(double **, double *, double *, const int&);    // Gaussan elimination of 5 diagonal matrix, gives q out of Aq=r

void calculate_matrix_A(double **, const int&, double *, double *, double *, const double&, const double&, const double&);

double calculate_volume(const int&, double *, const double&);

void check_volume(const int&, double&, double *, double *, const double&, const double&);

int check_second_time_derivative(const int&, double *, double *, double *, double& , double, int&);

//int check_newton_convergence(const int&, double *, double *, double *, double *, const double&, int *, int&, const int&, double&);



int check_newton_convergence(const int&, double *, double *, double *, double&, const double&, int&, int&, const int&, double&);


