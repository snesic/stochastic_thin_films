/*
 Routines to add noise, surface tension and vdw term to the Newton-Raphson functional
 */

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "random_no_generators.h"

using namespace std;


void create_noise(int, double *, double *, const double&, const double&, const double&, long *);

void create_funcional_st(int, double *, double *, double *h, int *, const double&, const double&, const double&);

void create_funcional_add_vdw_noise(int, double *, double *, double *, int *, const double&, const double&, const double&, const double&);
