/*

 */
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include "read_write_msg.h"

using namespace std;



void initial_condition(double *, int,double, double); // Droplet as: a  b + re^{(x-{x_mid})^4/r^4}

void initial_condition_sin(double *, int, double, double);   //Droplet as a  b + re^{(x-{x_mid})^4/r^4}

int call_initial(double *, double *, int, double, double, int, string);

void create_boundary_conditions(int *, int);
