/*
 Implicit algorithm to solve 2d stochastic lubrication equation (thin films).
 */

#include <iostream>
#include <string>
#include <math.h>

#include "read_write_msg.h"
#include "initial_conditions.h"
#include "newton-raphson_tools.h"
#include "create_forces.h"

using namespace std;

int simulation(int, int, double, double, int, double, double, double, double, double, double, long *, int, string);