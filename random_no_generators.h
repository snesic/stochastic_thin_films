/*
 *  
Random number generetors
Numerical Recipes in C
 *  Created by Svetozar Nesic on 26.2.10..
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __RANDOM_NO_GENERATORS_H_INCLUDED__
#define __RANDOM_NO_GENERATORS_H_INCLUDED__
#include <iostream>
#include <stdlib.h>
#include <math.h>


double ran3(long *);

double gauss(double, double, long *);			/* normal random variate generator */

#endif