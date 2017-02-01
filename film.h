#ifndef _FILM_H_INCLUDED__
#define _FILM_H_INCLUDED__
#include <iostream>
#include <algorithm>    // std::swap
#include <numeric>
#include <time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>



#include "read_write_msg.h"
#include "initial_conditions.h"
#include "newton-raphson_tools.h"
#define PI 3.141592653589793


using namespace std;




class film
{
public:
    
    double *h; double dx,dt,b;
    int L;
    
    film(int L=400, double dx=0.02, double dt=1e-5, double b=0.02, int ugao=0, double sigma=0);
    int simulation(double t, double& t0, long int *inic);
        
    
 private:
    
    
    double **a, *y, *q, *r, *r_fix, *h2, *h3, *y2, *y3, *noise, *old_sol;
    int *ind;
    double angle, A, B, sigma;
    double tol=1e-10;
    double theta=0.5;
    
    void set_h2_3();
    void set_y2_3();
    void create_noise(long *inic);
    void create_funcional_st(double *r1, double *h1, double *h13, double theta1);
    void create_funcional_add_vdw_noise();
    

};

#endif
