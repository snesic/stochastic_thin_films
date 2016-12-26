/*
 *  
 * Simulation of a droplet over a flat substrate decribed by the differential equation 

\frac{\partial h}{\partial t}=-\nabla [D(h)\nabla \nabla^2 h] + vdW !

where D(h)=h^3 is a diffusion. Standard sheme is used: D(h_i)=1/2(h_{i+1}+h_{i}). 

The equation is solved using IMPLICIT METHOD, \theta=1/2.  for solving the surface tension part. Gravity and vdW forces are implemented EXPLICITLY as well as the NOISE.

 *  Created by Svetozar Nesic on 26.2.10..
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include "read_write_msg.h"
#include "simulation.h"


#define PI 3.141592653589793

using namespace std;


int main(int argc, char * const argv[])
{   int L, t, t0, i_sim_0;
    double b, A, B, a_min, ugao;     // L System size; dx, dt - Time and Space lattice distances. b - Precursor layer.
    double tol;                  // tol - Tolerance of Newton-Raphson Method, theta - Implicit Mehtod.
    time_t vreme;
    string dir;
    stringstream check_dir, krajnje_vreme;
    struct stat sb;
    long int *inic = new long int;
    double dx, dt, theta=0.5, sigma;
    
    if(argc<3)  dir = "test_dir";
    else dir = argv[2];
    
    i_sim_0 = atoi(argv[1]);
    check_dir << dir << "/" << i_sim_0;
    
    if (stat(check_dir.str().c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
        { cout<< check_dir.str()<< " exists... The folowing simulation will be preformed:" << endl; }            // The simulation starts only if there is a directory!
    else {cout<< "This simulations is canceled as there is no directory to write down results!"<<endl; return 1;}    


L=3*266+1;
dx=0.01;
dt=1e-5;
t0=0;
t=2;
b=0.01;
    ugao=50;
    a_min=1;
    A=2*(1-cos(PI*ugao/180))*(a_min*b); B=A*(a_min*b); // vdW: The min of vdW potential determins B..
    sigma=1e-4;
tol=1e-10;

//*inic=time(NULL)+i_sim_0*500;
*inic=45;
vreme=time(NULL);


print_initial_parameters(sigma, dt, dx, L, b, A, a_min, t, t0, ugao, i_sim_0);
simulation(L, t, dx, dt, t0, b, A, B, sigma, theta, tol, inic, i_sim_0, check_dir.str());
  
    
vreme=time(NULL)-vreme;
krajnje_vreme.str("");
krajnje_vreme << trunc(vreme/3600) << " hours and " << trunc(fmod(vreme,3600)/60) << " min.";

/*  write_email_with_sim_data(L,t,dx,dt,b,A,sigma, check_dir.str(), krajnje_vreme.str(), i_sim_0);
    system("./send_email.sh");
*/
    
cout << "Calculation time is: " << vreme << endl;
cout << "which is " << krajnje_vreme.str() << endl;

return 0;
}

