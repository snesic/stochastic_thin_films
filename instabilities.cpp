/*
 *  
 * Simulation of a droplet over a flat substrate decribed by the differential equation 

\frac{\partial h}{\partial t}=-\nabla [D(h)\nabla \nabla^2 h] + vdW + gravity!

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
#include "film.h"

using namespace std;


int main(int argc, char * const argv[])
{   int L, i_sim_0, j_time;
    double b, ugao;     // L System size; dx, dt - Time and Space lattice distances. b - Precursor layer.
    
    time_t vreme;
    string dir;
    stringstream check_dir, krajnje_vreme;
    struct stat sb;
    long int *inic = new long int;
    double dx, dt, sigma;
    double t, t0, jdt, jdt_no;
    
    if(argc<3){ dir = "test_dir"; }
    else {dir = argv[2];}
    
    i_sim_0 = atoi(argv[1]);
    check_dir << dir << "/" << i_sim_0;
    
    if (stat(check_dir.str().c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
    { cout<< check_dir.str()<< " exists... The folowing simulation will be preformed:" << endl; }            // The simulation starts only if there is a directory!
    else {cout<< "This simulations is canceled as there is no directory to write down results!"<<endl; return 1;}
    
    
    L=3*266+1;
    dx=0.01;
    dt=1e-5;
    t0=0;
    t=7;
    b=0.01;
    ugao=50;
    sigma=1e-4;
    
    *inic=45; // time(NULL)+i_sim_0*500;
    vreme=time(NULL);
    
    
//    print_initial_parameters(sigma, dt, dx, L, b, A, a_min, t, t0, ugao, i_sim_0);

    film droplet(L, dx, dt, b, ugao, sigma);
    initial_condition_sin(droplet.h,L,dx,b);
    
    
    jdt=t0;
    j_time=0;
    jdt_no=t0+dt;

    do
    {   droplet.simulation(jdt_no, jdt, inic);
        cout << "Time step: " << jdt << " dt " << dt << " max:  " << droplet.h[266] << endl;
        write_in_File(droplet.h, L, dx, j_time, check_dir.str());
        jdt=jdt_no;
        if(jdt_no<1) jdt_no=jdt_no*10; else jdt_no=trunc(jdt+1);
        j_time++;
    }while(j_time<t);
    

    
    vreme=time(NULL)-vreme;
    krajnje_vreme.str("");
    krajnje_vreme << trunc(vreme/3600) << " hours and " << trunc(fmod(vreme,3600)/60) << " min.";
    cout << "Calculation time is: " << vreme << endl;
    cout << "which is " << krajnje_vreme.str() << endl;
    
    return 0;
}

