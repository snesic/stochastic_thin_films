


#include "simulation.h"


int simulation(int L, int t, double dx, double dt, int t0, double b, double A, double B, double sigma, double theta, double tol, long *inic, int i_sim, string dir)
{   int i, j_time, *ind;
    int no_newton_it, no_consecutive_neg, no_good_solutions;
    double jdt, jdt_no,dt_old;
    double *h, *y, *q, *r, *r_fix, *h2, *h3, *y2, *y3, *noise, *old_sol;
    double maxq, vol, vol_0;
    double **a;

    h     = new double [L];
    y     = new double [L];
    q     = new double [L];
    r     = new double [L];
    r_fix = new double [L];

    noise = new double [L+2];
    ind   = new int    [L+4];     // Determines boundary conditions

    h2    = new double [L];       // To make calculations faster -> h2 = h^2, h3=h^3, ...
    h3    = new double [L];
    y2    = new double [L];
    y3    = new double [L];
    old_sol = new double [L];
    
    a = new double *[L];
    for(i=0;i<L;i++)
        a[i]=new double [5];

    no_good_solutions=0;
    j_time=t0;
    jdt=t0;
    jdt_no=t0+dt;
    dt_old=dt;

    create_boundary_conditions(ind, L); // ind - index for h, includes zero-der. boundary conditions

    call_initial(h, old_sol, L, dx, b, t0, dir); // creates initial conditions -> returns the value of the average film height (if needed)
    
    vol_0 = calculate_volume(L, h, dx);
    vol=vol_0;

    
// Simulation begins
do
{
ponovo:
    
for(i=0;i<L;i++)
	{   y[i]=h[i];
        h2[i]=h[i]*h[i]; h3[i]=h[i]*h[i]*h[i];
    }
    
    create_noise(L, h3, noise, sigma, dt, dx, inic);
    create_funcional_st(L, r_fix, h, h3, ind, theta, dx,dt);
    create_funcional_add_vdw_noise(L, r_fix, h, noise, ind, A, B, dx, dt);

//--------Implicit iteration process - Newton-Raphson method----------------

    no_newton_it=0;
    do {
        no_newton_it++;
        for(i=0;i<L;i++)
            {  y2[i]=y[i]*y[i];
               y3[i]=y[i]*y[i]*y[i];
            }
        
        create_funcional_st(L, r, y, y3, ind, theta-1, dx, dt);
        
        for(i=0;i<L;i++)
            {  r[i]=r_fix[i]-r[i];
               q[i]=0;
            }
        
        calculate_matrix_A(a,L, y, y2, y3, theta, dx, dt);
        elimination(a,r,q,L); // Gauss elimination to calculate q from Aq=r;
 
        switch (check_newton_convergence(L, h, y, q, dt, jdt, no_good_solutions, no_consecutive_neg, no_newton_it, maxq))
            {
                case 0: break;
                case 1: goto ponovo;
                        break;
                case 2: goto kraj;
                        break;
            }
        
        } while(maxq>tol);

    no_consecutive_neg=0;

//---End of Newton-Raphson iterative method-->
//----Checking errors (vol conservation, second time derivative small enough)
 
    check_volume(L, vol_0, h, y, dx, jdt);
    if (check_second_time_derivative(L, h, y, old_sol, dt, dt_old, no_good_solutions) == 1 ) goto ponovo;
    

//------------The solution is good enough - write in file-------------------------

    if(jdt>=jdt_no)
    {
 		cout << "Time step: " << jdt <<  ". Volume " << vol << " dt "<< dt << " max:  " << h[266] << endl;
        write_in_File(h, ind, L, dx, A, B, j_time, dir);
        if(jdt_no<1) jdt_no=jdt_no*10; else jdt_no=trunc(jdt+1);
 		j_time++;
     }

    jdt=jdt+dt;
    dt_old=dt;

    if (no_good_solutions>100 && dt<1e-3 )   //Everything works fine => the time step increases
        {   no_good_solutions=0;
            dt=dt+0.001*dt;
        }

    for(i=0;i<L;i++) {old_sol[i]=h[i]; h[i]=y[i];}   //---New solution -> initial condition

    
} while (jdt<=t);   //----------------- End of the simulation - writing data -----------------//
    
kraj: 

    delete       h;
    delete       y;
    delete       q;
    delete       r;
    delete   r_fix;
    delete   noise;
    delete     ind;
    delete      h2;
    delete      h3;
    delete      y2;
    delete      y3;
    delete old_sol;
    delete[]     a;
    


    return 0;
}

