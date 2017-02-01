#include "film.h"


film::film(int L, double dx, double dt, double b, int ugao, double sigma): L(L), dx(dx), dt(dt), b(b), angle(ugao), sigma(sigma)
{
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
    for(int i=0;i<L;i++)
        a[i]=new double [5];

    
    double a_min = 1.0;
    
    A     = 2*(1-cos(PI*ugao/180))*(a_min*b);
    B     = A*(a_min*b); // vdW: The min of vdW potential determins B..
    tol   = 1e-10;
    theta = 0.5;
    create_boundary_conditions(ind, L); // ind - index for h, includes zero-der. boundary conditions
    
    
}

void film::set_h2_3()
{
    for(int i=0;i<L;i++)
    {
        y[i]=h[i];
        h2[i]=h[i]*h[i];
        h3[i]=h[i]*h2[i];
    }
    
}

void film::set_y2_3()
{
    for(int i=0;i<L;i++)
    {
        y2[i]=y[i]*y[i];
        y3[i]=y[i]*y2[i];
    }
}


void film::create_noise(long *inic)
{
    double dht=sqrt(dt*sigma);
    static boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
    generator(boost::mt19937(time(0)),boost::normal_distribution<>());

    
    for(int i=0;i<L;i++)
        noise[i+1]=dht/dx * sqrt(h3[i]/dx) * generator();
    
    // Boundary conditions for the noise:
    noise[0]=0;
    noise[1]=0;
    noise[2]=0;
    noise[L+1]=0;
    noise[L]=0;
    noise[L-1]=0;
}



void film::create_funcional_st(double *r1, double *h1, double *h13, double theta1)
{
    double dx4dt=dt/(dx*dx*dx*dx);
    
    for(int i=0;i<L;i++)
    {
        r1[i]=h1[i] - dx4dt*theta1*(
                       0.5*( h13[i]+h13[ind[i+3]] ) * (h1[ind[i+4]]-3*h1[ind[i+3]]+3*h1[i]-h1[ind[i+1]])
                      -0.5*( h13[ind[i+1]]+h13[i] ) * (h1[ind[i+3]]-3*h1[i]+3*h1[ind[i+1]]-h1[ind[i]]  )
                       );
    }
}



void film::create_funcional_add_vdw_noise()
{
    double dx2dt=dt/(dx*dx);
    
    for(int i=0;i<L;i++)
        r_fix[i] = r_fix[i] + dx2dt*(               // Multiply this by theta if vdw treated implicitly
                             ( 6.0*B/(h[ind[i+3]]+h[i]) - 2.0*A ) * (h[ind[i+3]] - h[i])
                            -( 6.0*B/(h[ind[i+1]]+h[i]) - 2.0*A ) * (h[i] - h[ind[i+1]])
                             )
        + ( noise[i+2] - noise[i] )/2.0;
}




int film::simulation(double t, double& t0, long int *inic)
{   int i;
    int no_newton_it, no_consecutive_neg, no_good_solutions;
    double jdt,dt_old;
    double maxq, vol_0;
    no_good_solutions=0;
    jdt = t0;
    dt_old = dt;

    vol_0 = calculate_volume(L, h, dx);
    
    do
    {
        // Simulation begins
    ponovo:
        
        if(dt<1e-15) {cout<< "Wrong parameters" << endl; goto kraj;}
        
        set_h2_3();
        
        create_noise(inic);
        create_funcional_st(r_fix, h, h3,theta);
        create_funcional_add_vdw_noise();
        
        //--------Implicit iteration process - Newton-Raphson method----------------
        
        no_newton_it=0;
        do {
            no_newton_it++;
            set_y2_3();
            
            create_funcional_st(r, y, y3,theta-1);
            
            for(i=0;i<L;i++)
            {   r[i]=r_fix[i]-r[i];
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
        
        jdt=jdt+dt;
        dt_old=dt;
        
        if (no_good_solutions>100 && dt<1e-3 )   //Everything works fine => the time step increases
        {
            no_good_solutions=0;
            dt=dt+0.001*dt;
        }

        //swap(h,old_sol);
        //copy ( y, y+L, h );

        for(i=0;i<L;i++) { old_sol[i]=h[i]; h[i]=y[i]; }   //---New solution -> initial condition
        
        
    } while (jdt<=t);   //----------------- End of the simulation - writing data -----------------//
    
     t0=jdt;
    
kraj:
    return 0;
}





