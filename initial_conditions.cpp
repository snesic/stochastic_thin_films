

#include "initial_conditions.h"
#define PI 3.141592653589793

void initial_condition(double *h, int L, double dx, double b)   // Initial condition: Droplet as a  b + re^{(x-{x_mid})^4/r^4}
{
  int i; double height_drop, x_mid, x_pos; 

  x_mid=L/2*dx; 
  height_drop=1;  

  for(i=0;i<L;i++)
    {    x_pos=i*dx; 
	     h[i]=b + height_drop * exp(-pow((x_pos-x_mid)/height_drop,4));
    }
}

void initial_condition_sin(double *h, int L, double dx, double b)   // Initial condition: Droplet as a  b + re^{(x-{x_mid})^4/r^4}
{
    int i; double A, B;
    
    B=0.1; A=0.01;
    for(i=0;i<L;i++)
        h[i]= B + A*cos(2*PI*i*dx/2.66);
}

int call_initial( double *h, double *old_sol, int  L, double dx, double b, int t0, string dir)
{  string file_name;                        // from were the init cond is called
    stringstream sim_name;
    
    if(t0==0)
    { initial_condition_sin(h,L,dx,b); initial_condition_sin(old_sol,L,dx,b);}
    else {  sim_name.str(""); sim_name<<t0;          // Read initial position from a file!!!
        file_name= dir +"/droplet_at_" + sim_name.str()+ ".txt";
        if(read_position(h,L,file_name)==1)
        {cout<< "No initial condition file, simulation aborted!" <<endl;
            return 1;}
        read_position(old_sol,L,file_name);
    }
    
return 0;
}


void create_boundary_conditions(int *ind, int L)
{ int i;
    for(i=0;i<L;i++)
        ind[i+2]=i;        //  ind[i] - to use boundary conditions, shifted for 2
    ind[0]=2; ind[1]=1; ind[L+2]=L-2; ind[L+3]=L-3;
}
