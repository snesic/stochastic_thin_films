
#include "create_forces.h"



void create_noise(int L, double *h3, double *noise, double sigma, double dt, double dx, long *inic)
{ int i;
  double dht=sqrt(dt*sigma);
    
    for(i=0;i<L;i++)
        noise[i+1]=dht/dx * sqrt(h3[i]/dx) * gauss(0,1,inic);

    // Boundary conditions for the noise:
    noise[0]=0;
    noise[1]=0;
    noise[2]=0;
    noise[L+1]=0;
    noise[L]=0;
    noise[L-1]=0;
}



void create_funcional_st(int L, double *r, double *h, double *h3, int *ind, double theta, double dx, double dt)
{ int i; double dx4dt=dt/(dx*dx*dx*dx);
    
for(i=0;i<L;i++)
    r[i]=h[i]           // Functional to be minimized, the part that doesn't change during the iterations
            - dx4dt*theta*(
                 0.5*( h3[i]+h3[ind[i+3]] ) * (h[ind[i+4]]-3*h[ind[i+3]]+3*h[i]-h[ind[i+1]])
               - 0.5*( h3[ind[i+1]]+h3[i] ) * (h[ind[i+3]]-3*h[i]+3*h[ind[i+1]]-h[ind[i]]  )
                           );
}




void create_funcional_add_vdw_noise(int L, double *r, double *h, double *noise, int *ind, double A, double B, double dx, double dt)
{ int i; double dx2dt=dt/(dx*dx);
    
    for(i=0;i<L;i++)
        r[i] = r[i] + dx2dt*(               // Multiply this by theta if vdw treated implicitly
                              ( 6*B/(h[ind[i+3]]+h[i]) - 2*A ) * (h[ind[i+3]] - h[i])
                             -( 6*B/(h[ind[i+1]]+h[i]) - 2*A ) * (h[i] - h[ind[i+1]])
                            )
                    + ( noise[i+2] - noise[i] )/2;
    
}


