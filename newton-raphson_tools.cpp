
#include "newton-raphson_tools.h"



void elimination(double ** a, double *r, double *q, const int& L)    // Gaussan elimination of 5 diagonal matrix, gives q out of Aq=r
{
  int i,j,l; double c;

  for(l=0;l<L-2;l++)
    for(j=0;j<2;j++)
     {
   	c=a[j+1+l][1-j]/a[l][2]; 
	a[j+1+l][1-j]=0;
   	for(i=1;i<4;i++)
      	    a[j+1+l][i-j]=a[j+1+l][i-j]-c*a[l][1+i];
    	    r[j+1+l]=r[j+1+l]-r[l]*c;
      }

  c=a[L-1][1] / a[L-2][2];
  a[L-1][1]=0;
  a[L-1][2]=a[L-1][2]-a[L-2][3]*c;

  r[L-1]=r[L-1]-r[L-2]*c;

  // calculating solution reverse 

   q[L-1]=r[L-1]/a[L-1][2];
   q[L-2]=(r[L-2] - q[L-1]*a[L-2][3])/a[L-2][2];

   for(i=L-3;i>=0;i--)
       q[i]=(r[i] - q[i+1]*a[i][3] - q[i+2]*a[i][4])/a[i][2];

}

void calculate_matrix_A(double **a, const int& L, double *y, double *y2, double *y3, const double& theta, const double& dx, const double& dt)
{ int i,j; double dx4dt = dt/(dx*dx*dx*dx);
    
    // B.C. for the coefficients left:
    a[0][2] =   -3*( y3[0] + y3[1] ) -   3*y2[0]*( y[2]-4*y[1]+3*y[0] );
    a[0][3] =    2*( y3[0] + y3[1] ) - 1.5*y2[1]*( y[2]-4*y[1]+3*y[0] );
    a[0][4] = -0.5*( y3[0] + y3[1] );
    
    a[1][1] =  0.5*( y3[1] + y3[2] ) + 1.5*( y3[0] + y3[1] ) + 1.5*y2[0]*( y[2]-4*y[1]+3*y[0] );
    a[1][2] = -1.5*y2[1]*( y[3]-4*y[2]+7*y[1]-4*y[0] )
    -1.5*( y3[0] + 2*y3[1] + y3[2]);
    a[1][3] = -1.5*y2[2]*( y[3]-3*y[2]+3*y[1]-y[0] )
    + 1.5*( y3[1] + y3[2])
    + 0.5*( y3[0] + y3[1] );
    a[1][4] = -0.5*( y3[1] + y3[2] );
    
    // B. C. for the newton_Raphson equation:
    
    a[0][3]=2*a[0][3];
    a[1][2]=a[1][2] + a[0][4];//-0.5*(y3[0]+y3[1]);
    a[0][4]=2*a[0][4];
    
    for(i=2;i<L-2;i++)
    {
        a[i][0] = -0.5*( y3[i] + y3[i-1] );
        a[i][1] =  0.5*( y3[i] + y3[i+1] ) + 1.5*y2[i-1]*( y[i+1]-3*y[i]+3*y[i-1]-y[i-2] )
        + 1.5*( y3[i] + y3[i-1] );
        a[i][2] = -1.5*( y3[i] + y3[i+1] ) - 1.5*y2[i] * ( y[i+2]-4*y[i+1]+6*y[i]-4*y[i-1]+y[i-2] )
        -1.5*( y3[i] + y3[i-1] );
        a[i][3] =  1.5*( y3[i] + y3[i+1] ) - 1.5*y2[i+1]*( y[i+2]-3*y[i+1]+3*y[i]-y[i-1] )
        + 0.5*( y3[i] + y3[i-1] );
        a[i][4] = -0.5*( y3[i] + y3[i+1] );
    }
    
    // B.C. for the coefficients right:
    a[L-2][0] = -0.5*( y3[L-3] + y3[L-2] );
    a[L-2][1] =  0.5*( y3[L-2] + y3[L-1] ) + 1.5*y2[L-3]*( y[L-1]-3*y[L-2]+3*y[L-3]-y[L-4] )
    +1.5*( y3[L-3] + y3[L-2] );
    a[L-2][2] = -1.5*( y3[L-2] + y3[L-1] ) - 1.5*y2[L-2] * ( -4*y[L-1]+7*y[L-2]-4*y[L-3]+y[L-4] )
    -1.5*( y3[L-2] + y3[L-3] );
    a[L-2][3] =  1.5*( y3[L-2] + y3[L-1] ) - 1.5*y2[L-1]*( -3*y[L-1]+4*y[L-2]-y[L-3] )
    +0.5*( y3[L-2] + y3[L-3] );
    
    a[L-1][0] = -0.5*( y3[L-1] + y3[L-2] );
    a[L-1][1] =    2*( y3[L-1] + y3[L-2] ) + 1.5*y2[L-2]*( 4*y[L-2]-3*y[L-1]-y[L-3] );
    a[L-1][2] = -  3*( y3[L-1] + y3[L-2] ) -   3*y2[L-1]*( y[L-3]+3*y[L-1]-4*y[L-2] );
    
    // B. C. for the newton_Raphson equation:
    
    a[L-1][1]=2*a[L-1][1];
    a[L-2][2]=a[L-2][2] + a[L-1][0]; //-0.5*(y3[L-1]+y3[L-2]);
    a[L-1][0]=2*a[L-1][0];
    
    for(i=0;i<L;i++)
    for(j=0;j<5;j++)
    a[i][j]=theta*a[i][j]*dx4dt*(-1);     // This should be changed (should be done before)
    for(i=0;i<L;i++) a[i][2]=1+a[i][2];
}

double calculate_volume(const int& L, double *h, const double& dx)
{ int i; double vol=0;
    
for(i=0;i<L;i++)
    vol = vol + h[i];
    
return vol*dx;
}

void check_volume(const int& L, double& vol_0, double *h, double *y, const double& dx, const double& jdt)
{
int i, maxq;
double vol;
    
    vol=calculate_volume(L,y,dx);

    if(fabs(vol-vol_0)/vol_0>0.03)
    {   cout<< "Problem!!!! No Volume conservation at " << jdt << "!! Initial vol:  " << vol_0 << "  and now:  " << vol << endl; vol_0=vol;
        vol=0; maxq=0;
        for(i=0;i<L;i++)
            if(fabs(h[i]-y[i])>vol)
                {vol=h[i]-y[i]; maxq=i;}
        cout<< vol << "  Biggest change at the position: " << maxq << "  and time:  "<< jdt <<  endl;
        vol=fabs(vol);
    }
}

int check_second_time_derivative(const int& L, double *h, double *y, double *old_sol, double& dtt, double dt_old, int& no_good_solutions)
{ double maxq, t_err, dt=dtt;
    int i;

    maxq = (2*dt/dt_old)*(y[0]*dt_old+old_sol[0]*dt-(dt_old+dt)*h[0])/(h[0]*(dt+dt_old));

    for(i=0;i<L;i++)
    {   t_err = (2*dt/dt_old)*(y[i]*dt_old+old_sol[i]*dt-(dt_old+dt)*h[i])/(h[i]*(dt+dt_old));
        if(t_err>maxq) maxq=t_err;
    }
    
    if (maxq<0.01)
    {   no_good_solutions++;
        return 0;
    }
    else{ if(maxq>0.1)
            {   no_good_solutions=0;
                dtt=dt-pow(10,-fabs(round(log10(dt)))-2);
                //dtt=floorf(dt * pow(10,fabs(round(log10(dt)))+1) *100 ) / (100*pow(10,fabs(round(log10(dt)))+1));
                cout << "Time derivative error, repeat iteration for smaller dt " << dtt << endl;
                return 1;
            }
        }
    
return 0;
}

int check_newton_convergence(const int& L, double *h, double *y, double *q, double& dt, const double& jdt, int& no_good_solutions, int& no_consecutive_neg, const int& no_newton_it, double& mq)
{   double maxq=0;
    int i;
    
    for(i=0;i<L;i++)
    {   if(maxq<fabs(q[i])) maxq=fabs(q[i]);
        
        if(y[i]+q[i]<0)
            {
            cout<< h[i] << "  at " << jdt << " , on the position " <<  i   << "  is negative!" << endl;
            no_consecutive_neg++;
            if(no_consecutive_neg>10) {cout << "Wrong parameters -> Simulation breaks!" << endl; return 2;}
            return 1;
            }
        else
        y[i]=y[i]+q[i];
    }
    
    if(no_newton_it>10)
        {   no_good_solutions=0;
            dt=dt-pow(10,-fabs(round(log10(dt)))-2);
            return 1;
        }
    
    mq=maxq;
    
return 0;
}

