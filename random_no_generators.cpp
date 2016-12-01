
#include "random_no_generators.h"

#define PI 3.141592653589793
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

using namespace std;


double ran3(long *idum)
{
static int inext,inextp;
static long ma[56]; 
static int iff=0;
long mj,mk;
int i,ii,k;


	if (*idum < 0 || iff == 0) { 

			iff=1;
			mj=labs(MSEED-labs(*idum)); 
			mj %= MBIG; 
			ma[55]=mj;
			mk=1;
			for (i=1;i<=54;i++) { 
					ii=(21*i) % 55; 
					ma[ii]=mk; 
					mk=mj-mk;
					if (mk < MZ) mk += MBIG;
				        mj=ma[ii];
			}

			for (k=1;k<=4;k++) 
				for(i=1;i<=55;i++) {
							ma[i] -= ma[1+(i+30) % 55];
							if (ma[i] < MZ) ma[i] += MBIG;
				    		}
	inext=0; 
	inextp=31; 
	*idum=1;
	}

	if (++inext == 56) inext=1; 
	if (++inextp == 56) inextp=1; 
	mj=ma[inext]-ma[inextp]; 
	if (mj < MZ) mj += MBIG; 
	ma[inext]=mj;
 
return mj*FAC;
}

double gauss(double m, double s, long *idum)			/* normal random variate generator */
{				     				   /* mean m, standard deviation s */
	
double ran3(long *idum);
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * ran3(idum) - 1.0;
			x2 = 2.0 * ran3(idum) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}
	if(!isfinite(m + y1 * s))  cout<< "Generator!!"<<endl;
	return( m + y1 * s );
}
