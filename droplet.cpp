#include"droplet.h"
#define PI 3.141592653589793

using namespace std;

int droplet::read_surface(string dir_n, int iter, int t)
{   string get_word;
    int i, *ind;
    string file_name;
    stringstream iteration_str,time_step;
    iteration=iter; time=t;
    dir_name=dir_n;
    i=0;
   
    time_step.str("");     iteration_str.str(""); 
    time_step<<time; iteration_str<< iteration; 
    file_name= dir_name + "/" + iteration_str.str() + "/droplet_at_" + time_step.str() + ".txt";
    ifstream File(file_name.c_str());

    if(File.is_open())
        while(!File.eof())
        {
            getline(File, get_word);
            i++;
        }
    else
    {
        cout << endl << "Failed to open file " << file_name; cout<< endl;
        return 1;
	}

    L=i-1;

    File.clear();
    File.seekg(0, ios::beg);

    h       = new double [L];
    h_prim  = new double [L];
    h_sec   = new double [L];
    h_thr   = new double [L];
    x       = new double [L];
    Sq      = new double [L];
    ind     = new int    [L+4];     // Determines boundary conditions

    i=0;
    if(File.is_open())
        while(!File.eof())
        {
            File>> x[i];
            File>> h[i];
            i++;
        }

    dx= x[1]-x[0];

//----------Calculating derivatives

    for(i=0;i<L;i++)
        ind[i+2]=i;        // ind[i] - to use boundary conditions, shifted for 2
    ind[0]=2; ind[1]=1; ind[L+2]=L-2; ind[L+3]=L-3;

    x_max=x[0]; h_max=h[0];   // Global maximum!
    for(i=1;i<L;i++)
    {	if(h[i]>h_max) {h_max=h[i]; x_max=x[i];}
        h_prim[i] = (h[ind[i+3]]-h[ind[i+1]])/(2*dx);
        h_sec[i]  = (h[ind[i+3]]+h[ind[i+1]]-2*h[i])/(dx*dx);
    	h_thr[i]  = (h[ind[i+4]] - 2*h[ind[i+3]] + 2*h[ind[i+1]] - h[ind[i]])/(8*dx*dx*dx);
    }
    if(h[L-1]>h_max) {h_max=h[L-1]; x_max=x[L-1];}

    droplets_number();
    drop_distances();
    delete     ind;
    
return 0;
}

//---------------------------------Volume-----------------------------------------------

double droplet::volume()
{int i; double vol;

    vol=0;
    for(i=0;i<L;i++)
        vol=vol+h[i];

return (vol*dx);
}

//-------------------------Number of drops----------------------------

int droplet::droplets_number()
{   int i, j, no_par_der;  int drop_no;
 
    drop_no=0;
    j=1;
    
    if(h[0]>mean_height() )
	{
        drop_no++;
        j++;
    }                                // Droplet at the left boundary

    for(i=20;i<L-30;i++)
        if(h_prim[i-1]>0 && h_prim[i]>=0 && h_prim[i+1]<0 && h[i]>mean_height())
            {   //Droples in the domain located using derivatives
                drop_no++;
                drop_no++;
                j++;
                i = i+40;
            }
 
    if(h[L-1]>mean_height())
        drop_no++;
	                           // Droplet at the left boundary
return drop_no; 
}


//-------------------------Distances between drops----------------------------

int droplet::drop_distances()
{   int i, j, no_par_der, drop_no;
 
    x_poz_max = new double[droplets_number()];
    h_poz_max = new double[droplets_number()];
    l_j       = new double[droplets_number()];
    
    for(i=0;i<droplets_number();i++)
    {
        x_poz_max[i]=0;
        h_poz_max[i]=0;
        l_j[i]=0;
    }

    j=0;
    if(h[0]>mean_height() )
		{
            x_poz_max[j] = x[0];
		 	h_poz_max[j] = h[0];
			j++;
		}                                // Droplet at the left boundary

    for(i=20;i<L-30;i++)
        if(h_prim[i-1]>0 && h_prim[i]>=0 && h_prim[i+1]<0 && h[i]>mean_height())
        {       x_poz_max[j] = x[i];
                h_poz_max[j] = h[i];        //Droples in the domain located using derivatives
                j++;
                i = i + 40;
        }

    if(h[L-1]>mean_height())
    {       x_poz_max[j] = x[L-1];
            h_poz_max[j] = h[L-1];
    }                                // Droplet at the left boundary

return 0; 
}


double droplet::calculate_l_i(int j_droplet, int i_droplet)
{
    double x_poz, vol, l_t;
    int i,j;

    vol=0;

    if (i_droplet==0) // The droplet is not just a half
    {
        x_poz=0;
        i=0;
        do
        {
     		vol=vol + (h[i]-b);
     		i++;
        }while (h[i]>0.02);         // Precursor is set to 0.01, this has to be changed if precursor changed
        
        vol=vol*dx;
        l_t=0;

        do
        {
            l_t=l_t + dx*pow(i*dx,2)*(h[i]-b);
            i++;
        }while (h[i]>0.02);
        
        l_t=l_t/vol;
        return sqrt(2*l_t);
        
    }else if (i_droplet==L-1) // The droplet is not just a half
        {
            x_poz=(L-1)*dx;
            i=L-1;
            vol=0;
            
            do
            {
                vol=vol + (h[i]-b);
                i--;
            }while (h[i]>0.02);
            
            vol=vol*dx;
            l_t=0;
            i=L-1;
            
            do
            {
                l_t=l_t + dx*pow(x_poz-i*dx,2)*(h[i]-b);
                i--;
            }while (h[i]>0.02);
            
            l_t=l_t/vol;
            
            return sqrt(2*l_t);
            
        } else
        {    // If the droplet is not just a half
            vol=0;
            i=1;
            vol=vol + (h[i_droplet] - b);
            do
            {
                vol=vol + ( h[i_droplet+i] + h[i_droplet-i] - 2*b);
                i++;
            }while (h[i_droplet+i]>0.02);
  	
            vol=vol*dx;
            i=1;
            x_poz=i_droplet*h[i_droplet];

            do
            {
                x_poz=x_poz+(i_droplet+i)*h[i_droplet+i] + (i_droplet-i)*h[i_droplet-i];
                i++;
            }while (h[i_droplet+i]>0.02);
         
            x_poz=dx*dx*x_poz/vol;

            l_t=dx*pow(((i_droplet)*dx-x_poz),2)*(h[i_droplet]-b);
            i=1;
            do
            {
                l_t=l_t + dx*pow(((i_droplet+i)*dx-x_poz),2)*(h[i_droplet+i]-b);
                l_t=l_t + dx*pow(((i_droplet-i)*dx-x_poz),2)*(h[i_droplet-i]-b);
                i++;
            }while (h[i+i_droplet]>0.02);
            l_t=l_t/vol;
            return sqrt(l_t);
        }
return 0;
}


int droplet::tanner()
{
  double x_poz;
  int i,j;

/*
i=0;
do       
j=int (x_poz_max[i]/dx);

          do 
                x_poz=x_poz+i*h[i];
          while (h[j]>0.02)
            x_poz=dx*dx*x_poz/(this->volume());

//    cout << x_poz_max[i+1]-x_poz_max[i] << "  " << i <<  endl;
/*
       l[i_f]=0;
            for(i=0;i<L;i++)
                if(h[i]>2*b)
					l_rel[i_f]=l_rel[i_f] + dx*pow(i*dx-x_poz,2)*(h[i]-b);
				l_rel[i_f]=l_rel[i_f]/volume;
        }
    
	for(i=0;i<br_fajlova_end-br_fajlova_0;i++)
    	     l=l+sqrt(l_rel[i])/(br_fajlova_end-br_fajlova_0);

} while(x_poz_max[i]!=0);
*/
return 0;
}

double droplet::mean_height()
{   int i; double mean;
    mean=0;
    for(i=1;i<L;i++)
        mean = mean+h[i];
return mean/L;
}

double droplet::min_height()
{   int i; double min;
    min=h[0];
    for(i=1;i<L;i++)
        if(min>h[i]) min=h[i];
return min; 
}

double droplet::max_height()
{   int i; double max;
    max=h[0];
    for(i=1;i<L;i++)
        if(max<h[i]) max=h[i];
return max; 
}


int droplet::fourier_trasform()
{
	int i;
	fftw_plan p;
	fftw_complex *izlaz;  
   double *h_p,m;

    h_p= new double [L];

    for(i=0;i<L;i++) Sq[i]=0;
 
    m=mean_height(); 
    for(i=0;i<L;i++) h_p[i]=h[i]-m;

//---------Fourier transform----------------------//
		
    izlaz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L);  // def complex sequence
    p= fftw_plan_dft_r2c_1d(L,h_p,izlaz,FFTW_ESTIMATE);  // power basis
    fftw_execute(p);
		
    Sq[0] = Sq[0] + ( izlaz[0][0]*izlaz[0][0]  + izlaz[0][1]*izlaz[0][1] );
    for(i=1;i<L/2;i++)
        Sq[i] = Sq[i] + ( izlaz[i][0]*izlaz[i][0] + izlaz[i][1]*izlaz[i][1] )/L;
    
    fftw_destroy_plan(p);
    free(izlaz);

//	    upisi_Sq(); 
delete h_p;
}



int droplet::fourier_trasform_cosine()
{   int i;
    fftw_plan p;
    fftw_complex *izlaz;
    double *h_p,m, *izlaz_cos;

    h_p		= new double [L]; 
    izlaz_cos= new double [L];

    for(i=0;i<L;i++) Sq[i]=0;
 
    m=mean_height(); 
    for(i=0;i<L;i++)
        h_p[i]=h[i]-m;

//---------Fourier transform----------------------//

    p= fftw_plan_r2r_1d(L,h_p,izlaz_cos,FFTW_REDFT10,FFTW_ESTIMATE);  // Cosine basis
    fftw_execute(p);
    Sq[0] = Sq[0] + izlaz_cos[0]*izlaz_cos[0];
    
    for(i=1;i<L;i++)
        Sq[i] = Sq[i] +  izlaz_cos[i]*izlaz_cos[i]/L; // Cosine Structure Factor
	
    fftw_destroy_plan(p);
    delete izlaz_cos;
	//    upisi_Sq();
return 0;
}


int droplet::upisi_Sq()
{
    ofstream dat;
    int i;
    string file_name;
    stringstream iteration_str,time_step;
    
    time_step.str("");     iteration_str.str("");
    time_step<<time; iteration_str<< iteration;
    file_name= dir_name + "/" + iteration_str.str() + "/droplet_sq_at_" + time_step.str() + ".txt";
    
    dat.open(file_name.c_str(), ios::out);
	    for(i=1;i<L/2;i++)
	        dat << 2*PI*i/L << " " << Sq[i] << endl;
    dat.close();

return 0;
}

void droplet::Delete_droplet()
{
delete this;
}


void histogram(double *niz, long int l, double max, double min, int korak, string dir_name, int time, int broj )
{
	double duzina=(max-min)/korak;
	long int i,p;
	ofstream dat;		
	double *histo_niz, norm;
	string file_name;
	norm=0;
	histo_niz=new double[korak];

    stringstream time_step;

    time_step.str("");
    time_step<<time;
 
		
    for (i=0;i<korak;i++) histo_niz[i]=0;

	for(i=0;i<l;i++)
		{
            p= (niz[i]-min)/duzina;
            histo_niz[p]++;
		}
			
//-----------------Norm---------------------------------//

	for(i=0;i<korak;i++)
		norm=norm+histo_niz[i]*duzina;
			
    if(broj==1)	file_name= dir_name + "/histogram_l_" + time_step.str() + ".txt";
    if(broj==2)	file_name= dir_name + "/histogram_x_" + time_step.str() + ".txt";
    if(broj==3)	file_name= dir_name + "/histogram_h_" + time_step.str() + ".txt";
	
    dat.open(file_name.c_str(),  ios::out);
	for(i=0;i<korak;i++)
		if(histo_niz[i]>0)
			dat << min + duzina/2 + i * duzina  << "\t" << histo_niz[i]/norm<<endl;	
	dat.close();
    delete histo_niz;
    
}




