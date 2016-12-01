#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <fftw3.h>

using namespace std;

class droplet{

private: 
double *x, *h_prim, *h_sec, *h_thr;
double x_max, h_max, *l; 
int iteration, time;
double dx, b;
string dir_name; 

public:
        droplet() {};
	~droplet() { 	delete Sq;
			delete x_poz_max;
			delete h_poz_max;
			delete l_j;
			delete x;
			delete h;
			delete h_prim;
			delete h_sec;
			delete h_thr; }

    double *Sq, *x_poz_max, *h_poz_max, *l_j, *h;
	int L;
    int    read_surface(string, int , int);
	double volume();
    int    droplets_number();
	int    drop_distances();
	int    tanner();
	double calculate_l_i(int, int);
	double mean_height();
	double min_height();
	double max_height();
	int    fourier_trasform();
	int    fourier_trasform_cosine();
    void   Delete_droplet();

protected:
    int    upisi_Sq();

};

void histogram(double *, long int , double, double, int, string, int, int);
