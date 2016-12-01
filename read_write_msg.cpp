
#include "read_write_msg.h"


int write_email_with_sim_data(int L,int t,double dx,double dt,double b,double A,double sigma, string dir, string calc_time, int i_sim_0)
{
	string pera; char ch;
	int i,j;
    
    ifstream text("text_email.txt");
    ofstream text_4_mail;
    
	// Make sure the file stream is good
	if(!text) {
		cout << endl << "Failed to open the email file ";
		return 1;
    }
    
    text_4_mail.open("text_email_f.txt", ios::out);
    
    while(!text.eof())
    {
       ch=text.get();
       if(ch==' ' || ch=='\n')
        {  do { text_4_mail<< ch; ch=text.get(); } while (ch==' ' || ch=='\n');
        } text.seekg (-1, ifstream::cur);
        text >> pera;
        text_4_mail<< pera;
        if (pera=="number") text_4_mail<< " " << i_sim_0;
        if (pera=="after") text_4_mail<< " " << calc_time;
        if (pera=="are:") { text_4_mail<< "\n" << "Grid length L = "<<L<<endl;
                            text_4_mail<< "Total time = " << t << endl;
                            text_4_mail<< "Grid distance = " << dx << endl;
                            text_4_mail<< "Time step = " << dt << endl;
                            text_4_mail<< "Precursor = " << b << endl;
                            text_4_mail<< "Hamaker constant = " << A << endl;
                            text_4_mail<< "Noise strenght = " << sigma << endl;
                            text_4_mail<< "The results are written in = " << dir;
                            }
    }
        
      return 0;
}

void print_initial_parameters(double sigma, double dt, double dx, int L, double b, double A, double a_min, int t, int t0, double ugao, double i_sim_0)
{
    cout<< "Simulation of a stochastic lubrication equation with noise!  sigma = " << sigma << endl;
    cout<< "with the folowing parameters: " << " dt = " << dt << "  dx= " << dx << "   L = " << L << " b = " << b << "  A = " << A << endl;
    cout<< "where the vcW force has a minimum at h = " << a_min << " * b " << endl;
    cout<< "Total time: " << t << "    and total length: "<< L*dx << endl;
    cout<< "and the initial time: " << t0 << " , contact angle    "<< ugao << endl;
    cout<< "Starting iteration: "<< i_sim_0 << endl;
}

int write_in_File(double *h, int *ind, int L, double dx, double A, double B, int j_time, string dir) // write data into files.
{
    int i;
    ofstream dat;
    string file_name;                        // Were results will be written
    stringstream sim_name;
    
    
    sim_name.str("");
    sim_name << j_time;
    
    file_name= dir + "/droplet_at_" + sim_name.str()+ ".txt";
    dat.open(file_name.c_str(), ios::out);
	dat.precision(10);
	dat<< scientific;
    for(i=0;i<L;i++)
        dat << scientific << i*dx << "  " << scientific << h[i] << endl;
    dat.close();

return 0;
}

int read_position(double *poz, int L, string filename)
{
	std::ifstream inFile(filename.c_str());
	int i,j; double vst; 
	
	
	// Make sure the file stream is good
	if(!inFile) {
		cout << endl << "Failed to open file " << filename << endl;
		return 1;
	}
	i=0; j=0;
	while(!inFile.eof() && i<L) 
	{  inFile >> vst;
		if(j==1) poz[i]=vst;
	   j++; 
	   if (j==2) {j=0; i++;}	
		
	}
    
return 0;
}

