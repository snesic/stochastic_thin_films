
#ifndef __READ_WRITE_MSG_H_INCLUDED__
#define __READ_WRITE_MSG_H_INCLUDED__

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;


int write_email_with_sim_data(int, int, double, double, double, double, double, string, string, int);

void print_initial_parameters(double, double, double, int, double, double, double, int, int, double, double);

int write_in_File(double *, int *, int, double, double , double ,int, string);

int read_position(double *, int, string);

#endif 
