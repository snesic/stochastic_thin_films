# stochastic_thin_films
2d Stochastic Thin Films - Implicit algorithm to solve 1d lubrication equation with thermal noise  

Use make to compile, or lub_eq.sh to compile and run. 

To run in terminal: [./run_thin_film sim_no dir_name]. The program only runs if there is a directory dir_name and subdirectory sim_no (sim_no must be integer), for example [test_dir/0].  



Files: 

lub_eq.sh   -> compile and run the program on UNIX systems 

makefile -> compile by make on UNIX systems


instabilities.cpp ---> main file, sets the parameters and calls simulation 

simulation.h   ----> implicit stochastic thin film solver. Solves the equation and writes down the solutions (thin film surface) at integer times from t0 to t.  

create_forces.h  ---> Creates forces (surface tension, vdw and thermal noise) for the simulation 

random_no_generators.h -> Gaussian random number generator

initial_conditions.h -> Initial conditions (if t0=0 [Flat film, randomly perturbed flat film, droplet], if 
                                            if t0>0 reads from file droplet_surface_t0.txt)

newton-raphson_tools.h -> routines for Newton-Raphson method

read_write_msg.h -> input output messages



Analysis of the data: 
droplet.h 
