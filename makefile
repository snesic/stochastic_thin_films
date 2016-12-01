CXX = g++
CXXFLAGS = -Wno-deprecated -w -O4 -Wall

DIFF = ./sdiff
PRE = ./
MAJOR = 1
MINOR = 0

%.o:           	%.cpp
		$(CXX) $(CXXFLAGS) -c $*.cpp

everything:    	instabilities

instabilities_obj = instabilities.o simulation.o create_forces.o newton-raphson_tools.o initial_conditions.o random_no_generators.o read_write_msg.o 
instabilities:    $(instabilities_obj)
		$(CXX) -o $@ $(instabilities_obj) -L. -lm $(CXXFLAGS)



instabilities.o: instabilities.cpp simulation.h read_write_msg.h

simulation.o:    simulation.cpp initial_conditions.h newton-raphson_tools.h create_forces.h read_write_msg.h

create_forces.o: create_forces.cpp random_no_generators.cpp 

initial_conditions.o: initial_conditions.cpp random_no_generators.h read_write_msg.h

random_no_generators.o: random_no_generators.cpp 

read_write_msg.o: read_write_msg.cpp 


instabilities.txx:   	instabilities
		$(PRE) instabilities > instabilities.txx
		$(DIFF)instabilities.txt instabilities.txx


