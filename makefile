CXX = g++ 
CXXFLAGS =  -Wno-deprecated -w -O4 -Wall

DIFF = ./sdiff
PRE = ./
MAJOR = 1
MINOR = 0

%.o:           	%.cpp
		$(CXX) $(CXXFLAGS) -c $*.cpp

everything:    	instabilities

instabilities_obj = instabilities.o film.o initial_conditions.o newton-raphson_tools.o read_write_msg.o

instabilities:    $(instabilities_obj)
		$(CXX) -o $@ $(instabilities_obj) -L. -lm $(CXXFLAGS)


instabilities.o: instabilities.cpp read_write_msg.h film.h

film.o: film.cpp initial_conditions.h newton-raphson_tools.h read_write_msg.h

initial_conditions.o: initial_conditions.cpp read_write_msg.h

read_write_msg.o: read_write_msg.cpp 

newton-raphson_tools.o: newton-raphson_tools.cpp



instabilities.txx:     instabilities
		$(PRE) instabilities > instabilities.txx
		$(DIFF)instabilities.txt instabilities.txx


