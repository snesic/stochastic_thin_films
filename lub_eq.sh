#!/bin/bash

g++ -c random_no_generators.cpp -Wno-deprecated -w -O4
g++ -c read_write_msg.cpp -Wno-deprecated -w -O4
g++ -c newton-raphson_tools.cpp -Wno-deprecated -w -O4
g++ -c initial_conditions.cpp -Wno-deprecated -w -O4
g++ -c create_forces.cpp -Wno-deprecated -w -O4
g++ -c simulation.cpp -Wno-deprecated -w -O4
g++ instabilities.cpp simulation.o create_forces.o newton-raphson_tools.o initial_conditions.o random_no_generators.o read_write_msg.o -o run_thin_film -Wno-deprecated -w -O4
rc=$?
echo $rc
if [[ $rc -eq 0 ]]; then
    echo Program starts:
    if [ ! -d "test_dir/0" ]; then
        mkdir "test_dir"
        mkdir "test_dir/0"
    fi

    ./run_thin_film 0 test_dir
exit $?
fi
exit $rc



#./proba
exit
