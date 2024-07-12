Compilation
Compile the program using the following command:
mpic++ w.cpp -o e -fopenmp
Explanation of the command:
•   mpic++: The MPI C++ compiler.
•   -fopenmp: Compiler flag to enable OpenMP.
•   w.cpp: The source file.
•   -o e: Specifies the output executable file name
Execution
If run on single device
Execution
For single device
mpiexec -n 10 ./e
For multiple device
mpiexec  -n 10 -f machinefile ./e
