-----------------------------------------------------------------------------------------------------------

This program is introduced in the paper: Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians (2023).

-----------------------------------------------------------------------------------------------------------

Instructions for the Permutation Matrix Representation Quantum Monte Carlo for spin-1/2 Hamiltonians:

1. Prepare the Hamiltonian input text file "H.txt".
   Each line contains "J q_1 sigma_1 q_2 sigma_2 ...", where J is a constant, q_i is a spin index, and sigma_i = X, Y, and Z correspond to the Pauli matrices. It is also possible to use 1, 2, and 3 instead of X, Y and Z.

2. Check the parameter values of the simulation in the header file "parameters.hpp" such as the number of Monte-Carlo updates and the inverse temperature.

3. Check the list of standard observables in the header file "parameters.hpp"

4. Optionally, prepare the input text files for additional custom observables in the same format as for the Hamiltonian.
   For example, if one uses the same text file for the Hamiltonian and for an observable, then the mean energy will be measured.
   For a custom observable to be computed correctly, it is necessary to account for the restrictions on its structure - see the paper for details.
   File names for observables match the pattern "O*.txt", i.e., they begin with the letter "O" and have the extension "txt".

5. Compile and run "prepare.cpp". This utility reads the input files and prepares the "hamiltonian.hpp" header file:

		g++ -O3 -std=c++11 -o prepare.bin prepare.cpp
		./prepare.bin H.txt $(ls O*.txt  2> /dev/null)

6. Compile and run PMRQMC:

		g++ -O3 -std=c++11 -o PMRQMC.bin PMRQMC.cpp
		./PMRQMC.bin

-----------------------------------------------------------------------------------------------------------

Instructions for parallel computing with MPI (optional):

1. Perform the steps 1 through 5 above.
   Note that the parameters "Tsteps" and "steps" are numbers of Monte Carlo updates performed by each MPI process rather than total numbers of Monte Carlo updates.

2. Compile "PMRQMC_mpi.cpp":

		mpicxx -O3 -o PMRQMC_mpi.bin PMRQMC_mpi.cpp

3. Run the MPI application "PMRQMC_mpi.bin" using mpirun or any MPI-compatible job scheduler.

-----------------------------------------------------------------------------------------------------------
