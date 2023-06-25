This program is introduced in the paper: Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians (2023).

Instructions for the Permutation Matrix Representation Quantum Monte Carlo for spin-1/2 Hamiltonians.

1. Prepare the Hamiltonian input text file "H.txt".
   Each line contains "J q_1 sigma_1 q_2 sigma_2 ...", where J is a constant, q_i is a spin index, and sigma_i = X, Y, and Z correspond to the Pauli matrices. It is also possible to use 1, 2 and 3 instead of X, Y and Z.

2. Check the parameter values of the simulation in the header file "parameters.h" such as the number of Monte-Carlo updates.

3. Check the list of standard observables in the header file "parameters.h"

4. Optionally, prepare the input text files for additional custom observables in the same format as for the Hamiltonian.
   For example, if one uses the same text file for the Hamiltonian and for an observable, then the mean energy will be measured.
   For a custom observable to be computed correctly, it is necessary to account for the restrictions on its structure - see the paper for details.
   File names for observables match the pattern "O*.txt", i.e., they begin with the letter "O" and have the extension "txt".

5. Run the bash script "run.sh".
