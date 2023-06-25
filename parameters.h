//
// This program is introduced in the paper:
// Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians (2023).
//
// This program is licensed under a Creative Commons Attribution 4.0 International License:
// http://creativecommons.org/licenses/by/4.0/
//

//
// Below are the parameter values:
//

#define Tsteps 1000000ull            // number of Monte-Carlo initial equilibration updates (use "ull" suffix for unsigned long long)
#define steps  10000000ull           // number of Monte-Carlo updates (use "ull" suffix for unsigned long long)
#define stepsPerMeasurement 100ull   // number of Monte-Carlo updates per measurement (use "ull" suffix for unsigned long long)
#define qmax     300                 // upper bound for the maximal length of the sequence of permutation operators
#define Nbins    250                 // number of bins for the error estimation via binning analysis
#define exhaustive_cycle_search      // comment this line for a more restrictive cycle search

//
// Below is the list of standard observables:
//

#define MEASURE_H                    // <H>             is measured when this line is not commented
#define MEASURE_H2                   // <H^2>           is measured when this line is not commented
#define MEASURE_Hdiag                // <H_{diag}>      is measured when this line is not commented
#define MEASURE_Hdiag2               // <H_{diag}^2>    is measured when this line is not commented
#define MEASURE_Hoffdiag             // <H_{offdiag}>   is measured when this line is not commented
#define MEASURE_Hoffdiag2            // <H_{offdiag}^2> is measured when this line is not commented
