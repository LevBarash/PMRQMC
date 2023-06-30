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

#define Tsteps 1000000               // number of Monte-Carlo initial equilibration updates
#define steps  10000000              // number of Monte-Carlo updates
#define stepsPerMeasurement 100      // number of Monte-Carlo updates per measurement
#define beta   1.0                   // inverse temperature

//
// Below is the list of standard observables:
//

#define MEASURE_H                    // <H>             is measured when this line is not commented
#define MEASURE_H2                   // <H^2>           is measured when this line is not commented
#define MEASURE_HDIAG                // <H_{diag}>      is measured when this line is not commented
#define MEASURE_HDIAG2               // <H_{diag}^2>    is measured when this line is not commented
#define MEASURE_HOFFDIAG             // <H_{offdiag}>   is measured when this line is not commented
#define MEASURE_HOFFDIAG2            // <H_{offdiag}^2> is measured when this line is not commented

//
// Below are the implementation parameters:
//

#define qmax     300                 // upper bound for the maximal length of the sequence of permutation operators
#define Nbins    250                 // number of bins for the error estimation via binning analysis
#define EXHAUSTIVE_CYCLE_SEARCH      // comment this line for a more restrictive cycle search
