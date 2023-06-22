#define Tsteps 1000000ull            // number of Monte-Carlo initial equilibration updates (use "ull" suffix for unsigned long long)
#define steps  10000000ull           // number of Monte-Carlo updates (use "ull" suffix for unsigned long long)
#define stepsPerMeasurement 100ull   // number of Monte-Carlo updates per measurement (use "ull" suffix for unsigned long long)
#define qmax     300                 // upper bound for the maximal length of the sequence of permutation operators
#define Nbins    250                 // number of bins for the error estimation via binning analysis
#define exhaustive_cycle_search      // comment this line for a more restrictive cycle search
