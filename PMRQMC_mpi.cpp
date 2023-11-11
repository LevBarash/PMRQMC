//
// This program implements Permutation Matrix Representation Quantum Monte Carlo for arbitrary spin-1/2 Hamiltonians.
//
// This program is introduced in the paper:
// Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians (2023).
//
// This program is licensed under a Creative Commons Attribution 4.0 International License:
// http://creativecommons.org/licenses/by/4.0/
//
// ExExFloat datatype and calculation of divided differences are described in the paper:
// L. Gupta, L. Barash, I. Hen, Calculating the divided differences of the exponential function by addition and removal of inputs, Computer Physics Communications 254, 107385 (2020)
//

#include"mainqmc.hpp"
#include<fstream>
#include<mpi.h>

double start_time, elapsed_time;

void compute(){
	start_time = MPI_Wtime(); 
	divdiff dd(q+4,500); d=&dd; init();
	for(step=0;step<Tsteps;step++) update();
	for(measurement_step=0;measurement_step<measurements;measurement_step++){
		for(step=0;step<stepsPerMeasurement;step++) update(); measure();
	}
	meanq /= measurements;
	elapsed_time = MPI_Wtime()-start_time;
}

void printout(){
	double Rsum[N_all_observables] = {0}; double sgn_sum = 0; int i,k,o=0;
	double over_bins_sum[N_all_observables] = {0}; double over_bins_sum_sgn = 0;
	double over_bins_sum_cov[N_all_observables] = {0}; double mean[N_all_observables]; double stdev[N_all_observables];
	for(i=0;i<Nbins;i++) sgn_sum += bin_mean_sgn[i]; sgn_sum /= Nbins;
	for(i=0;i<Nbins;i++) over_bins_sum_sgn += (bin_mean_sgn[i] - sgn_sum)*(bin_mean_sgn[i] - sgn_sum); over_bins_sum_sgn /= (Nbins*(Nbins-1));
	std::cout << std::setprecision(9);
	std::cout << "mean(sgn(W)) = " << sgn_sum << std::endl;
	std::cout << "std.dev.(sgn(W)) = " << sqrt(over_bins_sum_sgn) << std::endl;
	if(qmax_achieved) std::cout << std::endl << "Warning: qmax = " << qmax << " was achieved. The results may be incorrect. The qmax parameter should be increased." << std::endl;
	for(i=0;i<Ncycles;i++) if(!cycles_used[i]) std::cout << "Warning: cycle No. " << i << " of length " << cycle_len[i] << " was not used" << std::endl;
	std::cout << "mean(q) = " << meanq << std::endl;
	std::cout << "max(q) = "<< maxq << std::endl;
	for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
		std::cout << "Observable #" << ++o << ": "<< name_of_observable(k) << std::endl;
		for(i=0;i<Nbins;i++) Rsum[k] += bin_mean[k][i]; Rsum[k] /= Nbins;
		for(i=0;i<Nbins;i++) over_bins_sum[k] += (bin_mean[k][i] - Rsum[k])*(bin_mean[k][i] - Rsum[k]); over_bins_sum[k] /= (Nbins*(Nbins-1));
		for(i=0;i<Nbins;i++) over_bins_sum_cov[k] += (bin_mean[k][i] - Rsum[k])*(bin_mean_sgn[i] - sgn_sum); over_bins_sum_cov[k] /= (Nbins*(Nbins-1));
		mean[k] = Rsum[k]/sgn_sum*(1 + over_bins_sum_sgn/sgn_sum/sgn_sum) - over_bins_sum_cov[k]/sgn_sum/sgn_sum;
		stdev[k] = fabs(Rsum[k]/sgn_sum)*sqrt(over_bins_sum[k]/Rsum[k]/Rsum[k] + over_bins_sum_sgn/sgn_sum/sgn_sum - 2*over_bins_sum_cov[k]/Rsum[k]/sgn_sum);
		std::cout << "mean(O) = " << mean[k] << std::endl;
		std::cout << "std.dev.(O) = " << stdev[k] << std::endl;
	}
	std::cout << "Elapsed cpu time = " << elapsed_time << " seconds" << std::endl;
}


int    gathered_qmax_achieved;
double gathered_elapsed_time;
double gathered_meanq;
double gathered_maxq;
double gathered_bin_mean_sgn[Nbins];
double gathered_bin_mean[Nbins];

int main(int argc, char* argv[]){
	int i, k, mpi_rank, mpi_size; divdiff_init();
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	if(steps < Nbins*stepsPerMeasurement){
		std::cout << "Error: steps cannot be smaller than Nbins*stepsPerMeasurement." << std::endl;
		exit(1);
	}
	std::cout << "Starting calculation for MPI process No. " << mpi_rank << std::endl; fflush(stdout);
	compute();
	std::cout << "Calculation completed for MPI process No. " << mpi_rank
	          << ", elapsed time = " << elapsed_time << " seconds, RNG seed = " << rng_seed << std::endl; fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	if(mpi_rank == 0){
		std::cout << std::endl;
		std::cout << "Parameters: beta = " << beta << ", Tsteps = " << Tsteps << ", steps = " << steps << std::endl << std::endl;
		std::cout << "Number of MPI processes: " << mpi_size << std::endl;
		std::cout << std::endl << "Output of the MPI process No. 0:" << std::endl;
		printout();
		std::cout << std::endl;
	}
	double Rsum[N_all_observables] = {0}; double sgn_sum = 0; int o = 0;
	double over_bins_sum[N_all_observables] = {0}; double over_bins_sum_sgn = 0;
	double over_bins_sum_cov[N_all_observables] = {0}; double mean[N_all_observables]; double stdev[N_all_observables];
	MPI_Barrier(MPI_COMM_WORLD);
	if(mpi_rank == 0) std::cout << "Total number of MC updates = " << steps*(unsigned long long)mpi_size << std::endl;
	MPI_Reduce(&elapsed_time,&gathered_elapsed_time,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&meanq,&gathered_meanq,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); gathered_meanq /= mpi_size;
	if(mpi_rank == 0) std::cout << "Total mean(q) = " << gathered_meanq << std::endl;
	MPI_Reduce(&maxq,&gathered_maxq,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	if(mpi_rank == 0) std::cout << "Total max(q) = " << gathered_maxq << std::endl;
	MPI_Reduce(&qmax_achieved,&gathered_qmax_achieved,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
	if(mpi_rank == 0 && gathered_qmax_achieved) std::cout << "Warning: qmax = " << qmax << " was achieved by at least one of the MPI processes. The results may be incorrect. The qmax parameter should be increased." << std::endl;
	MPI_Reduce(bin_mean_sgn,gathered_bin_mean_sgn,Nbins,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	if(mpi_rank == 0){
		for(i=0;i<Nbins;i++) gathered_bin_mean_sgn[i] /= mpi_size;
		for(i=0;i<Nbins;i++) sgn_sum += gathered_bin_mean_sgn[i]; sgn_sum /= Nbins;
		for(i=0;i<Nbins;i++) over_bins_sum_sgn += (gathered_bin_mean_sgn[i] - sgn_sum)*(gathered_bin_mean_sgn[i] - sgn_sum); over_bins_sum_sgn /= (Nbins*(Nbins-1));
		std::cout << std::setprecision(9);
		std::cout << "Total mean(sgn(W)) = " << sgn_sum << std::endl;
		std::cout << "Total std.dev.(sgn(W)) = " << sqrt(over_bins_sum_sgn) << std::endl;
	}
	for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
		MPI_Reduce(bin_mean[k],gathered_bin_mean,Nbins,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		if(mpi_rank == 0){
			std::cout << "Total of observable #" << ++o << ": "<< name_of_observable(k) << std::endl;
			for(i=0;i<Nbins;i++) gathered_bin_mean[i] /= mpi_size;
			for(i=0;i<Nbins;i++) Rsum[k] += gathered_bin_mean[i]; Rsum[k] /= Nbins;
			for(i=0;i<Nbins;i++) over_bins_sum[k] += (gathered_bin_mean[i] - Rsum[k])*(gathered_bin_mean[i] - Rsum[k]); over_bins_sum[k] /= (Nbins*(Nbins-1));
			for(i=0;i<Nbins;i++) over_bins_sum_cov[k] += (gathered_bin_mean[i] - Rsum[k])*(gathered_bin_mean_sgn[i] - sgn_sum); over_bins_sum_cov[k] /= (Nbins*(Nbins-1));
			mean[k] = Rsum[k]/sgn_sum*(1 + over_bins_sum_sgn/sgn_sum/sgn_sum) - over_bins_sum_cov[k]/sgn_sum/sgn_sum;
			stdev[k] = fabs(Rsum[k]/sgn_sum)*sqrt(over_bins_sum[k]/Rsum[k]/Rsum[k] + over_bins_sum_sgn/sgn_sum/sgn_sum - 2*over_bins_sum_cov[k]/Rsum[k]/sgn_sum);
			std::cout << "Total mean(O) = " << mean[k] << std::endl;
			std::cout << "Total std.dev.(O) = " << stdev[k] << std::endl;
		}
	}
	if(mpi_rank == 0){
		std::cout << "Total elapsed cpu time = " << gathered_elapsed_time << " seconds" << std::endl;
		elapsed_time = MPI_Wtime()-start_time;
		std::cout << std::endl << "Wall-clock time = " << elapsed_time << " seconds" << std::endl;
	}
	MPI_Finalize();
	divdiff_clear_up();
	return 0;
}
