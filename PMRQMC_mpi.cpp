//
// This program implements Permutation Matrix Representation Quantum Monte Carlo for arbitrary spin-1/2 Hamiltonians.
//
// This program is introduced in the paper:
// Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians, Physical Review Research 6, 013281 (2024).
//
// This program is licensed under a Creative Commons Attribution 4.0 International License:
// http://creativecommons.org/licenses/by/4.0/
//
// ExExFloat datatype and calculation of divided differences are described in the paper:
// L. Gupta, L. Barash, I. Hen, Calculating the divided differences of the exponential function by addition and removal of inputs, Computer Physics Communications 254, 107385 (2020)
//

#include"mainQMC.hpp"
#include<mpi.h>

double start_time, elapsed_time;
double mean_O[N_all_observables]; double stdev_O[N_all_observables];
double sgn_mean, sgn_variance, sgn_stdev;
void signalHandler(int signum){	if(save_data_flag==0) save_data_flag = 1; }

void compute(){
	start_time = MPI_Wtime(); 
	if(!resume_calc) std::cout << "Starting calculation for MPI process No. " << mpi_rank << ", RNG seed = " << rng_seed << std::endl; fflush(stdout);
	if(TstepsFinished){
		if(step>0 && step<stepsPerMeasurement && measurement_step<measurements){
			for(;step<stepsPerMeasurement;step++) update(); measure(); measurement_step++;
		}
	} else{
		for(;step<Tsteps;step++) update(); TstepsFinished = 1;
	}
	for(;measurement_step<measurements;measurement_step++){
		for(step=0;step<stepsPerMeasurement;step++) update(); measure();
	}
#ifdef SAVE_COMPLETED_CALCULATION
	save_QMC_data(0);
#endif
	meanq /= measurements;
	elapsed_time = MPI_Wtime()-start_time;
	std::cout << "Calculation completed for MPI process No. " << mpi_rank
	          << ", elapsed time = " << elapsed_time << " seconds" << std::endl; fflush(stdout);
}

void process_single_run(){
	double Rsum[N_all_observables] = {0}; sgn_mean = 0; int i,k;
	double over_bins_sum[N_all_observables] = {0}; sgn_variance = 0;
	double over_bins_sum_cov[N_all_observables] = {0};
	for(i=0;i<Nbins;i++) sgn_mean += bin_mean_sgn[i]; sgn_mean /= Nbins;
	for(i=0;i<Nbins;i++) sgn_variance += (bin_mean_sgn[i] - sgn_mean)*(bin_mean_sgn[i] - sgn_mean); sgn_variance /= (Nbins*(Nbins-1));
	for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
		for(i=0;i<Nbins;i++) Rsum[k] += bin_mean[k][i]; Rsum[k] /= Nbins;
		for(i=0;i<Nbins;i++) over_bins_sum[k] += (bin_mean[k][i] - Rsum[k])*(bin_mean[k][i] - Rsum[k]); over_bins_sum[k] /= (Nbins*(Nbins-1));
		for(i=0;i<Nbins;i++) over_bins_sum_cov[k] += (bin_mean[k][i] - Rsum[k])*(bin_mean_sgn[i] - sgn_mean); over_bins_sum_cov[k] /= (Nbins*(Nbins-1));
		mean_O[k] = Rsum[k]/sgn_mean*(1 + sgn_variance/sgn_mean/sgn_mean) - over_bins_sum_cov[k]/sgn_mean/sgn_mean;
		stdev_O[k] = fabs(Rsum[k]/sgn_mean)*sqrt(over_bins_sum[k]/Rsum[k]/Rsum[k] + sgn_variance/sgn_mean/sgn_mean - 2*over_bins_sum_cov[k]/Rsum[k]/sgn_mean);
	}
}

void printout_single_run(){
	std::cout << std::setprecision(9); int i,k,o=0;
	std::cout << "mean(sgn(W)) = " << sgn_mean << std::endl;
	std::cout << "std.dev.(sgn(W)) = " << sqrt(sgn_variance) << std::endl;
	if(qmax_achieved) std::cout << std::endl << "Warning: qmax = " << qmax << " was achieved. The results may be incorrect. The qmax parameter should be increased." << std::endl;
	for(i=0;i<Ncycles;i++) if(!cycles_used[i]) std::cout << "Warning: cycle No. " << i << " of length " << cycle_len[i] << " was not used" << std::endl;
	std::cout << "mean(q) = " << meanq << std::endl;
	std::cout << "max(q) = "<< maxq << std::endl;
	for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
		std::cout << "Observable #" << ++o << ": "<< name_of_observable(k) << std::endl;
		std::cout << "mean(O) = " << mean_O[k] << std::endl;
		std::cout << "std.dev.(O) = " << stdev_O[k] << std::endl;
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
#ifdef SAVE_UNFINISHED_CALCULATION
	signal(SIGTERM,signalHandler);
#endif
	int i, k, o=0; divdiff_init();
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	if(steps < Nbins*stepsPerMeasurement){
		std::cout << "Error: steps cannot be smaller than Nbins*stepsPerMeasurement." << std::endl;
		MPI_Finalize(); exit(1);
	}
	if(N == 0){
		std::cout << "Error: no particles found. At least one particle must be described by the Hamiltonian." << std::endl;
		MPI_Finalize(); exit(1);
	}

	if(mpi_rank == 0) resume_calc = check_QMC_data();
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&resume_calc,1,MPI_INT,0,MPI_COMM_WORLD); init_rng();
	divdiff dd(q+4,500); d=&dd;
	if(resume_calc){ load_QMC_data(); init_basic(); } else init();
	MPI_Barrier(MPI_COMM_WORLD);
	compute(); process_single_run();
	MPI_Barrier(MPI_COMM_WORLD);
	if(mpi_rank == 0){
		std::cout << std::endl;
		std::cout << "Parameters: beta = " << beta << ", Tsteps = " << Tsteps << ", steps = " << steps << std::endl << std::endl;
		std::cout << "Number of MPI processes: " << mpi_size << std::endl;
		std::cout << std::endl << "Output of the MPI process No. 0:" << std::endl << std::endl;
		printout_single_run();
		std::cout << std::endl;
	}
	if(mpi_size>4){
		if(mpi_rank == 0) std::cout << "Testing thermalization" << std::endl << std::endl;
		MPI_Barrier(MPI_COMM_WORLD);
		double* gathered_mean = new double[mpi_size]; double gathered_stdev, mean_mean, std_mean; sgn_stdev = sqrt(sgn_variance);
		MPI_Gather(&sgn_mean,1,MPI_DOUBLE,gathered_mean,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Reduce(&sgn_stdev,&gathered_stdev,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		if(mpi_rank == 0){
			mean_mean = std_mean = 0; gathered_stdev /= mpi_size;
			for(i=0;i<mpi_size;i++) mean_mean += gathered_mean[i]; mean_mean /= mpi_size;
			for(i=0;i<mpi_size;i++) std_mean += (gathered_mean[i] - mean_mean)*(gathered_mean[i] - mean_mean);
			std_mean /= (mpi_size - 1); std_mean = sqrt(std_mean);
			// std::cout << "mean of std.dev.(sgn(W)) = " << gathered_stdev << ", std.dev. of mean(sgn(W)) = " << std_mean;
			// if(gathered_stdev >= 0.7 * std_mean) std::cout << ": test passed" << std::endl; else std::cout << ": test failed" << std::endl;
		}
		for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
			MPI_Gather(&mean_O[k],1,MPI_DOUBLE,gathered_mean,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Reduce(&stdev_O[k],&gathered_stdev,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			if(mpi_rank == 0){
				std::cout << "Observable #" << ++o << ": "<< name_of_observable(k);
				mean_mean = std_mean = 0; gathered_stdev /= mpi_size;
				for(i=0;i<mpi_size;i++) mean_mean += gathered_mean[i]; mean_mean /= mpi_size;
				for(i=0;i<mpi_size;i++) std_mean += (gathered_mean[i] - mean_mean)*(gathered_mean[i] - mean_mean);
				std_mean /= (mpi_size - 1); std_mean = sqrt(std_mean);
				std::cout << ", mean of std.dev.(O) = " << gathered_stdev << ", std.dev. of mean(O) = " << std_mean;
				if(gathered_stdev >= 0.7 * std_mean) std::cout << ": test passed" << std::endl; else std::cout << ": test failed" << std::endl;
			}
		}
		delete[] gathered_mean; if(mpi_rank == 0) std::cout << std::endl;
	}
	if(mpi_rank == 0) std::cout << "Collecting statistics and finalizing the calculation" << std::endl << std::endl;
	double Rsum[N_all_observables] = {0}; sgn_mean = 0; o = 0;
	double over_bins_sum[N_all_observables] = {0}; sgn_variance = 0;
	double over_bins_sum_cov[N_all_observables] = {0};
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
		for(i=0;i<Nbins;i++) sgn_mean += gathered_bin_mean_sgn[i]; sgn_mean /= Nbins;
		for(i=0;i<Nbins;i++) sgn_variance += (gathered_bin_mean_sgn[i] - sgn_mean)*(gathered_bin_mean_sgn[i] - sgn_mean); sgn_variance /= (Nbins*(Nbins-1));
		std::cout << std::setprecision(9);
		std::cout << "Total mean(sgn(W)) = " << sgn_mean << std::endl;
		std::cout << "Total std.dev.(sgn(W)) = " << sqrt(sgn_variance) << std::endl;
	}
	for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
		MPI_Reduce(bin_mean[k],gathered_bin_mean,Nbins,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		if(mpi_rank == 0){
			std::cout << "Total of observable #" << ++o << ": "<< name_of_observable(k) << std::endl;
			for(i=0;i<Nbins;i++) gathered_bin_mean[i] /= mpi_size;
			for(i=0;i<Nbins;i++) Rsum[k] += gathered_bin_mean[i]; Rsum[k] /= Nbins;
			for(i=0;i<Nbins;i++) over_bins_sum[k] += (gathered_bin_mean[i] - Rsum[k])*(gathered_bin_mean[i] - Rsum[k]); over_bins_sum[k] /= (Nbins*(Nbins-1));
			for(i=0;i<Nbins;i++) over_bins_sum_cov[k] += (gathered_bin_mean[i] - Rsum[k])*(gathered_bin_mean_sgn[i] - sgn_mean); over_bins_sum_cov[k] /= (Nbins*(Nbins-1));
			mean_O[k] = Rsum[k]/sgn_mean*(1 + sgn_variance/sgn_mean/sgn_mean) - over_bins_sum_cov[k]/sgn_mean/sgn_mean;
			stdev_O[k] = fabs(Rsum[k]/sgn_mean)*sqrt(over_bins_sum[k]/Rsum[k]/Rsum[k] + sgn_variance/sgn_mean/sgn_mean - 2*over_bins_sum_cov[k]/Rsum[k]/sgn_mean);
			std::cout << "Total mean(O) = " << mean_O[k] << std::endl;
			std::cout << "Total std.dev.(O) = " << stdev_O[k] << std::endl;
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
