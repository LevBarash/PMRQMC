//
// This utility analyses "qmc_data_*.dat" files in the current directory
// and outputs detailed summary of the completed calculations
//

#include<iostream>
#include<filesystem>
#include<vector>
#include<string>
#include<cctype>
#include"mainQMC.hpp"

namespace fs = std::filesystem;

double mean_derived_O[N_derived_observables], stdev_derived_O[N_derived_observables];
double jackknife_O[N_derived_observables], jackknife_sum[N_derived_observables], sgn_meanJ, sgn_varianceJ;
double sgn_mean, sgn_variance, sgn_stdev;
double elapsed_time;

std::vector<double> gathered_mean_sgn;
std::vector<double> gathered_mean_O[N_all_observables];
std::vector<double> gathered_mean_derived_O[N_derived_observables];
double gathered_stdev_sgn;
double gathered_stdev_O[N_all_observables];
double gathered_stdev_derived_O[N_derived_observables];
int    gathered_qmax_achieved;
double gathered_meanq;
double gathered_maxq;
double gathered_elapsed_time;
double gathered_bin_mean_sgn[Nbins];
double gathered_bin_mean[N_all_observables][Nbins];

void process_single_run(){
	double Rsum[N_all_observables] = {0}; sgn_mean = 0; int i,j,k,o;
	double over_bins_sum[N_all_observables] = {0}; sgn_variance = 0;
	double over_bins_sum_cov[N_all_observables] = {0};
	for(i=0;i<Nbins;i++) sgn_mean += bin_mean_sgn[i]; sgn_mean /= Nbins;
	for(i=0;i<Nbins;i++) sgn_variance += (bin_mean_sgn[i] - sgn_mean)*(bin_mean_sgn[i] - sgn_mean); sgn_variance /= (Nbins*(Nbins-1));
	for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
		for(i=0;i<Nbins;i++) Rsum[k] += bin_mean[k][i]; Rsum[k] /= Nbins;
		for(i=0;i<Nbins;i++) over_bins_sum[k] += (bin_mean[k][i] - Rsum[k])*(bin_mean[k][i] - Rsum[k]); over_bins_sum[k] /= (Nbins*(Nbins-1));
		for(i=0;i<Nbins;i++) over_bins_sum_cov[k] += (bin_mean[k][i] - Rsum[k])*(bin_mean_sgn[i] - sgn_mean); over_bins_sum_cov[k] /= (Nbins*(Nbins-1));
		mean_O_backup[k] = mean_O[k] = Rsum[k]/sgn_mean*(1 + sgn_variance/sgn_mean/sgn_mean) - over_bins_sum_cov[k]/sgn_mean/sgn_mean;
		stdev_O[k] = fabs(Rsum[k]/sgn_mean)*sqrt(over_bins_sum[k]/Rsum[k]/Rsum[k] + sgn_variance/sgn_mean/sgn_mean - 2*over_bins_sum_cov[k]/Rsum[k]/sgn_mean);
	}
	for(o=0;o<N_derived_observables;o++) if(valid_derived_observable(o)){
		mean_derived_O[o] = compute_derived_observable(o); jackknife_sum[o] = 0;
	}
	for(j=0;j<Nbins;j++){
		sgn_meanJ = sgn_varianceJ = 0;
		for(i=0;i<Nbins;i++) if(i!=j) sgn_meanJ += bin_mean_sgn[i]; sgn_meanJ /= (Nbins-1);
		for(i=0;i<Nbins;i++) if(i!=j) sgn_varianceJ += (bin_mean_sgn[i] - sgn_meanJ)*(bin_mean_sgn[i] - sgn_meanJ); sgn_varianceJ /= ((Nbins-1)*(Nbins-2));
		for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
			Rsum[k] = over_bins_sum_cov[k] = 0;
			for(i=0;i<Nbins;i++) if(i!=j) Rsum[k] += bin_mean[k][i]; Rsum[k] /= (Nbins-1);
			for(i=0;i<Nbins;i++) if(i!=j) over_bins_sum_cov[k] += (bin_mean[k][i] - Rsum[k])*(bin_mean_sgn[i] - sgn_meanJ); over_bins_sum_cov[k] /= ((Nbins-1)*(Nbins-2));
			mean_O[k] = Rsum[k]/sgn_meanJ*(1 + sgn_varianceJ/sgn_meanJ/sgn_meanJ) - over_bins_sum_cov[k]/sgn_meanJ/sgn_meanJ;
		}
		for(o=0;o<N_derived_observables;o++) if(valid_derived_observable(o)){
			jackknife_O[o] = compute_derived_observable(o);
			jackknife_sum[o] += (jackknife_O[o] - mean_derived_O[o])*(jackknife_O[o] - mean_derived_O[o]);
		}
	}
	for(o=0;o<N_derived_observables;o++) if(valid_derived_observable(o)) stdev_derived_O[o] = sqrt(jackknife_sum[o]*(Nbins-1)/Nbins);
	for(k=0;k<N_all_observables;k++) if(valid_observable[k]) mean_O[k] = mean_O_backup[k];
}

void process_data(int index){ int i,k;
	meanq /= measurements; process_single_run(); sgn_stdev = sqrt(sgn_variance);
	gathered_mean_sgn.push_back(sgn_mean);
	gathered_stdev_sgn += sgn_stdev;
	for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
		gathered_mean_O[k].push_back(mean_O[k]);
		gathered_stdev_O[k] += stdev_O[k];
	}
	for(k=0;k<N_derived_observables;k++) if(valid_derived_observable(k)){
		gathered_mean_derived_O[k].push_back(mean_derived_O[k]);
		gathered_stdev_derived_O[k] += stdev_derived_O[k];
	}
	gathered_meanq += meanq;
	gathered_elapsed_time += elapsed_time;
	if(gathered_maxq < maxq) gathered_maxq = maxq;
	if(gathered_qmax_achieved < qmax_achieved) gathered_qmax_achieved = qmax_achieved;
	for(i=0;i<Nbins;i++) gathered_bin_mean_sgn[i] += bin_mean_sgn[i];
        for(k=0;k<N_all_observables;k++) for(i=0;i<Nbins;i++) gathered_bin_mean[k][i] += bin_mean[k][i];
	mpi_size++;
}

#define fin(obj)  { finput.read((char *)& obj, sizeof(obj)); }

int is_completed(std::string filename){
	std::ifstream finput(filename,std::ios::binary); int r = 0;
	if(finput.good()){
		fin(rng); fin(bin_length_old);
		fin(cycle_len); fin(cycles_used); fin(cycles_used_backup); fin(cycle_min_len); fin(cycle_max_len);
		fin(found_cycles); fin(min_index); fin(max_index); fin(TstepsFinished);
		fin(in_bin_sum); fin(bin_mean); fin(in_bin_sum_sgn); fin(bin_mean_sgn); 
		fin(q); fin(qmax_achieved); fin(lattice); fin(z); fin(P); fin(P_in_cycles);
		fin(Sq); fin(Sq_backup); fin(Sq_subseq); fin(Sq_gaps);
		fin(Energies); fin(Energies_backup); fin(eoccupied); fin(currEnergy); fin(valid_observable);
		fin(currD); fin(old_currD); fin(currD_partial); fin(zero); fin(currWeight); fin(step); fin(measurement_step); 
		fin(rng_seed); fin(meanq); fin(maxq); fin(elapsed_time); if(finput.gcount()==0) elapsed_time = 0;
		finput.close();
		if(TstepsFinished && measurement_step*stepsPerMeasurement >= steps) r = 1;
	}
	return r;
}

int main(int argc, char* argv[]){
	if(steps < Nbins*stepsPerMeasurement){
		std::cout << "Error: steps cannot be smaller than Nbins*stepsPerMeasurement." << std::endl;
		exit(1);
	}
	for(const auto & entry : fs::directory_iterator(fs::current_path()))
		if(entry.is_regular_file() && entry.path().extension() == ".dat"){
			std::string filename = entry.path().filename().string();
			if(filename.find("qmc_data_") == 0){
		                std::string middle = filename.substr(9, filename.size() - 13);
				if(!middle.empty() && std::all_of(middle.begin(), middle.end(), ::isdigit))
					if(is_completed(filename)) process_data(std::stoi(middle));
			}
		}
	if(mpi_size == 0){
		std::cout << "Error: completed calculations not found" << std::endl; exit(1);
	}
	std::cout << std::endl; int i, j, k, o=0; std::cout << std::setprecision(9);
	std::cout << "Parameters: beta = " << beta << ", Tsteps = " << Tsteps << ", steps = " << steps << std::endl << std::endl;
	std::cout << "Number of MPI processes: " << mpi_size << std::endl << std::endl;
	if(mpi_size>4){
		std::cout << "Testing thermalization" << std::endl << std::endl; double mean_mean, std_mean;
		mean_mean = std_mean = 0; gathered_stdev_sgn /= mpi_size;
		for(i=0;i<mpi_size;i++) mean_mean += gathered_mean_sgn[i]; mean_mean /= mpi_size;
		for(i=0;i<mpi_size;i++) std_mean += (gathered_mean_sgn[i] - mean_mean)*(gathered_mean_sgn[i] - mean_mean);
		std_mean /= (mpi_size - 1); std_mean = sqrt(std_mean);
		// std::cout << "mean of std.dev.(sgn(W)) = " << gathered_stdev_sgn << ", std.dev. of mean(sgn(W)) = " << std_mean;
		// if(gathered_stdev_sgn >= 0.7 * std_mean) std::cout << ": test passed" << std::endl; else std::cout << ": test failed" << std::endl;
		for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
			std::cout << "Observable #" << ++o << ": "<< name_of_observable(k);
			mean_mean = std_mean = 0; gathered_stdev_O[k] /= mpi_size;
			for(i=0;i<mpi_size;i++) mean_mean += gathered_mean_O[k][i]; mean_mean /= mpi_size;
			for(i=0;i<mpi_size;i++) std_mean += (gathered_mean_O[k][i] - mean_mean)*(gathered_mean_O[k][i] - mean_mean);
			std_mean /= (mpi_size - 1); std_mean = sqrt(std_mean);
			std::cout << ", mean of std.dev.(O) = " << gathered_stdev_O[k] << ", std.dev. of mean(O) = " << std_mean;
			if(gathered_stdev_O[k] >= 0.7 * std_mean) std::cout << ": test passed" << std::endl; else std::cout << ": test failed" << std::endl;
		}
		for(k=0;k<N_derived_observables;k++) if(valid_derived_observable(k)){
			std::cout << "Derived observable: "<< name_of_derived_observable(k);
			mean_mean = std_mean = 0; gathered_stdev_derived_O[k] /= mpi_size;
			for(i=0;i<mpi_size;i++) mean_mean += gathered_mean_derived_O[k][i]; mean_mean /= mpi_size;
			for(i=0;i<mpi_size;i++) std_mean += (gathered_mean_derived_O[k][i] - mean_mean)*(gathered_mean_derived_O[k][i] - mean_mean);
			std_mean /= (mpi_size - 1); std_mean = sqrt(std_mean);
			std::cout << ", mean of std.dev.(O) = " << gathered_stdev_derived_O[k] << ", std.dev. of mean(O) = " << std_mean;
			if(gathered_stdev_derived_O[k] >= 0.7 * std_mean) std::cout << ": test passed" << std::endl; else std::cout << ": test failed" << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << "Collecting statistics and finalizing the calculation" << std::endl << std::endl;
	double Rsum[N_all_observables] = {0}; sgn_mean = 0; o = 0;
	double over_bins_sum[N_all_observables] = {0}; sgn_variance = 0;
	double over_bins_sum_cov[N_all_observables] = {0};
	std::cout << "Total number of MC updates = " << steps*(unsigned long long)mpi_size << std::endl;
	gathered_meanq /= mpi_size;
	std::cout << "Total mean(q) = " << gathered_meanq << std::endl;
	std::cout << "Total max(q) = " << gathered_maxq << std::endl;
	if(gathered_qmax_achieved) std::cout << "Warning: qmax = " << qmax << " was achieved by at least one of the MPI processes. The results may be incorrect. The qmax parameter should be increased." << std::endl;
	for(i=0;i<Nbins;i++) gathered_bin_mean_sgn[i] /= mpi_size;
	for(i=0;i<Nbins;i++) sgn_mean += gathered_bin_mean_sgn[i]; sgn_mean /= Nbins;
	for(i=0;i<Nbins;i++) sgn_variance += (gathered_bin_mean_sgn[i] - sgn_mean)*(gathered_bin_mean_sgn[i] - sgn_mean); sgn_variance /= (Nbins*(Nbins-1));
	std::cout << std::setprecision(9);
	std::cout << "Total mean(sgn(W)) = " << sgn_mean << std::endl;
	std::cout << "Total std.dev.(sgn(W)) = " << sqrt(sgn_variance) << std::endl;
	for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
		std::cout << "Total of observable #" << ++o << ": "<< name_of_observable(k) << std::endl;
		for(i=0;i<Nbins;i++) gathered_bin_mean[k][i] /= mpi_size;
		for(i=0;i<Nbins;i++) Rsum[k] += gathered_bin_mean[k][i]; Rsum[k] /= Nbins;
		for(i=0;i<Nbins;i++) over_bins_sum[k] += (gathered_bin_mean[k][i] - Rsum[k])*(gathered_bin_mean[k][i] - Rsum[k]); over_bins_sum[k] /= (Nbins*(Nbins-1));
		for(i=0;i<Nbins;i++) over_bins_sum_cov[k] += (gathered_bin_mean[k][i] - Rsum[k])*(gathered_bin_mean_sgn[i] - sgn_mean); over_bins_sum_cov[k] /= (Nbins*(Nbins-1));
		mean_O[k] = Rsum[k]/sgn_mean*(1 + sgn_variance/sgn_mean/sgn_mean) - over_bins_sum_cov[k]/sgn_mean/sgn_mean;
		stdev_O[k] = fabs(Rsum[k]/sgn_mean)*sqrt(over_bins_sum[k]/Rsum[k]/Rsum[k] + sgn_variance/sgn_mean/sgn_mean - 2*over_bins_sum_cov[k]/Rsum[k]/sgn_mean);
		std::cout << "Total mean(O) = " << mean_O[k] << std::endl;
		std::cout << "Total std.dev.(O) = " << stdev_O[k] << std::endl;
	}
	for(o=0;o<N_derived_observables;o++) if(valid_derived_observable(o)){
		mean_derived_O[o] = compute_derived_observable(o); jackknife_sum[o] = 0;
	}
	for(j=0;j<Nbins;j++){
		sgn_meanJ = sgn_varianceJ = 0;
		for(i=0;i<Nbins;i++) if(i!=j) sgn_meanJ += gathered_bin_mean_sgn[i]; sgn_meanJ /= (Nbins-1);
		for(i=0;i<Nbins;i++) if(i!=j) sgn_varianceJ += (gathered_bin_mean_sgn[i] - sgn_meanJ)*(gathered_bin_mean_sgn[i] - sgn_meanJ); sgn_varianceJ /= ((Nbins-1)*(Nbins-2));
		for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
			Rsum[k] = over_bins_sum_cov[k] = 0;
			for(i=0;i<Nbins;i++) if(i!=j) Rsum[k] += gathered_bin_mean[k][i]; Rsum[k] /= (Nbins-1);
			for(i=0;i<Nbins;i++) if(i!=j) over_bins_sum_cov[k] += (gathered_bin_mean[k][i] - Rsum[k])*(gathered_bin_mean_sgn[i] - sgn_meanJ); over_bins_sum_cov[k] /= ((Nbins-1)*(Nbins-2));
			mean_O[k] = Rsum[k]/sgn_meanJ*(1 + sgn_varianceJ/sgn_meanJ/sgn_meanJ) - over_bins_sum_cov[k]/sgn_meanJ/sgn_meanJ;
		}
		for(o=0;o<N_derived_observables;o++) if(valid_derived_observable(o)){
			jackknife_O[o] = compute_derived_observable(o);
			jackknife_sum[o] += (jackknife_O[o] - mean_derived_O[o])*(jackknife_O[o] - mean_derived_O[o]);
		}
	}
	for(o=0;o<N_derived_observables;o++) if(valid_derived_observable(o)){
		stdev_derived_O[o] = sqrt(jackknife_sum[o]*(Nbins-1)/Nbins);
		std::cout << "Total of derived observable: " << name_of_derived_observable(o) << std::endl;
		std::cout << "Total mean(O) = " << mean_derived_O[o] << std::endl;
		std::cout << "Total std.dev.(O) = " << stdev_derived_O[o] << std::endl;
	}
	std::cout << "Total elapsed cpu time = " << gathered_elapsed_time << std::endl;
	return 0;
}
