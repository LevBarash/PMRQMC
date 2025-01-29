//
// This utility analyses "qmc_data_*.dat" files in the current directory
// and outputs slurm muliple program configuration to finalize unfinished calculations
// without allocating resources to the completed calculations
//

#include<iostream>
#include<filesystem>
#include<vector>
#include<string>
#include<cctype>
#include"mainQMC.hpp"

namespace fs = std::filesystem;

#define fin(obj)  { finput.read((char *)& obj, sizeof(obj)); }

int is_unfinished(std::string filename, unsigned long long &mstep){
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
		fin(rng_seed); fin(meanq); fin(maxq);
		finput.close();
		if(!TstepsFinished || measurement_step*stepsPerMeasurement < steps){
			r = 1; mstep = TstepsFinished ? measurement_step : 0;
		}
	}
	return r;
}

int main(){
	if(steps < Nbins*stepsPerMeasurement){
		std::cout << "Error: steps cannot be smaller than Nbins*stepsPerMeasurement." << std::endl;
		exit(1);
	}
	std::vector<std::pair<int, unsigned long long>> threads; int core = 0, MaxIndex = 0; unsigned long long mstep;
	for(const auto & entry : fs::directory_iterator(fs::current_path()))
		if(entry.is_regular_file() && entry.path().extension() == ".dat"){
			std::string filename = entry.path().filename().string();
			if(filename.find("qmc_data_") == 0){
		                std::string middle = filename.substr(9, filename.size() - 13);
				if(!middle.empty() && std::all_of(middle.begin(), middle.end(), ::isdigit)){
					int index = std::stoi(middle); if(MaxIndex < index) MaxIndex = index;
					if(is_unfinished(filename,mstep)) threads.push_back({index,mstep});
				}
			}
		}
	std::sort(threads.begin(), threads.end(), [](const std::pair<int, unsigned long long>& a, const std::pair<int, unsigned long long>& b){ return a.first < b.first; });
	if(threads.size() == 0){
		std::cout << "Error: unfinished calculations not found" << std::endl; exit(1);
	}
	std::ofstream confFile("multiprog.conf");
	std::ofstream arrayFile("multiprog_array.slurm");
	std::ofstream slurmFile("multiprog.slurm");
	std::ofstream infoFile("multiprog_info.txt");
	arrayFile << "#!/bin/bash" << std::endl;
	arrayFile << "#SBATCH --partition=main" << std::endl;
	arrayFile << "#SBATCH --array=";
	infoFile << "List of unfinished calculations:" << std::endl << std::endl;
	for(size_t i=0;i<threads.size();i++){
		confFile << core << "\t./PMRQMC.bin " << threads[i].first << std::endl;
		infoFile << "Job " << core++ << " = MPI process No. " << threads[i].first << ": main steps completed: ";
		infoFile << threads[i].second*stepsPerMeasurement << " of " << steps << "." << std::endl;
		arrayFile << threads[i].first; if(i<threads.size()-1) arrayFile << ",";
	}
	arrayFile << std::endl; infoFile << std::endl;
	arrayFile << "#SBATCH --ntasks=1" << std::endl;
	arrayFile << "#SBATCH --time=2-0" << std::endl << std::endl;
	arrayFile << "srun ./PMRQMC.bin" << std::endl;
	slurmFile << "#!/bin/bash" << std::endl;
	slurmFile << "#SBATCH --partition=main" << std::endl;
	slurmFile << "#SBATCH --ntasks=" << threads.size() << std::endl;
	slurmFile << "#SBATCH --time=2-0" << std::endl << std::endl;
	slurmFile << "srun --output=\"multiprog_%t.out\" --multi-prog multiprog.conf" << std::endl;
	infoFile << "To complete unfinished calculations without allocating resources to those already completed, use one of the two Slurm scripts:" << std::endl<< "multiprog.slurm (to run the calculations as a single job) or multiprog_array.slurm (to run them as separate jobs)." << std::endl << std::endl;
	infoFile << "To perform postprocessing for the completed calculations, run datasummary utility" << std::endl << std::endl;
	confFile.close(); arrayFile.close(); slurmFile.close(); infoFile.close();
	return 0;
}
