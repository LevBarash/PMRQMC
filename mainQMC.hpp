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

#include<iostream>
#include<fstream>
#include<iomanip>
#include<complex>
#include<random>
#include<cstdlib>
#include<algorithm>
#include<csignal>
#include<bitset>
#include"divdiff.hpp"
#include"hamiltonian.hpp" // use a header file, which defines the Hamiltonian and the custom observables
#include"parameters.hpp"  // parameters of the simulation such as the number of Monte-Carlo updates

#define measurements (steps/stepsPerMeasurement)

#ifdef EXHAUSTIVE_CYCLE_SEARCH
	#define rmin 0                            // length r of sub-sequence is chosen randomly between rmin and rmax
	#define rmax cycle_max_len
	#define lmin r                            // cycle lengths must be between lmin and lmax
	#define lmax cycle_max_len
#else
	#define rmin (cycle_min_len-1)/2
	#define rmax (cycle_max_len+1)/2
	#define lmin 2*r-1
	#define lmax 2*r+1
#endif

static std::random_device rd;
static std::mt19937 rng;
static std::uniform_int_distribution<> dice2(0,1);
static std::uniform_int_distribution<> diceN(0,N-1);
static std::uniform_int_distribution<> diceNop(0,Nop-1);
static std::uniform_real_distribution<> val(0.0,1.0);
static std::geometric_distribution<> geometric_int(GAPS_GEOMETRIC_PARAMETER);

ExExFloat beta_pow_factorial[qmax]; // contains the values (-beta)^q / q!
double factorial[qmax]; // contains the values q!
int cycle_len[Ncycles];
int cycles_used[Ncycles];
int cycles_used_backup[Ncycles];
int cycle_min_len, cycle_max_len, found_cycles, min_index, max_index;

#ifndef MEASURE_CUSTOM_OBSERVABLES
#define Nobservables 0
#endif

const int N_all_observables = Nobservables + 8;
int valid_observable[N_all_observables];

unsigned long long bin_length = measurements / Nbins, bin_length_old;
double in_bin_sum[N_all_observables];
double bin_mean[N_all_observables][Nbins];
double in_bin_sum_sgn;
double bin_mean_sgn[Nbins];

int q;
int qmax_achieved=0;

divdiff *d, *dfs, *ds1, *ds2;

std::bitset<N> lattice;
std::bitset<N> z;
std::bitset<Nop> P;
std::bitset<Ncycles> P_in_cycles[Nop];

int Sq[qmax];	// list of q operator indices
int Sq_backup[qmax];
int Sq_subseq[qmax];
int Sq_gaps[qmax];
double Energies[qmax+1];
double Energies_backup[qmax+1];
int eoccupied[qmax+1];
double currEnergy;
std::complex<double> old_currD, currD;
std::complex<double> currD_partial[qmax];

#ifdef ABS_WEIGHTS
#define REALW    std::abs
#else
#define REALW    std::real
#endif

ExExFloat zero, currWeight; int TstepsFinished = 0;
unsigned long long step = 0;
unsigned long long measurement_step = 0;

unsigned int rng_seed;
double meanq = 0;
double maxq = 0;
double start_time;
int save_data_flag = 0, mpi_rank = 0, mpi_size = 0, resume_calc = 0;

#define fout(obj) { foutput.write((char *)& obj, sizeof(obj)); }
#define fin(obj)  { finput.read((char *)& obj, sizeof(obj)); }

void save_QMC_data(int printout = 1){
	if(printout) std::cout<<"SIGTERM signal detected. Saving unfinished calculation...";
	char fname[100]; if(mpi_size>0) sprintf(fname,"qmc_data_%d.dat",mpi_rank); else sprintf(fname,"qmc_data.dat");
	std::ofstream foutput(fname,std::ios::binary); double elapsed;
	fout(rng); fout(bin_length);
	fout(cycle_len); fout(cycles_used); fout(cycles_used_backup); fout(cycle_min_len); fout(cycle_max_len);
	fout(found_cycles); fout(min_index); fout(max_index); fout(TstepsFinished);
	fout(in_bin_sum); fout(bin_mean); fout(in_bin_sum_sgn); fout(bin_mean_sgn); 
	fout(q); fout(qmax_achieved); fout(lattice); fout(z); fout(P); fout(P_in_cycles);
	fout(Sq); fout(Sq_backup); fout(Sq_subseq); fout(Sq_gaps);
	fout(Energies); fout(Energies_backup); fout(eoccupied); fout(currEnergy); fout(valid_observable);
	fout(currD); fout(old_currD); fout(currD_partial); fout(zero); fout(currWeight); fout(step); fout(measurement_step); 
	fout(rng_seed); fout(meanq); fout(maxq);
#ifdef MPI_VERSION
	elapsed = MPI_Wtime() - start_time;
#else
	elapsed = (double)clock() / CLOCKS_PER_SEC - start_time;
#endif
	fout(elapsed);
	foutput.close(); if(printout) std::cout<<"done"<<std::endl; fflush(stdout);
}

int check_QMC_data(){
#ifdef RESUME_CALCULATION
	int i,r,g; char fname[100];
	if(mpi_rank==0 && mpi_size>0){
		r = 1;
		for(i=0;i<mpi_size;i++){
			sprintf(fname,"qmc_data_%d.dat",mpi_rank);
			std::ifstream finput(fname,std::ios::binary);
			g = finput.good(); if(g) finput.close(); r = r && g;
		}
	} else{
		if(mpi_size>0) sprintf(fname,"qmc_data_%d.dat",mpi_rank); else sprintf(fname,"qmc_data.dat");
		std::ifstream finput(fname,std::ios::binary);
		r = finput.good(); if(r) finput.close();
	}
#else
	int r = 0;
#endif
	return r;
}

void load_QMC_data(){
	char fname[100]; if(mpi_size>0) sprintf(fname,"qmc_data_%d.dat",mpi_rank); else sprintf(fname,"qmc_data.dat");
	std::ifstream finput(fname,std::ios::binary); double elapsed;
	if(finput.good()){
		fin(rng); fin(bin_length_old);
		fin(cycle_len); fin(cycles_used); fin(cycles_used_backup); fin(cycle_min_len); fin(cycle_max_len);
		fin(found_cycles); fin(min_index); fin(max_index); fin(TstepsFinished);
		fin(in_bin_sum); fin(bin_mean); fin(in_bin_sum_sgn); fin(bin_mean_sgn); 
		fin(q); fin(qmax_achieved); fin(lattice); fin(z); fin(P); fin(P_in_cycles);
		fin(Sq); fin(Sq_backup); fin(Sq_subseq); fin(Sq_gaps);
		fin(Energies); fin(Energies_backup); fin(eoccupied); fin(currEnergy); fin(valid_observable);
		fin(currD); fin(old_currD); fin(currD_partial); fin(zero); fin(currWeight); fin(step); fin(measurement_step); 
		fin(rng_seed); fin(meanq); fin(maxq); fin(elapsed); if(finput.gcount()==0) elapsed = 0;
		finput.close(); start_time -= elapsed;
		if(mpi_size > 0){
			if(TstepsFinished){
				std::cout<<"MPI process No. "<< mpi_rank <<" loaded unfinished calculation: all Tsteps completed, main steps completed: "
				    <<measurement_step*stepsPerMeasurement<<" of "<<steps<<"."<<std::endl;
			} else std::cout <<"MPI process No. "<< mpi_rank <<" loaded unfinished calculation: Tsteps completed: "<<step<<" of "<<Tsteps<<"."<<std::endl;
		} else{
			if(TstepsFinished){
				std::cout<<"Loaded unfinished calculation: all Tsteps completed, main steps completed: "
				    <<measurement_step*stepsPerMeasurement<<" of "<<steps<<"."<<std::endl;
			} else std::cout <<"Loaded unfinished calculation: Tsteps completed: "<<step<<" of "<<Tsteps<<"."<<std::endl;
		}
	} else{
		if(mpi_size>0) std::cout<<"MPI process No. "<<mpi_rank<<": error opening file "<<fname<<std::endl;
		else std::cout<<"Error opening file "<<fname<<std::endl; fflush(stdout);
#ifdef MPI_VERSION
		MPI_Abort(MPI_COMM_WORLD,1);
#else
		exit(1);
#endif
	}
	if(bin_length != bin_length_old){ // It is allowed to increase steps by an integer number of times for a completed calculation
		if(bin_length > 0 && bin_length % bin_length_old == 0){ // All other parameters should remain unchanged
			double sum; int i, j, o, m = bin_length / bin_length_old, curr_bins = Nbins/m; // merging bins together
			for(o=0;o<N_all_observables;o++) for(i=0;i<curr_bins;i++){
				sum = 0; for(j=0;j<m;j++) sum += bin_mean[o][m*i+j]; sum /= m;
				bin_mean[o][i] = sum;
			}
			for(i=0;i<curr_bins;i++){
				sum = 0; for(j=0;j<m;j++) sum += bin_mean_sgn[m*i+j]; sum /= m;
				bin_mean_sgn[i] = sum;
			}
			for(o=0;o<N_all_observables;o++){
				sum = 0; for(i=curr_bins*m;i<Nbins;i++) sum += bin_mean[o][i]*bin_length_old;
				in_bin_sum[o] += sum;
			}
			sum = 0; for(i=curr_bins*m;i<Nbins;i++) sum += bin_mean_sgn[i]*bin_length_old;
			in_bin_sum_sgn += sum;
		} else{
			std::cout << "Error: bin_length = " << bin_length <<" is not divisible by bin_length_old = " << bin_length_old << std::endl; fflush(stdout);
#ifdef MPI_VERSION
			MPI_Abort(MPI_COMM_WORLD,1);
#else
			exit(1);
#endif
		}
	}
}

double CalcEnergy(){ // calculate the energy <z | D_0 | z> of a given configuration of spins
	std::complex<double> sum = 0;
	for(int i=0;i<D0_size;i++) sum -= double(2*(int((D0_product[i] & (~lattice)).count())%2)-1) * D0_coeff[i];
	return sum.real();
}

std::complex<double> calc_d(int k){ // calculate d_k = <z | D_k | z> for the current configuration of spins
	std::complex<double> sum = 0;
	for(int i=0;i<D_size[k];i++) sum -= double(2*(int((D_product[k][i] & (~lattice)).count())%2)-1) * D_coeff[k][i];
	return sum;
}

void ApplyOperator(int k){
	lattice ^= P_matrix[k];
}

void GetEnergies(){
	currD = currD_partial[0] = 1;
	for(int i=0;i<q;i++){
		Energies[i] = CalcEnergy();
		ApplyOperator(Sq[i]);
		currD *= calc_d(Sq[i]);
		currD_partial[i+1] = currD;
	}
	currEnergy = Energies[q] = CalcEnergy();
}

ExExFloat GetWeight(){
	d->CurrentLength=0; GetEnergies();
	for(int i=0;i<=q;i++) d->AddElement(-beta*Energies[i]);
	return d->divdiffs[q] * beta_pow_factorial[q] * REALW(currD);
}

ExExFloat UpdateWeight(){
	int i, j, notfound, n=d->CurrentLength; double value;
	GetEnergies(); memset(eoccupied,0,(q+1)*sizeof(int));
	for(i=0;i<n;i++){
	        notfound = 1; value = d->z[i];
		for(j=0;j<=q;j++) if(eoccupied[j]==0 && value == -beta*Energies[j]){ notfound = 0; break; }
		if(notfound) break; eoccupied[j] = 1;
	}
	if(i==0) d->CurrentLength=0; else while(n>i){ d->RemoveElement(); n--; }
	j=0; while(i<=q){ while(eoccupied[j]) j++; d->AddElement(-beta*Energies[j++]); i++; }
        return d->divdiffs[q] * beta_pow_factorial[q] * REALW(currD);
}

ExExFloat UpdateWeightReplace(double removeEnergy, double addEnergy){
	if(removeEnergy != addEnergy){
		if(d->RemoveValue(-beta*removeEnergy)) d->AddElement(-beta*addEnergy); else{
			std::cout << "Error: energy not found" << std::endl; exit(1);
		}
	}
	return d->divdiffs[q] * beta_pow_factorial[q] * REALW(currD); // use this value only when the values of q and currD are correct
}

ExExFloat UpdateWeightDel(double removeEnergy1, double removeEnergy2){
	if(d->RemoveValue(-beta*removeEnergy1) && d->RemoveValue(-beta*removeEnergy2))
		return d->divdiffs[q] * beta_pow_factorial[q] * REALW(currD);  // use this value only when the values of q and currD are correct
	else{
		std::cout << "Error: energy not found" << std::endl; exit(1);
	}
}

ExExFloat UpdateWeightIns(double addEnergy1, double addEnergy2){
	d->AddElement(-beta*addEnergy1); d->AddElement(-beta*addEnergy2);
	return d->divdiffs[q] * beta_pow_factorial[q] * REALW(currD);    // use this value only when the values of q and currD are correct
}

int NoRepetitionCheck(int* sequence, int r){ // check for absence of repetitions in a sequence of length r
	int i,j,rep = 1;
	for(i=0;i<r && rep;i++) for(j=0;j<i;j++) if(sequence[j]==sequence[i]){ rep = 0; break;}
	return rep;
}

void PickSubsequence(int r){ // randomly picks a sequential sub-sequence of length r from Sq
	int i,m; m = int(val(rng)*(q-r+1)); // m is random integer between 0 and q-r
	for(i=0;i<r;i++) Sq_subseq[i] = Sq[i+m];
	min_index = m; max_index = m+r-1;
}

int FindCycles(int r){  // find all cycles of length between lmin and lmax, each containing all operators of the array Sq_subseq of length r.
	int i,k,sum; std::bitset<Ncycles> curr; curr.set();
	for(i=0;i<r;i++) curr &= P_in_cycles[Sq_subseq[i]];
#ifndef EXHAUSTIVE_CYCLE_SEARCH
	for(i=0;i<Ncycles;i++) if(curr[i]) if(cycle_len[i]<lmin || cycle_len[i]>lmax) curr.reset(i);
#endif
	found_cycles = curr.count();
	if(found_cycles > 0){
		k = int(val(rng)*found_cycles);
		i=sum=0; while(sum <= k) sum += curr[i++]; i--;
		return i;
	} else return -1;
}

void init_rng(){
#ifdef EXACTLY_REPRODUCIBLE
	rng_seed = 1000 + mpi_rank;
#else
	rng_seed = rd();
#endif
	rng.seed(rng_seed);
}

void init_basic(){
	double curr2=1; ExExFloat curr1; beta_pow_factorial[0] = curr1; factorial[0] = curr2;
	for(int k=1;k<qmax;k++){ curr1*=(-double(beta))/k; curr2*=k; beta_pow_factorial[k] = curr1; factorial[k] = curr2;}
	zero -= zero; currWeight = GetWeight();
}

void init(){
	int i,j; double curr2=1; ExExFloat curr1; beta_pow_factorial[0] = curr1; factorial[0] = curr2;
	for(q=1;q<qmax;q++){ curr1*=(-double(beta))/q; curr2*=q; beta_pow_factorial[q] = curr1; factorial[q] = curr2;}
	zero -= zero;
	lattice = 0; for(i=N-1;i>=0;i--) if(dice2(rng)) lattice.set(i); z = lattice; q=0;
	currWeight = GetWeight();
	for(i=0;i<Ncycles;i++) cycle_len[i] = cycles[i].count();
	cycle_min_len = 64; for(i=0;i<Ncycles;i++) cycle_min_len = min(cycle_min_len,cycle_len[i]);
	cycle_max_len = 0; for(i=0;i<Ncycles;i++) cycle_max_len = max(cycle_max_len,cycle_len[i]);
	for(i=0;i<Ncycles;i++) cycles_used[i] = 0;
	for(i=0;i<Nop;i++) for(j=0;j<Ncycles;j++) if(cycles[j].test(i)) P_in_cycles[i].set(j);
	for(i=0;i<N_all_observables;i++) in_bin_sum[i] = 0; in_bin_sum_sgn = 0;
	for(i=0;i<N_all_observables;i++) valid_observable[i] = 0;
#ifdef MEASURE_CUSTOM_OBSERVABLES
	for(i=0;i<Nobservables;i++) valid_observable[i] = 1;
#endif
#ifdef MEASURE_H
	valid_observable[Nobservables] = 1;
#endif
#ifdef MEASURE_H2
	valid_observable[Nobservables + 1] = 1;
#endif
#ifdef MEASURE_HDIAG
	valid_observable[Nobservables + 2] = 1;
#endif
#ifdef MEASURE_HDIAG2
	valid_observable[Nobservables + 3] = 1;
#endif
#ifdef MEASURE_HOFFDIAG
	valid_observable[Nobservables + 4] = 1;
#endif
#ifdef MEASURE_HOFFDIAG2
	valid_observable[Nobservables + 5] = 1;
#endif
#ifdef MEASURE_Z_MAGNETIZATION
	valid_observable[Nobservables + 6] = 1;
#endif
#ifdef MEASURE_Z_MAGNETIZATION2
	valid_observable[Nobservables + 7] = 1;
#endif
}

double Metropolis(ExExFloat newWeight){
	return min(1.0,fabs((newWeight/currWeight).get_double()));
}

void update(){
	if(save_data_flag){ save_QMC_data(); exit(0); };
	int i,m,p,r,u,oldq,cont; double oldE, oldE2, v = Nop>0 ? val(rng) : 1; ExExFloat newWeight; double Rfactor;
	if(v < 0.8){ // composite update
		Rfactor = 1; oldq = q; memcpy(Sq_backup,Sq,q*sizeof(int)); memcpy(cycles_used_backup,cycles_used,Ncycles*sizeof(int));
		newWeight = currWeight;
		do{
			cont = 0; v = val(rng);
			if(v < 0.25){ // attempt to swap Sq[m] and Sq[m+1]
				if(q>=2){
					m = int(val(rng)*(q-1)); // m is between 0 and (q-2)
					if(Sq[m]!=Sq[m+1]){
						oldE = Energies[m+1]; old_currD = currD;
						p = Sq[m]; Sq[m] = Sq[m+1]; Sq[m+1] = p;
						GetEnergies();
						newWeight = UpdateWeightReplace(oldE,Energies[m+1]);
						if(REALW(currD) == 0) cont = 1;
					}
				}
			} else if(v < 0.5){ // attempt to delete Sq[m] and Sq[m+1]
				if(q>=2){
					m = int(val(rng)*(q-1)); // m is between 0 and (q-2)
					if(Sq[m]==Sq[m+1]){
						oldE = Energies[m]; oldE2 = Energies[m+1]; old_currD = currD;
						for(i=m;i<q-2;i++) Sq[i] = Sq[i+2]; q-=2;
						GetEnergies(); Rfactor /= Nop;
						newWeight = UpdateWeightDel(oldE,oldE2);
						if(REALW(currD) == 0) cont = 1;
					}
				}
			} else if(v < 0.75){
				if(q+2<qmax){ // attempt to insert Sq[m] and Sq[m+1]
					m = int(val(rng)*(q+1)); // m is between 0 and q
					old_currD = currD; p = diceNop(rng);
					for(i=q-1;i>=m;i--) Sq[i+2] = Sq[i]; q+=2; Sq[m] = Sq[m+1] = p;
					GetEnergies(); Rfactor *= Nop;
					newWeight = UpdateWeightIns(Energies[m],Energies[m+1]);
					if(REALW(currD) == 0) cont = 1;
				} else qmax_achieved = 1;
			} else{ // attempting a fundamental cycle completion
				int j = 0, inv_pr; double wfactor;
				u = geometric_int(rng); // a random integer u is picked according to geometric distribution
				if(q >= u+rmin){
					inv_pr = min(rmax,q-u)-(rmin)+1;
					r = int(val(rng)*inv_pr) + (rmin);  // r is random integer between rmin and min(rmax,q-u)
					PickSubsequence(r+u); // indexes of the subsequence are min_index, min_index+1,..., max_index=min_index+r+u-1
					std::shuffle(Sq_subseq,Sq_subseq+r+u,rng);
					if(NoRepetitionCheck(Sq_subseq,r)){
						for(i=0;i<u;i++) Sq_gaps[i] = Sq_subseq[i+r];
						m = FindCycles(r);
						if(found_cycles > 0){ // cycles[m] is one of the found cycles, containing all the operators of Sq_subseq
							P = cycles[m]; for(i=0;i<r;i++) P.reset(Sq_subseq[i]);
							p = P.count(); // here, p is length of the complement sequence S'
							if(q+p-r < qmax){
								if(r<p)	     for(i=q-1;i>max_index;i--) Sq[i+p-r] = Sq[i]; // shift the values to the right
								else if(r>p) for(i=max_index+1;i<q;i++) Sq[i+p-r] = Sq[i]; // shift the values to the left
								for(i=0;i<p;i++){ while(!P.test(j)) j++; Sq_subseq[i] = Sq[min_index+i] = j++;}
								for(i=0;i<u;i++) Sq[min_index+p+i] = Sq_gaps[i]; // S' contains the remaining operators
								std::shuffle(Sq+min_index,Sq+min_index+p+u,rng);
								q += p-r; // the length q may have changed
								newWeight = UpdateWeight();
								wfactor = found_cycles; FindCycles(p); wfactor /= found_cycles;
								wfactor *= factorial[p]/factorial[r];
								wfactor *= inv_pr; inv_pr = min(rmax,q-u)-(rmin)+1; wfactor /= inv_pr;
								Rfactor *= wfactor;
								cycles_used[m] = 1;
								if(REALW(currD) == 0) cont = 1;
							} else qmax_achieved = 1;
						}
					}
				}
			}
		} while(cont &&
#ifdef HURRY_ON_SIGTERM
			!save_data_flag &&
#endif
			val(rng) > COMPOSITE_UPDATE_BREAK_PROBABILITY);
		if(REALW(currD) != 0 && val(rng) < Metropolis(newWeight*Rfactor)){
			currWeight = newWeight;
		} else{
			q = oldq;
			memcpy(Sq,Sq_backup,q*sizeof(int));
			memcpy(cycles_used,cycles_used_backup,Ncycles*sizeof(int));
			currWeight = UpdateWeight();
		}
	} else if(v < 0.9 && q>=2){ // attempting a block swap
		m = q==2 ? 0 : int(val(rng)*(q-1)); // m is between 0 and (q-2)
		oldE = currEnergy; for(i=0;i<=m;i++) ApplyOperator(Sq[i]);
		for(i=0;i<=m;i++)  Sq_backup[q-1-m+i] = Sq[i];
		for(i=m+1;i<q;i++) Sq_backup[i-m-1] = Sq[i];
		for(i=0;i<q;i++) { p = Sq[i]; Sq[i] = Sq_backup[i]; Sq_backup[i] = p;}
		memcpy(Energies_backup,Energies,(q+1)*sizeof(double));
		GetEnergies(); newWeight = UpdateWeightReplace(oldE,currEnergy);
		if(val(rng) < Metropolis(newWeight)){
			z = lattice; currWeight = newWeight;
		} else{ UpdateWeightReplace(currEnergy,oldE); currEnergy = oldE;
			lattice = z; memcpy(Sq,Sq_backup,q*sizeof(int));
			memcpy(Energies,Energies_backup,(q+1)*sizeof(double));
		}
	} else{ // flip of a random spin
		p = diceN(rng);	lattice.flip(p); newWeight = UpdateWeight();
		if(val(rng) < Metropolis(newWeight)) { z = lattice; currWeight = newWeight;}
			else { lattice.flip(p); currWeight = UpdateWeight();}
	}
}

#ifdef MEASURE_CUSTOM_OBSERVABLES

std::complex<double> calc_MD0(int n){ // calculate <z | MD_0 | z> for the current configuration of spins and observable n
	std::complex<double> sum = 0;
	for(int i=0;i<MD0_size[n];i++) sum -= double(2*(int((MD0_product[n][i] & (~lattice)).count())%2)-1) * MD0_coeff[n][i];
	return sum;
}

std::complex<double> calc_MD(int n, int k){ // calculate d_k = <z | MD_k | z> for the current configuration of spins and observable n
	std::complex<double> sum = 0;
	for(int i=0;i<MD_size[n][k];i++) sum -= double(2*(int((MD_product[n][k][i] & (~lattice)).count())%2)-1) * MD_coeff[n][k][i];
	return sum;
}

#endif

double measure_H(){
	double R = d->z[q]/(-beta);
	if(q > 0) R += (d->divdiffs[q-1]/d->divdiffs[q]).get_double()*q/(-beta);
	return R;
}

double measure_H2(){
	double R = (d->z[q]/(-beta))*(d->z[q]/(-beta));
	if(q>0) R += (d->z[q]/(-beta) + d->z[q-1]/(-beta))*(d->divdiffs[q-1]/d->divdiffs[q]).get_double()*q/(-beta);
	if(q>1) R += (d->divdiffs[q-2]/d->divdiffs[q]).get_double()*(q*(q-1))/(-beta)/(-beta);
	return R;
}

double measure_Hdiag(){
	return currEnergy;
}

double measure_Hdiag2(){
	return currEnergy*currEnergy;
}

double measure_Hoffdiag(){
	double R = 0;
	if(q > 0) R += (d->divdiffs[q-1]/d->divdiffs[q]).get_double()*q/(-beta);
	return R;
}

double measure_Hoffdiag2(){
	double R = (d->z[q]/(-beta))*(d->z[q]/(-beta)) + currEnergy*(currEnergy - 2*measure_H());
	if(q>0) R += (d->z[q]/(-beta) + d->z[q-1]/(-beta))*(d->divdiffs[q-1]/d->divdiffs[q]).get_double()*q/(-beta);
	if(q>1) R += (d->divdiffs[q-2]/d->divdiffs[q]).get_double()*(q*(q-1))/(-beta)/(-beta);
	return R;
}

double measure_Z_magnetization(){
	return (2.0*lattice.count() - N)/N; // this value is between -1 and 1
}

double measure_Z_magnetization2(){
	double m = measure_Z_magnetization();
	return m*m;
}

std::string name_of_observable(int n){
	std::string s;
	if(n < Nobservables){
#ifdef MEASURE_CUSTOM_OBSERVABLES
		s = Mnames[n];
#endif
	} else switch(n-Nobservables){
			case 0: s = "H";             break;
			case 1: s = "H^2";           break;
			case 2: s = "H_{diag}";      break;
			case 3: s = "H_{diag}^2";    break;
			case 4: s = "H_{offdiag}";   break;
			case 5: s = "H_{offdiag}^2"; break;
			case 6: s = "Z_magnetization"; break;
			case 7: s = "Z_magnetization^2"; break;
	}
	return s;
}

double measure_observable(int n){
	double R = 0;
	if(valid_observable[n]){
		if(n < Nobservables){
#ifdef MEASURE_CUSTOM_OBSERVABLES
			int i,k,len,cont;
			std::complex<double> T = calc_MD0(n);
			for(k=0;k<MNop[n];k++){
				P = MP[n][k]; len = P.count(); if(len>q) continue;
				if(!NoRepetitionCheck(Sq+(q-len),len)) continue;
				cont = 0; for(i=0;i<len;i++) if(!P.test(Sq[q-1-i])){ cont = 1; break;} if(cont) continue;
				T +=	(d->divdiffs[q-len]/d->divdiffs[q]).get_double() *
				        (beta_pow_factorial[q-len]/beta_pow_factorial[q]).get_double()/factorial[len] *
					(currD_partial[q-len]/currD) * calc_MD(n,k);
			}
#ifdef ABS_WEIGHTS
			R = std::abs(T)*currWeight.sgn()*cos(std::arg(T)+std::arg(currD));
#else
			R = std::real(currD*T)/std::real(currD); // we importance-sample Re(W_C A_C)/Re(W_C)
#endif
#endif
		} else  switch(n-Nobservables){
				case 0:	R = measure_H(); break;
				case 1:	R = measure_H2(); break;
				case 2:	R = measure_Hdiag(); break;
				case 3:	R = measure_Hdiag2(); break;
				case 4:	R = measure_Hoffdiag(); break;
				case 5:	R = measure_Hoffdiag2(); break;
				case 6: R = measure_Z_magnetization(); break;
				case 7: R = measure_Z_magnetization2(); break;
		}
	}
	return R;
}

void measure(){
	double R, sgn; int i;
	currWeight = GetWeight();
#ifdef ABS_WEIGHTS
	sgn = currWeight.sgn() * cos(std::arg(currD)); // arg(W) = arg(currD) + arg(currWeight), where arg(currWeight) = either 0 or Pi
#else
	sgn = currWeight.sgn();
#endif
	meanq += q; if(maxq < q) maxq = q; in_bin_sum_sgn += sgn;
	if((measurement_step+1) % bin_length == 0){
		in_bin_sum_sgn /= bin_length; bin_mean_sgn[measurement_step/bin_length] = in_bin_sum_sgn; in_bin_sum_sgn = 0;
	}
	for(i=0;i<N_all_observables;i++){
		R = measure_observable(i); in_bin_sum[i] += R*sgn;
		if((measurement_step+1) % bin_length == 0){
			in_bin_sum[i] /= bin_length; bin_mean[i][measurement_step/bin_length] = in_bin_sum[i]; in_bin_sum[i] = 0;
		}
	}
}

double mean_O[N_all_observables], stdev_O[N_all_observables], mean_O_backup[N_all_observables];

const int N_derived_observables = 2;  // we define number of derived observables

std::string name_of_derived_observable(int n){ // we define names of derived observables
	std::string s;
	switch(n){
		case 0 : s = "specific heat"; break;
		case 1 : s = "magnetic susceptibility"; break;
	}
	return s;
}

int valid_derived_observable(int n){ // we define which observables are needed for each derived observable
	int r = 0;
	switch(n){
		case 0: r = valid_observable[Nobservables+0] && valid_observable[Nobservables+1]; break;
		case 1: r = valid_observable[Nobservables+6] && valid_observable[Nobservables+7]; break;
	}
	return r;
}

double compute_derived_observable(int n){ // we compute the derived observables
	double R = 0;
	switch(n){
		case 0 : R = beta*beta*(mean_O[Nobservables+1] - mean_O[Nobservables+0]*mean_O[Nobservables+0]); break;
		case 1 : R = beta*(mean_O[Nobservables+7] - mean_O[Nobservables+6]*mean_O[Nobservables+6]); break;
	}
	return R;
}
