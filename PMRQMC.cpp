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

#include<iostream>
#include<iomanip>
#include<complex>
#include<random>
#include<cstdlib>
#include<algorithm>
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
static std::geometric_distribution<> geometric_int(0.8);

ExExFloat beta_pow_factorial[qmax]; // contains the values (-beta)^q / q!
double factorial[qmax]; // contains the values q!
int cycle_len[Ncycles];
int cycles_used[Ncycles];
int n_cycles[Nop+3]; // numbers of cycles of lengths 0,1,2,...,Nop+2, the last three values are always zeros.
int cycle_min_len, cycle_max_len, found_cycles, min_index, max_index;

#ifndef MEASURE_CUSTOM_OBSERVABLES
#define Nobservables 0
#endif

const int N_all_observables = Nobservables + 6;
int valid_observable[N_all_observables];

int bin_length = measurements / Nbins;
double in_bin_sum[N_all_observables];
double bin_mean[N_all_observables][Nbins];
double in_bin_sum_sgn;
double bin_mean_sgn[Nbins];

int q;
int qmax_achieved=0;

divdiff* d;

std::bitset<N> lattice;
std::bitset<N> z;
std::bitset<Nop> P;

int Sq[qmax];	// list of q operator indices
int Sq_backup[qmax];
int Sq_subseq[qmax];
int Sq_gaps[qmax];
double Energies[qmax+1];
double Energies_backup[qmax+1];
int EnergyIndexes[qmax];
int eoccupied[qmax];
double currEnergy;
std::complex<double> old_currD, currD;
std::complex<double> currD_partial[qmax];

ExExFloat one, currWeight;
unsigned long long step;
unsigned long long measurement_step;

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
	return d->divdiffs[q] * beta_pow_factorial[q] * currD.real();
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
        return d->divdiffs[q] * beta_pow_factorial[q] * currD.real();
}

ExExFloat UpdateWeightReplace(double removeEnergy, double addEnergy){
	if(removeEnergy != addEnergy) if(d->RemoveValue(-beta*removeEnergy)) d->AddElement(-beta*addEnergy); else{
		std::cout << "Error: energy not found" << std::endl; exit(1);
	}
	return d->divdiffs[q] * beta_pow_factorial[q] * currD.real(); // use this value only when the values of q and currD are correct
}

ExExFloat UpdateWeightDel(double removeEnergy1, double removeEnergy2){
	if(d->RemoveValue(-beta*removeEnergy1) && d->RemoveValue(-beta*removeEnergy2))
		return d->divdiffs[q] * beta_pow_factorial[q] * currD.real();  // use this value only when the values of q and currD are correct
	else{
		std::cout << "Error: energy not found" << std::endl; exit(1);
	}
}

ExExFloat UpdateWeightIns(double addEnergy1, double addEnergy2){
	d->AddElement(-beta*addEnergy1); d->AddElement(-beta*addEnergy2);
	return d->divdiffs[q] * beta_pow_factorial[q] * currD.real();    // use this value only when the values of q and currD are correct
}

int NoRepetitionCheck(int* sequence, int r){ // check for absence of repetitions in a sequence of length r
	int i,j,rep = 1;
	for(i=0;i<r && rep;i++) for(j=0;j<i;j++) if(sequence[j]==sequence[i]){ rep = 0; break;}
	return rep;
}

void PickSubsequence(int r){ // randomly picks a sequential sub-sequence of length r from Sq
	int i,j,m; m = int(val(rng)*(q-r+1)); // m is random integer between 0 and q-r
	for(i=0;i<r;i++) Sq_subseq[i] = Sq[i+m];
	min_index = m; max_index = m+r-1;
}

int FindCycles(int r){  // find all cycles of length between lmin and lmax, each containing all operators of the array Sq_subseq of length r.
	int i,j,not_contained;
	int found_cycle_list[Ncycles];
	found_cycles = 0;                                  // found_cycles contains the number of cycles found. it is global variable.
	for(i=0;i<Ncycles;i++){
		if(cycle_len[i]<lmin || cycle_len[i]>lmax) continue;
		not_contained = 0; for(j=0;j<r;j++) if(!cycles[i].test(Sq_subseq[j])){ not_contained = 1; break;}
		if(not_contained) continue;
		found_cycle_list[found_cycles++] = i;
	}
	return found_cycles>0 ? found_cycle_list[int(val(rng)*found_cycles)] : -1; // returns one of the found cycles chosen randomly.
}

void init(){
	int i,j,k,l; double curr2=1; ExExFloat curr1; beta_pow_factorial[0] = curr1; factorial[0] = curr2;
	for(q=1;q<qmax;q++){ curr1*=(-beta)/q; curr2*=q; beta_pow_factorial[q] = curr1; factorial[q] = curr2;}
	unsigned int rng_seed = rd(); rng.seed(rng_seed); std::cout << "RNG seed = " << rng_seed << std::endl;
	lattice = 0; for(i=N-1;i>=0;i--) if(dice2(rng)) lattice.set(i); z = lattice; q=0;
	currWeight = GetWeight();
	for(i=0;i<Ncycles;i++) cycle_len[i] = cycles[i].count();
	cycle_min_len = 64; for(i=0;i<Ncycles;i++) cycle_min_len = min(cycle_min_len,cycle_len[i]);
	cycle_max_len = 0; for(i=0;i<Ncycles;i++) cycle_max_len = max(cycle_max_len,cycle_len[i]);
	for(i=0;i<Ncycles;i++) cycles_used[i] = 0;
	for(i=0;i<Nop+3;i++) n_cycles[i] = 0; for(i=0;i<Ncycles;i++) n_cycles[cycle_len[i]]++;
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
}

double Metropolis(ExExFloat newWeight){
	return min(1.0,fabs((newWeight/currWeight).get_double()));
}

void update(){
	int i,m,p,r,u; double oldE, oldE2, v = Nop>0 ? val(rng) : 1; ExExFloat newWeight;
	if(v < 0.3){ //local move
		if(v<0.1 && q>=2){ // attempt to swap Sq[m] and Sq[m+1]
			m = int(val(rng)*(q-1)); // m is between 0 and (q-2)
			if(Sq[m]!=Sq[m+1]){
				oldE = Energies[m+1]; old_currD = currD;
				p = Sq[m]; Sq[m] = Sq[m+1]; Sq[m+1] = p; GetEnergies(); newWeight = UpdateWeightReplace(oldE,Energies[m+1]);
				if(val(rng) < Metropolis(newWeight)) currWeight = newWeight; else {
					Sq[m+1] = Sq[m]; Sq[m] = p; currD = old_currD;
					UpdateWeightReplace(Energies[m+1],oldE); Energies[m+1] = oldE;
				}
			}
		} else if(v<0.2){ // attempt to delete Sq[m] and Sq[m+1]
			if(q>=2){
				m = int(val(rng)*(q-1)); // m is between 0 and (q-2)
				if(Sq[m]==Sq[m+1]){
					oldE = Energies[m]; oldE2 = Energies[m+1]; old_currD = currD;
					memcpy(Sq_backup,Sq,q*sizeof(int)); memcpy(Energies_backup,Energies,(q+1)*sizeof(double));
					for(i=m;i<q-2;i++) Sq[i] = Sq[i+2]; q-=2;
					GetEnergies(); newWeight = UpdateWeightDel(oldE,oldE2);
					if(val(rng) < Metropolis(newWeight)/Nop) currWeight = newWeight; else{
						q+=2; memcpy(Sq,Sq_backup,q*sizeof(int)); memcpy(Energies,Energies_backup,(q+1)*sizeof(double));
						currD = old_currD; UpdateWeightIns(oldE,oldE2);
					}
				}
			}
		} else if(q+2<qmax){ // attempt to insert Sq[m] and Sq[m+1]
			m = int(val(rng)*(q+1)); // m is between 0 and q
			memcpy(Sq_backup,Sq,q*sizeof(int)); memcpy(Energies_backup,Energies,(q+1)*sizeof(double));
			old_currD = currD; p = diceNop(rng);
			for(i=q-1;i>=m;i--) Sq[i+2] = Sq[i]; q+=2; Sq[m] = Sq[m+1] = p;
			GetEnergies(); newWeight = UpdateWeightIns(Energies[m],Energies[m+1]);
			if(val(rng) < Metropolis(newWeight)) currWeight = newWeight; else{
				q-=2; memcpy(Sq,Sq_backup,q*sizeof(int)); memcpy(Energies,Energies_backup,(q+1)*sizeof(double));
				currD = old_currD; d->RemoveElement(); d->RemoveElement();
			}
		} else qmax_achieved = 1;
	} else if(v < 0.5){ // attempting a fundamental cycle completion
		int oldq, j = 0, s, t, inv_pr; double wfactor;
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
					memcpy(Sq_backup,Sq,q*sizeof(int));
					oldq = q; P = cycles[m]; for(i=0;i<r;i++) P.reset(Sq_subseq[i]);
					p = P.count(); // here, p is length of the complement sequence S'
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
					if(q < qmax && val(rng) < Metropolis(newWeight*wfactor)){
						currWeight = newWeight; cycles_used[m] = 1;
					} else{
						if(q>=qmax) qmax_achieved = 1;
						q = oldq; memcpy(Sq,Sq_backup,q*sizeof(int)); currWeight = UpdateWeight();
					}
				}
			}
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

double meanq = 0;
double maxq = 0;

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
	}
	return s;
}

double measure_observable(int n){
	double R = 0; int i,k,len,cont;
	if(valid_observable[n]) if(n < Nobservables){
#ifdef MEASURE_CUSTOM_OBSERVABLES
		std::complex<double> T = calc_MD0(n);
		for(k=0;k<MNop[n];k++){
			P = MP[n][k]; len = P.count(); if(len>q) continue;
			if(!NoRepetitionCheck(Sq+(q-len),len)) continue;
			cont = 0; for(i=0;i<len;i++) if(!P.test(Sq[q-1-i])){ cont = 1; break;} if(cont) continue;
			T +=	(d->divdiffs[q-len]/d->divdiffs[q]).get_double() *
			        (beta_pow_factorial[q-len]/beta_pow_factorial[q]).get_double()/factorial[len] *
				(currD_partial[q-len]/currD) * calc_MD(n,k);
		}
		R = (currD*T).real()/currD.real(); // we importance-sample Re(W_C A_C)/Re(W_C)
#endif
	} else  switch(n-Nobservables){
			case 0:	R = measure_H(); break;
			case 1:	R = measure_H2(); break;
			case 2:	R = measure_Hdiag(); break;
			case 3:	R = measure_Hdiag2(); break;
			case 4:	R = measure_Hoffdiag(); break;
			case 5:	R = measure_Hoffdiag2(); break;
	}
	return R;
}

void measure(){
	double R, sgn; int i;
	currWeight = GetWeight(); sgn = currWeight.sgn();
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

double get_cpu_time(){ return (double)clock() / CLOCKS_PER_SEC;}

int main(int argc, char* argv[]){
	double start_time = get_cpu_time();
	double Rsum[N_all_observables] = {0}; double sgn_sum = 0;
	double over_bins_sum[N_all_observables] = {0}; double over_bins_sum_sgn = 0;
	double over_bins_sum_cov[N_all_observables] = {0}; double mean[N_all_observables]; double stdev[N_all_observables];
	int i,k,o=0; divdiff_init(); divdiff dd(q+4,500); d=&dd; init();
	std::cout << "Parameters: beta = " << beta << ", Tsteps = " << Tsteps << ", steps = " << steps << std::endl;
	for(step=0;step<Tsteps;step++) update();
	for(measurement_step=0;measurement_step<measurements;measurement_step++){
		for(step=0;step<stepsPerMeasurement;step++) update(); measure();
	}
	for(i=0;i<Nbins;i++) sgn_sum += bin_mean_sgn[i]; sgn_sum /= Nbins;
	for(i=0;i<Nbins;i++) over_bins_sum_sgn += (bin_mean_sgn[i] - sgn_sum)*(bin_mean_sgn[i] - sgn_sum); over_bins_sum_sgn /= (Nbins*(Nbins-1));
	std::cout << std::setprecision(9);
	std::cout << "mean(sgn(W)) = " << sgn_sum << std::endl;
	std::cout << "std.dev.(sgn(W)) = " << sqrt(over_bins_sum_sgn) << std::endl;
	if(qmax_achieved) std::cout << "Warning: qmax = " << qmax << " was achieved" << std::endl;
	for(i=0;i<Ncycles;i++) if(!cycles_used[i]) std::cout << "Warning: cycle No. " << i << " of length " << cycle_len[i] << " was not used" << std::endl;
	std::cout << "mean(q) = " << meanq / measurements << std::endl;
	std::cout << "max(q) = "<< maxq << std::endl;
	for(k=0;k<N_all_observables;k++) if(valid_observable[k]){
		std::cout << "Observable #" << ++o << ": "<< name_of_observable(k) << std::endl;
		for(i=0;i<Nbins;i++) Rsum[k] += bin_mean[k][i]; Rsum[k] /= Nbins;
		for(i=0;i<Nbins;i++) over_bins_sum[k] += (bin_mean[k][i] - Rsum[k])*(bin_mean[k][i] - Rsum[k]); over_bins_sum[k] /= (Nbins*(Nbins-1));
		for(i=0;i<Nbins;i++) over_bins_sum_cov[k] += (bin_mean[k][i] - Rsum[k])*(bin_mean_sgn[i] - sgn_sum); over_bins_sum_cov[k] /= (Nbins*(Nbins-1));
		// std::cout << "mean(O*sgn(W)) = " << Rsum[k] << std::endl;
		// std::cout << "std.dev.(O*sgn(W)) = " << sqrt(over_bins_sum[k]) << std::endl;
		mean[k] = Rsum[k]/sgn_sum*(1 + over_bins_sum_sgn/sgn_sum/sgn_sum) - over_bins_sum_cov[k]/sgn_sum/sgn_sum;
		stdev[k] = fabs(Rsum[k]/sgn_sum)*sqrt(over_bins_sum[k]/Rsum[k]/Rsum[k] + over_bins_sum_sgn/sgn_sum/sgn_sum - 2*over_bins_sum_cov[k]/Rsum[k]/sgn_sum);
		std::cout << "mean(O) = " << mean[k] << std::endl;
		std::cout << "std.dev.(O) = " << stdev[k] << std::endl;
	}
	divdiff_clear_up();
	std::cout << "wall-clock cpu time = " << get_cpu_time()-start_time << " seconds" << std::endl;
	return 0;
}
