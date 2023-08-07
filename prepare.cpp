//
// This is an auxiliary program for Permutation Matrix Representation Quantum Monte Carlo for arbitrary spin-1/2 Hamiltonians.
//
// This program is introduced in the paper:
// Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians (2023).
//
// This program is licensed under a Creative Commons Attribution 4.0 International License:
// http://creativecommons.org/licenses/by/4.0/
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace std;

#define dbitset vector<char>

int no_qubit = 0; // number of qubits is a global variable!

// *****************************    Functions ******************************************* //
template<typename T>
void printMatrix(const vector<vector<T>>& matrix) {
    int m = matrix.size();
    int n = matrix[0].size();

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// This function is for binary addition (XOR)
vector<int> GF2_add(vector<int> vec1 , vector<int> vec2){
    vector<int> result;
    if (vec1.size() != vec2.size()){
        cerr << "The size of the added vectors are differen!" << endl;
    }

    for(int i = 0; i < vec1.size(); i++){
        result.push_back((vec1[i] + vec2[i])%2);
    }
    return result;
}

// Function to compute the mod 2 nullspace
vector<vector<int>> Null2(const vector<vector<int>>& matrix) {
    vector<vector<int>> nullspaceBasis;
    int numRows = matrix.size(); int numCols;    
    if(numRows==0) return nullspaceBasis; else numCols = matrix[0].size();
    vector<vector<int>> matrix_RE = matrix; // Copy the original matrix
    vector<int> nullspaceVector(numCols, 0);
    vector<int> marked_rows;

    // Gaussian Elimination (GF(2)) to obtain matrix_GE
    for(int j = 0; j < numCols; j++){
        for(int i=0; i < numRows; i++){
            if(matrix_RE[i][j] == 1){
                marked_rows.push_back(i);
                for(int k = 0; k < numCols; k++ ){
                    if(k != j){
                        if(matrix_RE[i][k] == 1){
                            // Adding column j to k!
                            for(int l = 0; l < matrix_RE.size(); l++){
                                matrix_RE[l][k] = matrix_RE[l][j] ^ matrix_RE[l][k];
                            }
                        }
                    }
                }
                break;
            }
        }
    }
    // sort the marked row
    std::sort(marked_rows.begin(), marked_rows.end());

    // Remove consecutive duplicates (taking only unique marked rows!)
    auto it = unique(marked_rows.begin(), marked_rows.end());
    marked_rows.erase(it, marked_rows.end());

    vector<vector<int>> marked_matrix;
    // Make the marked matrix:
    for(int i = 0 ; i < marked_rows.size() ; i++){
        marked_matrix.push_back(matrix_RE[marked_rows[i]]);
    }

    for(int i = 0; i < numRows; i++){
        auto iter = find(marked_rows.begin(), marked_rows.end(), i); 
        if(iter == marked_rows.end()){ // Finding unmarked (independent) rows
            // Step 1 - Then find the ones in that row!
            // 2 - Go over each one found in that row and find the corresponding
            // row which has that that index of one in that column

            // Step 1 : 
            vector<int> row_i = matrix_RE[i] , nullvector_j(numRows, 0);
            nullvector_j[i] = 1;
            for(int j = 0; j < numCols; j++){
                if(row_i[j] == 1){
                    // Go over the marked rows and extract the null space:
                    for(int k = 0; k < marked_matrix.size(); k++){
                        if(marked_matrix[k][j] == 1){
                            nullvector_j[marked_rows[k]] = 1;
                        } 
                    }
                }
            }
            nullspaceBasis.push_back(nullvector_j);
        }
    }
    return nullspaceBasis;
}

int Int_sum(vector<int> vec){
    int sum = 0;
    for(int i = 0; i < vec.size(); i++){
        sum += vec[i];
    }
    return sum;
}

// This function minimizes the size of the cycles (nullspace basis vectors with least 1s):
bool compare_null(const vector<int> & a, const vector<int> & b){
	return Int_sum(a) < Int_sum(b);
}

int cycle_minimize(vector<vector<int>>& null_eigs){ 
	int nullsize, k, m, null_k, changes_made = 0; vector<int> curr;
	nullsize = null_eigs.size();
	sort(null_eigs.begin(), null_eigs.end(), compare_null);
	for(k = nullsize-1; k > 0 ; k--){
		null_k = Int_sum(null_eigs[k]);
		for(m = 0 ; m < k; m++){
			curr = GF2_add( null_eigs[k] , null_eigs[m]);
			if(Int_sum(curr) < null_k){
				null_eigs[k] = curr; changes_made = 1;
				break;
			}
		}
	}
	return changes_made;
}

// This function minimizes the size of the cycles (nullspace basis vectors with least 1s):
// This Eig_minimize is specifically designed for operator permutations, and functions differently
//      than the one in Null_Generator, which is used to find closed cycles of minimum length 3.
int Eig_minimize(vector<vector<int>>& nullEigs){
    int changes_made = 0;
    int nullsize = nullEigs.size() , noops = nullEigs[0].size(); //noops is the number of operators

    vector<int> nulln(nullsize);
    for(int i = 0; i < nullsize; i++){
        nulln[i] = Int_sum(nullEigs[i]);
    }
    auto mincyc = min_element(nulln.begin(), nulln.end());

    vector<int> nullEigsMinind , nullEigsHighind, nullEigsOnesind; // the indices of eigenvectors ending with 1
    
    for(int i = 0; i < nullsize; i++){
        if(nulln[i] == *mincyc){
            nullEigsMinind.push_back(i);
        }
        else if(nullEigs[i][noops-1] == 1){
            nullEigsOnesind.push_back(i);
        }
        else{
            nullEigsHighind.push_back(i);
        }
    }
    int highsize = nullEigsHighind.size();
    int minsize = nullEigsMinind.size();
    
    for(int i=0; i < minsize; i++){
        int vecisum = Int_sum(nullEigs[nullEigsMinind[i]]);
        vector<int> bestVec = nullEigs[nullEigsMinind[i]]; 
        for(int j=0; j < minsize; j++){
            if(j!= i){
                vector<int> tryVecij = GF2_add(nullEigs[nullEigsMinind[i]] , nullEigs[nullEigsMinind[j]]);
                int tryvecijsum = Int_sum(tryVecij);
                if(tryvecijsum < vecisum){
                    vecisum = tryvecijsum;
                    bestVec = tryVecij;
                }
            }
        }
	if(nullEigs[nullEigsMinind[i]] != bestVec){ nullEigs[nullEigsMinind[i]] = bestVec; changes_made = 1;}
    }

    for(int k = 0; k < nullEigsOnesind.size(); k++){
        int oneskind = nullEigsOnesind[k];
        vector<int> highEigk = nullEigs[oneskind];
        vector<int> bestEigk = highEigk;
        int bestnullk = Int_sum(highEigk);
        for(int l = 0; l < nullsize; l++){
            if(l != oneskind){
                vector<int> highTomin = GF2_add(highEigk , nullEigs[l]);
                int lowk = Int_sum(highTomin);
                if(lowk < 3){
                    bestEigk = highTomin;
                    bestnullk = lowk;
                    break;
                }
                else if(lowk < bestnullk){
                    bestEigk = highTomin;
                    bestnullk = lowk;
                }
            }
        }
        if(nullEigs[oneskind] != bestEigk){ nullEigs[oneskind] = bestEigk; changes_made = 1;}
    }
    return changes_made;
}

// This function converts the array of integers into a corresponding binary string (used to convert the indices of Z into string of bitsets)
string int_to_str(vector<int> Z){
    string Z_string = "";
    if(Z.size() < 1){
        return "0";
    }
    std::sort(Z.begin() , Z.end());
    int Z_size = Z.size();
    int count = 1, ind_z_count = 0, max_z = Z[Z_size-1];

    while(ind_z_count < Z_size){
	if(Z[ind_z_count] < count){
		cout << endl << "Error: repeating spin indices are detected in a Pauli string" << endl; exit(1);
	} else if(Z[ind_z_count] == count){
            Z_string = "1" + Z_string;
            ind_z_count++;
        }
        else{
            Z_string = "0" + Z_string;
        }
        count++;
    }
    return Z_string;
}

// This function extracts the input file information into a vector of pairs!
vector<pair<complex<double>, vector<int>>> data_extract(const string& fileName){
    vector<pair<complex<double>, vector<int>>> data;

    ifstream inputFile(fileName);
    if (!inputFile) {
        cout << "Failed to open the input file!" << endl;
        return data;
    }

    string line;
    while (getline(inputFile, line)) {
        // The first non-empty element is the coefficient
        istringstream iss(line);
        pair<complex<double>,vector<int>> linedata;

        // Extracting the complex coefficient:
        double realpart, imagpart=0;
        char sign;
        string complexPart;
        iss >> complexPart; 

        istringstream complexIss(complexPart);
        complexIss >> realpart >> imagpart;
        linedata.first = complex<double> (realpart , imagpart);

        //Extracting the integer vectors of qubits and paulis:
        string token;
        vector<int> integers;
        while (iss >> token){
            // int num = std::stoi(token);
            int num = token == "X" || token == "x" ? 1 : 
                      token == "Y" || token == "y" ? 2 : 
                      token == "Z" || token == "z" ? 3 : std::stoi(token);
            integers.push_back(num);
        }
        linedata.second = integers;

        data.push_back(linedata);
        }
        inputFile.close();
    return data;
}

// Finding bitset in a vector of bitsets:
// The max bitset will be 64 (the maximum number of qubits), and we will cut off accordingly
// ..... when the number of qubits of the system is less than this number.
pair<bool , int> bit_is_in_set(dbitset bitstring, vector<dbitset> bitsetVector) {
    bool found = false;
    int found_indx = 0, ind_count = 0;
    pair<bool, int> output;
    for (const auto& bitset : bitsetVector) {
            if (bitset == bitstring) {
                found = true;
                found_indx = ind_count;
                break;
            }
            ind_count++;
        }
    output.first = found;
    output.second = found_indx;
    return output;
}

// This function downsizes the vector of bitsets from dbitset to bitset<no_qubit> to avoid 
// ..... having redundant zeros.
vector<vector<bool>> downsize_bitset(vector<dbitset> bitsetVector){
    extern int no_qubit;
    vector<vector<bool>> bitset_down;
    for(const auto& bitset : bitsetVector){
        vector<bool> bitset_vector;
        for(int i = 0; i < no_qubit; i++){
            bitset_vector.push_back(i < bitset.size() && bitset[i]);
        }
        reverse(bitset_vector.begin() , bitset_vector.end())
;        bitset_down.push_back(bitset_vector);
    }
    return bitset_down;
}

vector<vector<int>> bit_to_intvec(vector<vector<bool>> bitsetVec){
    extern int no_qubit;
    vector<vector<int>> bit_int;
    for(auto const& bitset : bitsetVec){
        vector<int> bit_int_i;
        for(int i=0; i < no_qubit; i++){
            bit_int_i.push_back(i < bitset.size() && bitset[i]);
        }
        reverse(bit_int_i.begin() , bit_int_i.end());
        bit_int.push_back(bit_int_i);
    }

    return bit_int;
}

// ------------- Functions to compute the Zs and Ps and coefficients ----------
struct BitsetComparator {
    bool operator()(const dbitset& lhs, const dbitset& rhs) const {
        int i = lhs.size()-1, j = rhs.size()-1, ans = 0;
        while(i >= 0 && lhs[i] == 0) i--;
        while(j >= 0 && rhs[j] == 0) j--;
        if(i < j) ans = 1; else if(i == j){
            while(i >= 0 && lhs[i] == rhs[i]) i--;
            if(i >= 0 && lhs[i] < rhs[i]) ans = 1;
	}
	return ans;
        // return lhs.to_ulong() < rhs.to_ulong();
    }
}; // This struct is for bitset comparison (in order to sort the bitsets)

typedef vector<complex<double>> Coeffs;
typedef vector<vector<int>> ZVecs;
struct PZdata {
    vector<dbitset> Ps;
    Coeffs coeffs;
    ZVecs Zs;
    ZVecs Z_track;
};

PZdata PZcomp(const vector<pair<complex<double>,vector<int>>>& data) {
    PZdata PZ_data;
    int l = data.size(),z_count = 0;
    extern int no_qubit;
    vector<dbitset> Ps;
    Coeffs coeffs;
    ZVecs Zs;
    ZVecs Z_track; //This vector maps the Zs to Ps: it is a many to one mapping!

    for (int i = 0; i<l; i++){
        complex<double> coeff_i = data[i].first;
        vector<int> zs_i; // For every line zs extracts the qubits on which a pauli Z acts! 
        vector<int> data_i = data[i].second; // Extracts the array of qubits and paulis for every line of input!
        dbitset bit_num; // This variable keeps track of the index of the permutation matrix we get for each line of data!

        for (size_t j = 0; j < data_i.size() / 2; j++) {
            // Format of the input file: The 1st, 3rd, 5th, ... indicate the qubits
            int qubit = data_i[2 * j];
            // Format of the input file: The 2nd, 4th, 6th, ... indicate the paulis
            int pauli_j = data_i[2 * j + 1];

            if (qubit <= 0){
                cout << endl << "Error: " << qubit << " is incorrect spin index. Spin indices must be positive integers." << endl; exit(1);
            }
            if (qubit > no_qubit) no_qubit = qubit;
            if (pauli_j == 1) {
		if(qubit > bit_num.size()) bit_num.resize(qubit);
                bit_num[qubit-1] = 1;
                // bit_num.set(qubit-1, true);
                // Add a 1 in the num-th position from the right of the bit string 
            } else if (pauli_j == 2) {
		if(qubit > bit_num.size()) bit_num.resize(qubit);
                //cpp_int perm_term = 1 << (qubit - 1);
                //num += perm_term;
                bit_num[qubit-1] = 1;
                // bit_num.set(qubit-1, true);
                coeff_i *= complex<double>(0, 1);
                zs_i.push_back(qubit);
            } else if (pauli_j == 3) {
                zs_i.push_back(qubit);
            }
        }
        coeffs.push_back(coeff_i);
        Zs.push_back(zs_i); // If zs is empty then no non-trivial diagonal components! (All identity operators)
            
        // Look for num in the previous list of Ps (permutations)
        pair<bool , int> bit_in_set = bit_is_in_set(bit_num , Ps);
        if(!bit_in_set.first){ 
            // num was not found in Ps, thus a new permutation matrix!
            Ps.push_back(bit_num);
            // We will add i-th element of coeffs and Zs to be associated with the current Ps!
            Z_track.push_back(vector<int> {i});
        }else {
            // num was found in Ps, and P_index will be the index of Ps that matched num.
            int P_index = bit_in_set.second;

            vector<int> z_indices = Z_track[P_index];

            bool z_found = false;
            for (int k = 0; k < z_indices.size(); k++){
                if (Zs[z_indices[k]] == zs_i){
                    coeffs[P_index] += coeff_i;
                    z_found = true;
                    break;
                }
            }
            // If the z array is new, we add it to the Z_track for the Ps associated with num! 
            if (!z_found){
                Z_track[P_index].push_back(i);
            }
        }
    }

    // Throw away the zero coefficients:
    vector<dbitset> Ps_kept;
    ZVecs Z_track_kept;

    for(int k = 0; k < Z_track.size(); k++){
        vector<int> ztrack_k;
        for(int l = 0; l < Z_track[k].size(); l++){
            int k_l_index = Z_track[k][l];
            complex<double> coeff_k_l = coeffs[k_l_index];
            if (abs(coeff_k_l) > 1e-8){
                //zero_coeffs_for_k.push_back(l);
                ztrack_k.push_back(k_l_index);
            }
        }
        if(ztrack_k.size() > 0){
            Z_track_kept.push_back(ztrack_k);
            Ps_kept.push_back(Ps[k]);
        }
    }

    // Sorting everything based on indices of Ps:
    vector<int> indices;
    for(int i = 0; i < Ps.size(); i++){
        indices.push_back(i);
    }
    sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        return BitsetComparator()(Ps_kept[a], Ps_kept[b]);
    });


    vector<dbitset> Ps_sorted;
    ZVecs Z_track_sorted;
    for(int i = 0; i < Ps_kept.size(); i++){
        Ps_sorted.push_back(Ps_kept[indices[i]]);
        Z_track_sorted.push_back(Z_track_kept[indices[i]]);
    }

    PZ_data.coeffs = coeffs; //Coeffs and Zs are kept as the original data
    PZ_data.Ps = Ps_sorted;
    PZ_data.Zs = Zs;
    PZ_data.Z_track = Z_track_sorted;

    return PZ_data;
}

int none(dbitset v){
    return std::all_of(v.begin(), v.end(), [](int i) { return i==0; });
}

ofstream output;

void output_complex(complex<double> z){
    if(z.imag() == 0) output << z.real(); else output << "{" << z.real() << "," << z.imag() << "}";
}

void main1(int argc , char* argv[]){
    string fileName(argv[1]);  // Reading the name of the input .txt file describing the Hamiltonian
    vector<pair<complex<double>, vector<int>>> data = data_extract(fileName);

    // Unpacking the data from the input file "fileName"
    PZdata PZ_data = PZcomp(data);
    vector<dbitset> Ps_bit = PZ_data.Ps;
    vector<vector<bool>> Ps = downsize_bitset(Ps_bit);
    vector<complex<double>> coefficients = PZ_data.coeffs;
    vector<vector<int>> Z_track = PZ_data.Z_track;
    vector<vector<int>> Zs = PZ_data.Zs;
    vector<string> Zs_string;
    bool D0_exists = false;
    //int no_qubit = PZ_data.no_qubit;

    // Converting the Z indices into string of bitsets!
    for(int i = 0; i < Zs.size();i++){
        Zs_string.push_back(int_to_str(Zs[i]));
    }

    // Convert Ps indices into vector of ints
    vector<vector<bool>> Ps_nontrivial = Ps;
    if(none(Ps_bit[0])){
        Ps_nontrivial.erase(Ps_nontrivial.begin());
        D0_exists = true;
    }

    // Minimizing the size of the fundamental cycless
    vector<vector<int>> Ps_binary = bit_to_intvec(Ps_nontrivial);
    vector<vector<int>> nullspace = Null2(Ps_binary);
    while(cycle_minimize(nullspace));
    int no_ps = Ps_binary.size();
    int nullity = nullspace.size();

    //string output_h = fileName.substr(0, fileName.find_last_of(".")) + ".h";
    // string output_h = "hamiltonian.hpp";
    // ofstream output(output_h);

    // ********************************************************************** //
    // -------------------- Creating the the .h file ------------------------ //
    if(output.is_open()){
        int actually_complex;
        output << "#define N        " << no_qubit << endl;
        output << "#define Nop      " << no_ps << endl;
        output << "#define Ncycles  " << nullity << endl;
        output << endl;
        
        // ---------------- Permutation matrices and cycles --------------------- //
        // The permutation bitsets
        output << "std::bitset<N> P_matrix[Nop] = {";
        for(int i = 0; i < no_ps; i++){
            output << "std::bitset<N>(\"";
            for(int j=0; j < Ps_nontrivial[i].size(); j++){
                output << Ps_nontrivial[i][j];
            }
            output << "\")";
            if(i < no_ps - 1){
                output << ", ";
            }
        }
        output << "};" << endl;

        // The cycle bitsets
        output << "std::bitset<Nop> cycles[Ncycles] = {";
        for(int i = 0; i < nullity; i++){
            output << "std::bitset<Nop>(\"";
            for(int j=no_ps-1; j >= 0; j--){
                output << nullspace[i][j];
            }
            output << "\")";
            if(i < nullity-1){
                output << ", ";
            }
        }
        output << "};" << endl;
        output << endl; 

        // ---------------------------------------------------------------------- //
        // -------------------------- Diagonal terms ---------------------------- //
        output << "const int D0_size = " << Z_track[0].size() << ";" << endl;

        actually_complex = 0;
        for(int i=0;i<Z_track[0].size();i++) if(coefficients[Z_track[0][i]].imag()!=0) actually_complex = 1;

        if(actually_complex) output << "std::complex<double> "; else output << "double ";
        output << "D0_coeff[D0_size] = {";
        for(int i=0;i<Z_track[0].size();i++){
            output_complex(coefficients[Z_track[0][i]]);
            if(i < Z_track[0].size() - 1) output << ", ";
        }
        output << "};" << endl;
        output << "std::bitset<N> D0_product[D0_size] = {";
        for(int i = 0; i < Z_track[0].size(); i++){
            vector<int> Zs_i = Zs[Z_track[0][i]];
            output << "std::bitset<N>(\"" << int_to_str(Zs_i) << "\")";
            if(i < Z_track[0].size() - 1){
                output << ", ";
            }
        }
        output << "};" << endl;
        output << endl; 

        // ------------------------------------------------------------------------ //
        // --------------------------- Off-Diagonal terms ------------------------- //
        // Finding D_maxsize:
        int D_max = 0 , D_start;
        vector<int> D_size;
        if(Ps_nontrivial.size() == Ps.size()) // In case of no diagonal term (i.e. D_0 = 0)
        {
            D_start = 0;
        }else{
            D_start = 1;
        }
        for(int i = D_start; i < Z_track.size(); i++){
            int z_size_i = Z_track[i].size();
            D_size.push_back(z_size_i);
            if(z_size_i > D_max){
                D_max = z_size_i;
            }
        }
        output << "const int D_maxsize = " << D_max << ";" << endl;
        output << "int D_size[Nop] = {";
        for(int i=0;i<no_ps;i++){
            output << D_size[i];
            if(i < no_ps - 1){
                output << ", ";
            }
        }
        output << "};" << endl;


        actually_complex = 0;
        for(int i = D_start; i < Z_track.size() ; i++)
           for(int j = 0; j < Z_track[i].size(); j++)
              if(coefficients[Z_track[i][j]].imag()!=0) actually_complex = 1;

        if(actually_complex) output << "std::complex<double> "; else output << "double ";
        output << "D_coeff[Nop][D_maxsize] = {";
        for(int i = D_start; i < Z_track.size() ; i++){
            output << "{";
            for(int j = 0; j < Z_track[i].size(); j++){
                output_complex(coefficients[Z_track[i][j]]);
                if(j < Z_track[i].size() - 1) output << ", ";
            }
            output << "}";
            if(i < Z_track.size()-1){
                output << ", ";
            }
        }
        output << "};" << endl;
        output << "std::bitset<N> D_product[Nop][D_maxsize] = {";
        for(int i = D_start; i < Z_track.size() ; i++){
            output << "{";
            for(int j = 0; j < Z_track[i].size(); j++){
                vector<int> z_ij = Zs[Z_track[i][j]];
                output << "std::bitset<N>(\"" << int_to_str(z_ij) << "\")";
                if(j < Z_track[i].size() - 1){
                    output << ", "; 
                }
            }
            output << "}";
            if(i < Z_track.size()-1){
                output << ", ";
            }
        }
        output << "};" << endl;

	if(no_ps == 0) output << endl << "#pragma GCC diagnostic ignored \"-Wdiv-by-zero\"" << endl;

        // output.close();
    }
}

template<typename T>
void PrintListPure(ofstream& output, T* list, int len){
	output << "{";
	for(int i=0;i<len;i++){
		output << list[i];
		if(i<len-1) output << ", ";
	}
	output << "}";
}

template<typename T>
void PrintList(ofstream& output, T* list, int len, string namelist){
	output << namelist << " = ";
        PrintListPure(output, list, len);
	output << ";" << endl;
}

void main2(int argc , char* argv[]){
    int Nobservables = argc - 2;
    vector<pair<complex<double>, vector<int>>> Op_data[Nobservables];

    string Ham_fileName(argv[1]);  // Reading the name of the input .txt file describing the Hamiltonian
    vector<pair<complex<double>, vector<int>>> Ham_data = data_extract(Ham_fileName);

    string Op_fileNames[Nobservables];
    for(int O=0;O<Nobservables;O++){
        Op_fileNames[O] = argv[O+2];  // Reading the name of the input .txt file describing the Hamiltonian
        Op_data[O] = data_extract(Op_fileNames[O]);
    }

    // Unpacking the Hamiltonian data from the input file 
    PZdata PZ_data_Ham = PZcomp(Ham_data);
    vector<dbitset> Ps_bit_Ham = PZ_data_Ham.Ps;
    vector<vector<bool>> Ps_Ham = downsize_bitset(Ps_bit_Ham);
    vector<complex<double>> coefficients_Ham = PZ_data_Ham.coeffs;
    vector<vector<int>> Z_track_Ham = PZ_data_Ham.Z_track;
    vector<vector<int>> Zs_Ham = PZ_data_Ham.Zs;
    vector<string> Zs_string_Ham;

    // Unpacking the Operator data from the input files
    PZdata PZ_data_Op[Nobservables];
    vector<dbitset> Ps_bit_Op[Nobservables];
    vector<vector<bool>> Ps_Op[Nobservables];
    vector<complex<double>> coefficients_Op[Nobservables];
    vector<vector<int>> Z_track_Op[Nobservables];
    vector<vector<int>> Zs_Op[Nobservables];
    vector<string> Zs_string_Op[Nobservables];
    bool D0_exists = false;
    bool P0_exists[Nobservables] = {false};

    for(int O=0;O<Nobservables;O++){
        PZ_data_Op[O] = PZcomp(Op_data[O]);
        Ps_bit_Op[O] = PZ_data_Op[O].Ps;
        Ps_Op[O] = downsize_bitset(Ps_bit_Op[O]);
        coefficients_Op[O] = PZ_data_Op[O].coeffs;
        Z_track_Op[O] = PZ_data_Op[O].Z_track;
        Zs_Op[O] = PZ_data_Op[O].Zs;
    }
    
    // Converting the Z indices into string of bitsets!
    for(int i = 0; i < Zs_Ham.size();i++){
        Zs_string_Ham.push_back(int_to_str(Zs_Ham[i]));
    }

    // Converting the Z indices into string of bitsets!
    for(int O=0;O<Nobservables;O++) for(int i = 0; i < Zs_Op[O].size();i++){
        Zs_string_Op[O].push_back(int_to_str(Zs_Op[O][i]));
    }

    // Convert Ps indices into vector of bools
    vector<vector<bool>> Ps_Ham_nontrivial = Ps_Ham;
    if(none(Ps_bit_Ham[0])){
        Ps_Ham_nontrivial.erase(Ps_Ham_nontrivial.begin());
        D0_exists = true;
    }

    // Convert Ps indices into vector of bools
    vector<vector<bool>> Ps_Op_nontrivial[Nobservables];
    for(int O=0;O<Nobservables;O++){
       Ps_Op_nontrivial[O]= Ps_Op[O];
       if(none(Ps_bit_Op[O][0])){
           Ps_Op_nontrivial[O].erase(Ps_Op_nontrivial[O].begin());
           P0_exists[O] = true;
       }
    }

    vector<vector<int>> Ps_binary_Ham = bit_to_intvec(Ps_Ham_nontrivial);
    vector<vector<int>> Ps_binary_Op[Nobservables];
    for(int O=0;O<Nobservables;O++) Ps_binary_Op[O] = bit_to_intvec(Ps_Op_nontrivial[O]);

    // Combine the permutations of Hamiltonians with each permutation from the operator and compute the nullspace
    int no_ps = Ps_Ham_nontrivial.size();
    int no_ops[Nobservables];
    for(int O=0;O<Nobservables;O++) no_ops[O] = Ps_Op_nontrivial[O].size();

    vector<vector<int>> Op_to_Ham[Nobservables];
    vector<int> zero_vector;
    for(int i = 0; i < no_ps; i++){
        zero_vector.push_back(0);
    }

    for(int O=0;O<Nobservables;O++) for(int k=0; k < Ps_binary_Op[O].size(); k++){
        vector<vector<int>> Ps_bin_total = Ps_binary_Ham, nullspace_k; 
        bool permutation_found = false;
        Ps_bin_total.push_back(Ps_binary_Op[O][k]);
        nullspace_k = Null2(Ps_bin_total);
        while(Eig_minimize(nullspace_k));
        int min_ops=100000 , min_index;
        for(int i = 0; i< nullspace_k.size(); i++){
            int sum_ops = 0;
            if(nullspace_k[i][no_ps] == 1){
                permutation_found = true;
                sum_ops = Int_sum(nullspace_k[i]);
                if(sum_ops < min_ops){
                    min_ops = sum_ops;
                    min_index = i;
                }
            }
        }
        if(permutation_found){
            Op_to_Ham[O].push_back(nullspace_k[min_index]);
        }
        else{
            Op_to_Ham[O].push_back(zero_vector);
        }
    }

    // string output_operator = "observables.hpp";
    // ofstream output(output_operator);

    vector<int> rel_perms[Nobservables];
    bool non_triv_offdiags_exists[Nobservables] = {false};

    for(int O=0;O<Nobservables;O++) for(int i = 0; i < no_ops[O]; i++){
        if(Int_sum(Op_to_Ham[O][i]) > 0){
            rel_perms[O].push_back(i);
        }
    }
    // Delete the Z_tracks of the trivial permutations:
    vector<vector<int>> Z_track_Op_kept[Nobservables];

    for(int O=0;O<Nobservables;O++){
       int counter = 0;
       for(int i = 0 ; i < no_ops[O]; i++){
           if(rel_perms[O][counter] == i){
               Z_track_Op_kept[O].push_back(Z_track_Op[O][i + int(P0_exists[O])]);
               counter++;
           }
       }
       no_ops[O] = rel_perms[O].size();
    }

    if(output.is_open()){
        int actually_complex;
	output << endl << endl;
	output << "#define MEASURE_CUSTOM_OBSERVABLES" << endl << endl;
        output << "const int Nobservables = " << Nobservables << ";" << endl;

        int no_ops_max = 0; for(int O=0; O<Nobservables; O++) no_ops_max = max(no_ops_max, no_ops[O]);

        output << "const int MNop_max = " << no_ops_max << ";" << endl;

        PrintList(output, no_ops, Nobservables, "int MNop[Nobservables]");
        output << endl;

        for(int O=0; O<Nobservables; O++) Op_fileNames[O] = "\"" + Op_fileNames[O]/*.substr(0,Op_fileNames[O].find_last_of("."))*/ + "\"";
        PrintList(output, Op_fileNames, Nobservables, "std::string Mnames[Nobservables]");

        output << "std::bitset<Nop> MP[Nobservables][MNop_max] = {";

        for(int O=0; O<Nobservables; O++){
               output << "{";
               non_triv_offdiags_exists[O] = true;
               for(int i = 0; i < no_ops[O]; i++){
                   output << "std::bitset<Nop>(\"";
                   for(int j=no_ps-1; j >= 0; j--){
                      output << Op_to_Ham[O][rel_perms[O][i]][j];
                   }
                   output << "\")";
                   if(i < no_ops[O]-1){
                       output << ", ";
                   }
               }
               output << "}";
               if(O<Nobservables-1) output << ",";
        }
	output << "};" << endl << endl;

        int MD0_size[Nobservables]; for(int O=0; O<Nobservables; O++) MD0_size[O] = Z_track_Op[O][0].size();
        int MD0_maxsize = 0; for(int O=0; O<Nobservables; O++) MD0_maxsize = max(MD0_maxsize, MD0_size[O]);

        output << "const int MD0_maxsize = " << MD0_maxsize << ";" << endl;

        PrintList(output, MD0_size, Nobservables, "int MD0_size[Nobservables]");

        actually_complex = 0;
        for(int O=0; O<Nobservables; O++)
           if(P0_exists[O]) for(int i=0;i<Z_track_Op[O][0].size();i++)
              if(coefficients_Op[O][Z_track_Op[O][0][i]].imag()!=0) actually_complex = 1;

        if(actually_complex) output << "std::complex<double> "; else output << "double ";
        output << "MD0_coeff[Nobservables][MD0_maxsize] = {";

        for(int O=0; O<Nobservables; O++){
            output << "{";
            if(P0_exists[O]) for(int i=0;i<Z_track_Op[O][0].size();i++){
                output_complex(coefficients_Op[O][Z_track_Op[O][0][i]]);
                if(i < Z_track_Op[O][0].size() - 1) output << ", ";
            }
            output << "}"; if(O<Nobservables-1) output << ",";
        }
       output << "};" << endl;

       output << "std::bitset<N> MD0_product[Nobservables][MD0_maxsize] = {";

       for(int O=0; O<Nobservables; O++){
            output << "{";
            if(P0_exists[O]) for(int i = 0; i < Z_track_Op[O][0].size(); i++){
                output << "std::bitset<N>(\"" << Zs_string_Op[O][Z_track_Op[O][0][i]] << "\")";
                if(i < Z_track_Op[O][0].size() - 1){
                    output << ", ";
                }
            }
            output << "}"; if(O<Nobservables-1) output << ",";
       }
       output << "};" << endl << endl;

       int MD_size[Nobservables][no_ops_max];
       int MD_maxsize[Nobservables] = {0};

       for(int O=0; O<Nobservables; O++) if(non_triv_offdiags_exists[O]){
            int D_max = 0;
            vector<int> D_size;
            for(int i = 0; i < no_ops[O]; i++){
                int z_size_i = Z_track_Op_kept[O][i].size();
                D_size.push_back(z_size_i);
                if(z_size_i > D_max){
                    D_max = z_size_i;
                }
            }
            MD_maxsize[O] = D_max;
            for(int i=0;i<no_ops[O];i++) MD_size[O][i] = D_size[i];
       } else MD_maxsize[O] = 0;

       int MD_Maxsize = 0; for(int O=0; O<Nobservables; O++) MD_Maxsize = max(MD_Maxsize, MD_maxsize[O]);

       output << "const int MD_Maxsize = " << MD_Maxsize << ";" << endl;

       PrintList(output, MD_maxsize, Nobservables, "int MD_maxsize[Nobservables]");

       output << "int MD_size[Nobservables][MNop_max] = {";
       for(int O=0; O<Nobservables; O++){ PrintListPure(output, MD_size[O], no_ops[O]); if(O<Nobservables-1) output<<", ";}
       output << "};" << endl;

       actually_complex = 0;
       for(int O=0; O<Nobservables; O++)
          if(non_triv_offdiags_exists[O]){
             for(int i = 0; i < no_ops[O]; i++)
                for(int j = 0; j < Z_track_Op_kept[O][i].size(); j++)
                   if(coefficients_Op[O][Z_track_Op_kept[O][i][j]].imag()!=0) actually_complex = 1;
          }

       if(actually_complex) output << "std::complex<double> "; else output << "double ";
       output << "MD_coeff[Nobservables][MNop_max][MD_Maxsize] = {";

       for(int O=0; O<Nobservables; O++){
         output << "{";
         if(non_triv_offdiags_exists[O]){
            for(int i = 0; i < no_ops[O]; i++){
                output << "{";
                for(int j = 0; j < Z_track_Op_kept[O][i].size(); j++){
                    output_complex(coefficients_Op[O][Z_track_Op_kept[O][i][j]]);
                    if(j < Z_track_Op_kept[O][i].size() - 1) output << ", ";
                }
                output << "}";
                if(i < no_ops[O]-1){
                    output << ", ";
                }
            }
         }
         output << "}"; if(O<Nobservables-1) output << ",";
       }
       output << "};" << endl;

       output << "std::bitset<N> MD_product[Nobservables][MNop_max][MD_Maxsize] = {";
       for(int O=0; O<Nobservables; O++){
            output << "{";
            for(int i = 0; i < Z_track_Op_kept[O].size() ; i++){
                output << "{";
                for(int j = 0; j < Z_track_Op_kept[O][i].size(); j++){
                    output << "std::bitset<N>(\"" << Zs_string_Op[O][Z_track_Op_kept[O][i][j]] << "\")";
                    if(j < Z_track_Op_kept[O][i].size() - 1){
                        output << ", "; 
                    }
                }
                output << "}";
                if(i < Z_track_Op_kept[O].size()-1){
                    output << ", ";
                }
            }
            output << "}"; if(O<Nobservables-1) output << ",";
       }
       output << "};" << endl;
       // output.close();
    }
}

int main(int argc , char* argv[]){
    if(argc < 2){ cout << "Usage: ./prepare.bin hamiltonian.txt observable_1.txt observable_2.txt ..." << endl; exit(1);}
    cout << "Preparing PMR for the Hamiltonian and computing the list of fundamental cycles...";
    output.open("hamiltonian.hpp");
    main1(argc, argv);
    cout << "done" << endl;
    if(argc >= 3){
       cout << "Preparing PMR for the observables and their representation via the permutation operators of the Hamiltonian...";
       main2(argc, argv);
       cout << "done" << endl;
    }
    output.close();
    return 0;
}
