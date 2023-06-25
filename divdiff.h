//
// These routines are introduced in the paper:
// L. Gupta, L. Barash, I. Hen, Calculating the divided differences of the exponential function by addition and removal of inputs, Computer Physics Communications 254, 107385 (2020)
//
// This program is licensed under a Creative Commons Attribution 4.0 International License:
// http://creativecommons.org/licenses/by/4.0/
//

#include<random>
#include<stdio.h>
#include<string.h>
#include<math.h>

double* invPowersOf2 {NULL};
const int maxexp = 100000;
const int extralen = 10;

template <typename T> T min (T a, T b) { return (b>=a)?a:b;}
template <typename T> T max (T a, T b) { return (b>=a)?b:a;}

class ExExFloat{
private:
	double mantissa;
	int exponent;
public:
	ExExFloat(){	mantissa = 0.5; exponent = 1;}
	void normalize(){ int tmp; mantissa = frexp(mantissa,&tmp); exponent += tmp;}
	void init_expmu(double mu){ double e = mu*1.4426950408889634; exponent = ceil(e); mantissa = pow(2.,e - ceil(e)); }
	void print(){
		double exp10, m;
		exp10 = (exponent*0.30102999566398114);
		m = mantissa*pow(10,exp10 - floor(exp10)); 
		exp10 = floor(exp10); if(fabs(m)<1){ exp10--; m*=10; }
		if((exp10<7)&&(exp10>-7)) printf("%.17f",m*pow(10,exp10)); else printf("%.17fe%.0f",m,exp10);
	}
	ExExFloat operator =(ExExFloat const &obj){ mantissa = obj.mantissa; exponent = obj.exponent;	return *this;}
	ExExFloat operator =(double const &obj){ mantissa = obj; exponent = 0; normalize(); return *this;}
	ExExFloat operator +(ExExFloat const &obj){ // important restriction: it is assumed here that each of the summands is not equal to zero
		ExExFloat res;
		if(obj.exponent >= exponent){
			res.mantissa = obj.mantissa + mantissa*invPowersOf2[obj.exponent - exponent];
			res.exponent = obj.exponent; res.normalize();
		} else{
			res.mantissa = mantissa + obj.mantissa*invPowersOf2[exponent - obj.exponent];
			res.exponent = exponent; res.normalize();
		}
		return res;
	}
	ExExFloat operator -(ExExFloat const &obj){ // important restriction: it is assumed here that each of the summands is not equal to zero
		ExExFloat res;
		if(obj.exponent >= exponent){
			res.mantissa = mantissa*invPowersOf2[obj.exponent - exponent] - obj.mantissa;
			res.exponent = obj.exponent; res.normalize();
		} else{
			res.mantissa = mantissa - obj.mantissa*invPowersOf2[exponent - obj.exponent];
			res.exponent = exponent; res.normalize();
		}
		return res;
	}
	ExExFloat operator +=(ExExFloat const &obj){ // important restriction: it is assumed here that each of the summands is not equal to zero
		if(obj.exponent >= exponent){
			mantissa = obj.mantissa + mantissa*invPowersOf2[obj.exponent - exponent];
			exponent = obj.exponent; normalize();
		} else{
			mantissa = mantissa + obj.mantissa*invPowersOf2[exponent - obj.exponent];
			exponent = exponent; normalize();
		}
		return *this;
	}
	ExExFloat operator -=(ExExFloat const &obj){ // important restriction: it is assumed here that each of the summands is not equal to zero
		if(obj.exponent >= exponent){
			mantissa = mantissa*invPowersOf2[obj.exponent - exponent] - obj.mantissa;
			exponent = obj.exponent; normalize();
		} else{
			mantissa = mantissa - obj.mantissa*invPowersOf2[exponent - obj.exponent];
			exponent = exponent; normalize();
		}
		return *this;
	}
	ExExFloat operator *(ExExFloat const &obj){
		ExExFloat res;	res.mantissa = mantissa * obj.mantissa;	res.exponent = exponent + obj.exponent;
		res.normalize(); return res;
	}
	ExExFloat operator /(ExExFloat const &obj){
		ExExFloat res;	res.mantissa = mantissa / obj.mantissa;
		res.exponent = exponent - obj.exponent;	res.normalize(); return res;
	}
	ExExFloat operator *(double const &obj){ ExExFloat res; res.mantissa = mantissa * obj; res.exponent = exponent; res.normalize(); return res; }
	ExExFloat operator /(double const &obj){ ExExFloat res; res.mantissa = mantissa / obj; res.exponent = exponent; res.normalize(); return res; }
	ExExFloat operator *=(ExExFloat const &obj){ mantissa *= obj.mantissa; exponent += obj.exponent; normalize(); return *this;}
	ExExFloat operator /=(ExExFloat const &obj){ mantissa /= obj.mantissa; exponent -= obj.exponent; normalize(); return *this;}
	ExExFloat operator *=(double const &obj){ mantissa *= obj; normalize(); return *this;}
	ExExFloat operator /=(double const &obj){ mantissa /= obj; normalize(); return *this;}
	int operator >=(double const &r){ // important restriction: it is assumed here that both values of mantissa are not negative
		if(r == 0) return (mantissa >= 0);
		else{
			ExExFloat R; R = r;
			if(exponent > R.exponent ) return 1;
			else if((exponent == R.exponent)&&(mantissa >= R.mantissa)) return 1;
			else return 0;
		}
	}
	int operator >=(ExExFloat const &r){ // important restriction: it is assumed here that both values of mantissa are not negative
		if(r.mantissa == 0) return (mantissa >= 0);
		else{
			if(exponent > r.exponent ) return 1;
			else if((exponent == r.exponent)&&(mantissa >= r.mantissa)) return 1;
			else return 0;
		}
	}
	double get_double(){ return ldexp(mantissa,exponent);}
	int sgn(){ return (mantissa > 0.0) - (mantissa < 0.0);}
	ExExFloat abs(){ExExFloat res; res.mantissa = fabs(mantissa); res.exponent = exponent; return res;}
	ExExFloat SqRt(){ ExExFloat res;
		if(exponent%2 == 0){ res.mantissa = sqrt(mantissa); res.exponent = exponent/2;} 
		else{ res.mantissa = sqrt(2*mantissa); res.exponent = (exponent-1)/2;}
		res.normalize(); return res;
	}
};

void divdiff_init(){
	invPowersOf2 = new double[maxexp];
	double curr=1; for(int i=0;i<maxexp;i++){ invPowersOf2[i] = curr; curr/=2; }
}

void divdiff_clear_up(){ delete[] invPowersOf2; }

double mean(double* z, int n){
	double sum=0; int i;
	for(i=0;i<n;i++) sum+=z[i];
	return sum/n;
}

double maxAbsDiff(double* z, int len){
	double zmax = z[0], zmin = z[0]; int i;
	for(i=1;i<len;i++) { zmin = min(zmin,z[i]); zmax = max(zmax,z[i]);}
	return fabs(zmax-zmin);
}

class divdiff{

protected:

double *z2;
ExExFloat *h, *ddd;
int s, maxlen = 10001, smax = 500; 
double mu; ExExFloat expmu;

public:

double *z;
ExExFloat *divdiffs;
int CurrentLength;

long long int getsize(){
	long long int sum=0;
	sum += sizeof(double)*maxlen;
	sum += 2*sizeof(ExExFloat)*(maxlen+extralen+1);
	sum += sizeof(ExExFloat)*maxlen*smax;
	sum += 2*sizeof(double*)+3*sizeof(ExExFloat*)+4*sizeof(int)+sizeof(double)+sizeof(ExExFloat);
	return sum;
}

divdiff(int maxlen_, int smax_){ // constructor
	maxlen = maxlen_; smax = smax_; AllocateMem();
	if(invPowersOf2==NULL){
		printf("Error: invPowersOf2 has not been initialized\n");
		exit(EXIT_FAILURE);
	}
}

divdiff(const divdiff& other){ // copy constructor
	maxlen=other.maxlen; smax=other.smax; AllocateMem();
	CurrentLength=other.CurrentLength; s=other.s; mu=other.mu; expmu=other.expmu;
	memcpy(z,other.z,maxlen*sizeof(double));
	memcpy(h,other.h,(maxlen+extralen+1)*sizeof(ExExFloat));
	memcpy(divdiffs,other.divdiffs,(maxlen+extralen+1)*sizeof(ExExFloat));
	memcpy(ddd,other.ddd,maxlen*smax*sizeof(ExExFloat));
}

divdiff& operator=(const divdiff& other){ // copy assignment operator
	FreeMem();
	maxlen=other.maxlen; smax=other.smax; AllocateMem();
	CurrentLength=other.CurrentLength; s=other.s; mu=other.mu; expmu=other.expmu;
	memcpy(z,other.z,maxlen*sizeof(double));
	memcpy(h,other.h,(maxlen+extralen+1)*sizeof(ExExFloat));
	memcpy(divdiffs,other.divdiffs,(maxlen+extralen+1)*sizeof(ExExFloat));
	memcpy(ddd,other.ddd,maxlen*smax*sizeof(ExExFloat));
	return *this;
}

~divdiff(){ // destructor
	FreeMem();
}

void AllocateMem(){
	z = new double[maxlen]; h = new ExExFloat[maxlen+extralen+1];
	divdiffs = new ExExFloat[maxlen+extralen+1]; ddd = new ExExFloat[maxlen*smax];
	CurrentLength = 0; s = 1;
}

void FreeMem(){
	delete[] z; delete[] h; delete[] divdiffs; delete[] ddd;
}

void PrintList(ExExFloat* list, int len, const char* namelist){
	int i,flag=1;
	printf("%s={",namelist);
	for(i=0;i<len;i++){
		list[i].print();
		if(i<len-1) printf(",");
	}
	printf("};\n");
}

void PrintList(double* list, int len, const char* namelist){
	int i,flag=1;
	printf("%s={",namelist);
	for(i=0;i<len;i++){
		printf("%.17g",list[i]);
		if(i<len-1) printf(",");
	}
	printf("};\n");
}

void PrintList(int* list, int len, const char* namelist){
	int i,flag=1;
	printf("%s={",namelist);
	for(i=0;i<len;i++){
		printf("%d",list[i]);
		if(i<len-1) printf(",");
	}
	printf("};\n");
}

void Backupz(int len){ z2 = new double[len]; memcpy(z2,z,len*sizeof(double));}
void Restorez(int len){ memcpy(z,z2,len*sizeof(double)); delete[] z2;}

int s_changed(){
	return fabs(z[CurrentLength-1]-mu)/3.5 > s;
}

void AddElement(double znew, int force_s = 0, double force_mu = 0){
	int j,k,n,l,N; ExExFloat curr; n = CurrentLength; l = n + 1; N = maxlen+extralen; z[n] = znew; CurrentLength++;
	if(CurrentLength==1){
		s = (force_s == 0) ? 1 : force_s;
		mu = (force_mu == 0) ? z[0] : force_mu;
		expmu.init_expmu(mu);
		h[0] = 1; for(k=1;k<=N;k++) h[k] = h[k-1]/s;
		if(mu != z[0]) for(k=N;k>0;k--) h[k-1] += h[k]*(z[0]-mu)/k;
		curr = expmu*h[0]; for(k=0;k<s-1;k++) { ddd[k*maxlen] = curr; curr*=h[0];}
		divdiffs[0].init_expmu(z[0]); // alternatively: divdiffs[0] = curr;
	} else if(s_changed()||(CurrentLength>=maxlen)) AddAll(CurrentLength, force_s);
	else{
		for(k=N;k>n;k--) h[k-1] += h[k]*(z[n]-mu)/k; curr = expmu*h[n];
		for(k=n;k>=1;k--) h[k-1] = (h[k-1]*n + h[k]*(z[n]-z[n-k]))/(n-k+1);
		for(k=0;k<s-1;k++){
			ddd[k*maxlen+n] = curr;
			curr = ddd[k*maxlen]*h[n]; for(j=1;j<=n;j++) curr += ddd[k*maxlen+j]*h[n-j];
		}
		divdiffs[n] = curr;
	}
}

void RemoveElement(){ 
	int k,n,N;
	if(CurrentLength>=1){ 
		n = CurrentLength - 1; N = maxlen+extralen;
		for(k=1;k<=n;k++) h[k-1] = (h[k-1]*(n-k+1) - h[k]*(z[n]-z[n-k]))/n;
		for(k=n+1;k<=N;k++) h[k-1] -= h[k]*(z[n]-mu)/k;
		CurrentLength--;
	}
}

int RemoveValue(double value, int force_s = 0, double force_mu = 0){ // remove from bulk
	int j,k,n = CurrentLength - 1,found = 0;
	for(k=n;k>=0;k--) if(z[k] == value){
		for(j=n;j>=k;j--) RemoveElement();
		for(j=k;j<n;j++) AddElement(z[j+1],force_s,force_mu);
		found = 1; break;
	}
	return found;
}

void AddAll(int len, int force_s = 0){ // input is taken from z, output is written to divdiffs, size of z should be not smaller than len
	int i,s; CurrentLength = 0;
	if(force_s == 0) s = (int)ceil(maxAbsDiff(z,len)/3.5); else s = force_s;
	if((s > smax)||(len >= maxlen)){
		i = maxlen; Backupz(i); FreeMem();
		if(s > smax) smax = max(smax*2,s);
		if(len >= maxlen) maxlen = max(maxlen*2,len);
		AllocateMem(); Restorez(i);
	}
	AddElement(z[0],s,mean(z,len));
	for(i=1;i<len;i++) AddElement(z[i]); // calculates the vector (d[z_0], 1! d[z_0,z_1], ..., n! d[z_0,z_1,...,z_n]).
}
};
