#include <vector>
#include <iostream>
#include "DenseMatrix.hpp"
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

ostream & operator<<(ostream & C , DenseMatrix & M){
	for(int i=0;i<M.size;i++){
		for(int j=0; j<M.size;j++){cout<<M.val[i*M.size+j]<<" ";}
			C<<endl;}
	return(C);
}

DenseMatrix::DenseMatrix(int n, char a){
	size=n;
	val.resize(n*n);
	for(int k=0; k<n;k++){for(int j=0; j<n ;j++){
		operator()(k,j)=exp(-((j-k)*(j-k))/double(n*n));
		}}
}

DenseMatrix::DenseMatrix(int n, double d, double low, double up){
	size=n;
	val.resize(n*n);
	for(int k=0; k<n;k++) operator()(k,k)=d;
	for(int k=1; k<n;k++) operator()(k,k-1)=up;
	for(int k=1; k<n;k++) operator()(k-1,k)=low;
};

void DenseMatrix::Load(const char* filename){
	ifstream F(filename);
	if(!F){cout<<"Erreur : Le fichier n'existe pas."<<endl;}
	string s;
	getline(F,s);
	F>>size;
	F>>size;
	val.resize(size*size);
	F>>s;
	for(int i=0;i<size*size;i++){
			F>>val[i];
		}
	}

vector<int> DenseMatrix::Argmax(){
	double max=0;
	int ligne=0;
	int col=0;
	vector<int> arg(2,0);
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			if(abs(operator()(i,j))>max) {max=abs(operator()(i,j)); ligne=i; col=j;}
		}}
	arg[0]=ligne;
	arg[1]=col;
	return(arg);	
	}


void DenseMatrix::decomp_LU(){
	DenseMatrix S(*this);
	for(int l=1; l<size;l++){						//Calcul du bloc L10
		operator()(l,0)=(1.0/S(0,0))*S(l,0);}
	for(int c=1; c<size; c++){						//Calcul du bloc U01
		operator()(0,c)=S(0,c);}
	for(int i=1; i<size; i++){					//Calcul S11 
		for(int j=1; j<size; j++){
			S(i,j)=S(i,j)-(1.0/S(0,0))*S(i,0)*S(0,j);
		}}
		
	for(int k=1; k<size; k++){
		operator()(k,k)=S(k,k);
		for(int l=k+1; l<size;l++){						//Calcul du bloc L(k+1,k)
			operator()(l,k)=(1.0/operator()(k,k))*S(l,k);}
		for(int c=k+1; c<size; c++){
			operator()(k,c)=S(k,c);
		}
		for(int i=k+1; i<size; i++){					//Calcul S11 
			for(int j=k+1; j<size; j++){
			S(i,j)=S(i,j)-(1.0/S(k,k))*S(i,k)*S(k,j);
		}}
	}
}

double FrobeniusNorm(const DenseMatrix& A){
	double s=0;
	double coef=0;
	for(int i=0; i<A.size; i++){
		for(int j=0; j<A.size; j++){
		coef=A(i,j);
		s+=coef*coef;}
	}	
	return sqrt(s);
}

void TestResoLine(int N){
	cout<<"Resolution du système Ax=b avec: \n\t-A la matrice du Laplacien (cf TP4) \n\t-b=(1,1,...,1)"<<endl;
	DenseMatrix A(N, 2.0, -1.0,-1.0);
	vector<double> b(N,1.0),x;
	x=ResoLU(A,b);
	cout<<x<<endl;
}

