#include <iostream>
#include <vector>
#include "SparseMatrix.hpp"
#include <cmath>

using namespace std;


double SparseMatrix::operator()(int i, int j) const{
	int pos0=row[i];
	int pos1=row[i+1];
	for(int k=pos0; k<pos1; k++){
		if(col[k]==j){return(val[k]);}}
	return(double(0));
}

SparseMatrix::SparseMatrix(const DenseMatrix & M){
	size=M.taille();
	vector<int> NewRow(1,0);
	vector<int> NewCol;
	vector<double> NewVal;
	int NN=0;
	for(int i=0; i<M.taille(); i++){
		for(int j=0; j<M.taille(); j++){
			if(M(i,j)!=0){
				NN++;
				NewCol.push_back(j);
				NewVal.push_back(M(i,j));
			}}
		NewRow.push_back(NN);
		}
	get_row()=NewRow;
	get_col()=NewCol;
	get_val()=NewVal;
	}


ostream & operator<<(ostream & C ,const SparseMatrix & M){
	for(int i=0;i<M.size;i++){
		for(int j=0; j<M.size;j++){cout<<M(i,j)<<" ";}
			C<<endl;}
	return(C);
}

double FrobeniusNorm(const SparseMatrix& A){
	double s=0;
	for(unsigned int i=0; i<A.val.size(); i++) s+=A.val[i]*A.val[i];
	return(sqrt(s));
}

vector<int> SparseMatrix::Argmax(){
	double max=abs(val[0]);
	int k=0;
	int ligne=0;
	int colonne=col[0];
	for(unsigned int i=0; i<val.size(); i++){
		if((int)i>=row[k]) k++;
		if(abs(val[i])>max){max=abs(val[i]);colonne=col[i];ligne=k-1;}		
}	
	while(row[ligne+1]==0) ligne++;
	vector<int> arg({ligne, colonne});
	return(arg);	
}

vector<double> MVProd(const SparseMatrix &A,const vector<double> &V){
	int N=V.size();
	vector<double> R(N);
	if(A.size!=N){cout<<"PRODUIT MAT VECT TAILLE ERROR"<<endl;}
	for(int i=0; i<A.size; i++){
		int pos0=A.row[i];
		int pos1=A.row[i+1];
		double s=0;
			for(int j=pos0; j<pos1; j++){
				s+=A.val[j]*V[A.col[j]];
			}
		R[i]=s;
	}
	return(R);
}

void SparseMatrix::Insert(int k, int j, double value){
	if(value !=0){ 											//Si value=0 : on ne fait rien
	bool TestPresence=false;
	int pos0=row[k];
	int pos1=row[k+1];
	int P=pos0;
	for(int i=pos0; i<pos1; i++){
		while(col[i]<=j && TestPresence==false){			//Les colonnes soit en ordre croissant
		P=i;
		if(col[i]==j ){val[i]+=value;TestPresence=true;}	//La valeur est déja dans la matrice
		else{break;}
	}}
	if(TestPresence==false){								//La valeur est nouvelle
		col.insert(col.begin()+P, j);
		val.insert(val.begin()+P, value);
		for(int i=k+1; i<=size; i++){row[i]++;};
	}}
}
