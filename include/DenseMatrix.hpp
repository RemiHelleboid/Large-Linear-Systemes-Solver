#include <vector>
#include <iostream>

using namespace std;

#ifndef DMATRIX
#define DMATRIX

template<class T>
ostream & operator<<(ostream & os,const vector<T> & V){  //Affichage des vecteurs
		if(V.size()==0) return(os);
	os<<'(';
	for(int i=0;i<(int)V.size()-1;i++) os<<V[i]<<", ";
	os<<V[V.size()-1]<<" )";
	return(os);
}

class DenseMatrix{
	private:
		int size;
		vector<double> val;
	public:
	//Constructeurs
		DenseMatrix(int n=0):size(n), val(vector<double>(n*n,0)){};		//Constructeur par defaut
		DenseMatrix(const DenseMatrix& M):size(M.size), val(M.val){};	//Constructeur par recopie
		DenseMatrix(int n, char a);  									//Matrice B
		DenseMatrix(int n, double d, double low, double up);  			//Matrice tribande
		DenseMatrix& operator=(const DenseMatrix &)=default;
	//Accesseur, Mutateurs, lecteurs, afficheurs
		int taille()const{return(size);};
		void Load(const char* filename);
		double get_val(int k)const{return val[k];}
		friend ostream & operator<<(ostream & C, DenseMatrix & M);
		double operator()(int i, int j) const{return(val[size*i+j]);};          
		double& operator()(int i, int j) {return(val[size*i+j]);};
	//Méthode
		vector<int> Argmax();											//(ligne, colone) qui maximise la matrice en valeur absolue
		void decomp_LU();
		friend vector<double> ResoLU(DenseMatrix &A, vector<double> &V);
		friend vector<double> inverse_triang_sup(DenseMatrix &A, vector<double> &V);
		friend vector<double> inverse_triang_inf(DenseMatrix &A, vector<double> &V);
		friend double FrobeniusNorm(const DenseMatrix& A);
	};

void TestResoLine(int N);

#endif
