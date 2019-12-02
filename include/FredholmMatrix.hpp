#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "SparseMatrix.hpp"
#include "DenseMatrix.hpp"
#include "redsvd.hpp"
#include "redsvdFile.hpp"



using namespace std;

#ifndef FREDI
#define FREDI

class FredholmMatrix{
	private:
		int size;
		SparseMatrix diag;
		vector<double> sigma;
		vector<vector<double>> lsv, rsv;
	public:
	//Constructeurs
		FredholmMatrix(int n0=0):size(n0), diag(SparseMatrix(n0)), sigma(vector<double>(0)){};
		FredholmMatrix(const FredholmMatrix &M):size(M.size), diag(M.diag), sigma(M.sigma), lsv(M.lsv), rsv(M.rsv){};
		FredholmMatrix& operator=(const FredholmMatrix &)=default;
		FredholmMatrix(const DenseMatrix &);
		FredholmMatrix(const SparseMatrix & D):size(D.taille()), diag(D){};
	//Accesseurs
		int taille()const{return(size);}
		vector<double> ValSing()const{return sigma;};
		SparseMatrix get_diag()const{return(diag);};
		void SVD ( const DenseMatrix & B, int r);
		void SVDEigen ( const DenseMatrix & B, int r);
		void CrossApprox(const DenseMatrix & B, int r);	
		void CrossApproxCompute(int r);
	//Methodes
		void Insert(int k, int j, double value);
		void Insert(double s, vector<double>& u, vector<double>& v);
		double operator()(int i, int j)const;
		double CoefSVD(int i, int j, int p)const;				//Coeffs A(i,j) de la SVD tronqué au rang p
	//Fonctions amicales
		friend ostream & operator<<(ostream & os,const FredholmMatrix & M);
		friend double FrobeniusNorm(const FredholmMatrix& A);
		friend double FHMDiff(const FredholmMatrix& F,const SparseMatrix &S, int p);  //Norme de Frobenius de F-S
		friend void MVProd(vector<double>& b, const FredholmMatrix& A, const vector<double>& x);
		friend void Solve(vector<double> &b, const FredholmMatrix& A, vector<double>& x);
		
	};

	
ostream & operator<<(ostream & os,const vector<double> & V);  //Affichage des vecteurs

void TestSVDEigen(int N);
void TestSVDCross(int N);



#endif



