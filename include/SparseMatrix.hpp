#include <iostream>
#include <vector>
#include "DenseMatrix.hpp"


using namespace std;

#ifndef SPARSEMAT
#define SPARSEMAT

class SparseMatrix{
	private:
		int size;
		vector<int> row;
		vector<int> col;
		vector<double> val;
	public:
	//Constructeurs
		SparseMatrix(int n=0):size(n),row(vector<int>(n+1)), col(vector<int>(0)), val(vector<double>(0,0.0)){}; //Const. nul
		SparseMatrix(const DenseMatrix &);
		SparseMatrix& operator=(const SparseMatrix &)=default;

	//Accesseurs, Mutateurs
		double operator()(int i, int j) const; 	         
		int taille()const {return size;};
		vector<int> get_row()const {return(row);};
		vector<int> get_col()const {return(col);};
		vector<double> get_val()const {return(val);};
		vector<int>& get_row() {return(row);};
		vector<int>& get_col() {return(col);};
		vector<double>& get_val() {return(val);};	
	//Méthodes	
		vector<int> Argmax();									//(ligne, colone) qui maximise la matrice en valeur absolue
		void Insert(int k, int j, double value);
	//Fonctions amicales
		friend ostream & operator<<(ostream & C,const SparseMatrix & M);
		friend vector<double> MVProd(const SparseMatrix &A,const vector<double> &V);
		friend double FrobeniusNorm(const SparseMatrix& A);

	};

#endif
