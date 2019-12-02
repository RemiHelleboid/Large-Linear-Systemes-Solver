#include <vector>
#include <iostream>
#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"
#include "FredholmMatrix.hpp"
using namespace std;
using namespace Eigen;



int main(){
	int N;
	cout<<endl<<"------------------------------------------------\n------------------------------------------------"<<endl;
	TestResoLine(5);
	cout<<"------------------------------------------------\n------------------------------------------------"<<endl<<endl;
	cout<<"Entrez la taille de la matrice pour la SVD : ";
	cin>>N;
	cout<<N<<'x'<<N<<endl;
	TestSVDEigen(N);	
	TestSVDCross(N);	

	return(0);
}
