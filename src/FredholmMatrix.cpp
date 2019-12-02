#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include "FredholmMatrix.hpp"
#include "DenseMatrix.hpp"

const float EPS = 0.00001f;

using namespace std;

void FredholmMatrix::Insert(double s, vector<double>& u, vector<double>& v){
	sigma.push_back(s);
	lsv.push_back(u);
	rsv.push_back(v);
}

double FredholmMatrix::operator()(int i, int j)const{
	double s=0;
	s+=diag(i,j);
	int r=sigma.size();
	for(int k=0; k<r; k++){
		s+=sigma[k]*lsv[k][i]*rsv[k][j];}
	return(s);
}

double FredholmMatrix::CoefSVD(int i, int j, int p)const{
	double s=0;
	int r=sigma.size();
	for(int k=0; k<p && k<r; k++){
		s+=sigma[k]*lsv[k][i]*rsv[k][j];}
	return(s);
}		



FredholmMatrix::FredholmMatrix(const DenseMatrix & M){
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
	diag.get_row()=NewRow;
	diag.get_col()=NewCol;
	diag.get_val()=NewVal;
}


void FredholmMatrix::Insert(int k, int j, double value){
	this->diag.Insert(k, j, value);
}



ostream & operator<<(ostream & os,const FredholmMatrix & M){
	for(int i=0;i<M.size;i++){
		for(int j=0; j<M.size;j++){os<<M(i,j)<<" ";}
			os<<endl;}
	return(os);
}

void MVProd(vector<double>& b, const FredholmMatrix& A, const vector<double>& x){
	if((int)b.size()!=A.size) cout<<"Produit matrice vecteur : Erreur de taille. b a été redimmensionné"<<endl;
	for(int i=0;i<A.size;i++){
		double s=0;
		for(int j=0; j<A.size;j++){
			s+=A(i,j)*x[j];}
		b[i]=s;
	}
}

vector<double> operator*(double p, vector<double> V){
	for(int i=0; i<(int)V.size(); i++) V[i]=p*V[i]; 
	return(V);
}

vector<double> inverse_triang_sup(DenseMatrix &A, vector<double> &V){
	int N=A.taille();
	vector<double> U(N);
	U[N-1]=V[N-1]/A(N-1,N-1);
	for(int i=N-2; i>=0; i--){
		double s=0;
		for(int j=N-1; j>i; j--){
			 s+=A(i,j)*U[j];	
		 }
		 U[i]=(V[i]-s)/A(i,i);
	 }
	 return(U);
}

vector<double> inverse_triang_inf(DenseMatrix &A, vector<double> &V){
	int N=A.taille();
	vector<double> U(N);			//Vecteur solution
	U[0]=V[0]/A(0,0);
	for(int i=0; i<N; i++){
		double s=0;
		for(int j=0; j<i; j++){
			 s+=A(i,j)*U[j];	
		 }
		 U[i]=(V[i]-s)/A(i,i);
	 }
	 return(U);
}


vector<double> ResoLU(DenseMatrix &A, vector<double> &V){
	DenseMatrix L(A.taille());
	DenseMatrix U(A.taille());
	A.decomp_LU();
	for(int i=0; i<A.taille(); i++){
		for(int j=i; j<A.taille(); j++){
			U(i,j)=A(i,j);}}
	for(int i=0; i<A.taille(); i++){
		for(int j=0; j<=i; j++){
			if(i==j){L(i,j)=1;}
			else{L(i,j)=A(i,j);}}}
	vector<double> Y;
	Y=inverse_triang_inf(L,V);
	vector<double> X;
	X=inverse_triang_sup(U,Y);
	return(X);	
}
	
void Solve(vector<double> &b, const FredholmMatrix& A, vector<double>& x){
	int N=A.taille();
	DenseMatrix M(N);
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++) M(i,j)=A(i,j);}
	b=ResoLU(M,x);
}


double FrobeniusNorm(const FredholmMatrix& A){
	double s=0;
	for(int i=0; i<A.size; i++){
		for(int j=0; j<A.size; j++) s+=A(i,j)*A(i,j);}	
	return sqrt(s);
}


void FredholmMatrix::SVD(const DenseMatrix & B, int r){
	int N=B.taille();
	Eigen::MatrixXf M(N,N);
	for(int i=0; i<N; i++){
		for(int j=0; j<N;j++){
			M(i,j)=B(i,j);
	}}
	REDSVD::RedSVD redsvd;
	redsvd.run(M, r);
	Eigen::MatrixXf U = redsvd.matrixU();
	Eigen::VectorXf S = redsvd.singularValues();
	Eigen::MatrixXf V = redsvd.matrixV();
	for(int i=1; i<r ;i++){
		double s=S[i];
		vector<double> u(N), v(N);
		for(int j=0; j<N; j++){
			u[j]=U(j,i);
			v[j]=V(j,i);}
			Insert(s,u,v);	
	cout<<S<<endl;
	Eigen::MatrixXf SV=U*S.asDiagonal()*V.transpose();
	cout<<"Norm "<<(SV-M).squaredNorm()<<endl;}
	
}

void FredholmMatrix::SVDEigen(const DenseMatrix & B, int r){
	int N=B.taille();
	Eigen::MatrixXd M(N,N);
	for(int i=0; i<N; i++){
		for(int j=0; j<N;j++){
			M(i,j)=B(i,j);
	}}
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::MatrixXd U = svd.matrixU();
	Eigen::VectorXd S = svd.singularValues();
	Eigen::MatrixXd V = svd.matrixV();
	for(int i=0; i<r ;i++){
		double s=S[i];
		vector<double> u(N), v(N);
		for(int j=0; j<N; j++){
			u[j]=U(j,i);
			v[j]=V(j,i);}
		Insert(s,u,v);	}

}

void FredholmMatrix::CrossApprox (const DenseMatrix & B, int r){
	DenseMatrix R(B);
	int N=R.taille();
	int k=0;
	int j=0;
	int D=sigma.size();
	for(int i=D; i<r; i++){
		vector<double> u(N), v(N);
		double s=0;
		j=R.Argmax()[0];
		k=R.Argmax()[1];
		for(int p=0; p<N; p++){
			u[p]=R(p,k);
			v[p]=R(j,p);
			if(R(j,k)!=0){s=1/R(j,k);}
			else{s=0;}
		}
			Insert(s,u,v);
		for(int a=0; a<N; a++){
			for(int b=0; b<N; b++){
				R(a,b)-=s*u[a]*v[b];
				}}
			}
}


void FredholmMatrix::CrossApproxCompute(int r){
	int N=size;
	int k=0;
	int j=0;
	int D=sigma.size();
	for(int i=D; i<=r; i++){
		vector<double> u(N), v(N);
		vector<int> ARG(2);
		double s=0;
		ARG=diag.Argmax();
		j=ARG[0];
		k=ARG[1];
		
		for(int p=0; p<N; p++){
			u[p]=diag(p,k);
			v[p]=diag(j,p);};
			if(diag(j,k)!=(double)0.0){s=1.0/(diag(j,k));}
			
			else{s=0;}				
			Insert(s,u,v);
		for(int a=0; a<N; a++){
			for(int b=0; b<N; b++){
				Insert(a,b,(-1)*s*u[a]*v[b]);
			}
		}
	}
}

double FHMDiff(const FredholmMatrix& F,const SparseMatrix &S, int p){
	double s=0;
	for(int i=0; i<S.taille(); i++){
		for(int j=0; j<S.taille(); j++) s+=pow(F.CoefSVD(i,j,p)-S(i,j), 2);}
		return(sqrt(s));
	}


void TestSVDCross(int N){
	cout<<"Test de l'algorithme d'approximation en croix "<<endl;
	ofstream F("dataCross.dat");
	DenseMatrix  M(N,'a');
	SparseMatrix S(M);
	FredholmMatrix FM1(S);
	for(int k=0; k<=15 && k<=N; k+=1){
		cout<<k;
		FM1.CrossApproxCompute(k);
		cout<<" DIFF	"<<FHMDiff(FM1,S,k)<<endl;	
		F<<k<<"	"<<(FHMDiff(FM1,S,k))<<endl;
}
	F.close();
}

void TestSVDEigen(int N){
	cout<<"Test de l'algorithme de la bibliothèque Eigen "<<endl;
	ofstream F("dataEigen.dat");
	DenseMatrix  M(N, 'a');
	FredholmMatrix FM1(N);
	SparseMatrix S(M);
	FM1.SVDEigen(M,N);
	for(int k=0; k<=15 && k<=N; k+=1){
		cout<<k<<" DIFF	"<<FHMDiff(FM1,S,k)<<endl;
		F<<k<<"	"<<(FHMDiff(FM1,S,k))<<endl;
	}
	F.close();
}
