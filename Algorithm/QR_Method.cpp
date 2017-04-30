#include <stdio.h>
#include <iostream>
#include "Matrix.cpp"
#include <cmath>
#include "QR.cpp"
#define VREle V_r->Elements
void SetDimention(Matrix* Origin, Matrix* Target){
	int n = Origin->getColumn();
	Target->setRow(n);
	Target->setColumn(n);
}

void MatrixSquared(Matrix* A, Matrix* Result){
	int n = A->getColumn();
	double sum = 0;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			for (int k = 0; k < n; k++){
				sum+=MEle[i][k]*MEle[k][j];
			}
			REle[i][j] = sum;
			sum = 0;
		}
	}
}

void Get2Eigens(double a11, double a12, double a21, double a22, double TwoEigens[][2]){
	double b = -1*(a11+a22);
	double c = a11*a22-a21*a12;
	double deta = b*b-4*c;

	if (deta == 0) {
		TwoEigens[0][0] = (-1*b+sqrt(deta))/2;
		TwoEigens[1][0] = 0;
		TwoEigens[0][1] = TwoEigens[0][0];
		TwoEigens[1][1] = 0;
	}
	else if (deta > 0){

		TwoEigens[0][0] = (-1*b+sqrt(deta))/2;
		TwoEigens[0][1] = (-1*b-sqrt(deta))/2;
		TwoEigens[1][0] = TwoEigens[1][1] = 0;
	}
	else {
		TwoEigens[0][0] = TwoEigens[0][1] = -1*b/2;
		TwoEigens[1][0] = sqrt(-1*deta)/2;
		TwoEigens[1][1] = -1*sqrt(-1*deta)/2;
	}
}
void Get2Eigens(Matrix* A, double TwoEigens[][2]){
	int n = A->getColumn();
	Get2Eigens(MEle[n-2][n-2], MEle[n-1][n-2], MEle[n-2][n-1], MEle[n-1][n-1], TwoEigens);
}
void QRMethod(Matrix* A, double Epthron)
{
	int IterationTimes = 1000000;
	int k = 1;
	int p = 0;
	int success = 1;
	int n = A->getColumn();
	int OriginalDimention = n;
	double t_r;
	double c_r;
	double h_r;
	double Norm;
	Vector* W_r = new Vector(n);
	Vector* P_r = new Vector(n);
	Vector* U_r = new Vector(n);
	Vector* V_r = new Vector(n);
	Vector* Q_r = new Vector(n);
	double s, t;
	double Eigen[2][n];
	double TwoEigens[2][2];
	Matrix* I = new Matrix(n, n);
	Matrix* R = new Matrix(n, n);
	Matrix* R_t = new Matrix(n, n);
	Matrix* A_t = new Matrix(n, n);

	// Quasi_Triangulation(A);

	// A->printMatrix();
	// cout << '\n';
	while (k < IterationTimes && success != 0){
		while(1){
			if (abs(MEle[A->getRow()-1][A->getColumn()-2]) <= Epthron) {
				Eigen[0][p] = MEle[A->getRow()-1][A->getColumn()-1];
				Eigen[1][p] = 0;
				p++;
				A->DownGrade();

				n = A->getColumn();
				if (n == 1) {
					Eigen[0][p] = MEle[n-1][n-1];
					Eigen[1][p] = 0;
					p++;
					//end all the loop
					success = 0;
					break;
				}
				else if (n == 2){
					Get2Eigens(A, TwoEigens);
					Eigen[0][p] = TwoEigens[0][0];
					Eigen[1][p] = TwoEigens[1][0];
					Eigen[0][p+1] = TwoEigens[0][1];
					Eigen[1][p+1] = TwoEigens[1][1];
					p+=2;
					//end all the loop
					success = 0;
					break;
				}
				else {
					continue;
				}

			}
			//5
			else if(abs(MEle[A->getRow()-2][A->getColumn()-3]) <= Epthron){
				Get2Eigens(MEle[A->getRow()-2][A->getColumn()-2], MEle[A->getRow()-2][A->getColumn()-1], MEle[A->getRow()-1][A->getColumn()-2], MEle[A->getRow()-1][A->getColumn()-1], TwoEigens);
				Eigen[0][p] = TwoEigens[0][0];
				Eigen[1][p] = TwoEigens[1][0];
				Eigen[0][p+1] = TwoEigens[0][1];
				Eigen[1][p+1] = TwoEigens[1][1];
				p+=2;
				A->DownGrade();
				A->DownGrade();
				
				n = A->getColumn();
				if (n == 1) {
					Eigen[0][p] = MEle[n-1][n-1];
					Eigen[1][p] = 0;
					p++;
					//end all the loop
					success = 0;
					break;
				}
				else if (n == 2){
					//Get2Eigens(A, TwoEigens);
					Eigen[0][p] = TwoEigens[0][0];
					Eigen[1][p] = TwoEigens[1][0];
					Eigen[0][p+1] = TwoEigens[0][1];
					Eigen[1][p+1] = TwoEigens[1][1];
					p+=2;
					//end all the loop
					success = 0;
					break;
				}
				else {
					continue;
				}

			}
			else{
				break;
			}

		}
		if (success == 0) {
			break;
		}
		n = A->getColumn();
		s = MEle[n-2][n-2] + MEle[n-1][n-1];
		t = MEle[n-2][n-2]*MEle[n-1][n-1]-MEle[n-1][n-2]*MEle[n-2][n-1];	
		SetDimention(A, R);
		SetDimention(A, I);
		W_r->setSize(n);
		P_r->setSize(n);
		U_r->setSize(n);
		V_r->setSize(n);
		Q_r->setSize(n);
		MatrixSquared(A, R);
		ScalarMulti(A, s);
		Identify(I);
		ScalarMulti(I, -1*t);
		MatrixSubstract(R, I);
		MatrixSubstract(R, A);
		ScalarMulti(A, 1/s);

		for (int r = 0; r < n-1; r++){
			Norm = 0.0;
			for (int i = r+1; i < n; i++){
				Norm+=MREle[i][r]*MREle[i][r];
			}
			if (Norm == 0.0) {	
				continue;
			}
			else{
				Norm+=MREle[r][r]*MREle[r][r];
				Norm = sqrt(Norm);

				c_r = -1*sgn(MREle[r][r])*Norm;
				h_r = c_r*(c_r - MREle[r][r]);
				//generate U
				for (int i = 0; i < n; i++){
					if (i >= r+1) {
						UREle[i] = MREle[i][r];
					}
					else if (i == r) {
						UREle[i] = MREle[i][r] - c_r;
					}
					else{
						UREle[i] = 0;
					}
				}
				
				Transpose(R, R_t);
				MatrixVectorMulti(R_t, U_r, V_r);
				ScalarMulti(V_r, 1/h_r);
				for (int i = 0; i < n; i++){
					for (int j = 0; j < n; j++){
						MREle[i][j]-=UREle[i]*VREle[j];
					}
				}

				Transpose(A, A_t);
				MatrixVectorMulti(A_t, U_r, P_r); 
				ScalarMulti(P_r, 1/h_r);
				MatrixVectorMulti(A, U_r, Q_r);
				ScalarMulti(Q_r, 1/h_r);

				t_r = DotProduct(P_r, U_r)/h_r;
				ScalarMulti(U_r, t_r);
				VectorSubstract(Q_r, U_r, W_r);
				ScalarMulti(U_r, 1/t_r);

				for (int i = 0; i < n; i++){
					for (int j = 0; j < n; j++){
						MEle[i][j]-=WEle[i]*UREle[j];
						MEle[i][j]-=UREle[i]*PEle[j];
					}
				}
			}
		}

		k++;
	}
	if (success == 0) {
		for (int i = 0; i < OriginalDimention; ++i){
			cout << Eigen[0][i] << "+i*" << Eigen[1][i] << '\n';
		}
	}
	else {
		cout << "cannot get all the eigens";
	}

}

int main5()
{
	double Epthron = 10e-12;
	Matrix* A = new Matrix(10,10);
	Matrix* A_r = new Matrix(10,10);
	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < 10; j++) {
			if (i == j){
				MEle[i][j] = 1.52*cos(i+1+1.2*(j+1));
			}
			else{
				MEle[i][j] = sin(0.5*(i+1)+0.2*(j+1));
			}
		}
	}

	// A->printMatrix();
	QRMethod(A, Epthron);
	//test(A, Epthron);

	return 0;
}