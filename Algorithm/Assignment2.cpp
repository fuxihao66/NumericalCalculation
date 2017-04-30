#include <stdio.h>
#include <iostream>
#include <cmath>
#include "QR_Method.cpp"
int main()
{
	double Epthron = 10e-12;
	Matrix* A = new Matrix(10,10);
	Matrix* A_r = new Matrix(10,10);
	Matrix* Q = new Matrix(10,10);
	Matrix* R = new Matrix(10,10);
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
	Recover(A, A_r);
	//A->printMatrix();
	Quasi_Triangulation(A);
	cout << "A^n-1 is:" << endl;
	A->printMatrix();
	QR_Division(A,Q,R);
	cout << "Q is:"<< endl;
	Q->printMatrix();
	cout << "R is:" << endl;
	R->printMatrix();
	cout << "all eigens are:" << endl;
	QRMethod(A, Epthron);
	
	return 0;
}