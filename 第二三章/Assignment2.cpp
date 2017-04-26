#include <stdio.h>
#include <iostream>
#include "Matrix.cpp"
#include "QR.cpp"
void Recover(Matrix* A_for_copying, Matrix* A_for_targeting){
	int row = A_for_targeting->getRow();
	int column = A_for_targeting->getColumn();
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < column; ++j)
		{
			A_for_targeting->Elements[i][j] = A_for_copying->Elements[i][j];
		}
	}
}
void Get2Eigens(Matrix* A, double* TwoEigens){

}
void GetRestEigen(Matrix* A, int p, double* EigenArray){
	int n = A->getColumn();
	double TwoEigens[2][2];
	if (n == 1) {
		EigenArray[0][p] = MEle[n-1][n-1];
		EigenArray[1][p] = 0;
		p++;
		break;
	}
	else if (n == 2){
		Get2Eigens(A, TwoEigens);
		EigenArray[0][p] = TwoEigens[0][0];
		EigenArray[1][p] = TwoEigens[1][0];
		EigenArray[0][p+1] = TwoEigens[0][1];
		EigenArray[1][p+1] = TwoEigens[1][1];
		p+=2;
		break;
	}
	else continue;
}
void QRMethod(Matrix* A)
{
	int IterationTimes;
	int k = 1;;
	int p;
	int n = A->getColumn();
	double s, t;
	double Epthron;
	double Eigen[2][n];
	Matrix* Q = new Matrix(n, n);
	Matrix* R = new Matrix(n, n);
	Quasi_Triangulation(A);
	while (k < IterationTimes){
		if (MEle[A->getRow()-1][A->getColumn()-2] <= Epthron) {
			Eigen[0][p] = MEle[A->getRow()-1][A->getColumn()-1];
			Eigen[1][p] = 0;
			p++;
			A->DownGrade();
			GetRestEigen(A);		
		}
		else if(MEle[A->getRow()-2][A->getColumn()-3] <= Epthron){
			//Get2Eigen;
			GetRestEigen();
		}

		//
		s = MEle[A->getRow()-2][A->getColumn()-2] + MEle[A->getRow()-1][A->getColumn()-1];
		t = MEle[A->getRow()-2][A->getColumn()-2]*MEle[A->getRow()-1][A->getColumn()-1]-MEle[A->getRow()-1][A->getColumn()-2]*MEle[A->getRow()-2][A->getColumn()-1];
		
		k++;
	}
	return 0;
}