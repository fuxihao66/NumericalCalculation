
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "Matrix.cpp"
#define AREle A_remain->Elements
using namespace std;
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
int main(int argc, char const *argv[])
{
	double Epthron = 10e-12;
	double c, b;
	double Miu[39];
	double a[501];
	double DET_A = 1;
	double lambda_s;
	double lambda_max;
	double lambda1;
	double lambda501;
	double cond;
	int s = 2;
	Matrix* A = new Matrix(5, 501); 
	//用于保护矩阵
	Matrix* A_remain = new Matrix(5, 501);
	A->setR(2);
	A->setS(2);

	for (int j = 0; j < 501; j++){
		a[j] = (1.64-0.024*(j+1))*sin(0.2*(j+1))-0.64*exp(0.1/(j+1));
	}
	b = 0.16;
	c = -0.064;
	//init the MATRIX A
	for (int j = 0; j < 501; j++){
		if (j > 1) {
			MEle[0][j]  = c;	
			AREle[0][j] = c;
		}
		if (j > 0){
			MEle[1][j]  = b;
			AREle[1][j] = b;
		}
		if (j < 500){
			MEle[3][j] 	= b;
			AREle[3][j] = b;
		}
		if (j < 499){
			MEle[4][j] 	= c;
			AREle[4][j] = c;
		}
		MEle[2][j] 	= a[j];
		AREle[2][j] = a[j];
	}

	//get the one with max module
	lambda_max = PowerMethod(A, Epthron);

	//translate to get the other
	A->Translate(abs(lambda_max));
	lambda1 = InversePowerMethod(A, Epthron) + abs(lambda_max);

	//compare lambda1 and lambda501
	if (lambda1 > lambda_max){
		lambda501 = lambda1;
		lambda1 = lambda_max;
	}
	else {
		lambda501 = lambda_max;
	}

	//get lambda_s
	Recover(A_remain, A);
	lambda_s = InversePowerMethod(A, Epthron);

	//A stands for L+U, thus the products of the diagonal is detA
	for (int i = 0; i < 501; ++i)
	{
		DET_A*=MEle[s][i];
	}

	cond = abs(lambda_max/lambda_s);

	//translate to get those close to Miu
	for (int i = 0; i < 39; ++i)
	{
		Miu[i] = lambda1 + (i+1)*(lambda501 - lambda1)/40;
		Recover(A_remain, A);
		//translate
		A->Translate(Miu[i]);

		cout << "lambda_miu" << i+1 << " = "<< setprecision(11) << scientific << InversePowerMethod(A, Epthron)+Miu[i] << endl;
	}

	cout << "lambda501 = " << lambda501 << endl;
	cout << "lambda1 = " << lambda1 << endl;
	cout << "lambda_s = " << lambda_s << endl;
	cout << "cond = " << cond << endl;
	cout << "det_A = " <<DET_A << endl;
	
	return 0;
}