int main(int argc, char const *argv[])
{
	double Epthron = 10e-12;
	double c, b;
	double Miu;
	double a[501];
	double DET_A;
	Matrix* A = new Matrix(5, 501); 

	Matrix* A_u = new Matrix(501, 501);


	A->setR(2);
	A->setS(2);
	//set value for A
	for (int j = 0; j < 501; j++){
		a[j] = (1.64-0.024*(j+1))*sin(0.2*(j+1))-0.64*exp(0.1/(j+1));
	}
	//init the MATRIX A
	for (int j = 0; j < 501; j++){
		if (j > 1) {
			MEle[0][j] = c;	
		}
		if (j > 0){
			MEle[1][j] = b;
		}
		if (j < 500){
			MEle[3][j] = b;
		}
		if (j < 499){
			MEle[4][j] = c;
		}
		MEle[3][j] = a[j];
	}
	
	//反幂法求lambda_s
	InversePowerMethod(A, Epthron);

	//幂法求模最大的
	PowerMethod(A, Epthron);

	//通过平移，求出另一个

	
	//cond使用 lambda1/lambda2
	//cond(A^-1) use the property of Eigenvalues(change from A's)

	//通过平移，求Miu
	// Miu = 


	//行列式 multiple all of the elements in diagonal
	DET_A = a[0];
	for (int i = 1; i < 501; i++){
		DET_A*=a[i];
	}

	return 0;
}