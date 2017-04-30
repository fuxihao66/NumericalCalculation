#include "Matrix.cpp"
#define YEle y->Elements
#define UEle u->Elements

double DotProduct(Vector* X, Vector* Y){
	double product = 0;
	for (int i = 0; i < X->getSize(); i++){
		product+=(X->Elements[i]*Y->Elements[i]);
	}
	return product;
}
double PowerMethod(Matrix* A, double Error){
	double Yita;
	double Beta;
	double FormerBeta;
	double LocalError;
	int column = A->getColumn();
	int row = A->getRow();
	//the s of Matrix_A
	int s = A->getS();
	Vector* y = new Vector(column);
	Vector* u = new Vector(column);

	//generate u0
	for (int i = 0; i < column; i++){
		UEle[i] = 0.1;
	}

	//begin loop, until the error is smaller the given one
	do {
		Yita = u->getNorm2();
		for (int i = 0; i < u->getSize(); i++){
			YEle[i] = UEle[i]/Yita;
		}

		//modify u
		if (row == column){
			for (int i = 0; i < row; i++){
				UEle[i] = 0;
				for (int j = 0; j < column; j++){
					UEle[i] += MEle[i][j]*YEle[j];
				}
			}
		}
		else {
			for (int i = 0; i < column; i++){
				UEle[i] = 0;
				for (int j = 0; j < column; j++){
					if (i-j+s >= 0 && i-j+s < row)
					{
						UEle[i]+=MEle[i-j+s][j]*YEle[j];
					}
				}
			}
		}

		FormerBeta = Beta;
		Beta = DotProduct(y, u);
		LocalError = abs(Beta - FormerBeta)/abs(Beta);
	} while (LocalError > Error);

	return Beta;
}
double InversePowerMethod(Matrix* A, double Error){
	double Beta;
	double FormerBeta;
	double LocalError;
	double Yita;
	int column = A->getColumn();
	int row = A->getRow();

	Vector* y = new Vector(column);
	Vector* u = new Vector(column);

	EquationGroup* E = new EquationGroup();
    E->initEquation(A, y);
	//do the decomposition
	if (row == column){
		E->TrianglarDecomposition();
	}
	else {
		E->TrianglarDecompStrip();
	}

	//generate u0
	for (int i = 0; i < column; i++){
		UEle[i] = 0.1;
	}
	//begin loop, until the error is smaller the given one
	do {
		Yita = u->getNorm2();
		for (int i = 0; i < column; i++){
			YEle[i] = UEle[i]/Yita;
		}

		//solve equation to get u
		E->initEquation(A, y);
		if (row == column){
			E->Substitution();
		}
		else {
			E->SubstitutionStrip();
		}

		for (int i = 0; i < column; i++){
			UEle[i] = E->getSolution()->Elements[i];
		}

		FormerBeta = Beta;
		Beta = DotProduct(y, u);
		LocalError = abs(Beta - FormerBeta)/abs(Beta);
	} while (LocalError > Error);

	return Beta;
}

int main(int argc, char const *argv[])
{
	int row, column;
	int s, r;
	double eigen;
	cin >> row;
	cin >> column;
	cin >> s;
	cin >> r;
	Matrix* A = new Matrix(row-s-r, column);
	A->printMatrix();
	//A->setValue(row, column, s, r);
	//eigen = PowerMethod(A, 10e-10);
	//cout << eigen;
	return 0;
}
