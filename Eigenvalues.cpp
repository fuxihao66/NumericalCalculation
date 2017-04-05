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
	double Beta;
	double FormerBeta;
	double LocalError;
	double n = A->getColumn();
	Vector* y = new Vector(n);
	Vector* u = new Vector(n);
	double Yita;
	int s = A->getS();
	//generate u0
	for (int i = 0; i < n; i++){
		UEle[i] = 0.1;
	}
	//begin loop, until the error is smaller the given one
	do {
		Yita = u->getNorm2();
		for (int i = 0; i < u->getSize(); i++){
			YEle[i] = UEle[i]/Yita;
		}

		//modify u
		if (A->getRow() == A->getColumn()){
			for (int i = 0; i < n; i++){
				UEle[i] = 0;
				for (int j = 0; j < ; j++){
					UEle[i] += MEle[i][j]*YEle[j];
				}
			}
		}
		else {
			for (int i = 0; i < n; i++){
				UEle[i] = 0;
				for (int j = 0; j < n; j++){
					if (i-j+s+1 >= 0 && i-j+s-1 < n)
					{
						UEle[i]+=MEle[i-j+s+1][j]*YEle[j];
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

	Vector* y = new Vector(A->getColumn());
	Vector* u = new Vector(A->getColumn());
	EquationGroup* E = new EquationGroup(A, y);

	//do the decomposition
	if (A->getRow() == A->getColumn()){
		E->TrianglarDecomposition();
	}
	else {
		E->TrianglarDecompStrip();
	}

	//generate u0
	for (int i = 0; i < u->getSize(); i++){
		UEle[i] = 0.1;
	}
	//begin loop, until the error is smaller the given one
	do {
		Yita = u->getNorm2();
		for (int i = 0; i < u->getSize(); i++){
			YEle[i] = UEle[i]/Yita;
		}

		//solve equation to get u
		E->initEquation(A, y);
		if (A->getRow() == A->getColumn()){
			// E->
		}
		else {
			E->SubstitutionStrip();
		}

		for (int i = 0; i < u->getSize(); i++){
			UEle[i] = E->getSolution()[i];
		}

		FormerBeta = Beta;
		Beta = DotProduct(y, u);
		LocalError = abs(Beta - FormerBeta)/abs(Beta);
	} while (LocalError > Error);
	
	return Beta;	
}