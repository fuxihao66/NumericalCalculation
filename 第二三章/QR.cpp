#define UEle U_r->Elements
#define REle Result->Elements
#define FEle Faciend->Elements 
#define SUEle Subtrahend->Elements
#define MIEle Minuend->Elements
#define M1Ele A1->Elements
#define WEle W_r->Elements
#define PEle P_r->Elements
#define QEle Q->Elements
void Identify(Matrix* A){
	int n = A->getColumn();
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; j++){
			if (i != j){
				MEle[i][j] = 0;
			}
			else{
				MEle[i][j] = 1;
			}
		}
	}
}
void MatrixVectorMulti(Matrix* A, Vector* Faciend, Vector* Result){
	int n = A->getColumn();
	for (int i = 0; i < n; ++i)
	{
		REle[i] = 0;
		for (int j = 0; j < n; j++){
			REle[i]+=MEle[i][j]*FEle[j]
		}
	}
}
void ScalarMulti(Vector* Faciend, double coefficient){
	int n = Faciend->getSize();
	for (int i = 0; i < n; i++){
		FEle[i]*=coefficient;
	}
}
void VectorSubstract(Vector* Subtrahend, Vector* Minuend, Vector* Result){
	int n = Subtrahend->getSize();
	for (int i = 0; i < n; i++){
		REle[i] = SUEle[i] - MIEle[i];
	}
}
int sgn(double x){
	if (x <= 0) {
		return -1;
	}
	else{
		return 1;
	}
}
void Transpose(Matrix* A, Matrix* A_T){
	int n = A->getColumn();
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			A_T->Elements[j][i] = MEle[i][j];
		}
	}
}
void Quasi_Triangulation(Matrix* A){
	int n = A->getColumn();
	double Norm;
	double c_r;
	double h_r;
	double t_r;
	Vector U_r = new Vector(n);
	Vector P_r = new Vector(n);
	Vector Q_r = new Vector(n);
	Vector W_r = new Vector(n);
	// Matrix_A 转置
	Matrix* A_T = new Matrix(n,n);
	Transpose(A,A_T);

 	for (int r = 0; r < n-2; r++)
	{
		Norm = 0.0;
		for (int i = r+2; i < n; ++i)
		{
			if (MEle[i][r] != 0)
			{
				for (int j = r+1; j < n; j++){
					Norm += MEle[i][r]*MEle[i][r];
				}
				Norm = sqrt(Norm);
				break;
			}

		}

		if (Norm == 0.0) {	
			continue;
		}
		else {
			c_r = -1*sgn(MEle[r+1][r])*Norm;
			h_r = c_r*(c_r - MEle[r+1][r]);
			//generate U
			for (int i = 0; i < n; i++){
				if (i >= r+2) {
					UEle[i] = MEle[i][r];
				}
				else if (i == r+1) {
					UEle[i] = MEle[i][r] - c_r;
				}
				else{
					UEle[i] = 0;
				}
			}

			MatrixVectorMulti(A_T, U_r, P_r);
			ScalarMulti(P_r, 1/h_r);

			MatrixVectorMulti(A, U_r, Q_r);
			ScalarMulti(Q_r, 1/h_r);

			t_r = DotProduct(P_r, U_r);
			t_r/=h_r;

			ScalarMulti(U_r, t_r);
			VectorSubstract(Q_r, U_r, W_r);
			ScalarMulti(U_r, 1/t_r);

			for (int i = 0; i < n; i++){
				for (int j = 0; j < n; j++){
					MEle-= WEle[i]*UEle[j];
					MEle-= UEle[i]*PEle[j];
				}
			}
		}

		
	}
}

void QR_Division(Matrix* A, Matrix* Q, Matrix* R){
	int n = A->getColumn();
	double Norm;
	double c_r;
	double h_r;
	Vector W_r = new Vector(n);
	Vector P_r = new Vector(n);
	Vector U_r = new Vector(n);
	Matrix* A_T = new Matrix(n, n);
	Transpose(A, A_T);
  	Identify(Q);
	for (int r = 0; r < n-1; r++){
		Norm = 0.0;
		for (int i = r+1; i < n; ++i)
		{
			if (MEle[i][r] != 0)
			{
				for (int j = r; j < n; j++){
					Norm += MEle[i][r]*MEle[i][r];
				}
				Norm = sqrt(Norm);
				break;
			}

		}

		//
		if (Norm == 0.0) {	
			continue;
		}
		else{

			c_r = -1*sgn(MEle[r][r])*Norm;
			h_r = c_r*(c_r - MEle[r][r]);
			//generate U
			for (int i = 0; i < n; i++){
				if (i >= r+1) {
					UEle[i] = MEle[i][r];
				}
				else if (i == r) {
					UEle[i] = MEle[i][r] - c_r;
				}
				else{
					UEle[i] = 0;
				}
			}

			MatrixVectorMulti(Q, U_r, W_r);
			ScalarMulti(W_r, 1/h_r);
			
			MatrixVectorMulti(A_T, U_r, P_r);
			ScalarMulti(P_r, 1/h_r);

			for (int i = 0; i < n; i++){
				for (int j = 0; j < n; j++){
					QEle[i][j]-= WEle[i]UEle[j];
					MEle[i][j]-= UEle[i]PEle[j];
				}
			}
		}
	}
}

