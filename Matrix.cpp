#include "stdio.h"
#include "stdlib.h"
#include "iostream"
#include "vector"
#define MEle A->Elements
#define VEle b->Elements
#define SEle Solution->Elements
using namespace std;

class Matrix{
private:	
	int row;
	int column;

public:
	double **Elements;
	Matrix(int row, int column){
		this->row = row;
		this->column = column;
		Elements = (double**)malloc(sizeof(double*)*row);
		for (int i = 0; i < row; ++i)
		{
			Elements[i] = (double*)malloc(sizeof(double)*column);
		}
	}
	void setValue(){
		double **p = Elements;
		cout << "Input Matrix:\n";
		for (int i = 0; i < this->row; ++i)
		{
			for (int j = 0; j < this->column; ++j)
			{
				cin >> p[i][j];
			}
		}
	}
	int getSize(){
		return this->row;
	}
};
class Vector(){
private:
	int size;
public:
	double* Elements;
	Vector(int size){
		this->size = size;
		Elements = (double*)malloc(sizeof(double)*size);
	}
	void setValue(){
		double *p = Elements;
		cout << "Input Vector:\n"
		for (int i = 0; i < this->size; ++i)
		{
			cin >> p[i];
		}
	}
}
class EquationGroup{
private:
	Matrix* A;
	Vector* b;
	Vector* Solution;
public:
	void initEquation(int row, int column){
		A = new Matrix(row, column);
		A->setValue();
		b = new Vector(row);
		b->setValue();
		Solution = new Vector(row);
	}
	void Gauss(){
		int n = A->getSize();
		int tempVar;
		for (int k = 0; i <= n-2; ++i)
		{
			if (MEle[k][k] == 0)
			{
				//抛出异常
			}
			for (int i = k+1; i <= n-1; ++i)
			{
				tempVar = MEle[i][k]/MEle[k][k];
				for (int j = k+1; i <= n-1; ++i)
				{
					MEle[i][j] = MEle[i][j] - tempVar*MEle[k][j];
				}
				VEle[i] = VEle[i] - tempVar*VEle[k];
			}
		}

		SEle[n-1] = VEle[n-1]/MEle[n-1][n-1];
		tempVar = 0;
		for (int i = n-2; i >= 0; --i)
		{
			for (int j = i+1; j <= n-1; ++j)
			{
				tempVar+=(MEle[i][j]*SEle[j]);
			}
			SEle[i] = (b[i]-tempVar)/MEle[i][i];
		}
	}
	void GaussPivot(){
		int n = A->getSize();
		int tempVar;
		for (int k = 0; i <= n-2; ++i)
		{
			//选行号
			int MaxRowIndex = k;
			int MaxNumber = MEle[k][k];
			for (int p = k+1; p <=n-1; ++p)
			{
				if (MEle[p][k] > MaxNumber)
				{
					MaxRowIndex = p;
				}
			}
			//交换
			for (p = k; p <= n-1; ++p)
			{
				tempVar = MEle[k][p];
				MEle[k][p] = MEle[MaxRowIndex][p];
				MEle[MaxRowIndex][p] = tempVar;
			}

			for (int i = k+1; i <= n-1; ++i)
			{
				tempVar = MEle[i][k]/MEle[k][k];
				for (int j = k+1; i <= n-1; ++i)
				{
					MEle[i][j] = MEle[i][j] - tempVar*MEle[k][j];
				}
				VEle[i] = VEle[i] - tempVar*VEle[k];
			}
		}

		SEle[n-1] = VEle[n-1]/MEle[n-1][n-1];
		tempVar = 0;
		for (int i = n-2; i >= 0; --i)
		{
			for (int j = i+1; j <= n-1; ++j)
			{
				tempVar+=(MEle[i][j]*SEle[j]);
			}
			SEle[i] = (b[i]-tempVar)/MEle[i][i];
		}
	}
	void TriangleDivision(){

	}
	void TriangleDivisionPivot(){

	}
};






int main(int argc, char const *argv[])
{
	int row, column;
	cout << "Input the number of row and column:\n";
	cin >> row;
	cin >> column;
	EquationGroup* E = new EquationGroup();
	E->initEquation(row, column);
	E->Gauss();
	return 0;
}