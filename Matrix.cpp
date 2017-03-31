#include "stdio.h"
#include "stdlib.h"
#include "iostream"
#include "float.h"
#include <cmath>
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
class Vector{
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
		cout << "Input Vector:\n";
		for (int i = 0; i < this->size; ++i)
		{
			cin >> p[i];
		}
	}
    void printValues(){
        for (int i = 0; i < this->size; i++){
            cout << Elements[i] << "/";
        }
    }
};
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

	Vector* getSolution(){
		return Solution;
	}
	void Gauss(){
		int n = A->getSize();
		double tempVar;
		for (int k = 0; k <= n-2; ++k)
		{
			if (MEle[k][k] == 0)
			{
				//抛出异常
			}
			for (int i = k+1; i <= n-1; ++i)
			{
				tempVar = MEle[i][k]/MEle[k][k];
				for (int j = k+1; j <= n-1; ++j)
				{
					MEle[i][j] = MEle[i][j] - tempVar*MEle[k][j];
				}
				VEle[i] = VEle[i] - tempVar*VEle[k];
			}
		}

		SEle[n-1] = VEle[n-1]/MEle[n-1][n-1];
		
		for (int i = n-2; i >= 0; --i)
		{
            tempVar = 0;
			for (int j = i+1; j <= n-1; ++j)
			{
				tempVar+=(MEle[i][j]*SEle[j]);
			}
			SEle[i] = (VEle[i]-tempVar)/MEle[i][i];
		}
	}
	void GaussPivot(){
		int n = A->getSize();
		double tempVar;
		for (int k = 0; k <= n-2; ++k)
		{
			//选行号
			int MaxRowIndex = k;
			int MaxNumber = MEle[k][k];
			for (int p = k+1; p <=n-1; ++p)
			{
				 if (abs(MEle[p][k]) > abs(MaxNumber))
				{
					MaxRowIndex = p;
                    MaxNumber = MEle[p][k];
				}
			}
			//交换
			for (int p = k; p <= n-1; ++p)
			{
				tempVar = MEle[k][p];
				MEle[k][p] = MEle[MaxRowIndex][p];
				MEle[MaxRowIndex][p] = tempVar;                
			}
			tempVar = VEle[k];
            VEle[k] = VEle[MaxRowIndex];
            VEle[MaxRowIndex] = tempVar;
			
			for (int i = k+1; i <= n-1; ++i)
			{
				tempVar = MEle[i][k]/MEle[k][k];
				for (int j = k+1; j <= n-1; ++j)
				{
					MEle[i][j] = MEle[i][j] - tempVar*MEle[k][j];
				}
				VEle[i] = VEle[i] - tempVar*VEle[k];
			}
		}

		SEle[n-1] = VEle[n-1]/MEle[n-1][n-1];
		
		for (int i = n-2; i >= 0; --i)
		{
            tempVar = 0;
			for (int j = i+1; j <= n-1; ++j)
			{
				tempVar+=(MEle[i][j]*SEle[j]);
			}
			SEle[i] = (VEle[i]-tempVar)/MEle[i][i];
		}
	}
    
	void TrianglarDecomposition(){
		int n = A->getSize();
		double u, l;
		double y[n];       
		for(int k = 0; k < n; k++){

			for(int j = k; j < n; j++){
				u = MEle[k][j];
				for(int t = 0; t <= k-1; t++){
					u-=MEle[k][t]*MEle[t][j];
				}

				MEle[k][j] = u;
			}
			for(int i = k+1; i < n; i++){
				l = MEle[i][k];
				for(int t = 0; t <= k-1; t++){
					l-=MEle[i][t]*MEle[t][k];
				}
				l/=MEle[k][k];

				MEle[i][k] = l;
			}

		}

		y[0] = VEle[0];
		for(int i = 1; i < n; i++){
			y[i] = VEle[i];
			for(int t = 0; t <= i-1; t++){
				y[i]-=MEle[i][t]*y[t];
			}
		}
		SEle[n-1] = y[n-1]/MEle[n-1][n-1];
		for(int i = n-2; i >= 0; i--){
			SEle[i] = y[i];
			for(int t = i+1; t < n; t++){
				SEle[i]-=MEle[i][t]*SEle[t];
			}
			SEle[i]/=MEle[i][i];
		}
	}
	void TrianglarDecompPivot(){
        int n = A->getSize();
		double u, l;
		double y[n];
		double MiddleValue[n];
		double temp;
		int MaxIndex[n];
		double MaxValue = 4.9406564584124654E-323;
		
		for(int k = 0; k < n; k++){
			
			//计算中间量
			//中间量直接存到MEle[i][k]里面
			for(int i = k; i < n; i++){
				MiddleValue[i] = MEle[i][k];
				for(int t = 0; t <= k-1; t++){
					MiddleValue[i]-=MEle[i][t]*MEle[t][k];
				}	
				if (abs(MiddleValue[i]) > abs(MaxValue)){
					MaxIndex[k] = i;
					MaxValue = MiddleValue[i];
				}
			}
			
			
			//交换
			if (MaxIndex[k] != k){
				for (int t = 0; t < n; t++){
					temp = MEle[k][t];
					MEle[k][t] = MEle[MaxIndex[k]][t];
					MEle[MaxIndex[k]][t] = temp;
				}
				temp = MiddleValue[k];
				MiddleValue[k] = MiddleValue[MaxIndex[k]];
				MiddleValue[MaxIndex[k]] = temp;
			}
			
			
			
			//求u  
			MEle[k][k] = MiddleValue[k];
			for(int j = k+1; j < n; j++){
				u = MEle[k][j];
				for(int t = 0; t <= k-1; t++){
					u-=MEle[k][t]*MEle[t][j];
				}
				MEle[k][j] = u;
			}
			for(int i = k+1; i < n; i++){
				MEle[i][k] = MiddleValue[i]/MEle[k][k];
			}

		}
		
		//求Qb
		for (int k = 0; k < n-1; k++){
			temp = VEle[k];
			VEle[k] = VEle[MaxIndex[k]];
			VEle[MaxIndex[k]] = VEle[k];
		}
		
		//回带过程
		y[0] = VEle[0];
		for(int i = 1; i < n; i++){
			y[i] = VEle[i];
			for(int t = 0; t <= i-1; t++){
				y[i]-=MEle[i][t]*y[t];
			}
		}
		SEle[n-1] = y[n-1]/MEle[n-1][n-1];
		for(int i = n-2; i >= 0; i--){
			SEle[i] = y[i];
			for(int t = i+1; t < n; t++){
				SEle[i]-=MEle[i][t]*SEle[t];
			}
			SEle[i]/=MEle[i][i];
		}
	}
	
	void Jacobi(double ErrorRequirement){
		int n = A->getSize();
		double Temp[n];
		double LocalError;
		for (int k = 0; k < n; k++){
			SEle[k] = 0;
		}
		while (LocalError < ErrorRequirement){
			for (int i = 0; i < n; i++){
				Temp[i] = VEle[i];
				for (int j = 0; j < n; j++){
					if (j != i){
						Temp[i]-=MEle[i][j]*SEle[i];
					}
				}
				Temp[i] = Temp[i]/MEle[i][i];
			}
			for (int i = 0; i < n; i++){
				SEle[i] = Temp[i];
			}
			
			//计算误差，如果满足则退出循环
		}
	}
	
	void GaussSeidel(double ErrorRequirement){
		int n = A->getSize();
		int LocalError;
		for (int k = 0; k < n; k++){
			SEle[k] = 0;
		}
		while (LocalError < ErrorRequirement){
			for (int i = 0; i < n; i++){
				SEle[i] = VEle[i]/MEle[i][i];
				for (int j = 0; j < i; j++){
					SEle[i]-=(MEle[i][j]/MEle[i][i]*SEle[j]);
				}
				for (int j = i+1; j < n; j++){
					SEle[i]-=(MEle[i][j]/MEle[i][i]*SEle[j]);
				}
			}
			//计算误差，如果满足则退出循环
		}
	}

	void SOR(double error){
		
	}
	
};

int main(int argc, char const *argv[])
{
	int row, column;
	Vector* Solution;
	cout << "Input the number of row and column:\n";
	cin >> row;
	cin >> column;
	EquationGroup* E = new EquationGroup();
	E->initEquation(row, column);

	E->TrianglarDecompPivot();
	Solution = E->getSolution();
    Solution->printValues();
	return 0;
}
