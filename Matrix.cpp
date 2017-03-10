#include "stdio.h"
#include "stdlib.h"
#include "iostream"
using namespace std;

class Matrix{
private:
	double **Elements;
	int size;
public:
	Matrix(int size){
		this->size = size;
		Elements = (double**)malloc(sizeof(double*)*size);
		for (int i = 0; i < size; ++i)
		{
			Elements[i] = (double*)malloc(sizeof(double)*size);
		}
	}
	void setValue(){
		double **p = Elements;
		cout << "Input:\n";
		for (int i = 0; i < this->size; ++i)
		{
			for (int j = 0; j < this->size; ++j)
			{
				cin >> p[i][j];
			}
		}
	}
	// void printValue(){
	// 	double **p = Elements;
	// 	for (int i = 0; i < this->size; ++i)
	// 	{
	// 		for (int j = 0; j < this->size; ++j)
	// 		{
	// 			cout << p[i][j];
	// 		}
	// 	}
	// }

	void Gauss(){
		
	}

	void TriangleDivision(){

	}
};

int main(int argc, char const *argv[])
{
	int size;
	cin >> size;
	Matrix *m1 = new Matrix(size);
	m1->setValue();


	return 0;
}