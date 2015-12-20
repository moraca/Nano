//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	MathMatrix.cpp
//OBJECTIVE:	A class for mathmatical operations of matrix
//AUTHOR:		Fei Han;
//E-MAIL:			fei.han@kaust.edu.sa	;
//====================================================================================

#include "MathMatrix.h"

//---------------------------------------------------------------------------
//Constructor

//to construct a matrix by (M*N) with zero values
MathMatrix::MathMatrix(const int& M, const int& N)
{
	vector<double> temp(N,0.0);
	element.assign(M,temp);
}
MathMatrix::MathMatrix(const vector<double>& vec2D)
{
	int N = int(vec2D.size());
	vector<double> temp(N);
	element.assign(1,temp);
	for (int i=0; i<N; i++)
	{
		element[0][i] = vec2D[i];
	}
}
MathMatrix::MathMatrix(const vector<vector<double> >& vec2D)
{
	element = vec2D;
}
//to construct a matrix by a (1*N) vector
MathMatrix::MathMatrix(const double *vec, const int& N)
{
	if (N<=0)
	{
		cout << "Error, the dimension of a 1D vector is less than 1, i.e., N<=0." << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(1,temp);
	for (int i=0; i<N; i++)
	{
		element[0][i] = vec[i];
	}
}
//to construct a matrix by a (M*N) vector
MathMatrix::MathMatrix(const double *vec, const int& M, const int& N)
{
	if (M<=0 || N<=0)
	{
		cout << "Error, one of dimensions of a 2D vector is less than 1, i.e., M<=0 or N<=0." << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(M,temp);
	if (N!=1)
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[M*i+j];
			}
		}
	}
	else 
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[i];
			}
		}
	}
}
//---------------------------------------------------------------------------
//A 1D array is assigned to a matrix
void MathMatrix::Evalu(const double *vec, const int& N)
{
	if (N<=0)
	{
		cout << "Error, the dimension of a 1D vector is less than 1, i.e., N<=0." << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(1,temp);
	for (int i=0; i<N; i++)
	{
		element[0][i] = vec[i];
	}
}
//---------------------------------------------------------------------------
//A 2D array is assigned to a matrix
void MathMatrix::Evalu(const double *vec, const int& M, const int& N)
{
	if (M<=0 || N<=0)
	{
		cout << "Error, one of dimensions of a 2D vector is less than 1, i.e., M<=0 or N<=0." << endl;
		exit (0);
	}
	vector<double> temp(N);
	element.assign(M,temp);
	if (N!=1)
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[M*i+j];
			}
		}
	}
	else 
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				element[i][j] = vec[i];
			}
		}
	}
}
//---------------------------------------------------------------------------
//To overload the assignment operator (a 1D vector is assigned to a matrix)
void MathMatrix::operator=(const vector<double>& vec)
{
	element.assign(1,vec);
}
//To overload the assignment operator (a 2D vector is assigned to a matrix)
void MathMatrix::operator=(const vector<vector<double> >& vec)
{
	element = vec;
}
//---------------------------------------------------------------------------
//To overload the multiplication operator (a matrix is multipled by a matrix)
MathMatrix MathMatrix::operator*(const MathMatrix& matrix)
{
	int L = 0;
	if (element[0].size() != matrix.element.size())
	{
		cout << "Error, it fails to implement multiplication between two matrices, the number of cloumn is not equal to the number of row!" << endl;
		exit(0);
	}
	else
	{
		L = int(matrix.element.size());
	}
	//The definition of a matrix: M rows and N column
	int M = int(element.size());
	int N = int(matrix.element[0].size());
	//A temporary matrix: Tem_mat
	MathMatrix Tem_mat(M,N);
	//An accumulator
	double conter;

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			conter = 0.0;
			for(int k=0; k<L; k++)
			{
				conter = conter + element[i][k] * matrix.element[k][j];
			}
			Tem_mat.element[i][j] = conter;
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//To overload the multiplication operator (a matrix is multipled by a real)
MathMatrix MathMatrix::operator*(const double& R)
{
	//The definition of a matrix: M rows and N column
	int M = int(element.size());
	int N = int(element[0].size());
	//A temporary matrix: Tem_mat
	MathMatrix Tem_mat(M,N);

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			Tem_mat.element[i][j] = element[i][j] * R;
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//To overload the addition operator (add a matrix to a matrix)
MathMatrix MathMatrix::operator+(const MathMatrix& matrix)
{
	if ((element.size() != matrix.element.size())||(element[0].size() != matrix.element[0].size()))
	{
		cout << "Error, it fails to implement addition between two matrices, the number of cloumn is not equal to the number of row!" << endl;
		exit(0);
	}
	//The definition of a matrix: M rows and N column
	int M = int(element.size());
	int N = int(element[0].size());
	//A temporary matrix: Tem_mat
	MathMatrix Tem_mat(M,N);

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			Tem_mat.element[i][j] = element[i][j] + matrix.element[i][j];
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//To overload the addition operator (add (1*1) matrix to a real)	
double MathMatrix::operator+(const double& R)
{
	//The definition of a matrix: M rows and N column
	int M = int(element.size());
	int N = int(element[0].size());
	if (M!=N||M==1)
	{
		cout << "Error, the dimensions of this matrix is not (1*1), so it cannot be added to a real!" << endl;
		exit(0);
	}
	//The definition of a materix
	double R2 = element[0][0] + R;

	return R2;
}
//---------------------------------------------------------------------------
//To overload the subtraction operator (subtract a matrix from a matrix)
MathMatrix MathMatrix::operator-(const MathMatrix& matrix)
{
	if ((element.size() != matrix.element.size())||(element[0].size() != matrix.element[0].size()))
	{
		cout << "Error, it fails to implement subtraction between two matrices, the number of cloumn is not equal to the number of row!" << endl;
		exit(0);
	}
	//The definition of a matrix: M rows and N column
	int M = int(element.size());
	int N = int(element[0].size());
	//A temporary matrix: Tem_mat
	MathMatrix Tem_mat(M,N);

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			Tem_mat.element[i][j] = element[i][j] - matrix.element[i][j];
		}
	}

	return Tem_mat;
}
//---------------------------------------------------------------------------
//to judge if a matrix is symmetric (-1: is not square matrix, 0: is a square matrix but nonsymmetric, 1: is a square matrix but symmetric
int MathMatrix::Symmetry(void)
{
	int M = int(element.size());
	int N = int(element[0].size());

	if (M==N)
	{
		for (int i=0; i<M; i++)
		{
			for (int j=i+1; j<M; j++)
			{
				if (fabs(element[i][j]-element[j][i])>Zero)
				{
					return 0;
				}
			}
		}
	}
	else
	{
		return -1;
	}
	return 1;
}
//---------------------------------------------------------------------------
//matrix inversion
MathMatrix MathMatrix::Inverse(void)
{
	vector<vector<double> > vec = element;
	//--------------------------------------------------
	//to judge if this matrix is a positive definite matrix

	//--------------------------------------------------
	//performs matrix inversion
	int stRank = int(vec.size());
	vector<double> b(stRank,0.0);
	for(int k=0; k<stRank; k++)
	{
		double w= vec[0][0];
		int m = stRank - k -1;
		for(int i=1; i<stRank; i++)
		{
			double g = vec[i][0];
			b[i] = g / w;
			if (i<=m)
			{
				b[i] = -b[i];
			}
			for(int j=1; j<=i; j++)
			{
				vec[i-1][j-1] = vec[i][j] + g * b[j];
			}
		}
		vec[stRank-1][stRank-1] = 1.0 / w;
		for(int i= 1; i<stRank; i++)
		{
			vec[stRank-1][i-1] =  b[i];
		}
	}
	for(int i=0; i<stRank-1; i++)
	{
		for(int j = i+1; j<stRank; j++)
		{
			vec[i][j] = vec[j][i];
		}
	}

	MathMatrix matrix = vec;
	return matrix;
}
//---------------------------------------------------------------------------
//matrix transposition
MathMatrix MathMatrix::Transpose(void)
{
	int M = int(element.size());
	int N = int(element[0].size());

	vector<vector<double> > vec;
	vector<double> temp(M);
	vec.assign(N,temp);

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<M; j++)
		{
			vec[i][j] = element[j][i];
		}
	}

	MathMatrix matrix = vec;
	return matrix;
}
//---------------------------------------------------------------------------
//To covert a column of matrix in a new matrix
MathMatrix MathMatrix::GetCol(int cn){
	MathMatrix matrix(int(element.size()),1);
	for(int i=0;i<int(element.size());i++)
		matrix.element[i][0]=element[i][cn];
	return matrix;
}                  
//---------------------------------------------------------------------------
//To calculate the number of row in a matrix
int MathMatrix::RowN(void)
{
	return int(element.size());
}
//---------------------------------------------------------------------------
//To calculate the number of column in a matrix
int MathMatrix::CalN(void)
{
	return int(element[0].size());
}
//============================================================================
