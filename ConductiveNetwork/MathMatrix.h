//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	MathMatrix.h
//OBJECTIVE:	A class for mathmatical operations of matrix
//AUTHOR:		Fei Han;
//E-MAIL:			fei.han@kaust.edu.sa	;
//====================================================================================

#ifndef MATHMATRIX_H
#define MATHMATRIX_H

#include<cmath>
#include<iomanip>
#include<iostream>
#include<vector>
#include<cstdlib>
using namespace std;

const double Zero = 1.0E-8;

//---------------------------------------------------------------------------
class MathMatrix
{
public:
	//-----------------------------------------------
	//Data Member
	vector<vector<double> > element;

	//-----------------------------------------------
	//Constructor
	MathMatrix(){};
	MathMatrix(const int&, const int&);							//to construct a matrix by (M*N) with zero values
	MathMatrix(const vector<double>&);						//to construct a vector in 1D
	MathMatrix(const vector<vector<double> >&);			//to construct a vector in 2D
	MathMatrix(const double *, const int&);					//to construct an array in 1D
	MathMatrix(const double *, const int&, const int&);	//to construct an array in 2D, the arguments are (&vec[0][0], m, n)
	//-----------------------------------------------
	//Member Functions
	void Evalu(const double *, const int&);						//a 1D array is assigned to a matrix
	void Evalu(const double *, const int& , const int&);	//a 2D array is assigned to a matrix, the arguments are (&vec[0][0], m, n)
	void operator=(const vector<double>&);					//to overload the assignment operator (a 1D vector is assigned to a matrix)
	void operator=(const vector<vector<double> >&);	//to overload the assignment operator (a 2D vector is assigned to a matrix)
	MathMatrix operator*(const MathMatrix&);				//to overload the multiplication operator (a matrix is multipled by a matrix)
	MathMatrix operator*(const double&);						//to overload the multiplication operator (a matrix is multipled by a real)
	MathMatrix operator+(const MathMatrix&);				//to overload the addition operator (add a matrix to a matrix)
	double operator+(const double&);								//to overload the addition operator (add (1*1) matrix to a real)	
	MathMatrix operator-(const MathMatrix&);				//to overload the subtraction operator (subtract a matrix from a matrix)
	MathMatrix Inverse(void);											//matrix inversion
	MathMatrix Transpose(void);										//matrix transposition
	int RowN(void);														//to calculate the number of row in a matrix
	int CalN(void);															//to calculate the number of column in a matrix
	int Symmetry(void);													//to judge if a matrix is symmetric
	//-----------------------------------------------------------------------------
	MathMatrix GetCol(int cn);										//to covert a column of matrix in a new matrix

	//overload the operator of output stream
	friend ostream& operator<<(ostream& o, const MathMatrix& matrix);
};
//overload the operator of output stream
inline ostream& operator<<(ostream& o,const MathMatrix& matrix)
{
	for (int i=0; i<int(matrix.element.size()); i++)
	{
		for(int j=0; j<int(matrix.element[i].size()); j++)
		{
			if (fabs(matrix.element[i][j])<1.0E-8)
			{
				o <<setw(8) << setfill(' ') << 0 << " ";
			}
			else
			{
				o <<setw(8) << setfill(' ') << matrix.element[i][j] << " ";
			}
		}
		o << endl;
	}
	return o;
}
//---------------------------------------------------------------------------
#endif
//===========================================================================
