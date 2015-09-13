//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Gauss.h
//OBJECTIVE:	A class definition for Gaussian quadratrure points
//AUTHOR:		Fei Han;
//E-MAIL:			fei.han@kaust.edu.sa	;
//====================================================================================

#ifndef GAUSS_H
#define GAUSS_H

#include<iostream>
#include<vector>
#include<cmath>
#include "Hns.h"
#include "Fem_3D.h"
using namespace hns;

//const double PI = 3.1415926535897932384626433832795;

//----------------------------------------------------------

class Gauss
{
public:
	//Data Member
	vector<Node> gauss; 
    vector<double> weight;

    //Constructor
	Gauss(int ipre=3){ precision=ipre; }
	
	//Member Functions
	//Generate 3D Gaussian quadrature points
	int Generate_gauss(ifstream &infile);

private:
	//Data Member
	int precision;

	//Generate a sequence of Gaussian quadrature points
	int Generate_gauss_array(vector<double> &gauss, vector<double> &weight)const;

	//读入信息一行，跳过注释行（以%开头）
	string Get_Line(ifstream &infile)const;
};//-----------------------------------------------------
#endif
//===========================================================================
