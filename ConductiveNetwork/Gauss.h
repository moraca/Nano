//===========================================================================
// Gauss.h
// 高斯积分点类头文件
// A class of Gauss integral points
//===========================================================================

#ifndef GAUSS_H
#define GAUSS_H

#include<iostream>
#include<vector>
#include<cmath>
#include "Hns.h"
#include "Fem_3D.h"
using namespace hns;

const double PI = 3.1415926535897932384626433832795;

//----------------------------------------------------------

class Gauss
{
public:
	//数据变量
	vector<Node> gauss; 
    vector<double> weight;

    //构造函数；
	Gauss(int ipre=3){ precision=ipre; }
	
	//成员函数
	//生成三维高斯节点
	int Generate_gauss(ifstream &infile);

private:

	//数据变量；
	int precision;

	//用于生成高斯点序列
	int Generate_gauss_array(vector<double> &gauss, vector<double> &weight)const;
	//读入信息一行，跳过注释行（以%开头）
	string Get_Line(ifstream &infile)const;
};//-----------------------------------------------------
#endif
//===========================================================================
