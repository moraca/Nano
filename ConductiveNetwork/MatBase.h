//===========================================================================
// MatBase.cpp
// 材料库类头文件
// A class of material database
//===========================================================================

#ifndef MATBASE_H
#define MATBASE_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<vector>
#include<string>
#include "MatPro.h"
#include "Gauss.h"
#include "Mesher.h"
#include "Hns.h"
using namespace hns;
//---------------------------------------------------------------------------
class MatBase
{
	public:

		//数据成员
		vector<MatPro> mats_vec;
		struct Decay_Para decay;

		//构造函数
		MatBase(){};

		//成员函数
		//生成材料数据库
		int Generate_matbase(ifstream &infile);

	private:
		//成员函数
		//输出所有材料的刚度矩阵
		void Print_stiffness_matrix(string print_name)const;
		void Print_stiffness_matrix(string print_name, MatPro &mat, int key)const;
		//预估局部模型等效刚度阵
		int Estimate_stiffness_long_range_interaction(const vector<Element> &elements, const vector<Node> &nodes, const double &grid_vol, 
								const vector<Node> &gauss, const vector<double> &weight, const double a[], MatPro &mat, const int &output_key)const;
		//读入信息一行，跳过注释行（以%开头）
		string Get_Line(ifstream &infile)const;
};
//-------------------------------------------------------
#endif
//===========================================================================
