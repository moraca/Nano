//===========================================================================
// MatPro.h
// 材料属性类头文件
// A class of material property
//===========================================================================

#ifndef MATPRO_H
#define MATPRO_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<cmath>
#include "Gauss.h"
#include "MathMatrix.h"
#include "Hns.h"

using namespace hns;

//---------------------------------------------------------------------------
class MatPro
{
	public:	
	
		//数据成员
		int type_val;   //三维材料性质类型，目前定义三种(0：各向同性; 1：横观各向同性; 2：正交各向异性)
		double E11, E22, E33, Nu12, Nu23, Nu13, G12, G23, G13;	//材料的弹性模量，泊松比，剪切模量
		double elas_matrix[6][6];		//Dij 弹性矩阵

		//构造函数；
		MatPro(int itype=0){ type_val = itype; }
		//---------------------------------------------
		//成员函数
		//设置各向同性材料参数的弹性模量、泊松比和剪切模量
		void set_ela_para(const double E, const double Nu);
		//设置横观各向同性材料的弹性模量、泊松比和剪切模量
		void set_ela_para(const double E1, const double E2, const double Nu1, const double Nu2, const double G2);
		//设置正交各向异性材料的弹性模量、泊松比和剪切模量
		void set_ela_para(const double iE11, const double iE22, const double iE33, const double iNu12, const double iNu23, const double iNu13, const double iG12, const double iG23, const double iG13);

		//设置各向同性材料参数的弹性模量、泊松比和剪切模量(用于长程力等效刚度)
		void set_ela_para(const double &E);
		//设置横观各向同性材料的弹性模量、泊松比和剪切模量(用于长程力等效刚度)(输入泊松比)
		void set_ela_para(const double &E1, const double &E2, const double &Nu2);
		//设置横观各向同性材料的弹性模量、泊松比和剪切模量(用于长程力等效刚度)(输入剪切模量)
		void set_ela_para_transverse(const double &E1, const double &E2, const double &G2);
		//设置正交各向异性材料的弹性模量、泊松比和剪切模量(用于长程力等效刚度)
		void set_ela_para(const double &iE1, const double &iE2, const double &iE3, const double &iNu12, const double &iNu23, const double &iNu13);
		//(局部模型)生成Dij弹性矩阵，此弹性矩阵对应的应变向量为（e11,e22,e33,2*e12,2*e23,2*e31）转置
		int Generate_local_elas_matrix();
		//(非局部模型)生成Dij弹性矩阵，此弹性矩阵对应的应变向量为（e11,e22,e33,2*e12,2*e23,2*e31）转置
		int Generate_nonlocal_elas_matrix();
		//根据Dij弹性矩阵，反求材料参数
		int Get_ele_para_by_ela_matrix();
		//根据理论公式反求系数，并验证一致性（由于要除以Horizon半径，所以Horizon半径应该大于1，避免误差过大）
		void Compare_coef_by_analysis_formula(const int mat_type, const struct Decay_Para &decay)const;
		//根据理论公式反求刚度矩阵，并验证一致性
		void Compare_matrix_by_analysis_formula(const struct Decay_Para &decay)const;

		void print()const;
};
//---------------------------------------------------------------------------
//记录衰减参数信息
struct Decay_Para
{	
	double R;			//衰减作用域半径(Horizon)
	double radius;		//衰减半径(公式参数)
	double acoe[6];	//各向异性系数
};
//------------------------------------------------
#endif
//===========================================================================
