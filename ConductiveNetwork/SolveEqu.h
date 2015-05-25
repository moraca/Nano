//===========================================================================
// SolveEqu.h
// 解线性方程组头文件
// A class to solving linear equations
//===========================================================================
#ifndef SOVLEEQU_H
#define SOVLEEQU_H

#include<assert.h>
#include<iostream>
#include<cmath>
#include<vector>
#include "Geometry_3D.h"
#include "Fem_3D.h"
#include "Hns.h"
using namespace hns;

//---------------------------------------------------------------------------
//定义解方程类
class SolveEqu
{
public:
	//紧缩存储刚度矩阵(紧缩存储节点对应关系)
	int izig(const vector<Node> &nodes, vector<long> &Iz, vector<int> &Ig)const;

	//处理周期边界条件约束
	void Periodical_boundary_constraints(const vector<Node> &nodes, const vector<int>* &peri_bnods, vector<long> &Iz, vector<int> &Ig, vector<double> &AK,  vector<double>* F)const;

	//设置固定位移约束
	int Fixed_displacement_constraints(ifstream &infile, const vector<Node> &nodes, int &bnod_num, vector<int> &ip, vector<double> &vp)const;
	//处理三维固定位移约束点的零位移值
	void Deal_with_displacement_zero_value(const int &bnod_num, const int &nod_size, const vector<long> &Iz, const vector<int> &Ig, const vector<int> &Ip,
																			 const vector<double> &vp, vector <double>* F, vector <double> &AK)const;
	
	//求解线性方程组
	void Solve_linear_equations(const int &bnod_num, const int &N, const vector<long> &Iz, const vector<int> &Ig, const vector<int> &ip, const vector<double> &vp, 
													  const vector<double> &A, const vector<double> &B, vector<double> &X)const;

	//根据周期边界条件附加其他节点的解
	void Add_periodical_solution(const vector<int>* &peri_bnods, const vector<Node> nodes, vector<vector<double> > &solution)const;

	//输出完全的刚度矩阵及右端项用于检测
	void Export_complete_matrix_equright(const vector<double> &A, const vector <double>* F,  const vector<Node> &nodes, const vector<long> &Iz, const vector<int> &Ig)const;

	//用完全的刚度矩阵及右端项检测周期边界条件的处理
	void Complete_matrix_equright_testing_periodical_constraints(const vector<int>* &peri_bnods, const vector<double> &A, const vector <double>* F, const vector<Node> &nodes, const vector<long> &Iz, const vector<int> &Ig)const;

	//用于输出周期边界上节点的结果用于比对
	void Compared_periodical_bounday_solutions(const int &loop_num, const vector<int>* &peri_bnods, vector<double> &solution)const;
private:
	//处理右端向量
	void deal_with_F(const int num[], const double dist[3], const vector<vector<int> > &relative_small, const vector<vector<int> > &relative_large,
								   const vector<vector<double> > &temp_AKii, const vector<vector<vector<double> > > &temp_AK, vector<double>* F)const;
	
	//处理AK[ii], AK[ij]和AK[jj]的部分
	void deal_with_AK_ii_ij_jj(const int num[], vector<vector<int> > &relative_small, vector<vector<int> > &relative_large, vector<vector<double> > &temp_AKii, vector<vector<vector<double> > > &temp_AK)const;

	//处理AK[ki]和AK[kj]的部分
	void deal_with_AK_ki_kj(const int num[], vector<vector<int> > &relative_small, vector<vector<int> > &relative_large, vector<vector<vector<double> > > &temp_AK)const;	
	
	//处理三维固定位移约束点的非零位移值
	void displacement_nonzero_value(const int &Kw, const int &Kg, const int &bnod_num, const vector<int> &Ip, const vector<double> &vp, vector<double> &W, vector<double> &G)const;

	//矩阵AK乘向量U
	void mabvm(const int &N, const int &N1, const vector<long> &Iz, const vector<int> &Ig, const vector<double> &AK, const vector<double> &U, vector<double> &V)const;

	//读入信息一行，跳过注释行（以%开头）
	string Get_Line(ifstream &infile)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
