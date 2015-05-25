//===========================================================================
// Global_Stiff_Matrix.h
// 求解总刚阵类头文件
// A Class of the Global Stiff Matrix
//===========================================================================
#ifndef Global_Stiff_Matrix_H
#define Global_Stiff_Matrix_H

#include "Fem_3D.h"
#include "MatPro.h"
#include "Gauss.h"
#include "Mesher.h"
#include "Geometry_3D.h"
#include "time.h"
#include<omp.h>  //openmp

#define CHUNKSIZE 1		//defines the chunk size as 1 contiguous iteration
//---------------------------------------------------------------------------
class Global_Stiff_Matrix
{
	public:
	//数据变量

	//构造函数
		Global_Stiff_Matrix(){};

	//成员函数
		//生成总体刚度矩阵
		int Gen_global_stiff_matrix(ifstream &infile, const vector<Point_3D> &cnps, const vector<Node> &nodes, const vector<Element> &elements, const vector<MatPro> &mats, const struct Decay_Para &decay, 
														 const struct RVE_Geo &cell, const struct CNT_Geo &cnts, const string &com_mod, const vector<long> &Iz, const vector<int> &Ig, vector<double> &total_matrix, vector<double>* nl_equright, double backup_ege[]);

		//计算等效应变能
		int Calculate_equivalent_energy_simple(const vector<long> &Iz, const vector<int> &Ig, const vector<double> &total_matrix, const vector<double>* nl_equright, const double backup_ege[], const vector<vector<double> > &disp_solu,
																			  double equivalent_energy[])const;
		int Calculate_equivalent_energy(const vector<Element> &elements, const vector<MatPro> &mats, const struct Decay_Para &decay, const vector<vector<double> > &disp_solu, const struct RVE_Geo &cell,
																const string &com_mod, double equivalent_energy[])const;

	private:

	//成员函数
		//生成RVE背景网格单元(接触力作用)单刚阵(六面体)（CNT增强界面层局部连续模型部分）
		void CNT_Reinforcement_Local_Continuum(double (*element_stiff_matrix)[24], const vector<MatPro> &mats, const double &Jacobi, const double (*gauss_dfx)[3][8], 
																					   const vector<double> &weight, const double (*gwc_left)[10])const;
		//决定单元的相关单元, 用于计算高斯点距纳米管的最小距离
		void Deter_relative_elements_min_dist(const vector<Node> &nodes, const vector<Element> &elements, const struct RVE_Geo &cell, const double &dist, vector<vector<int> > &relemd)const;
		//生成RVE背景网格单元(接触力作用)单刚阵
		void Generate_Contactforce_Elestiff(double (*element_stiff_matrix)[24], const vector<MatPro> &mats, const double &Jacobi, const double (*gauss_dfx)[3][8], const vector<double> &weight)const;
		//生成RVE背景网格单元(长程力作用)单刚阵
		void Generate_Longforce_Elestiff(double (*element_stiff_matrix1)[24], double (*element_stiff_matrix2)[24], const double (*gauss_ns)[8], const double (*gauss_nw)[8], const vector<double> &weight,
																	 const struct Decay_Para &decay, const double (*gauss_po)[3], const double elec_left[], const double elec_right[], const double (*gwc_left)[10], const double (*gwc_right)[10])const;
		//生成边界单元(长程力作用)单刚阵(六面体)等效方程右端向量
		void Generate_Longforce_Equivalent_Equright(double (*element_equright)[24], const double &Jacobi, const double (*gauss_nw)[8], const vector<double> &weight, const struct Decay_Para &decay,
																						  const double (*gauss_po)[3], const double elec_left[], const double elec_right[], const double (*gwc_left)[10], const double (*gwc_right)[10], const double peri_dist[], double ege[])const;
		//将单刚添加到总刚(行列单元节点相同)
		void Add_to_gsmatrix(const double (*element_stiff_matrix)[24], const vector<long> &Iz, const vector<int> &Ig, vector<double> &total_matrix, const Element &element)const;
		//将单刚添加到总刚(行列单元节点不同)
		void Add_to_gsmatrix(const double (*element_stiff_matrix)[24], const vector<long> &Iz, const vector<int> &Ig, vector<double> &total_matrix, const Element &ele_row, const Element &ele_col)const;
		//计算单元每个高斯点关于最近距离纳米管的权重值
		int Generate_element_gauss_weight_cnt(ifstream &infile, const double &dist, const struct RVE_Geo &cell, const struct CNT_Geo &cnts, const vector<Point_3D> &cnps, const vector<Element> &elements, 
																			   const vector<Node> &nodes, const vector<Node> &gauss, const vector<double> &weight, const double (*gauss_po)[3], const double (*ele_cent)[3], double (*gauss_wcnt)[10], double &Jacobi)const;
		//计算单元每个高斯点的数据
		void Generate_element_gauss_data(const vector<Element> &elements, const vector<Node> &nodes, const string &com_mod, const vector<Node> &gauss, const vector<double> &weight,
																	  double (*gauss_ns)[8], double (*gauss_nw)[8], double (*gauss_dfx)[3][8], double &Jacobi,  double (*gauss_po)[3], double (*ele_cent)[3])const;	
		//判断那些由于周期边界条件产生的相关单元的情况，并修改相关值
		void Verify_Periodical_Condition_Modify_Values(const double elec_left[], double elec_right[], double (*disp_right)[24], const struct RVE_Geo &cell, const struct Decay_Para &decay)const;
		//计算RVE背景网格单元(长程力作用)等效应变能
		void Calculate_Longforce_Energy(double elements_energy[], const double weight[], const double (*gauss_ns)[8], const struct Decay_Para &decay, const double (*gauss_po)[3], const double elec_left[],
																	const double elec_right[], const double (*gwc_left)[10], const double (*gwc_right)[10], const double (*disp_left)[24], const double (*disp_right)[24], const int &GS)const;
		//计算RVE背景网格单元(接触力作用)等效应变能
		void Calculate_Contactforce_Energy(double ele_energy[], const double (*disp_left)[24], const vector<MatPro> &mats, const double &Jacobi, const double (*gauss_dfx)[3][8], const double weight[], const int &GS)const;
		//读入信息一行，跳过注释行（以%开头）
		string Get_Line(ifstream &infile)const;

};
//---------------------------------------------------------------------------
#endif
//===========================================================================
