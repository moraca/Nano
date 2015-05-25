//===========================================================================
// Algorithm_FEM.h
// 算法及有限元实现模块类头文件
// A class to implement the special algorithm and finite element method
//===========================================================================
#ifndef AlGORITHM_H
#define AlGORITHM_H

#include<vector>
#include<iomanip>
#include<iostream>
#include<fstream>
#include<sstream>
#include "Fem_3D.h"
#include "MatBase.h"
#include "Mesher.h"
#include "Global_Stiff_Matrix.h"
#include "SolveEqu.h"

#include "limits.h"
#include "time.h"
#include "Hns.h"
using namespace hns;

const int NDIM=3;
const int IRN=1;					         // the number of element group
const int NGAP=1;

//---------------------------------------------------------------------------
//均匀化参数求解类
class Algorithm_FEM
{
	public:
		//数据变量
		SolveEqu *Solv;  //线性方程组求解类
		Global_Stiff_Matrix *Glosmat; //刚度矩阵组装类
		double equivalent_energy[9];  //单元等效能量
		vector<vector<double> > U_Solution;  //位移解向量
		string backup_file_name;  //备份位移解及等效应变能数据文件名
		string wr_mod;	//0 表示将写位移解数据，1表示将读位移解数据备份
		string com_mod;	//记录计算模式: Mod[0] 代表Local Continuum计算部分，Mod[1] 代表Nonlocal Continuum计算部分
		string homo_comp; //是否进行均匀化参数计算

		//构造函数
		Algorithm_FEM(){};
	
		//成员函数
		//求解均匀化参数;
		int Solve(ifstream &infile, vector<Node> &nodes, const vector<int>* peri_bnods, vector<Element> &elements,  const vector<MatPro> &mats, const struct Decay_Para &decay, const struct RVE_Geo &cell,
						const struct CNT_Geo &cnts, const vector<Point_3D> &cnps);

	private:
		//读入计算类型数据
		int Import_computational_data(ifstream &infile);
		//判断节点的相关节点和单元的相关单元
		int Deter_relative_nodes_elements(vector<Node> &nodes, vector<Element> &elements, const double &dist, const struct RVE_Geo &cell, string mod);
		//数据读入或者写入二进制文件，以备以后单独求解线性方程组时用到
		int Write_read_matrix_equright_data(const long &Total_Matrix_Size, vector<double> &total_matrix, vector<double>* nl_equright, double backup_ege[], string &maeq_mod);
		//读取位移解
		void read_u_solution(const string &output_file_name);
		//输出位移解及等效应变能
		void write_u_solution_equivalent_energy(const string &output_file_name)const;
		//读入信息一行，跳过注释行（以%开头）
		string Get_Line(ifstream &infile)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
