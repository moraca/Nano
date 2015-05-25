//========================================================
// Postprocessor.h
// 后处理过程头文件
// A class to implement preprocessor
//========================================================

#ifndef POSTPROCESSOR_H
#define POSTPROCESSOR_H


#include<iomanip>
#include<iostream>
#include<fstream>
#include<sstream>
#include"Fem_3D.h"
#include "MatBase.h"
#include "Geometry_3D.h"
#include "Global_Stiff_Matrix.h"
#include "time.h"
#include "Hns.h"
using namespace hns;

//---------------------------------------------------------------------------
class Postprocessor 
{
	public:
		//构造函数
		Postprocessor(){};  

		//成员函数
		int Treatment(ifstream &infile, vector<Node> &nodes, const vector<int>* peri_bnods, const vector<Element> &elements, const vector<MatPro> &mats, 
								const vector<vector<double> > &u_solution, const string &wr_mod, const string &bakcup_file_name, double equivalent_energy[], const int &samples_count); //执行后处理过程

	private:

		//数据变量
		vector<vector<double> > u;
		vector<vector<double> > ele_strain;
		vector<vector<double> > nod_strain;

		//成员函数
		//输出或读取位移解
		void read_u_solution_equivalent_energy(const string &output_file_name, double equivalent_energy[]);
		//计算单元形心处应变向量
		int Calculate_Elements_Strain(const vector<Node> &nodes, const vector<Element> &elements);
		//计算节点形心处应变向量
		int Calculate_Nodes_Strain(vector<Node> &nodes, const vector<int>* &peri_bnods, const vector<Element> &elements);
		//输出Tecplot可视化位移场云图
		int Export_Displacement_Contour(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		//输出网格变形图
		int Export_Deformed_Mesh(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		//输出网格变形图(多种材料)
		int Export_Deformed_Mesh_Multimats(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements, const vector<MatPro> &mats)const;
		//输出应变场云图
		int Export_Strain_Mesh(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		//输出应变场云图(标准网格, 本身没有变形)
		int Export_Strain_Contour(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		//计算并输出均匀化系数
		int Calculate_Export_Homogenized_Coefficients(const string &output_file_name, const vector<Element> &elements, const double equivalent_energy[], const int &samples_count)const;
		//读取输入文件一行，跳过注释行（以%开始）
        char* Get_Line(ifstream& input_stream, char* read_line);
  		//读入信息一行，跳过注释行（以%开头）
		string Get_Line(ifstream &infile)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
