//===========================================================================
// Mesher.h
// 三维网格剖分类头文件
// A class of 3D mesh generation
//===========================================================================

#ifndef MESHER_H
#define MESHER_H

#include<iomanip>
#include<iostream>
#include<fstream>
#include<sstream>
#include "Fem_3D.h"
#include "GeoNano.h"
#include "time.h"
#include "Hns.h"
using namespace hns;

//---------------------------------------------------------------------------
class Mesher 
{
	public:

		//数据成员
		vector<Node> nodes;				//节点类向量
		vector<Element> elements;		//六面体单元类向量
		vector<int> peri_bnods[2];		//周期边界节点编号一维向量([1]中记录)及开始节点位置([2]中记录: 0为开始点，1为中间点)

		//构造函数；
		Mesher(){};
		//网格生成
		int Generate_mesh(ifstream &infile, struct RVE_Geo &cell_geo, vector<Point_3D> &cnps);
		//生成格子用于预估局部模型等效刚度阵
		int Generate_grids_for_effective_stiffness(ifstream &infile, const double &decayR, double &grid_vol);
	
	private:
		double dx, dy, dz; //各个方向的剖分细度

		//成员函数
		//背景六面体网格生成
		int Brick_background_mesh(ifstream &infile, struct RVE_Geo &cell_geo);
		//赋予网格单元材料属性
		int Element_material_property(ifstream &infile);
		//记录单元的相关纳米管线段
		int Element_relative_cnts(vector<Point_3D> &cnps, const struct RVE_Geo &cell_geo);
		//输出检测单元是否匹配了正确的相关纳米管线段信息
		int Test_element_relative_cnts(vector<Point_3D> &cnps)const;
		//输出Tecplot可视化网格数据
		int Export_mesh_data_tecplot(string output_file_name)const;
		//根据给定的i j k以及最大i_max j_max k_max, 决定给定节点的位置（角点、边界线、边界面、内部）
		int Deter_node_type(const int i, const int j, const int k, const int i_max, const int j_max, const int k_max)const;
		//一维存储周期边界节点编号及起始点信息
		int Record_periodic_boundary_nodes_number(const struct RVE_Geo &cell_geo);
		//读入信息一行，跳过注释行（以%开头）
		string Get_Line(ifstream &infile)const;
};
//-------------------------------------------------------
#endif
//=====================================================
