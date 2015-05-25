//===========================================================================
// Fem_3D.h
// 三维有限单元类头文件
// A class of FEM
//===========================================================================

#ifndef FEM_3D_H
#define FEM_3D_H

#include<iostream>
#include<cmath>
#include<vector>
using namespace std;

//---------------------------------------------------------------------------
//定义节点类
class Node
{	
	public:
		double x, y, z; //节点的坐标值

		int type;	//0内点，1边界面点，2边界线点，3角点
		vector<int> relative_nods;	//相关节点向量（节点编号）
		vector<int> relative_eles;	//相关单元向量（单元编号）

		//构造函数
		Node(const double ix=0, const double iy=0, const double iz=0);
		
		//成员函数；
		double	distance_to(const Node& n); //距离
};
//---------------------------------------------------------------------------
//单元类头文件；
class Element
{
	public:	
		int type;	//表示单元的形状类型(三个数字xyz： x 维度；y单元的节点个数；z单元的形函数幂次)；
						//例如，121: 一维两节点(线段)线性形函数；231: 二维三节点(三角形)线性形函数；241: 二维四节点(四边形)线性形函数；
						//341: 三维四节点(四面体)线性形函数；361: 三维六节点(三棱柱)线性形函数；381：三维八节点(六面体)线性形函数
		int mat;	//表示单元的材料物性；
		vector<int> nodes_id;
		vector<int> relative_eles;	//相关单元向量（单元编号）
		vector<int> relative_cnts;	//相关纳米管线段 （纳米管上节点编号）
		
		//构造函数；
		Element(){};
		
		//成员函数；
};
//---------------------------------------------------------------------------
//定义六面体单元
class Hexahedron
{
	public:
        int nodesId[8];	 //eight nodes;
        int materialId;
        int facesId[6];

		//构造函数；
		Hexahedron(){};
        Hexahedron(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
		
		//成员函数；
};

#endif
//===========================================================================
