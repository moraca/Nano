//===========================================================================
// Preprocessor.h
// 前处理过程头文件：生成材料信息库，计算区域的几何建模并进行网格剖分
// A class to implement preprocessor
//===========================================================================

#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

#include<iomanip>
#include<iostream>
#include<fstream>
#include<sstream>
#include "Fem_3D.h"
#include "MatBase.h"
#include "GeoNano.h"
#include "Mesher.h"
#include "time.h"
#include "Hns.h"

using namespace hns;

//---------------------------------------------------------------------------
class Preprocessor 
{
	public:
		//数据成员
		MatBase *Matrial;			//前处理之材料信息
		GeoNano *Geonano;		//前处理之几何信息
		Mesher *Mesh;				//前处理之网格信息

		//构造函数
		Preprocessor(){};  

		//成员函数
		int Implement(ifstream &infile, const int &samples_count); //执行前处理过程

private:

		//成员函数
		//读取输入文件一行，跳过注释行（以%开始）
        char* Get_Line(ifstream& input_stream, char* read_line);
  		//读入信息一行，跳过注释行（以%开头）
		string Get_Line(ifstream &infile)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
