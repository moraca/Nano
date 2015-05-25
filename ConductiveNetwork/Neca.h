//===========================================================================
//主要计算过程
//===========================================================================

#ifndef NECA_H
#define NECA_H

#include "Hns.h"
#include "Preprocessor.h"
#include "Algorithm_FEM.h"
#include "Postprocessor.h"
#include "RNetwork.h"

#include "limits.h"
#include "time.h"
using namespace hns;

//---------------------------------------------------------------------------
class Neca
{
	public: 
		//数据成员
		Preprocessor *Prep;					//前处理
		Algorithm_FEM *Algofem;		//算法及有限元实现
		Postprocessor *Post; //后处理
        RNetwork *RNet;
    
		//构造函数
		Neca(){};

		//成员函数
		int Begin(ifstream &infile, const int &samples_count);

	private:

		//成员函数
		string Get_Line(ifstream &infile)const;				//读入信息一行，跳过注释行（以%开头）
        int ReadFromFile(vector<Point_3D> &point, vector<double> &radius, RVE_Geo &cell_geo, CNT_Geo &cnts_geo, Region_Geo &cnt_regions, vector<vector<long int> > &structure, vector<vector<int> > &sectioned_domain_cnt);
};
//------------------------------------------------
#endif
//===========================================================================
