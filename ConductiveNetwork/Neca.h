//===========================================================================
//��Ҫ�������
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
		//���ݳ�Ա
		Preprocessor *Prep;					//ǰ����
		Algorithm_FEM *Algofem;		//�㷨������Ԫʵ��
		Postprocessor *Post; //����
        RNetwork *RNet;
    
		//���캯��
		Neca(){};

		//��Ա����
		int Begin(ifstream &infile, const int &samples_count);

	private:

		//��Ա����
		string Get_Line(ifstream &infile)const;				//������Ϣһ�У�����ע���У���%��ͷ��
        int ReadFromFile(vector<Point_3D> &point, vector<double> &radius, RVE_Geo &cell_geo, CNT_Geo &cnts_geo, Region_Geo &cnt_regions, vector<vector<long int> > &structure, vector<vector<int> > &sectioned_domain_cnt);
};
//------------------------------------------------
#endif
//===========================================================================
