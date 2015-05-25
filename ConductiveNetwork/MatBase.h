//===========================================================================
// MatBase.cpp
// ���Ͽ���ͷ�ļ�
// A class of material database
//===========================================================================

#ifndef MATBASE_H
#define MATBASE_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<vector>
#include<string>
#include "MatPro.h"
#include "Gauss.h"
#include "Mesher.h"
#include "Hns.h"
using namespace hns;
//---------------------------------------------------------------------------
class MatBase
{
	public:

		//���ݳ�Ա
		vector<MatPro> mats_vec;
		struct Decay_Para decay;

		//���캯��
		MatBase(){};

		//��Ա����
		//���ɲ������ݿ�
		int Generate_matbase(ifstream &infile);

	private:
		//��Ա����
		//������в��ϵĸնȾ���
		void Print_stiffness_matrix(string print_name)const;
		void Print_stiffness_matrix(string print_name, MatPro &mat, int key)const;
		//Ԥ���ֲ�ģ�͵�Ч�ն���
		int Estimate_stiffness_long_range_interaction(const vector<Element> &elements, const vector<Node> &nodes, const double &grid_vol, 
								const vector<Node> &gauss, const vector<double> &weight, const double a[], MatPro &mat, const int &output_key)const;
		//������Ϣһ�У�����ע���У���%��ͷ��
		string Get_Line(ifstream &infile)const;
};
//-------------------------------------------------------
#endif
//===========================================================================
