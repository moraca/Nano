//===========================================================================
// MatPro.h
// ����������ͷ�ļ�
// A class of material property
//===========================================================================

#ifndef MATPRO_H
#define MATPRO_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<cmath>
#include "Gauss.h"
#include "MathMatrix.h"
#include "Hns.h"

using namespace hns;

//---------------------------------------------------------------------------
class MatPro
{
	public:	
	
		//���ݳ�Ա
		int type_val;   //��ά�����������ͣ�Ŀǰ��������(0������ͬ��; 1����۸���ͬ��; 2��������������)
		double E11, E22, E33, Nu12, Nu23, Nu13, G12, G23, G13;	//���ϵĵ���ģ�������ɱȣ�����ģ��
		double elas_matrix[6][6];		//Dij ���Ծ���

		//���캯����
		MatPro(int itype=0){ type_val = itype; }
		//---------------------------------------------
		//��Ա����
		//���ø���ͬ�Բ��ϲ����ĵ���ģ�������ɱȺͼ���ģ��
		void set_ela_para(const double E, const double Nu);
		//���ú�۸���ͬ�Բ��ϵĵ���ģ�������ɱȺͼ���ģ��
		void set_ela_para(const double E1, const double E2, const double Nu1, const double Nu2, const double G2);
		//���������������Բ��ϵĵ���ģ�������ɱȺͼ���ģ��
		void set_ela_para(const double iE11, const double iE22, const double iE33, const double iNu12, const double iNu23, const double iNu13, const double iG12, const double iG23, const double iG13);

		//���ø���ͬ�Բ��ϲ����ĵ���ģ�������ɱȺͼ���ģ��(���ڳ�������Ч�ն�)
		void set_ela_para(const double &E);
		//���ú�۸���ͬ�Բ��ϵĵ���ģ�������ɱȺͼ���ģ��(���ڳ�������Ч�ն�)(���벴�ɱ�)
		void set_ela_para(const double &E1, const double &E2, const double &Nu2);
		//���ú�۸���ͬ�Բ��ϵĵ���ģ�������ɱȺͼ���ģ��(���ڳ�������Ч�ն�)(�������ģ��)
		void set_ela_para_transverse(const double &E1, const double &E2, const double &G2);
		//���������������Բ��ϵĵ���ģ�������ɱȺͼ���ģ��(���ڳ�������Ч�ն�)
		void set_ela_para(const double &iE1, const double &iE2, const double &iE3, const double &iNu12, const double &iNu23, const double &iNu13);
		//(�ֲ�ģ��)����Dij���Ծ��󣬴˵��Ծ����Ӧ��Ӧ������Ϊ��e11,e22,e33,2*e12,2*e23,2*e31��ת��
		int Generate_local_elas_matrix();
		//(�Ǿֲ�ģ��)����Dij���Ծ��󣬴˵��Ծ����Ӧ��Ӧ������Ϊ��e11,e22,e33,2*e12,2*e23,2*e31��ת��
		int Generate_nonlocal_elas_matrix();
		//����Dij���Ծ��󣬷�����ϲ���
		int Get_ele_para_by_ela_matrix();
		//�������۹�ʽ����ϵ��������֤һ���ԣ�����Ҫ����Horizon�뾶������Horizon�뾶Ӧ�ô���1������������
		void Compare_coef_by_analysis_formula(const int mat_type, const struct Decay_Para &decay)const;
		//�������۹�ʽ����նȾ��󣬲���֤һ����
		void Compare_matrix_by_analysis_formula(const struct Decay_Para &decay)const;

		void print()const;
};
//---------------------------------------------------------------------------
//��¼˥��������Ϣ
struct Decay_Para
{	
	double R;			//˥��������뾶(Horizon)
	double radius;		//˥���뾶(��ʽ����)
	double acoe[6];	//��������ϵ��
};
//------------------------------------------------
#endif
//===========================================================================
