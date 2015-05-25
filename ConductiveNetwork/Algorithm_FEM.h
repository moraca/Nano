//===========================================================================
// Algorithm_FEM.h
// �㷨������Ԫʵ��ģ����ͷ�ļ�
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
//���Ȼ����������
class Algorithm_FEM
{
	public:
		//���ݱ���
		SolveEqu *Solv;  //���Է����������
		Global_Stiff_Matrix *Glosmat; //�նȾ�����װ��
		double equivalent_energy[9];  //��Ԫ��Ч����
		vector<vector<double> > U_Solution;  //λ�ƽ�����
		string backup_file_name;  //����λ�ƽ⼰��ЧӦ���������ļ���
		string wr_mod;	//0 ��ʾ��дλ�ƽ����ݣ�1��ʾ����λ�ƽ����ݱ���
		string com_mod;	//��¼����ģʽ: Mod[0] ����Local Continuum���㲿�֣�Mod[1] ����Nonlocal Continuum���㲿��
		string homo_comp; //�Ƿ���о��Ȼ���������

		//���캯��
		Algorithm_FEM(){};
	
		//��Ա����
		//�����Ȼ�����;
		int Solve(ifstream &infile, vector<Node> &nodes, const vector<int>* peri_bnods, vector<Element> &elements,  const vector<MatPro> &mats, const struct Decay_Para &decay, const struct RVE_Geo &cell,
						const struct CNT_Geo &cnts, const vector<Point_3D> &cnps);

	private:
		//���������������
		int Import_computational_data(ifstream &infile);
		//�жϽڵ����ؽڵ�͵�Ԫ����ص�Ԫ
		int Deter_relative_nodes_elements(vector<Node> &nodes, vector<Element> &elements, const double &dist, const struct RVE_Geo &cell, string mod);
		//���ݶ������д��������ļ����Ա��Ժ󵥶�������Է�����ʱ�õ�
		int Write_read_matrix_equright_data(const long &Total_Matrix_Size, vector<double> &total_matrix, vector<double>* nl_equright, double backup_ege[], string &maeq_mod);
		//��ȡλ�ƽ�
		void read_u_solution(const string &output_file_name);
		//���λ�ƽ⼰��ЧӦ����
		void write_u_solution_equivalent_energy(const string &output_file_name)const;
		//������Ϣһ�У�����ע���У���%��ͷ��
		string Get_Line(ifstream &infile)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
