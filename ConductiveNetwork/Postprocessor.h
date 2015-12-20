//========================================================
// Postprocessor.h
// �������ͷ�ļ�
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
		//���캯��
		Postprocessor(){};  

		//��Ա����
		int Treatment(ifstream &infile, vector<Node> &nodes, const vector<int>* peri_bnods, const vector<Element> &elements, const vector<MatPro> &mats, 
								const vector<vector<double> > &u_solution, const string &wr_mod, const string &bakcup_file_name, double equivalent_energy[], const int &samples_count); //ִ�к������

	private:

		//���ݱ���
		vector<vector<double> > u;
		vector<vector<double> > ele_strain;
		vector<vector<double> > nod_strain;

		//��Ա����
		//������ȡλ�ƽ�
		void read_u_solution_equivalent_energy(const string &output_file_name, double equivalent_energy[]);
		//���㵥Ԫ���Ĵ�Ӧ������
		int Calculate_Elements_Strain(const vector<Node> &nodes, const vector<Element> &elements);
		//����ڵ����Ĵ�Ӧ������
		int Calculate_Nodes_Strain(vector<Node> &nodes, const vector<int>* &peri_bnods, const vector<Element> &elements);
		//���Tecplot���ӻ�λ�Ƴ���ͼ
		int Export_Displacement_Contour(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		//����������ͼ
		int Export_Deformed_Mesh(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		//����������ͼ(���ֲ���)
		int Export_Deformed_Mesh_Multimats(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements, const vector<MatPro> &mats)const;
		//���Ӧ�䳡��ͼ
		int Export_Strain_Mesh(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		//���Ӧ�䳡��ͼ(��׼����, ����û�б���)
		int Export_Strain_Contour(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const;
		//���㲢������Ȼ�ϵ��
		int Calculate_Export_Homogenized_Coefficients(const string &output_file_name, const vector<Element> &elements, const double equivalent_energy[], const int &samples_count)const;
		//��ȡ�����ļ�һ�У�����ע���У���%��ʼ��
        char* Get_Line(ifstream& input_stream, char* read_line);
  		//������Ϣһ�У�����ע���У���%��ͷ��
		string Get_Line(ifstream &infile)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
