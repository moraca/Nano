//===========================================================================
// Mesher.h
// ��ά�����ʷ���ͷ�ļ�
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

		//���ݳ�Ա
		vector<Node> nodes;				//�ڵ�������
		vector<Element> elements;		//�����嵥Ԫ������
		vector<int> peri_bnods[2];		//���ڱ߽�ڵ���һά����([1]�м�¼)����ʼ�ڵ�λ��([2]�м�¼: 0Ϊ��ʼ�㣬1Ϊ�м��)

		//���캯����
		Mesher(){};
		//��������
		int Generate_mesh(ifstream &infile, struct RVE_Geo &cell_geo, vector<Point_3D> &cnps);
		//���ɸ�������Ԥ���ֲ�ģ�͵�Ч�ն���
		int Generate_grids_for_effective_stiffness(ifstream &infile, const double &decayR, double &grid_vol);
	
	private:
		double dx, dy, dz; //����������ʷ�ϸ��

		//��Ա����
		//������������������
		int Brick_background_mesh(ifstream &infile, struct RVE_Geo &cell_geo);
		//��������Ԫ��������
		int Element_material_property(ifstream &infile);
		//��¼��Ԫ��������׹��߶�
		int Element_relative_cnts(vector<Point_3D> &cnps, const struct RVE_Geo &cell_geo);
		//�����ⵥԪ�Ƿ�ƥ������ȷ��������׹��߶���Ϣ
		int Test_element_relative_cnts(vector<Point_3D> &cnps)const;
		//���Tecplot���ӻ���������
		int Export_mesh_data_tecplot(string output_file_name)const;
		//���ݸ�����i j k�Լ����i_max j_max k_max, ���������ڵ��λ�ã��ǵ㡢�߽��ߡ��߽��桢�ڲ���
		int Deter_node_type(const int i, const int j, const int k, const int i_max, const int j_max, const int k_max)const;
		//һά�洢���ڱ߽�ڵ��ż���ʼ����Ϣ
		int Record_periodic_boundary_nodes_number(const struct RVE_Geo &cell_geo);
		//������Ϣһ�У�����ע���У���%��ͷ��
		string Get_Line(ifstream &infile)const;
};
//-------------------------------------------------------
#endif
//=====================================================
