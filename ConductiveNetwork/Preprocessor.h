//===========================================================================
// Preprocessor.h
// ǰ�������ͷ�ļ������ɲ�����Ϣ�⣬��������ļ��ν�ģ�����������ʷ�
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
		//���ݳ�Ա
		MatBase *Matrial;			//ǰ����֮������Ϣ
		GeoNano *Geonano;		//ǰ����֮������Ϣ
		Mesher *Mesh;				//ǰ����֮������Ϣ

		//���캯��
		Preprocessor(){};  

		//��Ա����
		int Implement(ifstream &infile, const int &samples_count); //ִ��ǰ�������

private:

		//��Ա����
		//��ȡ�����ļ�һ�У�����ע���У���%��ʼ��
        char* Get_Line(ifstream& input_stream, char* read_line);
  		//������Ϣһ�У�����ע���У���%��ͷ��
		string Get_Line(ifstream &infile)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
