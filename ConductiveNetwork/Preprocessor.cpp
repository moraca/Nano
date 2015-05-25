//===========================================================================
// Preprocessor.cpp
// ǰ������̣����ɲ�����Ϣ�⣬��������ļ��ν�ģ�����������ʷ�
// Member functions in a class to implement preprocessor
//===========================================================================
#include "Preprocessor.h"

//--------------------------------------------------------------------------
//ִ��ǰ�������
int Preprocessor::Implement(ifstream &infile, const int &samples_count)
{
	clock_t ct0,ct1;  //���������¼ִ�еĿ�ʼ�ͽ���ʱ�䣬�ֱ���㲢������ɲ�����Ϣ�����ν�ģ�������ʷֵĺ�ʱ
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���ɲ��Ͽ���Ϣ	
	ct0 = clock();
	hout << "-_- ��ʼ��ȡ������Ϣ......"<<endl;
	Matrial = new MatBase;
	if(Matrial->Generate_matbase(infile)==0) return 0; //���ɲ������ݿ�
	ct1 = clock();
	hout << "    ��ȡ������Ϣ�ܺ�ʱ"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"�롣"<<endl;
	hout << "^_^ ��ȡ������Ϣ��ϣ�"<<endl<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���ν�ģ
	ct0 = clock();
	hout << "-_- ��ʼ���ν�ģ......"<<endl;
	hout << "GeoNano in Pre-processor begin......"<<endl;
	Geonano = new GeoNano;
	if(Geonano->Geometric_modeling(infile, samples_count)==0) return 0;       //���ν�ģ
	ct1 = clock();
	hout << "GeoNano time: "<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<" secs"<<endl;
	//hout << "    ��ģ�ܺ�ʱ"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"�롣"<<endl;
	hout << "^_^ ���ν�ģ��ϣ�"<<endl<<endl;
	hout << "GeoNano in Pre-processor end......"<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	/*/��������
	ct0 = clock();
	hout << "-_- ��ʼ�������� (Mesher start)......"<<endl;
	Mesh = new Mesher;
	if(Mesh->Generate_mesh(infile, Geonano->cell_geo, Geonano->cnps)== 0) return 0;  //���������屳������
	ct1 = clock();
	hout << "    ���������ܺ�ʱ"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"�롣"<<endl;
	hout << "^_^ ����������ϣ� ....... (Mesher end) ......"<<endl<<endl; //*/

	return 1;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//��ȡ�����ļ�һ�У�����ע���У���%��ʼ��
char* Preprocessor::Get_Line(ifstream& input_stream, char* read_line)
{
	//�����һ��char read_line[200]     
	input_stream.getline(read_line,200);
	//����ע����
	while(!input_stream.eof() && read_line[0] == '%')
	{
		input_stream.getline(read_line,200);
	}

	return read_line;
}
//����һ����Ϣ��������ע���У���"%"��ͷ����
string Preprocessor::Get_Line(ifstream &infile)const
{
	string s;
	//������Ϣһ��
	getline(infile,s);
	//����ע����     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===========================================================================
