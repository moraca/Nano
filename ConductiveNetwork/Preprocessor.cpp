//===========================================================================
// Preprocessor.cpp
// 前处理过程：生成材料信息库，计算区域的几何建模并进行网格剖分
// Member functions in a class to implement preprocessor
//===========================================================================
#include "Preprocessor.h"

//--------------------------------------------------------------------------
//执行前处理过程
int Preprocessor::Implement(ifstream &infile, const int &samples_count)
{
	clock_t ct0,ct1;  //定义变量记录执行的开始和结束时间，分别计算并输出生成材料信息、几何建模和网格剖分的耗时
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//生成材料库信息	
	ct0 = clock();
	hout << "-_- 开始读取材料信息......"<<endl;
	Matrial = new MatBase;
	if(Matrial->Generate_matbase(infile)==0) return 0; //生成材料数据库
	ct1 = clock();
	hout << "    读取材料信息总耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
	hout << "^_^ 读取材料信息完毕！"<<endl<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//几何建模
	ct0 = clock();
	hout << "-_- 开始几何建模......"<<endl;
	hout << "GeoNano in Pre-processor begin......"<<endl;
	Geonano = new GeoNano;
	if(Geonano->Geometric_modeling(infile, samples_count)==0) return 0;       //几何建模
	ct1 = clock();
	hout << "GeoNano time: "<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<" secs"<<endl;
	//hout << "    建模总耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
	hout << "^_^ 几何建模完毕！"<<endl<<endl;
	hout << "GeoNano in Pre-processor end......"<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	/*/生成网格
	ct0 = clock();
	hout << "-_- 开始生成网格 (Mesher start)......"<<endl;
	Mesh = new Mesher;
	if(Mesh->Generate_mesh(infile, Geonano->cell_geo, Geonano->cnps)== 0) return 0;  //生成六面体背景网格
	ct1 = clock();
	hout << "    生成网格总耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
	hout << "^_^ 网格生成完毕！ ....... (Mesher end) ......"<<endl<<endl; //*/

	return 1;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//读取输入文件一行，跳过注释行（以%开始）
char* Preprocessor::Get_Line(ifstream& input_stream, char* read_line)
{
	//读入的一行char read_line[200]     
	input_stream.getline(read_line,200);
	//跳过注释行
	while(!input_stream.eof() && read_line[0] == '%')
	{
		input_stream.getline(read_line,200);
	}

	return read_line;
}
//读入一行信息，并跳过注释行（以"%"开头）；
string Preprocessor::Get_Line(ifstream &infile)const
{
	string s;
	//读入信息一行
	getline(infile,s);
	//跳过注释行     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===========================================================================
