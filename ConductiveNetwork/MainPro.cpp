//===========================================================================
//Non-local Elasticity models for Composite Analysis (NECA)
//在整个工程的开始加这么一段程序的目的是为了以后做动态问题时
//需要多时间步循环运算做准备
//===========================================================================
#include "Neca.h"
#include <string>
#include <sys/time.h>
int main(int argc, char** argv)
{
	long double timer;
	struct  timeval  start, end;
	gettimeofday(&start, NULL);

	//将输入文件名记录在infile中
	string in_file;
	if(argc > 1)
	{
		in_file = argv[1];
	}
	else
	{
		cout << "请输入用户输入参数文本文件名： " << endl;
		in_file = "input.dat";
		cout << in_file << endl;
	};

	//将输入文件名记录在outfile中
	string out_file;
	if(argc > 2)
	{
		out_file = argv[2];
	}
	else
	{
		cout << "请输入程序运行记录及输出文本文件名： " << endl;
		out_file = "output.dat";
		cout << out_file << endl;
	};

	//开始计算

	//确定输出文件
	if( out_file.size() > 0 ) open_deffo_stream( (char*)out_file.c_str() );
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//定义计算样本
	Neca *Compute =  new  Neca;

	int samples_count = 0;	//执行样本计数器
	int samples_num=1;		//设定样本数，最少为1
	do									//执行每个样本
	{
		//-----------------------------------------------------------------------------------------------
		//计算样本数累计
		samples_count++;

		//-----------------------------------------------------------------------------------------------
		//打开用户输入参数文本文件
		ifstream infile;
		infile.open(in_file.c_str());
		if(!infile) { hout << "不能打开输入文件: "  << in_file << endl;  return 0; }

		//---------------------------------------------------------------------------------------------
		string s;		//读入信息一行
		getline(infile,s);	//跳过注释行  
		while(!infile.eof() && s.substr(0,1)=="%")	getline(infile,s);
		istringstream istr(s);		//读取计算样本数信息
		istr >> samples_num;
		if(samples_num<1) { hout << "样本数输入错误！请检查！ " << endl;  return 0; }

		//----------------------------------------------------------------------------------------------
		//计算样本
		if(Compute->Begin(infile, samples_count)==0) return 0;

		//----------------------------------------------------------------------------------------------
	}while(samples_count<samples_num);
	
	//删除计算样本
	delete Compute;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//输出真实时间
	gettimeofday( &end, NULL);
	timer = end.tv_sec-start.tv_sec+(end.tv_usec-start.tv_usec)/1000000.0; 
	hout << "真实计算总共耗时：" << timer << "秒" << endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//关闭输出文件
	close_deffo_stream();

	return 1;
}
//===========================================================================
