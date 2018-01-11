//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Hns.h
//OBJECTIVE:	Classes for Han's namespace definition
//AUTHOR:		Fei Han;
//E-MAIL:			fei.han@kaust.edu.sa	;
//====================================================================================

#ifndef HNS_H
#define HNS_H

#include<fstream>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<string>

using namespace std;

//---------------------------------------------------------------------------
class SetWP
{
public:
	SetWP& operator()(int w)
	{
		width = w;
		precision = 0;
		return *this;
	}
	SetWP& operator()(int w,int p)
	{
		width = w;
		precision = p;
		return *this;
	}
	int width;
	int precision;
};
//---------------------------------------------------------------------------
class foutstream 
{
public:
	//Constructor
	foutstream();

	//the definitions for operators
	foutstream& operator << (short n);
	foutstream& operator << (unsigned short n);
	foutstream& operator << (int n);
	foutstream& operator << (unsigned int n);
//	foutstream& operator << (size_t n);
	foutstream& operator << (long n);
	foutstream& operator << (unsigned long n);
	foutstream& operator << (float f);
	foutstream& operator << (double f);
	foutstream& operator << (long double f);
	foutstream& operator << (char c);
	foutstream& operator << (bool b);
	foutstream& operator << (string s);
	foutstream& operator << (foutstream& (*os)(foutstream& ys));
	foutstream& operator << (SetWP& sw);
	foutstream& operator << (const char* c);

	int width(int n);
	int precision(int n);
	int setf(int n);
	int setf(int n1,int n2); 

//	void open_deffo_stream(char* file,int mod=ios::app);		//app: continue to write
	void open_deffo_stream(char* file,int mod=ios::out);			//out: rewrite
	void close_deffo_stream( );

	string output_file;
private:
	ostream *out_stream; 
	ofstream of_stream;
	bool is_open_file;
};
//---------------------------------------------------------------------------
namespace hns
{
	extern foutstream hout;
//	void open_deffo_stream(char* file,int mod=ios::app);	//app: continue to write
	void open_deffo_stream(char* file,int mod=ios::out);		//out: rewrite

	void close_deffo_stream( );
	foutstream& endl(foutstream& os);
	extern SetWP setwp;
}

//---------------------------------------------------------------------------
#endif
//===========================================================================
