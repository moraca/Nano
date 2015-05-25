//===========================================================================
// Hns.cpp
// 韩非的名空间类头文件
// Member functions in classes of Fei Han's namespace definition
//===========================================================================
#include "Hns.h"
//---------------------------------------------------------------------------
//名空间成员定义
foutstream hns::hout;
SetWP hns::setwp;
//---------------------------------------------------------------------------
void hns::open_deffo_stream(char* file,int mod)
{
	hout.open_deffo_stream(file,mod);
}
//---------------------------------------------------------------------------
void hns::close_deffo_stream( )
{
	hout.close_deffo_stream( );
}                                                                       
//---------------------------------------------------------------------------
foutstream& hns::endl(foutstream& os)
{
	return (os << "\n");
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//文件输出流(foutstream)类成员定义
foutstream::foutstream()
{          
	is_open_file = false;
	out_stream   = &cout;
}
//---------------------------------------------------------------------------
void foutstream::open_deffo_stream(char* file,int mod)
{
	if( is_open_file ) close_deffo_stream();
	of_stream.open(file);
	if( !of_stream )
	{
		cout << "无法打开文件！" << file << std::endl;
		return ;
	}
	out_stream = &of_stream;
	output_file = string(file);
	is_open_file = true;
}
//---------------------------------------------------------------------------
void foutstream::close_deffo_stream( )
{
	of_stream.close();
	out_stream   = &cout;
	is_open_file = false;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( short n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( unsigned short n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( int n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( unsigned int n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
//foutstream& foutstream:: operator <<( size_t n )
//{
//	(*out_stream) << n ;
//	out_stream->flush();
//	return *this;
//}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( long n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( unsigned long n )
{
	(*out_stream) << n ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( float f )
{
	(*out_stream) << f ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( double f )
{
	(*out_stream) << f ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( long double f )
{
	(*out_stream) << f ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( char c )
{
	(*out_stream) << c ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( bool b )
{
	(*out_stream) << b ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( string s )
{
	(*out_stream) << s ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<( foutstream& (*os)(foutstream& yt) )
{
	(*os)(*this);
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator <<(const char* c )
{
	(*out_stream) << c ;
	out_stream->flush();
	return *this;
}
//---------------------------------------------------------------------------
foutstream& foutstream:: operator << (SetWP& sw)
{
	this->width(sw.width);
	if( sw.precision != 0 )
	{
		this->precision(sw.precision);
	}
	return *this;
}
//---------------------------------------------------------------------------
int foutstream:: width ( int n )
{
	return int(out_stream->width(n));
}
//---------------------------------------------------------------------------
int foutstream:: precision( int n )
{
	return int(out_stream->precision(n));
}
//---------------------------------------------------------------------------
int foutstream:: setf( int n )
{
//	return out_stream->setf(n);
	return 1;
}
//---------------------------------------------------------------------------
int foutstream:: setf( int n1, int n2 )
{
//	return out_stream->setf(n1,n2);
	return 1;
}
//===========================================================================
