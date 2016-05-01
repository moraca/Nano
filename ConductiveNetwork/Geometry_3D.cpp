//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Geometry_3D.cpp
//OBJECTIVE:	A class for basic geometic elements in 3D such as point, line and plane 
//AUTHOR:		Fei Han;
//E-MAIL:			fei.han@kaust.edu.sa	;
//====================================================================================

#include "Geometry_3D.h"

//---------------------------------------------------------------------------
//The member function for the 3D point class
Point_3D::Point_3D( double px, double py, double pz )
{
	x = px;
	y = py;
	z = pz;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator+( Point_3D &pt )
{
	Point_3D rp( x + pt.x, y + pt.y, z + pt.z );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator+( const Point_3D &pt )
{
	Point_3D rp( x + pt.x, y + pt.y, z + pt.z );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator+( double d )
{
	Point_3D rp( x + d, y + d, z + d );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator-( Point_3D &pt )
{
	Point_3D rp( x - pt.x, y - pt.y, z - pt.z );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator-( const Point_3D &pt )
{
	Point_3D rp( x - pt.x, y - pt.y, z - pt.z );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator-( double d )
{
	Point_3D rp( x - d, y - d, z - d );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator*( double d )
{
	Point_3D rp( x*d, y*d, z*d );
	return rp;
}
//---------------------------------------------------------------------------
Point_3D Point_3D::operator/( double d )
{
	Point_3D rp( x/d, y/d, z/d );
	return rp;
}
//---------------------------------------------------------------------------
bool Point_3D::operator==( Point_3D &pt )
{
	return (x==pt.x&&y==pt.y&&z==pt.z);
}
//---------------------------------------------------------------------------
bool Point_3D::operator!=( Point_3D &pt )
{
	return (x!=pt.x||y!=pt.y||z!=pt.z);
}
//---------------------------------------------------------------------------
double Point_3D::distance_to(const Point_3D &pt )const
{
	double rv2 = (x-pt.x)*(x-pt.x)+(y-pt.y)*(y-pt.y)+(z-pt.z)*(z-pt.z);
	return sqrt(rv2);
}
//---------------------------------------------------------------------------
double Point_3D::distance_to(const double &px, const double &py, const double &pz)const
{
	double rv2 = (x-px)*(x-px)+(y-py)*(y-py)+(z-pz)*(z-pz);
	return sqrt(rv2);
}
//---------------------------------------------------------------------------
Point_3D Point_3D::cross(Point_3D &point)
{
    double a = y*point.z - point.y*z;
    double b = -(x*point.z - point.x*z);
    double c = x*point.y - point.x*y;
    
    return Point_3D(a,b,c);
}
//---------------------------------------------------------------------------
double Point_3D::dot(Point_3D &point)
{
    return x*point.x + y*point.y + z*point.z;
}

//===========================================================================
//The member function for the 3D line class
//---------------------------------------------------------------------------
//Constructor
Line_3D::Line_3D(Point_3D p0, Point_3D p1)
{
	point[0] = p0;
	point[1] = p1;
	xm = p1.x - p0.x;
	yn = p1.y - p0.y;
	zl = p1.z - p0.z;
	len = length();
	if(len==0) virtual_line = false;
	else virtual_line = true;
}
//---------------------------------------------------------------------------   
//The length of segment
double Line_3D::length() 
{
	double dx = point[1].x-point[0].x;
	double dy = point[1].y-point[0].y;
	double dz = point[1].z-point[0].z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}
//---------------------------------------------------------------------------
//The distance from a point to a line
double Line_3D::distance_point_to_line(const Point_3D *point_temp)const
{
	double dis = 0;
	if(xm==0&&yn==0&&zl==0)
	{
		hout << "Attention: this segment is reduced to a point!" <<endl;
		double X = point_temp->x-point[0].x;
		double Y = point_temp->y-point[0].y;
		double Z = point_temp->z-point[0].z;
		dis = sqrt(X*X+Y*Y+Z*Z);
	}
	else
	{
		double X = yn*(point_temp->y-point[0].y)-zl*(point_temp->z-point[0].z);
		double Y = zl*(point_temp->z-point[0].z)-xm*(point_temp->x-point[0].x);
		double Z = xm*(point_temp->x-point[0].x)-yn*(point_temp->y-point[0].y);
		dis = sqrt(X*X+Y*Y+Z*Z)/sqrt(xm*xm+yn*yn+zl*zl);
	}
	return dis;
}
double Line_3D::distance_point_to_line(const Point_3D &point_temp)const
{
	double dis = 0;
	if(xm==0&&yn==0&&zl==0)
	{
		hout << "Attention: this segment is reduced to a point!" <<endl;
		double X = point_temp.x-point[0].x;
		double Y = point_temp.y-point[0].y;
		double Z = point_temp.z-point[0].z;
		dis = sqrt(X*X+Y*Y+Z*Z);
	}
	{
		double X = yn*(point_temp.y-point[0].y)-zl*(point_temp.z-point[0].z);
		double Y = zl*(point_temp.z-point[0].z)-xm*(point_temp.x-point[0].x);
		double Z = xm*(point_temp.x-point[0].x)-yn*(point_temp.y-point[0].y);
		dis = sqrt(X*X+Y*Y+Z*Z)/sqrt(xm*xm+yn*yn+zl*zl);
	}
	return dis;
}
double Line_3D::distance_point_to_line(const double dx, const double dy, const double dz)const
{
	double dis = 0;
	if(xm==0&&yn==0&&zl==0)
	{
		hout << "Attention: this segment is reduced to a point!" <<endl;
		double X = dx-point[0].x;
		double Y = dy-point[0].y;
		double Z = dz-point[0].z;
		dis = sqrt(X*X+Y*Y+Z*Z);
	}
	{
		double X = yn*(dy-point[0].y)-zl*(dz-point[0].z);
		double Y = zl*(dz-point[0].z)-xm*(dx-point[0].x);
		double Z = xm*(dx-point[0].x)-yn*(dy-point[0].y);
		dis = sqrt(X*X+Y*Y+Z*Z)/sqrt(xm*xm+yn*yn+zl*zl);
	}
	return dis;
}
//---------------------------------------------------------------------------
//To judge if a point is in a segment
int Line_3D::contain(const Point_3D &point_temp)const
{
	//to judge if the distance from a point to two endpoints is larger than the distance between endpoints
	if( fabs(point_temp.distance_to(point[0])+point_temp.distance_to(point[1])-point[0].distance_to(point[1]))>Zero ) return 0; 
	return 1;
}
//===========================================================================

//The member function for the 3D plane class
//---------------------------------------------------------------------------
//Constructor
Plane_3D::Plane_3D(double para[])
{
	for(int i=0; i<4; i++)	coef[i] = para[i];
	if(coef[0]==0.0&&coef[1]==0.0&&coef[2]==0.0) virtual_plane = false;
	else virtual_plane = true;
}
//---------------------------------------------------------------------------
//To judge if a point is contained in this plane
int Plane_3D::contain(const Point_3D &point_temp)const
{
	if( coef[0]*point_temp.x+coef[1]*point_temp.y+coef[2]*point_temp.z+coef[3]==0 ) return 1;	//in the plane
	return 0;	//out of the plane
}
//---------------------------------------------------------------------------
//To judge if a point is contained in this plane
int Plane_3D::contain(const double dx, const double dy, const double dz)const
{
	if( coef[0]*dx+coef[1]*dy+coef[2]*dz+coef[3]==0 ) return 1; //in the plane
	return 0;  //out of the plane
}
//===========================================================================
