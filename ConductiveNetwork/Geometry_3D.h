//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Geometry_3D.h
//OBJECTIVE:	A class for basic geometic elements in 3D such as point, line and plane 
//AUTHOR:		Fei Han;
//E-MAIL:			fei.han@kaust.edu.sa	;
//====================================================================================

#ifndef GEOMETRY_3D_H
#define GEOMETRY_3D_H

#include<cmath>
#include<stdlib.h>
#include<vector>
#include "Hns.h"
#include "MathMatrix.h"
using namespace hns;

const int MAX_INT = 65536; //2^16 for calculating a random number

//---------------------------------------------------------------------------
//The definition for points in 3D
class Point_3D
{
	public:
		//Data Member
		double x, y, z;
		int flag;
    
   		//Constructor
		Point_3D(){};
		Point_3D( double px, double py, double pz );
    
		//Member Functions
		Point_3D operator+( Point_3D &pt );
		Point_3D operator+( const Point_3D &pt );
		Point_3D operator+( double d );
		Point_3D operator-( double d );
		Point_3D operator-( Point_3D &pt );
		Point_3D operator-( const Point_3D &pt );
		Point_3D operator*( double d );
		Point_3D operator/( double d );
		bool operator==( Point_3D &pt );
		bool operator!=( Point_3D &pt );
		double distance_to(const Point_3D &pt)const;
		double distance_to(const double &px, const double &py,  const double &pz)const;
		Point_3D cross(Point_3D &point);
		double dot(Point_3D &point);
};
//---------------------------------------------------------------------------
//The definition for lines (segements) in 3D
class Line_3D
{
	public:
		//Data Member
		 Point_3D point[2];	//the coordinates of two endpoints of a segment
		 double len;				//the length of a segment
		 bool virtual_line;		//to mark if it is a virtual(false) segment (false: reduced to a point; true: a real segment)

		//Constructor
		Line_3D(Point_3D p0, Point_3D p1);
		
		//Member Functions
		double length();				//the length of a segment
        double distance_point_to_line(const Point_3D *point_temp)const;			//the distance from a point to a line
        double distance_point_to_line(const Point_3D &point_temp)const;		//the distance from a point to a line
        double distance_point_to_line(const double dx, const double dy, const double dz)const;  //the distance from a point to a line
		int contain(const Point_3D &point_temp)const;    //to judge if a point is in a segment

	private:			
		//Data Member
        double xm, yn, zl ;       //coefficients for a line: (x-x0)/xm=(y-y0)/yn=(z-z0)/zl
};
//---------------------------------------------------------------------------
//The definition for plane in 3D
class Plane_3D
{
	public:

		//Data Member
		double coef[4];			//four coefficients for an equation of plane, i.e., ax+by+cz+d=0
		bool virtual_plane;		///to mark if it is a virtual(false) plane (false: a virtual plane (normal vector is (0,0,0)); true: a real plane)

		//Constructor
		Plane_3D(){};
		Plane_3D(double para[]);
		
		//Member Functions
		int contain(const Point_3D &point_temp)const;										//to judge if a point is contained in this plane
		int contain(const double dx, const double dy, const double dz)const;		//to judge if a point is contained in this plane
};
//---------------------------------------------------------------------------
//Structural data for an ellipsoid
struct elliparam		
{	
	double x, y, z;		//the center point (x,y,z) of an ellipsoid
	double a, b, c;		//the long, middle and short axis of an ellipsoid
	double alpha1, alpha2, alpha3;	//[(alpha1,beta1,gamma1),(alpha2,beta2,gamma2),(alpha3,beta3,gamma3)] 
	double beta1, beta2, beta3;		//are 9 angles between three axes of ellipsoid (a,b,c) with three coordinate axes (ox,oy,oz).
	double gamma1, gamma2, gamma3;
};
//---------------------------------------------------------------------------
//Structural data for a cuboid
struct cuboid	
{	
	Point_3D poi_min;					//Define an origin point for a cubid
	double len_x, wid_y, hei_z;		//Define length, width and height for a cuboid
	double volume;							//Define the volume of a cuboid
};
#endif
//===========================================================================
