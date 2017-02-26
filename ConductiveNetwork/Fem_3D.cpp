//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Fem_3D.h
//OBJECTIVE:	Member functions in a class of finite element method (FEM)
//AUTHOR:		Fei Han;
//E-MAIL:			fei.han@kaust.edu.sa	;
//====================================================================================

#include "Fem_3D.h"

//=========================================
//Constructor function for node class
Node::Node(const double ix, const double iy, const double iz)
{
	x=ix;
	y=iy;
	z=iz;
}
//---------------------------------------------------------------------------
//To calculate the distance
double Node::distance_to(const Node& n)
{
	return(sqrt( (x-n.x)*(x-n.x)+(y-n.y)*(y-n.y)+(z-n.z)*(z-n.z)) );
}

//=========================================
//Constructor function for element class

//=========================================
//Constructor function for hexahedron class
Hexahedron::Hexahedron(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{
	nodesId[0] = n1; nodesId[1] = n2;
	nodesId[2] = n3; nodesId[3] = n4;
	nodesId[4] = n5; nodesId[5] = n6;
	nodesId[6] = n7; nodesId[7] = n8;
}
//=================================================================
