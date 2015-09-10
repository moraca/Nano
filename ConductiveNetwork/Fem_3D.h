//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Fem_3D.h
//OBJECTIVE:	A class definition for finite element method (FEM)
//AUTHOR:		Fei Han;
//E-MAIL:			fei.han@kaust.edu.sa	;
//====================================================================================
#ifndef FEM_3D_H
#define FEM_3D_H

#include<iostream>
#include<cmath>
#include<vector>
using namespace std;

//---------------------------------------------------------------------------
//The definition for node
class Node
{	
	public:
		double x, y, z;		//the coordinates of a node

		int type;				//the type of node: 0 inner node; 1 a node on the surface ; 2 a node on the boundary; 3 a node at the corner
		vector<int> relative_nods;		//a vector for the relative nodes of this node (using the number of node)
		vector<int> relative_eles;		//a vector for the relative elements of this node (using the number of element)

		//Constructor
		Node(const double ix=0, const double iy=0, const double iz=0);
		
		//Member Functions
		double	distance_to(const Node& n);		//to calculate the distance with another node
};
//---------------------------------------------------------------------------
//The definition for element
class Element
{
	public:	
		int type;	//the shape and order of an element (using 3 numbers: ijk, i denotes dimensions; j denotes the number of nodes in one element; k denotes the order of shape function)
						//for example, 121: one dimension two nodes (segment) linear shape function; 231: two dimension thress nodes (triangular) linear shape function; 241: two dimension four nodes (quadrilateral) linear shape function
						//341: three dimensions four nodes (tetrahedron) linear shape function; 361: thress dimensions six nodes (prism) linear shape function; 381: thress dimensions eight nodes (hexahedron) linear shape function
		int mat;		//the material property of element
		vector<int> nodes_id;
		vector<int> relative_eles;	//a vector for the relative elements of this element (using the number of element)
		vector<int> relative_cnts;	//a vector for the relative cnts of this element (using the number of points in nanotubes)
		
		//Constructor
		Element(){};
		
		//Member Functions
};
//---------------------------------------------------------------------------
//The definition for hexahedron element
class Hexahedron
{
	public:
        int nodesId[8];	 //eight nodes;
        int materialId;
        int facesId[6];

		//Constructor
		Hexahedron(){};
        Hexahedron(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
		
		//Member Functions
};

#endif
//===========================================================================
