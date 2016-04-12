//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	GCH.h
//OBJECTIVE:	A class for a hybrid Graphene nanoplatelet-Carbon Nanotubes particles
//AUTHOR:		Angel Mora;
//E-MAIL:			angel.mora@kaust.edu.sa	;
//====================================================================================


#ifndef GCH_H
#define GCH_H

#include<cmath>
#include<stdlib.h>
#include<vector>
#include "Hns.h"
#include "MathMatrix.h"
#include "Geometry_3D.h"

//---------------------------------------------------------------------------
//The definition for points in 3D
class GCH
{
public:
    //Data Member
    cuboid gnp;
    Point_3D center;
    MathMatrix rotation;
    vector<Point_3D> top_cnts, bottom_cnts, corners;
    vector<vector<int> > structure_top, structure_bottom;
    vector<int> cnts, cnts_top, cnts_bottom;
    int flag;
    
    //Constructor
    GCH(){};
    GCH(double len_x, double wid_y, double thick_z);
    
    //Member Functions
    double Distance_to(GCH hybrid);
    int Overlapping_g(GCH hybrid);

};

#endif
//===========================================================================