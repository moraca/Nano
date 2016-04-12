//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	GCH.cpp
//OBJECTIVE:	A class for a hybrid Graphene nanoplatelet-Carbon Nanotubes particles
//AUTHOR:		Angel Mora;
//E-MAIL:			angel.mora@kaust.edu.sa	;
//====================================================================================

#include "GCH.h"

//Constructor that initializes the graphene nanoplatelet geometry
GCH::GCH(double len_x, double wid_y, double thick_z)
{
    //Geometry of the GNP
    gnp.len_x = len_x;
    gnp.wid_y = wid_y;
    gnp.hei_z = thick_z;    
}
//This function determines the separation between two hybrid particles
double GCH::Distance_to(GCH hybrid){
    return 0;
}

//This function determines if the graphene of the two particles are overlapping
int GCH::Overlapping_g(GCH hybrid){
    return 1;
}