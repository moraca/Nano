//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Background_grid.h
//OBJECTIVE:	
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef BACKGROUND_GRID_H
#define BACKGROUND_GRID_H

#include <iostream>

#include "Input_Reader.h"

//-------------------------------------------------------
class Background_grid
{
public:
    //Data Member
    vector<vector<int> > sectioned_domain_cnt; //Shell sub-regions to identify CNTs that need to be trimmed
    vector<vector<long int> > structure; //Vector with global number of points in the same form as the points_in vector
    
    //Constructor
    Background_grid(){};
    
    //Member Functions
    int Generate_background_grids(const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, vector<Point_3D> &points_in);
    int Fill_structure_and_shell(const struct Geom_RVE &sample, vector<Point_3D> &points_in);
    int Add_to_shell(struct Geom_RVE sample, Point_3D point, vector<vector<int> > &sectioned_domain_cnt);
    int Find_shell(double x_in, double x_min, double len_x, double dx);
    
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================