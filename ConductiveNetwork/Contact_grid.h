//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Background_grid.h
//OBJECTIVE:
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef CONTACT_GRID_H
#define CONTACT_GRID_H

#include <iostream>

#include "Input_Reader.h"

//-------------------------------------------------------
class Contact_grid
{
public:
    //Data Member
    vector<vector< long int> > sectioned_domain; //
    
    //Constructor
    Contact_grid(){};
    
    //Member Functions
    int Generate_contact_grid(const vector<vector<long int> > &structure, const vector<int> &cnts_inside, vector<Point_3D> &points_in, const struct Geom_RVE &sample, struct Cutoff_dist &cutoffs, const struct Nanotube_Geo &cnts, double lx, double ly, double lz);
    int calculate_t(int a, int b, int c, int sx, int sy);
    
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================
