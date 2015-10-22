//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Background_vectors.h
//OBJECTIVE:	Using nested shells on the background to mark the CNTs for trimming faster in each observation window
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef BACKGROUND_VECTORS_H
#define BACKGROUND_VECTORS_H

#include<iomanip>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include "Input_Reader.h"

//-------------------------------------------------------
class Background_vectors
{
public:
    //Data Member
    
    //Constructor
    Background_vectors(){};
    
    //Member Functions
    int Generate_shells_and_structure(const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const vector<Point_3D> &points_out, vector<vector<int> > &shells_cnt);
    int Add_to_shell(const struct Geom_RVE &sample, const Point_3D &point, vector<vector<int> > &shells_cnt);
    int Find_shell(double x_in, double x_0, double len_x, double dx, double win_min_x, double win_max_x, vector<vector<int> > &shells_cnt);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================