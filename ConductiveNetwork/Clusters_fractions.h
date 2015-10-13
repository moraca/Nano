//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Clusters_fractions.h
//OBJECTIVE:
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef __Nanocode_clean__Clusters_fractions__
#define __Nanocode_clean__Clusters_fractions__

#include <iostream>

#include "Input_Reader.h"

//-------------------------------------------------------
class Clusters_fractions
{
public:
    //Data Member
    vector<vector<int> > sectioned_domain_cnt; //Shell sub-regions to identify CNTs that need to be trimmed
    vector<vector<long int> > structure; //Vector with global number of points in the same form as the points_in vector
    
    //Constructor
    Clusters_fractions(){};
    
    //Member Functions
    int Calculate_fractions(const vector<vector<long int> > &structure, const vector<int> &cnts_inside, const vector<vector<int> > &isolated, const vector<Point_3D> points_in, vector<double> &families_lengths, vector<double> &branches_lengths, vector<double> fractions);
    double CNT_length(const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, int CNT);
    void Append_1d_vector_to_file(const vector<double> &list, const string &filename);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================