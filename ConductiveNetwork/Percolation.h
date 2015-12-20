//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Percolation.h
//OBJECTIVE:	Determine which clusters percolate and in which directions
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef __Nanocode_clean__Percolation__
#define __Nanocode_clean__Percolation__

#include <iostream>

#include "Input_Reader.h"

//---------------------------------------------------------------------------
class Percolation
{
public:
    //Data Member
    vector<int> family; //This determines the family to which a cluster belongs to
    //0 for x-x; 1 for y-y; 2 for z-z; 3 for x-x and y-y; 4 for x-x and z-z; 5 for y-y and z-z; 6 for the three directions
    
    //Constructor
    Percolation(){};
    
    //Member Functions
    int Determine_percolating_clusters(const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const vector<vector<int> > &boundary_cnt, const vector<int> &labels, const vector<int>  &labels_labels, const vector<int>  &label_map, vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated, const int &window);
    int Cluster_CNT_percolation(const vector<vector<int> > &boundary_cnt, const vector<int> &labels, const vector<int>  &labels_labels, const vector<int>  &label_map, vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated);
    void Fill_percolation_flags_all_directions(const vector<vector<int> > &boundary_cnt, vector<vector<short int> > &perc_flag, const vector<int> &labels, const vector<int>  &labels_labels, const vector<int>  &label_map);
    void Fill_percolation_flags_single_direction(const vector<int> &boundary_vector, int boundary_number, vector<vector<short int> > &perc_flag, const vector<int> &labels, const vector<int>  &labels_labels, const vector<int>  &label_map);
    int Find_root(int L, const vector<int> &labels_labels);
    int Check_percolation_all_clusters(const vector<vector<int> > &boundary_cnt, vector<vector<short int> > &perc_flag, vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated);
    int Check_percolation_single_cluster(const vector<short int> &cluster_flag, int &family);
    int Single_CNT_percolation(const vector<vector<int> > &boundary_cnt, const vector<int> &labels, const vector<int>  &labels_labels, vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated, const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const int &window);
    void Fill_percolation_flags_all_directions_CNT_percoaltion(const vector<vector<int> > &boundary_cnt, vector<vector<short int> > &perc_flag, const vector<int> &labels_iso, int px, int py, int pz);
    void Fill_percolation_flags_single_direction_CNT_percoaltion(const vector<int> &boundary_vector, int boundary_number, vector<vector<short int> > &perc_flag, const vector<int> &labels_iso);
    int Check_percolation_CNTs(const vector<vector<int> > &boundary_cnt, vector<vector<short int> > &perc_flag, vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated);
    int Check_percolation_single_CNT(vector<short int> cluster_flag, int &family);
    
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
