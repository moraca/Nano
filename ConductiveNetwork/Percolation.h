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
    int Determine_percolated_clusters(const int &window, const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const struct GNP_Geo &gnps, const vector<vector<int> > &boundary_cnt, const vector<int> &labels_cnt, const vector<vector<int> > &boundary_gnp, const vector<int> &labels_gnp, vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated_cnt, vector<vector<int> > &clusters_gnp, vector<vector<int> > &isolated_gnp);
    int Cluster_percolation(const int &ignore, const vector<vector<int> > &boundary_1, const vector<int> &labels_1, vector<vector<int> > &clusters_1, vector<vector<int> > &isolated_1, const vector<vector<int> > &boundary_2, const vector<int> &labels_2, vector<vector<int> > &clusters_2, vector<vector<int> > &isolated_2);
    void Fill_percolation_flags_all_directions(const vector<vector<int> > &boundary, vector<vector<short int> > &perc_flag, const vector<int> &labels);
    void Fill_percolation_flags_single_direction(const vector<int> &boundary_vector, int boundary_number, vector<vector<short int> > &perc_flag, const vector<int> &labels);
    int Check_percolation_all_clusters(const int &ignore, vector<vector<short int> > &perc_flag, vector<vector<int> > &clusters_1, vector<vector<int> > &isolated_1, vector<vector<int> > &clusters_2, vector<vector<int> > &isolated_2, vector<int> &family);
    int Check_percolation_single_cluster(const vector<short int> &cluster_flag, int &family);
    int Single_particle_percolation(const int &window, const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const struct GNP_Geo &gnps, const vector<vector<int> > &boundary, const vector<int> &labels, vector<vector<int> > &clusters, vector<vector<int> > &isolated);
    int Determine_direction_flags(const int &window, const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const struct GNP_Geo &gnps, int &px, int &py, int &pz);
    int Determine_isolated_labels(const vector<vector<int> > &isolated, vector<int> &labels_iso);
    void Fill_single_particle_percolation_flags_all_directions(const vector<vector<int> > &boundary, vector<vector<short int> > &perc_flag, const vector<int> &labels_iso, const int &px, const int &py, const int &pz);
    void Fill_single_particle_percolation_flags_single_direction(const vector<int> &boundary_vector, const int &boundary_number, vector<vector<short int> > &perc_flag, const vector<int> &labels_iso);
    int Check_particles_for_percolation(vector<vector<short int> > &perc_flag, vector<vector<int> > &clusters, vector<vector<int> > &isolated);
    int Check_percolation_single_particle(const vector<short int> &cluster_flag, int &family);
    
    
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
