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
#include "GenNetwork.h"
#include "Printer.h"

//-------------------------------------------------------
class Clusters_fractions
{
public:
    //Data Member
    
    //Constructor
    Clusters_fractions(){};
    
    //Member Functions
    int Calculate_fractions(const struct Geom_RVE &sample, const vector<int> &cnts_inside, const vector<int> &gnps_inside, const vector<vector<long int> > &structure, const vector<Point_3D> points_in, const vector<double> &radii, const vector<vector<int> > &isolated, const vector<GCH> &hybrid_particles, const vector<vector<int> > &isolated_gnp, const vector<vector<long int> > &all_dead_indices, const vector<vector<long int> > &all_percolated_indices, const vector<vector<int> > &all_dead_gnps, const vector<vector<int> > &all_percolated_gnps);
    int Calculate_volumes_cnt(const vector<int> &cnts_inside, const vector<vector<long int> > &structure, const vector<Point_3D> points_in, const vector<double> &radii, const vector<vector<int> > &isolated, const vector<vector<long int> > &all_dead_indices, const vector<vector<long int> > &all_percolated_indices, vector<double> &volumes_cnt, vector<double> &dead_branches, double &cnt_volume);
    double CNT_cluster_volume(const vector<int> &cluster, const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const vector<double> &radii);
    double CNT_length(const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const int &CNT);
    double CNT_indices_volume(const vector<long int> &indices, const vector<Point_3D> &points_in, const vector<double> &radii);
    int Calculate_volumes_gnp(const struct Geom_RVE &sample, const vector<int> &gnps_inside, const vector<GCH> &hybrid_particles, const vector<vector<int> > &isolated_gnps, const vector<vector<int> > &all_dead_gnps, const vector<vector<int> > &all_percolated_gnps, vector<double> &volumes_gnp, vector<double> &dead_gnps, double &gnp_volume);
    double GNP_cluster_volume(const struct cuboid &gvcub, const vector<int> &cluster, const vector<GCH> &hybrid_particles);
    int Append_to_volume_and_fraction_files(const vector<double> &volumes_vec, const double &volume, const string &filename);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================