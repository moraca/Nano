//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Backbone_Network.h
//OBJECTIVE:	To determine the backbone network and dead branches in the percolation network
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef BACKBONE_NETWORK_H
#define BACKBONE_NETWORK_H

#include "Direct_Electrifying.h"
#include "Hoshen_Kopelman.h"
#include "Input_Reader.h"
#include "Printer.h"
#include <algorithm>


//-------------------------------------------------------
class Backbone_Network
{
public:
    //Data Member
    vector<vector<long int> > dead_indices;
    vector<vector<long int> > percolated_indices;
    vector<int> percolated_gnps;
    vector<int> dead_gnps;
    
    //Constructor
    Backbone_Network(){};
    
    //Member Functions
    int Determine_backbone_network(const int &family, const int &n_cluster, const int &R_flag, const int &tecplot_flag, Direct_Electrifying *DEA, Hoshen_Kopelman *HoKo, const Electric_para &electric_param, const Cutoff_dist &cutoffs, const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_in_gnp, const vector<GCH> &hybrid_particles, vector<double> &families_lengths, vector<double> &branches_lengths, vector<vector<long int> > &all_dead_indices, vector<vector<long int> > &all_indices, vector<vector<int> > &gnp_dead_indices, vector<vector<int> > &gnp_indices);
    double Zero_current(const int &n_cluster, const int &R_flag, Direct_Electrifying *DEA, Hoshen_Kopelman *HoKo, const struct Electric_para &electric_param, const Cutoff_dist &cutoffs, const vector<Point_3D> &point_list, const vector<Point_3D> &points_in_gnp, const vector<double> &radii, const vector<GCH> &hybrid_particles, vector<vector<double> > &currents_cnt, vector<vector<double> > &currents_gnp);
    double Voltage_difference(const long int &P1, const long int &P2, const vector<int> &LM_matrix, const vector<double> &voltages);
    int Find_dead_branches(const double &zero_cutoff, const vector<vector<double> > &currents_cnt, const vector<int> &cluster, vector<vector<long int> > &elements);
    int Find_dead_gnps(const double &zero_cutoff, const vector<vector<double> > &currents_gnp, const vector<int> &cluster_gch);
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    int Calculate_lengths(const int &family, const vector<Point_3D> &points_in, vector<double> &families_lengths, vector<double> &branches_lengths);
    double Segment_length(long int index1, long int index2, const vector<Point_3D> &points_in);
    void Add_indices_to_global_vectors(const int &family, vector<vector<long int> > &all_dead_indices, vector<vector<long int> > &all_indices);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================