//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Backbone_Network.h
//OBJECTIVE:	To determine the backbone network and dead branches in the percolation network
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef BACKBONE_NETWORK_H
#define BACKBONE_NETWORK_H

#include "Input_Reader.h"

//-------------------------------------------------------
class Backbone_Network
{
public:
    //Data Member
    vector<vector<long int> > dead_indices;
    vector<vector<long int> > percolated_indices;
    
    //Constructor
    Backbone_Network(){};
    
    //Member Functions
    int Determine_backbone_network(const int &family, const vector<int> &cluster, const vector<double> &voltages, const vector<int> &LM_matrix, const vector<vector<long int> > &elements, const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, vector<double> &families_lengths, vector<double> &branches_lengths);
    int Find_dead_branches(const vector<double> &voltages, const vector<int> &cluster, const vector<int> &LM_matrix, const vector<vector<long int> > &elements, const vector<vector<long int> > &structure);
    double Zero_voltage(const vector<double> &voltages, const vector<int> &cluster, const vector<int> &LM_matrix, const vector<vector<long int> > &elements);
    int Calculate_lengths(const int &family, const vector<Point_3D> &points_in, vector<double> &families_lengths, vector<double> &branches_lengths);
    double Segment_length(long int index1, long int index2, const vector<Point_3D> &points_in);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================