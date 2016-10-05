//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Electrical_analysis.h
//OBJECTIVE:	Extract the backbone and calculate the electrical resistance
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef ELECTRICAL_ANALYSIS_h
#define ELECTRICAL_ANALYSIS_h

#include "Input_Reader.h"
#include "Backbone_Network.h"
#include "Cutoff_Wins.h"
#include "Direct_Electrifying.h"
#include "Hoshen_Kopelman.h"
#include "Printer.h"

//-------------------------------------------------------
class Electrical_analysis
{
public:
    //Data Member
    vector<double> resistors;
    
    int Perform_analysis_on_clusters(const int &iteration, const vector<int> &family, Hoshen_Kopelman *HoKo, const vector<vector<long int> > &structure, const vector<vector<short int> > &boundary_flags, const vector<vector<int> > &boundary_cnt, const vector<Point_3D> &point_list, const vector<double> &radii, const struct Geom_RVE &geom_rve, const struct Electric_para &electric_param, const struct Cutoff_dist &cutoffs, vector<GCH> &hybrid_particles, vector<double> &families_lengths, vector<double> &branches_lengths, vector<vector<long int> > &all_dead_indices, vector<vector<long int> > &all_indices, vector<vector<int> > &gnp_dead_indices, vector<vector<int> > &gnp_indices);
    int Vector_of_directions(const int &family, vector<int> &directions);
    int Convert_index_to_structure(const vector<int> &cluster, const vector<vector<long int> > &indices, vector<vector<long int> > &structure, vector<int> &backbone_cnts);
    int Update_hybrids(const vector<int> &cluster_gch, const vector<vector<long int> > &structure, const vector<vector<long int> > &backbone_structure, vector<GCH> &hybrid_particles);
    int Calculate_parallel_resistor(const int &direction, Direct_Electrifying * DEA, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<int> > &boundary_cnt, const struct Electric_para &electric_param, vector<vector<double> > &paralel_resistors);
    double Current_of_element_in_boundary(const long int &P1, const long int &P2, const double &radius, Direct_Electrifying *DEA, const struct Electric_para &electric_param, const vector<Point_3D> &point_list);
    int Calculate_matrix_resistances(const double &matrix_resistivity, const struct Geom_RVE &geom_rve, vector<double> &matrix_resistances);
    int Calculate_resistances(const vector<double> &matrix_resistances, const vector<vector<double> > &paralel_resistors, vector<double> &resistors);

private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================