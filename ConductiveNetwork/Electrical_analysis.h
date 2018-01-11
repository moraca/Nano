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
    
    int Perform_analysis_on_clusters(const int &avoid_resistance_flag, const vector<int> &family, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &point_list_gnp, const struct Geom_RVE &window, const struct Electric_para &electric_param, const struct Cutoff_dist &cutoffs, vector<GCH> &hybrid_particles, vector<vector<long int> > &all_dead_indices, vector<vector<long int> > &all_indices, vector<vector<int> > &gnp_dead_indices, vector<vector<int> > &gnp_indices);
    int Electrical_analysis_along_each_percolated_direction (const int &R_flag, const int &n_cluser, const vector<int> &family, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, Backbone_Network *Backbonet, const struct Electric_para &electric_param, const struct Cutoff_dist &cutoffs, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<long int> > &contacts_initial, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &point_list_gnp, vector<GCH> &hybrid_particles, vector<vector<double> > &paralel_resistors);
    int Vector_of_directions(const int &family, vector<int> &directions);
    int Convert_index_to_structure(const vector<int> &cluster, const vector<vector<long int> > &indices, vector<vector<long int> > &structure, vector<int> &backbone_cnts);
    int Update_hybrids(const vector<int> &cluster_gch, const vector<vector<long int> > &structure, const vector<vector<long int> > &backbone_structure, vector<GCH> &hybrid_particles);
    int Update_vectors_for_hoko_cutwins(const int &n_cnts, const int &n_gnps, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, Hoshen_Kopelman *HoKo_Re, Cutoff_Wins *Cutwins_Re, const vector<vector<long int> > &contacts_initial, const vector<vector<long int> > &structure, const vector<int> &backbone_cnts, const vector<int> &percolated_gnps);
    int Calculate_parallel_resistor(const int &direction, Direct_Electrifying * DEA, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<int> > &clusters_cnt, const vector<vector<int> > &boundary_cnt, const vector<Point_3D> &point_list_gnp, const vector<GCH> &hybrid_particles, const vector<vector<int> > &clusters_gnp, const vector<vector<int> > &boundary_gnp, const vector<vector<short int> > &boundary_flags_gnp, const struct Electric_para &electric_param, vector<vector<double> > &paralel_resistors);
    double Current_of_element_in_boundary(const long int &P1, const long int &P2, const double &radius, Direct_Electrifying *DEA, const struct Electric_para &electric_param, const vector<Point_3D> &point_list);
    double Current_of_edges_in_boundary(const int &side, Direct_Electrifying *DEA, const struct Electric_para &electric_param, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<Point_3D> &point_list_gnp, const GCH &hybrid);
    int Calculate_matrix_resistances(const double &matrix_resistivity, const struct Geom_RVE &window, vector<double> &matrix_resistances);
    int Calculate_resistances_and_resistivities(const struct Geom_RVE &window, const vector<double> &matrix_resistances, const vector<vector<double> > &paralel_resistors, vector<double> &resistors, vector<double> &resistivities);

private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================