//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Direct_Electrifying.h
//OBJECTIVE:	The direct eletrifying algorithm (C.Y. Li and T.W. Chou, Int. J. Mod Phys C, 20, 2009, 423-33.)
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef DIRECTELECTRIFYING_H
#define DIRECTELECTRIFYING_H

#include "Input_Reader.h"
#include "Printer.h"
#include "Triangulation.h"
#include <algorithm>

//-------------------------------------------------------
class Direct_Electrifying
{
public:
    //Data Member
    vector<double> voltages;
    vector<double> resistances;
    vector<vector<long int> > elements; //This vector will store the elements. It is needed to trim the CNTs
    vector<vector<long int> > elements_tunnel; //This vector will store the tunnel elements. It is needed to calculate the zero-current cutoff
    vector<int> LM_matrix;//Local mapping matrix. It maps from point number to node number. It is also used to calculate the currents
    vector<vector<vector<int> > > boundary_node_map;
    int reserved_nodes; //This is the number of nodes assigned to the boudaries with prescribed voltage
    
    //Constructor
    Direct_Electrifying(){};
    
    //Member Functions
    int Calculate_voltage_field(const int &family, const int &R_flag, const vector<vector<long int> > &structure, vector<vector<long int> > &contacts_point, const vector<vector<short int> > &boundary_flags, const vector<Point_3D> &point_list, const vector<int> &cluster, const vector<int> &cluster_gch, const vector<double> &radii, const struct Electric_para &electric_param, const struct Cutoff_dist &cutoffs, vector<GCH> &hybrid_particles);
    void Initialize_boundary_node_map();
    int Get_LM_matrix(const vector<vector<long int> > &structure, const vector<vector<long int> > &contacts_point, const vector<vector<short int> > &boundary_flags, const vector<int> &cluster, const int &family, int &global_nodes, vector<int> &LM_matrix, vector<vector<long int> > &elements);
    int Get_global_nodes(const int &family);
    int Is_in_relevant_boundary(int family, short int boundary_node);
    void Add_point_to_LM_matrix(long int P, int family, const vector<vector<short int> > &boundary_flags, int &global_nodes, vector<int> &LM_matrix);
    int Get_boundary_node(const vector<short int> &boundary_flag, const int &family);
    int Fill_sparse_stiffness_matrix(const int &R_flag, const vector<vector<long int> > &structure, vector<vector<long int> > &contacts_point, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &cluster_gch, const vector<int> &cluster, const int &nodes, const double &d_vdw, const vector<int> &LM_matrix, const vector<vector<long int> > &elements, vector<GCH> &hybrid_particles, vector<vector<double> > &KEFT, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, const struct Electric_para &electric_param);
    int Fill_2d_matrices(const vector<vector<long int> > &elements, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &cluster, const vector<int> &LM_matrix, const struct Electric_para &electric_param, const double &d_vdw, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point);
    int Fill_2d_matrices_gch(const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &cluster_gch, const vector<int> &LM_matrix, const struct Electric_para &electric_param, vector<GCH> &hybrid_particles, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal);
    int Fill_2d_matrices_unit_resistors(const vector<vector<long int> > &elements, const vector<Point_3D> &point_list, const vector<int> &cluster, const vector<int> &LM_matrix, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point);
    int Fill_2d_matrices_gch_unit_resistors(const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<int> &cluster_gch, const vector<int> &LM_matrix, vector<GCH> &hybrid_particles, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal);
    void Check_for_other_elements(const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &LM_matrix, const struct Electric_para &electric_param, const double d_vdw, const long int &P1, const long int &node1, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point);
    void Check_for_other_unit_elements(const vector<Point_3D> &point_list, const vector<int> &LM_matrix, const long int &P1, const long int &node1, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point);
    void Add_elements_to_sparse_stiffness(const long int &node1, const long int &node2, const double &Re, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal);
    double Calculate_resistance_cnt(const vector<Point_3D> &point_list, const long int &P1, const long int &P2, const double &radius, const double &resistivity);
    double Calculate_resistance_tunnel(const vector<double> &radii, const struct Electric_para &electric_param, const Point_3D &P1, const Point_3D &P2, const double &d_vdw);
    double Calculate_resistance_gnp(const Point_3D &P1, const Point_3D &P2, const double &rad1, const double &rad2, const GCH &hybrid, const struct Electric_para &electric_param);
    void Remove_from_vector(long int num, vector<long int> &vec);
    int Check_repeated_col_ind_2d(const int &nodes, const vector<vector<long int> > &structure, vector<vector<long int> > &contacts_point, const vector<Point_3D> &point_list, const vector<int> &LM_matrix, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d);
    void Find_repeated_elements(const vector<long int> &vector_in, vector<long int> &elements, vector<int> &indices);
    int Export_matlab_sparse_matrix(const vector<vector<long int> > &col_ind_2d, const vector<vector<double> > &values_2d, const vector<double> &diagonal, const string &filename);
    void From_2d_to_1d_vectors(const vector<vector<long int> > &col_ind_2d, const vector<vector<double> > &values_2d, vector<vector<double> > &KEFT, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal);
    int Solve_DEA_equations_CG_SSS(const int &R_flag, long int nodes, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, const struct Electric_para &electric_param, vector<vector<double> > &KEFT);
    void Get_voltage_vector(const double &nodes, vector<double> &voltages);
    void Conjugate_gradient(long int nodes, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, vector<double> &R, vector<double> &P);
    double V_dot_v(const vector<double> &A, const vector<double> &B);
    void spM_V_SSS(const vector<double> &V, const vector<long int> &rowptr, const vector<long int> &colind, const vector<double> &diagonal, const vector<double> &values, vector<double> &R);
    void V_plus_aW(const vector<double> &W, const double &a, vector<double> &V);
    void W_plus_aV(const vector<double> &W, const double &a, vector<double> &V);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================