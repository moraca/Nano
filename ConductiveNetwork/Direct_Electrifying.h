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

//-------------------------------------------------------
class Direct_Electrifying
{
public:
    //Data Member
    vector<double> voltages;
    double resistance;
    vector<vector<long int> > elements; //This vector will store the elements. So it is needed to trim the CNTs
    vector<int> LM_matrix;//Local mapping matrix. It maps from point number to node number. It is also used to calculate the currents
    
    //Constructor
    Direct_Electrifying(){};
    
    //Member Functions
    int Calculate_voltage_field(const vector<vector<long int> > &structure, vector<vector<long int> > &contacts_point, const vector<vector<short int> > &boundary_flags, const vector<int> &cluster, const vector<double> &radii, int family, const struct Electric_para &electric_param);
    int Get_LM_matrix(const vector<vector<long int> > &structure, const vector<vector<long int> > &contacts_point, const vector<vector<short int> > &boundary_flags, const vector<int> &cluster, int family, int &global_nodes, vector<int> &LM_matrix, vector<vector<long int> > &elements);
    int Is_in_relevant_boundary(int family, short int boundary_node);
    void Fill_sparse_stiffness_matrix(const vector<vector<long int> > &structure, vector<vector<long int> > &contacts_point, const vector<int> &cluster, int nodes, const vector<int> &LM_matrix, const vector<vector<long int> > &elements, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, vector<double> &R);
    void Fill_2d_matrices(const vector<vector<long int> > &elements, const vector<int> &cluster, const vector<int> &LM_matrix, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point);
    void Check_for_other_elements(const vector<int> &LM_matrix, long int P1, long int node1, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point);
    void Add_elements_to_sparse_stiffness(const vector<vector<long int> > &contacts_point, const vector<int> &LM_matrix, long int P1, long int node1, long int k, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal);
    void Remove_from_vector(long int num, vector<long int> &vec);
    void From_2d_to_1d_vectors(const vector<vector<long int> > &col_ind_2d, const vector<vector<double> > &values_2d, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, vector<double> &R);
    void Solve_DEA_equations_CG_SSS(long int nodes, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, const struct Electric_para &electric_param, vector<double> &R);
    void Conjugate_gradient(long int nodes, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, vector<double> &R, vector<double> &P);
    double V_dot_v(vector<double> A, vector<double> B);
    void spM_V_SSS(vector<double> V, vector<long int> rowptr, vector<long int> colind, vector<double> diagonal, vector<double> values, vector<double> &R);
    void V_plus_aW(vector<double> W, double a, vector<double> &V);
    void W_plus_aV(vector<double> W, double a, vector<double> &V);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================