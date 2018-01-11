//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Direct_Eletrifying.cpp
//OBJECTIVE:	The direct eletrifying algorithm (C.Y. Li and T.W. Chou, Int. J. Mod Phys C, 20, 2009, 423-33.)
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Direct_Electrifying.h"

//Calculate the voltage values at contact points and endpoints
int Direct_Electrifying::Calculate_voltage_field(const int &family, const int &n_cluster, const int &R_flag, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &point_list_gnp, const struct Electric_para &electric_param, const struct Cutoff_dist &cutoffs, vector<GCH> &hybrid_particles)
{
    //First we need to prepare the matrices that the direct electrifying needs
    //The first matrix will be the local mapping (LM) matrix. This matrix maps from point number in the structure
    //to node number for the solver
    //The number of reserved nodes is calculated. These is the number of boundaries with prescribed voltage
    reserved_nodes = Get_global_nodes(family);
    //This variable is used to assing node numbers and after calling Get_LM_matrix, it will contain the number of nodes in the network
    int global_nodes = reserved_nodes;
    //Initialize the size of the LM matrix to be equal to the number of points
    LM_matrix.assign(HoKo->contacts_point.size(), -1);
    //Initialize the size of the LM matrix for GNP points to be equal to the number of GNP points
    LM_matrix_gnp.assign(point_list_gnp.size(), -1);
    
    //hout << "LM_matrix.size()="<<LM_matrix.size()<<" clusters_cnt.size()="<<HoKo->clusters_cnt.size()<<endl;
    //hout << "LM_matrix_gnp.size()="<<LM_matrix_gnp.size()<<" clusters_gch.size()="<<HoKo->clusters_gch.size()<<endl;
    //If there is a cluster wiht one CNT, then there is no need to do the DEA
    //This heppens when the CNT length is larger than one or more dimensions of the observation window
    if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluster].size() == 1) {
        return 1;
    }
    
    //Initialize the vector boundary_node_map
    Initialize_boundary_node_map();
    
    //Initialize the size of the elements matrix to be equal to the number of CNTs
    vector<long int> empty;
    elements.assign(structure.size(), empty);
    //Initialize the vector of GNP points to triangulate
    vector<vector<long int> > gnp_triangulation_points(structure_gnp.size(), empty);
    
    //Construct the LM matrices
    if (!HoKo->clusters_gch.size()) //If there are only CNTs
    {
        if (!Get_LM_matrix_cnts_only(family, n_cluster, HoKo, Cutwins, structure, global_nodes, LM_matrix, elements)) {
            hout << "Error in Calculate_voltage_field when calling Get_LM_matrix_cnts_only" << endl;
            return 0;
        }
        hout << "Get_LM_matrix_cnts_only"<<endl;
    }
    else //If there are mixed, hybrids or GNPs
    {
        if (!Get_LM_matrices(family, n_cluster, HoKo, Cutwins, structure, structure_gnp, global_nodes, LM_matrix, LM_matrix_gnp, elements, gnp_triangulation_points)) {
            hout << "Error in Calculate_voltage_field when calling Get_LM_matrices" << endl;
            return 0;
        }
        hout << "Get_LM_matrices"<<endl;
        
    }
    
    /*/
    Printer *P = new Printer;
    if (R_flag == 1) {
        P->Print_2d_vec(gnp_triangulation_points, "gnp_triangulation_points_R.txt");
        P->Print_2d_vec(Cutwins->boundary_gnp, "boundary_gnp_R.txt");
        P->Print_2d_vec(elements, "elements_R.txt");
        P->Print_1d_vec(LM_matrix_gnp, "LM_matrix_gnp_R.txt");
        P->Print_1d_vec(LM_matrix, "LM_matrix_R.txt");
        P->Print_1d_vec(HoKo->gnp_contacts, "gnp_contacts_R.txt");
    } else if (R_flag == 0) {
        P->Print_2d_vec(gnp_triangulation_points, "gnp_triangulation_points_unit.txt");
        P->Print_2d_vec(Cutwins->boundary_gnp, "boundary_gnp_unit.txt");
        P->Print_2d_vec(elements, "elements_unit.txt");
        P->Print_1d_vec(LM_matrix_gnp, "LM_matrix_gnp_unit.txt");
        P->Print_1d_vec(LM_matrix, "LM_matrix_unit.txt");
        P->Print_1d_vec(HoKo->gnp_contacts, "gnp_contacts_unit.txt");
    } else {
        hout << "Error in Fill_sparse_stiffness_matrix. Invalid R_flag:" << R_flag << ". Only 0 and 1 are valid flags." << endl;
        return 0;
    }
    delete P;//*/

    //hout << "global_nodes="<<global_nodes<<endl;
    //Variables for using the SSS for the sparse matrix
    vector<long int> col_ind, row_ptr;
    vector<double> values, diagonal;
    vector<vector<double> > KEFT;
    //With the LM matrix, now fill the sparse stiffness matrix
    if (!Fill_sparse_stiffness_matrix(R_flag, global_nodes, cutoffs.van_der_Waals_dist, n_cluster, HoKo, structure, elements, LM_matrix, point_list, radii, structure_gnp, gnp_triangulation_points, LM_matrix_gnp, point_list_gnp, hybrid_particles, KEFT, col_ind, row_ptr, values, diagonal, electric_param)) {
        hout << "Error in Calculate_voltage_field when calling Fill_sparse_stiffness_matrix" << endl;
        return 0;
    }
    
    //Clear vector to free memory
    gnp_triangulation_points.clear();
    
    //hout << "Solve_DEA_equations_CG_SSS"<<endl;
    //This is where the actual direct electrifying algorithm (DEA) takes place
    if (!Solve_DEA_equations_CG_SSS(R_flag, global_nodes, col_ind, row_ptr, values, diagonal, electric_param, KEFT)) {
        hout << "Error in Calculate_voltage_field when calling Solve_DEA_equations_CG_SSS" << endl;
        return 0;
    }
    
    return 1;
}
//
int Direct_Electrifying::Get_global_nodes(const int &family)
{
    if (family == 6) {
        //If family is 6, then a cluster percolates in the three directions. Hence we need 6 voltages: 0 to 5. 6 is the first available node
        return 6;
    } else if ( (3 <= family) && (family <= 5) ){
        //If family is 3, 4 or 5, then a cluster percolates in two directions. Hence we need 4 voltages: 0 to 3. 4 is the first available node
        return 4;
    } else {
        //If a family is 0, 1 or 2, then a cluster percolates in one direction. Hence we need 2 voltages: 0 to 1. 2 is the first available node
        return 2;
    }
}
//This function initializes the vector boundary_node_map. This can be initialized in the h-file using Xcode.
//Using the make file I get an error because it cannot be initialized in the h-file (for what I could understand)
//So instead of using an array, I use a vector and create this function to add the values
//The orginal initialization in the h-file was:
//
//The boundary_node_map is used to assign a reserved node number to the boudaries of the observation window depending on the
//number of directions in which a cluster percolates
//It has this form: boundary_node_map[family][direction][side]
//direction is 0 for x, 1 for y, or 2 for z
//side is 0 for x0,y0 or z0; or 1 for x1, y1 or z1
//The -1 are used to produce errors, so they can be used when debugging
//boundary_node_map[4][3][2] = {
//    {{ 0, 1}, { 2, 3}, {-1,-1}}, //family 3
//    {{ 0, 1}, {-1,-1}, { 2, 3}}, //family 4
//    {{-1,-1}, { 0, 1}, { 2, 3}}, //family 5
//    {{ 0, 1}, { 2, 3}, { 4, 5}}  //family 6
//};
void Direct_Electrifying::Initialize_boundary_node_map(){
    //Numerical elements
    //hout << "Numerical elements" << endl;
    vector<int> els_01;
    els_01.push_back(0);
    els_01.push_back(1);
    vector<int> els_23;
    els_23.push_back(2);
    els_23.push_back(3);
    vector<int> els_45;
    els_45.push_back(4);
    els_45.push_back(5);
    vector<int> els_n1n1;
    els_n1n1.push_back(-1);
    els_n1n1.push_back(-1);
    
    //Family vectors
    //hout << "Family vectors" << endl;
    vector<vector<int> > family_3;
    family_3.push_back(els_01);
    family_3.push_back(els_23);
    family_3.push_back(els_n1n1);
    vector<vector<int> > family_4;
    family_4.push_back(els_01);
    family_4.push_back(els_n1n1);
    family_4.push_back(els_23);
    vector<vector<int> > family_5;
    family_5.push_back(els_n1n1);
    family_5.push_back(els_01);
    family_5.push_back(els_23);
    vector<vector<int> > family_6;
    family_6.push_back(els_01);
    family_6.push_back(els_23);
    family_6.push_back(els_45);
    
    //Global vector for mapping
    //hout << "Global vector for mapping "<<boundary_node_map.size() << endl;
    boundary_node_map.clear();
    boundary_node_map.push_back(family_3);
    boundary_node_map.push_back(family_4);
    boundary_node_map.push_back(family_5);
    boundary_node_map.push_back(family_6);
    
    /*/
    hout << "Elements of global vector for mapping" << endl;
    for (int i = 0; i < (int)boundary_node_map.size(); i++) {
        for (int j = 0; j < (int)boundary_node_map[i].size(); j++) {
            for (int k = 0; k < (int)boundary_node_map[i][j].size(); k++) {
                hout << "boundary_node_map["<<i<<"]["<<j<<"]["<<k<<"]="<<boundary_node_map[i][j][k]<<endl;
            }
        }
    }//*/
}
//Build the LM matrix and the elements matrix
//By building the elements matrix in this step, I avoid to use the contacts_cnt_point vector that I used in previous versions
//Also, by building it at this stage I have the nodes in order
int Direct_Electrifying::Get_LM_matrices(const int &family, const int &n_cluster, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure, const vector<vector<long int> > &structure_gnp, int &global_nodes, vector<int> &LM_matrix, vector<int> &LM_matrix_gnp, vector<vector<long int> > &elements, vector<vector<long int> > &gnp_triangulation_points)
{
    //=========================================================
    //Mixed and GNP contacts pre-processing
    
    //Vector with flags indicating a CNT point has a contact with a GNP point
    vector<long int> cnt_contacts_point(HoKo->contacts_point.size(), 0);
    //Vector with flags indicating a GNP point has a contact with another GNP point
    vector<long int> gnp_contacts_point(LM_matrix_gnp.size(), 0);
    //Fill the vectors of flags
    if (!Fill_mixed_contact_flags(family, HoKo->mixed_contacts, Cutwins->boundary_flags_cnt, cnt_contacts_point, gnp_contacts_point)) {
        hout << "Error in Get_LM_matrix when calling Fill_mixed_contact_flags" << endl;
        return 0;
    }
    if (!Fill_gnp_contact_flags(HoKo->gnp_contacts, gnp_contacts_point)) {
        hout << "Error in Get_LM_matrix when calling Fill_gnp_contact_flags" << endl;
        return 0;
    }
    
    //=========================================================
    //LM matrix for CNT points
    
    //Initialize last CNT node as -1
    last_cnt_node = -1;
    
    //Check if there are any CNT clusters
    if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluster].size()) {
        
        //Scan all CNTs in the cluster
        for (int i = 0; i < (int)HoKo->clusters_cnt[n_cluster].size(); i++) {
            
            //Current CNT
            int CNT = HoKo->clusters_cnt[n_cluster][i];
            
            //Always check the first point as most likely boundary points are an endpoint of the CNT
            long int P = structure[CNT].front();
            Add_point_to_LM_matrix(P, family, Cutwins->boundary_flags_cnt, global_nodes, LM_matrix);
            elements[CNT].push_back(P);
            
            //Scan the rest of the CNT for contacts except the last point, starting on j=1 since j=0 was already taken care of
            for (int j = 1; j < (int)structure[CNT].size()-1; j++) {
                
                //Point number
                P = structure[CNT][j];
                //hout<<"P="<<P<<' ';
                
                //Check if the CNT point has contacts with other CNTs or GNPs
                if (HoKo->contacts_point[P].size() || cnt_contacts_point[P]) {
                    //If the point has contacts, then add the mapping in the LM_matix
                    Add_point_to_LM_matrix(P, family, Cutwins->boundary_flags_cnt, global_nodes, LM_matrix);
                    //Add the node to the corresponding CNT. This will be an element node
                    elements[CNT].push_back(P);
                }
            }
            
            //hout<<"elements["<<CNT<<"].size()="<<elements[CNT].size();
            //Always check the last point as most likely boundary points are an endpoint of the CNT
            P = structure[CNT].back();
            Add_point_to_LM_matrix(P, family, Cutwins->boundary_flags_cnt, global_nodes, LM_matrix);
            elements[CNT].push_back(P);
            
            //Check that the elements vector has a valid size
            if (elements[CNT].size() <= 1) {
                hout << "Error in Get_LM_matrix. The vector elements["<<CNT<<"] has size "<< elements[CNT].size();
                hout << " but it has to have at least two elements." << endl;
                hout << "\tCNT has "<<structure[CNT].size()<<" points"<<endl;
                return 0;
            }
            //hout<<endl;
        }
        
        //Update last CNT node, the last assigned CNT node is (global_nodes - 1)
        last_cnt_node = global_nodes - 1;
        
    }
    //=========================================================
    //LM matrix for GNP points
    
    //Check if there are any GNP clusters
    if (HoKo->clusters_gch.size() && HoKo->clusters_gch[n_cluster].size()) {
        
        //Scan all GNPs in the cluster
        for (int i = 0; i < (int)HoKo->clusters_gch[n_cluster].size(); i++) {
            
            //Current GNP
            int GNP = HoKo->clusters_gch[n_cluster][i];
            
            //Scan the points of the current GNP
            for (int j = 0; j < (int)structure_gnp[GNP].size(); j++) {
                
                //current point
                long int P = structure_gnp[GNP][j];
                
                //Check if the current point has any contacts
                if (gnp_contacts_point[P]) {
                    
                    //If the point has contacts, then add a new node
                    Add_point_to_LM_matrix(P, family, Cutwins->boundary_flags_gnp, global_nodes, LM_matrix_gnp);
                    
                    //Add the point to the GNP points that will be part of the triangulation
                    gnp_triangulation_points[GNP].push_back(P);
                }
                //check if the point is a boundary point and is in relevant boundary
                else if (Cutwins->boundary_flags_gnp[P].size() == 2 && Is_in_relevant_boundary(family, Cutwins->boundary_flags_gnp[P][0])) {
                    
                    //Only when a boundary GNP point is in a relevant boundary, a node is added
                    //If a boundary GNP point is at a non-relevan boundary, there is no need to add a node for that point
                    LM_matrix_gnp[P] = Get_boundary_node(Cutwins->boundary_flags_gnp[P], family);
                    
                    //Add the point to the GNP points that will be part of the triangulation
                    gnp_triangulation_points[GNP].push_back(P);
                    
                }
            }        
        }
    }
    
    return 1;
}
//Build the LM matrix and the elements matrix
//By building it at this stage I have the nodes in order
int Direct_Electrifying::Get_LM_matrix_cnts_only(const int &family, const int &n_cluster, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure, int &global_nodes, vector<int> &LM_matrix, vector<vector<long int> > &elements)
{
    //=========================================================
    //LM matrix for CNT points
    
    //Initialize last CNT node as -1
    last_cnt_node = -1;
    
    //Check if there are any CNT clusters
    if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluster].size()) {
        
        //Scan all CNTs in the cluster
        for (int i = 0; i < (int)HoKo->clusters_cnt[n_cluster].size(); i++) {
            
            //Current CNT
            int CNT = HoKo->clusters_cnt[n_cluster][i];
            
            //Always check the first point as most likely boundary points are an endpoint of the CNT
            long int P = structure[CNT].front();
            Add_point_to_LM_matrix(P, family, Cutwins->boundary_flags_cnt, global_nodes, LM_matrix);
            elements[CNT].push_back(P);
            
            //Scan the rest of the CNT for contacts except the last point, starting on j=1 since j=0 was already taken care of
            for (int j = 1; j < (int)structure[CNT].size()-1; j++) {
                
                //Point number
                P = structure[CNT][j];
                //hout<<"P="<<P<<' ';
                
                //Check if the CNT point has contacts with other CNTs
                if (HoKo->contacts_point[P].size()) {
                    //If the point has contacts, then add the mapping in the LM_matix
                    Add_point_to_LM_matrix(P, family, Cutwins->boundary_flags_cnt, global_nodes, LM_matrix);
                    //Add the node to the corresponding CNT. This will be an element node
                    elements[CNT].push_back(P);
                }
            }
            
            //hout<<"elements["<<CNT<<"].size()="<<elements[CNT].size();
            //Always check the last point as most likely boundary points are an endpoint of the CNT
            P = structure[CNT].back();
            Add_point_to_LM_matrix(P, family, Cutwins->boundary_flags_cnt, global_nodes, LM_matrix);
            elements[CNT].push_back(P);
            
            //Check that the elements vector has a valid size
            if (elements[CNT].size() <= 1) {
                hout << "Error in Get_LM_matrix. The vector elements["<<CNT<<"] has size "<< elements[CNT].size();
                hout << " but it has to have at least two elements." << endl;
                hout << "\tCNT has "<<structure[CNT].size()<<" points"<<endl;
                return 0;
            }
            //hout<<endl;
        }
        
        //Update last CNT node, the last assigned CNT node is (global_nodes - 1)
        last_cnt_node = global_nodes - 1;
        
    }
    
    return 1;
}
//This function finds all CNT and GNP points in contact with GNP points, then it fills the contact flags
int Direct_Electrifying::Fill_mixed_contact_flags(const int &family, const vector<contact_pair> &mixed_contacts, const vector<vector<short int> > &boundary_flags_cnt, vector<long int> &cnt_contacts_point, vector<long int> &gnp_contacts_point)
{
    //Scan all mixed contacts
    for (int i = 0; i < (int)mixed_contacts.size(); i++) {
        
        //If the CNT is inside the cluster, get the CNT point number of contact i
        long int P = mixed_contacts[i].point1;
        
        //If the CNT point is a boundary point in a relevant boundary, ignore the tunnel
        if (!(boundary_flags_cnt[P].size() && Is_in_relevant_boundary(family, boundary_flags_cnt[P][0]))) {
            
            //Set the value of the corresponding CNT contact point to 1
            cnt_contacts_point[P] = 1;
            
            //Get the GNP point number of contact i
            P = mixed_contacts[i].point2;
            
            //Set the value of the corresponding GNP contact point to 1
            gnp_contacts_point[P] = 1;

        }
    }
    
    return 1;
}
//This function finds all GNP points in contact with other GNP points, then it fills the contact flags
int Direct_Electrifying::Fill_gnp_contact_flags(const vector<contact_pair> &gnp_contacts, vector<long int> &gnp_contacts_point)
{
    //Scan all mixed contacts
    for (int i = 0; i < (int)gnp_contacts.size(); i++) {
        
        //Get the first GNP point number of contact i
        long int P = gnp_contacts[i].point1;
        
        //Set the value of the corresponding GNP contact point to 1
        gnp_contacts_point[P] = 1;
        
        //Get the second GNP point number of contact i
        P = gnp_contacts[i].point2;
        
        //Set the value of the corresponding GNP contact point to 1
        gnp_contacts_point[P] = 1;
        
    }
    
    return 1;
}
//
void Direct_Electrifying::Add_point_to_LM_matrix(long int P, int family, const vector<vector<short int> > &boundary_flags, int &global_nodes, vector<int> &LM_matrix)
{
    //check if the point is in a relevant boudary
    if ((boundary_flags[P].size()==2) && Is_in_relevant_boundary(family, boundary_flags[P][0])) {
        //If the point is in a relevant boundary add the reserved node number
        LM_matrix[P] = Get_boundary_node(boundary_flags[P], family);
    } else {
        //If the point is not in a boundary, then add a new node number to the point
        LM_matrix[P] = global_nodes;
        //Increase the number of nodes
        global_nodes++;
    }
    
}
//This function checks if a boundary point is in a relevant boundary depending on the family the cluster belongs to.
int Direct_Electrifying::Is_in_relevant_boundary(int family, short int boundary_node)
{
    if ( (boundary_node==0) && ((family==0)||(family==3)||(family==4)||(family==6)) )
        return 1;
    if ( (boundary_node==1) && ((family==1)||(family==3)||(family==5)||(family==6)) )
        return 1;
    if ( (boundary_node==2) && ((family==2)||(family>=4)) )
        return 1;
    else
        return 0;
}
//This function determines the boundary node number based on the family of the cluster and the location of the point
int Direct_Electrifying::Get_boundary_node(const vector<short int> &boundary_flag, const int &family)
{
    if ( (family < 0) || (family > 6)) {
        //The statement below will cause a segmentation fault error
        hout << "Error in Get_boundary_node. Family has an invalid value: " << family<<endl;
        return boundary_node_map[-1][boundary_flag[0]][boundary_flag[1]];
    } else if (family <= 2) {
        //If the family is 0, 1 or 2, then there is percolation in one direction and only two boundaries are assigned the same
        //node number. In this case we can use directly boundary_flag[1], which only has valaues 0 or 1
        return (int)boundary_flag[1];
    } else {
        //If the family is 3, 4, 5 or 6, we need to use the boundary_node_map vector
        //The boundary_node_map[0] corresponds to family three, hence I need boundary_node_map[family-3]
        return boundary_node_map[family-3][boundary_flag[0]][boundary_flag[1]];
    }
    
}
//This function creates the sparse stifness matrix that will be used to solve the sytem of equations
//The sparse version is more efficient computationally speaking
int Direct_Electrifying::Fill_sparse_stiffness_matrix(const int &R_flag, const int &nodes, const double &d_vdw, const int &n_cluster, Hoshen_Kopelman *HoKo, const vector<vector<long int> > &structure, const vector<vector<long int> > &elements, const vector<int> &LM_matrix, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<vector<long int> > &gnp_triangulation_points, const vector<int> &LM_matrix_gnp, const vector<Point_3D> &point_list_gnp, vector<GCH> &hybrid_particles, vector<vector<double> > &KEFT, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, const struct Electric_para &electric_param)
{
    //------------------------------------------------------------------------
    //Start with the 2D vectors
    
    //Variables for the 2D matrices. Initialize them with the proper number of rows (or columns)
    vector<long int> empty_long;
    vector<vector<long int> > col_ind_2d(nodes, empty_long);
    vector<double> empty_double;
    vector<vector<double> > values_2d(nodes, empty_double);
    //Set the diagonal vector to the proper size and initialize it with zeros
    diagonal.clear();
    diagonal.assign(nodes, 0);
    
    //vector<vector<long int> > contacts_point_tmp = HoKo->contacts_point;
    
    hout << "Fill_2d_matrices_cnts"<<endl;
    //Fill the 2D matrices with the contributions of the CNTs when there are CNT clusters
    if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluster].size()) {
        
        if (!Fill_2d_matrices_cnts(R_flag, elements, point_list, radii, HoKo->clusters_cnt[n_cluster], LM_matrix, electric_param, d_vdw, col_ind_2d, values_2d, diagonal, HoKo->contacts_point)){
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices" << endl;
            return 0;
        }
    }
    /*/
    if (R_flag == 0) {
        Printer *P = new Printer;
        P->Print_2d_vec(col_ind_2d, "col_ind_2d_unit_preTr.txt");
        delete P;
    } else {
        Printer *P = new Printer;
        P->Print_2d_vec(col_ind_2d, "col_ind_2d_R_preTr.txt");
        delete P;
    }
    //*/
    
    //hout << "Fill_2d_matrices_gnp"<<endl;
    //Fill the 2D matrices with the contributions of the GNPs when there are GNP clusters
    if (HoKo->clusters_gch.size() && HoKo->clusters_gch[n_cluster].size()) {
        
        //hout << "GNP resistors" << endl;
        if (!Fill_2d_matrices_gnp(R_flag, HoKo->clusters_gch[n_cluster], structure, point_list, radii, LM_matrix, point_list_gnp, gnp_triangulation_points, LM_matrix_gnp, electric_param, hybrid_particles, col_ind_2d, values_2d, diagonal)) {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_gnp" << endl;
            return 0;
        }
        
        //Fill the vector of flags for GNPs in the cluster
        vector<short int> gnps_inside_flags(hybrid_particles.size(),0);
        if (!Flags_gnps_inside(HoKo->clusters_gch[n_cluster], gnps_inside_flags)) {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Flags_gnps_inside" << endl;
            return 0;
        }
        
        /*if (R_flag == 0) {
            Printer *P = new Printer;
            P->Print_2d_vec(col_ind_2d, "col_ind_2d_unit_preGNPT.txt");
            delete P;
        } else {
            Printer *P = new Printer;
            P->Print_2d_vec(col_ind_2d, "col_ind_2d_R_preGNPT.txt");
            delete P;
        }
        
        hout << "GNP-GNP tunnels" << endl;//*/
        //Fill the 2D matrices with the contributions of the GNP-GNP tunnels
        if (!Fill_2d_matrices_gnp_gnp_tunnels(R_flag, HoKo->gnp_contacts, point_list_gnp, LM_matrix_gnp, HoKo->clusters_gch[n_cluster], gnps_inside_flags, electric_param, d_vdw, hybrid_particles, col_ind_2d, values_2d, diagonal)) {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_gnp_gnp_tunnels" << endl;
            return 0;
        }
        
        /*if (R_flag == 0) {
            Printer *P = new Printer;
            P->Print_2d_vec(col_ind_2d, "col_ind_2d_unit_preMixT.txt");
            delete P;
        } else {
            Printer *P = new Printer;
            P->Print_2d_vec(col_ind_2d, "col_ind_2d_R_preMixT.txt");
            delete P;
        }
        
        hout << "CNT-GNP tunnels" << endl;//*/
        //Fill the 2D matrices with the contributions of the CNT-GNP tunnels
        if (!Fill_2d_matrices_cnt_gnp_tunnels(R_flag, HoKo->mixed_contacts, point_list, radii, LM_matrix, point_list_gnp, LM_matrix_gnp, gnps_inside_flags, electric_param, d_vdw, hybrid_particles, col_ind_2d, values_2d, diagonal)) {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_cnt_gnp_tunnels" << endl;
            return 0;
        }
        
    }

    
    //------------------------------------------------------------------------
    //------------------------------------------------------------------------
    //Check that there are no repeated nodes on each column of the stiffness matrix
    //If there are, then check if they can be added or if there is an error
    /*/hout << "Check_repeated_col_ind_2d" << endl;
    if (R_flag == 0) {
        Printer *P = new Printer;
        P->Print_2d_vec(col_ind_2d, "col_ind_2d_unit.txt");
        delete P;
    } else {
        Printer *P = new Printer;
        P->Print_2d_vec(col_ind_2d, "col_ind_2d_R.txt");
        delete P;
    }
    if (!Check_repeated_col_ind_2d(nodes, structure, contacts_point_tmp, point_list, LM_matrix, point_list_gnp, LM_matrix_gnp, col_ind_2d, values_2d)) {
        hout << "Error in Fill_sparse_stiffness_matrix when calling Check_repeated_col_ind_2d. R_flag=" << R_flag << endl;
        return 0;
    }//*/
    //Free some memory before continuing
    //contacts_point_tmp.clear();
    
    /*/
    if (R_flag == 1) {
        Printer *P = new Printer;
        //Export stiffness matrix to check conditioning number
        Export_matlab_sparse_matrix(col_ind_2d, values_2d, diagonal, "Matrix_R.dat");
        P->Print_1d_vec(resistances, "resistances_R.txt");
        P->Print_2d_vec(values_2d, "values_2d_R.txt");
        P->Print_1d_vec(diagonal, "diagonal_R.txt");
        delete P;
    } else if (R_flag == 0) {
        Printer *P = new Printer;
        Export_matlab_sparse_matrix(col_ind_2d, values_2d, diagonal, "Matrix_unit.dat");
        P->Print_1d_vec(resistances, "resistances_unit.txt");
        P->Print_2d_vec(values_2d, "values_2d_unit.txt");
        P->Print_1d_vec(diagonal, "diagonal_unit.txt");
        delete P;
    } else {
        hout << "Error in Fill_sparse_stiffness_matrix. Invalid R_flag:" << R_flag << ". Only 0 and 1 are valid flags." << endl;
        return 0;
    }//*/
    
    //------------------------------------------------------------------------
    //------------------------------------------------------------------------
    //Convert from 2D vectors to 1D vectors
    
    //Initialize vectors
    values.clear();
    col_ind.clear();
    row_ptr.clear();
    //The first element of row_ptr is zero
    row_ptr.push_back(0);
    //empty_double.push_back(0);
    KEFT.clear();
    vector<double> zeros(reserved_nodes,0);
    KEFT.assign(nodes-reserved_nodes, zeros);
    
    hout << "From_2d_to_1d_vectors"<<endl;
    From_2d_to_1d_vectors(col_ind_2d, values_2d, KEFT, col_ind, row_ptr, values, diagonal);
    //hout << "From_2d_to_1d_vectors done"<<endl;
    
    return 1;
}
//This function adds the contributions of the CNT and junction resistors to the stiffness matrix
int Direct_Electrifying::Fill_2d_matrices_cnts(const int &R_flag, const vector<vector<long int> > &elements, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &cluster, const vector<int> &LM_matrix, const struct Electric_para &electric_param, const double &d_vdw, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point)
{
    //Variables
    int CNT;
    long int P1, P2, node1, node2;
    
    //Scan every CNT in the cluster
    for (long int i = 0; i < (long int)cluster.size(); i++) {
        
        //Current CNT in cluster
        CNT = cluster[i];
        
        //First check the special case of a CNT with three nodes, where the two endpoints are in the same boundary
        if (elements[CNT].size() == 3 && LM_matrix[elements[CNT].front()] == LM_matrix[elements[CNT].back()]) {
            
            //Find node numbers of the points
            P1 = elements[CNT][0];
            node1 = LM_matrix[P1];
            P2 = elements[CNT][1];
            node2 = LM_matrix[P2];
            long int P3 = elements[CNT][2];
            long int node3 = LM_matrix[P3];
            
            //Calculate resistance depending of the R_flag
            double Re;
            if (R_flag == 1) {
                Re = Calculate_resistance_cnt(point_list, P1, P2, radii[CNT], electric_param.resistivity_CF);
            } else if (R_flag == 0) {
                Re = 1.0;
            } else {
                hout << "Error in Fill_2d_matrices_cnts. Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
                return 0;
            }
            //resistances.push_back(Re);
            
            //Add to sparse stiffnes
            Add_elements_to_sparse_stiffness(node1, node2, Re, col_ind_2d, values_2d, diagonal);
            
            //Check if the first two nodes have any contacts and add the corresponding contributions to the
            //stiffness matrix
            //hout << "Check_for_other_elements first contact"<<endl;
            Check_for_other_elements(R_flag, point_list, radii, LM_matrix, electric_param, d_vdw, P1, node1, col_ind_2d, values_2d, diagonal, contacts_point);
            //hout << "Check_for_other_elements middle contact"<<endl;
            Check_for_other_elements(R_flag, point_list, radii, LM_matrix, electric_param, d_vdw, P2, node2, col_ind_2d, values_2d, diagonal, contacts_point);
            
            //Re-calcualte the resistance if needed
            if (R_flag == 1) {
                Re = Calculate_resistance_cnt(point_list, P2, P3, radii[CNT], electric_param.resistivity_CF);
            }
            
            //For the second segment add the contributions to the stiffness matrix without calling the function
            //since only existing values are modified, no values are added (pushed back) to the vectors
            diagonal[node2] += 1/Re;
            diagonal[node3] += 1/Re;
            
            //node3 is boundary node so node2>node3, thus -1/Re is added to values[node2].back()
            values_2d[node2].back() += -1/Re;
            
            //Check if the last node has any contacts and add the corresponding contributions to the
            //hout << "Check_for_other_elements third contact"<<endl;
            Check_for_other_elements(R_flag, point_list, radii, LM_matrix, electric_param, d_vdw, P3, node3, col_ind_2d, values_2d, diagonal, contacts_point);
            
        }
        //Else continue with the general case
        else {
            //hout << "CNT="<<CNT<<" elements[CNT].size()="<<elements[CNT].size()<<endl;
            for (long int j = 0; j < (long int)elements[CNT].size()-1; j++) {
                
                //Find node numbers of the first two points
                P1 = elements[CNT][j];
                node1 = LM_matrix[P1];
                P2 = elements[CNT][j+1];
                node2 = LM_matrix[P2];
                
                //Add the elements to the sparse vectors
                //hout << "P1=" << P1 <<" LM_matrix[P1]="<<LM_matrix[P1]<< " P2=" << P2<<" LM_matrix[P2]="<<LM_matrix[P2] << endl;
                //hout << "Add_elements_to_sparse_stiffness "<<j<<" node1="<<node1<<" node2="<<node2<<endl;
                
                //Sometimes two points at a boundary will be in contact, but since both are at a boundary
                //it makes no sense to have tunneling there
                //Thus ignore elements that have a boundary node
                //Equivalently only add elements when both nodes are different
                if (node1 != node2) {
                    double Re;
                    if (R_flag == 1) {
                        Re = Calculate_resistance_cnt(point_list, P1, P2, radii[CNT], electric_param.resistivity_CF);
                    } else if (R_flag == 0) {
                        Re = 1.0;
                    } else {
                        hout << "Error in Fill_2d_matrices_cnts. Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
                        return 0;
                    }
                    //resistances.push_back(Re);
                    Add_elements_to_sparse_stiffness(node1, node2, Re, col_ind_2d, values_2d, diagonal);
                    
                    //Check if the current node1 has any contacts and add the corresponding contributions to the
                    //stiffness matrix
                    //hout << "Check_for_other_elements nested loop "<<j<<endl;
                    Check_for_other_elements(R_flag, point_list, radii, LM_matrix, electric_param, d_vdw, P1, node1, col_ind_2d, values_2d, diagonal, contacts_point);
                }
            }
            //Check if the last node has any contacts and add the corresponding contributions to the
            //stiffness matrix
            P1 = elements[CNT].back();
            node1 = LM_matrix[P1];
            Check_for_other_elements(R_flag, point_list, radii, LM_matrix, electric_param, d_vdw, P1, node1, col_ind_2d, values_2d, diagonal, contacts_point);
            //hout << "Check_for_other_elements "<<endl;
        }

    }
    
    return 1;
}
//This function calculates the length of a CNT segment that corresponds to one element (resistor)
//Using the resisitivity as input parameter the resistance of the CNT segment is calculated
double Direct_Electrifying::Calculate_resistance_cnt(const vector<Point_3D> &point_list, const long int &P1, const long int &P2, const double &radius, const double &resistivity)
{
    //Variable to store the CNT length
    double length = 0;
    //Calculate the length of each CNT segment
    for (long int i = P1; i < P2; i++) {
        //Calculate the distance from point i to i+1 and add it to the total length
        length = length + point_list[i].distance_to(point_list[i+1]);
    }
    
    //Calculate resistance as R = rho*l/A; A = PI*r*r
    return resistivity*length/(PI*radius*radius);
}
//
void Direct_Electrifying::Add_elements_to_sparse_stiffness(const long int &node1, const long int &node2, const double &Re, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal)
{
    double Re_inv = 1/Re;
    //Add the diagonal elements of the stiffness matrix
    diagonal[node1] += Re_inv;
    diagonal[node2] += Re_inv;
    
    //Add the off diagonal elements of the stiffness matrix
    if (node1 > node2) {
        col_ind_2d[node1].push_back(node2);
        //This is the resistance between the two nodes
        values_2d[node1].push_back(-Re_inv);
    } else {
        col_ind_2d[node2].push_back(node1);
        //This is the resistance between the two nodes
        values_2d[node2].push_back(-Re_inv);
    }
    //hout << "Added ";
}
//Check if the current node1 has any contacts and add the corresponding contributions to the stiffness matrix
void Direct_Electrifying::Check_for_other_elements(const int &R_flag, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &LM_matrix, const struct Electric_para &electric_param, const double d_vdw, const long int &P1, const long int &node1, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point)
{
    if (contacts_point[P1].size()) {
        for (long int k = 0; k < (long int)contacts_point[P1].size(); k++) {
            long int P2 = contacts_point[P1][k];
            //hout <<" contact P2="<<P2;
            long int node2 = LM_matrix[P2];
            //When using the actual resistances, there are two cases when a tunnel resistor should be ignored:
            //Case 1 (Only for actual resistors):
            //This fuction applies the DEA on the backbone, thus some CNTs will be deleted and so their contacts
            //Instead of updating the vector of contacts, ignore the deleted contacts
            //The deleted contacts will have a node2 of -1
            //
            //Case 2 (for both actual and unit resistors):
            //Sometimes there is tunneling between a boundary point and a non-boundary point
            //If the CNT that contains the non-boundary point happens to be in contact with the
            //boundary (this is actually likely to happen), then there will be two resistors connecting
            //the same nodes: a CNT resistor and a tunnel resistor
            //A tunnel from a boundary to a CNT that is already in contact with the boundary does not
            //seem to be physically correct, so the tunnels between the boundary and a CNT are avoided
            //Boundary nodes have a node number below the variable reserved_nodes, so to include a tunnel
            //resistors, both of its nodes should not be reserved nodes
            //
            //Note that case 1 is solved with the condition of case 2, since -1 is less than the reserved nodes
            //so I use the same if-statement as with unit resistors
            if (node1 >= (long int)reserved_nodes && node2 >= (long int)reserved_nodes) {
                //Calculate tunnel resistance
                double Re = 0.0;
                if (R_flag == 1) {
                    
                    //Call the function for CNT-CNT tunnels
                    Re = Calculate_resistance_tunnel(radii[point_list[P1].flag], point_list[P1], radii[point_list[P2].flag], point_list[P2], electric_param, d_vdw);
                    
                } else if (R_flag == 0) {
                    
                    //Use unit resistance
                    Re = 1.0;
                    
                    //Add tunnel element (only used to find zero current in backbone)
                    vector<long int> empty;
                    elements_tunnel.push_back(empty);
                    elements_tunnel.back().push_back(P1);
                    elements_tunnel.back().push_back(P2);
                }
                //resistances.push_back(Re);
                //Add the elements to the sparse vectors
                Add_elements_to_sparse_stiffness(node1, node2, Re, col_ind_2d, values_2d, diagonal);
                //Remove the contac tha was used so it is not used again in the future
                Remove_from_vector(P1, contacts_point[P2]);
                //hout << "Removed ";
            }
        }
        //Remove all contacts of P1
        contacts_point[P1].clear();
    }
}
//This function calculates the electrical resistance due to tunneling using Hu et al. approach
//Only CNT-CNT tunnel resistances are calculated with this function
double Direct_Electrifying::Calculate_resistance_tunnel(const double &rad1, const Point_3D &P1, const double &rad2, const Point_3D &P2, const struct Electric_para &electric_param, const double &d_vdw)
{
    
    //Variable to store the tunnel area
    double A;
    
    //Calculate the separation between CNTs, which is the distance between points minus the radii
    double separation = P1.distance_to(P2);
    separation = separation  - rad1 - rad2;
    
    //The tunnel area is equal to the CNT cross-sectional area
    //The thinnest CNT is taken as it will limit the flow of electrons
    if (rad1 < rad2)
        //rad1 is the smallest of the two
        A = PI*rad1*rad1;
    else
        //rad2 is the smallest of the two
        A = PI*rad2*rad2;
    
    
    //==============================================================================
    //==============================================================================
    //Call the function that uses scaled input parameters to avoid numerical errors when reading the input file
    return Exponential_tunnel(electric_param, d_vdw, A, separation);
}
//This function calculates the electrical resistance due to tunneling using Hu et al. approach
//Only GNP-GNP and CNT-GNP tunnel resistances are calculated with this function
//Tunnel flags are as follows:
//1: GNP-GNP tunnel
//2: CNT-GNP tunnel
double Direct_Electrifying::Calculate_resistance_tunnel(const int &tunnel_flag, const GCH &hybrid1, const Point_3D &P1, const double &rad2, const Point_3D &P2, const struct Electric_para &electric_param, const double &d_vdw)
{
    
    //Variable to store the separation between particles
    double separation;
    
    //Variable to store the tunnel area
    double A;
    
    //Calculate the minimum possible distance between points and area depending on the type of tunnel
    //GNP-GNP tunnel
    if (tunnel_flag == 1) {
        
        //Calculate the separation between GNPs, which is the distance from P2 to GNP1
        separation = Calculate_distance_tunnel_point_gnp(hybrid1, P1, P2);
        
        //For a GNP-GNP tunnel, take a rectangle with lengths where rad1 is the thickness of
        //GNP1 and rad2 is the thickness of GNP2
        double rad1 = hybrid1.gnp.hei_z;
        A = rad1*rad2;
        
    }
    //CNT-GNP tunnel
    else if (tunnel_flag == 2) {
        
        //Calculate the separation between the GNP and CNT, which is the distance from P2 to the GNP1 minus the radius of CNT
        separation = Calculate_distance_tunnel_point_gnp(hybrid1, P1, P2) - rad2;
        
        //For a mixed contact, take the crossectional area of the CNT
        //Use the value of rad1
        A = PI*rad2*rad2;
        
    } else  {
        hout << "Error in Calculate_resistance_tunnel. Invalid tunnel flag: " << tunnel_flag <<". Valid flags are 1 and 2 only" << endl;
        return 0;
    }
    
    //==============================================================================
    //==============================================================================
    //Call the function that uses scaled input parameters to avoid numerical errors when reading the input file
    return Exponential_tunnel(electric_param, d_vdw, A, separation);
}
//This function calcualtes the distance of the tunnel between a point and a GNP point
//It is assumend that point1 is in the GNP and point2 is in another particle (either CNT or GNP)
double Direct_Electrifying::Calculate_distance_tunnel_point_gnp(const GCH &hybrid1, const Point_3D &P1, const Point_3D &P2)
{
    //Temporary object to take advantage of the demapping function
    Hoshen_Kopelman *HKtmp = new Hoshen_Kopelman;
    
    //Demap the points using as reference hybrid1
    Point_3D point1 = HKtmp->Demap_gnp_point(hybrid1, P1);
    Point_3D point2 = HKtmp->Demap_gnp_point(hybrid1, P2);
    
    //Delete temporary object
    delete HKtmp;
    
    //Flags to determine if the point is at a GNP surface, edge or vertex
    int fx = 0, fy = 0, fz = 0, flag = 0;
    
    //This point will be used as reference
    //If point1 is at a surface or edge, it can be used to generate new points on the surface or edge
    //If point1 is a vertex, this point will be that vertex
    Point_3D reference(0,0,0);
    
    //Detemine if point1 is at any boundary
    if ( abs(hybrid1.gnp.len_x/2 - point1.x)  < Zero) {
        fx = 1;
        reference.x = hybrid1.gnp.len_x/2;
    } else if ( abs(hybrid1.gnp.len_x/2 + point1.x) < Zero) {
        fx = -1;
        reference.x = -hybrid1.gnp.len_x/2;
    }
    if ( abs(hybrid1.gnp.wid_y/2 - point1.y)  < Zero) {
        fy = 1;
        reference.y = hybrid1.gnp.wid_y/2;
    } else if ( abs(hybrid1.gnp.wid_y/2 + point1.y) < Zero) {
        fy = -1;
        reference.y = -hybrid1.gnp.wid_y/2;
    }
    if ( abs(hybrid1.gnp.hei_z/2 - point1.z)  < Zero) {
        fz = 1;
        reference.z = hybrid1.gnp.hei_z/2;
    } else if ( abs(hybrid1.gnp.hei_z/2 + point1.z) < Zero) {
        fz = -1;
        reference.z = -hybrid1.gnp.hei_z/2;
    }
    
    //Calculate the total flag
    flag = abs(fx) + abs(fy) + abs(fz);
    
    //Detemine if point1 is at a GNP surface, edge or vertex based on the total flag
    if (flag == 1) {
        
        //point1 is in a surface, so I need three points to calculate the distance from point2 to the plane corresponding to the surface of point1
        //point1 is already in the surface, so I only need two more points
        //Procedure from: http://mathworld.wolfram.com/Point-PlaneDistance.html
        
        //The new points will start as opposite corners
        //To set them to the same plane they will have an equal coordinate
        //This coordinate will be given by the flags
        Point_3D x1(hybrid1.gnp.len_x/2,hybrid1.gnp.wid_y/2,hybrid1.gnp.hei_z/2);
        Point_3D x2(-hybrid1.gnp.len_x/2,-hybrid1.gnp.wid_y/2,-hybrid1.gnp.hei_z/2);
        
        //The non-zero flag determines the shared coordinate that puts x1 and x2 in the same plane, which is the same plane as point1
        if (fx) {
            x1.x = reference.x;
            x2.x = reference.x;
        }
        else if(fy) {
            x1.y = reference.y;
            x2.y = reference.y;
        }
        else if(fz) {
            x1.z = reference.z;
            x2.z = reference.z;
        }
        
        //Now that there are three points in a plane, calculate the distance from point2 to the plane
        
        //First, calculate a vector normal to the plane
        Point_3D x3 = point1 - x1; //For some reason the argument of cross() cannot be (point1 - x1)
        Point_3D normal = (x2 - x1).cross(x3);
        //Second, make unit normal vector
        normal = normal/sqrt(normal.dot(normal));
        //Third, calculate distance (the abs is used because this formula gives a signed distance)
        x3 = point2 - x1; //For some reason the argument of dot() cannot be (point2 - x1)
        return abs(normal.dot(x3));
        
    }
    //point1 is at an edge
    else if (flag == 2) {
        
        //point1 is on an edge, so I need two points to calculate the distance from point2 to the line corresponding to the edge of point1
        //point1 is already on the edge, so I only need one more point
        //Procedure from: http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        
        //Initialize the new point in a corner
        Point_3D x1(hybrid1.gnp.len_x/2,hybrid1.gnp.wid_y/2,hybrid1.gnp.hei_z/2);
        
        //Modify the coordinates that will set the point in the line
        if (fx) {
            x1.x = reference.x;
        }
        if(fy) {
            x1.y = reference.y;
        }
        if(fz) {
            x1.z = reference.z;
        }
        
        //Now that there are two points in a line, calculate the distance from point2 to the line
        Point_3D x3 = x1 - point2;
        //Calculate the cross product of the numerator
        x3 = (point1 - x1).cross(x3);
        //Calcualte the numerator
        double num = sqrt(x3.dot(x3));
        //Calculate the subtraction of the denominator
        x3 = point1 - x1;
        //calculate the denominator
        double den = sqrt(x3.dot(x3));
        
        //Return the absolute value of the division
        return abs(num/den);
        
    }
    //point1 is a vertex
    else if (flag == 3) {
        
        //point1 is a vertex, so calculate the distance from point1 to point2
        return point1.distance_to(point2);
    }
    
    
    hout << "Error. The GNP point is not on a surface, edge or vertex."<<endl;
    return -1;
}
double Direct_Electrifying::Exponential_tunnel(const struct Electric_para &electric_param, const double &d_vdw, const double &A, double &separation)
{
    return 1e7;
    //==============================================================================
    //==============================================================================
    //Pre-processing
    
    //Check if the tunnel distance is below the van der Waals distance
    //This happens when the penetrating model is used
    if (separation < d_vdw)
        //If the distance between the points is below the van der Waals distance, then set the separation equal to the van der Waals distance
        separation = d_vdw;
    
    //==============================================================================
    //==============================================================================
    //Input parameters are scaled to avoid numerical errors when reading, therefore the following steps are taken
    
    //Calculate quantity associted with a squared root
    double sqrt_tmp = sqrt(2*electric_param.e_mass*electric_param.lambda_barrier*electric_param.e_charge);
    
    //Calculate the exponential term
    double exp_tmp = exp(4000*PI*separation*sqrt_tmp/electric_param.h_plank);
    
    //Calculate term that multiplies the exponential
    double denominator_tmp = A*electric_param.e_charge*electric_param.e_charge*sqrt_tmp;
    double mult_tmp = 10*electric_param.h_plank*electric_param.h_plank*separation/denominator_tmp;
    
    //Calculate tunnel resistance
    return mult_tmp*exp_tmp;
}
//This function deletes an specified number from a vector
//If the number is not there, nothing happens
void Direct_Electrifying::Remove_from_vector(long int num, vector<long int> &vec)
{
    for (long int i = 0; i < (long int)vec.size(); i++)
        if (vec[i] == num) {
            vec.erase(vec.begin()+i);
            break;
        }
}
//Add contributions from the resistors formed by the triangulation on the GNP
int Direct_Electrifying::Fill_2d_matrices_gnp(const int &R_flag, const vector<int> &cluster_gnp, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &LM_matrix, const vector<Point_3D> &point_list_gnp, const vector<vector<long int> > &gnp_triangulation_points, const vector<int> &LM_matrix_gnp, const struct Electric_para &electric_param, vector<GCH> &hybrid_particles, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal)
{
    //Variables
    long int P1, P2;
    double Re;
    
    //Triangulation object
    Triangulation *delaunay = new Triangulation;
    
    //Scan every hybrid particle, perform the triangulations and add the elements to the stiffness matrix
    for (long int i = 0; i < (long int)cluster_gnp.size(); i++) {
        
        //current hybrid particle
        int hyb = cluster_gnp[i];
        
        //Pre-procesing on the points to triangulate
        //In case there is more than one point at the same boundary, only use one
        vector<int> tmp_cnts;
        if (!Same_boundary_node_preprocessing(structure, point_list, LM_matrix, LM_matrix_gnp, gnp_triangulation_points[hyb], hybrid_particles[hyb], tmp_cnts)) {
            hout << "Error in Fill_2d_matrices_gnp when calling Same_boundary_node_preprocessing" << endl;
            return 0;
        }
        
        //Perform triangulation
        if (!delaunay->Generate_3d_trangulation(last_cnt_node, point_list, structure, point_list_gnp, gnp_triangulation_points[hyb], hybrid_particles[hyb])) {
            hout << "Error in Fill_2d_matrices_gnp when calling delaunay->Generate_3d_trangulation" << endl;
            return 0;
        }
        //hout << "Triangulation " << i << ", " << gnp_triangulation_points[hyb].size() << " points in GNP " << hyb<<", " << hybrid_particles[hyb].triangulation.size() << " edges"<< endl;
        
        //Add back elements deleted in the pre-processing
        for (int j = 0; j < (int)tmp_cnts.size(); j++) {
            if (tmp_cnts[j] >= 0) {
                hybrid_particles[hyb].cnts_top.push_back(tmp_cnts[j]);
            } else {
                hybrid_particles[hyb].cnts_bottom.push_back(-tmp_cnts[j]);
            }
        }
        
        //Vector to store elements from the triangulation that need to be deleted
        vector<int> to_delete;
        
        //Add elements from the triangulation
        for (int j = 0; j < (int)hybrid_particles[hyb].triangulation.size(); j++) {
            
            //Get teh point numbers of the triangulation edge
            P1 = hybrid_particles[hyb].triangulation[j][0];
            P2 = hybrid_particles[hyb].triangulation[j][1];
            
            //Particle numbers
            int particle1, particle2;
            
            //Get the nodes
            int node1, node2;
            //Get the points
            Point_3D point1, point2;
            //Get the radii or thickness
            double rad1, rad2;
            //flags: CNT point (1) or GNP point (0)
            if (hybrid_particles[hyb].triangulation_flags[j][0]) {
                //Get the data from CNT
                node1 = LM_matrix[P1];
                point1 = point_list[P1];
                particle1 = point_list[P1].flag;
                rad1 = radii[particle1];
            } else {
                //Get the data from GNP
                node1 = LM_matrix_gnp[P1];
                point1 = point_list_gnp[P1];
                particle1 = point_list_gnp[P1].flag;
                rad1 = hybrid_particles[hyb].gnp.hei_z;
            }
            if (hybrid_particles[hyb].triangulation_flags[j][1]) {
                //Get the data from CNT
                node2 = LM_matrix[P2];
                point2 = point_list[P2];
                particle2 = point_list[P2].flag;
                rad2 = radii[particle2];
            } else {
                //Get the data from GNP
                node2 = LM_matrix_gnp[P2];
                point2 = point_list_gnp[P2];
                particle2 = point_list_gnp[P2].flag;
                rad2 = hybrid_particles[hyb].gnp.hei_z;
            }
            
            //Check if the egde should be ignored, added to all vectors (of the sparse matix) or just to the values and diagonal vectors
            
            //CASE 1: Check if the edge has to be ignored
            //Sometimes a CNT seed is also a boundary point, in that case the triangulation will result in having repeated elements in the stiffness matrix
            //This since two boundary nodes have the same node number, thus this adds a resistor connected to the same node on both ends
            //Because of this, tunneling is ignored on a bonudary point where voltage is applied
            //Thus ignore elements that have two points that are boundary nodes
            if (node1 == node2 && node1 < reserved_nodes && node2 < reserved_nodes) {
                
                //If the edge needs to be ignored, then this element needs to be deleted from the triangulation
                to_delete.push_back(j);
                hout << "TO DELETE: flag1="<<hybrid_particles[hyb].triangulation_flags[j][0]<<" P1="<<P1<<" node1="<<node1<<" flag2="<<hybrid_particles[hyb].triangulation_flags[j][1]<<" P2="<<P2<<" node2="<<node2<<endl;
            }
            else {
                
                //Whatever the case is, the resistance is needed
                //So first calculate the resistance of the edge and then check in which of the remaining cases the edge belongs to
                //Resistance of the "conduction band" in the GNP
                if (R_flag == 1) {
                    
                    //Calculate the flag for the type of triangulation resistor
                    if (!hybrid_particles[hyb].triangulation_flags[j][0] && !hybrid_particles[hyb].triangulation_flags[j][1]) {
                        
                        //flag = 1;
                        //Calculate the triangulation resistor considering a GNP-GNP edge
                        Re = Calculate_resistance_gnp(1, point1, point2, rad1, rad2, hybrid_particles[hyb], electric_param);
                    }
                    
                    else if (hybrid_particles[hyb].triangulation_flags[j][0] && hybrid_particles[hyb].triangulation_flags[j][1]) {
                        
                        //flag = 0;
                        //Calculate the triangulation resistor considering a CNT-CNT edge
                        Re = Calculate_resistance_gnp(0, point1, point2, rad1, rad2, hybrid_particles[hyb], electric_param);
                    }
                    
                    else if ( hybrid_particles[hyb].triangulation_flags[j][0] && !hybrid_particles[hyb].triangulation_flags[j][1]){
                        
                        //flag = 2;
                        //Calculate the triangulation resistor considering a GNP-CNT edge
                        //where the particle 1 is the CNT
                        Re = Calculate_resistance_gnp(2, point1, point2, rad1, rad2, hybrid_particles[hyb], electric_param);
                    }
                    
                    else {
                        
                        //flag = 2;
                        //Calculate the triangulation resistor considering a GNP-CNT edge
                        //where the particle 2 is the CNT
                        Re = Calculate_resistance_gnp(2, point2, point1, rad2, rad1, hybrid_particles[hyb], electric_param);
                    }
                
                } else if (R_flag == 0) {
                    
                    //Use unit resistance
                    Re = 1.0;
                    
                } else {
                    hout << "Error in Fill_2d_matrices_gnp. Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
                    return 0;
                }
                //resistances.push_back(Re);
                
                //CASE 2:
                //Sometimes a triangulation edge has a CNT seed connected to a boundary GNP point without any other node in between
                //In this case, the CNT resistor and the triangulation edge connect the same nodes so this results
                //in repeated elements in the stiffness matrix
                //Thus, check if this is happening
                if(Only_add_values(LM_matrix, hybrid_particles[hyb].triangulation_flags[j], node1, particle1, node2, particle2)){
                    
                    //Identify the CNT node
                    int cnt_node = 0, gnp_node = 0;
                    if (hybrid_particles[hyb].triangulation_flags[j][0]) {
                        cnt_node = node1;
                        gnp_node = node2;
                    }
                    else {
                        cnt_node = node2;
                        gnp_node = node1;
                    }
                    
                    //Add contributions to diagonal
                    diagonal[cnt_node] += 1/Re;
                    diagonal[gnp_node] += 1/Re;
                    
                    //In this situation, gnp_node > cnt_node as the CNT node is a boundary node
                    //thus -1/Re is added to values[gnp_node][k], where k is the same when
                    //col_ind_2d[gnp_node][k] == cnt_node
                    for (int k = 0; k < (int)col_ind_2d[gnp_node].size(); k++) {
                        if (col_ind_2d[gnp_node][k] == cnt_node) {
                            values_2d[gnp_node][k] += -1/Re;
                            break;
                        }
                    }
                    
                }
                //CASE 3:
                //Add the edge following the standard procedure
                else {
                    
                    //hout << "flag1="<<hybrid_particles[hyb].triangulation_flags[j][0]<<" P1="<<P1<<" node1="<<node1<<" flag2="<<hybrid_particles[hyb].triangulation_flags[j][1]<<" P2="<<P2<<" node2="<<node2<<" Re="<<Re<<endl;
                    Add_elements_to_sparse_stiffness(node1, node2, Re, col_ind_2d, values_2d, diagonal);
                }

                
            }
            
        }
        
        //Delete the elements that have boundary nodes
        for (int k = (int)to_delete.size()-1; k >= 0; k--) {
            int index = to_delete[k];
            hybrid_particles[hyb].triangulation.erase(hybrid_particles[hyb].triangulation.begin()+index);
        }
    }
    
    //hout << "Triangulation done" << endl;
    //Delete triangulation object
    delete delaunay;
    //hout << "Triangulation object deleted" << endl;
    
    return 1;
}
//
int Direct_Electrifying::Same_boundary_node_preprocessing(const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<int> &LM_matrix, const vector<int> &LM_matrix_gnp, const vector<long int> &gnp_triangulation_points, GCH &hybrid, vector<int> &tmp_cnts)
{
    
    //Scan all CNTs in the top surface, starting from the last one
    for (int j = (int)hybrid.cnts_top.size() - 1; j >= 0; j--) {
        
        //Current CNT
        int CNT = hybrid.cnts_top[j];
        
        //Current CNT seed point
        long int P = structure[CNT].front();
        
        //Get the node for the CNT seed
        int node = LM_matrix[P];
        
        //if the node is at a (relevant) boundary node, then add the CNT to the tmp_cnts vector and delete it from the hybrid
        if (node < reserved_nodes) {
            
            tmp_cnts.push_back(CNT);
            hybrid.cnts_top.erase(hybrid.cnts_top.begin()+j);
        }
    }
    
    //Scan all CNTs in the bottom surface, starting from the last one
    for (int j = (int)hybrid.cnts_bottom.size() - 1; j >= 0; j--) {
        
        //Current CNT
        int CNT = hybrid.cnts_bottom[j];
        
        //Current CNT seed point
        long int P = structure[CNT].front();
        
        //Get the node for the CNT seed
        int node = LM_matrix[P];
        
        //if the node is at a (relevant) boundary node, then add the CNT to the tmp_cnts vector and delete it from the hybrid
        if (node < reserved_nodes) {
            
            tmp_cnts.push_back(-CNT);
            hybrid.cnts_bottom.erase(hybrid.cnts_bottom.begin()+j);
        }
    }
    
    return 1;
}
//This function checks if the triangulation edge should be only added to the values and diagonal vectors
int Direct_Electrifying::Only_add_values(const vector<int> &LM_matrix, const vector<short int> &edge_flags, const long int &node1, const int &particle1, const long int &node2, const int &particle2)
{
    //Then check if only one point is CNT point and it has only two points in the elements vector
    //When the sum of the flags is 1, then there is one CNT point and one GNP point
        
        //Check which node of the triangulation is a CNT point, and if that CNT has only two nodes
        if ( (edge_flags[0] && elements[particle1].size() == 2) ) {
            
            //The CNT is particle 1, the second node of the CNT is the one that might be at the boundary and
            //causing the repeated element in the stiffness matrix (a CNT seed is always the first point of a CNT)
            long int node3 = LM_matrix[elements[particle1][1]];
            
            //Check if the node of the GNP point (node2) and the second CNT node (node3) are the same
            //In that case, ignore the edge
            if (node2 == node3) {
                return 1;
            }
        }
        else if ( (edge_flags[1] && elements[particle2].size() == 2) ) {
            
            //The CNT is particle 2, the second node of the CNT is the one that might be at the boundary and
            //causing the repeated element in the stiffness matrix (a CNT seed is always the first point of a CNT)
            long int node3 = LM_matrix[elements[particle2][1]];
            
            //Check if the node of the GNP point (node1) and the second CNT node (node3) are the same
            //In that case, ignore the edge
            if (node1 == node3) {
                return 1;
            }
        }
    
    //If none of the cases above, do not ignore the edge
    return 0;
}
//This function calculates the resistance that comes from the triangulation on the GNPs
//It uses the resistivities along the surface and along the thickness so the resistance
//can be calculated along any direction the triangulation edge may have
//Flag:
//0: CNT-CNT edge
//1: GNP-GNP edge
//2: CNT-GNP edge
double Direct_Electrifying::Calculate_resistance_gnp(const int &flag, const Point_3D &P1, const Point_3D &P2, const double &rad1, const double &rad2, const GCH &hybrid, const struct Electric_para &electric_param)
{
    //Calculate the distance between the points in contact
    double L = P1.distance_to(P2);
    
    //Unit vector in the direction from P1 to P2
    Point_3D u_direction = (P2 - P1)/L;
    
    //Multiply resistivity tensor by the u_direction vector
    Point_3D direction(u_direction.x*electric_param.resistivity_GNP_surf, u_direction.y*electric_param.resistivity_GNP_surf, u_direction.z*electric_param.resistivity_GNP_t);
    
    //Calculate cross sectional area of conduction band depending on the tye of particles
    double A = 0;
    if (flag == 0) {
        //Use CNT cross-sectional area using average of the two radii
        A = PI*(rad1 + rad2)*(rad1 + rad2)/4;
    } else if (flag == 1) {
        //Use a rectangle with sides equal to the thicknesses of GNPs
        A = rad1*rad2;
    } else if (flag == 2) {
        //Use cross-sectional area of CNT
        A = PI*rad1*rad1;
    } else {
        hout << "Error in Calculate_resistance_gnp. Invalid flag:" << flag << ". Valid flags are 0, 1 and 2 only." << endl;
        return 0;
    }
    
    
    //Calculate the angle of the actual cross sectional area
    Point_3D x_dir(1.0,0.0,0.0); //unit vector in the x direction, parallel to the surface of the conduction band
    x_dir.rotation(hybrid.rotation, hybrid.center); //rotate the unit vector so that it is the global coordinate system
    double cosA = abs(u_direction.dot(x_dir)); //cosine of the angle
    
    //The resistance is the magnitude of the vector direction multiplied by the distance between points
    // and divided by the actual cross sectional area
    return direction.distance_to(0, 0, 0)*L/(A*cosA*cosA);
}
//
int Direct_Electrifying::Flags_gnps_inside(const vector<int> &clusters_gnp, vector<short int> &gnps_inside_flags)
{
    //Scan GNPs in the cluster
    for (int i = 0; i < (int)clusters_gnp.size(); i++) {
        
        //Current GNP
        int GNP = clusters_gnp[i];
        
        //Assign the corresponding flag
        gnps_inside_flags[GNP] = 1;
    }
    
    return 1;
}
//Add contributions from the resistors formed by GNP-GNP tunnels
int Direct_Electrifying::Fill_2d_matrices_gnp_gnp_tunnels(const int &R_flag, const vector<contact_pair> &gnp_contacts, const vector<Point_3D> &point_list_gnp, const vector<int> &LM_matrix_gnp, const vector<int> &clusters_gnp, vector<short int> &gnps_inside_flags, const struct Electric_para &electric_param, const double &d_vdw, vector<GCH> &hybrid_particles, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal)
{
    //Scan all GNP contacts
    for (int i = 0; i < (int)gnp_contacts.size(); i++) {
        
        //First GNP
        int GNP1 = gnp_contacts[i].particle1;
        
        //Second GNP
        int GNP2 = gnp_contacts[i].particle2;
        
        //Check if the GNPs are inside the cluster, if not just ignore this tunnel
        if (gnps_inside_flags[GNP1] && gnps_inside_flags[GNP2]) {
            
            //If the first GNP is inside the cluster, then add the tunnel
            
            //First GNP point
            long int P1 = gnp_contacts[i].point1;
            int node1 = LM_matrix_gnp[P1];
            
            //Second GNP point
            long int P2 = gnp_contacts[i].point2;
            int node2 = LM_matrix_gnp[P2];
            
            //Debugging
            if (node1 == -1 || node2 == -1) {
                hout << "There is an invalid node number in a GNP-GNP tunnel (contact "<<i<<" of "<<gnp_contacts.size()<<"): "<< endl;
                hout <<"\tGNP1="<<GNP1<<" inside_flag="<<gnps_inside_flags[GNP1]<<" P1="<<P1<<" node1="<<node1<<endl;
                hout <<"\tGNP2="<<GNP2<<" inside_flag="<<gnps_inside_flags[GNP2]<<" P2="<<P2<<" node2="<<node2<<endl;
                return 0;
            }
            
            //Only add the tunnel when none of the nodes are boundary nodes
            //Ignore the tunnel if any of the nodes is a boundary node
            if (node1 >= reserved_nodes && node2 >= reserved_nodes) {
                
                //Calculate tunnel resistance
                double Re;
                
                //Check the flag
                if (R_flag == 1) {
                    
                    //Call the function for GNP tunnel and set the flag for GNP-GNP tunnel
                    Re = Calculate_resistance_tunnel(1, hybrid_particles[GNP1], point_list_gnp[P1], hybrid_particles[GNP2].gnp.hei_z, point_list_gnp[P2], electric_param, d_vdw);
                    
                } else if (R_flag == 0) {
                    
                    //Use unit resistance
                    Re = 1.0;
                    
                    //Add a new tunnel element (only used to find zero current in backbone)
                    vector<long int> empty;
                    elements_gnp_tunnel.push_back(empty);
                    elements_gnp_tunnel.back().push_back(P1);
                    elements_gnp_tunnel.back().push_back(P2);
                    
                } else {
                    
                    hout << "Error in Fill_2d_matrices_tunnels (GNP-GNP tunnel). Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
                    return 0;
                }
                
                //Add the elements to the sparse vectors
                //hout << "node1=" << node1 << " node2=" << node2 << " Re=" << Re <<endl;
                Add_elements_to_sparse_stiffness(node1, node2, Re, col_ind_2d, values_2d, diagonal);
                //resistances.push_back(Re);
                

            }
        }
    }
    
    return 1;
}
//Add contributions from the resistors formed by CNT-GNP tunnels
int Direct_Electrifying::Fill_2d_matrices_cnt_gnp_tunnels(const int &R_flag, const vector<contact_pair> &mixed_contacts, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &LM_matrix, const vector<Point_3D> &point_list_gnp, const vector<int> &LM_matrix_gnp, vector<short int> &gnps_inside_flags, const struct Electric_para &electric_param, const double &d_vdw, vector<GCH> &hybrid_particles, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal)
{
    //Scan all mixed contacts
    for (int i = 0; i < (int)mixed_contacts.size(); i++) {
        
        //Get GNP number
        int GNP = mixed_contacts[i].particle2;
        
        //Check if the GNP is inside the cluster
        if (gnps_inside_flags[GNP]) {
            
            //If the GNP is inside the cluster, then add the contribution of the tunnel
            
            //GNP point
            long int P2 = mixed_contacts[i].point2;
            int node2 = LM_matrix_gnp[P2];
            
            
            //Get CNT number
            int CNT = mixed_contacts[i].particle1;
            //CNT point
            long int P1 = mixed_contacts[i].point1;
            int node1 = LM_matrix[P1];
            
            //Debugging
            //hout <<i<<" CNT="<<CNT<<" P1="<<P1<<" node1="<<node1<<endl;
            //hout <<" GNP2="<<GNP<<" flag="<<gnps_inside_flags[GNP]<<" P2="<<P2<<" node2="<<node2<<endl;
            
            //Only add the tunnel if both nodes are different and if both nodes are not boundary nodes
            if (node1 != node2 && node1 >= reserved_nodes && node2 >= reserved_nodes) {
                //Calculate tunnel resistance
                double Re;
                if (R_flag == 1) {
                    
                    //Call the function for GNP tunnel and set the flag for CNT-GNP tunnel
                    Re = Calculate_resistance_tunnel(2, hybrid_particles[GNP], point_list_gnp[P2], radii[CNT], point_list[P1], electric_param, d_vdw);
                    
                } else if (R_flag == 0) {
                    
                    //Use unit ressitance
                    Re = 1.0;
                    
                    //Add a new tunnel element (only used to find zero current in backbone)
                    vector<long int> empty;
                    elements_mixed_tunnel.push_back(empty);
                    elements_mixed_tunnel.back().push_back(P1);
                    elements_mixed_tunnel.back().push_back(P2);
                    
                } else {
                    hout << "Error in Fill_2d_matrices_tunnels (CNT-GNP tunnel). Invalid resistor flag:" << R_flag << ". Valid flags are 0 and 1 only." << endl;
                    return 0;
                }
                
                //Add the elements to the sparse vectors
                Add_elements_to_sparse_stiffness(node1, node2, Re, col_ind_2d, values_2d, diagonal);
                
            }
            
        }

    }
    
    return 1;
}
//Test function to find repeated elements in the vector col_ind_2d
int Direct_Electrifying::Check_repeated_col_ind_2d(const int &nodes, const vector<vector<long int> > &structure, vector<vector<long int> > &contacts_point, const vector<Point_3D> &point_list, const vector<int> &LM_matrix, const vector<Point_3D> &point_list_gnp, const vector<int> &LM_matrix_gnp, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d)
{
    //Create inverse mapping for the LM_matrix, i.e., from node number to point number
    //Initialize the LM_inverse with the number of nodes
    vector<long int> LM_inverse(nodes, -1);
    
    //Vector of flags to choose CNT points or GNP points
    //1: CNT node
    //0: GNP node
    vector<short int> point_flag(nodes, 1);
    
    //Vectors of boundary points
    vector<long int> empty_long;
    vector<vector<long int> > gnp_boundary_points(reserved_nodes, empty_long), cnt_boundary_points(reserved_nodes, empty_long);
    
    //Fill the LM_inverse matrix
    for (long int i = 0; i < (long int)LM_matrix.size(); i++) {
        int current_node = LM_matrix[i];
        if (current_node != -1) {
            //Check if the node is a boundary node
            if (current_node < reserved_nodes) {
                cnt_boundary_points[current_node].push_back(i);
                LM_inverse[current_node] = -2;
            }
            else {
                LM_inverse[current_node] = i;
            }
        }
    }
    for (long int i = 0; i < (long int)LM_matrix_gnp.size(); i++) {
        int current_node = LM_matrix_gnp[i];
        if (current_node != -1) {
            //Check if the node is a boundary node
            if (current_node < reserved_nodes) {
                gnp_boundary_points[current_node].push_back(i);
                LM_inverse[current_node] = -2;
            }
            else {
                LM_inverse[current_node] = i;
            }
            point_flag[current_node] = 0;
        }
    }
    
    //Flag to print boundary points
    //This flag will be activated whenever there is repeatition involving boundary nodes
    int print_flag = 0;
    
    //Loop over the rows of the 2D sparse matrix
    for (int i = 0; i < (int)col_ind_2d.size(); i++) {
        
        //Find repeated elements in the vector col_ind_2d[i]
        vector<long int> repeated_elements;
        vector<int> indices;
        Find_repeated_elements(col_ind_2d[i], repeated_elements, indices);
        
        if (repeated_elements.size()) {
            
            //Check if node is a boundary node
            if (i < reserved_nodes) {
                
                //Flag to print once "boundary node i"
                int print_once_flag = 1;
                
                for (int j = 0; j < (int)repeated_elements.size(); j++) {
                    
                    //If there are repeated elements connecting two boundaries, these repeated elements are not important
                    //When transfroming 2D vectors to 1D vectors the sub-matrix that contains resistors connecting boundaries is removed
                    //This, because the sub-matrix used to solve the system of equations does not consider these elements
                    //Thus, only print the information of repeated elements when they are not boundary nodes
                    if (repeated_elements[j] >= reserved_nodes) {
                        
                        if (print_once_flag) {
                            
                            hout << "boundary node " << i <<endl;
                            hout << "repeated elements: "<<repeated_elements.size()<<endl;
                            
                            //Activate print flag
                            print_flag = 1;
                            
                            //Set the flag to 0
                            print_once_flag = 0;
                        }
                        
                        //store the repeated node in a variable to ease reading
                        long int node2 = repeated_elements[j];
                        
                        //Get the point number of the repeated node
                        long int P2 = LM_inverse[node2];
                        
                        int particle2;
                        string type2;
                        
                        //Get particle number and type of repeated point
                        if (point_flag[node2]) {
                            particle2 = point_list[P2].flag;
                            type2 = "CNT";
                        } else {
                            particle2 = point_list_gnp[P2].flag;
                            type2 = "GNP";
                        }
                        
                        //Print info of repeated point
                        hout << j << ") node " << node2 << " P2="<<P2<<' '<<type2<<"="<<particle2<<endl;
                        if (type2 == "GNP") {
                            hout <<"("<<point_list_gnp[P2].x<<", "<<point_list_gnp[P2].y<<", "<<point_list_gnp[P2].z<<')'<<endl;
                        } else {
                            hout <<"("<<point_list[P2].x<<", "<<point_list[P2].y<<", "<<point_list[P2].z<<')'<<endl;
                            hout <<" CNT[0]="<<structure[particle2].front()<<" n_points="<<structure[particle2].size();
                            hout <<" CNT[last]="<<structure[particle2].back()<< endl;
                            hout <<" nodes in CNT: ";
                            for (int nn = 0; nn < (int)elements[particle2].size(); nn++) {
                                hout << LM_matrix[elements[particle2][nn]] << ' ';
                            }
                            hout<<endl;
                        }
                    }
                }
                
                //Add a space to separate repeated contacts
                if (print_flag)
                    hout << endl;
                
            }
            //Not a boundary node
            else {
                
                //Get the point number of current node
                long int P1 = LM_inverse[i];
                
                int particle1;
                string type1;
                
                //Get particle number and type of current point
                if (point_flag[i]) {
                    particle1 = point_list[P1].flag;
                    type1 = "CNT";
                } else {
                    particle1 = point_list_gnp[P1].flag;
                    type1 = "GNP";
                }
                
                
                hout << "node " << i << " P1="<<P1<<' '<<type1<<"="<<particle1<<endl;
                if (type1 == "GNP") {
                    hout <<"("<<point_list_gnp[P1].x<<", "<<point_list_gnp[P1].y<<", "<<point_list_gnp[P1].z<<')'<<endl;
                } else {
                    hout <<"("<<point_list[P1].x<<", "<<point_list[P1].y<<", "<<point_list[P1].z<<')'<<endl;
                    hout <<" n_points="<<structure[particle1].size() <<endl;
                    hout <<" CNT[0]="<<structure[particle1].front()<<" CNT[last]="<<structure[particle1].back()<< endl;
                    hout <<" nodes in CNT: ";
                    for (int nn = 0; nn < (int)elements[particle1].size(); nn++) {
                        hout << LM_matrix[elements[particle1][nn]] << ' ';
                    }
                    hout<<endl;
                }
                
                hout << "repeated elements: "<<repeated_elements.size()<<endl;
                for (int j = 0; j < (int)repeated_elements.size(); j++) {
                    if (repeated_elements[j] < reserved_nodes) {
                        
                        //Activate print flag
                        print_flag = 1;
                        
                        //boundary node has repeated resistors with another boundary node
                        hout << j << ") boundary node " << repeated_elements[j] <<endl;
                    }
                    else {
                        
                        //store the repeated node in a variable to ease reading
                        long int node2 = repeated_elements[j];
                        
                        //Get the point number of the repeated node
                        long int P2 = LM_inverse[node2];
                        
                        int particle2;
                        string type2;
                        
                        //Get particle number and type of repeated point
                        if (point_flag[node2]) {
                            particle2 = point_list[P2].flag;
                            type2 = "CNT";
                        } else {
                            particle2 = point_list_gnp[P2].flag;
                            type2 = "GNP";
                        }
                        
                        //Print info of repeated point
                        hout << j << ") node " << node2 << " P2="<<P2<<' '<<type2<<"="<<particle2<<endl;
                        if (type2 == "GNP") {
                            hout <<"("<<point_list_gnp[P2].x<<", "<<point_list_gnp[P2].y<<", "<<point_list_gnp[P2].z<<')'<<endl;
                        } else {
                            hout <<"("<<point_list[P2].x<<", "<<point_list[P2].y<<", "<<point_list[P2].z<<')'<<endl;
                            hout <<" CNT[0]="<<structure[particle2].front()<<" n_points="<<structure[particle2].size();
                            hout <<" CNT[last]="<<structure[particle2].back()<< endl;
                            hout <<" nodes in CNT: ";
                            for (int nn = 0; nn < (int)elements[particle2].size(); nn++) {
                                hout << LM_matrix[elements[particle2][nn]] << ' ';
                            }
                            hout<<endl;
                        }
                    }
                }
                
                //Add a space to separate repeated contacts
                hout << endl;
                
            }
            
        }
        
    }
    
    //If the print flag was activated, print all points at the boundaries for both CNTs and GNPs
    if (print_flag) {
        
        //For reference, print out the number of CNTs
        hout << "Total CNTs: " << structure.size() << endl;
        
        hout << "GNP boundary points:" << endl;
        for (int i = 0; i<(int)gnp_boundary_points.size(); i++) {
            hout << "boundary " << i << ": " << endl;
            for (int j = 0; j < (int)gnp_boundary_points[i].size(); j++) {
                hout << gnp_boundary_points[i][j] << " ";
            }
            hout<<endl;
        }
        hout << "CNT boundary points:" << endl;
        for (int i = 0; i<(int)cnt_boundary_points.size(); i++) {
            hout << "boundary " << i << ": " << endl;
            for (int j = 0; j < (int)cnt_boundary_points[i].size(); j++) {
                hout << cnt_boundary_points[i][j] << " ";
            }
            hout<<endl;
        }
    }
    return 1;
}
//This function exports the 2D sparse vector to a format that matlab can use
int Direct_Electrifying::Export_matlab_sparse_matrix(const vector<vector<long int> > &col_ind_2d, const vector<vector<double> > &values_2d, const vector<double> &diagonal, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    //Skip the top rows of the matrix that correspond to the reserved nodes
    for (long int i = (long int)reserved_nodes; i < (long int)col_ind_2d.size(); i++) {
        //Matlab indices start in 1
        long int row = i-(long int)reserved_nodes+1;
        otec << row << '\t' << row << '\t' << diagonal[i] << endl;
        for (long int j = 0; j < (long int)col_ind_2d[i].size(); j++) {
            if (col_ind_2d[i][j] >= reserved_nodes) {
                //Matlab indices start in 1
                otec << row << '\t' << col_ind_2d[i][j]-reserved_nodes+1 << '\t' << values_2d[i][j] << endl;
                otec << col_ind_2d[i][j]-reserved_nodes+1 << '\t' << row << '\t' << values_2d[i][j] << endl;
            }
        }
    }
    //Close file
    otec.close();
    return 1;
}
//This function transforms the 2D vectors that contain the stiffness matrix into 1D vectors so they can be in the SSS format and
//make the matrix-vector multiplications faster
//For reference, the 2D stiffness matrix has this form:
//
//   | KE    KEF |
//   | KEFT  KF  |
//
//For the CG algorithm we need to extract KEFT and KF. We actually already have the lower left corner of the stiffness matrix
//
void Direct_Electrifying::From_2d_to_1d_vectors(vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<vector<double> > &KEFT, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal)
{
    //hout << "Fill 1D vectors" <<endl;
    //Skip the top rows of the matrix that correspond to the reserved nodes
    for (long int i = (long int)reserved_nodes; i < (long int)col_ind_2d.size(); i++) {
        //hout << "for(i)=" <<i << " col_ind_2d[i].size()="<<col_ind_2d[i].size()<< endl;
        for (long int j = 0; j < (long int)col_ind_2d[i].size(); j++) {
            //hout << "for(j)="<<j<<' ';
            if (col_ind_2d[i][j] >= reserved_nodes) {
                //hout << "if1 col_ind_2d["<<i<<"]["<<j<<"]="<<col_ind_2d[i][j]<<' ';
                //When the colum index is greater or equal to reserved_nodes, then it belogs to the lower-righ
                values.push_back(values_2d[i][j]);
                //The column numbers that I need are the numbers in the full matrix-reserved_nodes
                col_ind.push_back(col_ind_2d[i][j]-reserved_nodes);
            } else  {
                //hout << "if2 col_ind_2d["<<i<<"]["<<j<<"]="<<col_ind_2d[i][j]<<' ';
                //When the column index is less than reserved_nodes, I need to save the value of the resistance on the vector KEFT
                //so that it is used for the CG algorithm.
                //The column index is (of course) the column index of KEFT as 0<=col_ind_2d[i][j]<reserved_nodes in this else-statement
                //The row is the index iterator i-reserved_nodes
                long int col = col_ind_2d[i][j];
                KEFT[i-reserved_nodes][col] = values_2d[i][j];
            }
        }
        //hout << "row_ptr" << endl;
        row_ptr.push_back(values.size());
        
        //Free memory
        col_ind_2d[i].clear();
        values_2d[i].clear();
    }
    //Remove the elements of the diagonal that are not used
    //hout << "Remove the elements of the diagonal that are not used"<<endl;
    for (int i = 0; i < reserved_nodes; i++) {
        diagonal.erase(diagonal.begin());
    }
    
}
//This function solves the equation of the electric circuit as done in the Direct Electrifing Algorithm (DEA)
int Direct_Electrifying::Solve_DEA_equations_CG_SSS(const int &R_flag, long int nodes, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, const struct Electric_para &electric_param, vector<vector<double> > &KEFT)
{
    
    //=========================================
    // Set up variables for the Conjugate Gradient Algorithm
    
    //Voltage applied to the sample
    vector<double> VEF;
    //Initializa the prescribed voltage boundary conditions
    //The magnitude of the voltage depends on the R_flag
    if (R_flag == 1) {
        //If R_falg is 1, then use real voltage (from the input parameters)
        Get_voltage_vector(electric_param.applied_voltage, VEF);
    } else if (!R_flag) {
        //If R_falg is 0, then use the number of nodes as the magnitude for the voltage
        Get_voltage_vector((double)nodes, VEF);
    } else {
        hout << "Error in Solve_DEA_equations_CG_SSS. The R_flag has an invalid value: " << R_flag << endl;
        return 0;
    }
    hout << setwp(1,20) << "Maximum and minimum voltages = " << VEF.front() << ", " << VEF.back() << endl;
    
    //P is the search direction
    //R is the residual vector
    vector<double> P, R;
    //Initialize P and R
    P.assign(nodes-reserved_nodes,0);
    R.assign(nodes-reserved_nodes,0);
    
    //The residual vector is initialized with R = b - Ax0.
    //x0 is the initial guess. If we use x0 = 0 as initial guess then R = b
    //From the matrix equations b = - KEFT*VEF
    for (int i = 0; i < (int)KEFT.size(); i++) {
        
        for (int j = 1; j < (int)KEFT[i].size(); j++) {
            //the first element of VEF (i.e. VEF[0]) is zero, so it can be skipped to save computational time
            R[i] = R[i] - (KEFT[i][j]*VEF[j]);
        }
        
        //The search direction of the CG is initialized with the initial value of the residual
        P[i] = R[i];
    }
    
    //Clear vector to free memory
    KEFT.clear();
    
    //=========================================
    // Conjugate Gradient Algorithm
    hout << "Start CG" <<endl;
    Conjugate_gradient(nodes, col_ind, row_ptr, values, diagonal, R, P);
    
    //The known boundary conditions are added at the beginning of the solution
    for (int i = reserved_nodes-1; i >=0; i--) {
        voltages.insert(voltages.begin(), VEF[i]);
    }
    
    /*/Print the voltages
    Printer *Pr = new Printer;
    if (R_flag) {
        Pr->Print_1d_vec(voltages, "voltages_R.txt");
    } else {
        Pr->Print_1d_vec(voltages, "voltages_unit.txt");
    }
    delete Pr;//*/
    
    return 1;
}
//This function creates a voltage vector depending on the number of prescribed boundary conditios
void Direct_Electrifying::Get_voltage_vector(const double &nodes, vector<double> &voltages)
{
    //Clear the vector of voltages
    voltages.clear();
    for (int i = 0; i < reserved_nodes; i++) {
        voltages.push_back( nodes*((double)i) );
    }
}
//This function solves the system of equations using the CG gradient
//P is the search direction
//R is the residual vector
//The residual vector is initialized with R = b - Ax0.
//x0 is the initial guess. If we use x0 = 0 as initial guess then R = b
//Also, P = R = b
void Direct_Electrifying::Conjugate_gradient(long int nodes, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, vector<double> &R, vector<double> &P)
{
    //Preconditioner
    vector<double> M_inv(diagonal.size(),0.0);
    Jacobi_preconditioner(diagonal, M_inv);
    
    //Apply preconditoner
    vector<double> Y(R.size(),0.0);
    Apply_preconditioner(M_inv, R, P, Y);
    
    //Variables of the algorithm
    vector<double> AP; //Variable for the multiplication A*P
    AP.assign(nodes-reserved_nodes,0);
    double alpha, beta, rr0, rr;
    voltages.clear();
    voltages.assign(nodes-reserved_nodes, 0);
    
    //Maximum number of iterations for the CG
    long int max_iter = 10*nodes;
    //Iteration variable
    long int k;
    //Variable to check the status of the CG
    int test = 50000, test_inc = 50000;
    
    //Initial residual
    double R0 = 1.0E-10*sqrt(V_dot_v(R, R));
    //double R0 = Zero*sqrt(V_dot_v(R, R));
    //double R0 = 1.0E-10;
    //hout << "R0 = " << R0 << endl;
    
    //Predonditioned
    for (k = 1; k <= max_iter; k++) {
        //Calculate Ap
        spM_V_SSS(P, row_ptr, col_ind, diagonal, values, AP);
        //Calculate norm or residual of step k-1. Will be used later as convergence criteria and to calculate beta
        rr0 = V_dot_v(Y, R);
        //Step length
        alpha = rr0/(V_dot_v(P, AP));
        //Approximate solution
        //X = X + P*alpha;
        V_plus_aW(P, alpha, voltages);
        //Residual
        //R = R - AP*alpha;
        V_plus_aW(AP, -alpha, R);
        //Calculate norm or residual of step k. Used as convergence criteria and to calculate beta
        rr = V_dot_v(R, R);
        //Status update: print every test iterations
        if ( k == test){
            hout << "CG iteration " << k << endl;
            test = test + test_inc;
        }
        //Convergence criteria
        if (sqrt(rr) <= R0)
            break;
        //Update Y
        //Y = M_inv*R
        Componentwise_multiply(M_inv, R, Y);
        //Improvement of step
        beta = V_dot_v(Y, R)/rr0;
        //Search direction
        //P = Y + P*beta;
        W_plus_aV(Y, beta, P);
    }
    
    if (k >= max_iter)
        hout << "CG reached maximum number of iterations" << endl;
    hout << "CG iterations: " << k << endl;
    //hout << "RR = " << sqrt(rr) << endl;
}
//Calculate the Jacobi preconditioner M_inv
//M_inv[i] is simply 1/diagonal[i]
void Direct_Electrifying::Jacobi_preconditioner(const vector<double> &diagonal, vector<double> &M_inv)
{
    for (int i = 0; i < (int)diagonal.size(); i++) {
        //M_inv[i] is simply 1/diagonal[i]
        M_inv[i] = 1/diagonal[i];
    }
}
void Direct_Electrifying::Apply_preconditioner(const vector<double> &M_inv, const vector<double> &R, vector<double> &P, vector<double> &Y)
{
    //P = Y = M_inv*R
    for (int i = 0; i < (int)R.size(); i++) {
        P[i] = M_inv[i]*R[i];
        Y[i] = P[i];
    }
}
//This function performs a sparse-matrix-vector multiplication using the Symmetric Sparse Skyline (SSS) format
void Direct_Electrifying::spM_V_SSS(const vector<double> &V, const vector<long int> &rowptr, const vector<long int> &colind, const vector<double> &diagonal, const vector<double> &values, vector<double> &R)
{
    //Size of the system
    long int N = V.size();
    //Initialize result vector
    R.clear();
    R.assign(N,1);
    //SSS
    long int c;
    for (long int r = 0; r < N; r++) {
        R[r] = diagonal[r]*V[r];
        for (long int j = rowptr[r]; j < rowptr[r+1]; j++) {
            c = colind[j];
            R[r] = R[r] + values[j]*V[c];
            R[c] = R[c] + values[j]*V[r];
        }
    }
}
//This function calculates the dot product of two vectors
double Direct_Electrifying::V_dot_v(const vector<double> &A, const vector<double> &B)
{
    if (A.size() != B.size()){
        hout << "Vectors must have the same length. A.size()="<< A.size() << " B.size()=" << B.size() << endl;
        double tmp = 0;
        return 1/tmp;
    }
    
    //This variable will store the result of the dot product
    double dot_product = 0;
    
    for (int i = 0; i < (int)A.size(); i++) {
        dot_product = dot_product + A[i]*B[i];
    }
    
    return dot_product;
}
//This function solves the following operation: V = V + a*W
//where V and W are vectors and a is a scalar
void Direct_Electrifying::V_plus_aW(const vector<double> &W, const double &a, vector<double> &V)
{
    for (int i = 0; i < (int)V.size(); i++) {
        V[i] = V[i] + a*W[i];
    }
}
//This function solves the following operation: V = W + a*V
//where V and W are vectors and a is a scalar
void Direct_Electrifying::W_plus_aV(const vector<double> &W, const double &a, vector<double> &V)
{
    for (int i = 0; i < (int)V.size(); i++) {
        V[i] = W[i] + a*V[i];
    }
}
//Find the repeated elements in a vector
void Direct_Electrifying::Find_repeated_elements(const vector<long int> &vector_in, vector<long int> &elements, vector<int> &indices)
{
    for (int i = 0; i < (int)vector_in.size()-1; i++) {
        for (int j = i+1; j < (int)vector_in.size(); j++) {
            if (vector_in[i] == vector_in[j]) {
                elements.push_back(vector_in[j]);
                indices.push_back(i);
                indices.push_back(j);
            }
        }
    }
}
//This function multiplies two vectors componentwise
void Direct_Electrifying::Componentwise_multiply(const vector<double> &vector_in1, const vector<double> &vector_in2, vector<double> &vector_out)
{
    for (int i = 0; i < (int)vector_in1.size(); i++) {
        vector_out[i] = vector_in1[i]*vector_in2[i];
    }
}