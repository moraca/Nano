//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Direct_Eletrifying.cpp
//OBJECTIVE:	The direct eletrifying algorithm (C.Y. Li and T.W. Chou, Int. J. Mod Phys C, 20, 2009, 423-33.)
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Direct_Electrifying.h"

/*
 Input:
    vector<vector<long int> > structure
        Vector with the structure
    vector<vector<long int> > contacts_point
        Vector of point to point contacts.
    vector<vector<short int> > boundary_flags
        Vector of flags that indicates which CNT endpoints are on a boundary and in which boundary
    vector<int> cluster
        List of CNTs that belong to the percolated cluster were
    vector<double> radii
        List of radii. Using this vector allows for the code to be able to work with CNTs of different radii
    int family
        Percoalted family to which the cluster belongs to. This will help to set up the boundary conditions
    struct Electric_para electric_param
        structure that contains resistivity of CNTs and voltage applied to the sample
 
 Output (These two are class variables):
    vector<double> voltages
        Vector with the value of the voltages at every node (contact point)
    double resistances
        Electrical resistance of the network 
    vector<vector<long int> > elements;
        Vector with the elements for the DEA. This vector has the same size as structure. Each element is one segment of a CNT that is a resistor. Each elements[i] is a CNT and each elements[i][j] is a point. Every elements[i][0] is the first point of a CNT and every elements[i].back is the last point of that CNT. 
 */

//Calculate the voltage values at contact points and endpoints
int Direct_Electrifying::Calculate_voltage_field(const int &family, const int &R_flag, const vector<vector<long int> > &structure, vector<vector<long int> > &contacts_point, const vector<vector<short int> > &boundary_flags, const vector<Point_3D> &point_list, const vector<int> &cluster, const vector<int> &cluster_gch, const vector<double> &radii, const struct Electric_para &electric_param, const struct Cutoff_dist &cutoffs, vector<GCH> &hybrid_particles)
{
    //First we need to prepare the matrices that the direct electrifying needs
    //The first matrix will be the local mapping (LM) matrix. This matrix maps from point number in the structure
    //to node number for the solver
    //The number of reserved nodes is calculated. These is the number of boundaries with prescribed voltage
    reserved_nodes = Get_global_nodes(family);
    //This variable is used to assing node numbers and after calling Get_LM_matrix, it will contain the number of nodes in the network
    int global_nodes = reserved_nodes;
    //Initialize the size of the LM matrix to be equal to the number of points
    //In the case of hybrid particles, there are two more nodes per hybrid particle in the cluster--each new node corresponds to a center of a GNP surface
    LM_matrix.assign(contacts_point.size(), -1);
    
    //hout << "LM_matrix.size()="<<LM_matrix.size()<<" cluster.size()="<<cluster.size()<<endl;
    //If there is a cluster wiht one CNT, then there is no need to do the DEA
    //This heppens when the CNT length is larger than one or more dimensions of the observation window
    if (cluster.size() == 1) {
        return 1;
    }
    
    //Initialize the vector boundary_node_map
    Initialize_boundary_node_map();
    //Initialize the size of the elements matrix to be equal to the number of CNTs
    vector<long int> empty;
    elements.assign(structure.size(), empty);
    if(!Get_LM_matrix(structure, contacts_point, boundary_flags,cluster, family , global_nodes, LM_matrix, elements)){
        hout << "Error in Calculate_voltage_field when calling Get_LM_matrix" << endl;
        return 0;
    }
    
    //hout << "Fill_sparse_stiffness_matrix"<<endl;
    //Variables for using the SSS for the sparse matrix
    vector<long int> col_ind, row_ptr;
    vector<double> values, diagonal, R;
    vector<vector<double> > KEFT;
    //With the LM matrix, now fill the sparse stiffness matrix
    if (!Fill_sparse_stiffness_matrix(R_flag, structure, contacts_point, point_list, radii, cluster_gch, cluster, global_nodes, cutoffs.van_der_Waals_dist, LM_matrix, elements, hybrid_particles, KEFT, col_ind, row_ptr, values, diagonal, electric_param)){
        hout << "Error in Calculate_voltage_field when calling Fill_sparse_stiffness_matrix" << endl;
        return 0;
    }
    
    //hout << "Solve_DEA_equations_CG_SSS"<<endl;
    //This is where the actual direct electrifying algorithm (DEA) takes place
    if (!Solve_DEA_equations_CG_SSS(R_flag, global_nodes, col_ind, row_ptr, values, diagonal, electric_param, KEFT)) {
        hout << "Error in Calculate_voltage_field when calling Solve_DEA_equations_CG_SSS" << endl;
        return 0;
    }
    
    return 1;
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
int Direct_Electrifying::Get_LM_matrix(const vector<vector<long int> > &structure, const vector<vector<long int> > &contacts_point, const vector<vector<short int> > &boundary_flags, const vector<int> &cluster, const int &family, int &global_nodes, vector<int> &LM_matrix, vector<vector<long int> > &elements)
{
    //Variables
    int CNT;
    long int P;
    
    for (int i = 0; i < (int)cluster.size(); i++) {
        CNT = cluster[i];
        
        //Always check the first point as most likely boundary points are an endpoint of the CNT
        P = structure[CNT].front();
        Add_point_to_LM_matrix(P, family, boundary_flags, global_nodes, LM_matrix);
        elements[CNT].push_back(P);
        
        //Scan the rest of the CNT for contacts except the last point, starting on j=1 since j=0 was already taken care of
        for (int j = 1; j < (int)structure[CNT].size()-1; j++) {
            //Point number
            P = structure[CNT][j];
            //hout<<"P="<<P<<' ';
            if (contacts_point[P].size()) {
                //If the point has contacts, then add the elements to the LM_matix
                Add_point_to_LM_matrix(P, family, boundary_flags, global_nodes, LM_matrix);
                //Add the node to the corresponding CNT. This will be an element node
                elements[CNT].push_back(P);
            }
        }
        
        //hout<<"elements["<<CNT<<"].size()="<<elements[CNT].size();
        //Always check the last point as most likely boundary points are an endpoint of the CNT
        P = structure[CNT].back();
        Add_point_to_LM_matrix(P, family, boundary_flags, global_nodes, LM_matrix);
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
//This function determines the boundary node number based on the family of the cluster and the location of the point
int Direct_Electrifying::Get_boundary_node(const vector<short int> &boundary_flag, const int &family)
{
    if ( (family < 0) || (family > 6)) {
        //The statement below will cause a segmentation fault error
        hout << "Error in Get_boundary_node. Family has an invalid value: " << family<<endl;
        return boundary_node_map[-1][boundary_flag[0]][boundary_flag[1]];
    } else if (family <= 2) {
        //If the family is 0, 1 or 2, then there is percolation in une direction and only two boundaries are assigned the same
        //node number. In this case we can use directly boundary_flag[1], which only has valaues 0 or 1
        return (int)boundary_flag[1];
    } else {
        //If the family is 3, 4, 5 or 6, we need to use the boundary_node_map vector
        //The boundary_node_map[0] corresponds to family three, hence I need boundary_node_map[family-3]
        return boundary_node_map[family-3][boundary_flag[0]][boundary_flag[1]];
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
//This function creates the sparse stifness matrix that will be used to solve the sytem of equations
//The sparse version is more efficient computationally speaking
int Direct_Electrifying::Fill_sparse_stiffness_matrix(const int &R_flag, const vector<vector<long int> > &structure, vector<vector<long int> > &contacts_point, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &cluster_gch, const vector<int> &cluster, const int &nodes, const double &d_vdw, const vector<int> &LM_matrix, const vector<vector<long int> > &elements, vector<GCH> &hybrid_particles, vector<vector<double> > &KEFT, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, const struct Electric_para &electric_param)
{
    //------------------------------------------------------------------------
    //Start with the 2D vectors
    
    //Variables for the 2D matrices. Initialize them with the proper number of rows (or columns)
    vector<vector<long int> > col_ind_2d;
    vector<long int> empty_long;
    col_ind_2d.assign(nodes, empty_long);
    vector<vector<double> > values_2d;
    vector<double> empty_double;
    values_2d.assign(nodes, empty_double);
    //Set the diagonal vector to the proper size and initialize it with zeros
    diagonal.clear();
    diagonal.assign(nodes, 0);
    
    vector<vector<long int> > contacts_point_tmp = contacts_point;
    
    //Check the R_flag
    if (R_flag == 1) {
        //If R_flag == 1, then calculate actual resistances
        
        //Printer *P = new Printer;
        //hout << "Fill_2d_matrices"<<endl;
        //Fill the 2D matrices with the contributions of the CNTs
        if (!Fill_2d_matrices(elements, point_list, radii, cluster, LM_matrix, electric_param, d_vdw, col_ind_2d, values_2d, diagonal, contacts_point)){
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices" << endl;
            return 0;
        }
        /*/------------------------------------------------------------------------
        //------------------------------------------------------------------------
        //CHECK
        P->Print_2d_vec(col_ind_2d, "col_ind_2d_R_preT.txt");
        hout << "=========================================================================================" << endl;
        hout << "=========================================================================================" << endl;
        hout << "Check_repeated_col_ind_2d" << endl;
        if (!Check_repeated_col_ind_2d(nodes, structure, contacts_point_tmp, point_list, LM_matrix, col_ind_2d, values_2d)) {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Check_repeated_col_ind_2d " << R_flag << endl;
            return 0;
        }
        hout << "=========================================================================================" << endl;
        hout << "=========================================================================================" << endl;
        //CHECK
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        //*/
        
        //hout << "Fill_2d_matrices_gch"<<endl;
        //Fill the 2D matrices with the contributions of the GNPs
        if (!Fill_2d_matrices_gch(structure, point_list, radii, cluster_gch, LM_matrix, electric_param, hybrid_particles, col_ind_2d, values_2d, diagonal)) {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_gch" << endl;
            return 0;
        }

        //Sort resistances
        sort(resistances.begin(),resistances.end());
        //P->Print_1d_vec(resistances, "resistances_DEA_Re.txt");
        //P->Print_2d_vec(col_ind_2d, "col_ind_2d_R.txt");
        //P->Print_1d_vec(diagonal, "diagonal_R.txt");
        //P->Print_2d_vec(values_2d, "values_2d_R.txt");
        //delete P;
        Export_matlab_sparse_matrix(col_ind_2d, values_2d, diagonal, "Matrix_R.dat");
        
    } else if (R_flag == 0) {
        //If R_flag == 0, then calculate unit resistances
        
        //Printer *P = new Printer;
        //hout << "Fill_2d_matrices_unit_resistors"<<endl;
        //Fill the 2D matrices with the contributions of the CNTs
        if (!Fill_2d_matrices_unit_resistors(elements, point_list, cluster, LM_matrix, col_ind_2d, values_2d, diagonal, contacts_point)){
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices" << endl;
            return 0;
        }
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        /*/CHECK
        P->Print_2d_vec(col_ind_2d, "col_ind_2d_unit_preT.txt");
        P->Print_2d_vec(contacts_point_tmp, "contacts_point_tmp.txt");
        P->Print_2d_vec(elements, "elements.txt");
        P->Print_2d_vec(structure, "structure.txt");
        hout << "*****************************************************************************************" << endl;
        hout << "*****************************************************************************************" << endl;
        hout << "Check_repeated_col_ind_2d" << endl;
        if (!Check_repeated_col_ind_2d(nodes, structure, contacts_point_tmp, point_list, LM_matrix, col_ind_2d, values_2d)) {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Check_repeated_col_ind_2d " << R_flag << endl;
            return 0;
        }
        hout << "*****************************************************************************************" << endl;
        hout << "*****************************************************************************************" << endl;
        //CHECK
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        //*/

        //hout << "Fill_2d_matrices_gch_unit_resistors"<<endl;
        //Fill the 2D matrices with the contributions of the GNPs
        if (!Fill_2d_matrices_gch_unit_resistors(structure, point_list, cluster_gch, LM_matrix, hybrid_particles, col_ind_2d, values_2d, diagonal)) {
            hout << "Error in Fill_sparse_stiffness_matrix when calling Fill_2d_matrices_gch_unit_resistors" << endl;
            return 0;
        }
        //P->Print_2d_vec(col_ind_2d, "col_ind_2d_unit.txt");
        //P->Print_1d_vec(diagonal, "diagonal_unit.txt");
        ////P->Print_2d_vec(values_2d, "values_2d_unit.txt");
        //P->Print_2d_vec(elements_tunnel, "elements_tunnel.txt");
        //delete P;
        Export_matlab_sparse_matrix(col_ind_2d, values_2d, diagonal, "Matrix_unit.dat");
        
    } else {
        hout << "Error in Fill_sparse_stiffness_matrix. R_flag can only be 1 or 0. Actual value is " << R_flag << endl;
        return 0;
    }
    
    //------------------------------------------------------------------------
    //------------------------------------------------------------------------
    //Chech that there are no repeated nodes on each column of the stiffness matrix
    //If there are, then check if the can be added or if there is an error
    hout << "Check_repeated_col_ind_2d" << endl;
    if (!Check_repeated_col_ind_2d(nodes, structure, contacts_point_tmp, point_list, LM_matrix, col_ind_2d, values_2d)) {
        hout << "Error in Fill_sparse_stiffness_matrix when calling Check_repeated_col_ind_2d " << R_flag << endl;
        return 0;
    }
    contacts_point_tmp.clear();
    
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
    
    //hout << "From_2d_to_1d_vectors"<<endl;
    From_2d_to_1d_vectors(col_ind_2d, values_2d, KEFT, col_ind, row_ptr, values, diagonal);
    //hout << "From_2d_to_1d_vectors done"<<endl;
    
    return 1;    
}
//This function adds the contributions of the CNT and junction resistors to the stiffness matrix
int Direct_Electrifying::Fill_2d_matrices(const vector<vector<long int> > &elements, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &cluster, const vector<int> &LM_matrix, const struct Electric_para &electric_param, const double &d_vdw, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point)
{
    //Variables
    int CNT;
    long int P1, P2, node1, node2;
    
    //Scan every CNT in the cluster
    for (long int i = 0; i < (long int)cluster.size(); i++) {
        CNT = cluster[i];
        //hout << "CNT="<<CNT<<" elements[CNT].size()="<<elements[CNT].size()<<endl;
        for (long int j = 0; j < (long int)elements[CNT].size()-1; j++) {
            //Find node numbers of the first two elements
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
                double Re = Calculate_resistance_cnt(point_list, P1, P2, radii[CNT], electric_param.resistivity_CF);
                resistances.push_back(Re);
                Add_elements_to_sparse_stiffness(node1, node2, Re, col_ind_2d, values_2d, diagonal);
                
                //Check if the current node1 has any contacts and add the corresponding contributions to the
                //stiffness matrix
                //hout << "Check_for_other_elements nested loop "<<j<<endl;
                Check_for_other_elements(point_list, radii, LM_matrix, electric_param, d_vdw, P1, node1, col_ind_2d, values_2d, diagonal, contacts_point);
            }
        }
        //Check if the last node has any contacts and add the corresponding contributions to the
        //stiffness matrix
        P1 = elements[CNT].back();
        node1 = LM_matrix[P1];
        Check_for_other_elements(point_list, radii, LM_matrix, electric_param, d_vdw, P1, node1, col_ind_2d, values_2d, diagonal, contacts_point);
        //hout << "Check_for_other_elements "<<endl;
    }
    
    return 1;
}
//Add contributions from the resistors formed by the triangulation on the GNP
int Direct_Electrifying::Fill_2d_matrices_gch(const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &cluster_gch, const vector<int> &LM_matrix, const struct Electric_para &electric_param, vector<GCH> &hybrid_particles, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal)
{
    //Variables
    long int P1, P2, node1, node2;
    double Re;
    
    //Triangulation object
    Triangulation *delaunay = new Triangulation;
    
    //Scan every hybrid particle, perform the triangulations and add the elements to the stiffness matrix
    for (long int i = 0; i < (long int)cluster_gch.size(); i++) {
        //current hybrid particle
        int hyb = cluster_gch[i];
        
        //Perform triangulations
        if (!delaunay->Generate_3d_trangulation(point_list, structure, hybrid_particles[hyb])) {
            hout << "Error in Fill_2d_matrices_gch" << endl;
            return 0;
        }
        //hout << "Triangulation " << i << endl;
        
        //Vector to store elements from the triangulation that need to be deleted
        vector<int> to_delete;
        
        //Add elements from the top triangulation
        for (int j = 0; j < (int)hybrid_particles[hyb].triangulation.size(); j++) {
            P1 = hybrid_particles[hyb].triangulation[j][0];
            node1 = LM_matrix[P1];
            P2 = hybrid_particles[hyb].triangulation[j][1];
            node2 = LM_matrix[P2];
            //Sometimes a CNT seed is also a boundary point, in that case the triangulation will result
            //in having repeated elements in the stiffness matrix, since there is no tunneling on a bonudary point
            //Thus ignore elements that have a boundary node
            //Equivalently only add elements when both nodes are non-boundary nodes
            if (node1 >= (long int)reserved_nodes && node2 >= (long int)reserved_nodes) {
                //Resistance of the "conduction band" in the GNP
                resistances.push_back(Re);
                Re = Calculate_resistance_gnp(point_list[P1], point_list[P2], radii[point_list[P1].flag], radii[point_list[P2].flag], hybrid_particles[hyb], electric_param);
                Add_elements_to_sparse_stiffness(node1, node2, Re, col_ind_2d, values_2d, diagonal);
            } else {
                //If one node is a boundary node this element needs to be deleted from the triangulation
                to_delete.push_back(j);
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
//This function adds the contributions of the CNT and junction resistors to the stiffness matrix
//Unit resistors are used
int Direct_Electrifying::Fill_2d_matrices_unit_resistors(const vector<vector<long int> > &elements, const vector<Point_3D> &point_list, const vector<int> &cluster, const vector<int> &LM_matrix, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point)
{
    //Variables
    int CNT;
    long int P1, P2, node1, node2;
    
    //Scan every CNT in the cluster
    for (long int i = 0; i < (long int)cluster.size(); i++) {
        CNT = cluster[i];
        //hout << "CNT="<<CNT<<" elements[CNT].size()="<<elements[CNT].size()<<endl;
        for (long int j = 0; j < (long int)elements[CNT].size()-1; j++) {
            //Find node numbers of the first two elements
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
            //Equivalently only add elements when both nodes are non-boundary nodes
            if (node1 != node2) {
                Add_elements_to_sparse_stiffness(node1, node2, 1, col_ind_2d, values_2d, diagonal);
                
                //Check if the current node1 has any contacts and add the corresponding contributions to the
                //stiffness matrix
                //hout << "Check_for_other_elements nested loop "<<j<<endl;
                Check_for_other_unit_elements(point_list, LM_matrix, P1, node1, col_ind_2d, values_2d, diagonal, contacts_point);
            }
        }
        //Check if the last node has any contacts and add the corresponding contributions to the
        //stiffness matrix
        //hout << "Check_for_other_elements "<<endl;
        P1 = elements[CNT].back();
        node1 = LM_matrix[P1];
        Check_for_other_unit_elements(point_list, LM_matrix, P1, node1, col_ind_2d, values_2d, diagonal, contacts_point);
    }
    
    return 1;
}
//Add contributions from the resistors formed by the triangulation on the GNP
//Unit resistors are used
int Direct_Electrifying::Fill_2d_matrices_gch_unit_resistors(const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<int> &cluster_gch, const vector<int> &LM_matrix, vector<GCH> &hybrid_particles, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal)
{
    //Variables
    long int P1, P2, node1, node2;
    
    //Triangulation object
    Triangulation *delaunay = new Triangulation;
    
    //Scan every hybrid particle, perform the triangulations and add the elements to the stiffness matrix
    for (long int i = 0; i < (long int)cluster_gch.size(); i++) {
        //current hybrid particle
        int hyb = cluster_gch[i];
        
        //Perform triangulations
        if (!delaunay->Generate_3d_trangulation(point_list, structure, hybrid_particles[hyb])) {
            hout << "Error in Fill_2d_matrices_gch" << endl;
            return 0;
        }
        
        //Add elements from the top triangulation
        for (int j = 0; j < (int)hybrid_particles[hyb].triangulation.size(); j++) {
            P1 = hybrid_particles[hyb].triangulation[j][0];
            node1 = LM_matrix[P1];
            P2 = hybrid_particles[hyb].triangulation[j][1];
            node2 = LM_matrix[P2];
            //Sometimes a CNT seed is also a boundary point, in that case the triangulation will result
            //in having repeated elements in the stiffness matrix, since there is no tunneling on a bonudary point
            //Thus ignore elements that have a boundary node
            //Equivalently only add elements when both nodes are non-boundary nodes
            if (node1 >= (long int)reserved_nodes && node2 >= (long int)reserved_nodes) {
                Add_elements_to_sparse_stiffness(node1, node2, 1, col_ind_2d, values_2d, diagonal);
            }
        }
    }
    
    //Delete triangulation object
    delete delaunay;
    
    return 1;
}
//Check if the current node1 has any contacts and add the corresponding contributions to the stiffness matrix
void Direct_Electrifying::Check_for_other_elements(const vector<Point_3D> &point_list, const vector<double> &radii, const vector<int> &LM_matrix, const struct Electric_para &electric_param, const double d_vdw, const long int &P1, const long int &node1, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point)
{
    if (contacts_point[P1].size()) {
        for (long int k = 0; k < (long int)contacts_point[P1].size(); k++) {
            long int P2 = contacts_point[P1][k];
            //hout <<" contact P2="<<P2;
            long int node2 = LM_matrix[P2];
            //When using the actual resistances, there are two cases when a tunnel resistor should be ignored:
            //Case 1:
            //This fuction applies the DEA on the backbone, thus some CNTs will be deleted and so their contacts
            //Instead of updating the vector of contacts, ignore the deleted contacts
            //The deleted contacts will have a node2 of -1
            //
            //Case 2 (same as in unit resistors):
            //Sometimes there is tunneling between a boundary point and a non-boundary point
            //If the CNT that contains the non-boundary point happens to be in contact with the
            //boundary (this is actually likely to happen), then there will be two resistors connecting
            //the same nodes: a CNT resistor and a tunner resistor
            //A tunnel from a boundary to a CNT that is already in contact with the boundary does not
            //seem to be physically correct, so the tunnels between the boundary and a CNT are avoided
            //Boundary nodes have a node number below the variable reserved_nodes, so to include a tunnel
            //resistors, both of its nodes should not be reserved nodes
            //
            //Note that case 1 is solved with the condition of case 2, since -1 is less than the reserved nodes
            //so I use the same if-statement as with unit resistors
            if (node1 >= (long int)reserved_nodes && node2 >= (long int)reserved_nodes) {
                //Calculate tunnel resistance
                double Re = Calculate_resistance_tunnel(radii, electric_param, point_list[P1], point_list[P2], d_vdw);
                resistances.push_back(Re);
                //Add the elements to the sparse vectors
                Add_elements_to_sparse_stiffness(node1, node2, Re, col_ind_2d, values_2d, diagonal);
                //Remove the contac tha was used so it is not used again in the future
                Remove_from_vector(P1, contacts_point[P2]);
                //hout << "Removed ";
                //Add tunnel element
                vector<long int> empty;
                elements_tunnel.push_back(empty);
                elements_tunnel.back().push_back(P1);
                elements_tunnel.back().push_back(P2);
            }
        }
        //Remove all contacts of P1
        contacts_point[P1].clear();
    }
}
//Check if the current node1 has any contacts and add the corresponding contributions to the stiffness matrix
//This function considers all resistances to have value 1
void Direct_Electrifying::Check_for_other_unit_elements(const vector<Point_3D> &point_list, const vector<int> &LM_matrix, const long int &P1, const long int &node1, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point)
{
    if (contacts_point[P1].size()) {
        for (long int k = 0; k < (long int)contacts_point[P1].size(); k++) {
            long int P2 = contacts_point[P1][k];
            //hout <<" contact P2="<<P2;
            long int node2 = LM_matrix[P2];
            //Sometimes there is tunneling between a boundary point and a non-boundary point
            //If the CNT that contains the non-boundary point happens to be in contact with the
            //boundary (this is actually likely to happen), then there will be two resistors connecting
            //the same nodes: a CNT resistor and a tunner resistor
            //A tunnel from a boundary to a CNT that is already in contact with the boundary does not
            //seem to be physically correct, so the tunnels between the boundary and a CNT are avoided
            //Boundary nodes have a node number below the variable reserved_nodes, so to include a tunnel
            //resistors, both of its nodes should not be reserved nodes
            if (node1 >= (long int)reserved_nodes && node2 >= (long int)reserved_nodes) {
                //Add the elements to the sparse vectors
                Add_elements_to_sparse_stiffness(node1, node2, 1, col_ind_2d, values_2d, diagonal);
                //Remove the contac tha was used so it is not used again in the future
                Remove_from_vector(P1, contacts_point[P2]);
                //hout << "Removed ";
                //Add tunnel element
                vector<long int> empty;
                elements_tunnel.push_back(empty);
                elements_tunnel.back().push_back(P1);
                elements_tunnel.back().push_back(P2);
            }
        }
        //Remove all contacts of P1
        contacts_point[P1].clear();
    }
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
/*/This function calculates the electrical resistance due to tunneling using Li et al. approach
double Direct_Electrifying::Calculate_resistance_tunnel(const vector<double> &radii, const struct Electric_para &electric_param, const Point_3D &P1, const Point_3D &P2, const double &d_vdw)
{
    //Calculate the distance between the points in contact
    double distance = P1.distance_to(P2);
    
    //Calculate the minimum possible distance between points (r1 + r2 + d_vdw)
    int CNT1 = P1.flag, CNT2 = P2.flag;
    double min_dist = radii[CNT1] + radii[CNT2] + d_vdw;
    //Check if the distance between points is below the sum of radii plus van der Waals distance
    //This happens when the penetrating model is used
    if (distance < min_dist)
        //If the distance between the points is below min_dist, then set this distance equal to min_dist
        distance = min_dist;
    
    //Calculate the thickness of the film = distance - r1 - r2 - 2*d_vdw = distance - min_dist - d_vdw
    double film = distance - min_dist - d_vdw;
    //Check if it is a negative value then the film thickness is zero so the resistance due to tunnel is zero
    //this might happen due to using the penetrating model or numerical error
    if (film < Zero)
        return 0;
    
    //Approximate cross section of the tunnel as a rectangle formed by the CNTs
    double Ac = radii[CNT1]*radii[CNT2];
    
    //Approximate the voltage across the film, U
    double U = electric_param.e_charge*film/(Ac*electric_param.K_polymer*electric_param.e0_vacuum);
    
    //Calculate t1, t2, Dt
    double t1 = 6/(electric_param.K_polymer*electric_param.CNT_work_function);
    double Dt = film*(1 - 46/(3*electric_param.CNT_work_function*electric_param.K_polymer*film + 20 - 2*U*electric_param.K_polymer*film) );
    double t2 = Dt + t1;
    
    //Calculate phi
    double phi = electric_param.CNT_work_function - U*(t1+t2)/(2*film);
    double tmp = t2*(film-t1)/(t1*(film-t2));
    phi = phi - 5.75*log(tmp)/(electric_param.K_polymer*Dt);
    
    //Calculate the current density, J
    double J = phi*exp(-1.025*Dt*sqrt(phi)) - (phi+U)*exp(-1.025*Dt*sqrt(phi+U));
    J = 6.2e10*J/(Dt*Dt);
    
    //Calculate tunnel resistance
    return U/(Ac*J);
}//*/
//This function calculates the electrical resistance due to tunneling using Hu et al. approach
double Direct_Electrifying::Calculate_resistance_tunnel(const vector<double> &radii, const struct Electric_para &electric_param, const Point_3D &P1, const Point_3D &P2, const double &d_vdw)
{
    //Variable to store the separation between CNTs
    double separation;
    
    //Calculate the distance between the points in contact
    double distance = P1.distance_to(P2);
    
    //Calculate the minimum possible distance between points (r1 + r2 + d_vdw)
    double min_dist = radii[P1.flag] + radii[P2.flag] + d_vdw;
    //Check if the distance between points is below the sum of radii plus van der Waals distance
    //This happens when the penetrating model is used
    if (distance < min_dist)
        //If the distance between the points is below min_dist, then set the separation equal to the van der Waals distance
        separation = d_vdw;
    else
        //If the distance between the points is above min_dist, then calculate the separation between CNTs
        separation = distance - min_dist;
    
    //==============================================================================
    //==============================================================================
    //Input parameters are scaled to avoid numerical errors when reading, therefore the following steps are taken
    
    //Calculate quantity associted with a squared root
    double sqrt_tmp = sqrt(2*electric_param.e_mass*electric_param.lambda_barrier*electric_param.e_charge);
    
    //Calculate the exponential term
    double exp_tmp = exp(4000*PI*separation*sqrt_tmp/electric_param.h_plank);
    
    //Calculate term that multiplies the exponential
    double denominator_tmp = radii[P1.flag]*radii[P2.flag]*electric_param.e_charge*electric_param.e_charge*sqrt_tmp;
    double mult_tmp = 10*electric_param.h_plank*electric_param.h_plank*separation/denominator_tmp;
    
    //Calculate tunnel resistance
    return mult_tmp*exp_tmp;
}
//This function calculates the resistance that comes from the triangulation on the GNPs
//It uses the resistivities along the surface and along the thickness so the resistance
//can be calculated along any direction the triangulation edge may have
double Direct_Electrifying::Calculate_resistance_gnp(const Point_3D &P1, const Point_3D &P2, const double &rad1, const double &rad2, const GCH &hybrid, const struct Electric_para &electric_param)
{
    //Calculate the distance between the points in contact
    double L = P1.distance_to(P2);
    
    //Unit vector in the direction from P1 to P2
    Point_3D u_direction = (P2 - P1)/L;
    
    //Multiply resistivity tensor by the u_direction vector
    Point_3D direction(u_direction.x*electric_param.resistivity_GNP_surf, u_direction.y*electric_param.resistivity_GNP_surf, u_direction.z*electric_param.resistivity_GNP_t);
    
    //Calculate cross sectional area of conduction band parallel to the GNP surface
    double A = PI*(rad1 + rad2)*(rad1 + rad2)/4;
    
    //Calculate the angle of the actual cross sectional area
    Point_3D x_dir(1.0,0.0,0.0); //unit vector in the x direction, parallel to the surface of the conduction band
    x_dir.rotation(hybrid.rotation, hybrid.center); //rotate the unit vector so that it is the global coordinate system
    double cosA = abs(u_direction.dot(x_dir)); //cosine of the angle
    
    //The resistance is the magnitude of the vector direction multiplied by the distance between points
    // and divided by the actual cross sectional area
    return direction.distance_to(0, 0, 0)*L/(A*cosA*cosA);
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
//Test function to find repeated elements in the vector col_ind_2d
int Direct_Electrifying::Check_repeated_col_ind_2d(const int &nodes, const vector<vector<long int> > &structure, vector<vector<long int> > &contacts_point, const vector<Point_3D> &point_list, const vector<int> &LM_matrix, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d)
{
    //Create inverse mapping for the LM_matrix, i.e., from node number to point number
    //Initialize the LM_inverse with the number of nodes
    vector<long int> LM_inverse(nodes, -1);
    
    //Fill the LM_inverse matrix
    for (long int i = 0; i < (long int)LM_matrix.size(); i++) {
        int current_node = LM_matrix[i];
        if (current_node != -1) {
            LM_inverse[current_node] = i;
        }
    }
    
    //Loop over the rows of the 2D sparse matrix
    for (int i = 0; i < (int)col_ind_2d.size(); i++) {
        //Find repeated elements in the vector col_ind_2d[i]
        vector<long int> repeated_elements;
        vector<int> indices;
        Find_repeated_elements(col_ind_2d[i], repeated_elements, indices);
        //If there are repeated elements print them in the output file
        if (repeated_elements.size() == 1) {
            //Get the point number of the node with repeated indices
            long int P = LM_inverse[i];
            //Get the CNT of P
            int CNT = point_list[P].flag;
            //Get the node numbers of the endpoints of the CNT
            int node_front = LM_matrix[elements[CNT].front()];
            int node_back = LM_matrix[elements[CNT].back()];
            //Check if the two endpoints of the CNT are the same and they are boundary nodes
            if (node_back == node_front && node_front < reserved_nodes) {
                //If the two endpoints are in the same node and are boundary nodes, that means that the CNT has two contacts with the boundary and one contact in between
                //This is the only case in with two nodes have two resistors
                //So in this case add the quantities in the vector values
                values_2d[i][indices[0]] = values_2d[i][indices[0]] + values_2d[i][indices[1]];
                
                //Remove repeated element from the vectors
                //Te repeated element is the second index in the vector indices
                col_ind_2d[i].erase(col_ind_2d[i].begin()+indices[1]);
                values_2d[i].erase(values_2d[i].begin()+indices[1]);
                
            }
        } else if (repeated_elements.size()) {
            //If there is more than one element repeated then something else is going wrong and will need this output to help me understand
            hout << "col_ind_2d[" << i << "] (";
            hout <<"point="<<LM_inverse[i]<<" CNT="<<point_list[LM_inverse[i]].flag;
            hout <<" CNT[0]="<<structure[point_list[LM_inverse[i]].flag].front();
            hout <<" CNT[last]="<<structure[point_list[LM_inverse[i]].flag].back()<<")";
            hout << endl;
            hout <<"("<<point_list[LM_inverse[i]].x<<", "<<point_list[LM_inverse[i]].y<<", ";
            hout <<point_list[LM_inverse[i]].z<<')'<<endl;
            hout << "Total contacts: "<<contacts_point[LM_inverse[i]].size()<<endl;
            for (int k = 0; k < (int)contacts_point[LM_inverse[i]].size(); k++) {
                long int contact_tmp = contacts_point[LM_inverse[i]][k];
                int CNT_tmp = point_list[contact_tmp].flag;
                hout << "Contact "<<k<<": P="<<contact_tmp<<" node2="<<LM_matrix[contact_tmp]<<" CNT="<<CNT_tmp;
                hout << " CNT[0]="<<structure[CNT_tmp].front()<<" CNT[last]="<<structure[CNT_tmp].back()<<endl;
                hout <<"("<<point_list[contact_tmp].x<<", "<<point_list[contact_tmp].y<<", ";
                hout <<point_list[contact_tmp].z<<')'<<endl;
            }
            for (int j = 0; j < (int)repeated_elements.size(); j++) {
                long int P = LM_inverse[repeated_elements[j]];
                int CNT = point_list[P].flag;
                hout <<"Repeated node="<<repeated_elements[j]<<" point="<<P<<" CNT="<<CNT;
                hout <<" CNT[0]="<<structure[CNT].front()<<" CNT[last]="<<structure[CNT].back()<<" d="<<point_list[P].distance_to(point_list[LM_inverse[i]]);
                hout << endl;
                hout <<"("<<point_list[P].x<<", "<<point_list[P].y<<", "<<point_list[P].z<<')'<<endl;
                //Print out the contacts of repeated element
                hout << "Total contacts of repeated node: "<<contacts_point[P].size()<<endl;
                for (int k = 0; k < (int)contacts_point[P].size(); k++) {
                    long int contact_tmp = contacts_point[P][k];
                    int CNT_tmp = point_list[contact_tmp].flag;
                    hout << "Contact "<<k<<": P="<<contact_tmp<<" CNT="<<CNT_tmp;
                    hout << " CNT[0]="<<structure[CNT_tmp].front()<<endl;
                }
            }
            hout << endl;
        }
    }
    return 1;
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
void Direct_Electrifying::From_2d_to_1d_vectors(const vector<vector<long int> > &col_ind_2d, const vector<vector<double> > &values_2d, vector<vector<double> > &KEFT, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal)
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
    //From the matrix equations R = b = - KEFT*VEF
    for (long int i = 0; i < (long int)KEFT.size(); i++) {
        //the first element of VEF (i.e. VEF[0]) is zero, so it can be skipped to save computational time
        for (long int j = 1; j < (long int)KEFT[i].size(); j++) {
            R[i] = R[i] - (KEFT[i][j]*VEF[j]);
        }
        //The search direction of the CG is initialized with the first value of the residual
        P[i] = R[i];
    }

    //=========================================
    // Conjugate Gradient Algorithm
    Conjugate_gradient(nodes, col_ind, row_ptr, values, diagonal, R, P);
    
    //The known boundary conditions are added at the beginning of the solution
    for (int i = reserved_nodes-1; i >=0; i--) {
        voltages.insert(voltages.begin(), VEF[i]);
    }
    
    //Print the voltages
    Printer *Pr = new Printer;
    if (R_flag) {
        Pr->Print_1d_vec(voltages, "voltages_R.txt");
    } else {
        Pr->Print_1d_vec(voltages, "voltages_unit.txt");
    }
    delete Pr;
    
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
void Direct_Electrifying::Conjugate_gradient(long int nodes, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, vector<double> &R, vector<double> &P)
{
    //Variables of the algorithm
    vector<double> AP;
    AP.assign(nodes-reserved_nodes,0);
    double alpha, beta, rr0, rr;
    voltages.clear();
    voltages.assign(nodes-reserved_nodes, 0);
    
    //Maximum number of iterations for the CG
    long int max_iter = 100*nodes;
    //Iteration variable
    long int k;
    //Variable to check the status of the CG
    int test = 50000, test_inc = 50000;
    
    //Initial residual
    double R0 = 1.0E-12*sqrt(V_dot_v(R, R));
    //double R0 = Zero*sqrt(V_dot_v(R, R));
    //double R0 = 1.0E-10;
    //hout << "R0 = " << R0 << endl;
    
    for (k = 1; k <= max_iter; k++) {
        //Calculate Ap
        spM_V_SSS(P, row_ptr, col_ind, diagonal, values, AP);
        //Calculate norm or residual of step k-1. Will be used later as convergence criteria and to calculate beta
        rr0 = V_dot_v(R, R);
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
        //Status update: print every hundred iterations
        if ( k == test){
            hout << "CG iteration " << k << endl;
            test = test + test_inc;
        }
        //Convergence criteria
        if (sqrt(rr) <= R0)
            break;
        //Improvement of step
        beta = rr/rr0;
        //Search direction
        //P = R + P*beta;
        W_plus_aV(R, beta, P);
    }
    
    if (k >= max_iter)
        hout << "CG reached maximum number of iterations" << endl;
    hout << "CG iterations: " << k << endl;
    //hout << "RR = " << sqrt(rr) << endl;
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