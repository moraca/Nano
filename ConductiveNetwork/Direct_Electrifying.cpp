//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Direct_Eletrifying.h
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
int Direct_Electrifying::Calculate_voltage_field(const vector<vector<long int> > &structure, vector<vector<long int> > &contacts_point, const vector<vector<short int> > &boundary_flags, const vector<int> &cluster, const vector<double> &radii, int family, const struct Electric_para &electric_param)
{
    //First we need to prepare the matrices that the direct electrifying needs
    //The first matrix will be the local mapping (LM) matrix. This matrix maps from point number in the structure
    //to node number for the solver
    //The number of reserved nodes is calculated. These is the number of boundaries with prescribed voltage
    reserved_nodes = Get_global_nodes(family);
    //This variable is used to assing node numbers and after calling Get_LM_matrix, it will contain the number of nodes in the network
    int global_nodes = reserved_nodes;
    //Initialize the size of the LM matrix to be equal to the number of points
    LM_matrix.assign(contacts_point.size(), -1);
    //hout << "LM_matrix.size()="<<LM_matrix.size()<<" cluster.size()="<<cluster.size()<<endl;
    //Initialize the vector boundary_node_map
    Initialize_boundary_node_map();
    //Initialize the size of the elements matrix to be equal to the number of CNTs
    vector<long int> empty;
    elements.assign(structure.size(), empty);
    if(!Get_LM_matrix(structure, contacts_point, boundary_flags,cluster, family , global_nodes, LM_matrix, elements)){
        hout << "Error in Calculate_voltage_field" << endl;
        return 0;
    }
    
    //hout << "Fill_sparse_stiffness_matrix"<<endl;
    //Variables for using the SSS for the sparse matrix
    vector<long int> col_ind, row_ptr;
    vector<double> values, diagonal, R;
    vector<vector<double> > KEFT;
    //With the LM matrix, now fill the sparse stiffness matrix
    Fill_sparse_stiffness_matrix(structure, contacts_point, cluster, global_nodes, LM_matrix, elements, KEFT, col_ind, row_ptr, values, diagonal);
    
    //hout << "Solve_DEA_equations_CG_SSS"<<endl;
    //This is where the actual direct electrifying algorithm (DEA) takes place
    Solve_DEA_equations_CG_SSS(global_nodes, col_ind, row_ptr, values, diagonal, electric_param, KEFT);
    //hout << "end"<<endl;
    
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
int Direct_Electrifying::Get_LM_matrix(const vector<vector<long int> > &structure, const vector<vector<long int> > &contacts_point, const vector<vector<short int> > &boundary_flags, const vector<int> &cluster, int family, int &global_nodes, vector<int> &LM_matrix, vector<vector<long int> > &elements)
{
    //Variables
    int CNT;
    long int P;
    
    for (int i = 0; i < (int)cluster.size(); i++) {
        CNT = cluster[i];
        
        //Scan the CNT contacts
        for (int j = 0; j < (int)structure[CNT].size(); j++) {
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
        //Check if the endpoints of the CNT are already in the elements matrix, otherwise add them to the elements matrix and the LM_matrix
        if (elements[CNT].front() != structure[CNT].front()){
            P = structure[CNT].front();
            elements[CNT].insert(elements[CNT].begin(), P);
            Add_point_to_LM_matrix(P, family, boundary_flags, global_nodes, LM_matrix);
        }
        if (elements[CNT].back() != structure[CNT].back()){
            P = structure[CNT].back();
            elements[CNT].push_back(P);
            Add_point_to_LM_matrix(P, family, boundary_flags, global_nodes, LM_matrix);
        }
        
        //Check that the elements vector has a valid size
        if (elements[CNT].size() <= 1) {
            hout << "Error in Get_LM_matrix. The vector elements["<<CNT<<"] has size "<< elements[CNT].size();
            hout << " but it has to have at least two elements." << endl;
            hout << "\tCNT has "<<structure[CNT].size()<<" points"<<endl;
        }
        //hout<<endl;
    }
    
    return 1;
}

int Direct_Electrifying::Get_global_nodes(int family)
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
void Direct_Electrifying::Fill_sparse_stiffness_matrix(const vector<vector<long int> > &structure, vector<vector<long int> > &contacts_point, const vector<int> &cluster, int nodes, const vector<int> &LM_matrix, const vector<vector<long int> > &elements, vector<vector<double> > &KEFT, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal)
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
    
    
    //hout << "Fill_2d_matrices"<<endl;
    //Fill the 2D matrices
    Fill_2d_matrices(elements, cluster, LM_matrix, col_ind_2d, values_2d, diagonal, contacts_point);
    
    
    //------------------------------------------------------------------------
    //Convert from 2D vectors to 1D vectors
    
    //Initialize vectors
    values.clear();
    col_ind.clear();
    row_ptr.clear();
    //The two first elements of row_ptr are zero
    //row_ptr.push_back(0);
    row_ptr.push_back(0);
    //empty_double.push_back(0);
    KEFT.clear();
    vector<double> zeros(reserved_nodes,0);
    KEFT.assign(nodes-reserved_nodes, zeros);
    
    //hout << "From_2d_to_1d_vectors"<<endl;
    From_2d_to_1d_vectors(col_ind_2d, values_2d, KEFT, col_ind, row_ptr, values, diagonal);
    
    
}

void Direct_Electrifying::Fill_2d_matrices(const vector<vector<long int> > &elements, const vector<int> &cluster, const vector<int> &LM_matrix, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point)
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
            Add_elements_to_sparse_stiffness(node1, node2, col_ind_2d, values_2d, diagonal);
            
            //Check if the current node1 has any contacts and add the corresponding contributions to the
            //stiffness matrix
            //hout << "Check_for_other_elements nested loop "<<j<<endl;
            Check_for_other_elements(LM_matrix, P1, node1, col_ind_2d, values_2d, diagonal, contacts_point);
        }
        //Check if the last node has any contacts and add the corresponding contributions to the
        //stiffness matrix
        //hout << "Check_for_other_elements "<<endl;
        P1 = elements[CNT].back();
        node1 = LM_matrix[P1];
        Check_for_other_elements(LM_matrix, P1, node1, col_ind_2d, values_2d, diagonal, contacts_point);
    }
}

//Check if the current node1 has any contacts and add the corresponding contributions to the
//stiffness matrix
void Direct_Electrifying::Check_for_other_elements(const vector<int> &LM_matrix, long int P1, long int node1, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal, vector<vector<long int> > &contacts_point)
{
    if (contacts_point[P1].size()) {
        for (long int k = 0; k < (long int)contacts_point[P1].size(); k++) {
            long int P2 = contacts_point[P1][k];
            //hout <<" contact P2="<<P2;
            long int node2 = LM_matrix[P2];
            //Add the elements to the sparse vectors
            Add_elements_to_sparse_stiffness(node1, node2, col_ind_2d, values_2d, diagonal);
            //Remove the contac tha was used so it is not used again in the future
            Remove_from_vector(P1, contacts_point[P2]);
            //hout << "Removed ";
        }
        //Remove all contacts of P1
        contacts_point[P1].clear();
    }
    
}

void Direct_Electrifying::Add_elements_to_sparse_stiffness(long int node1, long int node2, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values_2d, vector<double> &diagonal)
{
    //Add the diagonal elements of the stiffness matrix
    diagonal[node1] += 1;
    diagonal[node2] += 1;
    
    //Add the off diagonal elements of the stiffness matrix
    if (node1 > node2) {
        col_ind_2d[node1].push_back(node2);
        //This is the resistance between the two nodes
        values_2d[node1].push_back(-1);
    } else {
        col_ind_2d[node2].push_back(node1);
        //This is the resistance between the two nodes
        values_2d[node2].push_back(-1);
    }
    //hout << "Added ";
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
    for (long int i = reserved_nodes; i < (long int)col_ind_2d.size(); i++) {
        //hout << "for(i)=" <<i << " col_ind_2d[i].size()="<<col_ind_2d[i].size()<< endl;
        for (long int j = 0; j < (long int)col_ind_2d[i].size(); j++) {
            //hout << "for(j)="<<j<<' ';
            if (col_ind_2d[i][j] >= reserved_nodes) {
                //hout << "if1 ";
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
void Direct_Electrifying::Solve_DEA_equations_CG_SSS(long int nodes, const vector<long int> &col_ind, const vector<long int> &row_ptr, const vector<double> &values, const vector<double> &diagonal, const struct Electric_para &electric_param, vector<vector<double> > &KEFT)
{
    
    //=========================================
    // Set up variables for the Conjugate Gradient Algorithm
    
    //Voltage applied to the sample
    vector<double> VEF;
    //Initializa the prescribed voltage boundary conditions
    Get_voltage_vector(nodes, VEF);
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
    Pr->Print_1d_vec(voltages, "voltages.txt");
    delete Pr;
}

//This function creates a voltage vector depending on the number of prescribed boundary conditios
void Direct_Electrifying::Get_voltage_vector(const long int &nodes, vector<double> &voltages)
{
    //Clear the vector of voltages
    voltages.clear();
    for (int i = 0; i < reserved_nodes; i++) {
        voltages.push_back( ((double)nodes)*((double)i) );
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
    long int max_iter = 10*nodes;
    //Iteration variable
    long int k;
    //Variable to check the status of the CG
    int test = 500, test_inc = 250;
    
    //Initial residual
    double R0 = 1.0E-12*sqrt(V_dot_v(R, R));
    //double R0 = Zero*sqrt(fabs(V_dot_v(R, R)));
    //double R0 = 1.0E-10;
    
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
}

//This function calculates the dot product of two vectors
double Direct_Electrifying::V_dot_v(vector<double> A, vector<double> B)
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
void Direct_Electrifying::spM_V_SSS(vector<double> V, vector<long int> rowptr, vector<long int> colind, vector<double> diagonal, vector<double> values, vector<double> &R)
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
void Direct_Electrifying::V_plus_aW(vector<double> W, double a, vector<double> &V)
{
    for (int i = 0; i < (int)V.size(); i++) {
        V[i] = V[i] + a*W[i];
    }
}

//This function solves the following operation: V = W + a*V
//where V and W are vectors and a is a scalar
void Direct_Electrifying::W_plus_aV(vector<double> W, double a, vector<double> &V)
{
    for (int i = 0; i < (int)V.size(); i++) {
        V[i] = W[i] + a*V[i];
    }    
}