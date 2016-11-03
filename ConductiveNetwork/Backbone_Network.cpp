//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Backbone_Network.cpp
//OBJECTIVE:	To determine the backbone network and dead branches in the percolation network
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Backbone_Network.h"

//This function scans all CNTs inside one cluster and saves the indices of the points for the dead branches and conducting segmenents
//This function takes advantage of the fact that CNTs are numbered consecutively
//The vector dead_indices can only have sizes 0, 2 or 4 in case the whole CNT conducts, or it has one dead branch or it has two dead branches
//respectively
//The vector percolated_indices can only have sizes 0 or 2 in case the CNT does not conduct or it has a conducting segment, respectively
int Backbone_Network::Determine_backbone_network(const int &family, const int &R_flag, const int &tecplot_flag, Direct_Electrifying *DEA, const Electric_para &electric_param, const vector<GCH> &hybrid_particles, const vector<int> &cluster,const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const vector<double> &radii, vector<int> &cluster_gch, vector<double> &families_lengths, vector<double> &branches_lengths, vector<vector<long int> > &all_dead_indices, vector<vector<long int> > &all_indices, vector<vector<int> > &gnp_dead_indices, vector<vector<int> > &gnp_indices)
{
    //First check if the percolated cluster has only one CNT, in that case there are no dead branches
    //If the size of the cluster is greater than one, then probably there are branches that need to be removed
    if (cluster.size() > 1) {
        //Define the 'Zero-cutoff' of the system
        double zero_cutoff = Zero_current(R_flag, DEA, electric_param, cluster, points_in, radii, hybrid_particles, cluster_gch);
        
        //Find the dead branches of each CNT in the cluster. Save the information on the vectors dead_indices and percolated_indices
        if (!Find_dead_branches_simplified(DEA->elements, cluster, zero_cutoff)) {
            hout << "Error in Determine_backbone_network when calling Find_dead_branches" << endl;
            return 0;
        }
        
        //Find the dead GNPs
        if (!Find_dead_gnps_simplified(zero_cutoff, cluster_gch, gnp_dead_indices[family], gnp_indices[family])) {
            hout << "Error in Determine_backbone_network when calling Find_dead_gnps" << endl;
            return 0;
        }
        
    } else {
        //When there is only one percolated CNT, then the percolated_indices only contains
        //the first and last point of the CNT and dead_indices remains empty
        int CNT = cluster.front();
        //This varibale is used to initialize the vectors below
        vector<long int> empty;
        percolated_indices.push_back(empty);
        percolated_indices.back().push_back(structure[CNT].front());
        percolated_indices.back().push_back(structure[CNT].back());
        //Initialize dead_indices as it needs to have the same size as percolated_indices
        dead_indices.push_back(empty);
    }
    
    //Use the dead_indices and percolated_indices to calculate the fraction of CNTs that belong to each family
    if (!Calculate_lengths(family, points_in, families_lengths, branches_lengths)) {
        hout << "Error in Determine_backbone_network when calling Calculate_lengths" << endl;
        return 0;
    }
    //hout << "Calculate_lengths "<<endl;
    
    //Finally, add the indices to the global vectors so that they can be exported as tecplot files
    Add_indices_to_global_vectors(family, all_dead_indices, all_indices);
    //hout << "Add_indices_to_global_vectors" << endl;
    
    return 1;
}
//This function determines the cutoff for "zero current"
//This step is necessary due to floating point errors
double Backbone_Network::Zero_current(const int &R_flag, Direct_Electrifying *DEA, const struct Electric_para &electric_param, const vector<int> &cluster, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<GCH> &hybrid_particles, const vector<int> &cluster_gch)
{
    //Variables used to calculate the current
    double I, Re;
    //Variable to store the cutoff for zero voltage
    double zero_cutoff;
    //Variables
    int CNT;
    long int P1, P2;
    //Vector to store all the currents
    vector<double> currents;
    //Initialize the vectors of element currents and resistances for the CNTs
    vector<double> empty_double;
    current_e.assign(cluster.size(), empty_double);
    resistance_r.assign(cluster.size(), empty_double);
    
    //Calculate all currents from CNTs
    for (long int i = 0; i < (long int)cluster.size(); i++) {
        CNT = cluster[i];
        for (long int j = 0; j < (long int)DEA->elements[CNT].size()-1; j++) {
            P1 = DEA->elements[CNT][j];
            P2 = DEA->elements[CNT][j+1];
            //Calculate the voltage difference
            //Temporarily the voltage difference will be the current value, unless the element resistance is not 1
            I = abs(Voltage_difference(P1, P2, DEA->LM_matrix, DEA->voltages));
            //Check the resistance flag to use the appropiate resistance
            //If R_flag is 0, the voltage difference is the current, so I remains the same
            if (R_flag) {
                //If R_flag is 1, then use the actual resistance
                
                //Calculate the Resistance of the CNT segment
                Re = DEA->Calculate_resistance_cnt(point_list, P1, P2, radii[CNT], electric_param.resistivity_CF);
                //Calculate current
                I = I/Re;
                //Add resistance to the vector of element resistances
                resistance_r[i].push_back(Re);
            }
            //Add current to vector of all currents
            currents.push_back(I);
            //Add current to the vector of element currents
            current_e[i].push_back(I);
        }
    }
    
    //Calculate all currents from tunnels
    for (int i = 0; i < (int)DEA->elements_tunnel.size(); i++) {
        //Tunnel elements have only two elements per vector
        P1 = DEA->elements_tunnel[i][0];
        P2 = DEA->elements_tunnel[i][1];
        //Calculate the voltage difference
        //Temporarily the voltage difference will be the current value, unless the element resistance is not 1
        I = abs(Voltage_difference(P1, P2, DEA->LM_matrix, DEA->voltages));
        //Check the resistance flag to use the appropiate resistance
        //If R_flag is 0, the voltage difference is the current, so I remains the same
        if (R_flag) {
            //If R_flag is 1, then use the actual resistance
            
            //Calculate the Resistance of the tunnel
            Re = DEA->Calculate_resistance_tunnel(radii, electric_param, point_list[P1], point_list[P2], 0.34);
            //Calculate current
            I = I/Re;
        }
        //Add current to vector of all currents
        currents.push_back(I);
    }
    
    //Initialize the vectors of element currents and resistances for the GNPs
    current_gnp.assign(cluster_gch.size(), empty_double);
    resistance_gnp.assign(cluster_gch.size(), empty_double);
    
    //Calculate currents from GNPs
    for (int i = 0; i < (int)cluster_gch.size(); i++) {
        int hyb = cluster_gch[i];
        //Scan all edges of the triangulation
        for (int j = 0; j < (int)hybrid_particles[hyb].triangulation.size(); j++) {
            //Calculate the voltage difference
            P1 = hybrid_particles[hyb].triangulation[j][0];
            P2 = hybrid_particles[hyb].triangulation[j][1];
            //Temporarily the voltage difference will be the current value, unless the element resistance is not 1
            I = abs(Voltage_difference(P1, P2, DEA->LM_matrix, DEA->voltages));
            //Check the resistance flag to use the appropiate resistance
            //If R_flag is 0, the voltage difference is the current, so I remains the same
            if (R_flag) {
                //If R_flag is 1, then use the actual resistance
                
                //Calculate resistance of the triangulatin edge
                Re = DEA->Calculate_resistance_gnp(point_list[P1], point_list[P2], radii[point_list[P1].flag], radii[point_list[P2].flag], hybrid_particles[hyb], electric_param);
                //Calculate current
                I = I/Re;
                //Add resistance to the vector of element resistances
                resistance_gnp[i].push_back(Re);
            } 
            //Add current to vector of all currents
            currents.push_back(I);
            //Add current to the vector of element currents
            current_gnp[i].push_back(I);
        }
    }
    
    //Sort currents
    sort(currents.begin(),currents.end());
    
    //The error cutoff seems to work well with a drop in 9 orders of magnitude of the current. So that is how the cutoff is set.
    //This idea comes from Li and Chou's paper of the DEA in which using a voltage of 1V, a drop in 9 orders of magnitude
    //in the current gave good results.
    zero_cutoff = currents.back()*1e-9;
    
    Printer *P = new Printer;
    if (R_flag) {
        P->Print_1d_vec(currents, "currents_R.txt");
    } else {
        P->Print_1d_vec(currents, "currents.txt");
    }
    delete P;
    
    return zero_cutoff;
}
//This function calculates the voltage difference between two nodes
double Backbone_Network::Voltage_difference(const long int &P1, const long int &P2, const vector<int> &LM_matrix, const vector<double> &voltages)
{
    //Get the node numbers
    int node1 = LM_matrix[P1];
    int node2 = LM_matrix[P2];
    //Calculate the voltage difference
    return voltages[node2] - voltages[node1];
}
//This function scans all CNTs and calculates the currents again, then it decides which CNTs are part of the backbone and which CNTs are not percolated
int Backbone_Network::Find_dead_branches_simplified(vector<vector<long int> > &elements, const vector<int> &cluster, const double &zero_cutoff)
{
    //This varibale is used to initialize the vectors below
    vector<long int> empty;
    //This variable will store the begining and end point of the non-conducting branches of each CNT in the cluster
    //Then each vector withing the vector can only have size 2 or 4
    dead_indices.assign(cluster.size(), empty);
    //This variable will store the begining and end point of the conducting segment of each CNT in the cluster
    //Then each vector withing the vector can only have size 2
    percolated_indices.assign(cluster.size(), empty);
    
    //Variables
    int CNT;
    //Variables used to calculate the current
    double I;
    
    //Scan all the CNTs that belong to the current cluster
    for (long int i = 0; i < (long int)cluster.size(); i++) {
        CNT = cluster[i];
        //The dead branches are in th extremes of a CNT, they cannot be in the middle.
        
        //Vector to store the elements that have a current above the cutoff
        vector<double> element_index;
        //vectors cluster and current_e have the same size
        for (int j = 0; j < (int)current_e[i].size(); j++) {
            //Get the element current
            I = current_e[i][j];
            if (I > zero_cutoff) {
                //If the current is above the zero cutoff, add the current j-index to the vector of element_index
                element_index.push_back(j);
            }
        }
        
        //Variable for the j-index
        int jj;
        //----------- Check for a conducting segment -----------
        if (element_index.size()) {
            //----------- Add conducting segment -----------
            //The first point of the conducting segment is given by the first j-index in element_index
            jj = element_index.front();
            percolated_indices[i].push_back(elements[CNT][jj]);
            //The second point of the conducting segment is given by the last j-index in element_index
            jj = element_index.back();
            //An element has two indices, and for the conducting segment, I need the second index of the last conducting element
            //Thus in this case I use jj+1
            percolated_indices[i].push_back(elements[CNT][jj+1]);
            
            //----------- Check for a dead branch on the front -----------
            //If the first j-index in element_index is NOT 0, then there is a dead brach at the front of the CNT
            if (element_index.front() != 0) {
                //The first index of the dead branch is the front of the CNT
                dead_indices[i].push_back(elements[CNT].front());
                //The second index is given by the first j-index in element_index
                jj = element_index.front();
                dead_indices[i].push_back(elements[CNT][jj]);
            }
            
            //----------- Check for a dead branch on the back -----------
            //If the last j-index in element_index is NOT current_e[i].size()-1, then there is a dead brach at the back of the CNT
            if (element_index.back() != (current_e[i].size()-1) ) {
                //The first index of the dead branch is given by the last j-index in element_index
                jj = element_index.back();
                //The last j-index is the first index of the last conducting element
                //An element has two indices, and for the dead branch, I need the second index of the last conducting element
                //Thus in this case I use jj+1
                dead_indices[i].push_back(elements[CNT][jj+1]);
                //The last index of the dead branch is the back of the CNT
                dead_indices[i].push_back(elements[CNT].back());
            }
            
        } else {
            //----------- Whole CNT is a dead branch -----------
            //If there is no conducting segment, the whole CNT is a dead brach
            //Thus, add the first and last points of the CNT to the dead indices
            dead_indices[i].push_back(elements[CNT].front());
            dead_indices[i].push_back(elements[CNT].back());
        }
        
        //Check that the vectors have the correct size. percolated_indices can only have size 0 or 2
        if ( (percolated_indices[i].size() != 2) && (percolated_indices[i].size() != 0) ) {
            hout << "Error in Find_dead_branches. The vector percolated_indices["<<i<<"] has size " << percolated_indices[i].size();
            hout << " but it can only have size 0 or 2." << endl;
            hout << "\tThe vector dead_indices["<<i<<"] has size " << dead_indices[i].size()<<endl;
            hout << "\tCNTs in cluster: " << cluster.size() << ", current CNT: "<<CNT;
            hout << ", elements[CNT].size()=" << elements[CNT].size() << endl;
            hout << "\tLast calculated current: "<<I<<", zero cutoff: "<<zero_cutoff<<endl;
            return 0;
        }
        //Check that the vectors have the correct size. dead_indices can only have size 0, 2 or 4
        if ( (dead_indices[i].size() != 4) && (dead_indices[i].size() != 2) && (dead_indices[i].size() != 0) ) {
            hout << "Error in Find_dead_branches. The vector dead_indices["<<i<<"] has size " << dead_indices[i].size();
            hout << " but it can only have size 0, 2 or 4." <<endl;
            hout << "\tThe vector percolated_indices["<<i<<"] has size " << percolated_indices[i].size()<<endl;
            hout << "\tCNTs in cluster: " << cluster.size() << ", current CNT: "<<CNT;
            hout << ", elements[CNT].size()=" << elements[CNT].size() << endl;
            hout << "\tLast calculated current: "<<I<<", zero cutoff: "<<zero_cutoff<<endl;
            return 0;
        }
        
    }
    
    //Printer *P = new Printer;
    //P->Print_2d_vec(percolated_indices, "percolated_indices_simplified.txt");
    //P->Print_2d_vec(dead_indices, "dead_indices_simplified.txt");
    //delete P;
    
    return 1;
}
//
int Backbone_Network::Find_dead_gnps_simplified(const double &zero_cutoff, vector<int> &cluster_gch, vector<int> &gnp_dead_indices, vector<int> &gnp_indices)
{
    //Variables
    int hyb;
    
    //Scan every GNP in the cluster
    //Scan backwards to avoid conflicts with the indices when deleting an element
    for (int i = 0; i < (int)cluster_gch.size() ; i++) {
        //Flag to determine if the GNP is part of the backbone or not
        int percolated = 0;
        
        //Current hybrid particle
        hyb = cluster_gch[i];
        
        //hout << "Top triangulation " << i << endl;
        //Scan every current in the resistors coming from the triangulation
        //current_gnp and cluster_gch have the same size
        for (int j = 0; j < (int)current_gnp[i].size(); j++) {
            //If at least one current (current_gnp[i][j]) is above the cutoff then the GNP is part of the backbone
            if (current_gnp[i][j] > zero_cutoff) {
                //Set the percolated flag to 1
                percolated = 1;
                
                //Add the GNP to the vector gnp_indices
                gnp_indices.push_back(hyb);
                //Add the GNP to the local vector of percolated GNPs
                percolated_gnps.push_back(hyb);
                //There is no need to calculate more currents so break the loop
                break;
            }
        }
        
        //Check if the percolated flag was set to 1
        if (!percolated) {
            //If not set to one, then add the GNP to the cluster of dead GNPs.
            gnp_dead_indices.push_back(hyb);            
        }
    }
    
    //Printer *P = new Printer;
    //P->Print_1d_vec(percolated_gnps, "percolated_gnps_simplified.txt");
    //delete P;
    
    return 1;
}
//This vector calculates the lengths of the CNTs that form part of the backbone and the dead branches
//It is assumed that the vectors are initialized and the quantities calculated for the current cluster are added
//To the corresponding element of the vector
//In both vectors, except for the last one, each element corresponds to a families in the following order:
//[0] - X
//[1] - Y
//[2] - Z
//[3] - XY
//[4] - XZ
//[5] - YZ
//[6] - XYZ
//[7] - NP (only families_lengths)
int Backbone_Network::Calculate_lengths(const int &family, const vector<Point_3D> &points_in, vector<double> &families_lengths, vector<double> &branches_lengths)
{
    //Variables to store lengths
    double backbone = 0, branches = 0;
    
    //Scan the percolated_indices vector to calculate the length of the backbone
    //Since the branches_lengths has the same length as percolated_indices, both vectors can be scanned in the same loop
    for (int i = 0; i < (int)percolated_indices.size(); i++) {
        //Check if there is a conducting segment, if there is add the corresponding length
        if (percolated_indices[i].size()){
            backbone = backbone + Segment_length(percolated_indices[i][0], percolated_indices[i][1], points_in);
        }
        //Check if there is one or two non-conducting segments
        if (dead_indices[i].size() == 2) {
            branches = branches + Segment_length(dead_indices[i][0], dead_indices[i][1], points_in);
        } else if (dead_indices[i].size() == 4){
            branches = branches + Segment_length(dead_indices[i][0], dead_indices[i][1], points_in);
            branches = branches + Segment_length(dead_indices[i][2], dead_indices[i][3], points_in);
        }
    }
    
    //Add the calculated lengths to the corresponding family
    families_lengths[family] = families_lengths[family] + backbone;
    branches_lengths[family] = branches_lengths[family] + branches;
    
    return 1;
}

//This function calculates the length of a segement of a CNT given by a sequence of points from (global number) index1 to index2
double Backbone_Network::Segment_length(long int index1, long int index2, const vector<Point_3D> &points_in)
{
    double length = 0;
    for (int long j = index1; j < index2; j++) {
        length = length + points_in[j].distance_to(points_in[j+1]);
    }
    return length;
}

//
void Backbone_Network::Add_indices_to_global_vectors(const int &family, vector<vector<long int> > &all_dead_indices, vector<vector<long int> > &all_indices)
{
    //Add indices for percolated segments
    for (int i = 0; i < (int)percolated_indices.size(); i++) {
        for (int j = 0; j < (int)percolated_indices[i].size(); j++) {
            all_indices[family].push_back(percolated_indices[i][j]);
        }
    }
    
    //Add indices for dead branches
    for (int i = 0; i < (int)dead_indices.size(); i++) {
        for (int j = 0; j < (int)dead_indices[i].size(); j++) {
            all_dead_indices[family].push_back(dead_indices[i][j]);
        }
    }
    
}