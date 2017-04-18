//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Electrical_analysis.cpp
//OBJECTIVE:	Extract the backbone and calculate the electrical resistance
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Electrical_analysis.h"

int Electrical_analysis::Perform_analysis_on_clusters(const int &iteration, const vector<int> &family, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &point_list_gnp, const struct Geom_RVE &geom_rve, const struct Electric_para &electric_param, const struct Cutoff_dist &cutoffs, vector<GCH> &hybrid_particles, vector<double> &families_lengths, vector<double> &branches_lengths, vector<vector<long int> > &all_dead_indices, vector<vector<long int> > &all_indices, vector<vector<int> > &gnp_dead_indices, vector<vector<int> > &gnp_indices)
{
    //Time variables
    time_t ct0, ct1;
    
    //Vector of parallel resistors
    //Each cluster will contribute with a resistor to each direction in which it percolates
    //So each cluster adds a parallel resistor on each percolated direction
    vector<double> empty_double;
    vector<vector<double> > paralel_resistors(3, empty_double);
    
    //Get the number of clusters
    int n_clusters = 0;
    //hout<<"clusters_cnt.size()="<<clusters_cnt.size()<<endl;
    //hout<<"clusters_gch()="<<clusters_gch.size()<<endl;
    if (HoKo->clusters_cnt.size()) {
        n_clusters = (int)HoKo->clusters_cnt.size();
        
    } else if (HoKo->clusters_gch.size()) {
        n_clusters = (int)HoKo->clusters_gch.size();
    }
    
    //Scan every percolated cluster
    for (int j = 0; j < n_clusters; j++) {
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Current iteration
        hout << "=============================" <<endl;
        hout << "Cluster " << j+1 << " of " << n_clusters <<", family " << family[j] << endl;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Direct Electrifying algorithm
        Direct_Electrifying *DEA = new Direct_Electrifying;
        
        //Resitor flag, set to 0 to use unit resistors
        int R_flag = 0;
        
        //Keep a copy
        vector<vector<long int> > contacts_initial(HoKo->contacts_point);
        
        ct0 = time(NULL);
        if (!DEA->Calculate_voltage_field(family[j], j, R_flag, HoKo, Cutwins, structure, point_list, radii, structure_gnp, point_list_gnp, electric_param, cutoffs, hybrid_particles)) {
            hout << "Error in Perform_analysis_on_cluster when calling DEA->Calculate_voltage_field using unit resistances" << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Calculate voltage field time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine the backbone and dead branches
        Backbone_Network *Backbonet = new Backbone_Network;
        ct0 = time(NULL);
        if (!Backbonet->Determine_backbone_network(family[j], j, R_flag, 1, DEA, HoKo, electric_param, cutoffs, structure, point_list, radii, structure_gnp, point_list_gnp, hybrid_particles, families_lengths, branches_lengths, all_dead_indices, all_indices, gnp_dead_indices, gnp_indices)) {
            hout << "Error in Perform_analysis_on_cluster when calling Determine_backbone_network" << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Determine backbone network time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Set now the R_flag to 1 to indicate that actual resistances will be used
        R_flag = 1;
        
        //Calculate the electrical resistances along each direction for current clusters
        if (!Electrical_analysis_along_each_percolated_direction(R_flag, j, family, HoKo, Cutwins, Backbonet, electric_param, cutoffs, structure, point_list, radii, contacts_initial, structure_gnp, point_list_gnp, hybrid_particles, paralel_resistors)) {
            hout << "Error in Perform_analysis_on_cluster when calling Electrical_analysis_along_each_percolated_direction" << endl;
            return 0;
        }
        
        //Delete objects to free memory
        delete DEA;
        delete Backbonet;
    }
    
    //Calculate the matrix resistances on each direction
    vector<double> matrix_resistances;
    if (!Calculate_matrix_resistances(electric_param.resistivity_matrix, geom_rve, matrix_resistances)) {
        hout << "Error in Perform_analysis_on_cluster when calling Calculate_matrix_resistances" << endl;
        return 0;
    }
    
    //Calculate the resistances and resistivities along each direction
    vector<double> resistivities;
    if (!Calculate_resistances_and_resistivities(geom_rve, matrix_resistances, paralel_resistors, resistors, resistivities)) {
        hout << "Error in Perform_analysis_on_cluster when calling Calculate_resistances_and_resistivities" << endl;
        return 0;
    }
    
    //Append resistors to a file
    Printer *P = new Printer;
    P->Append_1d_vec(resistors, "resistors.txt");
    P->Append_1d_vec(resistivities, "resistivities.txt");
    delete P;
    
    return 1;
}
//
int Electrical_analysis::Electrical_analysis_along_each_percolated_direction (const int &R_flag, const int &n_cluser, const vector<int> &family, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, Backbone_Network *Backbonet, const struct Electric_para &electric_param, const struct Cutoff_dist &cutoffs, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<long int> > &contacts_initial, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &point_list_gnp, vector<GCH> &hybrid_particles, vector<vector<double> > &paralel_resistors)
{
    //Time variables
    time_t ct0, ct1;
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Get the vector of directions
    ct0 = time(NULL);
    vector<int> directions;
    if (!Vector_of_directions(family[n_cluser], directions)) {
        hout << "Error in Perform_analysis_on_cluster when calling Vector_of_directions" << endl;
        return 0;
    }
    ct1 = time(NULL);
    hout << "Vector_of_directions time: "<<(int)(ct1-ct0)<<" secs."<<endl;
    
    //Set now the R_flag to 1 to indicate that actual resistances will be used
    //R_flag = 1;
    
    //Calculate the electrical resistance per direction
    for (int k = 0; k < (int)directions.size(); k++) {
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Direct Electrifying algorithm to calculate electrical resistance
        Direct_Electrifying *DEA_Re = new Direct_Electrifying;
        
        //Generate a structure vector
        vector<long int> empty_long;
        vector<vector<long int> > backbone_structure(structure.size(), empty_long);
        vector<int> backbone_cnts;
        //First check if there are CNTs in the cluster
        if (HoKo->clusters_cnt.size() && HoKo->clusters_cnt[n_cluser].size()) {
            if (!Convert_index_to_structure(HoKo->clusters_cnt[n_cluser], Backbonet->percolated_indices, backbone_structure, backbone_cnts)) {
                hout << "Error in Perform_analysis_on_cluster when calling Convert_index_to_structure" << endl;
                return 0;
            }
            hout << "Convert_index_to_structure" << endl;
        }
        
        
        //Update the hybrid particles that belong to the backbone
        if (!Update_hybrids(Backbonet->percolated_gnps, structure, backbone_structure, hybrid_particles)) {
            hout << "Error in Perform_analysis_on_cluster when calling Update_hybrids" << endl;
            return 0;
        }
        //hout << "Update_hybrids" << endl;
        
        //Temporary HoKo and Cutwins
        Hoshen_Kopelman *HoKo_Re = new Hoshen_Kopelman;
        Cutoff_Wins *Cutwins_Re = new Cutoff_Wins;
        if (!Update_vectors_for_hoko_cutwins((int)structure.size(), (int)structure_gnp.size(), HoKo, Cutwins, HoKo_Re, Cutwins_Re, contacts_initial, backbone_structure, backbone_cnts, Backbonet->percolated_gnps)) {
            hout << "Error in Perform_analysis_on_cluster when calling Update_vectors_for_hoko_cutwins" << endl;
            return 0;
        }
        //hout << "Update_vectors_for_hoko_cutwins" << endl;
        
        //Run a new DEA to obtain the new voltage field in the backbone using the actual resistances
        ct0 = time(NULL);
        if (!DEA_Re->Calculate_voltage_field(directions[k], 0, R_flag, HoKo_Re, Cutwins_Re, backbone_structure, point_list, radii, structure_gnp, point_list_gnp, electric_param, cutoffs, hybrid_particles)) {
            hout << "Error in Perform_analysis_on_cluster when calling DEA_Re->Calculate_voltage_field using actual resistances" << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Calculate voltage field on backbone network time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //With the new voltage field calculate the current going through a face and calculate the resistance along that direction
        if (!Calculate_parallel_resistor(directions[k], DEA_Re, point_list, radii, HoKo->clusters_cnt, Cutwins_Re->boundary_cnt, point_list_gnp, hybrid_particles, HoKo->clusters_gch, Cutwins_Re->boundary_gnp, Cutwins_Re->boundary_flags_gnp, electric_param, paralel_resistors)) {
            hout << "Error in Perform_analysis_on_cluster when calling Calculate_parallel_resistor" << endl;
            return 0;
        }
        //hout << "Calculate_parallel_resistor" << endl;
        
        //Delete objects to free memory
        delete DEA_Re;
        delete HoKo_Re;
        delete Cutwins_Re;
        
    }//*/
    
    return 1;
}
//This function generates a vector with the different directions in which the electrical resistance needs to be calculated
//e.g. if the family is 3 (i.e. XY), a vector with elements {0,1} is generated
int Electrical_analysis::Vector_of_directions(const int &family, vector<int> &directions)
{
    if (family == 6) {
        //If the family is 6 (XYZ); then this the resistance needs to be calculated in the three directions
        directions.assign(3, 0);
        directions[1] = 1;
        directions[2] = 2;
    } else if (family == 5) {
        //If the family is 5 (YZ); then this the resistance needs to be calculated in two directions: 1 and 2
        directions.push_back(1);
        directions.push_back(2);
    } else if (family == 4) {
        //If the family is 4 (XZ); then this the resistance needs to be calculated in two directions: 0 and 2
        directions.push_back(0);
        directions.push_back(2);
    } else if (family == 3) {
        //If the family is 5 (XY); then this the resistance needs to be calculated in two directions: 0 and 1
        directions.push_back(0);
        directions.push_back(1);
    } else if (family <= 2) {
        //If the family is 0, 1 or 2; then this is the only direction in which resistance needs to be calculated
        directions.push_back(family);
    } else {
        hout << "Invalid family: " << family << endl;
        return 0;
    }
    
    return 1;
}
//This function converts the data type index into data type structure
//The structure that is generated must had been initialized with the size of the original structure
int Electrical_analysis::Convert_index_to_structure(const vector<int> &cluster, const vector<vector<long int> > &indices, vector<vector<long int> > &backbone_structure, vector<int> &backbone_cnts)
{
    //The branches are given in pairs
    for (int i = 0; i < (int)indices.size(); i++) {
        //Check if the current index vector has any elements
        if (indices[i].size()) {
            //If the index vector has elements, then there is a conducting segment
            
            //CNT number of the conducting segment
            int CNT = cluster[i];
            //Add CNT to the cluster of backbone CNTs
            backbone_cnts.push_back(CNT);
            
            //Generate the backbone struture
            for (long int j = indices[i][0]; j <= indices[i][1]; j++) {
                backbone_structure[CNT].push_back(j);
            }
            
        }
    }
    return 1;
}
//This function updates the hybrid particles that are part of the backbone:
//     - Update the vector of CNTs attached to the GNP surface
//     - Clear the triangulation vectors
int Electrical_analysis::Update_hybrids(const vector<int> &cluster_gch, const vector<vector<long int> > &structure, const vector<vector<long int> > &backbone_structure, vector<GCH> &hybrid_particles)
{
    //Scan all hybrids in the cluster
    for (int i = 0; i < (int)cluster_gch.size(); i++) {
        //Current hybrid
        int hyb = cluster_gch[i];
        
        //---------------- TOP CNTs
        //Temporary vector initialized with the same elements as the cnts in the top surface of the GNP
        vector<int> top_tmp(hybrid_particles[hyb].cnts_top);
        //Clear the vector of CNTs on the top surface, so now CNTs are in the temporaty vector
        hybrid_particles[hyb].cnts_top.clear();
        //Scan all CNTs in the top surface of current hybrid
        for (int j = 0; j < (int) top_tmp.size(); j++) {
            //Current CNT
            int CNT = top_tmp[j];
            //Check if the initial points of the original structure and the backbone_structure are the same
            if (backbone_structure[CNT].size() && structure[CNT].front() == backbone_structure[CNT].front()) {
                //If initial points are the same, then add the CNT to the vetor of CNTs on the top surface
                hybrid_particles[hyb].cnts_top.push_back(CNT);
            }
        }
        
        //---------------- BOTTOM CNTs
        //Temporary vector initialized with the same elements as the cnts in the bottom surface of the GNP
        vector<int> bottom_tmp(hybrid_particles[hyb].cnts_bottom);
        //Clear the vector of CNTs on the bottom surface, so now CNTs are in the temporaty vector
        hybrid_particles[hyb].cnts_bottom.clear();
        //Scan all CNTs in the bottom surface of current hybrid
        for (int j = 0; j < (int) bottom_tmp.size(); j++) {
            //Current CNT
            int CNT = bottom_tmp[j];
            //Check if the initial points of the original structure and the backbone_structure are the same
            if (backbone_structure[CNT].size() && structure[CNT].front() == backbone_structure[CNT].front()) {
                //If initial points are the same, then add the CNT to the vetor of CNTs on the bottom surface
                hybrid_particles[hyb].cnts_bottom.push_back(CNT);
            }
        }

    }
    
    //When more than one percolated clusters are present, the triangulation for a GNP of one cluster
    //might be used when calculating the current of another cluster
    //Thus, clear all triangulations
    for (int i = 0; i < (int)hybrid_particles.size(); i++) {
        
        hybrid_particles[i].triangulation.clear();
        hybrid_particles[i].triangulation_flags.clear();
        
    }
    
    return 1;
}
//
int Electrical_analysis::Update_vectors_for_hoko_cutwins(const int &n_cnts, const int &n_gnps, Hoshen_Kopelman *HoKo, Cutoff_Wins *Cutwins, Hoshen_Kopelman *HoKo_Re, Cutoff_Wins *Cutwins_Re, const vector<vector<long int> > &contacts_initial, const vector<vector<long int> > &structure, const vector<int> &backbone_cnts, const vector<int> &percolated_gnps)
{
    //Update vectors of temporary HoKo and Cutwins
    HoKo_Re->contacts_point = contacts_initial;
    if (backbone_cnts.size())
        HoKo_Re->clusters_cnt.push_back(backbone_cnts);
    if (percolated_gnps.size())
        HoKo_Re->clusters_gch.push_back(percolated_gnps);
    Cutwins_Re->boundary_cnt = Cutwins->boundary_cnt;
    Cutwins_Re->boundary_gnp = Cutwins->boundary_gnp;
    Cutwins_Re->boundary_flags_cnt = Cutwins->boundary_flags_cnt;
    Cutwins_Re->boundary_flags_gnp = Cutwins->boundary_flags_gnp;
    
    //Vector flags to determine CNTs and GNPs in the cluster
    vector<short int> cnt_cluster_flags(n_cnts, 0), gnp_cluster_flags(n_gnps, 0);
    
    //Fill flags
    for (int i = 0; i < (int)backbone_cnts.size(); i++) {
        
        //Current CNT
        int CNT = backbone_cnts[i];
        
        //Set flag
        cnt_cluster_flags[CNT] = 1;
    }
    for (int i = 0; i < (int)percolated_gnps.size(); i++) {
        
        //Current GNP
        int GNP = percolated_gnps[i];
        
        //Set flag
        gnp_cluster_flags[GNP] = 1;
    }
    
    //Update gnp contact vectors
    for (int i = 0; i < (int)HoKo->gnp_contacts.size(); i++) {
        
        //Get the GNP numbers
        int GNP1 = HoKo->gnp_contacts[i].particle1;
        int GNP2 = HoKo->gnp_contacts[i].particle2;
        
        //Check if both GNPs are part of the cluster
        if (gnp_cluster_flags[GNP1] && gnp_cluster_flags[GNP2]) {
            
            //If both GNPs are in the cluster, then add the contact
            HoKo_Re->gnp_contacts.push_back(HoKo->gnp_contacts[i]);
        }
    }
    
    //Update mixed contact vectors
    for (int i = 0; i < (int)HoKo->mixed_contacts.size(); i++) {
        
        //Get the CNT and GNP numbers
        int CNT = HoKo->mixed_contacts[i].particle1;
        int GNP = HoKo->mixed_contacts[i].particle2;
        
        //Get the CNT point
        long int P_CNT = HoKo->mixed_contacts[i].point1;
        
        //Check if both GNP and CNT are part of the cluster and if the CNT point is still in the CNT
        if (cnt_cluster_flags[CNT] && gnp_cluster_flags[GNP] && P_CNT >= structure[CNT].front() && P_CNT <= structure[CNT].back()) {
            
            //If both GNPs are in the cluster, then add the contact
            HoKo_Re->mixed_contacts.push_back(HoKo->mixed_contacts[i]);
        }
    }
    
    /*/
    Printer *P =  new Printer;
    P->Print_1d_vec(HoKo->gnp_contacts, "gnp_contacts_unit.txt");
    P->Print_1d_vec(HoKo_Re->gnp_contacts, "gnp_contacts_R.txt");
    delete P;//*/
    
    return 1;
}
//Calculate the resistance along a direction
//Note: direction has values 0, 1 or 2 only to represent directions X, Y and Z, respectively
//The indices of boundary_cnt are as follows:
//0,1 for X boundaries
//2,3 for Y boundaries
//4,5 for Z boundaries
//Thus 2*direction will be 0, 2 or 4, i.e. the first boundary on each direction
//Then, 2*direction+1 will be 1, 3 or 5, i.e. the second boundary on each direction
int Electrical_analysis::Calculate_parallel_resistor(const int &direction, Direct_Electrifying * DEA, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<vector<int> > &clusters_cnt, const vector<vector<int> > &boundary_cnt, const vector<Point_3D> &point_list_gnp, const vector<GCH> &hybrid_particles, const vector<vector<int> > &clusters_gnp, const vector<vector<int> > &boundary_gnp, const vector<vector<short int> > &boundary_flags_gnp, const struct Electric_para &electric_param, vector<vector<double> > &paralel_resistors)
{
    //Currents from CNTs
    double I_total_cnt = 0;
    double I_total_cnt_check = 0;
    
    //Currents from GNPS
    double I_total_gnp = 0;
    double I_total_gnp_check = 0;

    
    //---------------- Currents through CNTs
    //Check if there are CNTs in the cluster
    if (clusters_cnt.size()) {
        
        //When there are no GNPs, if any of the two boundary vectors has no CNTs there is an error as there cannot be percolation in this direction
        if (!clusters_gnp.size() && ( !boundary_cnt[2*direction].size() || !boundary_cnt[2*direction+1].size() ) ) {
            hout << "One boundary vector along direction " << direction << " is empty, so there cannot be percolation along that direction."<< endl;
            hout << "\t boundary_cnt[" << 2*direction << "].size() = " << boundary_cnt[2*direction].size() << endl;
            hout << "\t boundary_cnt[" << 2*direction+1 << "].size() = " << boundary_cnt[2*direction+1].size() << endl;
            return 0;
        }
        
        long int P1, P2;
        
        //Scan all CNTs at the boundary
        for (int i = 0; i < (int)boundary_cnt[2*direction].size(); i++) {
            //Current CNT
            int CNT = boundary_cnt[2*direction][i];
            
            //Some CNTs on the boundary might not be part of the backbone or the geometric cluster
            //First check if there are any elements on the CNT, if there are no elements there, skip the CNT
            if (DEA->elements[CNT].size()) {
                //Check if the front and/or back of the CNT are in contact with the boundary
                
                //Get the points of the element at the front of the CNT
                P1 = DEA->elements[CNT].front();
                P2 = DEA->elements[CNT][1];
                //Add the current of the element at the front of the CNT
                I_total_cnt = I_total_cnt + Current_of_element_in_boundary(P1, P2, radii[CNT], DEA, electric_param, point_list);
                
                //Get the points of the element at the back of the CNT
                P1 = DEA->elements[CNT].back();
                int size = (int)DEA->elements[CNT].size();
                P2 = DEA->elements[CNT][size-2];
                //Add the current of the element at the front of the CNT
                I_total_cnt = I_total_cnt + Current_of_element_in_boundary(P1, P2, radii[CNT], DEA, electric_param, point_list);
                
            }
        }
        
        //=========================== CURRENT Check
        hout <<"//=========================== CURRENT Check"<<endl;
        //Scan all CNTs at the boundary
        for (int i = 0; i < (int)boundary_cnt[2*direction+1].size(); i++) {
            //Current CNT
            int CNT = boundary_cnt[2*direction+1][i];
            
            //Some CNTs on the boundary might not be part of the backbone or the geometric cluster
            //First check if there are any elements on the CNT, if there are no elements there, skip the CNT
            if (DEA->elements[CNT].size()) {
                //Check if the front and/or back of the CNT are in contact with the boundary
                
                //Get the points of the element at the front of the CNT
                P1 = DEA->elements[CNT].front();
                P2 = DEA->elements[CNT][1];
                //Add the current of the element at the front of the CNT
                I_total_cnt_check = I_total_cnt_check + Current_of_element_in_boundary(P1, P2, radii[CNT], DEA, electric_param, point_list);
                
                //Get the points of the element at the back of the CNT
                P1 = DEA->elements[CNT].back();
                int size = (int)DEA->elements[CNT].size();
                P2 = DEA->elements[CNT][size-2];
                //Add the current of the element at the front of the CNT
                I_total_cnt_check = I_total_cnt_check + Current_of_element_in_boundary(P1, P2, radii[CNT], DEA, electric_param, point_list);
                
            }
        }
        
        hout << "I_total_cnt="<<I_total_cnt<<" direction="<<direction<<endl;
        hout << "I_total_cnt_check="<<I_total_cnt_check<<" direction="<<direction<<endl;
        
    }
    
    //---------------- Currents through GNPs
    //Check if there are GNPs in the cluster
    if (clusters_gnp.size()) {
        
        //When there are no CNTs, if any of the two boundary vectors has no GNPs there is an error as there cannot be percolation in this direction
        if ( !clusters_cnt.size() && ( !boundary_gnp[2*direction].size() || !boundary_gnp[2*direction+1].size() )) {
            hout << "One boundary vector along direction " << direction << " is empty, so there cannot be percolation along that direction."<< endl;
            hout << "\t boundary_gnp[" << 2*direction << "].size() = " << boundary_gnp[2*direction].size() << endl;
            hout << "\t boundary_gnp[" << 2*direction+1 << "].size() = " << boundary_gnp[2*direction+1].size() << endl;
            return 0;
        }
        
        //hout << "Current GNP" << endl;
        //Scan all GNPs at the boundary
        for (int i = 0; i < (int)boundary_gnp[2*direction].size(); i++) {
            
            //Current GNP
            int GNP = boundary_gnp[2*direction][i];
            
            //Some GNPs on the boundary might not be part of the backbone or the geometric cluster
            //First check if there are any triangulation edges on the GNP, if there are no edges there, skip the GNP
            if (hybrid_particles[GNP].triangulation.size()) {
                
                hout << "GNP="<<GNP<<" dir="<<direction<<" boundary="<<2*direction<<endl;
                
                //Add the currents from the triangulation edges at the boundary
                I_total_gnp = I_total_gnp + Current_of_edges_in_boundary(0, DEA, electric_param, point_list, radii, point_list_gnp, hybrid_particles[GNP]);
                
            }
        }
        
        
        //=========================== CURRENT Check GNP
        hout <<"//=========================== CURRENT Check GNP"<<endl;
        //Scan all GNPs at the boundary
        for (int i = 0; i < (int)boundary_gnp[2*direction+1].size(); i++) {
            
            //Current GNP
            int GNP = boundary_gnp[2*direction+1][i];
            
            //Some GNPs on the boundary might not be part of the backbone or the geometric cluster
            //First check if there are any triangulation edges on the GNP, if there are no edges there, skip the GNP
            if (hybrid_particles[GNP].triangulation.size()) {
                
                hout << "GNP="<<GNP<<" dir="<<direction<<" boundary="<<2*direction+1<<endl;
                
                //Add the currents from the triangulation edges at the boundary
                I_total_gnp_check = I_total_gnp_check + Current_of_edges_in_boundary(1, DEA, electric_param, point_list, radii, point_list_gnp, hybrid_particles[GNP]);
            }
            
        }
        
        hout << "I_total_gnp="<<I_total_gnp<<" direction="<<direction<<endl;
        hout << "I_total_gnp_check="<<I_total_gnp_check<<" direction="<<direction<<endl;
    }
    

    //Calculate total currents
    double I_total = I_total_cnt + I_total_gnp;
    double I_total_check = I_total_cnt_check + I_total_gnp_check;
    hout << "I_total="<<I_total<<" direction="<<direction<<endl;
    hout << "I_total_check="<<I_total_check<<" direction="<<direction<<endl;
    
    //Calculate total resistance
    //To calculate the resistance in a single direction, the DEA was also run using a single direction
    //thus, there are only two boundary conditions and the voltage difference is equal to the input voltage
    double R_total = electric_param.applied_voltage/I_total;
    
    //Add resistor to vector of parallel resistors
    paralel_resistors[direction].push_back(R_total);
    
    return 1;
}
//This function checks if an element is at a boundary, and if so it calculates the current
//It is assumed that P1 is the point that is at either the front or the back of the CNT, only these points can be in contact with the boundary
double Electrical_analysis::Current_of_element_in_boundary(const long int &P1, const long int &P2, const double &radius, Direct_Electrifying *DEA, const struct Electric_para &electric_param, const vector<Point_3D> &point_list)
{
    //Check where is P1
    //hout <<"P1="<<P1<<" node1="<<DEA->LM_matrix[P1]<<" ("<<point_list[P1].x<<", "<<point_list[P1].y<<", ";
    //hout <<point_list[P1].z<<")"<<endl;
    //If P1 is node 0 or 1, then the current is calculated
    //This means it is on a valid boundary
    int node1 = DEA->LM_matrix[P1];
    int node2 = DEA->LM_matrix[P2];
    if (node1 <= 1) {
        
        hout <<"P1="<<P1<<" node1="<<DEA->LM_matrix[P1]<<" ("<<point_list[P1].x<<", "<<point_list[P1].y<<", "<<point_list[P1].z<<")"<<endl;//Get the node numbers
        hout << "P2="<<P2<<" node2="<<node2<<endl;
        //Calculate voltage difference on the element, node1 is at the boundary with voltage 0
        //so the voltage drop is from node2 to node1
        double V = DEA->voltages[node2] - DEA->voltages[node1];
        hout << "V1="<<DEA->voltages[node1]<<" V2="<<DEA->voltages[node2]<<"\nDV=" << V;
        //Calculate resistance of the element
        double Re;
        if (P1 > P2)
            //Check if the resistor is at the back of the CNT, since in that case the calculated resistance will be zero;
            //this because if P1 > P2, the function that calculates the resistance does not find points after P1
            //When the resistor is at the back of the elemnt P1 > P2, so in that case invert the point numbers
            Re = DEA->Calculate_resistance_cnt(point_list, P2, P1, radius, electric_param.resistivity_CF);
        else
            Re = DEA->Calculate_resistance_cnt(point_list, P1, P2, radius, electric_param.resistivity_CF);
        //Calculate current and add it to the total current
        hout << " Re=" << Re << " I=" << V/Re << endl;
        return V/Re;
    }
    
    //If P1 was not a the boundary then zero current is returned so zero current is added to the total current
    return 0.0;
}
//
double Electrical_analysis::Current_of_edges_in_boundary(const int &side, Direct_Electrifying *DEA, const struct Electric_para &electric_param, const vector<Point_3D> &point_list, const vector<double> &radii, const vector<Point_3D> &point_list_gnp, const GCH &hybrid)
{
    //Total curent in boundary edges
    double I_edges = 0.0;
    
    //Scan the triangulation points and look for the point at the right boundary
    for (int i = 0; i < (int)hybrid.triangulation.size(); i++) {
        
        //Get the points of the triangulation edge
        long int P1 = hybrid.triangulation[i][0];
        long int P2 = hybrid.triangulation[i][1];
        
        //Determne flag of triangulation resistor
        //0: CNT-CNT edge
        //1: GNP-GNP edge
        //2: CNT-GNP edge
        int flag = 0;
        
        //Get the nodes
        long int node1, node2;
        //Get the points
        Point_3D point1, point2;
        //Get the radii or thickness
        double rad1, rad2;
        //flags: CNT point (1) or GNP point (0)
        if (hybrid.triangulation_flags[i][0]) {
            //Get the data from CNT
            node1 = DEA->LM_matrix[P1];
            point1 = point_list[P1];
            rad1 = radii[point_list[P1].flag];
        } else {
            //Get the data from GNP
            node1 = DEA->LM_matrix_gnp[P1];
            point1 = point_list_gnp[P1];
            rad1 = hybrid.gnp.hei_z;
        }
        if (hybrid.triangulation_flags[i][1]) {
            //Get the data from CNT
            node2 = DEA->LM_matrix[P2];
            point2 = point_list[P2];
            rad2 = radii[point_list[P2].flag];
        } else {
            //Get the data from GNP
            node2 = DEA->LM_matrix_gnp[P2];
            point2 = point_list_gnp[P2];
            rad2 = hybrid.gnp.hei_z;
        }
        //hout << "flag1="<<hybrid.triangulation_flags[i][0]<<" P1="<<P1<<" node1="<<node1<<" flag2="<<hybrid.triangulation_flags[i][1]<<" P2="<<P2<<" node2="<<node2<<endl;
        //hout << " ("<<point1.x<<","<<point1.y<<","<<point1.z<<") ("<<point2.x<<","<<point2.y<<","<<point2.z<<")"<<endl;
        
        //Voltage
        double V = 0;
        
        //Check where is P1
        //If P1 is on the right side and is a GNP point, then the current is calculated as V(node2) - V(node1)
        if (node1 == side && !hybrid.triangulation_flags[i][0]) {
            
            hout << "flag1="<<hybrid.triangulation_flags[i][0]<<" P1="<<P1<<" node1="<<node1<<" flag2="<<hybrid.triangulation_flags[i][1]<<" P2="<<P2<<" node2="<<node2<<endl;
            hout << " ("<<point1.x<<","<<point1.y<<","<<point1.z<<") ("<<point2.x<<","<<point2.y<<","<<point2.z<<")"<<endl;
            
            //Calculate voltage difference on the triangulation edge, node1 is at the boundary with voltage 0
            //so the voltage drop is from node2 to node1
            V = DEA->voltages[node2] - DEA->voltages[node1];
            hout << "V1="<<DEA->voltages[node1]<<" V2="<<DEA->voltages[node2]<<" DV=" <<V<<endl;
            
            //Calculate the flag for the type of triangulation resistor
            if (!hybrid.triangulation_flags[i][1])
                flag = 1; //GNP-GNP edge
            else
                flag = 2; //CNT-GNP edge
            
            //Calculate resistance of the element
            //The CNT is particle 2, so it has to go first in the arguments
            double Re = DEA->Calculate_resistance_gnp(flag, point2, point1, rad2, rad1, hybrid,  electric_param);
            //Calculate current and add it to the total current
            hout << " Re=" << Re << " I=" << V/Re << endl;
            
            //Add to current
            I_edges = I_edges + V/Re;
            
        }
        //If P2 is on the right side, then the current is calculated as V(node1) - V(node2)
        else if (node2 == side && !hybrid.triangulation_flags[i][1]) {
            
            hout << "flag1="<<hybrid.triangulation_flags[i][0]<<" P1="<<P1<<" node1="<<node1<<" flag2="<<hybrid.triangulation_flags[i][1]<<" P2="<<P2<<" node2="<<node2<<endl;
            hout << " ("<<point1.x<<","<<point1.y<<","<<point1.z<<") ("<<point2.x<<","<<point2.y<<","<<point2.z<<")"<<endl;
            
            //Calculate voltage difference on the triangulation edge, node2 is at the boundary with voltage 0
            //so the voltage drop is from node1 to node2
            V = DEA->voltages[node1] - DEA->voltages[node2];
            hout << "V1="<<DEA->voltages[node2]<<" V2="<<DEA->voltages[node1]<<" DV="<<V<<endl;
            
            //Calculate the flag for the type of triangulation resistor
            if (!hybrid.triangulation_flags[i][0])
                flag = 1; //GNP-GNP edge
            else
                flag = 2; //CNT-GNP edge
            
            //Calculate resistance of the element
            double Re = DEA->Calculate_resistance_gnp(flag, point1, point2, rad1, rad2, hybrid,  electric_param);
            //Calculate current and add it to the total current
            hout << " Re=" << Re << " I=" << V/Re << endl;
            
            //Add to current
            I_edges = I_edges + V/Re;
            
        }
        
    }

    return I_edges;
}
//This function calculates the matrix resistance, depending on the direction of the applied voltage
int Electrical_analysis::Calculate_matrix_resistances(const double &matrix_resistivity, const struct Geom_RVE &geom_rve, vector<double> &matrix_resistances)
{
    //The resistance is calculated as rho*l/A
    double length, A;
    
    //Determine the values of A and length according to the direction in which the voltage is applied
    
    //------------------ Resistance along the x-direction
    length = geom_rve.len_x;
    A = geom_rve.hei_z*geom_rve.wid_y;
    matrix_resistances.push_back(matrix_resistivity*length/A);
    
    
    //------------------ Resistance along the y-direction
    length = geom_rve.wid_y;
    A = geom_rve.len_x*geom_rve.hei_z;
    matrix_resistances.push_back(matrix_resistivity*length/A);
    
    //------------------ Resistance along the z-direction
    length = geom_rve.hei_z;
    A = geom_rve.len_x*geom_rve.wid_y;
    matrix_resistances.push_back(matrix_resistivity*length/A);
    
    return 1;
}
//This function calculates the resistance on each direction from the vector of parallel resistors
int Electrical_analysis::Calculate_resistances_and_resistivities(const struct Geom_RVE &geom_rve, const vector<double> &matrix_resistances, const vector<vector<double> > &paralel_resistors, vector<double> &resistors, vector<double> &resistivities)
{
    //Scan each direction
    //Note that matrix_resistances and paralel_resistors have the same size
    for (int i = 0; i < (int)paralel_resistors.size(); i++) {
        if (!paralel_resistors[i].size()) {
            //If there are no resistors in direction i, then use the conductivity of the matrix
            resistors.push_back(matrix_resistances[i]);
        } else if (paralel_resistors[i].size() == 1) {
            //If this direction has only one resistor, then that is the resistance in that direction
            resistors.push_back(paralel_resistors[i].front());
        } else {
            //If there is more than one resistor, then use the equation to calculate the resistance of parallel resistors
            double R = 0;
            for (int j = 0; j < (int)paralel_resistors[i].size(); j++) {
                R = R + 1/paralel_resistors[i][j];
            }
            resistors.push_back(1/R);
        }
    }
    
    //Calculate resistivities along each direction
    double rho;
    //-----x-direction
    rho = resistors[0]*geom_rve.hei_z*geom_rve.wid_y/geom_rve.len_x;
    resistivities.push_back(rho);
    //-----y-direction
    rho = resistors[1]*geom_rve.hei_z*geom_rve.len_x/geom_rve.wid_y;
    resistivities.push_back(rho);
    //-----z-direction
    rho = resistors[2]*geom_rve.len_x*geom_rve.wid_y/geom_rve.hei_z;
    resistivities.push_back(rho);
    
    return 1;
}