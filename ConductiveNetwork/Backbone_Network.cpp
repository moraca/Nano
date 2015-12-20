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
int Backbone_Network::Determine_backbone_network(const int &family, const vector<int> &cluster, const vector<double> &voltages, const vector<int> &LM_matrix, const vector<vector<long int> > &elements, const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, vector<double> &families_lengths, vector<double> &branches_lengths, vector<vector<long int> > &all_dead_indices, vector<vector<long int> > &all_indices)
{
    //Find the dead branches of each CNT in the cluster. Save the information on the vectors dead_indices and percolated_indices
    if (!Find_dead_branches(voltages, cluster, LM_matrix, elements, structure)) {
        hout << "Error in Determine_backbone_network" << endl;
        return 0;
    }
    
    //Use the dead_indices and percolated_indices to calculate the fraction of CNTs that belong to each family
    if (!Calculate_lengths(family, points_in, families_lengths, branches_lengths)) {
        hout << "Error in Determine_backbone_network" << endl;
        return 0;
    }
    
    //Finally, add the indices to the global vectors so that they can be exported as tecplot files
    Add_indices_to_global_vectors(family, all_dead_indices, all_indices);
    
	return 1;
}


int Backbone_Network::Find_dead_branches(const vector<double> &voltages, const vector<int> &cluster, const vector<int> &LM_matrix, const vector<vector<long int> > &elements, const vector<vector<long int> > &structure)
{
    //Define the 'Zero-cutoff' of the system
    double zero_cutoff = Zero_voltage(voltages, cluster, LM_matrix, elements);
    
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
    long int P1, node1, P2, node2;
    //Variable for the current
    double I;
    
    //Scan all the CNTs that belong to the current cluster
    for (long int i = 0; i < (long int)cluster.size(); i++) {
        CNT = cluster[i];
        //The dead branches are in th extremes of a CNT, they cannot be in the middle.
        //Hence I just need to check the current of the first and last elements of the CNT
        
        //Check all elements from the begining until one that carries current is found (if any)
        for (long int j = 0; j < (long int)elements[CNT].size()-1; j++) {
            P1 = elements[CNT][j];
            node1 = LM_matrix[P1];
            P2 = elements[CNT][j+1];
            node2 = LM_matrix[P2];
            I = abs(voltages[node2] - voltages[node1]);
            //hout << "CNT="<<CNT<<" P1="<<P1<<" node1="<<node1<<" P2="<<P2<<" node2="<<node2<<" I="<<I<<endl;
            if (I < zero_cutoff){
                if ( !dead_indices[i].size() ){
                    //If the vector dead_indices[i] is empty, then push P1 as the first index
                    dead_indices[i].push_back(P1);
                    if (percolated_indices[i].size() == 1) {
                        //if the vector percolated_indices[i] has one element
                        //then the conducting segment ends and the non-conducting segment starts
                        //Hence, the conductive segment has last index P1
                        percolated_indices[i].push_back(P1);
                        //Also, the rest of the CNT is not-conducting so break the loop
                        break;
                    }
                } else if (dead_indices[i].size() == 2) {
                    //If the dead_indices[i] vector has size two, it means that there is a middle conducting segment and 2
                    //non-conducting segments on both ends of the CNT. The only way for dead_indices[i] to have size two is
                    //that it had size one and a conducting segment was found
                    //P1 will be the first index of the second segment
                    dead_indices[i].push_back(P1);
                    //P1 will be the second index of the conducting segement
                    percolated_indices[i].push_back(P1);
                    //If it reaches this part, then there is no point in checking the rest of the CNT, as for sure
                    //the rest of the CNT is not conducting. Hence, break the loop
                    break;
                }
            } else{
                if (dead_indices[i].size() == 1){
                    //If the current is greater than the cutoff and the dead_indices[i] has only one element,
                    //the non-conducting segment ends and the conducting segment starts
                    //Hence, the non-conductive segment at the front has as last index P1
                    dead_indices[i].push_back(P1);
                }
                //If the percolated_indices[i] is empty, then P1 is the first index
                if (!percolated_indices[i].size()) {
                    percolated_indices[i].push_back(P1);
                }
            }
        }
        
        //Take care of the last point of the CNT
        P2 = elements[CNT].back();
        //The last point of the CNT is the second index of the non-conducting segment only if I < zero_cutoff
        //otherwise it is the second index of the conducting segment
        if (I < zero_cutoff) {
            dead_indices[i].push_back(P2);
        } else {
            percolated_indices[i].push_back(P2);
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
        //Check that the vectors have the corrent size. dead_indices can only have size 0, 2 or 4
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
    
    return 1;
}

//This function determines the cutoff for "zero current"
//This step is necessary due to floating point errors
double Backbone_Network::Zero_voltage(const vector<double> &voltages, const vector<int> &cluster, const vector<int> &LM_matrix, const vector<vector<long int> > &elements)
{
    //Variable for the current
    double I;
    //Variable to store the cutoff for zero voltage
    double zero_cutoff;
    //Variables
    int CNT;
    long int P1, node1, P2, node2;
    //Vector to store the currents
    vector<double> currents;
    
    //First calculate all currents
    for (long int i = 0; i < (long int)cluster.size(); i++) {
        CNT = cluster[i];
        for (long int j = 0; j < (long int)elements[CNT].size()-1; j++) {
            P1 = elements[CNT][j];
            node1 = LM_matrix[P1];
            P2 = elements[CNT][j+1];
            node2 = LM_matrix[P2];
            I = abs(voltages[node2] - voltages[node1]);
            currents.push_back(I);
        }
    }
    
    //Sort currents
    sort(currents.begin(),currents.end());
    //Print1DVec(currents, "currents.txt");
    //vector<double> cutoffs;
    //Find the cutoff. 
    for(int i = (int) currents.size()-1; i >= 0 ; i--){
        //Temporarily store the drop in current in the variable zero_cutoff
        zero_cutoff = currents[i-1]/currents[i];
        //zero_cutoff = currents[i]/currents.back();
        //cutoffs.push_back(zero_cutoff);
        //hout << currents[i-1] << ' ' << currents[i] << ' ' << zero_cutoff << endl;
        //When there is a jump in the order of magnitude then that is the cutoff for zero current
        if (zero_cutoff < 1E-4){
            if (zero_cutoff < 1E-15){
                //In case the gap is from 1e-X to 0, where is 15 or more, then probably the zero current
                //should be a drop of ten orders of magnitude as it was before, otherwise, use the
                //gap to determine the zero cutoff
                zero_cutoff = currents.front()*1e-10;
                hout << "Gap in current = ("<<currents[i]<< ", "<<currents[i-1] << "), zero_cutoff = " << zero_cutoff << endl;
                break;
            } else {
                //A large enough drop in the order of magnitude of the current is taken to be the cutoff for "zero current"
                zero_cutoff = currents[i-1]*10;
                hout << "Gap in current = ("<<currents[i]<< ", "<<currents[i-1] << "), zero_cutoff = " << zero_cutoff << endl;
                break;
            }
        }
    }
    Printer *P = new Printer;
    P->Print_1d_vec(currents, "currents.txt");
    //P->Print_1d_vec(cutoffs, "cutoffs.txt");
    delete P;
    return zero_cutoff;
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