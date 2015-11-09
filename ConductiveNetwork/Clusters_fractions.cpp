//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Clusters_fractions.h
//OBJECTIVE:
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Clusters_fractions.h"
#include <iostream>

//This function calcultes the fractions of the percolated families. It also calculates the lengths of the non-percolated CNTs
//There is a check of the total length
int Clusters_fractions::Calculate_fractions(const vector<vector<long int> > &structure, const vector<int> &cnts_inside, const vector<vector<int> > &isolated, const vector<Point_3D> points_in, vector<double> &families_lengths, vector<double> &branches_lengths, vector<double> fractions)
{
    //First calculate the length of all CNTs
    int CNT;
    double total_length = 0;
    for (int i = 0; i < (int)cnts_inside.size(); i++) {
        CNT = cnts_inside[i];
        total_length = total_length + CNT_length(structure, points_in, CNT);
    }
    
    //Calculate the fractions of the percoalted families, i.e. from i=0 to i=6
    for (int i = 0; i <= 6; i++) {
        fractions[i] = families_lengths[i]/total_length;
    }
    
    //Now calculate the length of the isolated CNTs
    double isolated_length = 0;
    for (int i = 0; i < (int)isolated.size(); i++) {
        for (int j = 0; j < (int)isolated[i].size(); j++) {
            CNT = isolated[i][j];
            isolated_length = isolated_length + CNT_length(structure, points_in, CNT);
        }
    }
    
    //Add the length of isolated CNTs into the families_lengths vector
    double dead_branches = 0;
    for (int i = 0; i < (int)branches_lengths.size(); i++) {
        dead_branches = dead_branches + branches_lengths[i];
    }
    families_lengths[7] = isolated_length +  dead_branches;
    //Calculate the fraction of the NP family
    fractions[7] = families_lengths[7]/total_length;
    
    //Check that the families lengths and the total length are the same (or nearly the same)
    double check_length = 0;
    for (int i = 0; i < (int)families_lengths.size(); i++) {
        check_length = check_length + families_lengths[i];
    }
    families_lengths.push_back(check_length);
    families_lengths.push_back(total_length);
    
    //Append vectors to files
    Append_1d_vector_to_file(families_lengths, "clusters_lengths.txt");
    Append_1d_vector_to_file(branches_lengths, "dead_branches_lengths.txt");
    Append_1d_vector_to_file(fractions, "clusters_fractions.txt");
    
    //Calculate the fraction of the geometrically percolated clusters
    vector<double> geom_fractions(8,0);
    for (int i = 0; i < 7; i++) {
        geom_fractions[i] = (families_lengths[i] + branches_lengths[i])/total_length;
        //The total fractions of geometrically percolated clusters
        geom_fractions[7] = geom_fractions[7] + geom_fractions[i];
    }
    
    Append_1d_vector_to_file(geom_fractions, "geom_fractions.txt");
    
    if (abs((check_length - total_length)/total_length) >  Zero){
        hout << "Error in Calculate_fractions. The total length of the CNTs in the observation window does not match with "<<endl;
        hout << "the length of the CNTs in the percolated and non-percolated clusters. " << endl;
        hout << setwp(1,20) << "Length in clusters = " << check_length << endl;
        hout << "Length of CNTs inside the observation window = " << total_length << endl;
        return 0;
    }
    
    
    return 1;
}

//This function calculates the length of a segement of a CNT given by a sequence of points from (global number) index1 to index2
double Clusters_fractions::CNT_length(const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, int CNT)
{
    //Variable to store the length
    double length = 0;
    //Variables to store the point numbers
    long int P1, P2;
   for (int k = 0; k < (int)structure[CNT].size()-1; k++) {
        P1 = structure[CNT][k];
        P2 = structure[CNT][k+1];
        length = length + points_in[P1].distance_to(points_in[P2]);
    }
    return length;
}

//Print a vector of doubles with the specified filename
void Clusters_fractions::Append_1d_vector_to_file(const vector<double> &list, const string &filename)
{
    ofstream otec(filename.c_str(), std::ios_base::app);
    hout << "Appending to file: " << filename << endl;
    for (long int i=0; i < list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\t";
    }
    otec << "\n";
    otec.close();
}
