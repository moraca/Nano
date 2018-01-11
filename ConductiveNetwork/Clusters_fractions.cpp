//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Clusters_fractions.cpp
//OBJECTIVE:
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Clusters_fractions.h"
#include <iostream>

//This function calcultes the fractions of the percolated families. It also calculates the lengths of the non-percolated CNTs
//There is a check of the total length
int Clusters_fractions::Calculate_fractions(const struct Geom_RVE &sample, const vector<int> &cnts_inside, const vector<int> &gnps_inside, const vector<vector<long int> > &structure, const vector<Point_3D> points_in, const vector<double> &radii, const vector<vector<int> > &isolated, const vector<GCH> &hybrid_particles, const vector<vector<int> > &isolated_gnp, const vector<vector<long int> > &all_dead_indices, const vector<vector<long int> > &all_percolated_indices, const vector<vector<int> > &all_dead_gnps, const vector<vector<int> > &all_percolated_gnps)
{
    
    //Variables to store the volumes of particles
    vector<double> volumes_cnt(8,0);
    vector<double> volumes_gnp(8,0);
    
    //Variable to store the volumes of dead particles
    vector<double> dead_branches(7,0);
    vector<double> dead_gnps(7,0);
    
    //Variables to store total volumes of each particle
    double cnt_volume = 0.0;
    double gnp_volume = 0.0;
    
    //Variables to store total volumes
    double total_volume = 0.0;
    
    //Variables to store total volumes per family
    vector<double> total_volumes(8,0);
    vector<double> dead_particles(7,0);
    
    //=================================================================================
    //CNT volumes
    if (cnts_inside.size()) {
        
        //Calculate the volumes
        if (!Calculate_volumes_cnt(cnts_inside, structure, points_in, radii, isolated, all_dead_indices, all_percolated_indices, volumes_cnt, dead_branches, cnt_volume)) {
            hout << "Error in Calculate_fractions when calling Calculate_volumes_cnt" << endl;
            return 0;
        }
        
        //Append to files for CNT volumes and fractions
        if (!Append_to_volume_and_fraction_files(volumes_cnt, cnt_volume, "cnts_wrt_cnts")) {
            hout << "Error in Calculate_fractions when calling Append_to_volume_and_fraction_files" << endl;
            return 0;
        }
    }
    
    //=================================================================================
    //GNP volumes
    if (gnps_inside.size()) {
        
        //Calculate the volumes
        if (!Calculate_volumes_gnp(sample, gnps_inside, hybrid_particles, isolated_gnp, all_dead_gnps, all_percolated_gnps, volumes_gnp, dead_gnps, gnp_volume)) {
            hout << "Error in Calculate_fractions when calling Calculate_volumes_gnp" << endl;
            return 0;
        }
        
        //Append to files for GNP volumes and fractions
        if (!Append_to_volume_and_fraction_files(volumes_gnp, gnp_volume, "gnps_wrt_gnps")) {
            hout << "Error in Calculate_fractions when calling Append_to_volume_and_fraction_files" << endl;
            return 0;
        }
    }
    
    //=================================================================================
    //Total fractions
    
    //Total volume from particles volume
    total_volume = cnt_volume + gnp_volume;
    
    //Total volume from volumes vector
    double volume_check = 0.0;
    
    //Calculate volumes per family
    for (int i = 0; i < (int)total_volumes.size(); i++) {
        total_volumes[i] = volumes_cnt[i] + volumes_gnp[i];
        volume_check = volume_check + total_volumes[i];
    }
    for (int i = 0; i < (int)dead_particles.size(); i++) {
        dead_particles[i] = dead_branches[i] + dead_gnps[i];
    }
    
    
    //Append to files for total volumes and fractions
    if (!Append_to_volume_and_fraction_files(total_volumes, total_volume, "clusters")) {
        hout << "Error in Calculate_fractions when calling Append_to_volume_and_fraction_files" << endl;
        return 0;
    }
    if (!Append_to_volume_and_fraction_files(dead_branches, total_volume, "dead_particles")) {
        hout << "Error in Calculate_fractions when calling Append_to_volume_and_fraction_files" << endl;
        return 0;
    }
    
    //=================================================================================
    //Check that the two ways of calculating the volume give the same result
    if (abs((volume_check - total_volume)/total_volume) >  Zero){
        hout << "Error in Calculate_fractions. The total volume of the particles in the observation window does not match with "<<endl;
        hout << "the volume of the particles in the percolated and non-percolated clusters. " << endl;
        hout << setwp(1,20) << "Volume in clusters = " << volume_check << endl;
        hout << "Volume of particles inside the observation window = " << total_volume << endl;
        return 0;
    }
    
    return 1;
}
//
int Clusters_fractions::Calculate_volumes_cnt(const vector<int> &cnts_inside, const vector<vector<long int> > &structure, const vector<Point_3D> points_in, const vector<double> &radii, const vector<vector<int> > &isolated, const vector<vector<long int> > &all_dead_indices, const vector<vector<long int> > &all_percolated_indices, vector<double> &volumes_cnt, vector<double> &dead_branches, double &cnt_volume)
{
    
    //Calculate the volume of all CNTs inside the sample
    cnt_volume = CNT_cluster_volume(cnts_inside, structure, points_in, radii);
    
    //Variable to store the volume of dead branches
    double dead_volume = 0.0;
    
    //Calculate the volumes of percolated CNTs and dead branches, i.e. from i=0 to i=6
    for (int i = 0; i <= 6; i++) {
        
        //Volume of percolated CNTs in family i
        volumes_cnt[i] = CNT_indices_volume(all_percolated_indices[i], points_in, radii);
        
        //Volume of dead branches in family i
        dead_branches[i] = CNT_indices_volume(all_dead_indices[i], points_in, radii);
        
        //Add the total volume of dead branches
        dead_volume = dead_volume + dead_branches[i];
    }
    
    //Calculate the volume of isolated CNTs
    for (int i = 0; i < (int)isolated.size(); i++) {
        
        //Add the volume for each isolated cluster
        volumes_cnt[7] = volumes_cnt[7] + CNT_cluster_volume(isolated[i], structure, points_in, radii);
    }
    
    //Add the volume of isolated CNTs and dead branches
    volumes_cnt[7] = volumes_cnt[7] + dead_volume;
    
    return 1;
}
//Function that calculates the volume of a cluster of CNTs
double Clusters_fractions::CNT_cluster_volume(const vector<int> &cluster, const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const vector<double> &radii)
{
    //Calculate the volume of all CNTs
    double volume = 0;
    for (int i = 0; i < (int)cluster.size(); i++) {
        
        //Current CNT
        int CNT = cluster[i];
        
        //Add the volume of current CNT
        volume = volume + CNT_length(structure, points_in, CNT)*PI*radii[CNT]*radii[CNT];
    }
    
    return volume;
}
//This function calculates the length of a segement of a CNT given by a sequence of points from (global number) index1 to index2
double Clusters_fractions::CNT_length(const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const int &CNT)
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
//Function that calculates the volume of a cluster of CNTs given by indices
double Clusters_fractions::CNT_indices_volume(const vector<long int> &indices, const vector<Point_3D> &points_in, const vector<double> &radii)
{
    //Calculate the volume of all CNTs
    double volume = 0;
    for (int i = 0; i < (int)indices.size(); i=i+2) {
        
        //First point
        long int P1 = indices[i];
        
        //Second point
        long int P2 = indices[i+1];
        
        //Get CNT number
        int CNT = points_in[P1].flag;
        
        //Add the volume of CNT segment
        for (long int j = P1; j < P2; j++) {
            volume = volume + points_in[j].distance_to(points_in[j+1])*PI*radii[CNT]*radii[CNT];
        }
    }
    
    return volume;
}
//
int Clusters_fractions::Calculate_volumes_gnp(const struct Geom_RVE &sample, const vector<int> &gnps_inside, const vector<GCH> &hybrid_particles, const vector<vector<int> > &isolated_gnps, const vector<vector<int> > &all_dead_gnps, const vector<vector<int> > &all_percolated_gnps, vector<double> &volumes_gnp, vector<double> &dead_gnps, double &gnp_volume)
{
    //Generate cuboid to use function that approximates GNP volume
    struct cuboid gvcub;
    gvcub.poi_min = sample.origin;
    gvcub.len_x = sample.len_x;
    gvcub.wid_y = sample.wid_y;
    gvcub.hei_z = sample.hei_z;
    
    //Calculate the volume of all GNPs inside the sample
    gnp_volume = GNP_cluster_volume(gvcub, gnps_inside, hybrid_particles);
    
    //Variable to store the volume of dead GNPs
    double dead_volume = 0.0;
    
    //Calculate the volumes of the percolated and dead GNPs, i.e. from i=0 to i=6
    for (int i = 0; i <= 6; i++) {
        
        //Volume of percolated GNPs in family i
        volumes_gnp[i] = GNP_cluster_volume(gvcub, all_percolated_gnps[i], hybrid_particles);
        
        //Volume of dead GNPs in family i
        dead_gnps[i] = GNP_cluster_volume(gvcub, all_dead_gnps[i], hybrid_particles);
        
        //Add the total volume of dead GNPs
        dead_volume = dead_volume + dead_gnps[i];
    }
    
    //Calculate the volume of isolated GNPs
    for (int i = 0; i < (int)isolated_gnps.size(); i++) {
        
        //Add the volume for each isolated cluster
        volumes_gnp[7] = volumes_gnp[7] + GNP_cluster_volume(gvcub, isolated_gnps[i], hybrid_particles);
    }
    
    //Add the volume of isolated and dead GNPs
    volumes_gnp[7] = volumes_gnp[7] + dead_volume;
    
    return 1;
}
//
double Clusters_fractions::GNP_cluster_volume(const struct cuboid &gvcub, const vector<int> &cluster, const vector<GCH> &hybrid_particles)
{
    //Generate a temporary vector to take advantage of the function that approximates the volume
    //of a GNP in the GenNetwork class
    GenNetwork *GN_tmp = new GenNetwork;
    
    //Variable to store the total volume of GNPs
    double volume = 0.0;
    
    //Variable to store the volume of a single GNP
    double single_gnp = 0.0;
    
    //Calculate the volume of all GNPs inside the sample
    for (int i = 0; i < (int)cluster.size(); i++) {
        
        //Current GNP
        int GNP = cluster[i];
        
        //Approximate volume of GNP
        if (!GN_tmp->Approximate_gnp_volume(gvcub, hybrid_particles[GNP], single_gnp)) {
            return 0;
        }
        
        //Add volume of GNP to total
        volume = volume + single_gnp;
    }
    
    //Delete temporary object
    delete GN_tmp;
    
    return volume;
}
//
int Clusters_fractions::Append_to_volume_and_fraction_files(const vector<double> &volumes_vec, const double &volume, const string &filename)
{
    //Create printer object
    Printer *P = new Printer;
    
    //Create filename for volumes
    string filename_vol = filename;
    filename_vol = filename_vol.append("_volumes.txt");
    
    //Append volumes
    P->Append_1d_vec(volumes_vec, filename_vol);
    
    //Create filename for fractions
    string filename_frac = filename;
    filename_frac = filename_frac.append("_fractions.txt");
    
    //Calculate fractions
    vector<double> fractions;
    for (int i = 0; i < (int)volumes_vec.size(); i++) {
        fractions.push_back(volumes_vec[i]/volume);
    }
    
    //Append fractions
    P->Append_1d_vec(fractions, filename_frac);
    
    //Delete printer object
    delete P;
    
    return 1;
}
