//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Percolation.cpp
//OBJECTIVE:	Determine which clusters percolate and in which directions
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Percolation.h"

int Percolation::Determine_percolated_clusters(const int &window, const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const struct GNP_Geo &gnps, const vector<vector<int> > &boundary_cnt, const vector<int> &labels_cnt, const vector<vector<int> > &boundary_gnp, const vector<int> &labels_gnp, vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated_cnt, vector<vector<int> > &clusters_gnp, vector<vector<int> > &isolated_gnp)
{
    //flag to ignore variables for second particle
    int ignore_particle = 1;
    //Determine percoaltion in clusters made of only CNTs
    if (sample.particle_type == "CNT_wires") {
        //Check if any CNT cluster percolates
        //CNT vectors go first so GNP vectors are ignored
        if (!Cluster_percolation(ignore_particle, boundary_cnt, labels_cnt, clusters_cnt, isolated_cnt, boundary_gnp, labels_gnp, clusters_gnp, isolated_gnp)) {
            hout << "Error in Determine_percolating_clusters when calling Cluster_percolation for CNTs" << endl;
            return 0;
        }
        //Check if single CNTs percolate
        if (!Single_particle_percolation(window, sample, cnts, gnps, boundary_cnt, labels_cnt, clusters_cnt, isolated_cnt)) {
            hout << "Error in Determine_percolating_clusters when calling Single_particle_percolation for CNTs" << endl;
            return 0;
        }
    }
    //Determine percoaltion in clusters made of only GNPs
    else if (sample.particle_type == "GNP_cuboids") {
        //Check if any GNP cluster percolates
        //GNP vectors go first so CNT vectors are ignored
        if (!Cluster_percolation(ignore_particle, boundary_gnp, labels_gnp, clusters_gnp, isolated_gnp, boundary_cnt, labels_cnt, clusters_cnt, isolated_cnt)) {
            hout << "Error in Determine_percolating_clusters when calling Cluster_percolation for GNPs" << endl;
            return 0;
        }
    }
    //Determine percoaltion in clusters made of both particles
    else {
        //Do not ignore the second particle
        ignore_particle = 0;
        //Check if any cluster percolates
        if (!Cluster_percolation(ignore_particle, boundary_cnt, labels_cnt, clusters_cnt, isolated_cnt, boundary_gnp, labels_gnp, clusters_gnp, isolated_gnp)) {
            hout << "Error in Determine_percolating_clusters when calling Cluster_percolation for mixed fillers" << endl;
            return 0;
        }
    }
    
    return 1;
}
//This function determines if the clusters in the clusters_cnt variable actually percolate
int Percolation::Cluster_percolation(const int &ignore, const vector<vector<int> > &boundary_1, const vector<int> &labels_1, vector<vector<int> > &clusters_1, vector<vector<int> > &isolated_1, const vector<vector<int> > &boundary_2, const vector<int> &labels_2, vector<vector<int> > &clusters_2, vector<vector<int> > &isolated_2)
{
    //Check if there is any cluster at all. If the size of clusters_cnt is non-zero, then there are clusters
    //If the two particles are present, both clusters should have the same size
    if (clusters_1.size()) {
        //Initialize the flag vectors
        vector<short int> zeros(6,0);
        vector<vector<short int> > perc_flag(clusters_1.size(), zeros);
        //hout << "fill flags " << endl;
        //Scan all the boundary vectors and fill the flag vectors
        Fill_percolation_flags_all_directions(boundary_1, perc_flag, labels_1);
        //Check the vectors for the second particle when ignore flag is zero
        if (!ignore) {
            //hout << "fill flags 2" << endl;
            Fill_percolation_flags_all_directions(boundary_2, perc_flag, labels_2);
        }
        //hout << "all clusters " << endl;
        
        //Move all non-percolating clusters to the corresponding vector<vector>
        if (!Check_percolation_all_clusters(ignore, perc_flag, clusters_1, isolated_1, clusters_2, isolated_2, family)) {
            hout << "Error in Determine_percolating_clusters >> Check_percolation_all_clusters" << endl;
            return 0;
        }
        //hout << "clusters_1 end" << endl;
    }
    return 1;
}
//This function calls the Fill_percolation_flags_single_direction and usese three flags: px, py and pz.
//The flags are used in case one direction is not necessary to check for percolation. This is specially useful for single CNT percolation
void Percolation::Fill_percolation_flags_all_directions(const vector<vector<int> > &boundary, vector<vector<short int> > &perc_flag, const vector<int> &labels)
{
    //Scan all the boundary vectors and fill the flag vectors
    //X-direction
    //hout << "X-direction 0, boundary.size()="<<boundary.size()<<endl;
    Fill_percolation_flags_single_direction(boundary[0], 0, perc_flag, labels);
    //hout << "X-direction 1"<<endl;
    Fill_percolation_flags_single_direction(boundary[1], 1, perc_flag, labels);
    //Y-direction
    //hout << "Y-direction 0"<<endl;
    Fill_percolation_flags_single_direction(boundary[2], 2, perc_flag, labels);
    //hout << "Y-direction 1"<<endl;
    Fill_percolation_flags_single_direction(boundary[3], 3, perc_flag, labels);
    //Z-direction
    //hout << "Z-direction 0"<<endl;
    Fill_percolation_flags_single_direction(boundary[4], 4, perc_flag, labels);
    //hout << "Z-direction 1"<<endl;
    Fill_percolation_flags_single_direction(boundary[5], 5, perc_flag, labels);
    //hout << "all-direction"<<endl;
}

//This function fills the percolation_flags vector, which is usde to determine percolation
void Percolation::Fill_percolation_flags_single_direction(const vector<int> &boundary_vector, int boundary_number, vector<vector<short int> > &perc_flag, const vector<int> &labels)
{
    //Scan all particles at the boundary
    for (int i = 0; i < (int)boundary_vector.size(); i++) {
        //Current particle
        int particle = boundary_vector[i];
        //hout <<"particle="<<particle<<endl;
        //Label of the CNT, which is also the cluster number
        int L = labels[particle];
        //Chek if the boundary CNT has a valid label
        if (L != -1) {
            //Turn on the flag of the cluster L that corresponds to the boundary boundary_number
            perc_flag[L][boundary_number] = 1;
        }
    }
}
//Scan cluster by cluster and check whether it percolates or not
int Percolation::Check_percolation_all_clusters(const int &ignore, vector<vector<short int> > &perc_flag, vector<vector<int> > &clusters_1, vector<vector<int> > &isolated_1, vector<vector<int> > &clusters_2, vector<vector<int> > &isolated_2, vector<int> &family)
{
    //Move all non-percolating clusters to the corresponding vector<vector>
    //I start from the end of the vector to avoid issues with the index when removing an element
    //This variable will store the family number given by Check_percolation_single_cluster
    int fam;
    for (int i = (int)clusters_1.size() - 1; i >= 0 ; i--) {
        //hout <<"Check_clusters_percolation " << i << " size " << clusters.size() << endl;
        //hout << "percolation_flags.size()=" << percolation_flags.size() << endl;
        if (Check_percolation_single_cluster(perc_flag[i], fam)){
            //Add the family number to the vector of families
            family.insert(family.begin(), fam);
        } else {
            //percolation_flags and clusters have the same size
            
            //copy the non-percolated cluster to the vector of isolated clusters
            isolated_1.push_back(clusters_1[i]);
            //remove the non-percolating cluster
            clusters_1.erase(clusters_1.begin()+i);
            
            //If the vectors for the second particle are not ignored, do the same for that particle
            if (!ignore) {
                //copy the non-percolated cluster to the vector of isolated clusters
                isolated_2.push_back(clusters_2[i]);
                //remove the non-percolating cluster
                clusters_2.erase(clusters_2.begin()+i);
            }
        }
    }
    return 1;
}
//This function will check if there is percolation for a single cluster in x, y and/or z directions
int Percolation::Check_percolation_single_cluster(const vector<short int> &cluster_flag, int &family)
{
    //Falgs that will tell me the family the cluster belogns to
    int percolating[] = {0, 0, 0, 0, 0, 0, 0};
    
    //Boolean operations to find the different percolation directions
    percolating[0] = cluster_flag[0] && cluster_flag[1]; //x-x only
    percolating[1] = cluster_flag[2] && cluster_flag[3]; //y-y only
    percolating[2] = cluster_flag[4] && cluster_flag[5]; //z-z only
    percolating[3] = percolating[0] && percolating[1]; //x-x and y-y only
    percolating[4] = percolating[0] && percolating[2]; //x-x and z-z only
    percolating[5] = percolating[1] && percolating[2]; //y-y and z-z only
    percolating[6] = percolating[0] && percolating[1] && percolating[2]; //x-x, y-y and z-z
    
    //Scan the percolating array backwards
    for (int i = 6; i >=0 ; i--)
        if (percolating[i]){
            //If there is a non-zero percolating[i], then this cluster percolates and belongs to family i
            family = i;
            return 1;
        }
    //If all percolating[i] were 0, then there is no percolation in the cluster
    return 0;
}
//It assumed that this function is called when the CNTs have a lenth equal or grater than the dimensions of the sample
//In this function, the vector of isolated CNTs is scanned to look for percolated CNTs.
//A single CNT can only percolate on X, Y or Z.
int Percolation::Single_particle_percolation(const int &window, const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const struct GNP_Geo &gnps, const vector<vector<int> > &boundary, const vector<int> &labels, vector<vector<int> > &clusters, vector<vector<int> > &isolated)
{
    
    //These variables are flags to determine in which directions there might be percolation by a single CNT
    //hout << "0 ";
    int px = 0, py= 0, pz = 0;
    if (!Determine_direction_flags(window, sample, cnts, gnps, px, py, pz)) {
        //If in all three directions the particles are smaller than the dimensions of the observation window, then there is nothing more
        //to do in this function. This happens when px, px and py are zero, so the only way the sum is zero is when all three are zero.
        //Then, since Determine_direction_flags returns the sum of the flags, when the sum is zero, terminate the function
        return 1;
    }
    
    //This is similar to the labels generated by HK76. Since it is only needed when the particles are smaller than a dimension of
    //the observation window, it is created here and not during the HK76
    //The labels vector has the same size as the number of particles, so I use it to initialize labels_iso
    vector<int> labels_iso(labels.size(), -1);
    int single_particles = Determine_isolated_labels(isolated, labels_iso);
    
    //Check if there is any cluster at all. If single_CNTs is non-zero, then there are clusters
    if (single_particles) {
        //Initialize the flag vectors
        vector<short int> zeros(6,0);
        vector<vector<short int> > perc_flag(single_particles, zeros);
        //hout << "10 ";
        //Scan all the boundary vectors and fill the flag vectors
        Fill_single_particle_percolation_flags_all_directions(boundary, perc_flag, labels_iso, px, py, pz);
        //hout << "11 ";
        
        //The percolated CNTs will be stored in the clusters_tmp vector
        if (!Check_particles_for_percolation(perc_flag, clusters, isolated)) {
            hout << "Error in Single_CNT_percolation" << endl;
            return 0;
        }
    }
    return 1;
}
int Percolation::Determine_direction_flags(const int &window, const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const struct GNP_Geo &gnps, int &px, int &py, int &pz)
{
    //These are variables for the geometry of the observation window
    //Dimensions of the current observation window
    double w_x = sample.win_max_x - window*sample.win_delt_x;
    double w_y = sample.win_max_y - window*sample.win_delt_y;
    double w_z = sample.win_max_z - window*sample.win_delt_z;
    
    //Check CNT geometry if the sample has only CNTs
    if (sample.particle_type == "CNT_wires") {
        if (cnts.len_max >= w_x)
            px = 1;
        if (cnts.len_max >= w_y)
            py = 1;
        if (cnts.len_max >= w_z)
            pz = 1;
    }
    
    //Check GNP geometry if the sample has only GNPs
    else if (sample.particle_type == "GNP_cuboids") {
        if (gnps.len_max >= w_x)
            px = 1;
        if (gnps.len_max >= w_y)
            py = 1;
        if (gnps.len_max >= w_z)
            pz = 1;
    }
    
    //Return the sum of the percoaltion diraction flags
    return (px+py+pz);
}
int Percolation::Determine_isolated_labels(const vector<vector<int> > &isolated, vector<int> &labels_iso)
{
    //Number of single particles
    int single_particles = 0;
    //Fill label_map_iso
    for (int i = 0; i < (int)isolated.size(); i++) {
        //This function is called after the clusters are checked for percolation, so in the isolated vector there are single particles and clusters
        //As I am interested in single particles (isolated[i].size() == 1), when I find clusters (isolated[i].size() > 1) I break the loop
        if (isolated[i].size() > 1)
            break;
        //hout << "isolated["<<i<<"].size()="<<isolated[i].size()<<endl;
        int particle = isolated[i].front();
        labels_iso[particle] = i;
        single_particles++;
    }
    
    //Return the number single particle clusters as an indicator to proceed on checking for single particle percolation
    return single_particles;
}
//This function calls the Fill_percolation_flags_single_direction and usese three flags: px, py and pz.
//The flags are used in case one direction is not necessary to check for percolation. This is specially useful for single particle percolation
void Percolation::Fill_single_particle_percolation_flags_all_directions(const vector<vector<int> > &boundary, vector<vector<short int> > &perc_flag, const vector<int> &labels_iso, const int &px, const int &py, const int &pz)
{
    //Scan all the boundary vectors and fill the flag vectors
    //X-direction
    if (px) {
        Fill_single_particle_percolation_flags_single_direction(boundary[0], 0, perc_flag, labels_iso);
        Fill_single_particle_percolation_flags_single_direction(boundary[1], 1, perc_flag, labels_iso);
    }
    //Y-direction
    if (py) {
        Fill_single_particle_percolation_flags_single_direction(boundary[2], 2, perc_flag, labels_iso);
        Fill_single_particle_percolation_flags_single_direction(boundary[3], 3, perc_flag, labels_iso);
    }
    //Z-direction
    if (pz) {
        Fill_single_particle_percolation_flags_single_direction(boundary[4], 4, perc_flag, labels_iso);
        Fill_single_particle_percolation_flags_single_direction(boundary[5], 5, perc_flag, labels_iso);
    }
}
//This function fills the percolation_flags vector, which is usde to determine percolation
void Percolation::Fill_single_particle_percolation_flags_single_direction(const vector<int> &boundary_vector, const int &boundary_number, vector<vector<short int> > &perc_flag, const vector<int> &labels_iso)
{
    //Variables
    int L, particle;
    
    for (int i = 0; i < (int)boundary_vector.size(); i++) {
        //Current CNT
        particle = boundary_vector[i];
        //hout << "particle=" << particle << ' ';
        //Label of the particle, which also corresponds to its perc_flag number
        L = labels_iso[particle];
        //If the label is -1, then the particle is not an isolated cluster or it belongs to a cluster with more than one particle
        //Hence, if the label is different from -1 turn on the flag
        if (L != -1) {
            //Activate the flag corresponding to the boundary
            perc_flag[L][boundary_number] = 1;
        }
    }
}
int Percolation::Check_particles_for_percolation(vector<vector<short int> > &perc_flag, vector<vector<int> > &clusters, vector<vector<int> > &isolated)
{
    //Check if every cluster made of a single cluster percolates
    //I start from the end of the vector of percolated flags to avoid issues with the index when removing an element
    //This variable (fam) will store the family number given by Check_percolation_single_direction
    int fam, counter = 0;
    for (int i = (int)perc_flag.size() - 1; i >= 0 ; i--) {
        //hout <<"Check_clusters_percolation " << i << " isolated.size " << isolated.size() << endl;
        //hout << "perc_flag.size()=" << perc_flag.size() << endl;
        //If there is percolation, then add the isolated CNT to the clusters_cnt vector
        if (Check_percolation_single_particle(perc_flag[i], fam)){
            clusters.push_back(isolated[i]);
            //remove the non-percolating cluster
            isolated.erase(isolated.begin()+i);
            //Add the family number to the vector of families
            family.push_back(fam);
            //Increase the counter of percolated clusters
            counter++;
        }
    }
    hout << "Single particle percolated clusters: " << counter <<endl;
    return 1;
}
//This function will check if there is percolation for a single cluster consisting of a single particle in x, y or z directions
int Percolation::Check_percolation_single_particle(const vector<short int> &cluster_flag, int &family)
{
    //Check the three possible cases of single particle percolation
    //Falgs that will tell me the family the cluster belogns to
    int percolating[] = {0, 0, 0};
    
    //Boolean operations to find the different percolation directions
    percolating[0] = cluster_flag[0] && cluster_flag[1]; //x-x only
    percolating[1] = cluster_flag[2] && cluster_flag[3]; //y-y only
    percolating[2] = cluster_flag[4] && cluster_flag[5]; //z-z only
    
    //Scan the percolating array backwards
    for (int i = 0; i < 3 ; i++)
        if (percolating[i]){
            //If there is a non-zero percolating[i], then this cluster percolates and belongs to family i
            family = i;
            return 1;
        }
    //If all percolating[i] were 0, then there is no percolation in the cluster
    return 0;
}