//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Percolation.cpp
//OBJECTIVE:	Determine which clusters percolate and in which directions
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Percolation.h"


int Percolation::Determine_percolating_clusters(const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const vector<vector<int> > &boundary_cnt, const vector<int> &labels, const vector<int>  &labels_labels, const vector<int>  &label_map, vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated, const int &window)
{
    time_t ct0, ct1;
    //---------------------------------------------------------------------------------------------------------------------------
    //Check if the clusters of CNTs percolate
    hout << "Cluster percolation ";
    ct0 = time(NULL);
    if (!Cluster_CNT_percolation(boundary_cnt, labels, labels_labels, label_map, clusters_cnt, isolated)) {
        hout << "Error in Determine_percolating_clusters" << endl;
        return 0;
    }
    ct1 = time(NULL);
    hout << "time: " << (int)(ct1-ct0) <<" secs." << endl;
  
    //---------------------------------------------------------------------------------------------------------------------------
    //Check if the (single) CNTs that do not form clusters might percolate
    hout << "CNT percolation ";
    ct0 = time(NULL);
    if (!Single_CNT_percolation(boundary_cnt, labels, labels_labels, clusters_cnt, isolated, sample, cnts, window)) {
        hout << "Error in Single_CNT_percolation" << endl;
        return 0;
    }
    ct1 = time(NULL);
    hout << "time: " << (int)(ct1-ct0) <<" secs." << endl;
    
    //---------------------------------------------------------------------------------------------------------------------------
    hout << "percolated clusters = "<<  clusters_cnt.size() << endl;
    //Just a check. family and clusters_cnt MUST have the same size
    //hout << "Check"<<endl;
    if (family.size() != clusters_cnt.size()) {
        hout <<  "ERROR the vector of family number and clusters_cnt do not have the same size" << endl;
        hout << "fam size = " << family.size()<< "\t cluster size = "<<  clusters_cnt.size() << endl;
        return 0;
    }
    //Print2DVec(clusters_cnt_iso, "isolated.txt");
    return 1;
}

//This function determines if the clusters in the clusters_cnt variable actually percolate
int Percolation::Cluster_CNT_percolation(const vector<vector<int> > &boundary_cnt, const vector<int> &labels, const vector<int>  &labels_labels, const vector<int>  &label_map, vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated)
{
    //Check if there is any cluster at all. If the size of clusters_cnt is non-zero, then there are clusters
    if (clusters_cnt.size()) {
        //Initialize the flag vectors
        vector<short int> zeros(6,0);
        vector<vector<short int> > perc_flag(clusters_cnt.size(), zeros);
        //hout << "0 ";
        //Scan all the boundary vectors and fill the flag vectors
        Fill_percolation_flags_all_directions(boundary_cnt, perc_flag, labels, labels_labels, label_map);
        //hout << "1 ";
        
        //Move all non-percolating clusters to the corresponding vector<vector>
        if (!Check_percolation_all_clusters(boundary_cnt, perc_flag, clusters_cnt, isolated)) {
            hout << "Error in Determine_percolating_clusters >> Check_percolation_all_clusters" << endl;
            return 0;
        }
    }
    return 1;
}

//This function calls the Fill_percolation_flags_single_direction and usese three flags: px, py and pz.
//The flags are used in case one direction is not necessary to check for percolation. This is specially useful for single CNT percolation
void Percolation::Fill_percolation_flags_all_directions(const vector<vector<int> > &boundary_cnt, vector<vector<short int> > &perc_flag, const vector<int> &labels, const vector<int>  &labels_labels, const vector<int>  &label_map)
{
    //Scan all the boundary vectors and fill the flag vectors
    //X-direction
    Fill_percolation_flags_single_direction(boundary_cnt[0], 0, perc_flag, labels, labels_labels, label_map);
    Fill_percolation_flags_single_direction(boundary_cnt[1], 1, perc_flag, labels, labels_labels, label_map);
    //Y-direction
    Fill_percolation_flags_single_direction(boundary_cnt[2], 2, perc_flag, labels, labels_labels, label_map);
    Fill_percolation_flags_single_direction(boundary_cnt[3], 3, perc_flag, labels, labels_labels, label_map);
    //Z-direction
    Fill_percolation_flags_single_direction(boundary_cnt[4], 4, perc_flag, labels, labels_labels, label_map);
    Fill_percolation_flags_single_direction(boundary_cnt[5], 5, perc_flag, labels, labels_labels, label_map);
}

//This function fills the percolation_flags vector, which is usde to determine percolation
void Percolation::Fill_percolation_flags_single_direction(const vector<int> &boundary_vector, int boundary_number, vector<vector<short int> > &perc_flag, const vector<int> &labels, const vector<int>  &labels_labels, const vector<int>  &label_map)
{
    //Variables
    int L, root, CNT, cluster;
    
    for (int i = 0; i < (int)boundary_vector.size(); i++) {
        //Current CNT
        CNT = boundary_vector[i];
        //hout << "CNT=" << CNT << ' ';
        //Label of the CNT
        L = labels[CNT];
        //If the label is -1, then the CNT does not belong ot any cluster
        if (L != -1) {
            //Cluster number (root) of the CNT
            root = Find_root(L, labels_labels);
            //hout << "root=" << root << " L=" << L <<' ';
            //Find the proper cluster number by using the label map
            cluster = label_map[root];
            //hout << "cluster=" << cluster << ' ';
            //If the CNT is at a boundary but does not belong to any cluster, then it will map to a -1
            //Son only when the CNT maps to something different to -1 then it belongs to a cluster
            if (cluster != -1)
                //Turn on the flag corresponding to the boundary
                perc_flag[cluster][boundary_number] = 1;
        }
    }
}

//Find the root, i.e. the proper cluster number
int Percolation::Find_root(int L, const vector<int> &labels_labels)
{    
    while (labels_labels[L] <= 0){
        L = -labels_labels[L];
        //If labels_labels[L] = 0, then the root is zero, not necesarily L
        if (labels_labels[L] == 0)
            return 0;
    }
    
    return L;
}

int Percolation::Check_percolation_all_clusters(const vector<vector<int> > &boundary_cnt, vector<vector<short int> > &perc_flag, vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated)
{
    //Move all non-percolating clusters to the corresponding vector<vector>
    //I start from the end of the vector to avoid issues with the index when removing an element
    //This variable will store the family number given by Check_percolation_single_cluster
    int fam;
    for (int i = (int)clusters_cnt.size() - 1; i >= 0 ; i--) {
        //hout <<"Check_clusters_percolation " << i << " size " << clusters_cnt.size() << endl;
        //hout << "percolation_flags.size()=" << percolation_flags.size() << endl;
        if (Check_percolation_single_cluster(perc_flag[i], fam)){
            //Add the family number to the vector of families
            family.insert(family.begin(), fam);
        } else {
            //percolation_flags and clusters_cnt have the same size
            //hout << "NO\n";
            isolated.push_back(clusters_cnt[i]);
            //remove the non-percolating cluster
            clusters_cnt.erase(clusters_cnt.begin()+i);
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
int Percolation::Single_CNT_percolation(const vector<vector<int> > &boundary_cnt, const vector<int> &labels, const vector<int>  &labels_labels,vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated, const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const int &window)
{
    //These are variables for the geometry of the observation window
    //Dimensions of the current observation window
    double w_x = sample.win_max_x - window*sample.win_delt_x;
    double w_y = sample.win_max_y - window*sample.win_delt_y;
    double w_z = sample.win_max_z - window*sample.win_delt_z;
    //These variables are the coordinates of the lower corner of the observation window
    //double xmin = sample.origin.x + (sample.len_x - w_x)/2;
    //double ymin = sample.origin.y + (sample.wid_y - w_y)/2;
    //double zmin = sample.origin.z + (sample.hei_z - w_z)/2;
    
    //These variables are flags to determine in which directions there might be percolation by a single CNT
    int px = 0, py= 0, pz = 0;
    if (cnts.len_max >= w_x)
        px = 1;
    if (cnts.len_max >= w_y)
        py = 1;
    if (cnts.len_max >= w_z)
        pz = 1;
    
    //hout << "0 ";
    if (px+py+pz == 0) {
        //If in all three directions the CNTs are larger than the dimensions of the observation window, then there is nothing more
        //to do in this function. This happens when px, px and py are zero, so the only way the sum is zero is when all three are zero.
        //Then, I check the sum instead of the three conditions that each are zero.
        return 1;
    }
    
    //This is similar to the labels generated by HK76. Since it is only needed when the CNTs are smaller than a dimension of
    //the observation window, it is created here and not during the HK76
    vector<int>  labels_iso;
    //The labels vector has the same size as the number of CNTs, so I use it to initialize labels_iso
    labels_iso.assign(labels.size(), -1);
    
    /*/hout << "1 ";
    hout << "isolated.size()="<<isolated.size()<<endl;
    for (int i = 0; i < (int)isolated.size(); i++) {
        hout << "isolated["<<i<<"].size()="<<isolated[i].size()<<endl;
    }//*/
    //Number of single CNTs
    int single_CNTs = 0;
    //Fill label_map_iso
    for (int i = 0; i < (int)isolated.size(); i++) {
        //This function is called after the clusters are checked for percolation, so in the isolated vector there are single CNTs and clusters
        //As I am intereste in single CNTs (isolated[i].size() == 1), when I find clusters (isolated[i].size() > 1) I break the loop
        if (isolated[i].size() > 1)
            break;
        //hout << "isolated["<<i<<"].size()="<<isolated[i].size()<<endl;
        int CNT = isolated[i].front();
        labels_iso[CNT] = i;
        single_CNTs = i+1;
    }
    //hout << "2 ";
    
    //Check if there is any cluster at all. If single_CNTs is non-zero, then there are clusters
    if (single_CNTs) {
        //Initialize the flag vectors
        vector<short int> zeros(6,0);
        vector<vector<short int> > perc_flag(single_CNTs, zeros);
        //hout << "10 ";
        //Scan all the boundary vectors and fill the flag vectors
        Fill_percolation_flags_all_directions_CNT_percoaltion(boundary_cnt, perc_flag, labels_iso, px, py, pz);
        //hout << "11 ";
        
        //The percolated CNTs will be stored in the clusters_tmp vector
        vector<vector<int> > clusters_tmp;
        if (!Check_percolation_CNTs(boundary_cnt, perc_flag, clusters_tmp, isolated)) {
            hout << "Error in Single_CNT_percolation" << endl;
            return 0;
        }
        hout << "Single CNT clusters="<<clusters_tmp.size()<<endl;
        
        //Add the percolated CNTs to clusters_cnt vector
        for (int i = 0; i < (int)clusters_tmp.size(); i++) {
            clusters_cnt.push_back(clusters_tmp[i]);
        }
        //hout << "13 ";
        
    }
    
    return 1;
}
//This function calls the Fill_percolation_flags_single_direction and usese three flags: px, py and pz.
//The flags are used in case one direction is not necessary to check for percolation. This is specially useful for single CNT percolation
void Percolation::Fill_percolation_flags_all_directions_CNT_percoaltion(const vector<vector<int> > &boundary_cnt, vector<vector<short int> > &perc_flag, const vector<int> &labels_iso, int px, int py, int pz)
{
    //Scan all the boundary vectors and fill the flag vectors
    //X-direction
    if (px) {
        Fill_percolation_flags_single_direction_CNT_percoaltion(boundary_cnt[0], 0, perc_flag, labels_iso);
        Fill_percolation_flags_single_direction_CNT_percoaltion(boundary_cnt[1], 1, perc_flag, labels_iso);

    }
    //Y-direction
    if (py) {
        Fill_percolation_flags_single_direction_CNT_percoaltion(boundary_cnt[2], 2, perc_flag, labels_iso);
        Fill_percolation_flags_single_direction_CNT_percoaltion(boundary_cnt[3], 3, perc_flag, labels_iso);
    }
    //Z-direction
    if (pz) {
        Fill_percolation_flags_single_direction_CNT_percoaltion(boundary_cnt[4], 4, perc_flag, labels_iso);
        Fill_percolation_flags_single_direction_CNT_percoaltion(boundary_cnt[5], 5, perc_flag, labels_iso);
    }
}

//This function fills the percolation_flags vector, which is usde to determine percolation
void Percolation::Fill_percolation_flags_single_direction_CNT_percoaltion(const vector<int> &boundary_vector, int boundary_number, vector<vector<short int> > &perc_flag, const vector<int> &labels_iso)
{
    //Variables
    int L, CNT;
    
    for (int i = 0; i < (int)boundary_vector.size(); i++) {
        //Current CNT
        CNT = boundary_vector[i];
        //hout << "CNT=" << CNT << ' ';
        //Label of the CNT, which also corresponds to its perc_flag number
        L = labels_iso[CNT];
        //If the label is -1, then the CNT is an isolated cluster and need to check if the percolation flags need to be activated
        if (L != -1) {
            //Activate the flag corresponding to the boundary
            perc_flag[L][boundary_number] = 1;
        }
    }
}

int Percolation::Check_percolation_CNTs(const vector<vector<int> > &boundary_cnt, vector<vector<short int> > &perc_flag, vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated)
{
    //Check if every cluster made of a single cluster percolates
    //I start from the end of the vector of percolated flags to avoid issues with the index when removing an element
    //This variable (fam) will store the family number given by Check_percolation_single_direction
    int fam;
    for (int i = (int)perc_flag.size() - 1; i >= 0 ; i--) {
        //hout <<"Check_clusters_percolation " << i << " isolated.size " << isolated.size() << endl;
        //hout << "perc_flag.size()=" << perc_flag.size() << endl;
        //If there is percolation, then add the isolated CNT to the clusters_cnt vector
        if (Check_percolation_single_CNT(perc_flag[i], fam)){
            clusters_cnt.push_back(isolated[i]);
            //remove the non-percolating cluster
            isolated.erase(isolated.begin()+i);
            //Add the family number to the vector of families
            family.insert(family.begin(), fam);
        }
    }
    return 1;
}

//This function will check if there is percolation for a single cluster consisting of a single CNT in x, y or z directions
int Percolation::Check_percolation_single_CNT(vector<short int> cluster_flag, int &family)
{
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

