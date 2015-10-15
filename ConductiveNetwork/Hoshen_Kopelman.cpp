//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Hoshen_Kopelman.h
//OBJECTIVE:	The Hoshen_Kopelman Algorithm
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Hoshen_Kopelman.h"

/*
 
 This function implements the Hoshen-Kopelman algorithm and generates the vector of points that are in contact.
 
 There is the function Scan_sub_regions_then_delete that it is not used because it is too slow. I leave it here as a remark. That function first generates the 2D vector of contact points, then scans every 1D vector of contacts and deletes repeated number. This is way too slower than first checking if there are repeated contacts (with the funtion Check_repeated) and then adding the contact if it is a new one.

 
 Input:
    struct Cutoff_dist cutoffs
        Structure that contains the cutoff for tunneling
    vector<int> cnts_inside
        List of CNTs that are inside the observation window
    vector<vector<long int> > sectioned_domain
        List with all points grouped into sub-regions in order to reduce computacional cost of finding contacts
    vector<vector<long int> > structure
        Vector with the structure
    vector<Point_3D> points_in
        List of points
    vector<double> radii
        List of radii. Using this vector allows for the code to be able to work with CNTs of different radii
 
 Output (class variables):
    vector<vector<long int> > contacts_point
        Vector of point to point contacts. This is helpful for determining the resistor network on the direct electrifying algorithm
    vector<vector<int> > clusters_cnt
        Vector with clusters of CNT. The size of clusters_cnt is the number of clusters
    vector<vector<int> > isolated
        Vector with CNTs tha are isolated, i.e. form a cluster of 1 CNT. Each isolated[i] is a cluster of size 1. Later, non percolated clusters found in vector clusters_cnt are moved to isolated
 
 It also creates a connectivity vector of the points in contact. This connectivity vector is used to define the elements used in the direct electrifying algorithm.
 */

//To determinate nanotube clusters using Hoshen Kopelman Algorithm
int Hoshen_Kopelman::Determine_nanotube_clusters(const struct Cutoff_dist &cutoffs, const vector<int> &cnts_inside, const vector<vector<long int> > &sectioned_domain, const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const vector<double> &radii)
{
    //Label the CNTs and make the data structures for the direct electrifying algorithm
   if (!Scan_sub_regions(points_in, radii, cutoffs.tunneling_dist, sectioned_domain)){
        hout << "Error in Determinate_nanotube_clusters." <<endl;
        return 0;
    }
    
    //Make the clusters
    if (!Make_CNT_clusters(structure, points_in, cnts_inside)){
        hout << "Error in Determinate_nanotube_clusters." <<endl;
        return 0;
    }
    
	return 1;
}

//This function scans all the subregions to look for points close enough for tunneling to happen, i.e. points that are in contact.
//When points are in contact use the Hoshen-Kopelman algorithm
int Hoshen_Kopelman::Scan_sub_regions(const vector<Point_3D> &points_in, const vector<double> &radii, const double &tunnel, const vector<vector<long int> > &sectioned_domain)
{
    //These ints are just to store the global point number and the CNTs they belong to.
    //They are just intermediate variables and I only use them to make the code more readable
    long int P1, P2;
    int CNT1, CNT2;
    //Temporary vector of int's to store the contact pair
    vector<long int> empty;
    vector<int> empty_int;
    //the list of contacts has to be the same size as the list of points
    contacts_point.assign(points_in.size(), empty);
    //Varable to calculate the cutoff for tunneling
    double cutoff_t;
    //double separation;
    //hout << "tunnel=" << tunnel << endl;
    
    //Initialize the variables for the labels. The size of the vector labels has to be equal to the number of CNTs
    //It is initialized to -1 so if there is a bug in the code, there is going to be an error when using the -1 as an index
    labels.assign(radii.size(), -1);
    
    //Variable for an inner loop
    long int inner;
    
    for (long int i = 0; i < sectioned_domain.size(); i++) {
        inner = sectioned_domain[i].size();
        for (long int j = 0; j < inner-1; j++) {
            P1 = sectioned_domain[i][j];
            CNT1 = points_in[P1].flag;
            for (long int k = j+1; k<inner; k++) {
                P2 = sectioned_domain[i][k];
                CNT2 = points_in[P2].flag;
                //If distance below the cutoff and points belong to different CNT
                cutoff_t = radii[CNT1] + radii[CNT2] + tunnel;
                //hout <<"P1="<<P1<<" CNT1="<<CNT1<<" P2="<<P2<<" CNT2="<<CNT2;
                //hout <<" cutoff_t="<<cutoff_t;
                //hout <<" CNT1="<<CNT1<<" CNT2="<<CNT2;
                //hout<<" r1="<<radii[CNT1]<<" r2="<<radii[CNT2]<<endl;
                //First check if the CNTs are different. Only when the CNTs are different the distance between points is calculated
                //In this way calculation of all distances is avoided
                if ((CNT1!=CNT2)&&(points_in[P1].distance_to(points_in[P2]) <= cutoff_t)) {
                    if (!Check_repeated(contacts_point[P1], P2)) {
                        //Fill the vector of contacts contacts_point
                        contacts_point[P2].push_back(P1);
                        contacts_point[P1].push_back(P2);
                        //Add point contacts in CNT1
                    }
                    
                    //Here is where the actual HK76 algotihm takes place
                    if (!HK76(CNT1, CNT2)) {
                        hout << "Error in Contacts_and_HK76" << endl;
                        return 0;
                    }
                }
            }
        }
    }
    
    return 1;
}
//This function checks if the point Point is in the vector region
int Hoshen_Kopelman::Check_repeated(const vector<long int> &region, long int Point)
{
    for (long int i = 0; i < region.size(); i++) {
        if (Point == region[i]) {
            return 1;
        }
    }
    return 0;
}

//This function scans all the subregions to look for points close enough for tunneling to happen, i.e. points that are in contact.
//When points are in contact use the Hoshen-Kopelman algorithm
int Hoshen_Kopelman::Scan_sub_regions_then_delete(const vector<Point_3D> &points_in, const vector<double> &radii, const double &tunnel, const vector<vector<long int> > &sectioned_domain)
{
    //These ints are just to store the global point number and the CNTs they belong to.
    //They are just intermediate variables and I only use them to make the code more readable
    long int P1, P2;
    int CNT1, CNT2;
    //Temporary vector of int's to store the contact pair
    vector<long int> empty;
    vector<int> empty_int;
    //the list of contacts has to be the same size as the list of points
    contacts_point.assign(points_in.size(), empty);
    //Varable to calculate the cutoff for tunneling
    double cutoff_t;
    //double separation;
    //hout << "tunnel=" << tunnel << endl;
    
    //Initialize the variables for the labels. The size of the vector labels has to be equal to the number of CNTs
    //It is initialized to -1 so if there is a bug in the code, there is going to be an error when using the -1 as an index
    labels.assign(radii.size(), -1);
    
    //Variable for an inner loop
    long int inner;
    
    for (long int i = 0; i < sectioned_domain.size(); i++) {
        inner = sectioned_domain[i].size();
        for (long int j = 0; j < inner-1; j++) {
            P1 = sectioned_domain[i][j];
            CNT1 = points_in[P1].flag;
            for (long int k = j+1; k<inner; k++) {
                P2 = sectioned_domain[i][k];
                CNT2 = points_in[P2].flag;
                //If distance below the cutoff and points belong to different CNT
                cutoff_t = radii[CNT1] + radii[CNT2] + tunnel;
                //hout <<"P1="<<P1<<" CNT1="<<CNT1<<" P2="<<P2<<" CNT2="<<CNT2;
                //hout <<" cutoff_t="<<cutoff_t;
                //hout <<" CNT1="<<CNT1<<" CNT2="<<CNT2;
                //hout<<" r1="<<radii[CNT1]<<" r2="<<radii[CNT2]<<endl;
                //First check if the CNTs are different. Only when the CNTs are different the distance between points is calculated
                //In this way calculation of all distances is avoided
                if ((CNT1!=CNT2)&&(points_in[P1].distance_to(points_in[P2]) <= cutoff_t)) {
                    //Fill the vector of contacts contacts_point 
                    contacts_point[P2].push_back(P1);
                    contacts_point[P1].push_back(P2);
                    //Here is where the actual HK76 algotihm takes place
                    if (!HK76(CNT1, CNT2)) {
                        hout << "Error in Contacts_and_HK76" << endl;
                        return 0;                       
                    }
                }
            }
        }
    }
    
    //Delete repeated contacts
    hout << " Delete repeated contacts ";
    Delete_repeated_contacts();

    return 1;
}

//Function for the Hoshen-Kopelman (HK76) algorithm only
//It is assumed that this function is used only when a contact is found. Otherwise the results will be wrong and probably one cluster with all CNTs will be generated
int Hoshen_Kopelman::HK76(int CNT1, int CNT2) {
    //L is just a label variable and new_label will take the value of the newest cluster
    int L, new_label = 0;
    
    //If both labels of a CNT are -1, then both CNT will form a new cluster labeled with the current value in new_label
    //If only one of the CNTs has label -1, then use the same label as the other CNT
    //If CNTs have different label and both are not -1, then the labels need to merge
    //labels_labels[i] keeps the size of the label i. If two labels are merged, then labels_labels[i] will have a negative vaule equal to the negative of the smallest merged label
    if ( (labels[CNT1] == -1) && (labels[CNT2] == -1) ) {
        labels[CNT1] = new_label;
        labels[CNT2] = new_label;
        new_label++;
        labels_labels.push_back(2);
    } else if (labels[CNT1] == -1) {
        //the proper label is in labels[CNT2]
        L = labels[CNT2];
        if (labels_labels[L] <= 0) {
            //hout << "L=" << L << " LL[L]=" << labels_labels[L] << endl;
            L = Find_root(L);
            //hout << "L_root=" << L << " LL[L_root]=" << labels_labels[L] << endl;
        }
        labels[CNT1] = L;
        labels_labels[L] = labels_labels[L] + 1;
    } else if (labels[CNT2] == -1) {
        //the proper label is in labels[CNT1]
        L = labels[CNT1];
        if (labels_labels[L] <= 0) {
            //hout << "L=" << L << " LL[L]=" << labels_labels[L] << endl;
            L = Find_root(L);
            //hout << "L_root=" << L << " LL[L_root]=" << labels_labels[L] << endl;
        }
        labels[CNT2] = L;
        labels_labels[L]= labels_labels[L] + 1;
    } else if (labels[CNT1] != labels[CNT2]) {
        //Solve label "conflict" and merge labels
        if (!Merge_labels(Find_root(labels[CNT1]), Find_root(labels[CNT2]))){
            hout << "Error in HK76" << endl;
            return 0;
        }
    }
    return 1;
}

//Find in the label L is a root or it points to another label
//If it points to another label, find that label (i.e. the root)
int Hoshen_Kopelman::Find_root(int L) {
    
    while (labels_labels[L] <= 0){
        L = -labels_labels[L];
        //If labels_labels[L] = 0, then the root is zero, not necesarily L
        if (labels_labels[L] == 0)
            return 0;
    }
    
    return L;
}

//this function merges two clusters
int Hoshen_Kopelman::Merge_labels(int root1, int root2){
    if ( (labels_labels[root1] <= 0) || (labels_labels[root2] <= 0) ) {
        hout << "Error on merging clusters. Both labels are negative, at this point this should not happen."<<endl;
        hout << "root1=" << root1 << " LL[root1]=" << labels_labels[root1];
        hout << " root2=" << root2 << " LL[root2]=" << labels_labels[root2]<< endl;
        return 0;
    }
    if (root1 < root2) {
        //In this case the root is root1
        labels_labels[root1] = labels_labels[root1] + labels_labels[root2];
        labels_labels[root2] = -root1;
    } else if (root2 < root1) {
        //In this case the root is root2;
        labels_labels[root2] = labels_labels[root2] + labels_labels[root1];
        labels_labels[root1] = -root2;
    }
    //If root1 is equal to root2, the CNTs are already in the same cluster so there is nothing to do
    
    return 1;
}

//In this function all the contacts that are repeated due to the overlapping of subregions are discarded
void Hoshen_Kopelman::Delete_repeated_contacts()
{
    for (int i = 0; i < (int)contacts_point.size(); i++)
        Discard_repeated(contacts_point[i]);
}

//This Functions deletes the repeated items in the vector vec
void Hoshen_Kopelman::Discard_repeated(vector<long int> &vec)
{
    for (long int i = 0; i <  vec.size()-1; i++) {
        for (long int k = i+1; k < vec.size(); k++) {
            if (vec[i] == vec[k]){
                vec.erase(vec.begin() + k);
                k--;
            }
        }
    }
}

//In this function the clusters are made using the labels from the HK76
int Hoshen_Kopelman::Make_CNT_clusters(const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const vector<int> &cnts_inside)
{
    //This vector has a map of labels in order to know which ones are proper labels
    label_map.assign(labels_labels.size(),-1);
    //This variable will be used to count the number of cluster
    int counter = 0;
    
    //First scan all the labels of labels to find the proper labels and create the map of labels
    for (int i = 0; i < (int)labels_labels.size(); i++) {
        //if labels_labels[i] > 0, then i is a proper label
        if (labels_labels[i] > 0) {
            label_map[i] = counter;
            counter++;
        }
    }
    //Assing the correct size to the vector of clusters: the number of clusters is stored in the variable counter
    vector<int> empty;
    clusters_cnt.assign(counter, empty);
    
    //
    int root, L, n_cluster, CNT;
    //Now scan the vector of cnts_inside. Check the label of each CNT to make the clusters.
    //The vectos for the labels of labels and label map are used to group the CNTs
    for (int i = 0; i < (int)cnts_inside.size(); i++) {
        //Current CNT
        CNT = cnts_inside[i];
        //Store the label in the variable
        L = labels[CNT];
        //If a label[i] is -1, it means CNT_i is an isolated CNT. At this point "isolated" means that the CNT is not forming a cluster
        //although it might percolate but only when it's length is greater or equal than that of the observation window
        //Only when label[i] is diffenrent form -1, CNT_i belongs to a cluster
        if (L != -1){
            //Check if the label corresponds to a proper cluster. Otherwise find the root (i.e. proper cluster number)
            if (labels_labels[L] <= 0)
                root = Find_root(L);
            else
                root = L;
            n_cluster = label_map[root];
            //If the cluster number is negative, then there is an error with the label_map vector
            if (n_cluster < 0) {
                hout << "There is an error with the label_map vector. n_cluster = " << n_cluster ;
                hout << ", root=" << root << " LL[root]=" << labels_labels[root] << " L=" << L ;
                hout << " LL[L]=" << labels_labels[L] << " proper clusters=" << counter << endl;
                return 0;
            }
            clusters_cnt[n_cluster].push_back(CNT);
        }
    }
    
    return 1;
}


//-------------------------------------------------------------------------------------------------------------------------------------
