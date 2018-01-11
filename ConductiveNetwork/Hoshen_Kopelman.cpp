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
 
 Input:
    const struct Geom_RVE sample
        Sample information and geometry
    struct Cutoff_dist cutoffs
        Structure that contains the cutoff for tunneling
    vector<int> cnts_inside
        List of CNTs that are inside the observation window
    vector<vector<long int> > sectioned_domain
        List with all points grouped into sub-regions in order to reduce computacional cost of finding contacts
    vector<vector<long int> > structure
        Vector with the structure
    vector<Point_3D> points_in
        List of CNT points
    vector<double> radii
        List of radii. Using this vector allows for the code to be able to work with CNTs of different radii
    const vector<int> gnps_inside
        List of GNPs that are inside the observation window
    const vector<vector<long int> > sectioned_domain_gnp
        List with all points grouped into sub-regions in order to reduce computacional cost of finding contacts
    const vector<vector<long int> > structure_gnp
        Vector with the structure
    const vector<Point_3D> points_gnp
        List of GNP points
    const vector<GCH> hybrid_particles
        List of hybrid particles or GNPs
 
 Output (class variables):
    vector<int> labels
        Labels for HK76 (CNTs)
    vector<int> labels_labels
        Labels of labels for HK76 (CNTs)
    vector<int> labels_gnp
        Labels for HK76 (GNPs)
    vector<int> labels_labels_gnp
        Labels of labels for HK76 (GNPs)
    vector<vector<long int> > contacts_point
        Vector of point to point contacts. This is helpful for determining the resistor network on the direct electrifying algorithm
    vector<contact_pair> gnp_contacts
        List of GNP to GNP contacts
    vector<contact_pair> mixed_contacts
        List of CNT to GNP contacts
    vector<vector<int> > clusters_cnt
        Vector with clusters of CNT. The size of clusters_cnt is the number of clusters
    vector<vector<int> > isolated
        Vector with CNTs tha are isolated, i.e. form a cluster of 1 CNT. Each isolated[i] is a cluster of size 1. Later, non percolated clusters found in vector clusters_cnt are moved to isolated
    vector<vector<int> > clusters_gch
        Vector with clusters of GNP. The size of clusters_gch is the number of clusters
    vector<vector<int> > isolated_gch
        Vector with GNPs tha are isolated, i.e. form a cluster of 1 GNP. Each isolated_gch[i] is a cluster of size 1. Later, non percolated clusters found in vector clusters_gch are moved to isolated
 
 It also creates a connectivity vector of the points in contact. This connectivity vector is used to define the elements used in the direct electrifying algorithm.
 */

//To determinate hybrid particle clusters using Hoshen Kopelman Algorithm
int Hoshen_Kopelman::Determine_clusters(const struct Geom_RVE &sample, const struct Cutoff_dist &cutoffs, const vector<int> &cnts_inside, const vector<vector<long int> > &sectioned_domain, const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const vector<double> &radii, const vector<int> &gnps_inside, const vector<vector<long int> > &sectioned_domain_gnp, const vector<vector<int> > &sectioned_domain_hyb, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_gnp, const vector<GCH> &hybrid_particles)
{
    //There are three distinct cases to make clusters: CNTs, GNPs or having both (either hybrid or mixed
    //============================================================================================
    if (sample.particle_type == "CNT_wires") {
        //Label the CNTs and make the data structures for the direct electrifying algorithm
        if (!Scan_sub_regions_cnt(sample, points_in, gnps_inside, hybrid_particles, radii, cutoffs.tunneling_dist, sectioned_domain, structure)){
            hout << "Error in Determine_clusters when calling Scan_sub_regions_cnt." <<endl;
            return 0;
        }
        //Make the CNT clusters
        int n_clusters = (int)labels_labels.size();
        if (!Make_particle_clusters(n_clusters, cnts_inside, labels, isolated, clusters_cnt)){
            hout << "Error in Determine_clusters when calling Make_particle_clusters for CNT clusters." <<endl;
            return 0;
        }
        
    }
    //============================================================================================
    else if (sample.particle_type == "GNP_cuboids") {
        //Label the GNPs and make the data structures for the direct electrifying algorithm
        if (!Scan_sub_regions_gnp(points_gnp, hybrid_particles, cutoffs.tunneling_dist, sectioned_domain_gnp)){
            hout << "Error in Determine_clusters when calling Scan_sub_regions_gnp." <<endl;
            return 0;
        }
        //Make the GNP clusters
        int n_clusters = (int)labels_labels_gnp.size();
        if (!Make_particle_clusters(n_clusters, gnps_inside, labels_gnp, isolated_gch, clusters_gch)){
            hout << "Error in Determine_clusters when calling Make_particle_clusters for GNP clusters." <<endl;
            return 0;
        }
    }
    //============================================================================================
    else {
        
        //Time markers
        time_t ct0, ct1;
        
        ct0 = time(NULL);
        //Label the CNTs and make the data structures for the direct electrifying algorithm
        if (!Scan_sub_regions_cnt(sample, points_in, gnps_inside, hybrid_particles, radii, cutoffs.tunneling_dist, sectioned_domain, structure)){
            hout << "Error in Determine_clusters when calling Scan_sub_regions_cnt." <<endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Scan_sub_regions_cnt: " << (int)(ct1-ct0) <<" secs." << endl;
        
        ct0 = time(NULL);
        //Label the GNPs and make the data structures for the direct electrifying algorithm
        if (!Scan_sub_regions_gnp(points_gnp, hybrid_particles, cutoffs.tunneling_dist, sectioned_domain_gnp)){
            hout << "Error in Determine_clusters when calling Scan_sub_regions_gnp." <<endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Scan_sub_regions_gnp: " << (int)(ct1-ct0) <<" secs." << endl;
        
        //Vectors for the HK76 for clustering the CNT and GNP clusters
        vector<int> labels_mixed;
        vector<int> labels_labels_mixed;
        ct0 = time(NULL);
        //Search for mixed contacts
        if (!Scan_sub_regions_cnt_and_gnp(sample, cutoffs.tunneling_dist, points_in, radii, sectioned_domain, structure, points_gnp, hybrid_particles, gnps_inside, sectioned_domain_gnp, sectioned_domain_hyb, labels_mixed, labels_labels_mixed)) {
            hout << "Error in Determine_clusters when calling Scan_sub_regions_cnt_and_gnp" << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Scan_sub_regions_cnt_and_gnp: " << (int)(ct1-ct0) <<" secs." << endl;
        
        //Make the CNT clusters
        int n_clusters = (int)labels_labels_mixed.size();
        ct0 = time(NULL);
        if (!Make_particle_clusters(n_clusters, cnts_inside, labels, isolated, clusters_cnt)){
            hout << "Error in Determine_clusters when calling Make_particle_clusters for CNT clusters." <<endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Scan_sub_regions_cnt: " << (int)(ct1-ct0) <<" secs." << endl;
        
        ct0 = time(NULL);
        //Make the GNP clusters
        if (!Make_particle_clusters(n_clusters, gnps_inside, labels_gnp, isolated_gch, clusters_gch)){
            hout << "Error in Determine_clusters when calling Make_particle_clusters for GNP clusters." <<endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Scan_sub_regions_cnt: " << (int)(ct1-ct0) <<" secs." << endl;
        
    }
    
    return 1;
}
//This function scans all the subregions to look for points close enough for tunneling to happen, i.e. points that are in contact.
//When points are in contact use the Hoshen-Kopelman algorithm
int Hoshen_Kopelman::Scan_sub_regions_cnt(const struct Geom_RVE &sample, const vector<Point_3D> &points_in, const vector<int> &gnps_inside, const vector<GCH> &hybrid_particles, const vector<double> &radii, const double &tunnel, const vector<vector<long int> > &sectioned_domain, const vector<vector<long int> > &structure)
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
    
    //Initialize the variables for the labels. The size of the vector labels has to be equal to the number of CNTs
    //It is initialized to -1 so if there is a bug in the code, there is going to be an error when using the -1 as an index
    labels.assign(radii.size(), -1);
    
    //new_label will take the value of the newest cluster
    int new_label = 0;
    
    //Hybrid particle pre-processing
    //Set the flag for hybrid particles if the particle type is hybrid
    int hybrids_flag = sample.particle_type == "Hybrid_particles";
    //If there are hybrid particles the pre-processing functon is called
    if (hybrids_flag) {
        if (!Group_cnts_in_gnp(hybrid_particles, gnps_inside, new_label)) {
            hout << "Error in Scan_sub_regions when calling Group_cnts_in_gnp" << endl;
            return 0;
        }
    }
    
    //Scan every overlapping sub-region
    for (long int i = 0; i < (long int)sectioned_domain.size(); i++) {
        long int inner = (long int)sectioned_domain[i].size();
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
                //hout<<" r1="<<radii[CNT1]<<" r2="<<radii[CNT2]<<endl;
                //First check if the CNTs are different. Only when the CNTs are different the distance between points is calculated
                //In this way calculation of all distances is avoided
                if ((CNT1!=CNT2)&&(points_in[P1].distance_to(points_in[P2]) <= cutoff_t)) {
                    //If there are hybrid particles and both points are CNT seed points, ignore the contact
                    int ignore_flag = hybrids_flag && (structure[CNT1][0] == P1) && (structure[CNT2][0] == P2);
                    //Check if the contact has already been added
                    if (!ignore_flag && !Check_repeated(contacts_point[P1], P2)) {
                        //Fill the vector of contacts contacts_point
                        contacts_point[P2].push_back(P1);
                        contacts_point[P1].push_back(P2);
                    }
                    
                    //Here is where the actual HK76 algorithm takes place
                    if (!HK76(CNT1, CNT2, new_label, labels, labels_labels)) {
                        hout << "Error in Scan_sub_regions_cnt when calling HK76" << endl;
                        return 0;
                    }
                }
            }
        }
    }
    
    //Clean up the labels to find the proper labels, i.e. merged and consecutive labels starting at 0
    if (!Cleanup_labels(labels, labels_labels)) {
        hout << "Error in Scan_sub_regions_cnt when calling Cleanup_labels" << endl;
        return 0;
    }
    
    return 1;
}
//This function does the pre-processing of the hybrid particles
//A label is assigned to each CNT that belongs to a GNP
//Thus, there will be as many labels as hybrid particles. Then the HK76 will only solve label conflicts
int Hoshen_Kopelman::Group_cnts_in_gnp(const vector<GCH> &hybrid_particles, const vector<int> &gnps_inside, int &new_label)
{
    //Loop through all the hybrid particles inside the sample and then though all CNTs inside that particle and assign them the same label
    for (int i = 0; i < (int)gnps_inside.size(); i++) {
        //CNTs at the top surface
        int cnts_top = (int)hybrid_particles[i].cnts_top.size();
        
        for (int j = 0; j < cnts_top; j++) {
            int CNT = hybrid_particles[i].cnts_top[j];
            //The label to be assign is the iterator "i"
            labels[CNT] = i;
        }
        
        //CNTs at the bottom surface
        int cnts_bottom = (int)hybrid_particles[i].cnts_bottom.size();
        
        for (int j = 0; j < cnts_bottom; j++) {
            int CNT = hybrid_particles[i].cnts_bottom[j];
            //The label to be assign is the iterator "i"
            labels[CNT] = i;
        }
        
        //Update the labels_labels vector with the number of CNTs with label "i"
        labels_labels.push_back(cnts_top+cnts_bottom);
    }
    
    //Update new_label
    new_label = (int)labels_labels.size();
    
    return 1;
}
//This function checks if the point Point is in the vector region
int Hoshen_Kopelman::Check_repeated(const vector<long int> &contacts_vector, const long int &point)
{
    for (int i = 0; i < (int)contacts_vector.size(); i++) {
        if (point == contacts_vector[i]) {
            return 1;
        }
    }
    return 0;
}
//Function for the Hoshen-Kopelman (HK76) algorithm only
//It is assumed that this function is used only when a contact is found. Otherwise the results will be wrong and probably one cluster with all CNTs will be generated
int Hoshen_Kopelman::HK76(const int &CNT1, const int &CNT2, int &new_label, vector<int> &labels,  vector<int> &labels_labels) {
    //L is just a label variable and new_label will take the value of the newest cluster
    int L;
    
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
            L = Find_root(L, labels_labels);
            //hout << "L_root=" << L << " LL[L_root]=" << labels_labels[L] << endl;
        }
        labels[CNT1] = L;
        labels_labels[L] = labels_labels[L] + 1;
    } else if (labels[CNT2] == -1) {
        //the proper label is in labels[CNT1]
        L = labels[CNT1];
        if (labels_labels[L] <= 0) {
            //hout << "L=" << L << " LL[L]=" << labels_labels[L] << endl;
            L = Find_root(L, labels_labels);
            //hout << "L_root=" << L << " LL[L_root]=" << labels_labels[L] << endl;
        }
        labels[CNT2] = L;
        labels_labels[L]= labels_labels[L] + 1;
    } else if (labels[CNT1] != labels[CNT2]) {
        //Solve label "conflict" and merge labels
        int L1 = Find_root(labels[CNT1], labels_labels);
        int L2 = Find_root(labels[CNT2], labels_labels);
        if (!Merge_labels(L1, L2, labels_labels)){
            hout << "Error in HK76" << endl;
            return 0;
        }
    }
    return 1;
}

//Find in the label L is a root or it points to another label
//If it points to another label, find that label (i.e. the root)
int Hoshen_Kopelman::Find_root(const int &L, vector<int> &labels_labels)
{
    int L_root = L;
    while (labels_labels[L_root] <= 0){
        L_root = -labels_labels[L_root];
        //If labels_labels[L_root] = 0, then the root is zero, not necesarily L
        if (labels_labels[L_root] == 0)
            return 0;
    }
    
    return L_root;
}

//this function merges two clusters
int Hoshen_Kopelman::Merge_labels(const int &root1, const int &root2, vector<int> &labels_labels)
{
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
//This function scans all the subregions to look for points close enough for tunneling to happen, i.e. points that are in contact.
//When points are in contact use the Hoshen-Kopelman algorithm
int Hoshen_Kopelman::Scan_sub_regions_gnp(const vector<Point_3D> &points_gnp, const vector<GCH> &hybrid_particles, const double &tunnel, const vector<vector<long int> > &sectioned_domain_gnp)
{
    //These ints are just to store the global point number and the CNTs they belong to.
    //They are just intermediate variables and I only use them to make the code more readable
    long int P1, P2;
    int GNP1, GNP2;
    //Temporary vector of int's to store the contact pair
    vector<long int> empty;
    vector<int> empty_int;
    
    //Initialize the matrices that save the contacts between GNPs
    vector<vector<long int> > point_matrix;
    vector<vector<double> > distance_matrix;
    
    if (!Initialize_contact_matrices((int)hybrid_particles.size(), point_matrix, distance_matrix)) {
        hout << "Error in Scan_sub_regions_gnp when calling Initialize_contact_matrices" << endl;
        return 0;
    }
    
    //Initialize the variables for the labels. The size of the vector labels has to be equal to the number of CNTs
    //It is initialized to -1 so if there is a bug in the code, there is going to be an error when using the -1 as an index
    labels_gnp.assign(hybrid_particles.size(), -1);
    
    //new_label will take the value of the newest cluster
    int new_label = 0;
    
    //Scan every overlapping sub-region
    for (long int i = 0; i < (long int)sectioned_domain_gnp.size(); i++) {
        long int inner = (long int)sectioned_domain_gnp[i].size();
        for (long int j = 0; j < inner-1; j++) {
            P1 = sectioned_domain_gnp[i][j];
            GNP1 = points_gnp[P1].flag;
            for (long int k = j+1; k<inner; k++) {
                P2 = sectioned_domain_gnp[i][k];
                GNP2 = points_gnp[P2].flag;
                //If distance below the cutoff and points belong to different CNT
                
                //Check if the GNPs are different
                //Only when the GNPs are different it is worth to do the rest of computations
                if (GNP1 != GNP2) {
                    
                    //Map P2 to the local coordinates of GNP1
                    Point_3D demapped = Demap_gnp_point(hybrid_particles[GNP1], points_gnp[P2]);
                    
                    //Check if the point is inside the GNP bounding box that determines contact
                    if ( Judge_point_inside_bounding_box(hybrid_particles[GNP1].gnp, demapped, tunnel) ) {
                        
                        //If the point of GNP2 is inside the GNP1 bounding box, then the GNPs are in contact
                        
                        //Identify the GNP with the largest number, this will be the row
                        //The GNP with the lowest number will be the column
                        int row, column;
                        if (GNP1 < GNP2) {
                            row = GNP2;
                            column = GNP1;
                        } else {
                            row = GNP1;
                            column = GNP2;
                        }
                        
                        //Calculate the distance between the points in contact
                        double distance = points_gnp[P1].squared_distance_to(points_gnp[P2]);
                        
                        //hout <<"P1="<<P1<<" GNP1="<<GNP1<<" P2="<<P2<<" GNP2="<<GNP2;
                        //hout <<" row="<<row<<" col="<<column<<endl;
                        //hout <<"point_matrix.s="<<point_matrix.size()<<endl;
                        //hout <<"point_matrix[GNP1].s="<<point_matrix[GNP1].size()<<endl;
                        //hout <<"point_matrix[GNP2].s="<<point_matrix[GNP2].size()<<endl;
                        
                        //Check if the distance between points is smaller than the one in the distance_matrix
                        if (distance < distance_matrix[row][column]) {
                            
                            //If the distance is smaller than the one in the distance matrix, take this pair of points as the contact
                            point_matrix[GNP1][GNP2] = P1; //Point on GNP1 in contact with GNP2
                            point_matrix[GNP2][GNP1] = P2; //Point on GNP2 in contact with GNP1
                            
                            //Update the squared distance
                            distance_matrix[row][column] = distance;
                        }
                        
                        //Here is where the actual HK76 algorithm takes place
                        if (!HK76(GNP1, GNP2, new_label, labels_gnp, labels_labels_gnp)) {
                            hout << "Error in Scan_sub_regions_gnp when calling HK76" << endl;
                            return 0;
                        }
                    }
                }
            }
        }
    }
    
    //Create the vector of GNP contacts
    if (!Create_vector_of_gnp_contacts(point_matrix)) {
        hout << "Error in Scan_sub_regions_gnp when calling Create_vector_of_gnp_contacts" << endl;
        return 0;
    }
    
    //Clean up the labels to find the proper labels, i.e. merged and consecutive labels starting at 0
    if (!Cleanup_labels(labels_gnp, labels_labels_gnp)) {
        hout << "Error in Scan_sub_regions_gnp when calling Cleanup_labels" << endl;
        return 0;
    }
    return 1;
}
//This function initializes two contact matrices, one will keep the point numbers of the contacts adn the other the distances
//The point_matrix needs to be squared to keep the points of contact on each GNP
//The distance matrix can be lower triangular because distance from P1 to P2 is the same as the ditance from P2 to P1
int Hoshen_Kopelman::Initialize_contact_matrices(const int &n_GNPs, vector<vector<long int> > &point_matrix, vector<vector<double> > &distance_matrix)
{
    //vectors to increase the size of the matrices
    vector<long int> tmp_int(n_GNPs,-1);
    vector<double> tmp_double;
    
    for (int i = 0; i < n_GNPs; i++) {
        //Add the row for GNP i
        point_matrix.push_back(tmp_int);
        distance_matrix.push_back(tmp_double);
        
        //Add one column to the distance_matrix
        tmp_double.push_back(1.0);
    }
    
    return 1;
}
//This function maps a point to the local coordinates of a GNP
//It is asummed the center of coordinates is located at the GNP center
Point_3D Hoshen_Kopelman::Demap_gnp_point(const GCH &hybrid, const Point_3D &point_gnp2)
{
    //Point to store the result
    Point_3D demapped_point;
    
    //First subtract the displacement form the point in GNP2
    Point_3D tmp = point_gnp2 - hybrid.center;
    
    //Multiply the transpose of the rotation matrix by the point obtained in the previous step
    //The operation is done explicitly since a matrix-point multiplication is not defined
    demapped_point.x = hybrid.rotation.element[0][0]*tmp.x + hybrid.rotation.element[1][0]*tmp.y + hybrid.rotation.element[2][0]*tmp.z;
    
    demapped_point.y = hybrid.rotation.element[0][1]*tmp.x + hybrid.rotation.element[1][1]*tmp.y + hybrid.rotation.element[2][1]*tmp.z;
    
    demapped_point.z = hybrid.rotation.element[0][2]*tmp.x + hybrid.rotation.element[1][2]*tmp.y + hybrid.rotation.element[2][2]*tmp.z;

    return demapped_point;
}
//This function judges if a point is inside a GNP bounding box
int Hoshen_Kopelman::Judge_point_inside_bounding_box(const struct cuboid &gnp, const Point_3D &point, const double &extension)
{
    //If the point is outside the function returns 0: any of the operations will return 1; thus the OR wil return 1; then the NOT returns 0
    return !(point.x<-(gnp.len_x/2+extension)||point.x>(gnp.len_x/2+extension)||
            point.y<-(gnp.wid_y/2+extension)||point.y>(gnp.wid_y/2+extension)||
            point.z<-(gnp.hei_z/2+extension)||point.z>(gnp.hei_z/2+extension) );
}
//
int Hoshen_Kopelman::Create_vector_of_gnp_contacts(const vector<vector<long int> > &point_matrix)
{
    //Scan the distance matrix, this will be enough to include all contacts
    //The outer loop iterated over the row numbers
    for (int GNP1 = 0; GNP1 < (int)point_matrix.size(); GNP1++) {
        //The inner loop iterates over column numbers
        for (int GNP2 = 0; GNP2 < GNP1; GNP2++) {
            //Check if there is a contact bewteen GNP1 and GNP2
            if (point_matrix[GNP1][GNP2] != -1) {
                //create a new contact_pair
                struct contact_pair new_contact;
                new_contact.point1 = point_matrix[GNP1][GNP2];
                new_contact.particle1 = GNP1;
                new_contact.type1 = "GNP";
                new_contact.point2 = point_matrix[GNP2][GNP1];
                new_contact.particle2 = GNP2;
                new_contact.type2 = "GNP";
                //Add the contact pair to the vector of contacts
                gnp_contacts.push_back(new_contact);
            }
        }
    }
    return 1;
}
//In this function the clusters are made using the labels from the HK76
int Hoshen_Kopelman::Make_particle_clusters(const int &n_clusters, const vector<int> &particles_inside, vector<int> &labels, vector<vector<int> > &isolated, vector<vector<int> > &clusters_particles)
{
    //Assing the correct size to the vector of clusters: the number of clusters is the same as the number of proper labels, which is the same as the size of the vector labels_labels after clean up
    vector<int> empty;
    clusters_particles.assign(n_clusters, empty);
    
    //Now scan the vector of particles_inside. Check the label of each particle to make the clusters.
    for (int i = 0; i < (int)particles_inside.size(); i++) {
        //Current particle
        int particle = particles_inside[i];
        //Store the label in the variable
        int L = labels[particle];
        //If a label[i] is -1, it means particle_i is an isolated particle. At this point "isolated" means that the particle is not part of any cluster
        //Only when label[i] is diffenrent form -1, particle_i belongs to a cluster
        if (L != -1){
            //Since now the labels are consecutive and starting in 0, the label indicates the cluster number
            clusters_particles[L].push_back(particle);
        } else {
            isolated.push_back(empty);
            isolated.back().push_back(particle);
        }
    }
    
    return 1;
}
int Hoshen_Kopelman::Cleanup_labels(vector<int> &labels, vector<int> &labels_labels)
{
    //This vector has a map of labels in order to know which ones are root labels and renumbers them to a proper label
    vector<int> label_map(labels_labels.size(),-1);
    //This variable will be used to count the number of cluster
    int counter = 0;
    
    //Temporary vector of labels of labels initialized to be equal to the original vector
    vector<int> labels_labels_tmp(labels_labels);
    //Clear the original vector of labels
    labels_labels.clear();
    
    //First scan all the labels of labels to find the root labels
    //Create a map from root label to proper label
    //Proper labels are consecutive labels from 0 to N-1, where N is the number of clusters
    //i.e., these are merged and renumbered labels
    for (int i = 0; i < (int)labels_labels_tmp.size(); i++) {
        //if labels_labels_tmp[i] > 0, then i is a proper label
        if (labels_labels_tmp[i] > 0) {
            label_map[i] = counter;
            //Now the vector of labels of labels has only positive integers
            labels_labels.push_back(labels_labels_tmp[i]);
            counter++;
        }
    }
    
    //Scan all labels and change labels to keep only the proper labels
    for (int i = 0; i < (int)labels.size(); i++) {
        //Check that labels[i] is a valid label
        if (labels[i] != -1) {
            //Find the root label of labels[i]
            int root = Find_root(labels[i], labels_labels_tmp);
            //Find the proper label
            int proper = label_map[root];
            if (proper == -1) {
                hout << "Invalid proper label " << proper << ". Valid labels are positive integers or 0"<<endl;
                return 0;
            }
            //Change the current label by the proper label
            labels[i] = proper;
        }
    }
    
    return 1;
}
//This function scans all the subregions to look for points close enough for tunneling to happen, i.e. points that are in contact.
//When points are in contact use the Hoshen-Kopelman algorithm
int Hoshen_Kopelman::Scan_sub_regions_cnt_and_gnp(const struct Geom_RVE &sample, const double &tunnel, const vector<Point_3D> &points_in, const vector<double> &radii, const vector<vector<long int> > &sectioned_domain, const vector<vector<long int> > &structure, const vector<Point_3D> &points_gnp, const vector<GCH> &hybrid_particles, const vector<int> &gnps_inside, const vector<vector<long int> > &sectioned_domain_gnp, const vector<vector<int> > &sectioned_domain_hyb, vector<int> &labels_mixed, vector<int> &labels_labels_mixed)
{
    //These ints are just to store the global point number and the CNTs they belong to.
    //They are just intermediate variables and I only use them to make the code more readable
    long int P1;
    int CNT1, GNP2;
    
    //Temporary vector of int's to store the contact pair
    vector<long int> empty;
    vector<int> empty_int;
    //Contact matrix that assigns the contact number to a given GNP
    vector<vector<int> > gnp_contact_matrix( hybrid_particles.size(), empty_int);
    
    //Number of CNTs, used to calculate the particle number of GNPs
    int n_cnts = (int)labels.size();
    
    //Initialize mixed labels with the values of CNT and GNP labels
    if (!Initialize_mixed_labels(labels_mixed, labels_labels_mixed)) {
        hout << "Error in Scan_sub_regions_cnt_and_gnp when calling Initialize_mixed_labels" << endl;
        return 0;
    }
    
    //Vector to determine the GNP number of a CNT
    vector<int> cnt_gnp_numbers(radii.size(), -1);
    if (!Fill_cnt_gnp_numbers(hybrid_particles, cnt_gnp_numbers)) {
        hout << "Error in Scan_sub_regions_cnt_and_gnp when calling Fill_cnt_gnp_numbers" << endl;
        return 0;
    }
    
    //new_label will take the value of the newest cluster
    //when having mixed particles, there is a non-zero number of clusters, thus new_label is initialized with the size of labels_labels_mixed,
    //which is the number of CNT clusters plus GNP clusters
    int new_label = (int)labels_labels_mixed.size();
    
    //If hybrid particles are used, cluster together GNPs with their CNTs
    if (sample.particle_type == "Hybrid_particles") {
        if (!Cluster_gnps_and_cnts(hybrid_particles, gnps_inside, labels_mixed, labels_labels_mixed, new_label)) {
            hout << "Error in Scan_sub_regions_cnt_and_gnp when calling Cluster_gnps_and_cnts" << endl;
            return 0;
        }
    }
    
    //hout << "gnp_labels="<<gnp_labels<<" cnt_labels="<<cnt_labels<<endl;
    //hout <<"L_cnt="<<labels.size()<<" LL_cnt="<<labels_labels.size()<<"L_gnp="<<labels_gnp.size()<<" LL_gnp="<<labels_labels_gnp.size()<<endl;
    //Scan every overlapping sub-region and compare each point on CNT sub-regions with the GNP in the same sub-region
    //This is possible because both vectors have the same size
    for (long int i = 0; i < (long int)sectioned_domain.size(); i++) {
        
        //number of elements in sectioned_domain[i]
        long int inner1 = (long int)sectioned_domain[i].size();
        
        //Scan CNT points in subregion[i]
        for (long int j = 0; j < inner1; j++) {
            
            //Current point in the CNT sub-regions
            P1 = sectioned_domain[i][j];
            CNT1 = points_in[P1].flag;
            
            //Number of elements in sectioned_domain_hyb[i]
            long int inner2 = (long int)sectioned_domain_hyb[i].size();
            
            for (long int k = 0; k < inner2; k++) {
                
                //Current GNP in the sub-regions
                GNP2 = sectioned_domain_hyb[i][k];
                
                //First check if the CNT is not attached to the GNP. This eliminates adding resistors between a CNT seed, or close to a CNT seed, and a GNP point
                if (cnt_gnp_numbers[CNT1] != GNP2) {
                    
                    //Map the CNT point P to the local coordinates of the GNP
                    Point_3D demapped = Demap_gnp_point(hybrid_particles[GNP2], points_in[P1]);
                    
                    //Calculate the extension of the GNP bounding box, which is the same as the cutoff for tunneling
                    double extension = tunnel + radii[CNT1];
                    
                    //Check if the CNT point is inside the GNP bounding box that determines contact
                    if ( Judge_point_inside_bounding_box(hybrid_particles[GNP2].gnp, demapped, extension) ) {
                        
                        //If the CNT point is inside the boundaing box, then there us a contact
                        //Check if the contact is repeated or it is equivalent to an existing one
                        //if (!Check_repeated_or_equivalent(P1, P2, points_in, points_gnp, mixed_contacts, gnp_contact_matrix[GNP2])) {
                        if (!Check_repeated_or_equivalent_mixed_contact(P1, mixed_contacts, gnp_contact_matrix[GNP2])) {
                            
                            //If not repeated or equivalent, create a new contact
                            struct contact_pair tmp;
                            tmp.point1 = P1;
                            tmp.particle1 = CNT1;
                            tmp.type1 = "CNT";
                            //Temporaryly store the subregion where the contact was located in the variable point2
                            tmp.point2 = i;
                            tmp.particle2 = GNP2;
                            tmp.type2 = "GNP";
                            
                            //Add contact to vector of mixed contacts
                            mixed_contacts.push_back(tmp);
                            //Add contact number to matrix of GNP contacts
                            gnp_contact_matrix[GNP2].push_back((int)mixed_contacts.size()-1);
                        }
                        
                        //Get the particle number of the GNP for the mixed labels
                        int particle2 = GNP2+n_cnts;
                        
                        //Here is where the actual HK76 algorithm takes place
                        if (!HK76(CNT1, particle2, new_label, labels_mixed, labels_labels_mixed)) {
                            hout << "Error in Scan_sub_regions_cnt_and_gnp when calling HK76" << endl;
                            return 0;
                        }
                    }
                    
                }
            }
        }
    }
    
    //Delete repeated and equivalent mixed contacts
    if (!Remove_equivalent_mixed_contacs(structure, mixed_contacts, gnp_contact_matrix)) {
        hout << "Error in Scan_sub_regions_cnt_and_gnp when calling Remove_equivalent_mixed_contacs" << endl;
        return 0;
    }
    
    //Add the GNP point of the mixed contacts found
    if (!Add_gnp_point_to_contact(sectioned_domain_gnp, points_in, points_gnp, mixed_contacts)) {
        hout << "Error in Scan_sub_regions_cnt_and_gnp when calling Add_gnp_point_to_contact" << endl;
        return 0;
    }
    
    //Cleanup the mixed labels
    if (!Cleanup_labels(labels_mixed, labels_labels_mixed)) {
        hout << "Error in Scan_sub_regions_cnt_and_gnp when calling Cleanup_labels" << endl;
        return 0;
    }
    
    //Merge the two sets of labels
    if (!Merge_interparticle_labels(labels_mixed)) {
        hout << "Error in Scan_sub_regions_cnt_and_gnp when calling Merge_interparticle_labels" << endl;
        return 0;
    }
    
    return 1;
}
int Hoshen_Kopelman::Initialize_mixed_labels(vector<int> &labels_mixed, vector<int> &labels_labels_mixed)
{
    //labels_mixed will have a size equal to the number of CNTs plus GNPs
    int n_cnts = (int)labels.size();
    int n_gnps = (int)labels_gnp.size();
    labels_mixed.assign((n_cnts+n_gnps), -1);
    
    //Assign the CNT labels
    for (int i = 0; i < n_cnts; i++) {
        labels_mixed[i] = labels[i];
    }
    
    //Add the CNT labels of labels
    for (int i = 0; i < (int)labels_labels.size(); i++) {
        labels_labels_mixed.push_back(labels_labels[i]);
    }
    
    //Calculate the number of CNT clusters
    int cnt_clusters = (int)labels_labels.size();
    //Renumber and assign the GNP labels
    for (int j = 0; j < n_gnps; j++) {
        //Since the GNP labes are adjusted, check if the GNP has a valid label
        //If the GNP has the -1 label, then nothing needs to be done as the mixed label already has that value
        if (labels_gnp[j] != -1) {
            labels_gnp[j] = labels_gnp[j] + cnt_clusters;
            labels_mixed[j+n_cnts] = labels_gnp[j];
        }
    }
    
    //Add the GNP labels of labels
    for (int j = 0; j < (int)labels_labels_gnp.size(); j++) {
        labels_labels_mixed.push_back(labels_labels_gnp[j]);
    }
    
    return 1;
}
//Function that fills the GNP numbers for every CNT
//-1 indicates the CNT is not attached to any GNP
int Hoshen_Kopelman::Fill_cnt_gnp_numbers(const vector<GCH> &hybrid_particles, vector<int> &cnt_gnp_numbers) {
    
    //Scan all hybrid particles
    for (int i = 0; i < (int)hybrid_particles.size(); i++) {
        
        //Scan all CNTs on the top surface of current hybrid
        for (int j = 0; j < (int)hybrid_particles[i].cnts_top.size(); j++) {
            
            //Current CNT
            int CNT = hybrid_particles[i].cnts_top[j];
            
            //Update GNP number in vector
            cnt_gnp_numbers[CNT] = i;
        }
        
        //Scan all CNTs on the top surface of current hybrid
        for (int j = 0; j < (int)hybrid_particles[i].cnts_bottom.size(); j++) {
            
            //Current CNT
            int CNT = hybrid_particles[i].cnts_bottom[j];
            
            //Update GNP number in vector
            cnt_gnp_numbers[CNT] = i;
        }
    }
    
    return 1;
}
//function that assigns the same cluster to CNTs and their GNPs only for hybrid particles
int Hoshen_Kopelman::Cluster_gnps_and_cnts(const vector<GCH> &hybrid_particles, const vector<int> &gnps_inside, vector<int> &labels_mixed, vector<int> &labels_labels_mixed, int &new_label)
{
    //Number of CNTs
    int n_cnts = (int)labels.size();
    
    //
    
    //Scan vector of GNPs inside the sample
    for (int i = 0; i < (int)gnps_inside.size(); i++) {
        
        //Current GNP
        int GNP = gnps_inside[i];
        
        //Scan all CNTs on top surface
        for (int j = 0; j < (int)hybrid_particles[GNP].cnts_top.size(); j++) {
            
            //Get CNTnumber
            int CNT = hybrid_particles[GNP].cnts_top[j];
            
            //Get the particle number of the GNP for the mixed labels
            int particle2 = GNP+n_cnts;
            if (!HK76(CNT, particle2, new_label, labels_mixed, labels_labels_mixed)) {
                hout << "Error in Scan_sub_regions_cnt_and_gnp when calling HK76" << endl;
                return 0;
            }
        }
        
        //Scan all CNTs on bottom surface
        for (int j = 0; j < (int)hybrid_particles[GNP].cnts_bottom.size(); j++) {
            
            //Get CNTnumber
            int CNT = hybrid_particles[GNP].cnts_bottom[j];
            
            //Get the particle number of the GNP for the mixed labels
            int particle2 = GNP+n_cnts;
            
            //Here is where the actual HK76 algorithm takes place
            if (!HK76(CNT, particle2, new_label, labels_mixed, labels_labels_mixed)) {
                hout << "Error in Scan_sub_regions_cnt_and_gnp when calling HK76" << endl;
                return 0;
            }
        }
    }
    
    return 1;
}
//This function checks if a contact has alraedy been created
//It is assumed that the CNT is the first particle
int Hoshen_Kopelman::Check_repeated_or_equivalent_mixed_contact(const long int &point_cnt, const vector<contact_pair> &contacts, const vector<int> &gnp_contact_vector)
{
    //Iterate over the contacts of the GNP and check if it the contact has already been created
    for (int i = 0; i < (int)gnp_contact_vector.size(); i++) {
        
        //Current contact
        int cont_pair = gnp_contact_vector[i];
        
        //Check if the CNT point is the same
        if (point_cnt == contacts[cont_pair].point1) {
            
            //The CNt point is the same, thus the contact is repeated or equivalent and terminate the function with 1
            //There is no neeed to check the GNP point, since another contact between this CNT point
            //and another GNP point from the same GNP is possible
            //If I also check that the GNP point is the same, then the contacts would be different i
            return 1;
        }
    }
    
    //The contact was not found, so return 0
    return 0;
}
int Hoshen_Kopelman::Remove_equivalent_mixed_contacs(const vector<vector<long int> > &structure, vector<contact_pair> &contacts, vector<vector<int> > &gnp_contact_matrix)
{
    //Vector that will store all contacts to be deleted
    vector<int> to_delete;
    
    //Scan the contacts that each GNP has
    for (int i = 0; i < (int)gnp_contact_matrix.size(); i++) {
        
        //i is the current GNP
        
        //Group the contacts of the current GNP by CNTs
        vector<vector<int> > grouped_contacts;
        if(!Group_mixed_contacts_by_cnt(structure, contacts, gnp_contact_matrix[i], grouped_contacts)) {
            hout << "Error in Remove_equivalent_mixed_contacs when calling Group_mixed_contacts_by_cnt" <<endl;
            return 0;
        }
        
        //Group into CNT segments, select one contact per segment and delete the rest
        if (!Group_and_merge_consecutive_contacts(gnp_contact_matrix[i], grouped_contacts, contacts, to_delete)) {
            hout << "Error in Remove_equivalent_mixed_contacs when calling Group_and_merge_consecutive_contacts" <<endl;
            return 0;
        }
        
    }
    
    //Sort the vector to_delete
    sort(to_delete.begin(), to_delete.end());
    
    //Delete all contacts indicated by the to_delete vector, starting on the last element (which is the largest index)
    for (int i = (int)to_delete.size()-1; i >=0; i--) {
        contacts.erase(contacts.begin()+to_delete[i]);
    }
    
    return 1;
}
int Hoshen_Kopelman::Add_gnp_point_to_contact(const vector<vector<long int> > &sectioned_domain_gnp, const vector<Point_3D> &points_cnt, const vector<Point_3D> &points_gnp, vector<contact_pair> &contacts)
{
    
    //Scan the mixed contacts found
    for (int i = 0; i < (int)contacts.size(); i++) {
        
        //Get the sub-region number
        long int subregion = contacts[i].point2;
        
        //Get the GNP number
        int GNP = contacts[i].particle2;
        
        //Set the minimum distance equal to 1 micron, which is a lot larger than any contact
        double dist_min = 1;
        
        //Scan the points in the subregion and find the closest to point 1
        for (int j = 0; j < (int)sectioned_domain_gnp[subregion].size(); j++) {
            
            //Get the GNP point in the sectioned domain
            long int P_GNP = sectioned_domain_gnp[subregion][j];
            
            //Check the point is in the GNP
            if (points_gnp[P_GNP].flag == GNP) {
                
                //Get the CNT point
                long int P_CNT = contacts[i].point1;
                
                //If the GNP point is in in the GNP, then check if it is the closes to point1
                //Use squared distance to redue computations
                double dist2 = points_cnt[P_CNT].squared_distance_to(points_gnp[P_GNP]);
                if (dist2 < dist_min) {
                    
                    //Update GNP point in contact
                    contacts[i].point2 = P_GNP;
                    
                    //Update distance
                    dist_min = dist2;
                }
            }
        }
        
    }
    
    
    return 1;
}
//This function creates vectors with equal length as that of a CNT in the structure vector
//but instead of point numbers it has contact numbers
//Thus this vectors can be used to group several consecutive mixed contacts and then merge them into a single mixed contact
int Hoshen_Kopelman::Group_mixed_contacts_by_cnt(const vector<vector<long int> > &structure, const vector<contact_pair> &contacts, const vector<int> &gnp_contact_vector, vector<vector<int> > &grouped_contacts)
//int Hoshen_Kopelman::Group_mixed_contacts_by_cnt(const vector<Point_3D> &points_in, const vector<vector<long int> > &structure, const vector<contact_pair> &contacts, const vector<int> &gnp_contact_vector, vector<vector<int> > &grouped_contacts)
{
    //Mapping vector
    vector<int> cnt_map(structure.size(), -1);
    
    //Counter to map the CNTs
    int counter = 0;
    
    //Scan all the contacts in the current contact vector
    for (int i = 0; i < (int)gnp_contact_vector.size(); i++) {
        
        //get the current contact number
        int cont = gnp_contact_vector[i];
        
        //get the CNT number of contact i
        int CNT = contacts[cont].particle1;
        
        //Check if the CNT has already been mapped
        if (cnt_map[CNT] == -1) {
            
            //Map the CNT
            cnt_map[CNT] = counter;
            
            //Create a vector of -1 with the same size as the CNT in the structure
            //An extra -1 is added to the vector
            //In case the last point has a valid contact, the -1 makes easier to find the end of
            //a segment of consecutive contacts since this ensures all CNT will have a -1 at the end
            vector <int> cnt_contacts(structure[CNT].size()+1,-1);
            
            //Add this vector to the vector of grouped contacts
            grouped_contacts.push_back(cnt_contacts);
            
            //Increase the counter to map the next unmapped CNT
            counter++;
        }
        
        //get the mapped CNT number
        int CNT_mapped = cnt_map[CNT];
        
        //get the CNT point number
        long int P1 = contacts[cont].point1;
        
        //Get the positon on the structure[CNT] vector
        long int pos = P1 - structure[CNT].front();
        
        //Add the contact to the corresponding position on the grouped_contacts vector
        grouped_contacts[CNT_mapped][pos] = cont;
        
    }
    
    return 1;
}
//
int Hoshen_Kopelman::Group_and_merge_consecutive_contacts(const vector<int> &gnp_contact_vector, vector<vector<int> > &grouped_contacts, vector<contact_pair> &contacts, vector<int> &to_delete)
{
    //Scan every CNT in contact with the GNP
    for (int i = 0; i < (int)grouped_contacts.size(); i++) {
        
        //Flags that determine which point in the segment of consecutive contacts I am looking for
        int initial_flag = 1; //I will start looking for the initial point of a segment of consecutive contacts
        int final_flag = 0;
        
        //variables to store the initial and final contacts in a segment
        int initial = 0, final = 0;
        
        //Scan every contact in the CNT
        for (int j = 0; j < (int) grouped_contacts[i].size(); j++) {
            
            //If am looking fot the initial point check if the current contact is valid (!= -1)
            if (initial_flag && grouped_contacts[i][j] != -1) {
                
                //When an initial valid contact is found save the current contact as the initial contact of the segment
                initial = j;
                
                //Change the flag to zero as I will not be looking for an initial point
                initial_flag = 0;
                
                //Make sure the flag for the final point is 1
                final_flag = 1;
            }
            
            //If am looking fot the initial point check if the current contact is valid (!= -1)
            if (final_flag && grouped_contacts[i][j] == -1) {
                
                //When a final valid contact is found save the previous contact as the final contact of the segment
                final = j-1;
                
                //Change the flag to zero as I will not be looking for a final point
                final_flag = 0;
                
                //Make sure the flag for the initial point is 1
                initial_flag = 1;
                
                //A segment has been found, if the difference between final and initial is 0
                //thre is nothing to do as the contact will stay
                //but if the difference is 1 or more, then some contacts will be deleted
                int difference = final - initial;
                if ( difference > 0) {
                    
                    //find the middle contact
                    int middle = initial + (int)(difference/2);
                    
                    //Set the initial contact equal to the middle one
                    int initial_contact = grouped_contacts[i][initial];
                    int middle_contact = grouped_contacts[i][middle];
                    contacts[initial_contact] = contacts[middle_contact];
                    contacts[initial_contact] = contacts[middle_contact];
                    
                    
                    //Now, add all other contacts to the vector to_delete
                    for (int ii = initial+1; ii <= final; ii++) {
                        
                        to_delete.push_back(grouped_contacts[i][ii]);
                    }
                    
                }
                
            }
            
        }
        
    }
    
    return 1;
}
//This function checks if a contact has alraedy been created
//It is assumed that the CNT is the first particle
int Hoshen_Kopelman::Check_repeated_or_equivalent(const long int &point_cnt, const long int &point_gnp, const vector<Point_3D> &points_in, const vector<Point_3D> &points_gnp, vector<contact_pair> &contacts, vector<int> &gnp_contact_vector)
{
    //Get the CNT number
    int CNT = points_in[point_cnt].flag;
    
    //define the range of points to eliminate long tunnels
    int range = 10;
    
    //Iterate over the contacts of the GNP and check if it the contact has already been created
    for (int i = 0; i < (int)gnp_contact_vector.size(); i++) {
        
        //Current contact
        int cont_pair = gnp_contact_vector[i];
        
        //First check if the points are the same
        if (point_cnt == contacts[cont_pair].point1 && point_gnp == contacts[cont_pair].point2) {
            
            //The contact is repeated, thus terminate the function
            return 1;
        }
        //Check if the CNT is the same as in the current contact pair (CNTs are always particle1)
        //in case there is an equivalent contact
        else if (CNT == contacts[cont_pair].particle1) {
            
            //If the contact was not the same, now we check if the contact is equivalent
            
            //Compute the point numbers in the range to determine if the contact already exists
            long int start_range = contacts[cont_pair].point1-range, end_range = contacts[cont_pair].point1+range;
            
            //Check if point_cnt is inside the range
            if (start_range <= point_cnt && point_cnt <= end_range) {
                
                //If the CNT point is inside the range, then the contact is equivalent
                //Choose the pair of points whith the shortest squared distance (to reduce computations)
                double new_distance = points_in[point_cnt].squared_distance_to(points_gnp[point_gnp]);
                double previous_distance = points_in[contacts[cont_pair].point1].squared_distance_to(points_gnp[contacts[cont_pair].point2]);
                
                if ( new_distance < previous_distance) {
                    
                    //Update the points in contact
                    contacts[cont_pair].point1 = point_cnt;
                    contacts[cont_pair].point2 = point_gnp;
                }
                
                //The contact is equivalent, so return 1
                return 1;
            }
        }//*/
    }
    
    //An equivalent contact was not found, so return 0
    return 0;
}
//Renumber the two sets of labels to merge the clusters of CNTs and GNPs
int Hoshen_Kopelman::Merge_interparticle_labels(vector<int> &labels_mixed)
{
    //Iterate over the CNT labels
    for (int i = 0; i < (int)labels.size(); i++) {
        //Assign mixed label i to CNT i
        labels[i] = labels_mixed[i];
    }
    
    //Iterate over the GNP labels, which start after the last CNT label
    //i.e., at the mixed label index equal to the number of CNT labels
    int n_cnts = (int)labels.size();
    for (int i = n_cnts; i < (int)labels_mixed.size(); i++) {
        //Assign mixed label i to GNP i-n_cnts
        labels_gnp[i-n_cnts] = labels_mixed[i];
    }
    
    return 1;
}
//-------------------------------------------------------------------------------------------------------------------------------------