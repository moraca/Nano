//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Hoshen_Kopelman.h
//OBJECTIVE:	The Hoshen_Kopelman Algorithm
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef HOSHENKOPELMAN_H
#define HOSHENKOPELMAN_H

#include "Input_Reader.h"
#include "Printer.h"

//-------------------------------------------------------
class Hoshen_Kopelman
{
public:
    //Variables
    //Labels for HK76 (CNTs)
    vector<int> labels, labels_labels;
    //Labels for HK76 (GNPs)
    vector<int> labels_gnp, labels_labels_gnp;
    //Contact vectors
    vector<vector<long int> > contacts_point;
    vector<contact_pair> gnp_contacts;
    vector<contact_pair> mixed_contacts;
    //Cluster vectors to be used by other classes
    vector<vector<int> > clusters_cnt;
    vector<vector<int> > isolated;
    //Cluster vectors for hybrid particles
    vector<vector<int> > clusters_gch;
    vector<vector<int> > isolated_gch;
    //Data Member
    
    //Constructor
    Hoshen_Kopelman(){};
    
    //Member Functions
    //Hybrid particle
    int Determine_clusters(const struct Geom_RVE &sample, const struct Cutoff_dist &cutoffs, const vector<int> &cnts_inside, const vector<vector<long int> > &sectioned_domain, const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const vector<double> &radii, const vector<int> &gnps_inside, const vector<vector<long int> > &sectioned_domain_gnp, const vector<vector<long int> > &structure_gnp, const vector<Point_3D> &points_gnp, const vector<GCH> &hybrid_particles);
    int Scan_sub_regions_cnt(const struct Geom_RVE &sample, const vector<Point_3D> &points_in, const vector<int> &gnps_inside, const vector<GCH> &hybrid_particles, const vector<double> &radii, const double &tunnel, const vector<vector<long int> > &sectioned_domain, const vector<vector<long int> > &structure);
    int Group_cnts_in_gnp(const vector<GCH> &hybrid_particles, const vector<int> &gnps_inside, int &new_label);
    int Check_repeated(const vector<long int> &contacts_vector, const long int &point);
    int HK76(const int &CNT1, const int &CNT2, int &new_label, vector<int> &labels, vector<int> &labels_labels);
    int Find_root(const int &L, vector<int> &labels_labels);
    int Merge_labels(const int &root1, const int &root2, vector<int> &labels_labels);
    int Scan_sub_regions_gnp(const vector<Point_3D> &points_gnp, const vector<GCH> &hybrid_particles, const double &tunnel, const vector<vector<long int> > &sectioned_domain_gnp);
    int Initialize_contact_matrices(const int &n_GNPs,vector<vector<long int> > &point_matrix, vector<vector<double> > &distance_matrix);
    int Create_vector_of_gnp_contacts(const vector<vector<long int> > &point_matrix);
    int Make_particle_clusters(const int &n_clusters, const vector<int> &particles_inside, vector<int> &labels, vector<vector<int> > &isolated, vector<vector<int> > &clusters_particles);
    int Cleanup_labels(vector<int> &labels, vector<int> &labels_labels);
    int Scan_sub_regions_cnt_and_gnp(const struct Geom_RVE &sample, const double &tunnel, const vector<Point_3D> &points_in, const vector<double> &radii, const vector<vector<long int> > &sectioned_domain, const vector<Point_3D> &points_gnp, const vector<GCH> &hybrid_particles, const vector<int> &gnps_inside, const vector<vector<long int> > &sectioned_domain_gnp, vector<int> &labels_mixed, vector<int> &labels_labels_mixed);
    int Initialize_mixed_labels(vector<int> &labels_mixed, vector<int> &labels_labels_mixed);
    int Fill_cnt_gnp_numbers(const vector<GCH> &hybrid_particles, vector<int> &cnt_gnp_numbers);
    int Cluster_gnps_and_cnts(const vector<GCH> &hybrid_particles, const vector<int> &gnps_inside, vector<int> &labels_mixed, vector<int> &labels_labels_mixed, int &new_label);
    int Check_repeated_or_equivalent(const long int &point_cnt, const long int &point_gnp, const double &distance, const vector<Point_3D> &points_in, const vector<Point_3D> &points_gnp, vector<contact_pair> &contacts, vector<int> &gnp_contact_vector);
    void Update_contact(const long int &point1, const int &particle1, const long int &point2, const int &particle2, struct contact_pair &contact);
    int Merge_interparticle_labels(vector<int> &labels_mixed);

private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================