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

//-------------------------------------------------------
class Hoshen_Kopelman
{
public:
    //Variables
    //Labels for HK76
    vector<int> labels, labels_labels, label_map;
    //VEctors to be used by other classes
    vector<vector<long int> > contacts_point;
    vector<vector<int> > clusters_cnt;
    vector<vector<int> > isolated;
    //Data Member
    
    //Constructor
    Hoshen_Kopelman(){};
    
    //Member Functions
    int Determine_nanotube_clusters(vector<vector<long int> > structure, vector<vector<long int> > sectioned_domain, vector<Point_3D> points_in, vector<int> cnts_inside, vector<double> radii, struct Cutoff_dist cutoffs);
    int Scan_sub_regions(vector<Point_3D> points_in, vector<double> radii, double tunnel, vector<vector<long int> > sectioned_domain);
    int HK76(int CNT1, int CNT2);
    int Find_root(int L);
    int Merge_labels(int root1, int root2);
    void Delete_repeated_contacts();
    void Discard_repeated(vector<long int> &vec);
    int Make_CNT_clusters(vector<vector<long int> > structure, vector<Point_3D> points_in, vector<int> cnts_inside);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================