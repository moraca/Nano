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
    //Data Member
    
    //Constructor
    Hoshen_Kopelman(){};
    
    //Member Functions
    int Determine_nanotube_clusters(vector<vector<long int> > structure, vector<vector<long int> > sectioned_domain, vector<Point_3D> points_in, vector<int> cnts_inside, vector<double> radii, double tunnel, vector<vector<long int> > &contacts_point, vector<vector<long int> > &contacts_cnt_point, vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated);
    int Contacts_and_HK76(vector<Point_3D> points_in, vector<double> radii, double tunnel, vector<vector<long int> > &contacts_point, vector<vector<long int> > &contacts_cnt_point, vector<vector<long int> > sectioned_domain);
    void Fill_contact_vectors(long int P1, long int P2, int CNT1, int CNT2, vector<vector<long int> > &contacts_point, vector<vector<long int> > &contacts_cnt_point);
    int Check_repeated(vector<long int> region, long int Point);
    int HK76(int CNT1, int CNT2);
    int Find_root(int L);
    int Merge_labels(int root1, int root2);
    int Make_CNT_clusters(vector<vector<long int> > structure, vector<Point_3D> points_in, vector<int> cnts_inside, vector<vector<int> > &clusters_cnt, vector<vector<int> > &isolated);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================