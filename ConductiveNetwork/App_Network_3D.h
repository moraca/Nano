//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	App_Network_3D.h
//OBJECTIVE:	Create a 3D nanotube netwok
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef APPNETWORK3D_H
#define APPNETWORK3D_H

#include "time.h"
#include "Hns.h"
using namespace hns;

#include "Input_Reader.h"
#include "GenNetwork.h"
#include "Cutoff_Wins.h"
#include "Hoshen_Kopelman.h"
#include "Direct_Electrifying.h"
#include "Backbone_Network.h"
#include "Background_vectors.h"
#include "Contact_grid.h"
#include "Percolation.h"
#include "Clusters_fractions.h"
#include "Tecplot_Export.h"

//---------------------------------------------------------------------------
class App_Network_3D
{
public:
    //Data Member
    //vector<Point_3D> cnps;			//Define 3D point verctor of nanotuber points
    //vector<double> cnts_radius;		//Define the radius of every nanotube in the network
    
    //Constructor
    App_Network_3D(){};
    
    //Member Functions
    int Create_conductive_network_3D(Input *Init)const;
    int Export_tecplot_files(const int &iter, const struct Geom_RVE &sample, const vector<Point_3D> &points_in, const vector<double> &radii, const vector<vector<long int> > &structure, const vector<vector<int> > &isolated, vector<vector<long int> > &all_dead_indices, const vector<vector<long int> > &all_indices)const;
    int Export_tecplot_files_for_clusters(const string &type, const int &iter, const struct Geom_RVE &sample, const vector<Point_3D> &points_in, const vector<double> &radii, const vector<vector<long int> > &structure, const vector<vector<int> > &clusters_cnt, const vector<vector<int> > &isolated)const;
    int Convert_index_to_structure(const vector<long int> &indices, vector<vector<long int> > &structure)const;
    int Convert_cluster_to_structure(const vector<int> &cluster, const vector<vector<long int> > &structure_in, vector<vector<long int> > &structure_out)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================

