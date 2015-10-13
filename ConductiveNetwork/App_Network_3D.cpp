//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	App_Network_3D.cpp
//OBJECTIVE:	Create conductive nanotube network in 3D separated by backbone paths, dead branches and isolated clusters
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "App_Network_3D.h"

//Generate 3D conductive nanotube network separated by backbone paths, dead branches and isolated clusters
int App_Network_3D::Create_conductive_network_3D(Input *Init)const
{
	//Time markers for total simulation
	time_t ct0, ct1;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Network Generation with overlapping
	ct0 = time(NULL);
	
	vector<vector<Point_3D> > cnts_points;	//define two-dimensional vector of three-dimensional points for storing the CNT network
    vector<double> cnts_radius;						//define the radius of each nanotube in the network
    vector<Point_3D> points;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Network Generation with overlapping
	hout << "-_- To generate networks with overlapping......"<<endl;
	GenNetwork *Genet = new GenNetwork;
	ct0 = time(NULL);
	if(Genet->Generate_geometric_networks(Init->geom_rve, Init->cluster_geo, Init->nanotube_geo, cnts_points, cnts_radius)==0) return 0;
	ct1 = time(NULL);
	hout << "Network generation time: "<<(int)(ct1-ct0)<<" secs."<<endl;
    
    //Checking the angle between two segments in one nanotube (if less than PI/2, provide an alarm)
    if(Genet->CNTs_quality_testing(cnts_points)==0) return 0;
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //A new class of Tecplot_Export
    Tecplot_Export *Tecexpt = new Tecplot_Export;
    
    struct cuboid cub;														//Generate a cuboid for RVE
    cub.poi_min = Init->geom_rve.ex_origin;
    cub.len_x = Init->geom_rve.ex_len;
    cub.wid_y = Init->geom_rve.ey_wid;
    cub.hei_z = Init->geom_rve.ez_hei;
    
    //The geometric structure of CNT network (by threads in Tecplot)
    ct0 = time(NULL);
    if(Tecexpt->Export_network_threads(cub, cnts_points)==0) return 0;
    //The geometric structure of CNT network (by tetrahedron meshes in Tecplot) //Attention: little parts of nanotube volumes out of the cuboid
    if(Tecexpt->Export_cnt_network_meshes(cub, cnts_points, cnts_radius)==0) return 0;
    ct1 = time(NULL);
    hout << "Export_network_threads time: "<<(int)(ct1-ct0)<<" secs."<<endl;

    
    //-----------------------------------------------------------------------------------------------------------------------------------------
	ct0 = time(NULL);
    Background_grid *Bckg = new Background_grid;
    //From this function I get the intenal variables sectioned_domain_cnt and structure
    if (Bckg->Generate_background_grids(Init->geom_rve, Init->nanotube_geo, points) == 0) return 0;
	ct1 = time(NULL);
	hout << "Generate_background_grids time: "<<(int)(ct1-ct0)<<" secs."<<endl;
    
	for(int i=0; i<=Init->geom_rve.cut_num; i++)
	{
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine the local networks in cutoff windons
        Cutoff_Wins *Cutwins = new Cutoff_Wins;
        //From this function I get the internal variables cnts_inside and boundary_cnt
        if(Cutwins->Extract_observation_window(Bckg->sectioned_domain_cnt, i, Init->geom_rve, Init->nanotube_geo, Bckg->structure, cnts_radius, points)==0) return 0;
        
        //Calculate the dimensions of the cunrrent observation window
        //In this way the observation window goes from a large to a small one
        double lx = Init->geom_rve.len_x - i*Init->geom_rve.win_delt_x;
        double ly = Init->geom_rve.wid_y - i*Init->geom_rve.win_delt_y;
        double lz = Init->geom_rve.hei_z - i*Init->geom_rve.win_delt_z;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine the local networks in cutoff windons
        Contact_grid *Contacts = new Contact_grid;
        if (Contacts->Generate_contact_grid(Bckg->structure, Cutwins->cnts_inside, points, Init->geom_rve, Init->cutoff_dist, Init->nanotube_geo, lx, ly, lz)==0) return 0;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
		//Hoshen-Kopelman algorithm
		Hoshen_Kopelman *HoKo = new Hoshen_Kopelman;
		if(HoKo->Determine_nanotube_clusters(Bckg->structure, Contacts->sectioned_domain, points, Cutwins->cnts_inside, cnts_radius, Init->cutoff_dist)==0) return 0;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine percolation
        Percolation *Perc = new Percolation;
        if (Perc->Determine_percolating_clusters(Cutwins->boundary_cnt, HoKo->labels, HoKo->labels_labels, HoKo->label_map, HoKo->clusters_cnt, HoKo->isolated, Init->geom_rve, Init->nanotube_geo)==0) return 0;
        
        //These vectors are used to store the fractions of the different families in the current observation window
        //families_lengths has 7 elements because of the 6 percolated families and the non-percoalted CNTs, the same is true for fractions
        vector<double> families_lengths(7,0);
        vector<double> fractions(7,0);
        vector<double> branches_lengths(6,0);
        
        //Loop over the different clusters so that the direct electrifying algorithm is apllied on each cluster
        if (HoKo->clusters_cnt.size()) {
            for (int i = 0; i < (int)HoKo->clusters_cnt.size(); i++) {
                //-----------------------------------------------------------------------------------------------------------------------------------------
                //Direct Electrifying algorithm
                Direct_Electrifying *DEA = new Direct_Electrifying;
                if(DEA->Calculate_voltage_field(Bckg->structure, HoKo->contacts_point, Cutwins->boundary_flags, HoKo->clusters_cnt[i], cnts_radius, Perc->family[i], Init->electric_para)==0) return 0;
                
                //-----------------------------------------------------------------------------------------------------------------------------------------
                //Determine the backbone and dead branckes
                Backbone_Network *Backbonet = new Backbone_Network;
                if(Backbonet->Determine_backbone_network(Perc->family[i], HoKo->clusters_cnt[i], DEA->voltages, DEA->LM_matrix, DEA->elements, Bckg->structure, points, families_lengths, branches_lengths)==0) return 0;

            }
        } else {
            hout << "There are no percolating clusters" << endl;
        }
        
        //Calculate the fractions of CNTs that belong to each family and save them to a file
        Clusters_fractions *Fracs = new Clusters_fractions;
        if (Fracs->Calculate_fractions(Bckg->structure, Cutwins->cnts_inside, HoKo->isolated, points, families_lengths, branches_lengths, fractions)==0) return 0;
	}

	return 1;
}
//===========================================================================
