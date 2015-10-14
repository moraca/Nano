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
	
    vector<double> cnts_radius;						//Define the radius of each nanotube in the network
    vector<Point_3D> cnts_point;					//Define the set of cnt point in a 1D vector
    vector<vector<long int> > cnts_structure;		//The global number of points in the cnts
    vector<vector<int> > shells_cnt;                //Shell sub-regions to make the triming faster
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Network Generation with overlapping
    hout << "-_- To generate nanotube network......" << endl;
    ct0 = time(NULL);
    GenNetwork *Genet = new GenNetwork;
    if(Genet->Generate_nanotube_networks(Init->geom_rve, Init->cluster_geo, Init->nanotube_geo, cnts_point, cnts_radius, cnts_structure)==0) return 0;
    delete Genet;
    ct1 = time(NULL);
    hout << "Nanotube network generation time: " << (int)(ct1-ct0) <<" secs." << endl;
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
	ct0 = time(NULL);
    //
    Background_vectors *Bckg = new Background_vectors;
    Geom_RVE geo = Init->geom_rve;
    Nanotube_Geo nano = Init->nanotube_geo;
    if (Bckg->Generate_shells_and_structure(Init->geom_rve, Init->nanotube_geo, cnts_point, shells_cnt)==0) return 0;
	ct1 = time(NULL);
	hout << "Generate shells and structure time: "<<(int)(ct1-ct0)<<" secs."<<endl;//*/
    
	for(int i=0; i<=Init->geom_rve.cut_num; i++)
	{
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine the local networks in cutoff windons
        Cutoff_Wins *Cutwins = new Cutoff_Wins;
        //From this function I get the internal variables cnts_inside and boundary_cnt
        if(Cutwins->Extract_observation_window(shells_cnt, i, Init->geom_rve, Init->nanotube_geo, cnts_structure, cnts_radius, cnts_point)==0) return 0;
        
        //Calculate the dimensions of the cunrrent observation window
        //In this way the observation window goes from a large to a small one
        double lx = Init->geom_rve.len_x - i*Init->geom_rve.win_delt_x;
        double ly = Init->geom_rve.wid_y - i*Init->geom_rve.win_delt_y;
        double lz = Init->geom_rve.hei_z - i*Init->geom_rve.win_delt_z;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine the local networks in cutoff windons
        Contact_grid *Contacts = new Contact_grid;
        if (Contacts->Generate_contact_grid(cnts_structure, Cutwins->cnts_inside, cnts_point, Init->geom_rve, Init->cutoff_dist, Init->nanotube_geo, lx, ly, lz)==0) return 0;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
		//Hoshen-Kopelman algorithm
		Hoshen_Kopelman *HoKo = new Hoshen_Kopelman;
		if(HoKo->Determine_nanotube_clusters(cnts_structure, Contacts->sectioned_domain, cnts_point, Cutwins->cnts_inside, cnts_radius, Init->cutoff_dist)==0) return 0;
        
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
                if(DEA->Calculate_voltage_field(cnts_structure, HoKo->contacts_point, Cutwins->boundary_flags, HoKo->clusters_cnt[i], cnts_radius, Perc->family[i], Init->electric_para)==0) return 0;
                
                //-----------------------------------------------------------------------------------------------------------------------------------------
                //Determine the backbone and dead branckes
                Backbone_Network *Backbonet = new Backbone_Network;
                if(Backbonet->Determine_backbone_network(Perc->family[i], HoKo->clusters_cnt[i], DEA->voltages, DEA->LM_matrix, DEA->elements,cnts_structure, cnts_point, families_lengths, branches_lengths)==0) return 0;

            }
        } else {
            hout << "There are no percolating clusters" << endl;
        }
        
        //Calculate the fractions of CNTs that belong to each family and save them to a file
        Clusters_fractions *Fracs = new Clusters_fractions;
        if (Fracs->Calculate_fractions(cnts_structure, Cutwins->cnts_inside, HoKo->isolated, cnts_point, families_lengths, branches_lengths, fractions)==0) return 0;
	}//*/

	return 1;
}
//===========================================================================
