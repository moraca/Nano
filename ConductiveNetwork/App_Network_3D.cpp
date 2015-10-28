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
    Background_vectors *Bckg = new Background_vectors;
    Geom_RVE geo = Init->geom_rve;
    Nanotube_Geo nano = Init->nanotube_geo;
    if (Bckg->Generate_shells_and_structure(Init->geom_rve, Init->nanotube_geo, cnts_point, shells_cnt)==0) return 0;
	ct1 = time(NULL);
	hout << "Generate shells and structure time: "<<(int)(ct1-ct0)<<" secs."<<endl;//*/
    
	for(int i=0; i<=Init->geom_rve.cut_num; i++)
	{
        hout << "======================================================"<<endl;
        hout << "Iteration " << i << endl;
        time_t it0, it1;
        it0 = time(NULL);
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine the local networks in cutoff windons
        Cutoff_Wins *Cutwins = new Cutoff_Wins;
        //From this function I get the internal variables cnts_inside and boundary_cnt
        ct0 = time(NULL);
        if(Cutwins->Extract_observation_window(Init->geom_rve, Init->nanotube_geo, cnts_structure, cnts_radius, cnts_point, shells_cnt, i)==0) return 0;
        ct1 = time(NULL);
        hout << "Extract observation window time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine the local networks inside the cutoff windows
        Contact_grid *Contacts = new Contact_grid;
        ct0 = time(NULL);
        if (Contacts->Generate_contact_grid(Init->geom_rve, Init->cutoff_dist, Init->nanotube_geo, Cutwins->cnts_inside, cnts_structure, cnts_point, i)==0) return 0;
        ct1 = time(NULL);
        hout << "Generate contact grid time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
		//Hoshen-Kopelman algorithm
		Hoshen_Kopelman *HoKo = new Hoshen_Kopelman;
        ct0 = time(NULL);
		if(HoKo->Determine_nanotube_clusters(Init->cutoff_dist, Cutwins->cnts_inside, Contacts->sectioned_domain, cnts_structure, cnts_point, cnts_radius)==0) return 0;
        ct1 = time(NULL);
        hout << "Determine nanotube clusters time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine percolation
        Percolation *Perc = new Percolation;
        ct0 = time(NULL);
        if (Perc->Determine_percolating_clusters(Init->geom_rve, Init->nanotube_geo, Cutwins->boundary_cnt, HoKo->labels, HoKo->labels_labels, HoKo->label_map, HoKo->clusters_cnt, HoKo->isolated)==0) return 0;
        ct1 = time(NULL);
        hout << "Determine percolating clusters time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //These vectors are used to store the fractions of the different families in the current observation window
        //families_lengths has 8 elements because of the 7 percolated families and the non-percoalted CNTs, the same is true for fractions
        vector<double> families_lengths(8,0);
        vector<double> fractions(8,0);
        vector<double> branches_lengths(7,0);

        //Loop over the different clusters so that the direct electrifying algorithm is apllied on each cluster
        if (HoKo->clusters_cnt.size()) {
            //hout << "clusters_cnt.size()="<<HoKo->clusters_cnt.size()<<endl;
            for (int i = 0; i < (int)HoKo->clusters_cnt.size(); i++) {
                //-----------------------------------------------------------------------------------------------------------------------------------------
                //Direct Electrifying algorithm
                Direct_Electrifying *DEA = new Direct_Electrifying;
                ct0 = time(NULL);
                if(DEA->Calculate_voltage_field(cnts_structure, HoKo->contacts_point, Cutwins->boundary_flags, HoKo->clusters_cnt[i], cnts_radius, Perc->family[i], Init->electric_para)==0) return 0;
                ct1 = time(NULL);
                hout << "Calculate voltage field time: "<<(int)(ct1-ct0)<<" secs."<<endl;
                
                //-----------------------------------------------------------------------------------------------------------------------------------------
                //Determine the backbone and dead branckes
                Backbone_Network *Backbonet = new Backbone_Network;
                ct0 = time(NULL);
                if(Backbonet->Determine_backbone_network(Perc->family[i], HoKo->clusters_cnt[i], DEA->voltages, DEA->LM_matrix, DEA->elements,cnts_structure, cnts_point, families_lengths, branches_lengths)==0) return 0;
                ct1 = time(NULL);
                hout << "Determine backbone network time: "<<(int)(ct1-ct0)<<" secs."<<endl;

            }
        } else {
            hout << "There are no percolating clusters" << endl;
        }
        
        //Calculate the fractions of CNTs that belong to each family and save them to a file
        Clusters_fractions *Fracs = new Clusters_fractions;
        ct0 = time(NULL);
        if (Fracs->Calculate_fractions(cnts_structure, Cutwins->cnts_inside, HoKo->isolated, cnts_point, families_lengths, branches_lengths, fractions)==0) return 0;
        ct1 = time(NULL);
        hout << "Calculate fractions time: "<<(int)(ct1-ct0)<<" secs."<<endl;

        it1 = time(NULL);
        hout << "Iteration "<<i<<" time: "<<(int)(it1-it0)<<" secs."<<endl;
    }

	return 1;
}
//===========================================================================
