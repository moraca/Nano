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
    vector<Point_3D> gnps_point;          //Define the set of GNP point in a 1D vector
    vector<vector<long int> > cnts_structure;		//The global number of points in the cnts
    vector<vector<long int> > gnps_structure;		//The global number of points in the gnps
    vector<vector<int> > shells_cnt;                //Shell sub-regions to make the triming faster
    vector<GCH> hybrid_particles;                   //Define the set of GNP hybrid particles
    vector<vector<int> > shells_gnps;                //Shell sub-regions to make the deletion of GNPs faster
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Network Generation with overlapping
    hout << "-_- To generate nanotube network......" << endl;
    ct0 = time(NULL);
    GenNetwork *Genet = new GenNetwork;
    if(Genet->Generate_nanofiller_network(Init->geom_rve, Init->cluster_geo, Init->nanotube_geo, Init->gnp_geo, Init->cutoff_dist, cnts_point, gnps_point, hybrid_particles, cnts_radius, cnts_structure, gnps_structure)==0) return 0;
    delete Genet;
    ct1 = time(NULL);
    hout << "Nanotube network generation time: " << (int)(ct1-ct0) <<" secs." << endl;
    //Printer *P = new Printer;
    //P->Print_1d_vec(gnps_point, "gnps_point_00.txt");
    //delete P;
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
	ct0 = time(NULL);
    Background_vectors *Bckg = new Background_vectors;
    Geom_RVE geo = Init->geom_rve;
    Nanotube_Geo nano = Init->nanotube_geo;
    if (Bckg->Generate_shells(Init->geom_rve, Init->nanotube_geo, cnts_point, hybrid_particles, shells_cnt, shells_gnps)==0) return 0;
	ct1 = time(NULL);
	hout << "Generate shells and structure time: "<<(int)(ct1-ct0)<<" secs."<<endl;//*/
    
	for(int i=0; i<=Init->geom_rve.cut_num; i++)
	{
        hout << "======================================================"<<endl;
        hout << "Iteration " << i << endl;
        time_t it0, it1;
        it0 = time(NULL);
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //These vectors are used to export tecplot files
        vector<long int> empty;
        vector<vector<long int> > all_dead_indices(7,empty);
        vector<vector<long int> > all_indices(7,empty);
        vector<int> empty_int;
        vector<vector<int> > gnp_indices(7,empty_int);
        vector<vector<int> > gnp_dead_indices(7,empty_int);
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine the local networks in cutoff windons
        Cutoff_Wins *Cutwins = new Cutoff_Wins;
        //From this function I get the internal variables cnts_inside and boundary_cnt
        ct0 = time(NULL);
        if(Cutwins->Extract_observation_window(i, Init->geom_rve, Init->nanotube_geo, Init->gnp_geo, hybrid_particles, cnts_structure, gnps_structure, cnts_radius, cnts_point, gnps_point, shells_cnt, shells_gnps)==0) return 0;
        ct1 = time(NULL);
        hout << "Extract observation window time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine the local networks inside the cutoff windows
        Contact_grid *Contacts = new Contact_grid;
        ct0 = time(NULL);
        if (Contacts->Generate_contact_grid(i, Init->geom_rve, Init->cutoff_dist, Init->nanotube_geo, Cutwins->cnts_inside, cnts_point, cnts_structure, Cutwins->gnps_inside, gnps_point, gnps_structure)==0) return 0;
        ct1 = time(NULL);
        hout << "Generate contact grid time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
		//Hoshen-Kopelman algorithm
		Hoshen_Kopelman *HoKo = new Hoshen_Kopelman;
        ct0 = time(NULL);
        if (HoKo->Determine_clusters(Init->geom_rve, Init->cutoff_dist, Cutwins->cnts_inside, Contacts->sectioned_domain, cnts_structure, cnts_point, cnts_radius, Cutwins->gnps_inside, Contacts->sectioned_domain_gnps, gnps_structure, gnps_point, hybrid_particles)==0) return 0;
        ct1 = time(NULL);
        hout << "Determine nanotube clusters time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        /*/Save cluters into a file
        ct0 = time(NULL);
        if (Export_tecplot_files_for_clusters("Cluster", i, Init->geom_rve, cnts_point, cnts_radius, cnts_structure, HoKo->clusters_cnt, HoKo->isolated, hybrid_particles, HoKo->clusters_gch, HoKo->isolated_gch)==0) return 0;
        ct1 = time(NULL);
        hout << "Export tecplot files time: "<<(int)(ct1-ct0)<<" secs."<<endl;//*/
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine percolation
        Percolation *Perc = new Percolation;
        ct0 = time(NULL);
        if (Perc->Determine_percolated_clusters(i, Init->geom_rve, Init->nanotube_geo, Init->gnp_geo, Cutwins->boundary_cnt, HoKo->labels, Cutwins->boundary_gnp, HoKo->labels_gnp, HoKo->clusters_cnt, HoKo->isolated, HoKo->clusters_gch, HoKo->isolated_gch) == 0) return 0;
        ct1 = time(NULL);
        hout << "Determine percolating clusters time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        /*/-----------------------------------------------------------------------------------------------------------------------------------------
        //Save cluters into a file
        ct0 = time(NULL);
        if (Export_tecplot_files_for_clusters("Percolated", i, Init->geom_rve, cnts_point, cnts_radius, cnts_structure, HoKo->clusters_cnt, HoKo->isolated, hybrid_particles, HoKo->clusters_gch, HoKo->isolated_gch)==0) return 0;
        ct1 = time(NULL);
        hout << "Export tecplot files time: "<<(int)(ct1-ct0)<<" secs."<<endl;//*/
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //These vectors are used to store the fractions of the different families in the current observation window
        //families_lengths has 8 elements because of the 7 percolated families and the non-percoalted CNTs, the same is true for fractions
        vector<double> families_lengths(8,0);
        vector<double> fractions(8,0);
        vector<double> branches_lengths(7,0);

        //Loop over the different clusters so that the direct electrifying algorithm is aplied on each cluster
        if (HoKo->clusters_cnt.size() || HoKo->clusters_gch.size()) {
            //Perform the electrical analysis to obtain the backbone and calculate the electrical resistance
            ct0 = time(NULL);
            Electrical_analysis *Electric_A = new Electrical_analysis;
            if (Electric_A->Perform_analysis_on_clusters( i, Perc->family, HoKo, Cutwins, cnts_structure, cnts_point, cnts_radius, gnps_structure, gnps_point, Init->geom_rve, Init->electric_para, Init->cutoff_dist, hybrid_particles, families_lengths, branches_lengths, all_dead_indices, all_indices, gnp_dead_indices, gnp_indices)==0) return 0;
            ct1 = time(NULL);
            hout << "Perform electrical analysis time: "<<(int)(ct1-ct0)<<" secs."<<endl;
            delete Electric_A;
            
        } else {
            hout << "There are no percolated clusters" << endl;
        }
        
        /*/Calculate the fractions of CNTs that belong to each family and save them to a file
        Clusters_fractions *Fracs = new Clusters_fractions;
        ct0 = time(NULL);
        if (Fracs->Calculate_fractions(cnts_structure, Cutwins->cnts_inside, HoKo->isolated, cnts_point, families_lengths, branches_lengths, fractions)==0) return 0;
        ct1 = time(NULL);
        hout << "Calculate fractions time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        delete Fracs;//*/
        
        /*/Export 3D tecplot files
        ct0 = time(NULL);
        //if (Export_tecplot_files(i, Init->geom_rve, cnts_point, cnts_radius, hybrid_particles, HoKo->isolated_gch, cnts_structure, HoKo->isolated, all_dead_indices, all_indices, gnp_dead_indices, gnp_indices)==0) return 0;
        ct1 = time(NULL);
        hout << "Export tecplot files time: "<<(int)(ct1-ct0)<<" secs."<<endl;//*/
        
        /*/Export 2D tecplot files for triangulations
        ct0 = time(NULL);
        if (Export_triangulation_tecplot_files(i, Init->geom_rve, cnts_point, hybrid_particles)==0) return 0;
        ct1 = time(NULL);
        hout << "Export triangulation tecplot files time: "<<(int)(ct1-ct0)<<" secs."<<endl;//*/
        
        it1 = time(NULL);
        hout << "Iteration "<<i<<" time: "<<(int)(it1-it0)<<" secs."<<endl;
        
        //Delete objects to free memory
        delete Cutwins;
        delete Contacts;
        delete HoKo;
        delete Perc;
        
    }
    

	return 1;
}
//Export tecplot files
int App_Network_3D::Export_tecplot_files(const int &iter, const struct Geom_RVE &sample, const vector<Point_3D> &points_in, const vector<double> &radii, const vector<GCH> &hybrid_particles, const vector<int> &isolated_gch, const vector<vector<long int> > &structure, const vector<vector<int> > &isolated, vector<vector<long int> > &all_dead_indices, const vector<vector<long int> > &all_indices, vector<vector<int> > &gnp_dead_indices, const vector<vector<int> > &gnp_indices)const
{
    //These vectors will be used to create a structure-type vector to export to tecplot files
    vector<vector<long int> > percolated_tmp, dead_tmp;
    
    //These vectors will be used to create vectors of hybrid particles for percoalted and non-percolated clusters
    vector<GCH> percolated_particles, isolated_particles;
    
    //Arrange the clusters of hybrid particles according the the family of the cluster they belong to
    //This will make the data structures match
    
    //Tecplot export object
    Tecplot_Export *tec360 = new Tecplot_Export;
    
    //Geometry of observation window saved into a cuboid
    struct cuboid cub;
    //Dimensions of the current observation window
    cub.len_x = sample.win_max_x - iter*sample.win_delt_x;
    cub.wid_y = sample.win_max_y - iter*sample.win_delt_y;
    cub.hei_z = sample.win_max_z - iter*sample.win_delt_z;
    //These variables are the coordinates of the lower corner of the observation window
    cub.poi_min.x = sample.origin.x + (sample.len_x - cub.len_x)/2;
    cub.poi_min.y = sample.origin.y + (sample.wid_y - cub.wid_y)/2;
    cub.poi_min.z = sample.origin.z + (sample.hei_z - cub.hei_z)/2;

    //Names to be used to save visualization files
    vector<string> filenames;
    Initialize_filenames(filenames);
    
    /*/Export pecolated clusters and their dead branches
    for (int i = 0; i < 7; i++){
        //Check if the family is non empty. If it is non empty then a visualization file can be created
        if (all_indices[i].size()) {
            //Generate structure-type vectors
            if (!Convert_index_to_structure(all_indices[i], percolated_tmp)) {
                hout << "Error in Export_tecplot_files when converting percolated indices to structure"<<endl;
                return 0;
            }
            //Generate Tecplot files
            if ( !(tec360->Export_cnt_network_meshes(cub, points_in, radii, hybrid_particles, gnp_indices[i], percolated_tmp, filenames[i])))  {
                hout << "Error in Export_tecplot_files while translating and exporting directional clusters" <<endl;
                return 0;
            }
            //Clear the variable for the next iteration
            percolated_tmp.clear();
        }
        
        //It is possible that a CNT spans from one boundary to the other and has no contacts. In this case
        //The CNT is a cluster itself and has no dead branches, so I need to check separately if all_dead_indices[i]
        //is non empty. That is, the fact that a cluster percolates does not mean that there will be dead branches
        if (all_dead_indices[i].size()) {
            //Generate structure-type vectors
            if (!Convert_index_to_structure(all_dead_indices[i], dead_tmp)) {
                hout << "Error in Export_tecplot_files when converting dead indices to structure"<<endl;
                return 0;
            }
            //Generate Tecplot files
            if ( !(tec360->Export_cnt_network_meshes(cub, points_in, radii, hybrid_particles, gnp_dead_indices[i], dead_tmp, filenames[i+7])) ) {
                hout << "Error in Export_tecplot_files while translating and exporting directional clusters" <<endl;
                return 0;
            }
            //Clear the variable for the next iteration
            dead_tmp.clear();
        }
        
    }//*/
    
    //Create a structure vector for isolated CNTs
    vector<vector<long int> > iso_structure;
    for (int i = 0; i < (int)isolated.size(); i++) {
        for (int j = 0 ; j < (int)isolated[i].size(); j++) {
            int CNT = isolated[i][j];
            iso_structure.push_back(structure[CNT]);
        }
    }
    /*/Export the isolated CNTs
    if ( !(tec360->Export_cnt_network_meshes(cub, points_in, radii, hybrid_particles, isolated_gch, iso_structure, "SingleZone_isolated.dat")) ) {
        hout << "Error in Export_tecplot_files while translating and exporting directional clusters" <<endl;
        return 0;
    }//*/
    
    //Variables to use the command line
    int s;
    char command[100];
    //Move the visualization files to a new folder
    s = sprintf(command, "mkdir iter_%.4d", iter);
    s = system(command);
    s = sprintf(command, "mv Single*.dat iter_%.4d", iter);
    s = system(command);
    
    //delete tecplot object
    delete tec360;
    
    return 1;
    
}
//Filenames for the clusters
void App_Network_3D::Initialize_filenames(vector<string> &filenames)const
{
    //Save a separate file for each percolating direction
    filenames.push_back("SingleZone_00_xx.dat");
    filenames.push_back("SingleZone_01_yy.dat");
    filenames.push_back("SingleZone_02_zz.dat");
    filenames.push_back("SingleZone_03_xx_yy.dat");
    filenames.push_back("SingleZone_04_xx_zz.dat");
    filenames.push_back("SingleZone_05_yy_zz.dat");
    filenames.push_back("SingleZone_06_xx_yy_zz.dat");
    //Save a separate file for the dead branches of each percolating direction
    filenames.push_back("SingleZoneDead_00_xx.dat");
    filenames.push_back("SingleZoneDead_01_yy.dat");
    filenames.push_back("SingleZoneDead_02_zz.dat");
    filenames.push_back("SingleZoneDead_03_xx_yy.dat");
    filenames.push_back("SingleZoneDead_04_xx_zz.dat");
    filenames.push_back("SingleZoneDead_05_yy_zz.dat");
    filenames.push_back("SingleZoneDead_06_xx_yy_zz.dat");
}
//Export tecplot files
int App_Network_3D::Export_tecplot_files_for_clusters(const string &type, const int &iter, const struct Geom_RVE &sample, const vector<Point_3D> &points_in, const vector<double> &radii, const vector<vector<long int> > &structure, const vector<vector<int> > &clusters_cnt, const vector<vector<int> > &isolated, const vector<GCH> &hybrid_particles, const vector<vector<int> > &clusters_gch, const vector<vector<int> > &isolated_gch)const
{
    
    //Tecplot export object
    Tecplot_Export *tec360 = new Tecplot_Export;
    
    //Geometry of observation window saved into a cuboid
    struct cuboid cub;
    //Dimensions of the current observation window
    cub.len_x = sample.win_max_x - iter*sample.win_delt_x;
    cub.wid_y = sample.win_max_y - iter*sample.win_delt_y;
    cub.hei_z = sample.win_max_z - iter*sample.win_delt_z;
    //These variables are the coordinates of the lower corner of the observation window
    cub.poi_min.x = sample.origin.x + (sample.len_x - cub.len_x)/2;
    cub.poi_min.y = sample.origin.y + (sample.wid_y - cub.wid_y)/2;
    cub.poi_min.z = sample.origin.z + (sample.hei_z - cub.hei_z)/2;
    
    //Get the number of clusters
    int n_clusters = 0;
    if (clusters_cnt.size()) {
        n_clusters = (int)clusters_cnt.size();
        
    } else if (clusters_gch.size()) {
        n_clusters = (int)clusters_gch.size();
    }
    
    //Loop over the clusters
    for (int i = 0; i < n_clusters; i++) {
        //This vector will be used to create a structure-type vector
        vector<vector<long int> > structure_tmp;
        
        //Convert cluster into structure
        if (clusters_cnt.size() && clusters_cnt[i].size()) {
            
            if (!Convert_cluster_to_structure(clusters_cnt[i], structure, structure_tmp)) {
                hout << "Error in Export_tecplot_files while converting cluster to structure." << endl;
                return 0;
            }
        }
        
        //Create string variable to store filename
        string filename, zone_name;
        
        //Create filename
        ostringstream number;
        number << i;
        filename = filename.append(type);
        filename = filename.append("_");
        filename = filename.append(number.str());
        zone_name = filename;
        filename = filename.append(".dat");
        
        //Export cluster to a file
        if (!tec360->Export_cnt_network_meshes(cub, i, clusters_gch, clusters_cnt, structure_tmp, points_in, radii, hybrid_particles, filename, zone_name)) {
            hout << "Error in Export_tecplot_files while translating and exporting directional clusters. filename =" << filename <<endl;
            return 0;
        }
    }
    

    //Export isolated CNTs, if any
    if (isolated.size()) {
        
        //Create a structure vector for isolated CNTs
        vector<vector<long int> > iso_structure;
        for (int i = 0; i < (int)isolated.size(); i++) {
            for (int j = 0 ; j < (int)isolated[i].size(); j++) {
                int CNT = isolated[i][j];
                iso_structure.push_back(structure[CNT]);
            }
        }
        
        //Generate tecplot files
        if ( !(tec360->Export_cnt_network_meshes(cub, points_in, radii, iso_structure, "Isolated_cnts.dat")) ) {
            hout << "Error in Export_tecplot_files while translating and exporting directional clusters. filename = Cluster_isolated.dat" <<endl;
            return 0;
        }
    }
    
    //Export isolated GNPs, if any
    if (isolated_gch.size()) {
        
        //Create a 1D vector of isolated GNPs
        vector<int> isolated_gch_1d;
        for (int i = 0; i < (int)isolated_gch.size(); i++) {
            for (int j = 0; j < (int)isolated_gch[i].size(); j++) {
                isolated_gch_1d.push_back(isolated_gch[i][j]);
            }
        }
        
        //If there are GNPs in the 1D vector, call the Tecplot function
        if (isolated_gch_1d.size()) {
            
            ofstream otec("Isolated_gnps.dat");
            otec << "TITLE = GNP_Meshes_Singlezone" << endl;
            otec << "VARIABLES = X, Y, Z" << endl;
            
            
            
            string zone_name = "GNPs_iso";
            if(!tec360->Export_randomly_oriented_gnps(otec, hybrid_particles, isolated_gch_1d, zone_name)) {
                hout << "Error in Export_tecplot_files while calling Export_randomly_oriented_gnps" << endl;
                return 0;
            }
            
            otec.close();
        }

    }

    
    //Variables to use the command line
    int s;
    char command[100];
    //Move the visualization files to a new folder
    
    //Check if clusters or percolated
    if (type == "Cluster") {
        s = sprintf(command, "mkdir clusters_%.4d", iter);
        s = system(command);
        s = sprintf(command, "mv Isolated*.dat Cluster*.dat clusters_%.4d", iter);
        system(command);
    } else if (type == "Percolated") {
        s = sprintf(command, "mkdir percolated_%.4d", iter);
        s = system(command);
        s = sprintf(command, "mv Isolated*.dat Percolated*.dat percolated_%.4d", iter);
        s = system(command);
    }
    
    
    //delete tecplot object
    delete tec360;
    
    return 1;
    
}
//Export tecplot files
int App_Network_3D::Export_triangulation_tecplot_files(const int &iter, const struct Geom_RVE &sample, const vector<Point_3D> &points_in, const vector<GCH> &hybrid_particles)const
{
    //Tecplot export object
    Tecplot_Export *tec360 = new Tecplot_Export;
    
    //Ecan every hybrid particle
    for (int i = 0; i < (int)hybrid_particles.size(); i++) {
        //Create string variable to store filename
        string filename_top, filename_bottom;
        
        //Create filename
        ostringstream number;
        number << i;
        filename_top = filename_top.append("Triangulation_");
        filename_top = filename_top.append(number.str());
        filename_top = filename_top.append(".dat");
        if (!tec360->Export_triangulation_network_threads(hybrid_particles[i], points_in, hybrid_particles[i].triangulation, filename_top)) {
            hout << "Error in Export_triangulation_tecplot_files" << endl;
            return 0;
        }
    }
    
    //Variables to use the command line
    int s;
    char command[100];
    //Move the visualization files to a new folder
    s = sprintf(command, "mkdir iter_%.4d", iter);
    s = system(command);
    s = sprintf(command, "mv Triangulation_*.dat iter_%.4d", iter);
    s = system(command);
    
    //delete tecplot object
    delete tec360;
    
    return 1;
    
}
//This function converts the data type index into data type structure
int App_Network_3D::Convert_index_to_structure(const vector<long int> &indices, vector<vector<long int> > &structure)const
{
    //Empty vector
    vector<long int> empty;
    //The branches are given in pairs
    for (int i = 0; i < (int)indices.size(); i=i+2) {
        structure.push_back(empty);
        for (long int j = indices[i]; j <= indices[i+1]; j++) {
            structure.back().push_back(j);
        }
    }
    return 1;
}
//This function converts the data type cluster (set of CNTs) into data type structure
int App_Network_3D::Convert_cluster_to_structure(const vector<int> &cluster, const vector<vector<long int> > &structure_in, vector<vector<long int> > &structure_out)const
{
    //Empty vector
    vector<int> empty;
    //The branches are given in pairs
    for (int i = 0; i < (int)cluster.size(); i++) {
        int CNT = cluster[i];
        structure_out.push_back(structure_in[CNT]);
    }
    return 1;
}
//===========================================================================
