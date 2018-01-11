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
    hout << "-_- Generate nanofiller network......" << endl;
    ct0 = time(NULL);
    GenNetwork *Genet = new GenNetwork;
    if(!Genet->Generate_nanofiller_network(Init->geom_rve, Init->cluster_geo, Init->nanotube_geo, Init->gnp_geo, Init->cutoff_dist, Init->tec360_flags, cnts_point, gnps_point, hybrid_particles, cnts_radius, cnts_structure, gnps_structure)) return 0;
    delete Genet;
    ct1 = time(NULL);
    hout << "Nanotube network generation time: " << (int)(ct1-ct0) <<" secs." << endl;
    //Printer *P = new Printer;
    //P->Print_1d_vec(gnps_point, "gnps_point_00.txt");
    //delete P;
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    ct0 = time(NULL);
    Background_vectors *Bckg = new Background_vectors;
    if (!Bckg->Generate_shells(Init->geom_rve, Init->nanotube_geo, cnts_point, hybrid_particles, shells_cnt, shells_gnps)) {
        hout << "Error when generating shells" << endl;
        return 0;
    }
    delete Bckg;
    ct1 = time(NULL);
    hout << "Generate shells and structure time: "<<(int)(ct1-ct0)<<" secs."<<endl;//*/
    
    //Geometry of observation window saved into a cuboid
    struct Geom_RVE window;
    
    for(int i=0; i<=Init->geom_rve.cut_num; i++)
    {
        hout << "============================================================================"<<endl;
        hout << "============================================================================"<<endl;
        hout << "Iteration " << i+1 << " of " << Init->geom_rve.cut_num+1 << endl;
        time_t it0, it1;
        it0 = time(NULL);
        
        //Update observation window geometry
        //Dimensions of the current observation window
        window.len_x = Init->geom_rve.win_max_x - i*Init->geom_rve.win_delt_x;
        window.wid_y = Init->geom_rve.win_max_y - i*Init->geom_rve.win_delt_y;
        window.hei_z = Init->geom_rve.win_max_z - i*Init->geom_rve.win_delt_z;
        //These variables are the coordinates of the lower corner of the observation window
        window.origin.x = Init->geom_rve.origin.x + (Init->geom_rve.len_x - window.len_x)/2;
        window.origin.y = Init->geom_rve.origin.y + (Init->geom_rve.wid_y - window.wid_y)/2;
        window.origin.z = Init->geom_rve.origin.z + (Init->geom_rve.hei_z - window.hei_z)/2;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //These vectors are used to export tecplot files
        vector<long int> empty;
        vector<vector<long int> > all_dead_indices(7,empty);
        vector<vector<long int> > all_percolated_indices(7,empty);
        vector<int> empty_int;
        vector<vector<int> > all_percolated_gnp(7,empty_int);
        vector<vector<int> > all_dead_gnps(7,empty_int);
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine the local networks in cutoff windows
        Cutoff_Wins *Cutwins = new Cutoff_Wins;
        //From this function I get the internal variables cnts_inside and boundary_cnt
        ct0 = time(NULL);
        if(!Cutwins->Extract_observation_window(i, Init->geom_rve, Init->nanotube_geo, Init->gnp_geo, hybrid_particles, cnts_structure, gnps_structure, cnts_radius, cnts_point, gnps_point, shells_cnt, shells_gnps)) {
            hout << "Error when extracting observation window "<< i+1 << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Extract observation window time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine the local networks inside the cutoff windows
        Contact_grid *Contacts = new Contact_grid;
        ct0 = time(NULL);
        if (Contacts->Generate_contact_grid(i, Init->geom_rve, Init->cutoff_dist, Init->nanotube_geo, Cutwins->cnts_inside, cnts_point, cnts_structure, Cutwins->gnps_inside, gnps_point, gnps_structure)==0) {
            hout << "Error when generating contact grid" << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Generate contact grid time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Hoshen-Kopelman algorithm
        Hoshen_Kopelman *HoKo = new Hoshen_Kopelman;
        ct0 = time(NULL);
        if (!HoKo->Determine_clusters(Init->geom_rve, Init->cutoff_dist, Cutwins->cnts_inside, Contacts->sectioned_domain, cnts_structure, cnts_point, cnts_radius, Cutwins->gnps_inside, Contacts->sectioned_domain_gnps, Contacts->sectioned_domain_hyb, gnps_structure, gnps_point, hybrid_particles)) {
            hout << "Error when determining clusters" << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Determine nanotube clusters time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //Contacts are not needed anymore, so delete the object
        delete Contacts;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Check if Tecplot visualization files were requested for clusters
        if (Init->tec360_flags.clusters) {
            ct0 = time(NULL);
            if (!Export_tecplot_files_for_clusters("Cluster", i, Init->tec360_flags.clusters, Init->geom_rve, cnts_point, cnts_radius, cnts_structure, HoKo->clusters_cnt, HoKo->isolated, hybrid_particles, HoKo->clusters_gch, HoKo->isolated_gch)) {
                hout << "Error when exporting Tecplot files for clusters" << endl;
                return 0;
            }
            ct1 = time(NULL);
            hout << "Export tecplot files for clusters time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        }
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Determine percolation
        Percolation *Perc = new Percolation;
        ct0 = time(NULL);
        if (!Perc->Determine_percolated_clusters(i, Init->geom_rve, Init->nanotube_geo, Init->gnp_geo, Cutwins->boundary_cnt, HoKo->labels, Cutwins->boundary_gnp, HoKo->labels_gnp, HoKo->clusters_cnt, HoKo->isolated, HoKo->clusters_gch, HoKo->isolated_gch)) {
            hout << "Error when determining percoalted clusters" << endl;
            return 0;
        }
        //Create vector to delete object and free memmory
        vector<int> family(Perc->family);
        delete Perc;
        ct1 = time(NULL);
        hout << "Determine percolating clusters time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //Labels form HK76 are not needed anymore
        HoKo->labels.clear();
        HoKo->labels_labels.clear();
        HoKo->labels_gnp.clear();
        HoKo->labels_labels_gnp.clear();
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Check if Tecplot visualization files were requested for percolated clusters
        if (Init->tec360_flags.percolated_clusters) {
            ct0 = time(NULL);
            if (!Export_tecplot_files_for_clusters("Percolated", i, Init->tec360_flags.percolated_clusters, Init->geom_rve, cnts_point, cnts_radius, cnts_structure, HoKo->clusters_cnt, HoKo->isolated, hybrid_particles, HoKo->clusters_gch, HoKo->isolated_gch)) {
                hout << "Error when exporting Tecplot files for percoalted clusters" << endl;
                return 0;
            }
            ct1 = time(NULL);
            hout << "Export tecplot files for percolated clusters time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        }
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //These vectors are used to store the fractions of the different families in the current observation window
        //families_lengths has 8 elements because of the 7 percolated families and the non-percoalted CNTs, the same is true for fractions
        vector<double> families_lengths(8,0);
        vector<double> fractions(8,0);
        vector<double> branches_lengths(7,0);
        
        //Loop over the different clusters so that the direct electrifying algorithm is aplied on each cluster
        Electrical_analysis *Electric_A = new Electrical_analysis;
        if (HoKo->clusters_cnt.size() || HoKo->clusters_gch.size()) {
            
            //Perform the electrical analysis to obtain the backbone and calculate the electrical resistance
            ct0 = time(NULL);
            if (!Electric_A->Perform_analysis_on_clusters(Init->simu_para.avoid_resistance, family, HoKo, Cutwins, cnts_structure, cnts_point, cnts_radius, gnps_structure, gnps_point, window, Init->electric_para, Init->cutoff_dist, hybrid_particles, all_dead_indices, all_percolated_indices, all_dead_gnps, all_percolated_gnp)) {
                hout << "Error when performing electrical analysis" << endl;
                return 0;
            }
            ct1 = time(NULL);
            hout << "Perform electrical analysis time: "<<(int)(ct1-ct0)<<" secs."<<endl;
            
        } else {
            
            //Calculate the resistances and resistivities of the polymer matrix
            vector<double> matrix_resistances;
            if (!Electric_A->Calculate_matrix_resistances(Init->electric_para.resistivity_matrix, window, matrix_resistances)) {
                hout << "Error when calculating matrix resistances for a sample without percolated clusters" << endl;
                return 0;
            }
            
            //Calculate the resistances and resistivities along each direction
            vector<double> resistivities;
            //paralel_resistors is initialized with three empty vectors
            vector<vector<double> > paralel_resistors(3,resistivities);
            if (!Electric_A->Calculate_resistances_and_resistivities(window, matrix_resistances, paralel_resistors, Electric_A->resistors, resistivities)) {
                hout << "Error when calculating matrix resistivities for a sample without percolated clusters" << endl;
                return 0;
            }
            
            //Append resistors to a file
            Printer *P = new Printer;
            P->Append_1d_vec(Electric_A->resistors, "resistors.txt");
            P->Append_1d_vec(resistivities, "resistivities.txt");
            delete P;
            
            hout << "There are no percolated clusters" << endl;
        }
        delete Electric_A;
        
        //Calculate the fractions of CNTs that belong to each family and save them to a file
        Clusters_fractions *Fracs = new Clusters_fractions;
        ct0 = time(NULL);
        if (!Fracs->Calculate_fractions(Init->geom_rve, Cutwins->cnts_inside, Cutwins->gnps_inside, cnts_structure, cnts_point, cnts_radius, HoKo->isolated, hybrid_particles, HoKo->isolated_gch, all_dead_indices, all_percolated_indices, all_dead_gnps, all_percolated_gnp)) {
            hout << "Error when calculating clusters fractions" << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Calculate fractions time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //Delete objects to free memory
        delete Fracs;
        delete Cutwins;

        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Check if Tecplot visualization files were requested for the backbone
        if (Init->tec360_flags.backbone) {
            ct0 = time(NULL);
            if (!Export_tecplot_files(i, Init->tec360_flags.backbone, Init->geom_rve, cnts_point, cnts_radius, hybrid_particles, HoKo->isolated_gch, cnts_structure, HoKo->isolated, all_dead_indices, all_percolated_indices, all_dead_gnps, all_percolated_gnp)) {
                hout << "Error when exporting tecplot files for backbone networks" << endl;
                return 0;
            }
            ct1 = time(NULL);
            hout << "Export tecplot files for backbone time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        }
        
        //Delete objects to free memory
        delete HoKo;
        
        //-----------------------------------------------------------------------------------------------------------------------------------------
        //Check if Tecplot visualization files were requested for the triangulations
        if (Init->tec360_flags.triangulations) {
            ct0 = time(NULL);
            if (!Export_triangulation_tecplot_files(i, Init->geom_rve, cnts_point, gnps_point, hybrid_particles)) {
                hout << "Error when exporting tecplot files for GNP triangulations" << endl;
                return 0;
            }
            ct1 = time(NULL);
            hout << "Export tecplot files for triangulations time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        }
        
        it1 = time(NULL);
        hout << "Iteration "<<i+1<<" time: "<<(int)(it1-it0)<<" secs."<<endl;

        
    }
    
    
    return 1;
}
//Export tecplot files
int App_Network_3D::Export_tecplot_files_for_clusters(const string &type, const int &iter, const int &tecplot_flag, const struct Geom_RVE &sample, const vector<Point_3D> &points_in, const vector<double> &radii, const vector<vector<long int> > &structure, const vector<vector<int> > &clusters_cnt, const vector<vector<int> > &isolated, const vector<GCH> &hybrid_particles, const vector<vector<int> > &clusters_gch, const vector<vector<int> > &isolated_gch)const
{
    
    //Tecplot export object
    Tecplot_Export *tec360 = new Tecplot_Export;
    
    //Geometry of observation window saved into a cuboid
    struct cuboid window;
    //Dimensions of the current observation window
    window.len_x = sample.win_max_x - iter*sample.win_delt_x;
    window.wid_y = sample.win_max_y - iter*sample.win_delt_y;
    window.hei_z = sample.win_max_z - iter*sample.win_delt_z;
    //These variables are the coordinates of the lower corner of the observation window
    window.poi_min.x = sample.origin.x + (sample.len_x - window.len_x)/2;
    window.poi_min.y = sample.origin.y + (sample.wid_y - window.wid_y)/2;
    window.poi_min.z = sample.origin.z + (sample.hei_z - window.hei_z)/2;
    
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
        
        //Check if meshes or wires (3D lines) were requested
        if (tecplot_flag == 1 || tecplot_flag == 3) {
            
            //Customize the filename for wires
            string filename_wires = "Wires_";
            filename_wires.append(filename);
            
            //Wires (3D lines) were requested
            if (!tec360->Export_network_3dlines(window, i, clusters_gch, clusters_cnt, structure_tmp, points_in, hybrid_particles, filename_wires, zone_name)) {
                hout << "Error in Export_tecplot_files_for_clusters while exporting clusters wires (3D lines). filename =" << filename <<endl;
                return 0;
            }
        }
        if (tecplot_flag == 2 || tecplot_flag == 3) {
            
            //Customize the filename for wires
            string filename_mesh = "SingleZone_";
            filename_mesh.append(filename);
            
            //Meshes were requested
            if (!tec360->Export_network_meshes(window, i, clusters_gch, clusters_cnt, structure_tmp, points_in, radii, hybrid_particles, filename_mesh, zone_name)) {
                hout << "Error in Export_tecplot_files_for_clusters while exporting clusters meshes. filename =" << filename <<endl;
                return 0;
            }
        }
        if (tecplot_flag < 0 || tecplot_flag > 3) {
            hout << "Error in Export_tecplot_files_for_clusters while exporting clusters. Invalid flag: " << tecplot_flag <<endl;
            hout << "Valid flags to export clusters are 1 or 2 only." << endl;
            return 0;
        }
        
    }
    
    //delete tecplot object
    delete tec360;
    
    //Export isolated particles
    if (!Export_isolated_particles(tecplot_flag, window, points_in, radii, structure, isolated, hybrid_particles, isolated_gch)) {
        hout << "Error in Export_tecplot_files while calling Export_isolated_particles" <<endl;
        return 0;
    }
    
    //Variables to use the command line
    int s;
    char command[100];
    //Move the visualization files to a new folder
    
    //Check if clusters or percolated
    if (type == "Cluster") {
        s = sprintf(command, "mkdir clusters_%.4d", iter);
        s = system(command);
        
        //Check if wires and/or meshes were requested
        if (tecplot_flag == 1 || tecplot_flag == 3) {
            s = sprintf(command, "mv Wires*.dat clusters_%.4d", iter);
            s = system(command);
        }
        if (tecplot_flag == 2 || tecplot_flag == 3) {
            s = sprintf(command, "mv Single*.dat clusters_%.4d", iter);
            s = system(command);
        }
    } else if (type == "Percolated") {
        s = sprintf(command, "mkdir percolated_%.4d", iter);
        s = system(command);
        
        //Check if wires and/or meshes were requested
        if (tecplot_flag == 1 || tecplot_flag == 3) {
            s = sprintf(command, "mv Wires*.dat percolated_%.4d", iter);
            s = system(command);
        }
        if (tecplot_flag == 2 || tecplot_flag == 3) {
            s = sprintf(command, "mv Single*.dat percolated_%.4d", iter);
            s = system(command);
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
//Export only the isolated particles
int App_Network_3D::Export_isolated_particles(const int &tecplot_flag, const struct cuboid &cub, const vector<Point_3D> &points_in, const vector<double> &radii, const vector<vector<long int> > &structure, const vector<vector<int> > &isolated, const vector<GCH> &hybrid_particles, const vector<vector<int> > &isolated_gch)const
{
    
    //Create a 1D vector of isolated CNTs
    vector<int> empty_int;
    vector<vector<int> > isolated_cnt_1d(1,empty_int);
    
    //Create a 1D vector of isolated GNPs
    vector<vector<int> > isolated_gch_1d(1,empty_int);
    
    //Create cluster of isolated CNTs, if any
    if (isolated.size()) {
        
        //Scan all isolated CNTs
        for (int i = 0; i < (int)isolated.size(); i++) {
            for (int j = 0 ; j < (int)isolated[i].size(); j++) {
                
                //Current CNT
                int CNT = isolated[i][j];
                
                //Add the CNT to make a single cluster of isolated CNTs
                isolated_cnt_1d[0].push_back(CNT);
            }
        }
        
    }
    
    //Create cluster of isolated isolated GNPs, if any
    if (isolated_gch.size()) {
        
        //Scan all isolated GNPs
        for (int i = 0; i < (int)isolated_gch.size(); i++) {
            for (int j = 0; j < (int)isolated_gch[i].size(); j++) {
                
                //Current GNP
                int GNP = isolated_gch[i][j];
                
                //Add the GNP to make a single cluster of isolated GNPs
                isolated_gch_1d[0].push_back(GNP);
            }
        }
        
    }
    
    //Define the filename and family name
    string family = "Isolated";
    
    //Tecplot export object
    Tecplot_Export *tec360 = new Tecplot_Export;
    
    //Structure vector for tecplot function
    vector<vector<long int> > structure_tmp;
    if (!Convert_cluster_to_structure(isolated_cnt_1d[0], structure, structure_tmp)) {
        hout << "Error in Export_tecplot_files while converting cluster to structure." << endl;
        return 0;
    }
    
    //Check if meshes or wires (3D lines) were requested
    if (tecplot_flag == 1 || tecplot_flag == 3) {
        
        //Define the filename and family name
        string filename = "Wires_Isolated_particles.dat";
        
        //Wires (3D lines) were requested
        if (!tec360->Export_network_3dlines(cub, 0, isolated_gch_1d, isolated_cnt_1d, structure_tmp, points_in, hybrid_particles, filename, family)) {
            hout << "Error in Export_tecplot_files while exporting isolated particles wires (3D lines)." <<endl;
            return 0;
        }
    }
    if (tecplot_flag == 2 || tecplot_flag == 3) {
        
        //Define the filename and family name
        string filename = "SingleZone_Isolated_particles.dat";
        
        //Meshes were requested
        if (!tec360->Export_network_meshes(cub, 0, isolated_gch_1d, isolated_cnt_1d, structure_tmp, points_in, radii, hybrid_particles, filename, family)) {
            hout << "Error in Export_tecplot_files while exporting isolated particles meshes." <<endl;
            return 0;
        }
    }
    if (tecplot_flag < 0 || tecplot_flag > 3){
        hout << "Error in Export_isolated_particles while exporting tecplot files for isolated particles. Invalid flag: " << tecplot_flag <<endl;
        hout << "Valid flags to export isolated particles are 1 or 2 only." << endl;
        return 0;
    }
    
    delete tec360;
    
    return 1;
}
//Export tecplot files
int App_Network_3D::Export_tecplot_files(const int &iter, const int &tecplot_flag, const struct Geom_RVE &sample, const vector<Point_3D> &cnts_points, const vector<double> &cnts_radius, const vector<GCH> &hybrid_particles, const vector<vector<int> > &isolated_gch, const vector<vector<long int> > &structure, const vector<vector<int> > &isolated, vector<vector<long int> > &all_dead_indices, const vector<vector<long int> > &all_percolated_indices, const vector<vector<int> > &all_dead_gnps, const vector<vector<int> > &all_percolated_gnp)const
{
    //These vectors will be used to create a structure-type vector to export to tecplot files
    vector<vector<long int> > percolated_structure_tmp, dead_structure_tmp;
    
    //Cluster vectors for CNTs
    vector<int> empty_int;
    vector<vector<int> > clusters_cnt_tmp(7, empty_int);
    
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
    vector<string> filenames_mesh, filenames_wires, family_names;
    Initialize_filenames(filenames_mesh, filenames_wires, family_names);
    
    //Export pecolated clusters and their dead branches
    for (int i = 0; i < 7; i++){
        //Check if the family is non empty. If it is non empty then a visualization file can be created
        if (all_percolated_indices[i].size()) {
            
            //Generate structure-type vectors
            if (!Convert_index_to_structure(all_percolated_indices[i], percolated_structure_tmp,clusters_cnt_tmp[i])) {
                hout << "Error in Export_tecplot_files when converting percolated indices to structure"<<endl;
                return 0;
            }
            
            //Generate Tecplot files, check if meshes or wires (3D lines) were requested
            if (tecplot_flag == 1 || tecplot_flag == 3) {
                
                //Wires (3D lines) were requested
                if ( !tec360->Export_network_3dlines(cub,i, all_percolated_gnp, clusters_cnt_tmp,percolated_structure_tmp, cnts_points, hybrid_particles, filenames_wires[i], family_names[i]) )  {
                    hout << "Error in Export_tecplot_files when calling Export_network_meshes for backbone wires (3D lines)" <<endl;
                    return 0;
                }
            }
            if (tecplot_flag == 2 || tecplot_flag == 3) {
                
                //Meshes were requested
                if ( !tec360->Export_network_meshes(cub,i, all_percolated_gnp, clusters_cnt_tmp,percolated_structure_tmp, cnts_points, cnts_radius, hybrid_particles, filenames_mesh[i], family_names[i]) )  {
                    hout << "Error in Export_tecplot_files when calling Export_network_meshes for backbone meshes" <<endl;
                    return 0;
                }
            }
            if (tecplot_flag < 0 || tecplot_flag > 3) {
                hout << "Error in Export_tecplot_files while exporting tecplot files for backbone. Invalid flag: " << tecplot_flag <<endl;
                hout << "Valid flags to export backbone and dead branches are 1 or 2 only." << endl;
                return 0;
            }
            
            //Clear the variable for the next iteration
            percolated_structure_tmp.clear();
        }
        
        //It is possible that a CNT spans from one boundary to the other and is not in contact with other CNTs.
        //In this case, the CNT is a cluster itself and has no dead branches,
        //so I need to check separately if all_dead_indices[i] is non empty,
        //i.e., the fact that a cluster percolates does not mean that there will be dead branches
        if (all_dead_indices[i].size()) {
            
            //Generate structure-type vectors
            if (!Convert_index_to_structure(all_dead_indices[i], dead_structure_tmp,clusters_cnt_tmp[i])) {
                hout << "Error in Export_tecplot_files when converting dead indices to structure"<<endl;
                return 0;
            }
            
            //Generate Tecplot files, check if meshes or wires (3D lines) were requested
            if (tecplot_flag == 1 || tecplot_flag == 3) {
                
                //Wires (3D lines) were requested
                if ( !tec360->Export_network_3dlines(cub,i, all_dead_gnps, clusters_cnt_tmp, dead_structure_tmp, cnts_points, hybrid_particles, filenames_wires[i+7], family_names[i+7]) )  {
                    hout << "Error in Export_tecplot_files when calling Export_network_meshes for dead branches wires (3D lines)" <<endl;
                    return 0;
                }
            }
            if (tecplot_flag == 2 || tecplot_flag == 3) {
                
                //Meshes were requested
                if ( !tec360->Export_network_meshes(cub,i, all_dead_gnps,clusters_cnt_tmp, dead_structure_tmp, cnts_points, cnts_radius, hybrid_particles, filenames_mesh[i+7], family_names[i+7]) ) {
                    hout << "Error in Export_tecplot_files when calling Export_network_meshes for dead branches meshes" <<endl;
                    return 0;
                }
            }
            if (tecplot_flag < 0 || tecplot_flag > 3) {
                hout << "Error in Export_tecplot_files while exporting tecplot files for dead branches. Invalid flag: " << tecplot_flag <<endl;
                hout << "Valid flags to export backbone and dead branches are 1 or 2 only." << endl;
                return 0;
            }
            
            //Clear the variable for the next iteration
            dead_structure_tmp.clear();
        }
        
    }
    
    //delete tecplot object
    delete tec360;
    
    //Export isolated particles
    if (!Export_isolated_particles(tecplot_flag, cub, cnts_points, cnts_radius, structure, isolated, hybrid_particles, isolated_gch)) {
        hout << "Error in Export_tecplot_files while calling Export_isolated_particles" <<endl;
        return 0;
    }
    
    //Variables to use the command line
    int s;
    char command[100];
    //Move the visualization files to a new folder
    s = sprintf(command, "mkdir iter_%.4d", iter);
    s = system(command);
    
    //Check if wires and/or meshes were requested
    if (tecplot_flag == 1 || tecplot_flag == 3) {
        s = sprintf(command, "mv Wires*.dat iter_%.4d", iter);
        s = system(command);
    }
    if (tecplot_flag == 2 || tecplot_flag == 3) {
        s = sprintf(command, "mv Single*.dat iter_%.4d", iter);
        s = system(command);
    }
    
    
    return 1;
    
}
//Filenames for the clusters
void App_Network_3D::Initialize_filenames(vector<string> &filenames_mesh, vector<string> &filenames_wires, vector<string> &family_names)const
{
    //Mesh files
    //Save a separate file for each percolating direction
    filenames_mesh.push_back("SingleZone_00_X.dat");
    filenames_mesh.push_back("SingleZone_01_Y.dat");
    filenames_mesh.push_back("SingleZone_02_Z.dat");
    filenames_mesh.push_back("SingleZone_03_XY.dat");
    filenames_mesh.push_back("SingleZone_04_XZ.dat");
    filenames_mesh.push_back("SingleZone_05_YZ.dat");
    filenames_mesh.push_back("SingleZone_06_XYZ.dat");
    //Save a separate file for the dead branches of each percolating direction
    filenames_mesh.push_back("SingleZoneDead_00_X.dat");
    filenames_mesh.push_back("SingleZoneDead_01_Y.dat");
    filenames_mesh.push_back("SingleZoneDead_02_Z.dat");
    filenames_mesh.push_back("SingleZoneDead_03_XY.dat");
    filenames_mesh.push_back("SingleZoneDead_04_XZ.dat");
    filenames_mesh.push_back("SingleZoneDead_05_YZ.dat");
    filenames_mesh.push_back("SingleZoneDead_06_XYZ.dat");
    
    //3D Line files
    //Save a separate file for each percolating direction
    filenames_wires.push_back("Wires_00_X.dat");
    filenames_wires.push_back("Wires_01_Y.dat");
    filenames_wires.push_back("Wires_02_Z.dat");
    filenames_wires.push_back("Wires_03_XY.dat");
    filenames_wires.push_back("Wires_04_XZ.dat");
    filenames_wires.push_back("Wires_05_YZ.dat");
    filenames_wires.push_back("Wires_06_XYZ.dat");
    //Save a separate file for the dead branches of each percolating direction
    filenames_wires.push_back("WiresDead_00_X.dat");
    filenames_wires.push_back("WiresDead_01_Y.dat");
    filenames_wires.push_back("WiresDead_02_Z.dat");
    filenames_wires.push_back("WiresDead_03_XY.dat");
    filenames_wires.push_back("WiresDead_04_XZ.dat");
    filenames_wires.push_back("WiresDead_05_YZ.dat");
    filenames_wires.push_back("WiresDead_06_XYZ.dat");
    
    //Family names
    family_names.push_back("X");
    family_names.push_back("Y");
    family_names.push_back("Z");
    family_names.push_back("XY");
    family_names.push_back("XZ");
    family_names.push_back("YZ");
    family_names.push_back("XYZ");
    family_names.push_back("X_dead");
    family_names.push_back("Y_dead");
    family_names.push_back("Z_dead");
    family_names.push_back("XY_dead");
    family_names.push_back("XZ_dead");
    family_names.push_back("YZ_dead");
    family_names.push_back("XYZ_dead");
    
}
//Export tecplot files
int App_Network_3D::Export_triangulation_tecplot_files(const int &iter, const struct Geom_RVE &sample, const vector<Point_3D> &cnts_point, const vector<Point_3D> &gnps_point, const vector<GCH> &hybrid_particles)const
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
        if (!tec360->Export_triangulation_network_3dlines(hybrid_particles[i], cnts_point, gnps_point, hybrid_particles[i].triangulation, filename_top)) {
            hout << "Error in Export_triangulation_tecplot_files" << endl;
            return 0;
        }
    }
    //delete tecplot object
    delete tec360;
    
    //Variables to use the command line
    int s;
    char command[100];
    //Move the visualization files to a new folder
    s = sprintf(command, "mkdir iter_%.4d", iter);
    s = system(command);
    s = sprintf(command, "mv Triangulation_*.dat iter_%.4d", iter);
    s = system(command);
    
    
    return 1;
    
}
//This function converts the data type index into data type structure
int App_Network_3D::Convert_index_to_structure(const vector<long int> &indices, vector<vector<long int> > &structure, vector<int> &cluster)const
{
    //Empty vector
    vector<long int> empty;
    
    //Clear the cluster vector
    cluster.clear();
    
    //The branches are given in pairs
    for (int i = 0; i < (int)indices.size(); i=i+2) {
        
        //Add the CNT number to the cluster, which is the size of the structure
        cluster.push_back((int)structure.size());
        
        //Add a new CNT to the structure
        structure.push_back(empty);
        
        //Fill the points of the CNT
        for (long int j = indices[i]; j <= indices[i+1]; j++) {
            structure.back().push_back(j);
        }
    }
    
    
    
    return 1;
}
//===========================================================================
