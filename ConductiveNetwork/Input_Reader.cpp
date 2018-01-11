//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Input_Reader.cpp
//OBJECTIVE:	Reading the input data
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa;  angel.mora@kaust.edu.sa
//====================================================================================

#include "Input_Reader.h"

//---------------------------------------------------------------------------
//Read data
int Input::Read_Infile(ifstream &infile)
{

	cout << "Reading input file..." << endl;
	hout << "Reading input file..." << endl;

	while(!infile.eof())
	{
		istringstream istr(Get_Line(infile));
		if(infile.eof()) break;
		string str_temp;
		istr >> str_temp;

		if(str_temp.empty()) continue;  //skip over the space or Enter key after every input item
		else if(str_temp=="Application_Name") { if(Read_application_name(app_name, infile)==0) return 0; }
		else if(str_temp=="Simulation_Parameters")	{ if(Read_simulation_parameters(simu_para, infile)==0) return 0; }
		else if(str_temp=="RVE_Geometry")	{ if(Read_rve_geometry(geom_rve, infile)==0) return 0; }
		else if(str_temp=="Nanotube_Geometry")	{ if(Read_nanotube_geo_parameters(nanotube_geo, infile)==0) return 0; }
		else if(str_temp=="Cluster_Geometry")	{ if(Read_cluster_geo_parameters(cluster_geo, infile)==0) return 0; }
		else if(str_temp=="Cutoff_Distances")	{ if(Read_cutoff_distances(cutoff_dist, infile)==0) return 0; }
        else if(str_temp=="GNP_Geometry")	{ if(Read_gnp_geo_parameters(gnp_geo, infile)==0) return 0; }
		else if(str_temp=="Electrical_Parameters")	{ if(Read_electrical_parameters(electric_para, infile)==0) return 0; }
        else if(str_temp=="Tecplot_flags")	{ if(Read_tecplot_flags(tec360_flags, infile)==0) return 0; }
		else
		{ 
			cout << "Error: the keywords \"" << str_temp << "\" is not defined!" << endl; 
			hout << "Error: the keywords \"" << str_temp << "\" is not defined!" << endl; 
			return 0; 
		}

		//the real volume of cnts in the RVE
		if(geom_rve.mark&&nanotube_geo.mark)
		{
			//Calculate the real volume of nanotubes
			nanotube_geo.real_volume = nanotube_geo.volume_fraction * geom_rve.volume;
			//Define the extented RVE with one length of nanotube
			geom_rve.ex_origin.x = geom_rve.origin.x - nanotube_geo.len_max;
			geom_rve.ex_origin.y = geom_rve.origin.y - nanotube_geo.len_max;
			geom_rve.ex_origin.z = geom_rve.origin.z - nanotube_geo.len_max;
			geom_rve.ex_len = geom_rve.len_x + 2*nanotube_geo.len_max;
			geom_rve.ey_wid = geom_rve.wid_y+ 2*nanotube_geo.len_max;
			geom_rve.ez_hei = geom_rve.hei_z + 2*nanotube_geo.len_max;
		}
	}

	cout << "Reading the keywords is finished!" << endl;
	hout << "Reading the keywords is finished!" << endl;

	if(!app_name.mark) { cout << "Attention: \"Application_Name\" will use default parameters!" << endl; hout << "Attention: \"Application_Name\" will use default parameters!" << endl; }
	if(!simu_para.mark) { cout << "Attention: \"Simulation_Parameters\" will use default parameters!" << endl; hout << "Attention: \"Simulation_Parameters\" will use default parameters!" << endl; }
	if(!geom_rve.mark) { cout << "Attention: \"RVE_Geometry\" will use default parameters!" << endl; hout << "Attention: \"RVE_Geometry\" will use default parameters!" << endl; }
	if(!nanotube_geo.mark) {	cout << "Attention: \"Nanotube_Geometry\" will use default parameters!" << endl; hout << "Attention: \"Nanotube_Geometry\" will use default parameters!" << endl; }
	if(!cluster_geo.mark) { cout << "Attention: \"Cluster_Geometry\" will use default parameters!" << endl; hout << "Attention: \"Cluster_Geometry\" will use default parameters!" << endl; }
	if(!cutoff_dist.mark) {	cout << "Attention: \"Cutoff_Distances\" will use default parameters!" << endl; hout << "Attention: \"Cutoff_Distances\" will use default parameters!" << endl; }
    if(!gnp_geo.mark) {	cout << "Attention: \"GNP_Geometry\" will use default parameters!" << endl; hout << "Attention: \"GNP_Geometry\" will use default parameters!" << endl; }
	if(!electric_para.mark) {	cout << "Attention: \"Electrical_Parameters\" will use default parameters!" << endl; hout << "Attention: \"Electrical_Parameters\" will use default parameters!" << endl; }

	return 1;
}
//---------------------------------------------------------------------------
//Initialize data
int Input::Data_Initialization()
{
	//Initialize name of simulation
	app_name.keywords = "Application_Name";
	app_name.mark = false;
	app_name.str = "App_Electrical_Network_3D";

	//Initialize paramters of simulation
	simu_para.keywords = "Simulation_Parameters";
	simu_para.mark = false;
	simu_para.simu_name = "Test";
	simu_para.sample_num = 1;
	simu_para.create_read_network = "Create_Network";
    simu_para.avoid_resistance = 0;
    

	//Initialize the geometric parameters of the RVE
	geom_rve.keywords = "RVE_Geometry";
	geom_rve.mark = false;
    geom_rve.particle_type = "CNT_wires";
	geom_rve.origin.x = 0.0;
	geom_rve.origin.y = 0.0;
	geom_rve.origin.z = 0.0;
	geom_rve.origin.flag = 0;
	geom_rve.len_x = 1.0;
	geom_rve.wid_y = 1.0;
	geom_rve.hei_z = 1.0;
	geom_rve.ex_origin.x = 0.0;
	geom_rve.ex_origin.y = 0.0;
	geom_rve.ex_origin.z = 0.0;
	geom_rve.ex_origin.flag = 0;
	geom_rve.ex_len = 1.0;
	geom_rve.ey_wid = 1.0;
	geom_rve.ez_hei = 1.0;
	geom_rve.volume = geom_rve.len_x*geom_rve.wid_y*geom_rve.hei_z;
	geom_rve.density = 1.0;
	geom_rve.gs_minx = 1.0;
	geom_rve.gs_miny = 1.0;
	geom_rve.gs_minz = 1.0;
	geom_rve.win_max_x = 1.0;
	geom_rve.win_max_y = 1.0;
	geom_rve.win_max_z = 1.0;
	geom_rve.win_min_x = 1.0;
	geom_rve.win_min_y = 1.0;
	geom_rve.win_min_z = 1.0;
	geom_rve.win_delt_x = 1.0;
	geom_rve.win_delt_y = 1.0;
	geom_rve.win_delt_z = 1.0;

	//Initialize the geometric paramters of nanotubes
	nanotube_geo.keywords = "Nanotube_Geometry";
	nanotube_geo.mark = false;
	nanotube_geo.dir_distrib_type = "random";
	nanotube_geo.ini_pha = 0.0;
	nanotube_geo.ini_sita = 0.0;
	nanotube_geo.angle_max = 1.5707963267948966;
	nanotube_geo.step_length = 0.01;
	nanotube_geo.len_distrib_type = "uniform";
	nanotube_geo.len_min = 1.0;
	nanotube_geo.len_max = 1.0;
	nanotube_geo.rad_distrib_type = "uniform";
	nanotube_geo.rad_min = 0.005;
	nanotube_geo.rad_max = 0.005;
	nanotube_geo.criterion = "vol";
	nanotube_geo.volume_fraction = 0.0;
	nanotube_geo.accum_mode = 0;
	nanotube_geo.real_volume = 0.0;
	nanotube_geo.weight_fraction = 0.0;
	nanotube_geo.real_weight = 0.0;
	nanotube_geo.matrix_density = 1.06;
	nanotube_geo.linear_density = 5.8E-5;

	//Initialize the geometric paramters of nanotube clusters
	cluster_geo.keywords = "Cluster_Geometry";
	cluster_geo.mark = false;
	cluster_geo.print_key = 1;
	cluster_geo.vol_fra_criterion = 0;
	cluster_geo.amin = 0.0;
	cluster_geo.amax = 0.0;
	cluster_geo.bmin = 0.0;
	cluster_geo.cmin = 0.0;
	cluster_geo.growth_probability = 0.0;
	cluster_geo.volf_clust = 0.0;
	cluster_geo.cnt_real_volume = 0.0;

	//Initialize cutoff distances
	cutoff_dist.keywords = "Cutoff_Distances";
	cutoff_dist.mark = false;
	cutoff_dist.tunneling_dist = 0.0018;
	cutoff_dist.van_der_Waals_dist = 0.00034;

	//Initialize electrical parameters
	electric_para.keywords = "Electrical_Parameters";
	electric_para.mark = false;
	electric_para.applied_voltage = 1.0;
	electric_para.resistivity_CF = 0.001;
    
    //Initialize tecplot flags (do not export anything)
    tec360_flags.keywords = "Tecplot_flags";
    tec360_flags.mark = false;
    tec360_flags.generated_cnts = 0;
    tec360_flags.generated_gnps = 0;
    tec360_flags.clusters = 0;
    tec360_flags.percolated_clusters = 0;
    tec360_flags.backbone = 0;
    tec360_flags.triangulations = 0;

	cout << "^_^ Data initialization achieves" <<endl<<endl;
	hout << "^_^ Data initialization achieves" <<endl<<endl;

	return 1;
}
//---------------------------------------------------------------------------
//Reading the name of application case
int Input::Read_application_name(struct App_name &app_name, ifstream &infile)
{
	if(app_name.mark)
	{
		cout << "Attention: \"" << app_name.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << app_name.keywords << "\" has been input!" << endl;
		return 0;
	}
	else app_name.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> app_name.str;

	return 1;
}
//---------------------------------------------------------------------------
//Reading the parameters of simulation
int Input::Read_simulation_parameters(struct Simu_para &simu_para, ifstream &infile)
{
	if(simu_para.mark)
	{
		cout << "Attention: \"" << simu_para.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << simu_para.keywords << "\" has been input!" << endl;
		return 0;
	}
	else simu_para.mark = true;

	istringstream istr0(Get_Line(infile));
	istr0 >> simu_para.simu_name;			//Read the name of simulation

	istringstream istr1(Get_Line(infile));
	istr1 >> simu_para.sample_num;			//Read the number of samples
	if(simu_para.sample_num<1)	 {	hout << "Error: the number of samples less than 1." << endl; return 0; }
    
	istringstream istr2(Get_Line(infile));
	istr2 >> simu_para.create_read_network;		//Read keyword for creating a new network or reading a network from a file
	if(simu_para.create_read_network!="Create_Network"&&simu_para.create_read_network!="Read_Network"&&simu_para.create_read_network!="Network_From_Seeds")
	{ hout << "Error: Invalid keyword for 'create_read_network'. Valid options are 'Create_Network', 'Read_Network', or 'Network_From_Seeds'. Input was:" << simu_para.create_read_network << endl; return 0; }
    
    //Flag to avoid calculating the resistance of the network
    //This can be useful when only a geometric study or visualizations are needed
    // 1: Avoid calculating the resistance of the network
    // 0: Calculate the resistance of the network
    istringstream istr3(Get_Line(infile));
    istr3 >> simu_para.avoid_resistance;

	return 1;
}
//---------------------------------------------------------------------------
//Reading geometric information of the RVE
int Input::Read_rve_geometry(struct Geom_RVE &geom_rve, ifstream &infile)
{
	if(geom_rve.mark)
	{
		cout << "Attention: \"" << geom_rve.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << geom_rve.keywords << "\" has been input!" << endl;
		return 0;
	}
	else geom_rve.mark = true;

    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Define the type of particles inside the RVE: CNT_wires for CNTs or Hybrid_particles for GNP-CNT hybrid particles
    istringstream istr_perticle_type(Get_Line(infile));
    istr_perticle_type >> geom_rve.particle_type;
    if (geom_rve.particle_type != "CNT_wires" && geom_rve.particle_type != "GNP_cuboids" && geom_rve.particle_type != "Hybrid_particles" && geom_rve.particle_type != "GNP_CNT_mix") {
        cout << "Error: the type of particles shoud be one of the following: CNT_wires, GNP_cuboids, Hybrid_particles or GNP_CNT_mix." << endl;
        hout << "Error: the type of particles shoud be one of the following: CNT_wires, GNP_cuboids, Hybrid_particles or GNP_CNT_mix." << endl;
        return 0;
    }
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //If the network is generated from some seeds, then read them
    if (simu_para.create_read_network=="Network_From_Seeds") {
        
        //Different types of fillers need different number of seed
        //Thus, check the filler type and read the necessary number of seeds
        if (geom_rve.particle_type=="CNT_wires") {
            
            //7 seeds are required for CNT_wires
            geom_rve.network_seeds.assign(7, 0);
            
            //Read the line with the seeds
            istringstream istr_network_seeds(Get_Line(infile));
            //Add seeds to the vector
            for (int i = 0; i < 7; i++) {
                istr_network_seeds >> geom_rve.network_seeds[i];
                //hout << geom_rve.network_seeds[i] << endl;
            }
        }
    }
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the domain of RVE: the lower-left corner point of RVE and the length, width and height of RVE
	istringstream istr0(Get_Line(infile));
	istr0 >> geom_rve.origin.x >> geom_rve.origin.y >> geom_rve.origin.z;
	istr0 >> geom_rve.len_x >> geom_rve.wid_y >> geom_rve.hei_z;
	if(geom_rve.len_x<=0||geom_rve.wid_y<=0||geom_rve.hei_z<=0)
	{
		cout << "Error: the sizes of RVE should be positive!" << endl;
		hout << "Error: the sizes of RVE should be positive!" << endl;
		return 0;
	}
	geom_rve.volume = geom_rve.len_x*geom_rve.wid_y*geom_rve.hei_z;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the size range of the observation window and descrement by every step in x, y and z directions
	istringstream istr1(Get_Line(infile));
	istr1 >> geom_rve.win_max_x >> geom_rve.win_max_y >> geom_rve.win_max_z;
	istringstream istr2(Get_Line(infile));
	istr2 >> geom_rve.win_delt_x >> geom_rve.win_delt_y >> geom_rve.win_delt_z;
	istringstream istr3(Get_Line(infile));
	istr3 >> geom_rve.win_min_x >> geom_rve.win_min_y >> geom_rve.win_min_z;

	if(geom_rve.win_max_x<=0.0||geom_rve.win_max_y<=0.0||geom_rve.win_max_z<=0.0||
	   geom_rve.win_max_x>geom_rve.len_x||geom_rve.win_max_y>geom_rve.wid_y||geom_rve.win_max_y>geom_rve.hei_z)
	{
		cout << "Error: the win_max in each direction of RVE should be positive and must be smaller than the size of RVE." << endl;
		hout << "Error: the win_max in each direction of RVE should be positive and must be smaller than the size of RVE." << endl;
		return 0;
	}
	if(geom_rve.win_min_x<=0.0||geom_rve.win_min_y<=0.0||geom_rve.win_min_z<=0.0||
	   geom_rve.win_min_x>geom_rve.win_max_x||geom_rve.win_min_y>geom_rve.win_max_y||geom_rve.win_min_z>geom_rve.win_max_z)
	{
		cout << "Error: the win_min in each direction of RVE should be positive and must be smaller than max." << endl;
		hout << "Error: the win_min in each direction of RVE should be positive and must be smaller than max." << endl;
		return 0;
	}
	if(geom_rve.win_delt_x<=0.0||geom_rve.win_delt_y<=0.0||geom_rve.win_delt_z<=0.0)
	{
		cout << "Error: the win_delt in each direction of RVE should be positive." << endl;
		hout << "Error: the win_delt in each direction of RVE should be positive." << endl;
		return 0;
	}

	//Details: +Zero for reducing the error of division
	int num[3] = {	(int)((geom_rve.win_max_x-geom_rve.win_min_x + Zero)/geom_rve.win_delt_x), 
							(int)((geom_rve.win_max_y-geom_rve.win_min_y + Zero)/geom_rve.win_delt_y), 
							(int)((geom_rve.win_max_z-geom_rve.win_min_z + Zero)/geom_rve.win_delt_z)	};

	if(num[0]!=num[1]||num[0]!=num[2])
	{
		cout << "Error: the number of steps is different on each direction: " << endl;
        cout << "\tSteps on x = " << num[0] << endl;
        cout << "\tSteps on y = " << num[1] << endl;
        cout << "\tSteps on z = " << num[2] << endl;
		hout << "Error: the number of steps is different on each direction: " << endl;
        hout << "\tSteps on x = " << num[0] << endl;
        hout << "\tSteps on y = " << num[1] << endl;
        hout << "\tSteps on z = " << num[2] << endl;
		return 0;
	}
    else
    {
        //This code might be able to make the observation windows work when
        //the step size is non-integer multiple of (win_max-win_min)
        //However somewhere there is some kind of memory error
        //---AMC Aug 4, 2017
        /*/
        //Check that the last step is the smallest observation window
        int test_x = abs(geom_rve.win_max_x-num[0]*geom_rve.win_delt_x - geom_rve.win_min_x) > Zero;
        int test_y = abs(geom_rve.win_max_y-num[1]*geom_rve.win_delt_y - geom_rve.win_min_y) > Zero;
        int test_z = abs(geom_rve.win_max_z-num[2]*geom_rve.win_delt_z - geom_rve.win_min_z) > Zero;
        
        
        if (test_x || test_y || test_z) {
            //If along any direction the last step results in an observation window smaller than the minimum, then add one more
            num[0] = num[0] + 1;
            hout << "Increased number of observation windows" <<endl;
            return 0;
        }//*/
        geom_rve.cut_num = num[0];
        
    }

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the minimum size for background grids (looking for contact points)
	istringstream istr4(Get_Line(infile));
	istr4 >> geom_rve.gs_minx >> geom_rve.gs_miny >> geom_rve.gs_minz;
	if(geom_rve.gs_minx<=0||geom_rve.gs_miny<=0||geom_rve.gs_minz<=0)
	{
		cout << "Error: the number of segments in each direction of RVE should be positive!" << endl;
		hout << "Error: the number of segments in each direction of RVE should be positive" << endl;
		return 0;
	}
	else if((int)(geom_rve.win_max_x/geom_rve.gs_minx)>500||
			  (int)(geom_rve.win_max_y/geom_rve.gs_miny)>500||
			  (int)(geom_rve.win_max_z/geom_rve.gs_minz)>500)
	{
		cout << "Error: the number of divisions in one of boundary is too big (>500), which leads to the memory problem!" << endl;
		hout << "Error: the number of divisions in one of boundary is too big (>500), which leads to the memory problem!" << endl;
		return 0;	
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading the geometric parameters of nanotubes
int Input::Read_nanotube_geo_parameters(struct Nanotube_Geo &nanotube_geo, ifstream &infile)
{
	if(nanotube_geo.mark)
	{
		cout << "Attention: \"" << nanotube_geo.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << nanotube_geo.keywords << "\" has been input!" << endl;
		return 0;
	}
	else nanotube_geo.mark = true;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the initial growth direction type (random or specific) in a RVE
	istringstream istr_initial_direction(Get_Line(infile));
	istr_initial_direction >> nanotube_geo.dir_distrib_type;
	if(nanotube_geo.dir_distrib_type!="random"&&nanotube_geo.dir_distrib_type!="specific"){ hout << "Error: the direction distribution type must be either random or specific." << endl;	return 0; }
	if(nanotube_geo.dir_distrib_type=="specific")
	{
		istr_initial_direction  >> nanotube_geo.ini_sita >> nanotube_geo.ini_pha;
		if(nanotube_geo.ini_sita<0||nanotube_geo.ini_sita>PI||nanotube_geo.ini_pha<0||nanotube_geo.ini_pha>=2*PI)
		{
			hout << "Error: the specified angle is not in the acceptable range of (0, 2PI)." << endl;
			return 0;
		}
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the normal distribution range [-omega, omega] of the growth direction
	istringstream istr_angle_range(Get_Line(infile));
	istr_angle_range >> nanotube_geo.angle_max;
	if(nanotube_geo.angle_max>0.5*PI){ hout << "Error: the specified angle is not in the acceptable range of (-PI/2, PI/2)." << endl;	 return 0; }

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the step length (unit: micromether) of nanotube growth
	istringstream istr_step_len(Get_Line(infile));
	istr_step_len >> nanotube_geo.step_length;
	if(nanotube_geo.step_length<=0||
	   nanotube_geo.step_length>=0.25*geom_rve.len_x||
	   nanotube_geo.step_length>=0.25*geom_rve.wid_y||
       nanotube_geo.step_length>=0.25*geom_rve.hei_z)
	{ hout << "Error: the step length must be positive and 0.25 times lesser than the dimension of the RVE box." << endl;	return 0; }

	//-----------------------------------------------------------------------------------------------------------------------------------------
    //Define the distribution type (uniform or normal) of the length (unit: micromether) of nanotubes and the length range (min, max) of nanotubes in a RVE
    istringstream istr_cnt_len(Get_Line(infile));
    istr_cnt_len >> nanotube_geo.len_distrib_type;
    if(nanotube_geo.len_distrib_type!="uniform"&&nanotube_geo.len_distrib_type!="normal"){ hout << "Error: the distribution of the length should be either normal or uniform." << endl;	return 0; }
    istr_cnt_len >> nanotube_geo.len_min >> nanotube_geo.len_max;
    if(nanotube_geo.len_min<0||nanotube_geo.len_max<0||nanotube_geo.len_max<nanotube_geo.len_min){ hout << "Error: the length must be non-negative and min must be smaller than max." << endl; return 0; }

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the distribution type (uniform or normal) of the radius (unit: micromether) of nanotubes and the radius range (min, max) of nanotubes in a RVE
    istringstream istr_cnt_rad(Get_Line(infile));
    istr_cnt_rad >> nanotube_geo.rad_distrib_type;
    if(nanotube_geo.rad_distrib_type!="uniform"&&nanotube_geo.rad_distrib_type!="normal"){ hout << "Error: the distribution of the radius should be either normal or uniform." << endl;	return 0; }
    istr_cnt_rad >> nanotube_geo.rad_min >> nanotube_geo.rad_max;
    if(nanotube_geo.rad_min<0||nanotube_geo.rad_max<0||nanotube_geo.rad_max<nanotube_geo.rad_min||
	   nanotube_geo.rad_min>3*nanotube_geo.step_length||nanotube_geo.rad_max>0.05*nanotube_geo.len_min)
	{ hout << "Error: the radius must be non-negative, min must be smaller than max, min must be smaller than 3*step_length and max must be smaller than 0.05*len_min." << endl; return 0; }

    //---------------------------------------------------------------------------------------
    //-------- AMC
    //
    //Define the density of the CNTs (unit: gm/cm3)
    istringstream istr_density(Get_Line(infile));
    istr_density >> nanotube_geo.density;
    if (nanotube_geo.density <= 0) {
        hout << "Error: the density has to be a value greater than zero. Input value:" << nanotube_geo.density;
        hout << "." << endl;
        return 0;
    }
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
	//Define the volume or weight fraction of nanotubes in the RVE
	istringstream istr_cnt_vol(Get_Line(infile));
	istr_cnt_vol >> nanotube_geo.criterion;
	if(nanotube_geo.criterion=="vol")
	{
		istr_cnt_vol >> nanotube_geo.volume_fraction;
		if(nanotube_geo.volume_fraction>1||nanotube_geo.volume_fraction<0){ hout << "Error: the volume fraction must be between 0 and 1." << endl; return 0; }
		hout << "    The CNT volume fraction is "<< nanotube_geo.volume_fraction << endl;
        
		istr_cnt_vol >> nanotube_geo.accum_mode;
		if(nanotube_geo.accum_mode<0&&nanotube_geo.accum_mode>2){ hout <<"Error: the mode of accumulation should be between 0 and 2." << endl; return 0; }
        
		//The total volume of the nanotube network
		nanotube_geo.real_volume = nanotube_geo.volume_fraction*geom_rve.volume;
	}
	else if(nanotube_geo.criterion=="wt")
	{
		istr_cnt_vol >> nanotube_geo.weight_fraction;
		if(nanotube_geo.weight_fraction>1||nanotube_geo.weight_fraction<0){ hout << "Error: the volume fraction must be between 0 and 1." << endl; return 0; }
		hout << "    The weight fraction is " << nanotube_geo.weight_fraction << endl;
        
		istr_cnt_vol >> nanotube_geo.accum_mode;
		if(nanotube_geo.accum_mode<0&&nanotube_geo.accum_mode>2){ hout <<"Error: the mode of accumulation should be between 0 and 2." << endl; return 0;  }
        
		istr_cnt_vol >> nanotube_geo.linear_density;		//Read the linear density of a nanotube
		if(nanotube_geo.linear_density<0){ hout << "Error: the linear density of a nanotube should be non-nagetive." << endl; return 0; }
		
		istr_cnt_vol >> geom_rve.density;	 //Read the density of RVE. Here we ignore the volume of nantubes, so the density of RVE actually approximates to the density of matix
		if(geom_rve.density<0){ hout << "Error: the density of RVE should be non-nagetive." << endl; return 0; }
		if(nanotube_geo.linear_density>=nanotube_geo.matrix_density){ hout << "Error: the density of matrix or the linear density of a nanotube is wrong." << endl; return 0; }
        
		//The real weight of nanotubes
		nanotube_geo.real_weight = nanotube_geo.weight_fraction*geom_rve.volume*geom_rve.density;
	}
	else { hout << "Error: the criterian of generation is neither 'vol' nor 'wt'." << endl; return 0; }
    


	return 1;
}
//---------------------------------------------------------------------------
//Reading weighting function information
int Input::Read_cluster_geo_parameters(struct Cluster_Geo &cluster_geo, ifstream &infile)
{
	if(cluster_geo.mark)
	{
		cout << "Attention: \"" << cluster_geo.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << cluster_geo.keywords << "\" has been input!" << endl;
		return 0;
	}
	else cluster_geo.mark = true;

	istringstream istr_clust_para(Get_Line(infile));
	istr_clust_para >> cluster_geo.vol_fra_criterion;
	if(cluster_geo.vol_fra_criterion>1||cluster_geo.vol_fra_criterion<0){ hout << "Error: the volume fraction of clusters must be between 0 and 1." << endl; return 0; }
	if(cluster_geo.vol_fra_criterion!=0)  //0 means that it doesn't need to generate nanotube clusters
	{
		istr_clust_para >> cluster_geo.amin >> cluster_geo.amax;
		if(cluster_geo.amin<0||cluster_geo.amax<0||cluster_geo.amin>cluster_geo.amax){ hout << "Error: the length of a long axis must be non-negative and min must be smaller than max."<< endl; return 0; }
		
		istr_clust_para >> cluster_geo.bmin >> cluster_geo.cmin;
		if(cluster_geo.bmin<0||cluster_geo.bmin>cluster_geo.amin||
		   cluster_geo.cmin<0||cluster_geo.cmin>cluster_geo.bmin||
           cluster_geo.cmin>cluster_geo.amin)
		{ 
			hout << "Error: the lengths of middle and short axis must be non-negative;" << endl;
			hout << "or error: min of middle and short axes must be smaller than min of long axis;" << endl;
			hout << "or error: min of short axis must be smaller than min of middle axis." << endl;
			return 0; 
		}
		
		istr_clust_para >> cluster_geo.growth_probability;
		if(cluster_geo.growth_probability<0||cluster_geo.growth_probability>1){ hout << "Error: the growth probability of nanotubes in clusters must be between 0 and 1." << endl; return 0; }
		
		istr_clust_para >> cluster_geo.volf_clust;
		if(cluster_geo.volf_clust<0||cluster_geo.volf_clust>1) { hout << "Error: the volume fraction of nanotubes in clusters must be between 0 and 1." << endl; return 0; }

		istr_clust_para >> cluster_geo.print_key;
		if(cluster_geo.print_key!=0&&cluster_geo.print_key!=1&&cluster_geo.print_key!=2) { hout << "Error: the print_key of cluster_geo is not 0, 1 and 2." << endl; return 0; }

		//clear ellipsoid vector
		cluster_geo.ellips.clear();
	}

	return 1;
}
//---------------------------------------------------------------------------
//Reading material properties of elements
int Input::Read_cutoff_distances(struct Cutoff_dist &cutoff_dist, ifstream &infile)
{
	if(cutoff_dist.mark)
	{
		cout << "Attention: \"" << cutoff_dist.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << cutoff_dist.keywords << "\" has been input!" << endl;
		return 0;
	}
	else cutoff_dist.mark = true;

	istringstream istr0(Get_Line(infile));
	istr0 >> cutoff_dist.van_der_Waals_dist;

	istringstream istr1(Get_Line(infile));
	istr1 >> cutoff_dist.tunneling_dist;

	return 1;
}
//---------------------------------------------------------------------------
//Reading electrical properties of materials
int Input::Read_electrical_parameters(struct Electric_para &electric_para, ifstream &infile)
{
    if(electric_para.mark)
    {
        cout << "Attention: \"" << electric_para.keywords << "\" has been input!" << endl;
        hout << "Attention: \"" << electric_para.keywords << "\" has been input!" << endl;
        return 0;
    }
    else electric_para.mark = true;
    
    istringstream istr0(Get_Line(infile));
    istr0 >> electric_para.applied_voltage;
    
    //Carbon fiber resistivity
    istringstream istr1(Get_Line(infile));
    istr1 >> electric_para.resistivity_CF;
    
    //GNP resistivities
    istringstream istr2(Get_Line(infile));
    istr2 >> electric_para.sheet_resitance_GNP >> electric_para.resistivity_GNP_t >> electric_para.resistivity_GNP_surf;
    //hout << electric_para.resistivity_GNP_t <<" ,"<<electric_para.resistivity_GNP_surf<<endl;
    
    //Polymer matrix resistivity
    istringstream istr_rho(Get_Line(infile));
    istr_rho >> electric_para.resistivity_matrix;
    //hout << "electric_para.resistivity_matrix = " << electric_para.resistivity_matrix << endl;
    //Electrical constants: electron charge (C), permitivity of vaccum (F/m), CNT work function (V), dielectric constant of polymer
    //istringstream istr3(Get_Line(infile));
    //istr3 >> electric_para.e_charge >> electric_para.e0_vacuum >> electric_para.CNT_work_function >> electric_para.K_polymer;
    
    //Electrical constants:Planck’s constant (m2kg/s), electron charge (C), electron mass (Kg), height barrier (eV)
    istringstream istr4(Get_Line(infile));
    istr4 >> electric_para.h_plank >> electric_para.e_charge >> electric_para.e_mass >> electric_para.lambda_barrier;
    //hout << electric_para.h_plank <<", "<< electric_para.e_charge <<", "<< electric_para.e_mass <<", "<< electric_para.lambda_barrier << endl;
    
    return 1;
}
//---------------------------------------------------------------------------
//Reading crack information
int Input::Read_tecplot_flags(struct Tecplot_flags &tec360_flags, ifstream &infile)
{
	if(tec360_flags.mark)
	{
		cout << "Attention: \"" << tec360_flags.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << tec360_flags.keywords << "\" has been input!" << endl;
		return 0;
	}
	else tec360_flags.mark = true;
    
    //Flag to export generated CNTs:
    // 1: wires
    // 2: singlezone meshes
    // 3: multizone meshes
    //Do not export if flag set to 0
	istringstream istr0(Get_Line(infile));
	istr0 >> tec360_flags.generated_cnts;
    
    //Flag to export generated GNPs if set to 1
    //Do not export if flag set to 0
    istringstream istr1(Get_Line(infile));
    istr1 >> tec360_flags.generated_gnps;
    
    //Flag to export clusters as obtained from the Hoshen-Kopelman algorithm:
    // 1: wires (3D lines)
    // 2: meshes
    // 3: wires (3D lines) & meshes
    //Do not export if flag set to 0
    istringstream istr2(Get_Line(infile));
    istr2 >> tec360_flags.clusters;

    //Flag to export percolated clusters:
    // 1: wires (3D lines)
    // 2: meshes
    // 3: wires (3D lines) & meshes
    //Do not export if flag set to 0
    istringstream istr3(Get_Line(infile));
    istr3 >> tec360_flags.percolated_clusters;
    
    //Flag to export the backbone:
    // 1: wires (3D lines)
    // 2: meshes
    // 3: wires (3D lines) & meshes
    //Do not export if flag set to 0
    istringstream istr4(Get_Line(infile));
    istr4 >> tec360_flags.backbone;
    
    //Flag to export triangulations if set to 1
    //Do not export if flag set to 0
    istringstream istr5(Get_Line(infile));
    istr5 >> tec360_flags.triangulations;
    
	return 1;
}
//---------------------------------------------------------------------------------------
//-------- AMC
//---------------------------------------------------------------------------
//Reading the geometric parameters of GNPs
int Input::Read_gnp_geo_parameters(struct GNP_Geo &gnp_geo, ifstream &infile)
{
    if(gnp_geo.mark)
    {
        cout << "Attention: \"" << gnp_geo.keywords << "\" has been input!" << endl;
        hout << "Attention: \"" << gnp_geo.keywords << "\" has been input!" << endl;
        return 0;
    }
    else gnp_geo.mark = true;

    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Define if the generation of the CNTs on the GNP surface should be either parallel or independent
    istringstream istr_gnp_grow_type(Get_Line(infile));
    istr_gnp_grow_type >> gnp_geo.growth_type;
    if (gnp_geo.growth_type != "parallel" && gnp_geo.growth_type != "independent") {
        hout << "Error: the growth type of the CNTs on the GNP surface should be either parallel or independent." << endl;
        hout << "Input was: " << gnp_geo.growth_type << endl;
        return 0;
    }
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Define the GNP orientation type (random or specific) in a RVE
    istringstream istr_initial_orientation(Get_Line(infile));
    istr_initial_orientation >> gnp_geo.orient_distrib_type;
    if(gnp_geo.orient_distrib_type!="random"&&gnp_geo.orient_distrib_type!="specific"){ hout << "Error: the GNP orientation type must be either random or specific." << endl;	return 0; }
    if(gnp_geo.orient_distrib_type=="specific")
    {
        istr_initial_orientation  >> gnp_geo.ini_sita >> gnp_geo.ini_pha;
        if(gnp_geo.ini_sita<0||gnp_geo.ini_sita>PI||gnp_geo.ini_pha<0||gnp_geo.ini_pha>=2*PI)
        {
            hout << "Error: the specified angle is not in the acceptable range of (0, 2PI)." << endl;
            return 0;
        }
    }
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Define the step length (unit: micromether) of nanotube growth
    istringstream istr_step_len(Get_Line(infile));
    istr_step_len >> gnp_geo.discr_step_length;
    if(gnp_geo.discr_step_length<=0||
       gnp_geo.discr_step_length>=0.25*geom_rve.len_x||
       gnp_geo.discr_step_length>=0.25*geom_rve.wid_y||
       gnp_geo.discr_step_length>=0.25*geom_rve.hei_z)
    { hout << "Error: the step length must be positive and 0.25 times lesser than the dimension of the RVE box." << endl;	return 0; }
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Define the distribution type (uniform or normal) of the length and width (unit: micrometer) of GNPs and the range (min, max) of GNPs in a RVE
    istringstream istr_gnp_len(Get_Line(infile));
    istr_gnp_len >> gnp_geo.size_distrib_type;
    if(gnp_geo.size_distrib_type!="uniform"&&gnp_geo.size_distrib_type!="normal"){ hout << "Error: the distribution of the length should be either normal or uniform." << endl;	return 0; }
    istr_gnp_len >> gnp_geo.len_min >> gnp_geo.len_max;
    if(gnp_geo.len_min<0||gnp_geo.len_max<0||gnp_geo.len_max<gnp_geo.len_min){ hout << "Error: the length must be non-negative and min must be smaller than max." << endl; return 0; }
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Define the distribution type (uniform or normal) of the thickness (unit: micrometer) of GNPs and the thickness range (min, max) of GNPs in a RVE
    istringstream istr_gnp_thick(Get_Line(infile));
    istr_gnp_thick >> gnp_geo.thick_distrib_type;
    if(gnp_geo.thick_distrib_type!="uniform"&&gnp_geo.thick_distrib_type!="normal"){ hout << "Error: the distribution of the radius should be either normal or uniform." << endl;	return 0; }
    istr_gnp_thick >> gnp_geo.t_min >> gnp_geo.t_max;
    if(gnp_geo.t_min<0||gnp_geo.t_max<0||gnp_geo.t_max<gnp_geo.t_min)
    { hout << "Error: the GNP thickness must be non-negative, min must be smaller than max." << endl; return 0; }
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Define the mass ratio of GNP/nanotubes in the RVE
    istringstream istr_mass_ratio(Get_Line(infile));
    istr_mass_ratio >> gnp_geo.mass_ratio;
    if(gnp_geo.mass_ratio <= 0 && geom_rve.particle_type != "Hybrid_particles")
    { hout << "Error: the GNP/CNT mass ratio must be a value greater than zero." << endl; return 0; }
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Define the density of the CNTs (unit: gm/cm3)
    istringstream istr_density(Get_Line(infile));
    istr_density >> gnp_geo.density;
    if (gnp_geo.density <= 0) {
        hout << "Error: the density has to be a value greater than zero. Input value:" << nanotube_geo.density;
        hout << "." << endl;
        return 0;
    }
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Define the volume or weight fraction of GNPs in the RVE
    istringstream istr_gnp_vol(Get_Line(infile));
    istr_gnp_vol >> gnp_geo.criterion;
    if(gnp_geo.criterion=="vol")
    {
        //Read the input value to keep the flow of the reading file unaltered
        istr_gnp_vol >> gnp_geo.volume_fraction;
        //Check the particle type
        //If mixed CNTs+GNPs or hybrid particles are used, then the fraction of GNPs needs to be calculated from the CNT fraction
        //Then the CNT fraction needs to be adjusted
        if (geom_rve.particle_type == "GNP_CNT_mix" ) {
            
            //Calculate the GNP volume fraction
            gnp_geo.volume_fraction = nanotube_geo.density*nanotube_geo.volume_fraction/(gnp_geo.mass_ratio*gnp_geo.density + nanotube_geo.density);
            
            //Adjust the CNT volume fraction
            nanotube_geo.volume_fraction = nanotube_geo.volume_fraction - gnp_geo.volume_fraction;
            //Adjust the total volume of the nanotube network
            nanotube_geo.real_volume = nanotube_geo.volume_fraction*geom_rve.volume;
            
            hout << "    Given the CNT/GNP mass ratio of " << gnp_geo.mass_ratio << ", the CNT volume fraction was adjusted to " << nanotube_geo.volume_fraction << endl;
            hout << "    The GNP volume fraction is " << gnp_geo.volume_fraction <<endl;
        }
        else if(geom_rve.particle_type == "Hybrid_particles") {
            
            //If the mass ratio is positive proceed as with mixed particles
            if (gnp_geo.mass_ratio > Zero) {
                
                //Calculate the GNP volume fraction
                gnp_geo.volume_fraction = nanotube_geo.density*nanotube_geo.volume_fraction/(gnp_geo.mass_ratio*gnp_geo.density + nanotube_geo.density);
                
                //Adjust the CNT volume fraction
                //nanotube_geo.volume_fraction = nanotube_geo.volume_fraction - gnp_geo.volume_fraction;
                
                //Adjust the total volume of the nanotube network
                //nanotube_geo.real_volume = nanotube_geo.volume_fraction*geom_rve.volume;
                
                hout << "    Given the CNT/GNP mass ratio of " << gnp_geo.mass_ratio << ", the CNT volume fraction is " << (nanotube_geo.volume_fraction - gnp_geo.volume_fraction) << endl;
                hout << "    The GNP volume fraction is " << gnp_geo.volume_fraction <<endl;
            }
            //If the mass ratio is negative, then the absolute value is used as CNT density and only calculate the fraction of GNPs
            else {
                //Calculate the mass ratio from the CNT density
                //Use the smallest CNT geometry and largest GNP thickness to minimize mas ratio M
                //A minimum M will result in a maximum vlolume of GNP
                double M = 2*nanotube_geo.density*PI*nanotube_geo.rad_min*nanotube_geo.rad_min*nanotube_geo.len_min*abs(gnp_geo.mass_ratio);
                M = M/(gnp_geo.density*gnp_geo.t_max);
                
                hout << "Calculated mass ratio is " << M << endl;
                
                //Calculate the GNP volume fraction
                gnp_geo.volume_fraction = nanotube_geo.volume_fraction/M;
            }
            
        }
        else {
            
            //Check that the volume fraction is between 0 and 1
            if(gnp_geo.volume_fraction>1||gnp_geo.volume_fraction<0) {
                hout << "Error: the volume fraction must be between 0 and 1." << endl;
                return 0;
            }
            
        }
        
        //The real volume of GNPs
        gnp_geo.real_volume = gnp_geo.volume_fraction*geom_rve.volume;
    }
    else if(gnp_geo.criterion=="wt")
    {
        //Read the input value to keep the flow of the reading file unaltered
        istr_gnp_vol >> gnp_geo.volume_fraction;
        //Check the particle type
        //If mixed CNTs+GNPs or hybrid particles are used, then the fraction of GNPs needs to be calculated from the CNT fraction
        //Then the CNT fraction needs to be adjusted
        if (geom_rve.particle_type == "GNP_CNT_mix" || geom_rve.particle_type == "Hybrid_particles") {
            //Calculate the GNP weight fraction
            gnp_geo.weight_fraction = nanotube_geo.weight_fraction/(1+gnp_geo.mass_ratio);
            
            //Adjust the CNT weight fraction
            nanotube_geo.weight_fraction = nanotube_geo.weight_fraction - gnp_geo.weight_fraction;
            //Adjust the real weight of nanotubes
            nanotube_geo.real_weight = nanotube_geo.weight_fraction*geom_rve.volume*geom_rve.density;
            
            hout << "    Given the CNT/GNP mass ratio of " << gnp_geo.mass_ratio << ", the CNT weight fraction was adjusted to " << nanotube_geo.weight_fraction << endl;
            hout << "    The GNP weight fraction is " << gnp_geo.weight_fraction <<endl;
        } else {
            if(gnp_geo.weight_fraction>1||gnp_geo.weight_fraction<0){ hout << "Error: the volume fraction must be between 0 and 1." << endl; return 0; }
            hout << "    The weight fraction is " << gnp_geo.weight_fraction << endl;
        }
        
        //The real weight of GNPs
        gnp_geo.real_weight = gnp_geo.weight_fraction*geom_rve.volume*geom_rve.density;
    }
    else {
        hout << "Error: the criterion of generation is neither 'vol' nor 'wt'. Input was:"<< gnp_geo.criterion << endl;
        return 0;
    }
    
    
    return 1;
}
//---------------------------------------------------------------------------
//Read the input data in a whole line (to skip over the comment line starting with a '%')
string Input::Get_Line(ifstream &infile)const
{
	string s;
	//Read the input data in a whole line
	getline(infile,s);
	//to skip over the comment line starting with a '%'
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===========================================================================
