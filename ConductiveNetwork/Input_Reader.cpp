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
		else if(str_temp=="Electrical_Parameters")	{ if(Read_electrical_paramters(electric_para, infile)==0) return 0; }
		else 
		{ 
			cout << "Error: the keywords \"" << str_temp << "\" is not defined!" << endl; 
			hout << "Error: the keywords \"" << str_temp << "\" is not defined!" << endl; 
			return 0; 
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

	//Initialize the geometric parameters of the RVE
	geom_rve.keywords = "RVE_Geometry";
	geom_rve.mark = false;
	geom_rve.origin.x = 0.0;
	geom_rve.origin.y = 0.0;
	geom_rve.origin.z = 0.0;
	geom_rve.origin.flag = 0;
	geom_rve.len_x = 1.0;
	geom_rve.wid_y = 1.0;
	geom_rve.hei_z = 1.0;
	geom_rve.volume = geom_rve.len_x*geom_rve.wid_y*geom_rve.hei_z;
	geom_rve.Nx = 1;
	geom_rve.Ny = 1;
	geom_rve.Nz = 1;
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
	cluster_geo.vol_fra_criterion = 0;
	cluster_geo.amin = 0.0;
	cluster_geo.amax = 0.0;
	cluster_geo.bmin = 0.0;
	cluster_geo.cmin = 0.0;
	cluster_geo.growth_probability = 0.0;
	cluster_geo.wt_fra_cluster = 0.0;
	cluster_geo.cnt_real_weight = 0.0;
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

	istringstream istr2(Get_Line(infile));		
	istr2 >> simu_para.create_read_network;		//Read a signal to show if create a new network or read a previouse network from a file

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

	istringstream istr(Get_Line(infile));
	istr >> geom_rve.origin.x >> geom_rve.origin.y >> geom_rve.origin.z;
	istr >> geom_rve.len_x >> geom_rve.wid_y >> geom_rve.hei_z;
	if(geom_rve.len_x<0||geom_rve.wid_y<0||geom_rve.hei_z<0)
	{
		cout << "Error: the sizes of RVE should be positive!" << endl;
		hout << "Error: the sizes of RVE should be positive!" << endl;
		return 0;
	}
	geom_rve.volume = geom_rve.len_x*geom_rve.wid_y*geom_rve.hei_z;

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

	istringstream istr(Get_Line(infile));

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

	istringstream istr(Get_Line(infile));

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
	cutoff_dist.van_der_Waals_dist;

	istringstream istr1(Get_Line(infile));
	cutoff_dist.tunneling_dist;

	return 1;
}
//---------------------------------------------------------------------------
//Reading crack information
int Input::Read_electrical_paramters(struct Electric_para &electric_para, ifstream &infile)
{
	if(electric_para.mark)
	{
		cout << "Attention: \"" << electric_para.keywords << "\" has been input!" << endl;
		hout << "Attention: \"" << electric_para.keywords << "\" has been input!" << endl;
		return 0;
	}
	else electric_para.mark = true;

	istringstream istr0(Get_Line(infile));
	electric_para.applied_voltage;		

	istringstream istr1(Get_Line(infile));
	electric_para.resistivity_CF;

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
