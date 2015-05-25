//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Input_Reader.h
//OBJECTIVE:	Reading the input data
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa;  angel.mora@kaust.edu.sa
//====================================================================================

#ifndef INPUTREADER_H
#define INPUTREADER_H

#include<iomanip>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include "time.h"
#include "Hns.h"
using namespace hns;

#include "Geometry_3D.h"

//---------------------------------------------------------------------------
//Name of application case
struct App_name{
			string keywords;
			bool mark;
			string str;
		};
//Name of simulation
struct Simu_para{
			string keywords;
			bool mark;
			string str;
			int sample_num;
		};
//The geometry of the RVE
struct Geom_RVE{
			string keywords;
			bool mark;
			Point_3D origin;
			double len_x, wid_y, hei_z;
			double volume;
			int Nx, Ny, Nz;  //Define 'Nx Ny Nz' which are the number of segments in each direction by which the RVE is going to be divided (for looking for penetrating nanotubes)
			//Define the size range of the observation window and descrement by every step in x, y and z directions
			double win_max_x, win_max_y, win_max_z;
			double win_min_x, win_min_y, win_min_z;
			double win_delt_x, win_delt_y, win_delt_z;
		};
//The nanotube parameters in a network
struct Nanotube_para{
			string keywords;
			bool mark;
			string criterion;					//Define the volume or weight fraction of nanotubes in the RVE: vol, wt, nwt
			string dir_dist_type;			//Define the initial growth direction type (random or specific) in a RVE
			string len_dist_type;			//Define the distribution type (uniform or normal) of the length (unit: micromether) of nanotubes
			string rad_dist_type;			//Define the distribution type (uniform or normal) of the radius (unit: micromether) of nanotubes
			double step_length;				//Define the step length (unit: micromether) of nanotube growth
			double ini_sita, ini_pha;			//Define initial direction for 'specific' type in the spherical coordinates
			double angle_max;				//Define the angle 'omega' for the normal distribution range [-omega, omega] of the growth direction
			double len_min, len_max;		//Define the length range (min, max) of nanotubes
			double rad_min, rad_max;		//Define the radius range (min,max) of nanotubes
			double volume_fraction;		//Define the volume fraction of nanotubes
			double real_volume;				//Define the real volume of nanotubes
			double weight_fraction;			//Define the weight fraction of nanotubes
			double real_weight;				//Define the real weight of nanotubes
			double linear_density;			//Define the linear density of nanotubes
		};
//---------------------------------------------------------------------------
//记录纳米管团簇椭球参数信息
struct Clust_Geo
{
	double wt_fra_cluster;				//纳米管团簇中纳米管的重量分数
	double vol_fra_criterion;			//纳米管团簇椭球体积分数（界限）
	double amin;								//纳米管团簇椭球长轴取值范围最小值
	double amax;								//纳米管团簇椭球长轴取值范围最大值
	double bmin;								//纳米管团簇椭球中轴取值范围最小值
	double cmin;								//纳米管团簇椭球短轴取值范围最小值
	double growth_probability;		//纳米管在团簇椭球中生长的概率
	double real_volume_fraction;	//纳米管团簇椭球的实际体积
	double cnt_real_weight;			//纳米管团簇椭球中所包含的纳米管实际比重
	double cnt_real_volume;			//纳米管团簇椭球中所包含的纳米管实际体积
};
//The cutoff distances
struct Cutoff_dist{
			string keywords;
			bool mark;
			double van_der_Waals_dist;
			double tunneling_dist;
		};



//---------------------------------------------------------------------------
class Input 
{
	public:
		//Data members
		struct App_name app_name;
		struct Simu_name simu_name;
		struct Stif_loc stif_loc;
		struct Stif_nonloc stif_nonloc;
		struct Peri_para peri_para;
		struct Geom_RVE geom_rve;
		struct Grid_size grid_size, nonloc_gsize;
		struct Weight_func weight_func;
		struct Ele_prop ele_prop;
		struct Crack cracks;
		struct Damage damages;
		struct Model_Discret mod_disc;
		struct Iterative iter;
		struct RW_mod rw_mod;
		struct Gauss_Point gauss, nonloc_gau;
		struct Load load;
		struct Displace displace;

		//Constructor
		Input(){};  

		//Member functions
		int Data_Initialization();		//Initialize data
		int Read_Infile(ifstream &infile);		//Read data
		string Get_Line(ifstream &infile)const;									//读入信息一行，跳过注释行（以%开头）

private:
		//Member functions
		int Read_application_name(struct App_name &app_name, ifstream &infile);
		int Read_simulation_name(struct Simu_name &simu_name, ifstream &infile);
		int Read_local_stiffness(struct Stif_loc &stif_loc, ifstream &infile);
		int Read_local_stiffness_2D(struct Stif_loc &stif_loc, ifstream &infile);
		int Read_nonlocal_stiffness(struct Stif_nonloc &stif_nonloc, ifstream &infile);
		int Read_nonlocal_stiffness_2D(struct Stif_nonloc &stif_nonloc, ifstream &infile);
		int Read_peridynamic_parameters(struct Peri_para &peri_para, ifstream &infile);
		int Read_rve_geometry(struct Geom_RVE &geom_rve, ifstream &infile);
		int Read_grid_size(struct Grid_size &grid_size, ifstream &infile);
		int Read_weight_function(struct Weight_func &weight_func, ifstream &infile);
		int Read_element_properties(struct Ele_prop &ele_prop, ifstream &infile);
		int Read_crack(struct Crack &cracks, ifstream &infile);
		int Read_damage(struct Damage &damages, ifstream &infile);
		int Read_model_discret(struct Model_Discret &mod_disc, ifstream &infile);
		int Read_iterative(struct Iterative &iter, ifstream &infile);
		int Read_rw_mod(struct RW_mod &rw_mod, ifstream &infile);
		int Read_gauss(struct Gauss_Point &gauss, ifstream &infile);
		int Read_load(struct Load &load, ifstream &infile);
		int Read_load_2D(struct Load &load, ifstream &infile);
		int Read_displacement(struct Displace &displace, ifstream &infile);
		int Read_displacement_2D(struct Displace &displace, ifstream &infile);
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
