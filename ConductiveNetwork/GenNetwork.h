//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	GenNetwork.h
//OBJECTIVE:	To generate networks with overlapping
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef GENNETWORK_H
#define GENNETWORK_H

#include<iomanip>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include "Input_Reader.h"
#include "Geometry_3D.h"
#include "MathMatrix.h"
#include "Fem_3D.h"
#include "Gauss.h"
#include "Hns.h"
using namespace hns;

//-------------------------------------------------------
class GenNetwork
{
	public:
		//Data Member
		
		//Constructor
		GenNetwork(){};

		//Member Functions
		int Generate_geometric_networks(const struct Geom_RVE &geom_rve, struct Cluster_Geo &clust_geo, struct Nanotube_Geo &nanotube_geo)const;

	private:
		//Data Member

		//Generate a number of ellipsoids
		int Get_ellip_clusters(const struct cuboid &cub, struct Cluster_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle)const;
		//Generate a number of sperical clusters in regular arrangement
		int Get_spherical_clusters_regular_arrangement(const struct cuboid &cub, struct Cluster_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle)const;
		//Print the ellipsoid surfaces by grids
		void Export_cluster_ellipsoids_mesh(const struct cuboid &cub, const vector<struct elliparam> &ellips)const;
		//Export the data of ellipsoid surfaces
		void Export_cluster_ellipsoids_data(const vector<struct elliparam> &ellips, const double &ellip_ratio)const;
		//Randomly generate a seed (original point) of a CNT in the RVE
		int Get_seed_point(const struct cuboid &cub, int &seed, Point_3D &point)const;
		//Generate a random value through a probability distribution function
		int Get_random_value(const string &dist_type, const double &min, const double &max, int &seed, double &value)const;
		//Randomly generate a direction in the spherical coordinates as the original direction of CNT segments
		int Get_uniform_direction(const struct Nanotube_Geo &nanotube_geo, int &seed_sita, int &seed_pha, double &cnt_sita, double &cnt_pha)const;
		//Transform angles into matrix
		MathMatrix Get_transformation_matrix(const double &sita, const double &pha)const;
};
//-------------------------------------------------------
#endif
//===========================================================================