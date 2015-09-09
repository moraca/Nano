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
		int Generate_geometric_networks(const struct Geom_RVE &geom_rve, struct Cluster_Geo &cluster_geo)const;

	private:
		//Data Member

		//Generate a number of ellipsoids
		int Get_ellip_clusters(const struct Geom_RVE &cell, struct Cluster_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle)const;
		//生成特定位置的圆球团簇序列(避免随机情况时大数量椭球无法生成)
//		int Get_specific_sphere_clusters(const struct RVE_Geo &cell, struct Clust_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle, vector<struct elliparam> &ellips, const int &export_mod)const;   
};
//-------------------------------------------------------
#endif
//===========================================================================