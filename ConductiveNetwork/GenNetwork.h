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
		int Generate_geometric_networks(const struct Geom_RVE &geom_rve, const struct Cluster_Geo &clust_geo, const struct Nanotube_Geo &nanotube_geo, 
															vector<vector<Point_3D> > &cnts_points,  vector<double> &cnts_radius)const;
		//Checking the angle between two segments in one nanotube (if less than PI/2, provide an alarm)
		int CNTs_quality_testing(const vector<vector<Point_3D> > &cnts_points)const;
		//Generate the nodes and tetrahedron elements of nanotubes (No const following this function because a sum operation on two Point_3D points inside)
		int Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius);

	private:
		//Data Member

		//Generate a number of ellipsoids
		int Get_ellip_clusters(const struct cuboid &cub, const struct Cluster_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle)const;
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
		//Randomly generate a direction in the spherical coordinates, to have the positive Z-axis to be a central axis
		//Then, the radial angle, sita, obeys a normal distribution (sita \in fabs[(-omega,+omega)]) and the zonal angle, pha, obeys a uniform distribution (pha \in (0,2PI))
		int Get_normal_direction(const double &omega, int &seed_sita, int &seed_pha, double &cnt_sita, double &cnt_pha)const;
		//To calculate the coordinates of the new CNT point (transformation of coordinates)
		Point_3D Get_new_point(MathMatrix &Matrix, const double &Rad)const;
		//To judge if a point is included in a RVE
		int Judge_RVE_including_point(const struct cuboid &cub, const Point_3D &point)const;
		//Calculate all intersection points between the new segment and surfaces of RVE
		//(using a parametric equatio:  the parameter 0<t<1, and sort all intersection points from the smaller t to the greater t)  
		int Get_intersecting_point_RVE_surface(const struct cuboid &cub, const Point_3D &point0, const Point_3D &point1, vector<Point_3D> &ipoi_vec)const;
		//To calculate the effective portion (length) which falls into the given region (RVE)
		double Effective_length_given_region(const struct cuboid &cub, const Point_3D last_point, const Point_3D new_point)const;
		//Calculate the angles of a verctor in the spherical coordinate
		int Get_angles_vector_in_spherial_coordinates(const Point_3D &normal, double &sita, double &pha)const;
		//Calculate a group of equidistant points along the circumference which is on the plane defined by the center point of the circle and the normal vector
		int Get_points_circle_in_plane(const Point_3D &center, const double &sita, const double &pha, const double &radius, const int &num_sec, vector<Node> &nod_temp)const;
		//Calculate a group of projected points (which are on the plane with the center point of the circle and the normal vector) 
		//which are projected from a group of points on the previous circumference and projected along the direction of line_vec
		int Get_projected_points_in_plane(const Point_3D &center, const Point_3D &normal, const Point_3D &line, const int &num_sec, vector<Node> &nod_temp)const;

};
//-------------------------------------------------------
#endif
//===========================================================================