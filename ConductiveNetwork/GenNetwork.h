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
#include <random> 
#include "Input_Reader.h"
#include "Geometry_3D.h"
#include "MathMatrix.h"
#include "Fem_3D.h"
#include "Hns.h"
#include "Tecplot_Export.h"

using namespace hns;

//-------------------------------------------------------
class GenNetwork
{
public:
    //Data Member
    
    //Constructor
    GenNetwork(){};
    
    //Member Functions
    int Generate_nanofiller_network(const struct Geom_RVE &geom_rve, const struct Cluster_Geo &clust_geo, const struct Nanotube_Geo &nanotube_geo, const struct GNP_Geo &gnp_geo, const struct Cutoff_dist &cutoffs, const struct Tecplot_flags &tec360_flags, vector<Point_3D> &cpoints, vector<Point_3D> &gnps_points, vector<GCH> &hybrid_particles, vector<double> &cnts_radius, vector<vector<long int> > &cstructures, vector<vector<long int> > &gstructures)const;
    //Generate the nodes and tetrahedron elements of nanotubes (No const following this function because a sum operation on two Point_3D points inside)
    int Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius);
    //Generate the nodes and tetrahedron elements of nanotubes (No const following this function because a sum operation on two Point_3D points inside). This function uses a 1D point vector and a 2D structure vector that references the point vector
    int Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<Point_3D> &cnts_points, const vector<double> &cnts_radius, const vector<vector<long int> > &structure);    
    
//private:
    //Data Member
    
    //Member Functions
    
    //Generate a network defined by points and connections
    int Generate_cnt_network_threads_mt(const struct Geom_RVE &geom_rve, const struct Cluster_Geo &clust_geo, const struct Nanotube_Geo &nanotube_geo, const struct Cutoff_dist &cutoffs, vector<vector<Point_3D> > &cnts_points,  vector<double> &cnts_radius)const;
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    int Get_random_value_mt(const string &dist_type, mt19937 &engine, uniform_real_distribution<double> &dist, const double &min, const double &max, double &value)const;
    int Get_seed_point_mt(const struct cuboid &cub, Point_3D &point, mt19937 &engine_x, mt19937 &engine_y, mt19937 &engine_z, uniform_real_distribution<double> &dist)const;
    int Get_initial_direction_mt(const string &dir_distrib_type, const double &ini_sita, const double &ini_pha, mt19937 &engine_inital_direction, uniform_real_distribution<double> &dist_initial, MathMatrix &rotation)const;
    int Get_normal_direction_mt(const double &omega, double &cnt_sita, double &cnt_pha, mt19937 &engine_sita, mt19937 &engine_pha, uniform_real_distribution<double> &dist)const;
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    int Generate_gnp_network_mt(const struct GNP_Geo &gnp_geo, const struct Geom_RVE &geom_rve, const struct Cutoff_dist &cutoffs, vector<vector<Point_3D> > &gnps_points, vector<GCH> &hybrid_praticles, double &carbon_vol, double &carbon_weight)const;
    int Generate_cnt_network_threads_over_gnps_mt(const struct GNP_Geo &gnp_geo, const struct Geom_RVE &geom_rve, const struct Nanotube_Geo &nanotube_geo, const struct Cutoff_dist &cutoffs, vector<vector<Point_3D> > &cnts_points, vector<vector<Point_3D> > &gnps_points, vector<GCH> &hybrid_praticles, vector<double> &cnts_radius)const;
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // HYBRID PARTICLE
    int Generate_gnp(const struct GNP_Geo &gnp_geo, cuboid &gnp, mt19937 &engine_l, mt19937 &engine_t, uniform_real_distribution<double> &dist)const;
    int Handle_gnp_penetrations(const struct cuboid &gvcub, const vector<GCH> &hybrid_particles, GCH &hybrid, Point_3D &gnp_poi)const;
    int Check_penetration_and_move_gnp(const struct cuboid &gvcub, const vector<GCH> &hybrid_particles, GCH &hybrid, Point_3D &gnp_poi)const;
    int Calculate_new_gnp_location(const vector<GCH> &hybrid_particles, const vector<double> &cutoff, const vector<double> &distances, const vector<int> &indices, Point_3D &new_location)const;
    int GNP_inside_composite(const struct cuboid &gvcub, const GCH &hybrid)const;
    int Calculate_number_of_CNTs_on_GNP_surface(const struct Nanotube_Geo &nanotube_geo, const struct GNP_Geo &gnp_geo, const cuboid &gnp, const double &cnt_rad, const double &cnt_length, int &ns)const;
    int Generate_cnt_seeds(const struct Nanotube_Geo &nanotube_geo, const cuboid &gnp, vector<Point_3D> &seeds, vector<double> &radii, const double &d_vdw, const int &n_cnts, const double &z_coord, mt19937 &engine_x, mt19937 &engine_y, mt19937 &engine_rand, uniform_real_distribution<double> &dist);
    int Get_seed_point_2d_mt(const struct cuboid &cub, const double &cnt_rad, Point_3D &point, mt19937 &engine_x, mt19937 &engine_y, uniform_real_distribution<double> &dist)const;
    int Check_seed_overlapping(const vector<Point_3D> &seeds, const double &cnt_rad, const double &d_vdw, Point_3D &point)const;
    int Grow_cnts_on_gnp_surface(const struct Geom_RVE &geom_rve, const struct cuboid &excub, const GCH &hybrid, const struct Nanotube_Geo &nanotube_geo, const struct Cutoff_dist &cutoffs, vector<vector<Point_3D> > &cnts_points, vector<vector<Point_3D> > &gnps_points, vector<vector<Point_3D> > &particle, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain, vector<vector<int> > &global_coord_hybrid, vector<vector<long int> > &sect_domain_hybrid, vector<double> &cnts_radius, const vector<int> &n_subregions, const double &cnt_rad, const double  &cnt_length, const double &overlap_max_cutoff, const int &penetrating_model_flag, const int &ns, int &point_overlap_count, int &point_overlap_count_unique, int &cnt_reject_count, mt19937 &engine_x, mt19937 &engine_y, mt19937 &engine_sita, mt19937 &engine_pha, uniform_real_distribution<double> &dist)const;
    int Grow_single_cnt_on_gnp(const struct Geom_RVE &geom_rve, const struct cuboid &excub, const struct Nanotube_Geo &nanotube_geo, const struct Cutoff_dist &cutoffs, vector<vector<Point_3D> > &cnts_points, vector<vector<Point_3D> > &gnps_points, vector<vector<Point_3D> > &particle, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain, const vector<vector<int> > &global_coord_hybrid, const vector<vector<long int> > &sect_domain_hybrid, vector<Point_3D> &new_cnt, vector<double> &cnts_radius, const vector<int> &n_subregions, const MathMatrix &seed_multiplier, const double &cnt_rad, const double  &cnt_length, const int &penetrating_model_flag, int &point_overlap_count, int &point_overlap_count_unique, int &cnt_reject_count, mt19937 &engine_sita, mt19937 &engine_pha, uniform_real_distribution<double> &dist)const;
    int Grow_CNTs_on_GNP_surface_parallel(const struct Geom_RVE &geom_rve, const struct cuboid &excub, const struct cuboid &gnp, const struct Nanotube_Geo &nanotube_geo, const struct Cutoff_dist &cutoffs, vector<vector<Point_3D> > &cnts_points, vector<vector<Point_3D> > &gnps_points, vector<vector<Point_3D> > &particle, vector<Point_3D> &gnp_discrete, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain, vector<vector<int> > &global_coord_hybrid, vector<vector<long int> > &sect_domain_hybrid, vector<double> &cnts_radius, const vector<int> &n_subregions, const MathMatrix &multiplier, const Point_3D &gnp_poi, const double &cnt_rad, const double  &cnt_length, const int &penetrating_model_flag, const int &ns, int &point_overlap_count, int &point_overlap_count_unique, int &cnt_reject_count, mt19937 &engine_x, mt19937 &engine_y, mt19937 &engine_sita, mt19937 &engine_pha, uniform_real_distribution<double> &dist)const;
    int Copy_CNT(const vector<Point_3D> &initial_cnt, vector<Point_3D> &new_cnt)const;
    int Add_CNTs_to_hybrid_structure(const struct Geom_RVE &geom_rve, const vector<vector<Point_3D> > &particle, vector<Point_3D> &new_cnt, vector<vector<int> > &global_coord_hybrid, vector<vector<long int> > &sect_domain_hybrid, const vector<int> &n_subregions, const double &overlap_max_cutoff)const;
    int Add_CNTs_to_global_structure(const  struct Geom_RVE &geom_rve, const vector<vector<Point_3D> > &particle, vector<vector<Point_3D> > &cnts_points, vector<vector<int> > &global_coordinates, vector<vector<long int> > &sectioned_domain, vector<int> &n_subregions, vector<double> &cnts_radius, double &cnt_rad, double &overlap_max_cutoff, int &penetrating_model_flag)const;
    int Calculate_generated_volume(const struct GNP_Geo &gnp_geo, const struct cuboid &gvcub, const GCH &hybrid, const vector<vector<Point_3D> > &particle, const double &step_vol_para, const double &step_wei_para, double &vol_sum_cnt, double &wt_sum_cnt, double &vol_sum_gnp, double &wt_sum_gnp, double &vol_sum_tot, double &wt_sum_tot)const;
    int Calculate_generated_gnp_volume(const struct GNP_Geo &gnp_geo, const struct cuboid &gvcub, const GCH &hybrid, double &vol_sum, double &wt_sum)const;
    int Is_close_to_boundaries(const struct cuboid &gvcub, const GCH &hybrid)const;
    int Approximate_gnp_volume(const struct cuboid &gvcub, const GCH &hybrid, double &gnp_vol)const;
    int Discretize_gnp(const GCH &hybrid, const double &step, vector<Point_3D> &gnp_discrete)const;
    int Check_penetration_in_gch(const struct Geom_RVE &geom_rve, const struct Nanotube_Geo &nanotube_geo, const vector<vector<Point_3D> > &cnts, const vector<vector<Point_3D> > &gnps_points, const vector<vector<Point_3D> > &particle, const vector<Point_3D> &gnp_discrete, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain, const vector<vector<int> > &global_coord_hybrid, const vector<vector<long int> > &sect_domain_hybrid, const vector<double> &radii, const vector<Point_3D> &cnt_new, const vector<int> &n_subregions, const double &cnt_rad, const double &d_vdw, int &point_overlap_count, int &point_overlap_count_unique, Point_3D &point)const;
    int Get_penetrating_points(const vector<vector<Point_3D> > &cnts, const vector<vector<Point_3D> > &gnps_points,  const vector<vector<int> > &global_coordinates, const vector<long int> &subregion_vec, const vector<double> &radii, const double &cnt_rad, const double &d_vdw, const int &hybrid_flag, Point_3D &point, vector<Point_3D> &affected_points, vector<double> &cutoffs_p, vector<double> &distances)const;
    int Is_below_cutoff(const Point_3D &point_overlap, const Point_3D &point, const double &cutoff_p, vector<Point_3D> &affected_points, vector<double> &cutoffs_p, vector<double> &distances)const;
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //This functions initializes the vector sub-regions
    void Initialize_subregions(const struct Geom_RVE &geom_rve, vector<int> &nsubregions, vector<vector<long int> > &sectioned_domain)const;
    //Check if the current CNT is penetrating another CNT, i.e. is the new point is overlapping other point
    int Check_penetration(const struct Geom_RVE &geom_rve, const struct Nanotube_Geo &nanotube_geo, const vector<vector<Point_3D> > &cnts, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain, const vector<double> &radii, const vector<Point_3D> &cnt_new, const vector<int> &n_subregions, const double &cnt_rad, const double &d_vdW, int &point_overlap_count, int &point_overlap_count_unique, Point_3D &point)const;
    int Get_subregion(const struct Geom_RVE &geom_rve, const vector<int> &n_subregions, const Point_3D &point)const;
    void Get_penetrating_points(const vector<vector<Point_3D> > &cnts, const vector<vector<int> > &global_coordinates, const vector<long int> &subregion_vec, const vector<double> &radii, const double &cnt_rad, const double &d_vdw, Point_3D &point, vector<Point_3D> &affected_points, vector<double> &cutoffs_p, vector<double> &distances)const;
    void Move_point(const struct Geom_RVE &geom_rve, const struct Nanotube_Geo &nanotube_geo, const vector<Point_3D> &cnt_new, Point_3D &point, vector<double> &cutoffs, vector<double> &distances, vector<Point_3D> &affected_points)const;
    int Check_points_in_same_position(vector<double> &cutoffs, vector<double> &distances, vector<Point_3D> &affected_points)const;
    void Overlapping_points_same_position(const struct Geom_RVE &geom_rve, const struct Nanotube_Geo &nanotube_geo, const vector<Point_3D> &cnt_new, Point_3D &point)const;
    void One_overlapping_point(const vector<double> &cutoffs, const vector<double> &distances, const vector<Point_3D> &affected_points, Point_3D &point)const;
    void Two_overlapping_points(const vector<double> &cutoffs, const vector<Point_3D> &affected_points, Point_3D &point)const;
    void Three_or_more_overlapping_points(const vector<double> &cutoffs, const vector<double> &distances, const vector<Point_3D> &affected_points, Point_3D &point)const;
    int Check_segment_orientation(const Point_3D &point, const vector<Point_3D> &cnt_new)const;
    double Segment_angle_discriminant(const Point_3D &first, const Point_3D &second, const Point_3D &third)const;
    void Add_to_overlapping_regions(const struct Geom_RVE &geom_rve, double overlap_max_cutoff, Point_3D point, long int global_num, const vector<int> &n_subregions, vector<vector<long int> > &sectioned_domain)const;
    int Calculate_t(int a, int b, int c, int sx, int sy)const;
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //Checking the angle between two segments in one nanotube (if less than PI/2, provide an alarm)
    int CNTs_quality_testing(const vector<vector<Point_3D> > &cnts_points)const;
    //Generate a number of ellipsoids
    int Get_ellip_clusters(const struct cuboid &cub, const struct Cluster_Geo &clust_geo)const;
    //Generate a number of sperical clusters in regular arrangement
    int Get_spherical_clusters_regular_arrangement(const struct cuboid &cub, struct Cluster_Geo &clust_geo)const;
    //Print the ellipsoid surfaces by grids
    void Export_cluster_ellipsoids_mesh(const struct cuboid &cub, const vector<struct elliparam> &ellips)const;
    //Export the data of ellipsoid surfaces
    void Export_cluster_ellipsoids_data(const vector<struct elliparam> &ellips, const double &ellip_ratio)const;
    //Transform angles into matrix
    MathMatrix Get_transformation_matrix(const double &sita, const double &pha)const;
    //To calculate the coordinates of the new CNT point (transformation of coordinates)
    Point_3D Get_new_point(MathMatrix &Matrix, const double &Rad)const;
    //To judge if a point is included in a RVE
    int Judge_RVE_including_point(const struct cuboid &cub, const Point_3D &point)const;
    //To judge if a point is included in a RVE
    int Judge_RVE_including_point(const struct Geom_RVE &geom_rve, const Point_3D &point)const;
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
    //Transform the 2D cnts_points into 1D cpoints and 2D cstructuers
    int Transform_points(const string &type, vector<vector<Point_3D> > &cnts_points, vector<Point_3D> &cpoints, vector<vector<long int> > &cstructures)const;
};
//-------------------------------------------------------
#endif
//===========================================================================