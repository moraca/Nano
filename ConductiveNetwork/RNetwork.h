//
//  RNetwork.h
//  Nanocode
//
//  Created by Angel Mora Cordova on 3/24/14.
//  Copyright (c) 2014 Angel Mora. All rights reserved.
//

#ifndef __Nanocode__RNetwork__
#define __Nanocode__RNetwork__

#include<iomanip>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include "Geometry_3D.h"
#include "MathMatrix.h"
#include "GeoNano.h"
#include "Gauss.h"
#include "Hns.h"
#include "Fem_3D.h"
#include <algorithm>
#include<omp.h>  //openmp
#include <string>



#define CHUNKSIZE 1

using namespace hns;

class RNetwork{
public:
    RVE_Geo box_geometry; //Geometry of the box inside the sample
    double lx, ly, lz, xmin, ymin, zmin, xmax, ymax, zmax;//variables for the geometry of the sample.
            //By making them global I reduce computations required to acces them
    int secx, secy, secz; //Number of regions on each direction
    long int points0; //Initial number of points
    double tunnel, cutoff;
    double magnitude;
    double resistivity;
    double window0, window_max, dwindow; //Size of the observation window and step size for increasing the window size
    vector<int> cnts_inside;
    vector<char> cnts_inside_flags;
    //Here the cnts that are in the boundaries will be stored
    vector<int> bbdyx1_cnt, bbdyx2_cnt, bbdyy1_cnt, bbdyy2_cnt, bbdyz1_cnt, bbdyz2_cnt; //Boundary vectors.
    vector<int> family; //This determines the family to which a cluster belongs to
    //0 for x-x; 1 for y-y; 2 for z-z; 3 for x-x and y-y; 4 for x-x and z-z; 5 for y-y and z-z; 6 for the three directions
    vector<double> clusters_lengths, dead_branches_lengths, clusters_fractions; //In these vector the lengts of all the CNTs in each cluster will be stored and their fractions
    vector<Point_3D> sphere_c; //List with all centers of enclosing spheres
    //The flag is used to store the number of the closest sphere
    vector<double> sphere_r, closest_distance; //List with all radii of enclosing spheres and distance to the closest sphere
    vector<vector<int> > clusters_cnt, directional_clusters, isolated, dead_branches; //vectors for percolating and non percolating CNTs
    vector<vector<long int> > sectioned_domain; //List with all points numbered from 0 to cnts_t-1
    vector<vector<long int> > contacts, contacts_cnt_point; //Contact vectors that store point numbers
    vector<vector<int> > contacts_cnt; //Contact vector that stores CNT numbers
    vector<vector<int> > boundary_flags; //This vector will help find points on the boundary
    vector<vector<short int> > boundary_flags_cnt, percolation_flags; //These vectors will help find the percolation clusters
    
    RNetwork(){};
    
    int Construct(ifstream &infile, string structure_type, const int &samples_count, struct CNT_Geo cnts_geo, vector<Point_3D> points_in, const struct RVE_Geo &cell_geo, struct Region_Geo cnt_regions,vector<double> cnts_radius, vector<vector<long int> > structure, vector<vector<int> > sectioned_domain_cnt);
    int Read_parameters(ifstream &infile, const struct RVE_Geo cell_geo, CNT_Geo cnts_geo);
    int Print_input_data_files(vector<vector<int> > sectioned_domain_cnt, vector<Point_3D> point_list, vector<double> cnt_radius, RVE_Geo cell_geo, CNT_Geo cnts_geo, Region_Geo cnt_regions);
    int Analysis(double window, const struct RVE_Geo cell_geo, struct CNT_Geo cnts_geo, struct Region_Geo cnt_regions, vector<Point_3D> points_in, vector<double> cnts_radius, vector<vector<int> > sectioned_domain_cnt, vector<vector<long int> > structure);
    int Region_of_interest(vector<Point_3D> &points_in, struct Region_Geo cnt_regions, vector<vector<long int> > &structure, vector<vector<int> > sectioned_domain_cnt, vector<double> &radii);
    void Get_boundary_cnts(int &sx0, int &sx1, int &sy0, int &sy1, int &sz0, int &sz1, struct Region_Geo cnt_regions, vector<vector<int> > sectioned_domain_cnt);
    void Update_cnts_inside_flags(vector<vector<int> > sectioned_domain_cnt, long int t);
    void Discard_repeated(vector<int> &vec);
    void Discard_repeated(vector<long int> &vec);
    int Locate_and_trim_boundary_cnts(vector<Point_3D> &points_in,vector<vector<long int> > &structure, vector<double> &radii);
    string Where_is(double x, double y, double z);
    int New_boundary_point(vector<Point_3D> &points_in, vector<vector<long int> > &structure, long int insidePoint, long int outsidePoint, int CNT, string currentLocation);
    int Projected_Point(double x, double y, double z, double x1, double y1, double z1, double &xp, double &yp, double &zp );
    void Projection(double a, double x0_0, double x0_1, double x1_0, double x1_1, double x2_0, double x2_1, double &x1_p, double &x2_p);
    void Trim_CNT(vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<int> &cnts_inside, vector<double> &radii, long int boundary, int CNT);
    void Add_to_boundary_vectors(Point_3D point3d, long int point);
    void Add_CNT_to_boundary(vector<int> &boundary, int CNT);
    void Remove_from_vector(int num, vector<int> &vec);
    void Remove_from_vector(long int num, vector<long int> &vec);
    int Assign_region(const struct RVE_Geo cell_geo, vector<Point_3D> points_in, vector<vector<long int> > structure);
    long int calculate_t(long int a, long int b, long int c, int sx, int sy);
    int Check_contacts(vector<Point_3D> points_in, vector<double> radii, double tunnel, long int cnts_t);
    int Check_repeated(vector<long int> region, long int Point);
    int Check_repeated(vector<int> region, int Point);
    int Delete_contact(int contact, vector<vector<int> > &contacts_vector);
    int Delete_contact(long int contact, vector<vector<long int> > &contacts_vector);
    void Clear_vectors();
    int Make_CNT_clusters(struct CNT_Geo cnts_geo, vector<Point_3D> points_in, vector<vector<long int> > structure);
    int Single_cluster(int cnt_seed, int start, vector<int> &vec, vector<vector<int> > &contacts_vector);
    int Single_cluster(int cnt_seed, int start, vector<int> &vec, vector<vector<int> > &contacts_vector, vector<short int> &cluster_flag);
    int Check_clusters_percolation(vector<Point_3D> points_in, vector<vector<long int> > structure);
    int Check_percolation_single_cluster(vector<int> cluster, int &family);
    int Check_percolation_single_cluster(vector<short int> cluster_flag, int &family);
    int Intersection(vector<int> vec1, vector<int> vec2);
    int Find_spheres(vector<Point_3D> points_in, vector<vector<long int> > structure);
    int Split_cnts(vector<Point_3D> &points_in, vector<double> &cnts_radius, vector<vector<long int> > &structure);
    int Get_LM_matrix(int family, vector<int> cluster, vector<vector<long int> > structure, vector<vector<long int> > &contacts_cnt_point, vector<int> &LM_matrix);
    int Sort_vector(vector<long int> &contacts, vector<long int> cnt);
    int Is_in_relevant_boundary(int family, int boundary_node);
    void Fill_sparse_stiffness_matrix(long int nodes, vector<int> cluster, vector<vector<long int> > structure, vector<int> LM_matrix, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, MathMatrix &R);
    void Add_to_sparse_stiffness(long int node1, long int node2, vector<vector<long int> > &col_ind_2d, vector<vector<double> > &values);
    MathMatrix Solve_DEA_equations_CG_SSS(long int nodes, vector<long int> col_ind, vector<long int> row_ptr, vector<double> values, vector<double> diagonal, MathMatrix R);
    MathMatrix spM_V_SSS(MathMatrix V, vector<long int> rowptr, vector<long int> colind, vector<double> diagonal, vector<double> values);
    double V_dot_v(MathMatrix v1, MathMatrix v2);
    int Remove_currentless_CNTs(MathMatrix voltages, vector<int> &cluster, vector<int> LM_matrix, vector<int> &dead, vector<double> &cnts_radius, vector<vector<long int> > &structure);
    int Trim_currenles_CNT(vector<vector<long int> > &structure, vector<double> &cnts_radius, vector<int> &cluster, vector<int> &dead, int CNT, long int P1, long int P2);
    int Clusters_length(vector<Point_3D> points, vector<vector<long int> > structure);
    int Export_visualization_files(const struct RVE_Geo &cell_geo, vector<Point_3D> points_in, vector<double> cnts_radius, vector<vector<long int> > structure);
    void Update_box_geometry();
    //===================================================================================================
    //===================================================================================================
    //===================================================================================================
    void Print1DVec(const vector<Point_3D> &list, const string &filename);
    void Print1DVec(const vector<char> &list, const string &filename);
    void Print1DVec(const vector<int> &list, const string &filename);
    void Print1DVec(const vector<double> &list, const string &filename);
    void Append1DVec(const vector<double> &list, const string &filename);
    void Print1DVec(const vector<long int> &list, const string &filename);
    void Print2DVec(const vector<vector<int> > &num_mat, const string &filename);
    void Print2DVec(const vector<vector<long int> > &num_mat, const string &filename);
    void Print2DVec(const vector<vector<double> > &num_mat, const string &filename);
    string Get_Line(ifstream &infile)const;
    //===================================================================================================
    //===================================================================================================
    //===================================================================================================
    int Export_cnt_networks_meshes(const struct RVE_Geo &cell, vector<vector<long int> > structure, vector<int> cluster, vector<Point_3D> points_in, vector<double> cnts_radius, const string &filename);
    int Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, vector<vector<long int> > structure, vector<int> cluster, vector<Point_3D> points_in, vector<double> cnts_radius);
    int Export_cnts_meshes_singlezone(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const string &filename)const;
    int Get_angles_vector_in_spherial_coordinates(const Point_3D &normal, double &sita, double &pha)const;
    int Get_points_circle_in_plane(const Point_3D &center, const double &trans_sita, const double &trans_pha, const double &radius, const int &num_sec, vector<Node> &nod_temp)const;
    int Get_projected_points_in_plane(const Point_3D &center, const Point_3D &normal, const Point_3D &line, const int &num_sec, vector<Node> &nod_temp)const;
    MathMatrix Get_transformation_matrix(const double &sita, const double &pha)const;
};


#endif /* defined(__Nanocode__RNetwork__) */
