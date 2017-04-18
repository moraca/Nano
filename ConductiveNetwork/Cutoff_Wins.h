//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Cutoff_Wins.h
//OBJECTIVE:	To cutoff the windows
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef CUTOFFWINS_H
#define CUTOFFWINS_H

#include "Input_Reader.h"
#include "Background_vectors.h"
#include "Printer.h"

//-------------------------------------------------------
class Cutoff_Wins
{
public:
    //Data Member
    vector<int> cnts_inside; //List of CNTs inside the observation window
    vector<int> gnps_inside; //List of GNPs inside the observation window
    vector<vector<short int> > boundary_flags_cnt; //This vector will help find points on the boundary. This is used in the direct electrifying algorithm
    vector<vector<short int> > boundary_flags_gnp; //This vector will help find points on the boundary. This is used in the direct electrifying algorithm
    vector<vector<int> > boundary_cnt; //Boundary vector. It is used to determine percolation
    vector<vector<int> > boundary_gnp; //Boundary vector. It is used to determine percolation
    double xmin, ymin, zmin;
    double w_x, w_y, w_z;

    
    //Constructor
    Cutoff_Wins(){};
    
    //Member Functions
    int Extract_observation_window(const int &window, const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const struct GNP_Geo &gnps, vector<GCH> &hybrid_particles, vector<vector<long int> > &structure, vector<vector<long int> > &structure_gnp, vector<double> &radii, vector<Point_3D> &points_in, vector<Point_3D> &points_gnp, vector<vector<int> > &shells_cnt, vector<vector<int> > &shells_gnp);
    int Set_global_variables_for_geometry(const struct Geom_RVE &sample, const int &window);
    int Save_seeds(const vector<GCH> &hybrid_particles, const vector<vector<long int> > &structure, vector<long int> &seeds);
    int Compare_seeds(vector<GCH> &hybrid_particles, const vector<vector<long int> > &structure, const vector<long int> &seeds);
    int Trim_boundary_cnts(vector<vector<int> > &shells_cnt, const int &window, const struct Geom_RVE &sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<double> &radii);
    int First_index(vector<Point_3D> &points_in, vector<long int> &structure_CNT, int &new_CNT, int &index1);
    int Second_index(vector<Point_3D> &points_in, vector<long int> &structure_CNT, int &new_CNT, int &index2);
    string Where_is(Point_3D point);
    int New_boundary_point(struct Geom_RVE sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, long int insidePoint, long int outsidePoint, int CNT, string currentLocation);
    int Substitute_boundary_point(vector<Point_3D> &points_in, long int global_i, long int global_o);
    int Get_intersecting_point_RVE_surface(const Point_3D &point0, const Point_3D &point1, vector<Point_3D> &ipoi_vec);
    void Add_to_boundary_vectors(Point_3D point3d, long int point, int new_CNT);
    void Add_CNT_to_boundary(vector<int> &boundary, int CNT, long int point, short int flag1, short int flag2);
    int Change_repeated_seed(const int &CNT_original, int &index1_previous, int &index2_previous, int &index1_current, vector<vector<long int> > &structure, vector<Point_3D> &points_in);
    int Fill_cnts_inside(const vector<vector<long int> > &structure);
    int Trim_boundary_gnps(const struct GNP_Geo &gnps, const vector<GCH> &hybrid_particles, const vector<int> &shell_gnp, vector<Point_3D> &points_gnp, vector<vector<long int> > &structure_gnp);
    int Is_close_to_boundaries(const GCH &hybrid);
    int Remove_gnp_points_outside(const struct GNP_Geo &gnps, const GCH &hybrid, vector<Point_3D> &points_gnp, vector<long int> &gnp_discrete);
    int Find_inside_point(const GCH &hybrid, Point_3D &inside_point);
    int Find_projection_in_boundary(const Point_3D &inside, Point_3D &outside, const int &iterator, vector<vector<int> > &boundaries);
    int Add_GNPs_to_boundary(const vector<long int> &gnp_discrete, const vector<vector<int> > &boundaries, vector<Point_3D> &points_gnp);
    int Find_average_boundary_point(const vector<long int> &gnp_discrete, const vector<int> &boundary, vector<Point_3D> &points_gnp);
    int Update_discretization(vector<long int> &gnp_discrete, vector<long int> &gnp_discrete_in);
    int Fill_gnps_inside(const vector<vector<long int> > &structure_gnp);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================