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

//-------------------------------------------------------
class Cutoff_Wins
{
public:
    //Data Member
    vector<int> cnts_inside; //List of CNTs inside the observation window
    vector<vector<short int> > boundary_flags; //This vector will help find points on the boundary. This is used in the direct electrifying algorithm
    vector<vector<int> > boundary_cnt; //Boundary vector. It is used to determine percolation
    vector<vector<int> > sectioned_domain;
    double xmin, ymin, zmin;
    double w_x, w_y, w_z;

    
    //Constructor
    Cutoff_Wins(){};
    
    //Member Functions
    int Extract_observation_window(struct Geom_RVE sample, struct Nanotube_Geo cnts, vector<vector<long int> > &structure, vector<double> &radii, vector<Point_3D> &points_in, vector<vector<int> > &shells_cnt, int window);
    int Trim_boundary_cnts(vector<vector<int> > &shells_cnt, int window, struct Geom_RVE sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<double> &radii);
    string Where_is(Point_3D point);
    int New_boundary_point(struct Geom_RVE sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, long int insidePoint, long int outsidePoint, int CNT, string currentLocation);
    int Substitute_boundary_point(vector<Point_3D> &points_in, long int global_i, long int global_o);
    int Get_intersecting_point_RVE_surface(Point_3D &point0, Point_3D &point1, vector<Point_3D> &ipoi_vec);
    void Trim_CNT(vector<vector<int> > &shells_cnt, struct Geom_RVE sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<double> &radii, long int boundary, int CNT);
    void Add_to_boundary_vectors(struct Geom_RVE sample, Point_3D point3d, long int point);
    void Add_CNT_to_boundary(vector<int> &boundary, int CNT, long int point, short int flag1, short int flag2);
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================