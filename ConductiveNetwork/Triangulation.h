//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Triangulation.h
//OBJECTIVE:	Generates a Delaunay triangulation from a given a list of points
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include "Input_Reader.h"
#include "GenNetwork.h"

//-------------------------------------------------------
class Triangulation
{
public:
    //Data Member
    
    //Constructor
    Triangulation(){};
    
    //3D triangulation
    int Generate_3d_trangulation(const long int &last_cnt_node, const vector<Point_3D> &point_list, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list_gnp, const vector<long int> &gnp_junctions, GCH &hybrid);
    int Points_to_triangulate(const long int &last_cnt_node, const vector<Point_3D> &point_list, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list_gnp, const vector<long int> &gnp_junctions, vector<Point_3D> &points_out_3d, vector<long int> &points_out, vector<short int> &points_out_flags, GCH &hybrid);
    int Cnt_points_to_triangulate(const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<int> &cnt_list, vector<long int> &points_out, vector<Point_3D> &points_out_3d);
    int Gnp_points_to_triangulate(const int &n_cnt_points, const vector<Point_3D> &point_list_gnp, const vector<long int> &gnp_junctions, vector<Point_3D> &points_out_3d, vector<long int> &points_out, vector<short int> &points_out_flags);
    int Generate_trivial_triangulation(const vector<long int> &points_out, const vector<short int> &points_out_flags, GCH &hybrid);
    int Generate_supertriangle(const GCH &hybrid, vector<Point_3D> &vertices);
    int Bowyer_watson(const vector<Point_3D> &points_t_3d, vector<Point_3D> &vertices_s, const vector<long int> &points_t, const vector<short int> &points_t_flags, const vector<int> &vertices_t, vector<vector<int> > &triangles, GCH &hybrid);
    int Find_bad_triangles(const vector<Point_3D> &points_t, vector<Point_3D> &vertices_s, const int &point, vector<vector<int> > &triangles, vector<vector<int> > &bad_triangles_edges);
    int Is_in_circumcircle(const int &point, const vector<Point_3D> &points_t, vector<Point_3D> &vertices_s, const vector<int> &triangle);
    Point_3D Calculate_circumcenter(const vector<Point_3D> &points_t, vector<Point_3D> &vertices_s, const vector<int> &triangle);
    Point_3D Get_point(const int &index, const vector<Point_3D> &points_t, vector<Point_3D> &vertices);
    void Add_triangle_as_edges(const vector<int> &triangle, vector<vector<int> > &edges);
    int Add_new_triangles(const int &point, vector<vector<int> > &triangles, vector<vector<int> > &bad_triangles_edges);
    int Is_same_edge(const vector<int> &edge1, const vector<int> &edge2);
    int Final_triangulation_edges(const vector<vector<int> > &triangles, const vector<long int> &points_t, const vector<short int> &points_t_flags, const vector<int> &vertices_t, GCH &hybrid);
    int Valid_edges(const vector<int> &triangle, vector<vector<int> > &edges);   
};
//-------------------------------------------------------
#endif
//===========================================================================