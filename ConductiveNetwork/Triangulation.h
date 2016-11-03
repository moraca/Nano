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
    int Generate_3d_trangulation(const vector<Point_3D> &point_list, const vector<vector<long int> > &structure, GCH &hybrid);
    //2D triangulation
    int Generate_trangulations(const vector<Point_3D> &point_list, const vector<vector<long int> > &structure, const long int &GNP_top, const long int &GNP_bottom, GCH &hybrid);
    int Points_to_triangulate(const vector<vector<long int> > &structure, vector<long int> &points_out, GCH &hybrid);
    int Points_to_triangulate_single_surface(const vector<vector<long int> > &structure, const vector<int> &cnt_list, vector<long int> &points_out);
    int Generate_supertriangle(const GCH &hybrid, vector<Point_3D> &vertices);
    int Bowyer_watson(const vector<Point_3D> &point_list, vector<Point_3D> &vertices, const vector<long int> &points_t, vector<vector<long int> > &triangles, GCH &hybrid);
    int Find_bad_triangles(const vector<Point_3D> &point_list, vector<Point_3D> &vertices, const long int &point, vector<vector<long int> > &triangles, vector<vector<long int> > &bad_triangles_edges);
    int Is_in_circumcircle(const vector<Point_3D> &point_list, vector<Point_3D> &vertices, const vector<long int> &triangle, long int point);
    Point_3D Get_point(const vector<Point_3D> &point_list, vector<Point_3D> &vertices, long int index);
    Point_3D Calculate_circumcenter(const vector<Point_3D> &point_list, vector<Point_3D> &vertices, const vector<long int> &triangle);
    void Add_triangle_as_edges(const vector<long int> &triangle, vector<vector<long int> > &edges);
    int Add_new_triangles(const long int &point, vector<vector<long int> > &triangles, vector<vector<long int> > &bad_triangles_edges);
    int Is_same_edge(const vector<long int> &edge1, const vector<long int> &edge2);
    int Final_triangulation_edges(const vector<vector<long int> > &triangles, vector<vector<long int> > &edges);
    int Valid_edges(const vector<long int> &triangle, vector<vector<long int> > &edges);    
};
//-------------------------------------------------------
#endif
//===========================================================================