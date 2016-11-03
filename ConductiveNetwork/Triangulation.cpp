//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Triangulation.cpp
//OBJECTIVE:	Generates a Delaunay triangulation from a given a list of points
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Triangulation.h"

//Function to make the 3D triangulation
int Triangulation::Generate_3d_trangulation(const vector<Point_3D> &point_list, const vector<vector<long int> > &structure, GCH &hybrid)
{
    //Make sure the triangulation vectors of the hybrid are empty
    hybrid.triangulation.clear();
    
    //Vector of points to be triangulated
    vector<long int> points_t;
    //Get all points from bottom and top surfaces
    if (!Points_to_triangulate(structure, points_t, hybrid)) {
        hout << "Error in Generate_3d_trangulation when calling Points_to_triangulate" << endl;
        return 0;
    }
    
    //------------------------------------------------------------------------------------------------------------------
    //SUPER TRIANGLE
    vector<Point_3D> vertices;
    if (!Generate_supertriangle(hybrid, vertices)) {
        hout << "Error in Generate_3d_trangulation when calling Generate_supertetrahedron" << endl;
        return 0;
    }
    //Vector of triangles
    vector<vector<long int> > triangles;
    vector<long int> tmp;
    //The initial triangulation consists only of the supertriangle
    tmp.push_back(-1);tmp.push_back(-2);tmp.push_back(-3);
    triangles.push_back(tmp);
    
    //Add points to the triangulation sequentially using the Bowyer-Watson algorithm
    if (!Bowyer_watson(point_list, vertices, points_t, triangles, hybrid)) {
        hout << "Error in Generate_3d_trangulation when calling Bowyer_watson" << endl;
        return 0;
    }
    
    return 1;
}
//This function gets all points to be triangulated
//It first adds the CNT seeds from the top surface, then the CNT seeds of the bottoms surface
int Triangulation::Points_to_triangulate(const vector<vector<long int> > &structure, vector<long int> &points_out, GCH &hybrid)
{
    //triangulate top surface
    if (!Points_to_triangulate_single_surface(structure, hybrid.cnts_top, points_out)) {
        hout << "Error in Points_to_triangulate when calling Points_to_triangulate_single_surface (top)" << endl;
        return 0;
    }
    
    //triangulate bottom surface
    if (!Points_to_triangulate_single_surface(structure, hybrid.cnts_bottom, points_out)) {
        hout << "Error in Points_to_triangulate when calling Points_to_triangulate_single_surface (bottom)" << endl;
        return 0;
    }
    
    return 1;
}
//Gather in a single vector all points to be triangulated: all initial points of CNTs attached to GNP
int Triangulation::Points_to_triangulate_single_surface(const vector<vector<long int> > &structure, const vector<int> &cnt_list, vector<long int> &points_out)
{
    //Add intial points of CNTs attached to the GNP
    for (int i = 0; i < (int)cnt_list.size(); i++) {
        //current CNT
        int CNT = cnt_list[i];
        
        //Initial point of current CNT
        long int P = structure[CNT].front();
        
        //Add to vector of points
        points_out.push_back(P);
    }
    
    return 1;
}
//Generate supertriangle for a GNP surface
int Triangulation::Generate_supertriangle(const GCH &hybrid, vector<Point_3D> &vertices)
{
    //dimensions of the GNP
    double lx = hybrid.gnp.len_x;
    double ly = hybrid.gnp.wid_y;
    
    //Length of the equilateral supertriangle
    double leq = 2*ly/sqrt(3) + lx;
    
    //Point to store the vertices of the super triangle on the plane
    Point_3D vertex(0.0,0.0,hybrid.gnp.hei_z/2);
    //Rotated vertex point
    Point_3D vertex_rot;
    
    //Base vertex A
    vertex.y = 0.5*(ly + sqrt(3)*lx); //x is (still) zero
    //Rotate vertex (the center of the GNP is the origin of coordinates)
    vertex_rot = vertex.rotation(hybrid.rotation, hybrid.center);
    //Add to vector of vertices
    vertices.push_back(vertex_rot);
    
    //Base vertex B
    vertex.x = -0.5*leq;
    vertex.y = -0.5*ly;
    vertex_rot = vertex.rotation(hybrid.rotation, hybrid.center);
    //Add to vector of vertices
    vertices.push_back(vertex_rot);
    
    //Base vertex C
    //x coordinate has the same value as in B, but inverted sign
    vertex.x = -vertex.x; //y coordinate remains the same as in B
    //Rotate vertex (the center of the GNP is the origin of coordinates)
    vertex_rot = vertex.rotation(hybrid.rotation, hybrid.center);
    //Add to vector of vertices
    vertices.push_back(vertex_rot);
    
    return 1;
}
//This is the implementation of the Bowyer Watson algorithm
int Triangulation::Bowyer_watson(const vector<Point_3D> &point_list, vector<Point_3D> &vertices, const vector<long int> &points_t, vector<vector<long int> > &triangles, GCH &hybrid)
{
    //Add each point to the triangulation
    for (int i = 0; i < (int)points_t.size(); i++) {
        //Vector for bad triangles
        vector<vector<long int> > bad_triangles_edges;
        
        //Find the triangles whose circumcircle contains the current point and remove them from the triangulation
        if (!Find_bad_triangles(point_list, vertices, points_t[i], triangles, bad_triangles_edges)) {
            hout << "Error in Bowyer_watson when calling Find_bad_triangles" << endl;
            return 0;
        }
        
        //Find the edges that will form new triangles and add them to the triangulation
        if (!Add_new_triangles(points_t[i], triangles, bad_triangles_edges)) {
            hout << "Error in Bowyer_watson when calling Add_triangle_as_edges" << endl;
            return 0;
        }
        
    }
    
    //Triangulation is finished, so find the edges that make up the final triangulation:
    // -Remove those edges that contain a vertex of the supertriangle
    // -Remove repeated edges
    if ( !Final_triangulation_edges(triangles, hybrid.triangulation) ) {
        hout << "Error in Bowyer_watson when calling Final_triangulation_edges" <<endl;
        return 0;
}

    return 1;
    
}
//This function finds the bad triangles in a triangulation, i.e., the triangles that contain a point in their circumcircle
//The bad triangles are removed from the triangulation and added as edges
int Triangulation::Find_bad_triangles(const vector<Point_3D> &point_list, vector<Point_3D> &vertices, const long int &point, vector<vector<long int> > &triangles, vector<vector<long int> > &bad_triangles_edges)
{
    //Scan all triangles
    for (int j = (int)triangles.size()-1; j >= 0 ; j--) {
        //If point is inside the circumcircle then add the triangle to bad_triangles
        if (Is_in_circumcircle(point_list, vertices, triangles[j], point)) {
            //Add current triangle to bad_triangles as edges
            Add_triangle_as_edges(triangles[j], bad_triangles_edges);
            //Remove triangle from triangulation
            triangles.erase(triangles.begin()+j);
        }
    }
    
    return 1;
}
//This function checks if a point is inside the circumcircle of a triangle
int Triangulation::Is_in_circumcircle(const vector<Point_3D> &point_list, vector<Point_3D> &vertices, const vector<long int> &triangle, long int point)
{
    //Calculate circumcenter of triangle
    Point_3D C = Calculate_circumcenter(point_list, vertices, triangle);
    
    //Calculate squared radius of circumcircle (i.e. squared distance from center to any point)
    //Comparing squared distances will save time by avoiding calculating squared roots
    Point_3D P = Get_point(point_list, vertices, triangle[0]);
    double rad2 = C.squared_distance_to(P);
    
    //Squared distance from circumcenter to point
    P = Get_point(point_list, vertices, point);
    double dist = C.squared_distance_to(P);
    
    //If the squared distance between the point and C is smaller or equal than the squared radius,
    //then the point is inside the circumcircle so return 1
    return (dist <= rad2);
}
//There are some negative indices to refer to the vertices vectors while all positive indices refer to the point_list vector
//so this function deals with that and returns the proper point
Point_3D Triangulation::Get_point(const vector<Point_3D> &point_list, vector<Point_3D> &vertices, long int index)
{
    //If the index is equal or greater to zero, then the index remains the same
    if (index >= 0)
        return point_list[index];
    //If the index is negative, add one and change the sign to get the correct index in vertices
    else
        return vertices[-1-index];
}
//Given a triangle, this function calculates its circumcenter
Point_3D Triangulation::Calculate_circumcenter(const vector<Point_3D> &point_list, vector<Point_3D> &vertices, const vector<long int> &triangle)
{
    //Variables for calculations
    Point_3D P, Q, R;
    
    //Triangle points
    Point_3D P1 = Get_point(point_list, vertices, triangle[0]);
    Point_3D P2 = Get_point(point_list, vertices, triangle[1]);
    Point_3D P3 = Get_point(point_list, vertices, triangle[2]);
    
    //Calculate points P, Q and R = PxQ
    P = P1 - P3;
    Q = P2 - P3;
    R = P.cross(Q);
    
    //Calculate squared norms
    double Pn = P.dot(P);
    double Qn = Q.dot(Q);
    double Rn = R.dot(R);
    
    //Variable to store the circumcenter
    Point_3D C;
    
    //Calculate circumcircle
    C = (Q*Pn - P*Qn).cross(R)/(2*Rn);
    C = C + P3;
    
    return C;
}
//This function adds a triangle to the vector edges as the three edges that make up the triangle
void Triangulation::Add_triangle_as_edges(const vector<long int> &triangle, vector<vector<long int> > &edges)
{
    //Empty vector to add elements to edges
    vector<long int> empty;
    
    //Add edges to last vector of edges
    //Edge 1,2
    edges.push_back(empty);
    edges.back().push_back(triangle[0]);
    edges.back().push_back(triangle[1]);
    //Edge 1,3
    edges.push_back(empty);
    edges.back().push_back(triangle[0]);
    edges.back().push_back(triangle[2]);
    //Edge 2,3
    edges.push_back(empty);
    edges.back().push_back(triangle[1]);
    edges.back().push_back(triangle[2]);
    
}
//This function finds the edges that will make the new triangles in the triangulation
//then they new triangles are added to the triangulation with the edges found and a given point
int Triangulation::Add_new_triangles(const long int &point, vector<vector<long int> > &triangles, vector<vector<long int> > &bad_triangles_edges)
{
    //Scan each edge in bad triangles (except the last one) and compare it with the rest
    for (int j = 0; j < (int)bad_triangles_edges.size()-1; j++) {
        //Falg that indicates if an edge is shared (repeated)
        int shared = 0;
        for (int k=j+1; k < (int)bad_triangles_edges.size(); k++) {
            //Check if the two edges are the same
            if (Is_same_edge(bad_triangles_edges[j], bad_triangles_edges[k])) {
                //Change the value of the shared flag
                shared = 1;
                //Delete the repeated edge
                bad_triangles_edges.erase(bad_triangles_edges.begin()+k);
                //There is no need to continue checking so break the loop
                break;
            }
        }
        
        //Check if the edge was not shared
        if (!shared) {
            //If the edge is not shared, the two points of the edge and the current point make a new triangle
            triangles.push_back(bad_triangles_edges[j]);
            triangles.back().push_back(point);
        }
        
        //Check if it reached the last element in the iteration, that is if j==bad_triangles_edges.size()-2
        if ( j==((int)bad_triangles_edges.size()-2) ) {
            //If it reached this index value, it means that the last edge is was not repeated but was not scanned
            //Thus add it as a new triangle
            triangles.push_back(bad_triangles_edges.back());
            triangles.back().push_back(point);
        }
    }
    
    return 1;
}
//This function compares two edges and decides if they are equal or not
int Triangulation::Is_same_edge(const vector<long int> &edge1, const vector<long int> &edge2)
{
    return ( (edge1[0]==edge2[0] && edge1[1]==edge2[1]) || (edge1[0]==edge2[1] && edge1[1]==edge2[0]) );
}
//This function removes the edges that have vertices of the super triangle and
//generates the vector of edges that make up the triangulation
int Triangulation::Final_triangulation_edges(const vector<vector<long int> > &triangles, vector<vector<long int> > &edges)
{
    //Scan all triangles looking for valid edges
    for (int i = 0; i < (int)triangles.size(); i++) {
        //Check if there are any valid edges and add them to the vector of edges
        if (!Valid_edges(triangles[i], edges)) {
            hout << "Error in Final_triangulation_edges" << endl;
            return 0;
        }
    }
    
    //Delete repeated edges
    for (int i = 0; i < (int)edges.size()-1; i++) {
        for (int j = i+1; j < (int)edges.size(); j++) {
            //If the edge is repeated, then delete it
            if (Is_same_edge(edges[i], edges[j])) {
                edges.erase(edges.begin()+j);
                //break the loop, an edge cannot be shared by more than two triangles
                //break;
                j--; //I found that I had repeated edges and adding this seemed to help
            }
        }
    }
    
    return 1;
}
//Make edges using only valid edges, i.e., excluding the vertices of the supertriangle
int Triangulation::Valid_edges(const vector<long int> &triangle, vector<vector<long int> > &edges)
{
    //Temporary vector
    vector<long int> edge_tmp;
    
    //scan the vertices in the triangle
    for (int i = 0; i < (int)triangle.size(); i++) {
        //if a vertex is -1,-2 or -3 it is a supertriangle vertex so ignore it
        if (triangle[i] >= 0)
            edge_tmp.push_back(triangle[i]);
    }
    
    if (edge_tmp.size()==2)
        //In this case there is only one valid edge
        edges.push_back(edge_tmp);
    else if (edge_tmp.size()==3)
        //In this case the three vertices are valid so use the function that adds a triangle as edges
        Add_triangle_as_edges(edge_tmp, edges);
    
    return 1;
}
