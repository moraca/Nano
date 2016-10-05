//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Triangulation.cpp
//OBJECTIVE:	Generates a Delaunay triangulation from a given a list of points
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Triangulation.h"

//Function of the class that generates two triangulations, one for the top and one bottom surface of a GNP
int Triangulation::Generate_trangulations(const vector<Point_3D> &point_list, const vector<vector<long int> > &structure, const long int &GNP_top, const long int &GNP_bottom, GCH &hybrid)
{
    
    //Make sure the triangulation vectors of the hybrid are empty
    hybrid.triangulation_top.clear();
    hybrid.triangulation_bottom.clear();
    
    //Make the triangulation for the top surface
    if (!Single_surface_triangulation(1, GNP_top, point_list, structure, hybrid)) {
        hout << "Error in Generate_trangulations when calculating top triangulation" << endl;
        return 0;
    }
    
    //Make the triangulation for the bottom surface
    if (!Single_surface_triangulation(-1, GNP_bottom, point_list, structure, hybrid)) {
        hout << "Error in Generate_trangulations when calculating bottom triangulation" << endl;
        return 0;
    }    
    
    return 1;
}
//This function make a delaunay triangulation using as a starting triangulation the triangle given by the vector vertices
//This triangulation includes all initial points of CNTs on a GNP surface and the center of the GNP
int Triangulation::Single_surface_triangulation(const int &flag, const long int &GNP_center, const vector<Point_3D> &point_list, const vector<vector<long int> > &structure, GCH &hybrid)
{
    
    //Vector of points to be triangulated
    vector<long int> points_t;
    //Make a list of points to triangulate
    if (!Points_to_triangulate(flag, structure, points_t, hybrid)) {
        hout << "Error in Single_surface_triangulation when calling Points_to_triangulate" << endl;
        return 0;
    }
    
    //Variable for the vertices of the super triangle
    vector<Point_3D> vertices;
    //Add the center of the GNP to the vector vertices
    //vertices.push_back(hybrid.center);
    //Calculate the vertices of the super triangle for top surface and add them to the vertices vector
    if (!Generate_supertriangle(flag, hybrid, vertices)) {
        hout << "Error in Single_surface_triangulation when calling Single_surface_triangulation" << endl;
        return 0;
    }
    
    //Vector of triangles
    vector<vector<long int> > triangles;
    //Initialice triangulation with the triangulation given by the vertices of the supertriangle and the center of the GNP
    if (!Initialice_triangulation(triangles)) {
        hout << "Error in Single_surface_triangulation when calling Initialice_triangulation" << endl;
        return 0;
    }
    
    //Add points to the triangulation sequentially using the Bowyer-Watson algorithm
    if (!Bowyer_watson(flag, GNP_center, point_list, vertices, points_t, triangles, hybrid)) {
        hout << "Error in Single_surface_triangulation when calling Bowyer_watson" << endl;
        return 0;
    }
    
    return 1;
}
//This function gets all points to be triangulated depending on the value of an input flag:
//-1: Only bottom points are added
// 1: Only top points are added
// 0: Both bottom and top points are added
int Triangulation::Points_to_triangulate(const int &flag, const vector<vector<long int> > &structure, vector<long int> &points_out, GCH &hybrid)
{
    //If the flag is 0 or 1, the bottom points will be added
    if (flag >= 0) { //triangulate top surface
        if (!Points_to_triangulate_single_surface(structure, hybrid.cnts_top, points_out)) {
            hout << "Error in Points_to_triangulate when calling Points_to_triangulate_single_surface (top)" << endl;
            return 0;
        }
    }
    //If the flag is 0 or -1, the bottom points will be added
    if (flag <= 0) {//triangulate bottom surface
        if (!Points_to_triangulate_single_surface(structure, hybrid.cnts_bottom, points_out)) {
            hout << "Error in Points_to_triangulate when calling Points_to_triangulate_single_surface (bottom)" << endl;
            return 0;
        }
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
//flag=1 is for top surface
//flag=-1 is for bottom surface
int Triangulation::Generate_supertriangle(const int &flag, const GCH &hybrid, vector<Point_3D> &vertices)
{
    //dimensions of the GNP
    double lx = hybrid.gnp.len_x;
    double ly = hybrid.gnp.wid_y;
    
    //Length of the equilateral supertriangle
    double leq = 2*ly/sqrt(3) + lx;
    
    //z-coordinate
    double z = ((double)flag)*hybrid.gnp.hei_z/2;
    //Point to store the vertices of the super triangle on the plane
    Point_3D vertex(0.0,0.0,z);
    //Rotated vertex point
    Point_3D vertex_rot;
    
    //GNP surface center
    //Rotate vertex (the center of the GNP is the origin of coordinates)
    vertex_rot = vertex.rotation(hybrid.rotation, hybrid.center);
    //Add to vector of vertices
    vertices.push_back(vertex_rot);
    
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
//This function initialices the triangulation with the super triangle and the center point of the GNP
//This resutls in three triangles, the negative number represent:
//-1 is the center of the GNP
//-2,-3,-4 are the vertices of the super triangle
//These negative numbers correspond to the indices in the vertices vector
int Triangulation::Initialice_triangulation(vector<vector<long int> > &triangles)
{
    //Initial vector, it has the smallest point number; this is the center of GNP with vertex number -1
    vector<long int> initial;
    initial.push_back(-1);
    
    //Triangle 1
    triangles.push_back(initial);
    triangles.back().push_back(-2);
    triangles.back().push_back(-3);
    
    //Triangle 2
    triangles.push_back(initial);
    triangles.back().push_back(-2);
    triangles.back().push_back(-4);
    
    //Triangle 3
    triangles.push_back(initial);
    triangles.back().push_back(-3);
    triangles.back().push_back(-4);
    
    return 1;
}
//This is the implementation of the Bowyer Watson algorithm
int Triangulation::Bowyer_watson(const int &flag, const long int &GNP_center, const vector<Point_3D> &point_list, vector<Point_3D> &vertices, const vector<long int> &points_t, vector<vector<long int> > &triangles, GCH &hybrid)
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
    if ( (flag>=0) && (Final_triangulation_edges(flag, GNP_center, triangles, hybrid.triangulation_top)) )
        return 1;
    else if ( (flag==-1) && (Final_triangulation_edges(flag, GNP_center, triangles, hybrid.triangulation_bottom)) )
        return 1;
    else
        return 0;
    
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
//It also substitutes the encoded index number (-1) of the GNP center by the one used by the DEA
int Triangulation::Final_triangulation_edges(const int &flag, const long int &GNP_center, const vector<vector<long int> > &triangles, vector<vector<long int> > &edges)
{
    //Scan all triangles looking for valid edges
    for (int i = 0; i < (int)triangles.size(); i++) {
        //Check if there are any valid edges and add them to the vector of edges
        if (!Valid_edges(flag, GNP_center, triangles[i], edges)) {
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
//It also substitutes the encoded index number (-1) of the GNP center by the one used by the DEA
int Triangulation::Valid_edges(const int &flag, const long int &GNP_center, const vector<long int> &triangle, vector<vector<long int> > &edges)
{
    //Temporary vector
    vector<long int> edge_tmp;
    
    //scan the vertices in the triangle
    for (int i = 0; i < (int)triangle.size(); i++) {
        //if a vertex is -2,-3 or -4 it is a supertriangle vertex so ignore it
        if (triangle[i] >= 0)
            edge_tmp.push_back(triangle[i]);
        //Substitute the -1 by the GNP number given by the DEA
        else if ( flag!=0 && triangle[i]==-1)
            edge_tmp.push_back(GNP_center);
    }
    
    if (edge_tmp.size()==2)
        //In this case there is only one valid edge
        edges.push_back(edge_tmp);
    else if (edge_tmp.size()==3)
        //In this case the three vertices are valid so use the function that adds a triangle as edges
        Add_triangle_as_edges(edge_tmp, edges);
    
    return 1;
}
//===========================================================================================================================
//===========================================================================================================================
//===========================================================================================================================
//Function to make the 3D triangulation
int Triangulation::Generate_3d_trangulation(const vector<Point_3D> &point_list, const vector<vector<long int> > &structure, GCH &hybrid)
{
    //Make sure the triangulation vectors of the hybrid are empty
    hybrid.triangulation_top.clear();
    hybrid.triangulation_bottom.clear();
    
    //Vector of points to be triangulated
    vector<long int> points_t;
    //Get all points from bottom and top surfaces
    if (!Points_to_triangulate(0, structure, points_t, hybrid)) {
        hout << "Error in Single_surface_triangulation when calling Points_to_triangulate" << endl;
        return 0;
    }
    
    /*/------------------------------------------------------------------------------------------------------------------
    //USING A SUPER TETRAHEDRON
    //Generate supertetrahedron
    vector<Point_3D> vertices;
    if (!Generate_supertetrahedron(hybrid, vertices)) {
        hout << "Error in Single_surface_triangulation when calling Generate_supertetrahedron" << endl;
        return 0;
    }
    //Vector of triangles
    vector<vector<long int> > triangles;
    //Initialice triangulation with the triangulation given by the vertices of the supertetrahedron
    if (!Initialice_triangulation(triangles)) {
        hout << "Error in Single_surface_triangulation when calling Initialice_triangulation" << endl;
        return 0;
    }//*/
    
    //------------------------------------------------------------------------------------------------------------------
    //USING A SUPER TRIANGLE
    vector<Point_3D> vertices;
    if (!Generate_supertriangle(1, hybrid, vertices)) {
        hout << "Error in Single_surface_triangulation when calling Generate_supertetrahedron" << endl;
        return 0;
    }
    //Remove GNP center as it is not used in the 3D triangulation
    vertices.erase(vertices.begin());
    //Vector of triangles
    vector<vector<long int> > triangles;
    vector<long int> tmp;
    //The initial triangulation consists only of the supertriangle
    tmp.push_back(-1);tmp.push_back(-2);tmp.push_back(-3);
    triangles.push_back(tmp);//*/
    
    //Add points to the triangulation sequentially using the Bowyer-Watson algorithm
    if (!Bowyer_watson(0, 0, point_list, vertices, points_t, triangles, hybrid)) {
        hout << "Error in Single_surface_triangulation when calling Bowyer_watson" << endl;
        return 0;
    }
    
    return 1;
}
//Generate the initial triangulation
int Triangulation::Generate_supertetrahedron(const GCH &hybrid, vector<Point_3D> &vertices)
{
    //dimensions of the GNP
    double lx = hybrid.gnp.len_x;
    double ly = hybrid.gnp.wid_y;
    
    //Length of the equilateral supertriangle
    double leq = 2*ly/sqrt(3) + lx;
    
    //Point to store the vertices of the super triangle on the plane
    Point_3D vertex;
    //Rotated vertex point
    Point_3D vertex_rot;
    //Variable to store the unit vectors
    Point_3D unit;
    
    //Heigh vertex D
    //x coordinate has the same value as in B, but inverted sign
    vertex.x = 0.0;
    vertex.y = sqrt(3)*leq/6 - ly/2;
    vertex.z = sqrt(2/3)*leq;
    //Rotate vertex (the center of the GNP is the origin of coordinates)
    vertex_rot = vertex.rotation(hybrid.rotation, hybrid.center);
    //Add to vector of vertices
    vertices.push_back(vertex_rot);
    
    //Base vertex A
    vertex.y = 0.5*(ly + sqrt(3)*lx); //x is (still) zero
    vertex.z = hybrid.gnp.hei_z/2; //Vertices A, B and C have this z-coordinate value (initially)
    //Rotate vertex (the center of the GNP is the origin of coordinates)
    vertex_rot = vertex.rotation(hybrid.rotation, hybrid.center);
    //Unit vector DA
    unit = (vertex_rot - vertices.front())/(vertex_rot.distance_to(vertices.front()));
    //Add new location to vector of vertices
    vertices.push_back(vertices.front()+(unit*2*leq));
    
    //Base vertex B
    vertex.x = -0.5*leq;
    vertex.y = -0.5*ly;
    vertex_rot = vertex.rotation(hybrid.rotation, hybrid.center);
    //Unit vector DB
    unit = (vertex_rot - vertices.front())/(vertex_rot.distance_to(vertices.front()));
    //Add new location to vector of vertices
    vertices.push_back(vertices.front()+(unit*2*leq));
    
    //Base vertex C
    //x coordinate has the same value as in B, but inverted sign
    vertex.x = -vertex.x; //y coordinate remains the same as in B
    //Rotate vertex (the center of the GNP is the origin of coordinates)
    vertex_rot = vertex.rotation(hybrid.rotation, hybrid.center);
    //Unit vector DB
    unit = (vertex_rot - vertices.front())/(vertex_rot.distance_to(vertices.front()));
    //Add new location to vector of vertices
    vertices.push_back(vertices.front()+(unit*2*leq));
    
    return 1;
}