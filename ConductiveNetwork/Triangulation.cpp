//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Triangulation.cpp
//OBJECTIVE:	Generates a Delaunay triangulation from a given a list of points
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Triangulation.h"


//Function to make the 3D triangulation
int Triangulation::Generate_3d_trangulation(const long int &last_cnt_node, const vector<Point_3D> &point_list, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list_gnp, const vector<long int> &gnp_junctions, GCH &hybrid)
{
    
    //Make sure the triangulation vector of the hybrid is empty
    hybrid.triangulation.clear();
    //Make sure the triangulation flags vector of the hybrid is empty
    hybrid.triangulation_flags.clear();
    
    //------------------------------------------------------------------------------------------------------------------
    //Vectors of points to be triangulated
    vector<long int> points_t;
    //flags: CNT point (1) or a GNP point (0)
    vector<short int> points_t_flags;
    vector<Point_3D> points_t_3d;
    //Get all points from the CNT seeds (bottom and top surfaces) and GNP junctions
    if (!Points_to_triangulate(last_cnt_node, point_list, structure, point_list_gnp, gnp_junctions, points_t_3d, points_t, points_t_flags, hybrid)) {
        hout << "Error in Generate_3d_trangulation when calling Points_to_triangulate" << endl;
        return 0;
    }
    
    //Check if the number of points to triangulate is 3 or less, which would be a trivial case
    if (points_t.size() <= 3) {
        
        //If there are 3 points or less, make a trivial triangulation
        if (!Generate_trivial_triangulation(points_t, points_t_flags, hybrid)) {
            hout << "Error in Generate_3d_trangulation when calling Generate_trivial_triangulation" << endl;
            return 0;
        }
    }
    //If more than 3 points, make the triangulation using the Bowyer-Watson algorithm
    else {
        
        //Create a vector of consecutive numbers that will represent the points in points_t
        vector<int> vertices_t(points_t.size(),0);
        for (int i = 1; i < (int)vertices_t.size(); i++) {
            vertices_t[i] = i;
        }
        
        //------------------------------------------------------------------------------------------------------------------
        //SUPER TRIANGLE
        vector<Point_3D> vertices_s;
        if (!Generate_supertriangle(hybrid, vertices_s)) {
            hout << "Error in Generate_3d_trangulation when calling Generate_supertetrahedron" << endl;
            return 0;
        }
        //Vector of triangles
        vector<vector<int> > triangles;
        vector<int> tmp;
        //The initial triangulation consists only of the supertriangle
        tmp.push_back(-1);tmp.push_back(-2);tmp.push_back(-3);
        triangles.push_back(tmp);
        
        //------------------------------------------------------------------------------------------------------------------
        //Add points to the triangulation sequentially using the Bowyer-Watson algorithm
        if (!Bowyer_watson(points_t_3d, vertices_s, points_t, points_t_flags, vertices_t, triangles, hybrid)) {
            hout << "Error in Generate_3d_trangulation when calling Bowyer_watson" << endl;
            return 0;
        }
        
    }
    
    return 1;
}
//This function gets all points to be triangulated
//It first adds the CNT seeds from the top surface, then the CNT seeds of the bottoms surface, and finally the GNP junction points
int Triangulation::Points_to_triangulate(const long int &last_cnt_node, const vector<Point_3D> &point_list, const vector<vector<long int> > &structure, const vector<Point_3D> &point_list_gnp, const vector<long int> &gnp_junctions, vector<Point_3D> &points_out_3d, vector<long int> &points_out, vector<short int> &points_out_flags, GCH &hybrid)
{
    
    //Number of CNT points to triangulate
    int n_cnt_points = 0;
    
    //triangulation points for top surface
    n_cnt_points = Cnt_points_to_triangulate(structure, point_list, hybrid.cnts_top, points_out, points_out_3d);
    
    //triangulation points for bottom surface
    n_cnt_points += Cnt_points_to_triangulate(structure, point_list, hybrid.cnts_bottom, points_out, points_out_3d);
    
    //GNP triangulation points
    if (!Gnp_points_to_triangulate(n_cnt_points, point_list_gnp, gnp_junctions, points_out_3d, points_out, points_out_flags)) {
        hout << "Error in Points_to_triangulate when calling Gnp_points_to_triangulate" << endl;
        return 0;
    }
    
    return 1;
}
//Gather in a single vector all points to be triangulated: all initial points of CNTs attached to GNP
int Triangulation::Cnt_points_to_triangulate(const vector<vector<long int> > &structure, const vector<Point_3D> &point_list, const vector<int> &cnt_list, vector<long int> &points_out, vector<Point_3D> &points_out_3d)
{
    //Number of points to triangulate
    int n_points = 0;
    
    //Add intial points of CNTs attached to the GNP
    for (int i = 0; i < (int)cnt_list.size(); i++) {
        //current CNT
        int CNT = cnt_list[i];
        
        //Initial point of current CNT
        long int P = structure[CNT].front();
        
        //Add to vector of points
        points_out.push_back(P);
        
        //Add to vector of 3D points
        points_out_3d.push_back(point_list[P]);
        
        //Increase the count of points
        n_points++;
    }
    
    return n_points;
}
//flags: CNT point (1) or a GNP point (0)
int Triangulation::Gnp_points_to_triangulate(const int &n_cnt_points, const vector<Point_3D> &point_list_gnp, const vector<long int> &gnp_junctions, vector<Point_3D> &points_out_3d, vector<long int> &points_out, vector<short int> &points_out_flags)
{
    
    //Number of points to triangulate
    int n_points = 0;
    
    //Add GNP junction points
    for (int i = 0; i < (int)gnp_junctions.size(); i++) {
        
        //Current point
        long int P = gnp_junctions[i];
        
        //Add to vector of points
        points_out.push_back(P);
        
        //Add to vector of 3D points
        points_out_3d.push_back(point_list_gnp[P]);
        
        //Increase the count of points
        n_points++;
    }
    
    //Fill the vector of flags, they determine if a triangulation node is a CNT point (1) or a GNP point (0)
    //Initialize the vector with CNT flags
    points_out_flags.assign(n_cnt_points+n_points, 1);
    //Fill the GNP flags
    for (int i = n_cnt_points; i < n_cnt_points+n_points; i++) {
        points_out_flags[i] = 0;
    }
    
    
    return 1;
}
//This functions directly calculates generates a triangulation for the trivial cases (3 points or less)
int Triangulation::Generate_trivial_triangulation(const vector<long int> &points_out, const vector<short int> &points_out_flags, GCH &hybrid)
{
    //Get the number of points
    int n_case = (int)points_out.size();
    
    //Empty vectors
    vector<long int> empty_long;
    vector<short int> empty_short;
    
    //Check the trivial case
    if (n_case == 3) {
        
        //Generate a triangulation with three edges
        hybrid.triangulation.push_back(empty_long);
        hybrid.triangulation.back().push_back(points_out[0]);
        hybrid.triangulation.back().push_back(points_out[1]);
        hybrid.triangulation_flags.push_back(empty_short);
        hybrid.triangulation_flags.back().push_back(points_out_flags[0]);
        hybrid.triangulation_flags.back().push_back(points_out_flags[1]);
        
        hybrid.triangulation.push_back(empty_long);
        hybrid.triangulation.back().push_back(points_out[2]);
        hybrid.triangulation.back().push_back(points_out[1]);
        hybrid.triangulation_flags.push_back(empty_short);
        hybrid.triangulation_flags.back().push_back(points_out_flags[2]);
        hybrid.triangulation_flags.back().push_back(points_out_flags[1]);
        
        hybrid.triangulation.push_back(empty_long);
        hybrid.triangulation.back().push_back(points_out[0]);
        hybrid.triangulation.back().push_back(points_out[2]);
        hybrid.triangulation_flags.push_back(empty_short);
        hybrid.triangulation_flags.back().push_back(points_out_flags[0]);
        hybrid.triangulation_flags.back().push_back(points_out_flags[2]);

    }
    else if (n_case == 2) {
        
        //Generate a triangulation with one edge
        hybrid.triangulation.push_back(empty_long);
        hybrid.triangulation.back().push_back(points_out[0]);
        hybrid.triangulation.back().push_back(points_out[1]);
        hybrid.triangulation_flags.push_back(empty_short);
        hybrid.triangulation_flags.back().push_back(points_out_flags[0]);
        hybrid.triangulation_flags.back().push_back(points_out_flags[1]);
        
    }
    
    //If there is one point or no points then there is no need to make a triangulation
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
    vertex.y = 0.5*(ly + sqrt(3)*lx) + 1; //x is (still) zero
    //Rotate vertex (the center of the GNP is the origin of coordinates)
    vertex_rot = vertex.rotation(hybrid.rotation, hybrid.center);
    //Add to vector of vertices
    vertices.push_back(vertex_rot);
    
    //Base vertex B
    vertex.x = -0.5*leq - 1;
    vertex.y = -0.5*ly - 1;
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
int Triangulation::Bowyer_watson(const vector<Point_3D> &points_t_3d, vector<Point_3D> &vertices_s, const vector<long int> &points_t, const vector<short int> &points_t_flags, const vector<int> &vertices_t, vector<vector<int> > &triangles, GCH &hybrid)
{
    //Add each point to the triangulation
    for (int i = 0; i < (int)vertices_t.size(); i++) {
        //Vector for bad triangles
        vector<vector<int> > bad_triangles_edges;
        
        //Find the triangles whose circumcircle contains the current point and remove them from the triangulation
        if (!Find_bad_triangles(points_t_3d, vertices_s, vertices_t[i], triangles, bad_triangles_edges)) {
            hout << "Error in Bowyer_watson when calling Find_bad_triangles" << endl;
            return 0;
        }
        
        //Find the edges that will form new triangles and add them to the triangulation
        if (!Add_new_triangles(vertices_t[i], triangles, bad_triangles_edges)) {
            hout << "Error in Bowyer_watson when calling Add_triangle_as_edges" << endl;
            return 0;
        }
        
    }
    
    //Triangulation is finished, so find the edges that make up the final triangulation:
    // -Remove those edges that contain a vertex of the supertriangle
    // -Remove repeated edges
    if ( !Final_triangulation_edges(triangles, points_t, points_t_flags, vertices_t, hybrid) ) {
        hout << "Error in Bowyer_watson when calling Final_triangulation_edges" <<endl;
        return 0;
    }
    
    return 1;
}
//This function finds the bad triangles in a triangulation, i.e., the triangles that contain a point in their circumcircle
//The bad triangles are removed from the triangulation and added as edges
int Triangulation::Find_bad_triangles(const vector<Point_3D> &points_t, vector<Point_3D> &vertices_s, const int &point, vector<vector<int> > &triangles, vector<vector<int> > &bad_triangles_edges)
{
    //Scan all triangles
    for (int j = (int)triangles.size()-1; j >= 0 ; j--) {
        //If point is inside the circumcircle then add the triangle to bad_triangles
        if (Is_in_circumcircle(point, points_t, vertices_s, triangles[j])) {
            //Add current triangle to bad_triangles as edges
            Add_triangle_as_edges(triangles[j], bad_triangles_edges);
            //Remove triangle from triangulation
            triangles.erase(triangles.begin()+j);
        }
    }
    
    return 1;
}
//This function checks if a point is inside the circumcircle of a triangle
int Triangulation::Is_in_circumcircle(const int &point, const vector<Point_3D> &points_t, vector<Point_3D> &vertices_s, const vector<int> &triangle)
{
    //Calculate circumcenter of triangle
    Point_3D C = Calculate_circumcenter(points_t, vertices_s, triangle);
    
    //Calculate squared radius of circumcircle (i.e. squared distance from center to any point)
    //Comparing squared distances will save time by avoiding calculating squared roots
    Point_3D P = Get_point(triangle[0], points_t, vertices_s);
    double rad2 = C.squared_distance_to(P);
    
    //Squared distance from circumcenter to point
    P = Get_point(point, points_t, vertices_s);
    double dist = C.squared_distance_to(P);
    
    //If the squared distance between the point and C is smaller or equal than the squared radius,
    //then the point is inside the circumcircle so return 1
    return ( dist - rad2 <= Zero);
}
//Given a triangle, this function calculates its circumcenter
Point_3D Triangulation::Calculate_circumcenter(const vector<Point_3D> &points_t, vector<Point_3D> &vertices_s, const vector<int> &triangle)
{
    //Variables for calculations
    Point_3D P, Q, R;
    
    //Triangle points
    Point_3D P1 = Get_point(triangle[0], points_t, vertices_s);
    Point_3D P2 = Get_point(triangle[1], points_t, vertices_s);
    Point_3D P3 = Get_point(triangle[2], points_t, vertices_s);
    
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
//There are some negative indices to refer to the vertices vectors while all positive indices refer to the point_list vector
//so this function deals with that and returns the proper point
Point_3D Triangulation::Get_point(const int &index, const vector<Point_3D> &points_t, vector<Point_3D> &vertices)
{
    //If the index is equal or greater to zero, then the index remains the same
    if (index >= 0)
        return points_t[index];
    //If the index is negative, add one and change the sign to get the correct index in vertices
    else
        return vertices[-1-index];
}
//This function adds a triangle to the vector edges as the three edges that make up the triangle
void Triangulation::Add_triangle_as_edges(const vector<int> &triangle, vector<vector<int> > &edges)
{
    //Empty vector to add elements to edges
    vector<int> empty;
    
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
//then the new triangles are added to the triangulation with the edges found and a given point
int Triangulation::Add_new_triangles(const int &point, vector<vector<int> > &triangles, vector<vector<int> > &bad_triangles_edges)
{
    //Vector that counts the repetitions of each edge
    vector<int> repetitions_count(bad_triangles_edges.size(), 0);
    
    //Scan each edge in bad triangles and compare it with the rest
    for (int j = 0; j < (int)bad_triangles_edges.size(); j++) {
        
        //Scan all edges after edge j
        for (int k=j+1; k < (int)bad_triangles_edges.size(); k++) {
            
            //Check if the two edges are the same
            if (Is_same_edge(bad_triangles_edges[j], bad_triangles_edges[k])) {
                
                //If there are the same edge, increase the count for edges j and k
                repetitions_count[j] = repetitions_count[j] + 1;
                repetitions_count[k] = repetitions_count[k] + 1;
                
            }
        }
    }
    
    //Scan the vector of repetitions_count and add a new triangle if an edge has 0 repetitions
    for (int i = 0; i < (int)repetitions_count.size(); i++) {
        
        //Check if the count is zero
        if (!repetitions_count[i]) {
            
            //Add a new triangle using the two points of the edge and the current point
            triangles.push_back(bad_triangles_edges[i]);
            triangles.back().push_back(point);
        }
    }
    
    return 1;
}
//This function compares two edges and decides if they are equal or not
int Triangulation::Is_same_edge(const vector<int> &edge1, const vector<int> &edge2)
{
    return ( (edge1[0]==edge2[0] && edge1[1]==edge2[1]) || (edge1[0]==edge2[1] && edge1[1]==edge2[0]) );
}
//This function removes the edges that have vertices of the super triangle and
//generates the vector of edges that make up the triangulation
int Triangulation::Final_triangulation_edges(const vector<vector<int> > &triangles, const vector<long int> &points_t, const vector<short int> &points_t_flags, const vector<int> &vertices_t, GCH &hybrid)
{
    //Vector to store the triangulation edges using the consecutive numbering for vertices
    vector<vector<int> > edges;
    
    //Scan all triangles looking for valid edges
    for (int i = 0; i < (int)triangles.size(); i++) {
        //Check if there are any valid edges and add them to the vector of edges
        if (!Valid_edges(triangles[i], edges)) {
            hout << "Error in Final_triangulation_edges" << endl;
            return 0;
        }
    }
    /*/
    for (int i = 0; i < (int)edges.size(); i++) {
        hout << "edge" << i << ": ";
        for (int j = 0; j < (int)edges[i].size(); j++) {
            hout << edges[i][j] << ' ';
        }
        hout << endl;
    }//*/
    
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
    /*/
    hout << "Delete repeated edges" << endl;
    for (int i = 0; i < (int)edges.size(); i++) {
        hout << "edge" << i << ": ";
        for (int j = 0; j < (int)edges[i].size(); j++) {
            hout << edges[i][j] << ' ';
        }
        hout << endl;
    }//*/
    
    //Add the edges and flags to the hybrid using the point numbers
    
    //Empty vectors
    vector<long int> empty_long;
    vector<short int> empty_short;
    //Add the triangulation points and flags to the hybrid
    for (int i = 0; i < (int)edges.size(); i++) {
        
        //Get vertices
        int V1 = edges[i][0];
        int V2 = edges[i][1];
        
        //Get point numbers
        long int P1 = points_t[V1];
        long int P2 = points_t[V2];
        
        //Add triangulation edge to hybrid
        hybrid.triangulation.push_back(empty_long);
        hybrid.triangulation.back().push_back(P1);
        hybrid.triangulation.back().push_back(P2);
        
        //Get flags
        short int F1 = points_t_flags[V1];
        short int F2 = points_t_flags[V2];
        
        //Add triangulation flags to hybrid
        hybrid.triangulation_flags.push_back(empty_short);
        hybrid.triangulation_flags.back().push_back(F1);
        hybrid.triangulation_flags.back().push_back(F2);
    }
    
    return 1;
}
//Make edges using only valid edges, i.e., excluding the vertices of the supertriangle
int Triangulation::Valid_edges(const vector<int> &triangle, vector<vector<int> > &edges)
{
    //Temporary vector
    vector<int> edge_tmp;
    
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