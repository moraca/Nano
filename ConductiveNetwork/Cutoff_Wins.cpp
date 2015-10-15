//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Cutoff_Wins.cpp
//OBJECTIVE:	To cutoff the windows
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Cutoff_Wins.h"

/*
 Input:
    vector<vector<int> > shells_cnt
        List with all CNTs grouped into sub-regions in order to reduce computational cost of finding the CNTs intersected by the boundaries of the observation window
    int window
        Current sub-region that contains the CNTs to be trimmed
    struct Geom_RVE sample
        Geometry of the generated sample
    struct Nanotube_Geo cnts
        Geometry of the CNTs
    vector<vector<long int> > structure
        Vector with the structure. Since CNTs might get trimmed, this vector will be modified to delete points and add more CNTs and points.
    vector<double> radii
        List of radii. Using this vector allows for the code to be able to work with CNTs of different radii. As CNTs might be added, the radii of the new CNTs need also to be included in the vector
    vector<Point_3D> points_in
        Vector with the point coordinates. Since points at the boundaries of the observation window are going to be added, this vector needs to be modified
 
 Output (These three are class variables):
    vector<int> cnts_inside
        List with the CNTs that are inside the observation window and that will be used for future computations
    vector<vector<int> > boundary_cnt
        Vectors that contains the CNTs on each of the six boundaries. These are needed to determine percolation
            boundary_cnt[0] corresponds to x0
            boundary_cnt[1] corresponds to x1
            boundary_cnt[2] corresponds to y0
            boundary_cnt[3] corresponds to y1
            boundary_cnt[4] corresponds to z0
            boundary_cnt[5] corresponds to z1
 
 Modified inputs:
    vector<vector<long int> > structure
    vector<double> radii
    vector<Point_3D> points_in
 
 */

//This function removes the points that are outside the observation window.
//The vector cnts_inside is created so only the CNTs inside the obseration window are considered in other functions
int Cutoff_Wins::Extract_observation_window(struct Geom_RVE sample, struct Nanotube_Geo cnts, vector<vector<long int> > &structure, vector<double> &radii, vector<Point_3D> &points_in, vector<vector<int> > &shells_cnt, int window)
{
    //These are variables for the geometry of the observation window
    //Dimensions of the current observation window
    w_x = sample.win_max_x - window*sample.win_delt_x;
    w_y = sample.win_max_y - window*sample.win_delt_y;
    w_z = sample.win_max_z - window*sample.win_delt_z;
    //These variables are the coordinates of the lower corner of the observation window
    xmin = sample.origin.x + (sample.len_x - w_x)/2;
    ymin = sample.origin.y + (sample.wid_y - w_y)/2;
    zmin = sample.origin.z + (sample.hei_z - w_z)/2;
    hout<<"xmin="<<xmin<<" ymin="<<ymin<<" zmin="<<zmin<<endl;
    hout<<"w_x="<<w_x<<" w_y="<<w_y<<" w_z="<<w_z<<endl;

    //hout << "5 ";
    //Scan every Nanotube that is the boundary region. Delete and trim CNTs when needed.
    if (!Trim_boundary_cnts(shells_cnt, window, sample, points_in, structure, radii)){
        hout << "Error in Locate_and_trim_boundary_cnts (initial)" << endl;
        return 0;
    }
    
    //hout << "6 ";
    //Fill the vector cnts_inside
    for (int i = 0; i < (int)structure.size(); i++) {
        if (structure[i].size())
            cnts_inside.push_back(i);
    }
    //hout << "9 ";
	return 1;
}

int Cutoff_Wins::Trim_boundary_cnts(vector<vector<int> > &shells_cnt, int window, struct Geom_RVE sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<double> &radii)
{
    //These variables will help me locate the point with respect with the box
    string currentPoint, nextPoint;
    //Variables for current and next points and current CNT
    long int P1, P2;
    int CNT;
    //Empty vector to increase size of other vectors
    //Initialize the vector of boundary_flags with empty vectors
    vector<short int> empty_short;
    boundary_flags.assign(points_in.size(), empty_short);
    //Initialize the vector of boundary_flags with empty vectors
    vector<int> empty_int;
    boundary_cnt.assign(6, empty_int);
    //hout << "cnts_inside.size() = "<<cnts_inside.size()<<endl;
    for (long int i = 0; i < (long int)shells_cnt[window].size(); i++) {
        CNT = shells_cnt[window][i];
        //hout << "Check0 " <<  CNT << ' ' << structure[CNT].size() << ' ' << endl;
        //hout << " all=" << structure.size() << ' ' ;
        P1 = structure[CNT][0];
        currentPoint = Where_is(points_in[P1]);
        //hout << "Check0.1 " ;
        for (long int j = 1; j <= structure[CNT].size(); j++) {
            //hout << "Check1 " << i << ' ' << j << " CNT=" << CNT << ' ' << "P1=" << P1 << " s[CNT].s=" << structure[CNT].size() << ' ' << currentPoint << ' ';
            //Handle the last point
            if (j == structure[CNT].size()) {
                nextPoint = "nothing";
                //hout << " nothing "<<endl;
            } else {
                P2 = structure[CNT][j];
                nextPoint = Where_is(points_in[P2]);
                //hout << nextPoint <<endl;
            }
            //hout << " (" << points_in[P1].x << ", " << points_in[P1].y << ", " << points_in[P1].z << ") ";
            //hout << "Check2 " << endl;//*/
            
            //Check where is the current point and where the next point
            //First I check if the point is inside the box. If it is, what follows makes that assumption
            if (currentPoint == "inside") {
                if (nextPoint == "outside") {
                    //hout << "Check3 ";
                    //If the next point is outside, then it will be substituted by the projection at the boundary.
                    //Make the projection and add the new point to the point_list and the the structure vector
                    //j-1 is the current point
                    //j is the next point
                    if (!Substitute_boundary_point(points_in, structure[CNT][j-1], structure[CNT][j])) {
                        hout  << "Error in Locate_and_trim_boundary_cnts. currentPoint=" << currentPoint << " nextPoint=" << nextPoint << endl;
                        return 0;
                    }
                    //Trim the CNT from the projected boundary point, which now is in position j
                    Trim_CNT(shells_cnt, sample, points_in, structure, radii, j, CNT);
                    //Now the position of nextPoint is for the boundary point, so I need to update nextPoint
                    nextPoint = "boundary";
                    //hout << "Check5 ";
                }
            } else if (currentPoint == "outside") {
                if (nextPoint == "inside") {
                    //hout << "Check6 ";
                    //When the next point is inside, take the intersection and add the new point
                    //to the point_list and to the structure vector
                    if (!Substitute_boundary_point(points_in, structure[CNT][j], structure[CNT][j-1])) {
                        hout  << "Error in Locate_and_trim_boundary_cnts. currentPoint=" << currentPoint << " nextPoint=" << nextPoint << endl;
                        return 0;
                    }
                    //Now the position of currentPoint is for "boundary" and nextPoint remains as "inside", so I do not need to update positions
                    //hout << "Check7 ";
                } else {
                    //hout << "Check11 ";
                    //Remove the current point from the structure vector
                    structure[CNT].erase(structure[CNT].begin()+j-1);
                    //Update the iterator j. This is needed beacuase one element was deleted, so all positions shift one place
                    j--;
                }
                //hout << "Check12 ";
            } else if (currentPoint == "boundary") {
                //hout << "Check13 ";
                //If the next point is outside, then the CNT might need to be trimmed
                if (nextPoint == "outside") {
                    //hout << "Check14 ";
                    //If the boundary point is the first point of the CNT, then this point needs to be deleted and
                    //treated as if it was outside
                    if (j==1){
                        //hout << "Check14.1 ";
                        structure[CNT].erase(structure[CNT].begin());
                        j--;
                    } else {
                        //hout << "Check14.2 ";
                        //If the boundary is not the first point, then a segment of the CNT is inside the observation window so trim the CNT
                        Trim_CNT(shells_cnt, sample, points_in, structure, radii, j-1, CNT);
                        //Add the current point to the corresponding boundary vector
                        long int point_number = structure[CNT][j-1];
                        Add_to_boundary_vectors(sample, points_in[point_number], point_number);
                    }
                } else if ((nextPoint == "nothing")&&(j==1)) {
                    //hout << "Check13.1 ";
                    //If this is the last and only point, then delete it
                    structure[CNT].erase(structure[CNT].begin());
                    j--;
                } else {
                    //hout << "Check13.2 ";
                    //Add the current point to the corresponding boundary vector
                    long int point_number = structure[CNT][j-1];
                    Add_to_boundary_vectors(sample, points_in[point_number], point_number);
                }
                //hout << "Check15 ";
            } else {
                hout << "Current point has an invalid position. The only posibilities are: inside, outside, ";
                hout << "boundary or nothing. Current location is: " << currentPoint << endl;
                return 0;
            }
            //update current point
            currentPoint = nextPoint;
            P1 = P2;
            //hout << "Check16 " << endl;
        }
        //hout << "Check17 " << endl;
    }
    
    return 1;
}

string Cutoff_Wins::Where_is(Point_3D point)
{
    double x = point.x;
    double y = point.y;
    double z = point.z;
    //If the point is outside the observation window then any these conditions needs to be true
    if ((x < xmin)||(x > xmin+w_x)||(y < ymin)||(y > ymin+w_y)||(z < zmin)||(z > zmin+w_z))
        return "outside";
    //If the point it's at a boundary of the observation window, then only one of these conditions needs to be true,
    //provided that it is not outside
    else if ((abs(x - xmin) < Zero)||(abs(x - (xmin+w_x)) < Zero)||(abs(y - ymin) < Zero)||(abs(y - (ymin+w_y)) < Zero)||(abs(z - zmin) < Zero)||(abs(z - (zmin+w_z)) < Zero))
        return "boundary";
    //If the point is not outside the observation window nor at the boundary, then it's inside
    else
        return "inside";//*/
}

//When two consecutive points are found to be one outside and one inside, this function substitutes the outside point by the intersection
//of the segment between the two points with the boundaries of the observation window
int Cutoff_Wins::Substitute_boundary_point(vector<Point_3D> &points_in, long int global_i, long int global_o)
{
    vector<Point_3D> ipoi_vec;
    //The point at the boundary is in ipoi_vec[0]
    //hout<<"global_i="<<global_i<<" global_o="<<global_o<<endl;
    //hout<<"points_in[global_i]=("<<points_in[global_i].x<<','<<points_in[global_i].y<<','<<points_in[global_i].z<<')'<<endl;
    //hout<<"points_in[global_o]=("<<points_in[global_o].x<<','<<points_in[global_o].y<<','<<points_in[global_o].z<<')'<<endl;
    //hout<<"xmin="<<xmin<<" ymin="<<ymin<<" zmin="<<zmin<<endl;
    //hout<<"w_x="<<w_x<<" w_y="<<w_y<<" w_z="<<w_z<<endl;
    if (!Get_intersecting_point_RVE_surface(points_in[global_i], points_in[global_o], ipoi_vec)) {
        hout<< "Error in Substitute_boundary_point." <<endl;
        return 0;
    }
    //Substitute the outside point by the one at the boundary
    points_in[global_o].x = ipoi_vec[0].x;
    points_in[global_o].y = ipoi_vec[0].y;
    points_in[global_o].z = ipoi_vec[0].z;
        
    return 1;
}

int Cutoff_Wins::Get_intersecting_point_RVE_surface(Point_3D &point0, Point_3D &point1, vector<Point_3D> &ipoi_vec)
{
    double t_temp[6];
    //The planes (surfaces of RVE) perpendicular to X axis
    t_temp[0] = (xmin - point0.x)/(point1.x - point0.x);
    t_temp[1] = (xmin + w_x - point0.x)/(point1.x - point0.x);
    //The planes (surfaces of RVE) perpendicular to Y axis
    t_temp[2] = (ymin - point0.y)/(point1.y - point0.y);
    t_temp[3] = (ymin + w_y - point0.y)/(point1.y - point0.y);
    //The planes (surfaces of RVE) perpendicular to Z axis
    t_temp[4] = (zmin - point0.z)/(point1.z - point0.z);
    t_temp[5] = (zmin + w_z - point0.z)/(point1.z - point0.z);
    
    vector<double> t_ratio;
    for(int i=0; i<6; i++)
    {
        if(t_temp[i]>=0&&t_temp[i]<1)
        {
            //Binary insertion sort
            int left = 0;
            int right = (int)t_ratio.size()-1;
            while(right>=left)
            {
                int middle = (left + right)/2;
                if(fabs(t_ratio[middle] - t_temp[i])<Zero) goto T_Value_Same; //the case with same values
                else if(t_ratio[middle] > t_temp[i]) right = middle - 1;
                else left = middle + 1;
            }
            t_ratio.insert(t_ratio.begin()+left, t_temp[i]);	//insertion
        T_Value_Same: ;
        }
    }
    
    if((int)t_ratio.size()<1||(int)t_ratio.size()>3)
    {
        hout << "Error, the number of intersection points between the segement and the surfaces of RVE is " << (int)t_ratio.size() << ", less than one or more than three!" << endl;
        return 0;
    }
    
    Point_3D point_temp;
    for(int i=0; i<(int)t_ratio.size(); i++)
    {
        point_temp.x = point0.x+(point1.x-point0.x)*t_ratio[i];
        point_temp.y = point0.y+(point1.y-point0.y)*t_ratio[i];
        point_temp.z = point0.z+(point1.z-point0.z)*t_ratio[i];
        point_temp.flag = 1;		//a temporary point
        
        //---------------------------------------------------------------------------
        //Error correction
        if(fabs(point_temp.x-xmin)<Zero) point_temp.x = xmin;
        else if(fabs(point_temp.x-xmin-w_x)<Zero) point_temp.x = xmin + w_x;
        
        if(fabs(point_temp.y-ymin)<Zero) point_temp.y = ymin;
        else if(fabs(point_temp.y-ymin-w_y)<Zero) point_temp.y = ymin + w_y;
        
        if(fabs(point_temp.z-zmin)<Zero) point_temp.z = zmin;
        else if(fabs(point_temp.z-zmin-w_z)<Zero) point_temp.z = zmin + w_z;
        
        //---------------------------------------------------------------------------
        //Insert a new point
        ipoi_vec.push_back(point_temp);
    }
    
    return 1;
}

void Cutoff_Wins::Trim_CNT(vector<vector<int> > &shells_cnt, struct Geom_RVE sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<double> &radii, long int boundary, int CNT)
{
    //This bg variable is used to add the new CNT into the corresponding shell-sub-region
    Background_vectors *bg = new Background_vectors;
    
    //Here check which point is the current point
    vector<long int> empty;
    structure.push_back(empty);
    //The CNT will bre trimed from the point after the boundary point and up to the last point
    structure.back().insert(structure.back().begin(),structure[CNT].begin()+boundary+1, structure[CNT].end());
    //Erase form CNT
    structure[CNT].erase(structure[CNT].begin()+boundary+1, structure[CNT].end());
    
    //Update the CNT number of the points that were moved and add the new CNT to the corresponding shell or shells
    long int P;
    int new_CNT = ((int)structure.size()) - 1;
    for (long int i = 0; i < structure.back().size(); i++) {
        P = structure.back()[i];
        points_in[P].flag = new_CNT;
        bg->Add_to_shell(sample, points_in[P], shells_cnt);
    }
    
    //Update the radii vector
    //The new CNT is just a segment of the old one so they should have the same radius
    radii.push_back(radii[CNT]);
}

//Add the corrent point to the corrsponding boundary vector.
//The boundary vectors are used in the direct electrifying algorithm to find the nodes with known boundary conditions
void Cutoff_Wins::Add_to_boundary_vectors(struct Geom_RVE sample, Point_3D point3d, long int point)
{
    //Add point and CNT to the boundary vector
    double x = point3d.x;
    double y = point3d.y;
    double z = point3d.z;
   int CNT = point3d.flag;
    if ( abs(x - xmin) < Zero){
        Add_CNT_to_boundary(boundary_cnt[0], CNT, point,0,0);
    } else if ( abs(x - (xmin+w_x)) < Zero ){
        Add_CNT_to_boundary(boundary_cnt[1], CNT, point,0,1);
    } else if ( abs(y - ymin) < Zero ){
        Add_CNT_to_boundary(boundary_cnt[2], CNT, point,1,0);
    } else if ( abs(y - (ymin+w_y)) < Zero ){
        Add_CNT_to_boundary(boundary_cnt[3], CNT, point,1,1);
    } else if ( abs(z - zmin) < Zero ) {
        Add_CNT_to_boundary(boundary_cnt[4], CNT, point,2,0);
    } else if ( abs(z - (zmin+w_z)) < Zero ) {
        Add_CNT_to_boundary(boundary_cnt[5], CNT, point,2,1);
    }
}

//This function adds a CNT to the corresponding boundary vector.
//The falgs are used in the direct electrifying algorithm:
//flag1: indicates the direction 0 is x, 1 is y, 2 is z
//flag2: indicates which boundary 0 is for x0, y0 or z0; 1 is for x1, y1 or z1
void Cutoff_Wins::Add_CNT_to_boundary(vector<int> &boundary, int CNT, long int point, short int flag1, short int flag2)
{
    if (!boundary.size()) {
        boundary.push_back(CNT);
    } else if(boundary.back() != CNT){
        boundary.push_back(CNT);
    }
    boundary_flags[point].push_back(flag1);
    boundary_flags[point].push_back(flag2);
}