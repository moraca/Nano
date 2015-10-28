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
    hout<<"Observation window geometry:"<<endl;
    hout<<"xmin="<<xmin<<" ymin="<<ymin<<" zmin="<<zmin<<endl;
    hout<<"w_x="<<w_x<<" w_y="<<w_y<<" w_z="<<w_z<<endl;
    
    //Vector check for debugging. Comment or delete after debugging
    //vector<vector<long int> > structure_check(structure);
    //Print the points
    //Printer *P = new Printer;
   // P->Print_1d_vec(points_in, "points_in.txt");

    //hout << "5 ";
    //Scan every Nanotube that is the boundary region. Delete and trim CNTs when needed.
    if (!Trim_boundary_cnts(shells_cnt, window, sample, points_in, structure, radii)){
        hout << "Error in Locate_and_trim_boundary_cnts (initial)" << endl;
        return 0;
    }
    
    //hout << "6 ";
    //Fill the vector cnts_inside
    int flag = 0;
    for (int i = 0; i < (int)structure.size(); i++) {
        if (structure[i].size()) {
            if (structure[i].size() == 1) {
                hout << "Error in Extract_observation_window. A CNT ("<<i<<") has only one point. A CNT must have at least 2 points."<<endl;
                long int P = structure[i][0];
                hout << "\tP=("<<points_in[P].x<<", "<<points_in[P].y<<", "<<points_in[P].z<<")"<<endl;
                //hout << "\tThere were "<< structure_check.size()<<" CNTs before trimming. ";
                //hout << "CNT "<<i<<" had "<<structure_check[i].size()<<" points before trimming"<<endl;
                flag = 1;
            } else {
                cnts_inside.push_back(i);
            }
        }
    }
    if (flag) {
        return 0;
    }
    //hout << "9 ";
	return 1;
}

int Cutoff_Wins::Trim_boundary_cnts(vector<vector<int> > &shells_cnt, int window, struct Geom_RVE sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<double> &radii)
{
    //These variables will help me locate the point with respect with the box
    string currentPoint;
    //Variables for current and next points and current CNT
    long int P1;
    int CNT;
    //Empty vector to increase size of other vectors
    //Initialize the vector of boundary_flags with empty vectors
    vector<short int> empty_short;
    boundary_flags.assign(points_in.size(), empty_short);
    //Initialize the vector of boundary_flags with empty vectors
    vector<int> empty_int;
    boundary_cnt.assign(6, empty_int);
    //hout << "cnts_inside.size() = "<<cnts_inside.size()<<endl;
    //This varibale is used to initialize the vectors below
    vector<long int> empty_long;
    for (long int i = 0; i < (long int)shells_cnt[window].size(); i++) {
        CNT = shells_cnt[window][i];
        //hout << "Check0 " <<  CNT << ' ' << structure[CNT].size() << ' ' << endl;
        //hout << " all=" << structure.size() << ' ' ;
        //hout << "Check0.1 " ;
        
        //Vector to store the indices of the segments
        vector<int> branches_indices_CNT;
        
        //Scan all points in a given CNT
        for (int j = 0; j < (int)structure[CNT].size(); j++) {
            P1 = structure[CNT][j];
            currentPoint = Where_is(points_in[P1]);
            if ((currentPoint == "inside") ) {
                //If a point is "inside", save the local point number only if the current size of branches_indices[CNT] is even
                //When branches_indices[CNT].size() is even, that means one of two cases:
                //1) It is empty, so it is the first time it reaches this part, so a new segement needs to be started.
                //   Thus, add the local point number. In the following iterations, if the points are still inside,
                //   no point will be added since now the size of the vector is odd, i.e. 1
                //2) It is not empty and its size was made even by entering the "outside case", then it went back to the
                //   "inside" case. Hence, a new segment needs to be added. In the following iterations, if the points are still inside,
                //   no point will be added since now the size of the vector is odd.
                if (branches_indices_CNT.size()%2 == 0) {
                    branches_indices_CNT.push_back(j);
                }
            } else {
                if (branches_indices_CNT.size()%2 == 1) {
                    branches_indices_CNT.push_back(j);
                }
            }
        }
        
        //If the last two segments of a CNT are "outside"-"inside", then the vector branches_indices[CNT] will have odd size
        //in such case, the last point needs to be added
        if (branches_indices_CNT.size()%2 == 1) {
            int s = (int)structure[CNT].size()-1;
            branches_indices_CNT.push_back(s);
        }
        
        //After the segments have been defined, it is time to trim the CNT
        if (branches_indices_CNT.size() == 0) {
            //If there are no indices, that means all the CNT is outside, so delete all points of that CNT
            structure[CNT].clear();
        } else {
            //Scan each segment to handle boundary points
            for (int k = 0; k < (int)branches_indices_CNT.size(); k=k+2) {
                int index1 = branches_indices_CNT[k];
                int index2 = branches_indices_CNT[k+1];
                //Beginning of segement
                //Check if the first index is the first point of the CNT, otherwise add a boundary point
                if (index1 != 0) {
                    if (!First_index(points_in, structure[CNT], index1)){
                        hout << "Error in Trim_boundary_cnts2. branches_indices_CNT["<<k<<"]="<<branches_indices_CNT[k];
                        return 0;
                    }
                    //branches_indices[CNT][k] might be modified so I need to update it
                    branches_indices_CNT[k] = index1;
                }
                //End of segment
                //Check if the second index is the last point of the CNT, otherwise add a boundary point
                if (index2 != (int)structure[CNT].size()-1) {
                    if(!Second_index(points_in, structure[CNT], index2)){
                        hout << "Error in Trim_boundary_cnts2. branches_indices_CNT["<<k+1<<"]="<<branches_indices_CNT[k+1];
                        return 0;
                    }
                    //branches_indices[CNT][k+1] might be modified so I need to update it
                    branches_indices_CNT[k+1] = index2;
                }
                if (index1 == index2) {
                    hout << "Error in Trim_boundary_cnts2. index1 = index2 = "<<index1<<" on CNT "<<CNT<<endl;
                    hout << "points_in[structure[CNT][index2-1]] is "<<Where_is(points_in[structure[CNT][index2-1]]);
                    hout << " points_in[structure[CNT][index2]] is "<<Where_is(points_in[structure[CNT][index2]]) << endl;
                    return 0;
                }
                //If there are more than one segments, add the extra segments as new CNTs
                //A minimum of two segments means that k will have values 0 and 2, so whenever k is 2 or more there are multiple segments
                if (k >=2) {
                    //Add a new CNT
                    structure.push_back(empty_long);
                    
                    //This bg variable is used to add the new CNT into the corresponding shell-sub-region
                    Background_vectors *bg = new Background_vectors;
                    
                    //Add the points of the segment to the new CNT
                    //At the same time, update the CNT number of the points in the new CNT and add the CNT to the corresponding shell or shells
                    int new_CNT = ((int)structure.size()) - 1;
                    for (int kk = index1; kk <= index2; kk++) {
                        long int P = structure[CNT][kk];
                        structure.back().push_back(P);
                        points_in[P].flag = new_CNT;
                        bg->Add_to_shell(sample, points_in[P], shells_cnt);
                    }
                    
                    //Update the radii vector
                    //The new CNT is just a segment of the old one, so they should have the same radius
                    radii.push_back(radii[CNT]);
                }
            }
            //At this point all indices are inclusive of the boundary points, and these boundary points have been added into the
            //points_in vector by substituting outside points.
            
            //Move the points that are inside to the front of the CNT
            //If index1 is zero, the points of the first segment are already at the front of the CNT, so there is nothing to do
            int index1 = branches_indices_CNT.front();
            int index2 = branches_indices_CNT[1];
            if (index1 != 0) {
                for (int kk = index1; kk <= index2 ; kk++) {
                    structure[CNT][kk-index1] = structure[CNT][kk];
                }
            }
            //Remove the points that are outside or belong to other CNTs
            while ((int)structure[CNT].size() > (index2-index1+1)) {
                structure[CNT].pop_back();
            }
        }
        
    }

    return 1;
}

int Cutoff_Wins::First_index(vector<Point_3D> &points_in, vector<long int> &structure_CNT, int &index1)
{
    long int global_i = structure_CNT[index1];
    long int global_o = global_i-1;

    //Check if the outside point is in the boundary. This actually happened in some simulations so it is useful to check
    //So if the outside point is actually at the boundary, there is nothing to do.
    //Only when the outside point is not at the boundary, then we proceed to calculate the projection to the boundary
    if ( Where_is(points_in[global_o]) != "boundary") {
        if (!Substitute_boundary_point(points_in, global_i, global_o)){
            hout << "Error in First_index. global_i="<<global_i<<" global_o="<<global_o<<" structure_CNT.size()="<<structure_CNT.size();
            hout <<" index1="<<index1<<endl;
            hout <<"\tP_i=("<<points_in[global_i].x<<", "<<points_in[global_i].y<<", "<<points_in[global_i].z<<") P_o=(";
            hout <<points_in[global_o].x<<", "<<points_in[global_o].y<<", "<<points_in[global_o].z<<")"<<endl;
            return 0;
        }
        //Now, the outside poin has the coordinates of the boundary point
        Add_to_boundary_vectors(points_in[global_o], global_o);
    }
    //The first index is always inside, so:
    //    if the previous point is outside, the previous point needs to be included as it will bow be the boundary point
    //    if the previous point is boundary, it has to be added
    //Hence, independently of where the previous point is, it has to be included, so it will be the new index
    index1--;
    return 1;
}

int Cutoff_Wins::Second_index(vector<Point_3D> &points_in, vector<long int> &structure_CNT, int &index2)
{
    //The second index can be either inside or outside, so first I need to find out where it is and proceed accordingly
    if (Where_is(points_in[ structure_CNT[index2] ]) == "outside"){
        long int global_o = structure_CNT[index2];
        long int global_i = global_o-1;
        
        //Check if what is supposed to be the inside point is actually in the boundary.
        //If it hapens that the inside point is actually at the boundary, then this point has to be index 2
        if ( Where_is(points_in[global_i]) == "boundary") {
            index2--;
        } else{
            //Since the index2 point is outside, then the previous point is inside.
            //We then we proceed to calculate the projection to the boundary
            if (!Substitute_boundary_point(points_in, global_i, global_o)){
                hout << "Error in Second_index. global_i="<<global_i<<" global_o="<<global_o<<" structure_CNT.size()="<<structure_CNT.size();
                hout <<" index2="<<index2<<endl;
                hout <<"\tP_i=("<<points_in[global_i].x<<", "<<points_in[global_i].y<<", "<<points_in[global_i].z<<") P_o=(";
                hout <<points_in[global_o].x<<", "<<points_in[global_o].y<<", "<<points_in[global_o].z<<")"<<endl;
                return 0;
            }
            //Now, the outside point, i.e. index2, has the coordinates of the boundary point so I need to kep it unchanged
            Add_to_boundary_vectors(points_in[global_o], global_o);
        }
    }
    return 1;
}

//This function checks in which of these three location a point is placed:
//outside the observation window
//inside the observation window
//at the boundary of the observation window
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

//Add the corrent point to the corrsponding boundary vector.
//The boundary vectors are used in the direct electrifying algorithm to find the nodes with known boundary conditions
void Cutoff_Wins::Add_to_boundary_vectors(Point_3D point3d, long int point)
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