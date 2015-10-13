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
    vector<vector<int> > sectioned_domain_cnt
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
int Cutoff_Wins::Extract_observation_window(vector<vector<int> > &sectioned_domain_cnt, int window, struct Geom_RVE sample, struct Nanotube_Geo cnts, vector<vector<long int> > &structure, vector<double> &radii, vector<Point_3D> &points_in)
{    
    //hout << "6 ";
    //Scan every Nanotube that is the boundary region. Delete and trim CNTs when needed.
    if (!Trim_boundary_cnts(sectioned_domain_cnt, window, sample, points_in, structure, radii)){
        hout << "Error in Locate_and_trim_boundary_cnts (initial)" << endl;
        return 0;
    }
    
    //Fill the vector cnts_inside
    for (int i = 0; i < (int)structure.size(); i++) {
        if (structure[i].size())
            cnts_inside.push_back(i);
    }
    //hout << "9 ";
	return 1;
}

int Cutoff_Wins::Trim_boundary_cnts(vector<vector<int> > &sectioned_domain_cnt, int window, struct Geom_RVE sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<double> &radii)
{
    //These variables will help me locate the point with respect with the box
    string currentPoint, nextPoint;
    //Variables for current and next points and current CNT
    long int P1, P2;
    int CNT;
    //Empty vector to increase size of other vectors
    vector<short int> empty;
    //Resize the vector of boundary_flags with empty vector
    boundary_flags.assign(points_in.size(), empty);
    //hout << "cnts_inside.size() = "<<cnts_inside.size()<<endl;
    for (long int i = 0; i < sectioned_domain_cnt[window].size(); i++) {
        CNT = sectioned_domain_cnt[window][i];
        //hout << "Check0 " <<  CNT << ' ' << structure[CNT].size() << ' ' << endl;
        //hout << " all=" << structure.size() << ' ' ;
        P1 = structure[CNT][0];
        currentPoint = Where_is(sample, points_in[CNT]);
        //hout << "Check0.1 " ;
        for (long int j = 1; j <= structure[CNT].size(); j++) {
            //hout << "Check1 " << i << ' ' << j << " CNT=" << CNT << ' ' << "P1=" << P1 << " s[CNT].s=" << structure[CNT].size() << ' ' << currentPoint << ' ';
            //Handle the last point
            if (j == structure[CNT].size()) {
                nextPoint = "nothing";
                //hout << " nothing ";
            } else {
                P2 = structure[CNT][j];
                nextPoint = Where_is(sample, points_in[P2]);
                //hout << nextPoint << ' ';
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
                    if (!Substitute_boundary_point(sample, points_in, structure, j-1, j, CNT, currentPoint)) {
                        hout  << "Error in Locate_and_trim_boundary_cnts. currentPoint=" << currentPoint << " nextPoint=" << nextPoint << endl;
                        return 0;
                    }
                    //Trim the CNT from the projected boundary point, which now is in position j
                    Trim_CNT(sectioned_domain_cnt, sample, points_in, structure, radii, j, CNT);
                    //Now the position of nextPoint is for the boundary point, so I need to update nextPoint
                    nextPoint = "boundary";
                    //hout << "Check5 ";
                }
            } else if (currentPoint == "outside") {
                if (nextPoint == "inside") {
                    //hout << "Check6 ";
                    //When the next point is inside, take the intersection and add the new point
                    //to the point_list and to the structure vector
                    if (!Substitute_boundary_point(sample, points_in, structure, j, j-1, CNT, currentPoint)) {
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
                        structure[CNT].erase(structure[CNT].begin());
                        j--;
                    } else {
                        //If the boundary is not the first point, then a segment of the CNT is inside the observation window so trim the CNT
                        Trim_CNT(sectioned_domain_cnt, sample, points_in, structure, radii, j-1, CNT);
                        //Add the current point to the corresponding boundary vector
                        long int point_number = structure[CNT][j-1];
                        Add_to_boundary_vectors(sample, points_in[point_number], point_number);
                    }
                } else if ((nextPoint == "nothing")&&(j==1)) {
                    //If this is the last and only point, then delete it
                    structure[CNT].erase(structure[CNT].begin());
                    j--;
                } else {
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
    }
    
    return 1;
}

string Cutoff_Wins::Where_is(struct Geom_RVE sample, Point_3D point)
{
    double x = point.x;
    double y = point.y;
    double z = point.z;
    double xmin = sample.origin.x;
    double ymin = sample.origin.y;
    double zmin = sample.origin.z;
    double lx = sample.len_x;
    double ly = sample.wid_y;
    double lz = sample.hei_z;
    //If the point is outside the box then any these conditions needs to be true
    if ((x < xmin)||(x > xmin+lx)||(y < ymin)||(y > ymin+ly)||(z < zmin)||(z > zmin+lz))
        return "outside";
    //If the point it's at a boundary of the box, then only one of these conditions needs to be true,
    //provided that it is not outside
    else if ((abs(x - xmin) < Zero)||(abs(x - (xmin+lx)) < Zero)||(abs(y - ymin) < Zero)||(abs(y - (ymin+ly)) < Zero)||(abs(z - zmin) < Zero)||(abs(z - (zmin+lz)) < Zero))
        return "boundary";
    //If the point is not outside nor at the boundary, then it's inside
    else
        return "inside";//*/
}

//When two consecutive points are found to be one outside and one inside, this function substitutes the outside point by the intersection
//of the segment between the two points with the boundaries of the observation window
int Cutoff_Wins::Substitute_boundary_point(struct Geom_RVE sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, long int insidePoint, long int outsidePoint, int CNT, string currentLocation)
{
    //To make things easier, the inside point will be x,y,z and the outside point x1,y1,z1
    //Global coordinate
    long int global_i = structure[CNT][insidePoint];
    long int global_o = structure[CNT][outsidePoint];
    
    vector<Point_3D> ipoi_vec;
    //The point at the boundary is in ipoi_vec[0]
    if (!Get_intersecting_point_RVE_surface(sample, points_in[global_i], points_in[global_o], ipoi_vec)) {
        hout<< "Error in Substitute_boundary_point." <<endl;
        return 0;
    }
    //Substitute the outside point by the one at the boundary
    points_in[global_o].x = ipoi_vec[0].x;
    points_in[global_o].y = ipoi_vec[0].y;
    points_in[global_o].z = ipoi_vec[0].z;
        
    return 1;
}

int Cutoff_Wins::Get_intersecting_point_RVE_surface(struct Geom_RVE &sample, Point_3D &point0, Point_3D &point1, vector<Point_3D> &ipoi_vec)
{
	double t_temp[6];
	//”ÎX÷·¥π÷±∆Ω√Ê
	t_temp[0] = (sample.origin.x - point0.x)/(point1.x - point0.x);
	t_temp[1] = (sample.origin.x + sample.len_x - point0.x)/(point1.x - point0.x);
	//”ÎY÷·¥π÷±∆Ω√Ê
	t_temp[2] = (sample.origin.y - point0.y)/(point1.y - point0.y);
	t_temp[3] = (sample.origin.y + sample.wid_y - point0.y)/(point1.y - point0.y);
	//”ÎZ÷·¥π÷±∆Ω√Ê
	t_temp[4] = (sample.origin.z - point0.z)/(point1.z - point0.z);
	t_temp[5] = (sample.origin.z + sample.hei_z - point0.z)/(point1.z - point0.z);
    
	vector<double> t_ratio;
	for(int i=0; i<6; i++)
	{
		if(t_temp[i]>=0&&t_temp[i]<1)
		{
			//∂˛∑÷∑®≤Â»Î
			int left = 0;
			int right = (int)t_ratio.size()-1;
			while(right>=left)
			{
				int middle = (left + right)/2;
				if(fabs(t_ratio[middle] - t_temp[i])<Zero) goto T_Value_Same; //÷µœ‡Õ¨µƒ«Èøˆ
				else if(t_ratio[middle] > t_temp[i]) right = middle - 1;
				else left = middle + 1;
			}
			t_ratio.insert(t_ratio.begin()+left, t_temp[i]);	//≤Â»Î
        T_Value_Same: ;
		}
	}
    
	if((int)t_ratio.size()<1||(int)t_ratio.size()>3)
	{
		hout << "Get_intersecting_point_RVE_surface t_ratio.size()=" << (int)t_ratio.size() << ", …Ÿ”⁄“ª∏ˆªÚ’ﬂ∂‡”⁄»˝∏ˆ£° «ÎºÏ≤È£° " << endl;
		return 0;
	}
    
	Point_3D point_temp;
	for(int i=0; i<(int)t_ratio.size(); i++)
	{
		point_temp = point0+(point1-point0)*t_ratio[i];
		point_temp.flag = 1;		//÷–º‰µ„
		
		//---------------------------------------------------------------------------
		//ŒÛ≤Ó–ﬁ’˝
		if(fabs(point_temp.x-sample.origin.x)<Zero) point_temp.x = sample.origin.x;
		else if(fabs(point_temp.x-sample.origin.x-sample.len_x)<Zero) point_temp.x = sample.origin.x + sample.len_x;
        
		if(fabs(point_temp.y-sample.origin.y)<Zero) point_temp.y = sample.origin.y;
		else if(fabs(point_temp.y-sample.origin.y-sample.wid_y)<Zero) point_temp.y = sample.origin.y + sample.wid_y;
        
		if(fabs(point_temp.z-sample.origin.z)<Zero) point_temp.z = sample.origin.z;
		else if(fabs(point_temp.z-sample.origin.z-sample.hei_z)<Zero) point_temp.z = sample.origin.z + sample.hei_z;
        
		//---------------------------------------------------------------------------
		//≤Â»Î–¬µ„
		ipoi_vec.push_back(point_temp);
	}
    
	return 1;
}

void Cutoff_Wins::Trim_CNT(vector<vector<int> > &sectioned_domain_cnt, struct Geom_RVE sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<double> &radii, long int boundary, int CNT)
{
    //This bg variable is used to add the new CNT into the corresponding shell-sub-region
    Background_grid *bg = new Background_grid;
    
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
        bg->Add_to_shell(sample, points_in[P], sectioned_domain_cnt);
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
    double xmin = sample.origin.x;
    double ymin = sample.origin.y;
    double zmin = sample.origin.z;
    double lx = sample.len_x;
    double ly = sample.wid_y;
    double lz = sample.hei_z;
   int CNT = point3d.flag;
    if ( abs(x - xmin) < Zero){
        Add_CNT_to_boundary(boundary_cnt[0], CNT, point,0,0);
    } else if ( abs(x - (xmin+lx)) < Zero ){
        Add_CNT_to_boundary(boundary_cnt[1], CNT, point,0,1);
    } else if ( abs(y - ymin) < Zero ){
        Add_CNT_to_boundary(boundary_cnt[2], CNT, point,1,0);
    } else if ( abs(y - (ymin+ly)) < Zero ){
        Add_CNT_to_boundary(boundary_cnt[3], CNT, point,1,1);
    } else if ( abs(z - zmin) < Zero ) {
        Add_CNT_to_boundary(boundary_cnt[4], CNT, point,2,0);
    } else if ( abs(z - (zmin+lz)) < Zero ) {
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