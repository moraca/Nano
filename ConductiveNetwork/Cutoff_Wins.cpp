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
    vector<int> cnts_inside
        List with the GNPs that are inside the observation window and that will be used for future computations
    vector<vector<int> > boundary_cnt
        Vectors that contains the CNTs on each of the six boundaries. These are needed to determine percolation
            boundary_cnt[0] corresponds to x0
            boundary_cnt[1] corresponds to x1
            boundary_cnt[2] corresponds to y0
            boundary_cnt[3] corresponds to y1
            boundary_cnt[4] corresponds to z0
            boundary_cnt[5] corresponds to z1
    vector<vector<int> > boundary_cnt
        Vectors that contains the CNTs on each of the six boundaries. These are needed to determine percolation
            boundary_gnp[0] corresponds to x0
            boundary_gnp[1] corresponds to x1
            boundary_gnp[2] corresponds to y0
            boundary_gnp[3] corresponds to y1
            boundary_gnp[4] corresponds to z0
            boundary_gnp[5] corresponds to z1
 
 Modified inputs:
    vector<vector<long int> > structure
    vector<double> radii
    vector<Point_3D> points_in
 
 */

//This function removes the points that are outside the observation window.
//The vector cnts_inside is created so only the CNTs inside the obseration window are considered in other functions
int Cutoff_Wins::Extract_observation_window(const int &window, const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const struct GNP_Geo &gnps, vector<GCH> &hybrid_particles, vector<vector<long int> > &structure, vector<vector<long int> > &structure_gnp, vector<double> &radii, vector<Point_3D> &points_in, vector<Point_3D> &points_gnp, vector<vector<int> > &shells_cnt, vector<vector<int> > &shells_gnp)
{
    if (!Set_global_variables_for_geometry(sample, window)) {
        hout << "Error in Extract_observation_window when calling Set_global_variables_for_geometry" << endl;
        return 0;
    }
    //Output the current window geometry
    hout<<"Observation window geometry:"<<endl;
    hout<<"xmin="<<xmin<<" ymin="<<ymin<<" zmin="<<zmin<<endl;
    hout<<"w_x="<<w_x<<" w_y="<<w_y<<" w_z="<<w_z<<endl;
    
    //Vector to save initial seeds
    vector<long int> seeds;
    //Save the initial points of the CNTs that are attached to the GNP
    //Check if particle type is the hybrid
    //hout << "Save seeds" << endl;
    if (sample.particle_type == "Hybrid_particles") {
        if (!Save_seeds(hybrid_particles, structure, seeds)) {
            hout << "Error in Extract_observation_window when calling Save_seeds" << endl;
            return 0;
        }
    }
    
    //Trim the CNTs if there are CNTs in the structure
    //Check if the generated structure has CNTs, this happens when the particle type is not GNPs
    if (sample.particle_type != "GNP_cuboids") {
        
        //Scan every Nanotube that is the boundary region. Delete and trim CNTs when needed.
        if (!Trim_boundary_cnts(shells_cnt, window, sample, points_in, structure, radii)){
            hout << "Error in Extract_observation_window when calling Trim_boundary_cnts" << endl;
            return 0;
        }
        
        //Fill the vector cnts_inside
        if (!Fill_cnts_inside(structure)) {
            hout << "Error in Extract_observation_window when calling Fill_cnts_inside" << endl;
            return 0;
        }
    }
    
    //Update the CNTs attached to the GNP
    //Check if particle type is hybrid
    if (sample.particle_type == "Hybrid_particles") {
        //hout << "Compare_seeds" << endl;
        //Compare the initial points of the CNTs attached to the GNPs
        //If they are different, that means that the CNT is not attached to the GNP anymore
        if (!Compare_seeds(hybrid_particles, structure, seeds)) {
            hout << "Error in Extract_observation_window when calling Compare_seeds" << endl;
            return 0;
        }
    }
    
    //Remove GNPs that are outside the observation window
    //Check if the structure has GNPs, this happens when the particle type is not CNT
    if (sample.particle_type != "CNT_wires") {
        //Scan every GNP that is the boundary region. Delete GNPs and GNP points when needed.
        if (!Trim_boundary_gnps(gnps, hybrid_particles, shells_gnp[window], points_gnp, structure_gnp)){
            hout << "Error in Extract_observation_window when calling Trim_boundary_gnps" << endl;
            return 0;
        }
    
        //Fill the vector gnps_inside
        if (!Fill_gnps_inside(structure_gnp)) {
            hout << "Error in Extract_observation_window when calling Fill_gnps_inside" << endl;
            return 0;
        }
    }
    
    return 1;
}
//This function sets global variables
int Cutoff_Wins::Set_global_variables_for_geometry(const struct Geom_RVE &sample, const int &window)
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
    
    return 1;
}

//This function scans all hybrid particles and saves the intial points of its CNTs
int Cutoff_Wins::Save_seeds(const vector<GCH> &hybrid_particles, const vector<vector<long int> > &structure, vector<long int> &seeds)
{
    //Initialize the seeds vector with the same size as structure
    seeds.clear();
    seeds.assign(structure.size(), -1);
    
    //Loop over the hybrid particles
    for (int i = 0; i < (int)hybrid_particles.size(); i++) {
        //variable to store the CNT number
        int CNT;
        
        //Scan top CNTs
        for (int j = 0; j < (int)hybrid_particles[i].cnts_top.size(); j++) {
            CNT = hybrid_particles[i].cnts_top[j];
            seeds[CNT] = structure[CNT].front();
        }
        
        //Scan bottom CNTs
        for (int j = 0; j < (int)hybrid_particles[i].cnts_bottom.size(); j++) {
            CNT = hybrid_particles[i].cnts_bottom[j];
            seeds[CNT] = structure[CNT].front();
        }
    }
    
    return 1;
}
//This function scans all hybrid particles and saves the intial points of its CNTs
//If no CNTs have their initial point on the GNP, then it is removed (as this means the GNP has no CNTS attached)
int Cutoff_Wins::Compare_seeds(vector<GCH> &hybrid_particles, const vector<vector<long int> > &structure, const vector<long int> &seeds)
{
    //Loop over the hybrid particles
    for (int i = (int)hybrid_particles.size()-1; i >= 0; i--) {
        //variable to store the CNT number
        int CNT;
        
        //Temporary variables
        vector<int> top_tmp, bottom_tmp;
        //Scan top CNTs
        for (int j = 0; j < (int)hybrid_particles[i].cnts_top.size(); j++) {
            CNT = hybrid_particles[i].cnts_top[j];
            //Check if seeds are still the same. If so, save the CNT number in the temporary variable
            if (structure[CNT].size() && seeds[CNT] == structure[CNT].front()) {
                top_tmp.push_back(CNT);
            }
        }
        //Update the vector of CNTs at the top surface of the GNP
        hybrid_particles[i].cnts_top = top_tmp;
        
        //Scan bottom CNTs
        for (int j = 0; j < (int)hybrid_particles[i].cnts_bottom.size(); j++) {
            CNT = hybrid_particles[i].cnts_bottom[j];
            //Check if seeds are still the same. If so, save the CNT number in the temporary variable
            if (structure[CNT].size() && seeds[CNT] == structure[CNT].front()) {
                bottom_tmp.push_back(CNT);
            }
        }
        //Update the vector of CNTs at the bottom surface of the GNP
        hybrid_particles[i].cnts_bottom = bottom_tmp;
        
    }
    
    return 1;
}

int Cutoff_Wins::Trim_boundary_cnts(vector<vector<int> > &shells_cnt, const int &window, const struct Geom_RVE &sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<double> &radii)
{
    //These variables will help me find the point location with respect to the box
    string currentPoint;
    //Variables for current and next points and current CNT
    long int P1;
    int CNT;
    //Empty vector to increase size of other vectors
    //Initialize the vector of boundary_flags with empty vectors
    vector<short int> empty_short;
    boundary_flags_cnt.assign(points_in.size(), empty_short);
    //Initialize the vector of boundary_flags with empty vectors
    vector<int> empty_int;
    boundary_cnt.assign(6, empty_int);
    //hout << "cnts_inside.size() = "<<cnts_inside.size()<<endl;
    //This varibale is used to initialize the vectors below
    vector<long int> empty_long;
    //Scan all CNTs in the current shell
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
                //   "inside" case. Hence, a new segment needs to be added by adding the local point number.
                //   In the following iterations, if the points are still inside,
                //   no point will be added since now the size of the vector is odd.
                if (branches_indices_CNT.size()%2 == 0) {
                    branches_indices_CNT.push_back(j);
                }
            } else {
                //If a point is not "inside", save the local point number only if the current size of branches_indices[CNT] is odd
                //When branches_indices[CNT].size() is odd, there was a change from an inside point to a non-inside point
                //Thus, add the local point number. In the following iterations, if the points are still not inside,
                //no point will be added since now the size of the vector is even
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
            int new_CNT = CNT;
            for (int k = 0; k < (int)branches_indices_CNT.size(); k=k+2) {
                int index1 = branches_indices_CNT[k];
                int index2 = branches_indices_CNT[k+1];
                
                //Beginning of segement
                if (!First_index(points_in, structure[CNT], new_CNT, index1)){
                    hout << "Error in Trim_boundary_cnts. branches_indices_CNT["<<k<<"]="<<branches_indices_CNT[k];
                    return 0;
                }
                //branches_indices[CNT][k] might be modified so I need to update it
                branches_indices_CNT[k] = index1;
                
                //End of segment
                if(!Second_index(points_in, structure[CNT], new_CNT, index2)){
                    hout << "Error in Trim_boundary_cnts. branches_indices_CNT["<<k+1<<"]="<<branches_indices_CNT[k+1];
                    return 0;
                }
                //branches_indices[CNT][k+1] might be modified so I need to update it
                branches_indices_CNT[k+1] = index2;
                
                //Just a check. I was getting CNTs with one point, so the first and second indices were equal
                //This should not happen anymore so I'll leave it just in case and for debugging
                if (index1 == index2) {
                    hout << "Error in Trim_boundary_cnts. index1 = index2 = "<<index1<<" on CNT "<<CNT<<endl;
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
                    for (int kk = index1; kk <= index2; kk++) {
                        long int P = structure[CNT][kk];
                        structure.back().push_back(P);
                        points_in[P].flag = new_CNT;
                        bg->Add_to_shell(sample, points_in[P], shells_cnt);
                    }
                    //Delete the temporary Background_vectors object
                    delete bg;
                    
                    //Update the radii vector
                    //The new CNT is just a segment of the old one, so they should have the same radius
                    radii.push_back(radii[CNT]);
                }
                //Update the new_CNT number, when there is one segment, the new value of this variable is not used
                //If there are two or more sements, then the new value of this variable is used
                new_CNT = (int)structure.size();
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

int Cutoff_Wins::First_index(vector<Point_3D> &points_in, vector<long int> &structure_CNT, int &new_CNT, int &index1)
{
    //Check if the first index is the first point of the CNT, otherwise add a boundary point
    if (index1 == 0) {
        //If the first index is just the initial point of a CNT just check if it is a boundary point
        //In such case, add it to the boundary vectors
        long int P = structure_CNT.front();
        if ( Where_is(points_in[P]) == "boundary") {
            Add_to_boundary_vectors(points_in[P], P, new_CNT);
        }
    } else {
        //If the first index is not the first point of the CNT, then I need to add a boundary point
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
            //Now, the outside point has the coordinates of the boundary point
        }
        //The first index is always inside, so:
        //    if the previous point is outside, the previous point needs to be included as it will now be the boundary point
        //    if the previous point is boundary, it has to be added
        //Hence, independently of where the previous point is, it has to be included, so it will be the new index
        index1--;
        //Since the first index is always inside, and since in this "else"-statement we already know that it is not the
        //first point of the CNT, then for sure a boundary point needs to be added. Whether because the previous index
        //is a boundary point or because we added a boundary point to the previous index
        Add_to_boundary_vectors(points_in[global_o], global_o, new_CNT);
    }
    return 1;
}

int Cutoff_Wins::Second_index(vector<Point_3D> &points_in, vector<long int> &structure_CNT, int &new_CNT, int &index2)
{
    //String to store the location of index2 (this location is used more than once so this avoids calling Where_is multiple times)
    string index2_location = Where_is(points_in[ structure_CNT[index2] ]);
    
    //Find out where the second index is and proceed accordingly
    if (index2_location == "outside"){
        //If the second index is onutside, then for sure I need to add a boundary point
        long int global_o = structure_CNT[index2];
        long int global_i = global_o-1;
        
        //Check if what is supposed to be the inside point is actually in the boundary.
        //If it hapens that the inside point is actually at the boundary, then this point has to be index2
        //and that's all, nothing more to do
        if ( Where_is(points_in[global_i]) == "boundary") {
            index2--;
        } else{
            //In this case, since the index2 point is outside and the previous point is inside,
            //we then calculate the projection to the boundary
            if (!Substitute_boundary_point(points_in, global_i, global_o)){
                hout << "Error in Second_index. global_i="<<global_i<<" global_o="<<global_o<<" structure_CNT.size()="<<structure_CNT.size();
                hout <<" index2="<<index2<<endl;
                hout <<"\tP_i=("<<points_in[global_i].x<<", "<<points_in[global_i].y<<", "<<points_in[global_i].z<<") P_o=(";
                hout <<points_in[global_o].x<<", "<<points_in[global_o].y<<", "<<points_in[global_o].z<<")"<<endl;
                return 0;
            }
            //Now, the outside point, i.e. index2, has the coordinates of the boundary point so I need to kep it unchanged
            Add_to_boundary_vectors(points_in[global_o], global_o, new_CNT);
        }
    } else if (index2_location == "boundary") {
        //If second index is at a boundary, just add it to the boundary vectors;
        //there is no need to find a projection with the boundary
        long int P = structure_CNT[index2];
        Add_to_boundary_vectors(points_in[P], P, new_CNT);
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
        return "inside";
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

int Cutoff_Wins::Get_intersecting_point_RVE_surface(const Point_3D &point0, const Point_3D &point1, vector<Point_3D> &ipoi_vec)
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
void Cutoff_Wins::Add_to_boundary_vectors(Point_3D point3d, long int point, int new_CNT)
{
    //Add point and CNT to the boundary vector
    double x = point3d.x;
    double y = point3d.y;
    double z = point3d.z;
    if ( abs(x - xmin) < Zero){
        Add_CNT_to_boundary(boundary_cnt[0], new_CNT, point,0,0);
    } else if ( abs(x - (xmin+w_x)) < Zero ){
        Add_CNT_to_boundary(boundary_cnt[1], new_CNT, point,0,1);
    } else if ( abs(y - ymin) < Zero ){
        Add_CNT_to_boundary(boundary_cnt[2], new_CNT, point,1,0);
    } else if ( abs(y - (ymin+w_y)) < Zero ){
        Add_CNT_to_boundary(boundary_cnt[3], new_CNT, point,1,1);
    } else if ( abs(z - zmin) < Zero ) {
        Add_CNT_to_boundary(boundary_cnt[4], new_CNT, point,2,0);
    } else if ( abs(z - (zmin+w_z)) < Zero ) {
        Add_CNT_to_boundary(boundary_cnt[5], new_CNT, point,2,1);
    }
}

//This function adds a CNT to the corresponding boundary vector.
//The flags are used in the direct electrifying algorithm:
//flag1: indicates the direction 0 is x, 1 is y, 2 is z
//flag2: indicates which boundary 0 is for x0, y0 or z0; 1 is for x1, y1 or z1
void Cutoff_Wins::Add_CNT_to_boundary(vector<int> &boundary, int CNT, long int point, short int flag1, short int flag2)
{
    if (!boundary.size()) {
        //If the boundary vector is empty, then just add the CNT
        boundary.push_back(CNT);
    } else if(boundary.back() != CNT){
        //If the boundary vector is not empty, add the CNT only if it has not been added
        boundary.push_back(CNT);
    }
    //If only one point of the CNT is outside, but this point is not one of the end points,
    //then we have two CNTs that will share a boundary point
    //This will cause the vector boundary_flags[point] to have 4 elements, which causes problems
    //when assigning node numbers in the LM matrix
    //So if the flags are only added when vector boundary_flags[point] is empty
    //The repetition of the point can be safely ignored since the two CNTs will be at the same boundary
    //Thus element boundary_flags[point][0] will be the same as boundary_flags[point][2]
    //and boundary_flags[point][1] will be the same as boundary_flags[point][3]
    if (!boundary_flags_cnt[point].size()) {
        boundary_flags_cnt[point].push_back(flag1);
        boundary_flags_cnt[point].push_back(flag2);

    }
}
//Fill the vector cnts_inside
int Cutoff_Wins::Fill_cnts_inside(const vector<vector<long int> > &structure)
{
    //Flag to determine invalid CNTs
    int flag = 0;
    //Scan all CNTs in the structure
    for (int i = 0; i < (int)structure.size(); i++) {
        //A CNT needs at least two points
        if (structure[i].size() > 1) {
            cnts_inside.push_back(i);
        } else if (structure[i].size() == 1) {
            hout<<"Error in Extract_observation_window. CNT "<<i<<" has only one point. A CNT must have at least 2 points."<<endl;
            flag = 1;
        }
    }
    
    //If the flag was set, then there were CNTs with one point
    //The function is not terminated at the first CNTs with one point found so that all these CNTs can be displayed in the output file
    if (flag)
        return 0;
    else
        return 1;
}
//
int Cutoff_Wins::Trim_boundary_gnps(const struct GNP_Geo &gnps, const vector<GCH> &hybrid_particles, const vector<int> &shell_gnp, vector<Point_3D> &points_gnp, vector<vector<long int> > &structure_gnp)
{
    //Empty vector to increase size of other vectors
    //Initialize the vector of boundary_flags with empty vectors
    vector<short int> empty_short;
    boundary_flags_gnp.assign(points_gnp.size(), empty_short);
    //Initialize the vector of boundary_flags with empty vectors
    vector<int> empty_int;
    boundary_gnp.assign(6, empty_int);
    
    //Scan all particles inside the current shell
    for (int i = 0; i < (int)shell_gnp.size(); i++) {
        //Current GNP number
        int gnp_n = shell_gnp[i];
        //Check if the center is outside the observation window or is close enough to the boundary of the window that part of it is outside
        if (Where_is(hybrid_particles[gnp_n].center) == "outside" || Is_close_to_boundaries(hybrid_particles[gnp_n])) {
            //If it outside or close enough check which points are still inside the window
            //and remove the ones outside
            if (!Remove_gnp_points_outside(gnps, points_gnp, structure_gnp[gnp_n])) {
                hout << "Error in Trim_boundary_gnps when calling Remove_gnp_points_outside" << endl;
                return 0;
            }
        }
    }
    return 1;
}
//This function determines if a GNP is close enough to the boundaries
int Cutoff_Wins::Is_close_to_boundaries(const GCH &hybrid)
{
    //Calculate half the diagonal of the GNP, i.e. the radius of the sphere that encloses the GNP
    double sum = hybrid.gnp.len_x*hybrid.gnp.len_x + hybrid.gnp.wid_y*hybrid.gnp.wid_y + hybrid.gnp.hei_z*hybrid.gnp.hei_z;
    double radius = sqrt(sum)/2;
    
    if (hybrid.center.x < xmin + radius || hybrid.center.x > xmin + w_x - radius) {
        return 1;
    } else if (hybrid.center.y < ymin + radius || hybrid.center.y > ymin + w_y - radius) {
        return 1;
    } else if (hybrid.center.z < zmin + radius || hybrid.center.z > zmin + w_z - radius) {
        return 1;
    }
    
    //If none of the cases above, then the GNP is not close to the boundaries
    return 0;
}
//Function that removes points from the discretization of GNPs
//Only points outside the observation window are removed
int Cutoff_Wins::Remove_gnp_points_outside(const struct GNP_Geo &gnps, vector<Point_3D> &points_gnp, vector<long int> &gnp_discrete)
{
    //vector to store projection points on each boundary
    vector<int> empty;
    vector<vector<int> > boundaries(6, empty);
    
    //Square of the discretization step length
    //This will reduce computations when comparing against squared distance instead of actual distance, in this was a calculation of a squared root is omited
    //A 10% is added for numerical error
    double min_discretization = (gnps.discr_step_length*gnps.discr_step_length)*1.1;
    
    //Vector that stores the points inside or at the boundary from the discretization
    vector<long int> gnp_discrete_in;
    
    //Get the first point in the discretization
    long int previousP = gnp_discrete.front();
    string previousPoint = Where_is(points_gnp[previousP]);
    
    //Scan all points from the given GNP
    //start on the previous to last point in the since elements will be deleted and the last point is the first previous point
    for (int i = 1; i < (int) gnp_discrete.size(); i++) {
        //Get point number
        long int P = gnp_discrete[i];
        
        //Get location of P
        string currentPoint = Where_is(points_gnp[P]);
        
        //Check if current and previous points make an inside-outside sequence
        if (!Find_inside_outside_sequence(i, P, previousP, min_discretization, currentPoint, previousPoint, points_gnp, gnp_discrete_in, boundaries)) {
             hout << "Error in Remove_gnp_points_outside when calling Find_inside_outside_sequence" << endl;
            return 0;
        }
        
        //Update previous variables
        previousPoint = currentPoint;
        previousP = P;
    }
    
    //If there are points at any boundary:
    // 1) Get the average point at each boundary
    // 2) Substitute a boundary point by the average point
    // 3) Fill boundary vectors
    if (!Add_GNPs_to_boundary(gnp_discrete, boundaries, points_gnp)) {
        hout << "Error in Remove_gnp_points_outside when calling Add_GNPs_to_boundary" << endl;
        return 0;
    }
    
    //From the discretization, only keep the points inside or at the boundary
    //Those points are in the vector gnp_discrete_in
    if (!Update_discretization(gnp_discrete, gnp_discrete_in)) {
        hout << "Error in Remove_gnp_points_outside when calling Update_discretization" << endl;
        return 0;
    }

    
    return 1;
}
//Function that determines if a sequence inside-outside is present between two consecutive points
//If present, it calls the function that finds the projection in the boundary
int Cutoff_Wins::Find_inside_outside_sequence(const int &iterator, const long int &currentP, const long int &previousP, const double &min_discretization, string &currentPoint, string &previousPoint, vector<Point_3D> &points_gnp, vector<long int> &gnp_discrete_in, vector<vector<int> > &boundaries) {
    
    //Check if the current point is outside the observation window
    if (currentPoint == "outside") {
        
        //Check if the previous point was inside
        //If the previous point was inside, then a boundary point might need to be added
        //Two sonecutive points in the discretization could be on different sides of the GNP, so the separation needs to be below or equal to the discretization step
        //To save time compare using the squared distance to avoid calculating an squared root
        //Because of numerical error a 10% is added to the limit (see declaration of variable min_discretization)
        if (previousPoint == "inside" && (points_gnp[currentP].squared_distance_to(points_gnp[previousP]) < min_discretization) ) {
            
            //Find the boundary point
            if (!Find_projection_in_boundary(points_gnp[previousP], points_gnp[currentP], iterator, boundaries)) {
                hout << "Error in Find_inside_outside_sequence when calling Find_projection_in_boundary" << endl;
                return 0;
            }
            
            //Now the outside point is at the boundary, so update the location
            currentPoint = "boundary";
            
            //The current point now stays in the discretization, so add it to the gnp_discrete_in vector
            gnp_discrete_in.push_back(currentP);
        }
    }
    //If the current point is inside, then check if the previous point was outside
    else if (currentPoint == "inside") {
        
        //If current point was inside but previous point was outside, then a boundary point is added
        //Also check if the two points are not too far from each other
        if (previousPoint == "outside" && (points_gnp[currentP].squared_distance_to(points_gnp[previousP]) < min_discretization) ) {
            
            //Find the boundary point
            if (!Find_projection_in_boundary(points_gnp[currentP], points_gnp[previousP], iterator, boundaries)) {
                hout << "Error in Find_inside_outside_sequence when calling Find_projection_in_boundary" << endl;
                return 0;
            }
            
            //The previous point now stays in the discretization, so add it to the gnp_discrete_in vector
            gnp_discrete_in.push_back(previousP);

        }
        
        //The current point stays in the discretization, so add it to the gnp_discrete_in vector
        gnp_discrete_in.push_back(currentP);
        
    }
    //If the current point is at the boundary check if it has already been added
    else if (currentPoint == "boundary") {
        
        //If current point was not added is because this point happens to be at the boundary (and was not calculated from a projection)
        if (gnp_discrete_in.back() != currentP) {
            gnp_discrete_in.push_back(currentP);
        }
    }
    //If location is something else, there is an error
    else {
        hout << "Error. Invalid location: " << currentPoint << endl;
        return 0;
    }
    
    return 1;
}
//Given two points in a sequence inside-outside, this function finds the projection on the boundary
//It also adds the point the boundary vector or vectors
int Cutoff_Wins::Find_projection_in_boundary(const Point_3D &inside, Point_3D &outside, const int &iterator, vector<vector<int> > &boundaries) {
    
    //Vector that stores the projection on the boundary
    vector<Point_3D> ipoi_vec;
    
    //Get projection on boundary
    if (!Get_intersecting_point_RVE_surface(outside, inside, ipoi_vec)) {
        hout << "Error in Find_projection_in_boundary when calling Get_intersecting_point_RVE_surface" << endl;
        return 0;
    }
    
    //Substitute the outside point by the boundary point
    outside.x = ipoi_vec[0].x;
    outside.y = ipoi_vec[0].y;
    outside.z = ipoi_vec[0].z;
    
    //Add point and CNT to the boundary vector
    if ( abs(outside.x - xmin) < Zero ) {
        boundaries[0].push_back(iterator);
    } else if ( abs(outside.x - (xmin+w_x)) < Zero ) {
        boundaries[1].push_back(iterator);
    } else if ( abs(outside.y - ymin) < Zero ){
        boundaries[2].push_back(iterator);
    } else if ( abs(outside.y - (ymin+w_y)) < Zero ) {
        boundaries[3].push_back(iterator);
    } else if ( abs(outside.z - zmin) < Zero ){
        boundaries[4].push_back(iterator);
    } else if ( abs(outside.z - (zmin+w_z)) < Zero ) {
        boundaries[5].push_back(iterator);
    }
    
    return 1;
}
//This function updates the discretization by keeping only those points inside the sample or at the boundaries
int Cutoff_Wins::Update_discretization(vector<long int> &gnp_discrete, vector<long int> &gnp_discrete_in) {
    
    //resize the gnp_discrete vector
    gnp_discrete.resize(gnp_discrete_in.size());
    
    //Copy the elements of gnp_discrete_in into gnp_discrete
    for (int i = 0; i < (int)gnp_discrete_in.size(); i++) {
        gnp_discrete[i] = gnp_discrete_in[i];
    }
    
    return 1;
}
//This function add a GNP and GNP point to the boundary vector
int Cutoff_Wins::Add_GNPs_to_boundary(const vector<long int> &gnp_discrete, const vector<vector<int> > &boundaries, vector<Point_3D> &points_gnp) {
    
    //Scan boundaries vector
    for (int i = 0; i < (int)boundaries.size(); i++) {
        
        //Check if there are any points at the boundary
        if (boundaries[i].size()) {
            //Find average point of current boundary
            if (!Find_average_boundary_point(gnp_discrete, boundaries[i], points_gnp)) {
                hout << "Error in Add_GNPs_to_boundary when calling Find_average_boundary_point at iteration " << i << endl;
                return 0;
            }
            
            //Get the first point at the boundary, which is the average point
            long int P0 = gnp_discrete[ boundaries[i][0] ];
            
            //Get GNP number, use the avergage point (all points have the same GNP number)
            int GNP = points_gnp[P0].flag;
            
            //Add GNP number to the corresponding boundary vector
            boundary_gnp[i].push_back(GNP);
            
            //Calculate the first flag for the boundary flags
            short int flag1 = (short int)(i/2);
            
            //Calculate the second flag for the boundary flags
            short int flag2 = (short int) i%2;
            
            //Add flags to boundary_flags_gnp
            boundary_flags_gnp[P0].push_back(flag1);
            boundary_flags_gnp[P0].push_back(flag2);
            //hout<<"boundary_gnp["<<i<<"].push_back("<<GNP<<"), "<< "P0="<<P0<<", flag1="<<flag1<<", flag2="<<flag2<<endl;
        }
    }
    
    return 1;
}
//Given a set of GNP points at the boundary, find the average point
int Cutoff_Wins::Find_average_boundary_point(const vector<long int> &gnp_discrete, const vector<int> &boundary, vector<Point_3D> &points_gnp) {
    
    //Number of points in boundary
    int n_points = (int)boundary.size();
    
    //The average point will be stored in the first point indicated by the boundary indices
    long int P0 = gnp_discrete[boundary[0]];
    
    //Save the GNP number
    //This is done beacause the sum of two points returns a point with uninitialized flag, thus the GNP number is lost
    int GNP = points_gnp[P0].flag;
    
    //Iterate over points in boundary, start in the second point as the first will store the average
    for (int i = 1; i < n_points; i++) {
        //The boundary vector has indices of the points in the vector of discrete GNPs
        int index = boundary[i];
        
        //Get the point number
        long int P = gnp_discrete[index];
        
        //Add boundary point to average
        //Use the first point to store the average point
        points_gnp[P0] = points_gnp[P0] + points_gnp[P];
    }
    
    //Calculate average
    points_gnp[P0] = points_gnp[P0]/n_points;
    
    //Update the GNP number
    points_gnp[P0].flag = GNP;
    
    return 1;
}
//Function that fills the vector gnps_inside
int Cutoff_Wins::Fill_gnps_inside(const vector<vector<long int> > &structure_gnp)
{
    //Scan all GNP discretizations
    for (int i = 0; i < (int)structure_gnp.size(); i++) {
        //Check if there are any points on the GNP
        if (structure_gnp[i].size()) {
            //If there are still points in the discretization, that means that at least a fraction of the GNP is inside the observation window
            gnps_inside.push_back(i);
        }
    }
    return 1;
}