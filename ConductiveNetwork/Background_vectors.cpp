//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Background_vectors.cpp
//OBJECTIVE:	Using nested shells on the background to mark the CNTs for trimming faster in each observation window
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Background_vectors.h"

/*
 
 This class generates the vector:
 shells_cnt: is used to trim the CNTs. The CNTs are grouped according to the sub-regions, then need to be trimmed
 The shells will be created as follows: The smallest observation window (cuboid) will be one shell sub-region and will be used for the last element in the vector.
 Then, the next shell-subregion will be the volume of the next observation window minus the volume of the first one (the smallest observation window).
 The following shells will have the same form: the volume of the observation window minus the volume of the previous one.
 The last shell region will be the boundary layer. This will be used for the first element of the vector
 
 Input:
    struct Geom_RVE sample
        Geometry of the generated sample
    struct Nanotube_Geo cnts
        Geometry of the CNTs
    vector<Point_3D> points_in
        List of points
 
 Output:
    vector<vector<int> > shells_cnt
 
 Modified inputs:
 
 */

int Background_vectors::Generate_shells(const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const vector<Point_3D> &points_in, const vector<GCH> &hybrid_particles, vector<vector<int> > &shells_cnt, vector<vector<int> > &shells_gnps)
{
    //Initialize the shells_cnt vector
    //Empty vector to initialize the shells_cnt vector
    vector<int> empty;
    //sample.cut_num is the number of increments from the smallest to the largest observation widows
    //Then, there are sample.cut_num+1 observation windows
    //Then, there are sample.cut_num+2 shell sub-regions if we include the boundary layer.
    //Hence the shells_cnt vector will have sample.cut_num+2 elements
    shells_cnt.assign(sample.cut_num+2, empty);
    //hout << "sample.cut_num+2="<<sample.cut_num+2<<endl;
    
    //Scan all points to determine in which sub-regions the CNTs are located and construct the structure vector
    for (long int i = 0; i < (long int)points_in.size(); i++) {
        //Add to the corresponding shell
        if (!Add_to_shell(sample, points_in[i], shells_cnt)) {
            hout << "Error in Generate_shells_and_structure."<< endl;
            return 0;
        }
    }
    
    //Initialize the shells_gnps vector with the same size as the shells_cnt
    shells_gnps.assign(sample.cut_num+2, empty);
    
    //Scan all Hybrid particles to determine to which shells they belong
    for (int i = 0; i < (int)hybrid_particles.size(); i++) {
        //Add to the corresponding shell
        if (!Add_to_shells(sample, hybrid_particles[i], shells_gnps)) {
            hout << "Error in Generate_shells_and_structure."<< endl;
            return 0;
        }    }
    
    return 1;
}

//This function finds the corresponding shell sub-region where the point is located. Then it adds the CNT number to the
//corresponding vector in the 2D vector shells_cnt
int Background_vectors::Add_to_shell(const struct Geom_RVE &sample, const Point_3D &point, vector<vector<int> > &shells_cnt)
{
    //Number of shells
    int num_shells = (int)shells_cnt.size();
    
    //Find the shell that corresponds to the point
    int shell = Find_minimum_shell(sample, point, num_shells);
            
    //Finally add the CNT on the corresponding sectioned domain
    //hout << "shell="<<shell<<endl;
    //Add the CNT number to the shell sub-region, only if it is empty or if the last CNT is not the current CNT
    if ( (!shells_cnt[shell].size()) || (shells_cnt[shell].back() != point.flag) ){
        	shells_cnt[shell].push_back(point.flag);
    }
    
    return 1;
}

//This function finds the shell to which one point belongs to
//So it uses three times the function that finds the shell to which one coordinate belongs to
int Background_vectors::Find_minimum_shell(const struct Geom_RVE &sample, const Point_3D &point, const int &num_shells)
{
    //Find the shell based on the x coordinate
    //hout << "x_in=";
    int shell_x = Find_shell(point.x, sample.origin.x, sample.len_x, sample.win_delt_x, sample.win_min_x, sample.win_max_x, num_shells);
    //Find the shell based on the y coordinate
    //hout << "y_in=";
    int shell_y = Find_shell(point.y, sample.origin.y, sample.wid_y, sample.win_delt_y, sample.win_min_y, sample.win_max_y, num_shells);
    //Find the shell based on the z coordinate
    //hout << "z_in=";
    int shell_z = Find_shell(point.z, sample.origin.z, sample.hei_z, sample.win_delt_z, sample.win_min_z, sample.win_max_z, num_shells);
    
    //The shell to which the CNT of point is the outer-most shell from the three coordinates, that is, the minimum shell number overall
    int shell = shell_x;
    //With this if-statement I have shell = min(shell_x, shell_y)
    if (shell_y < shell)
        shell = shell_y;
    //With this if-statement I have shell = min( shell_z, min(shell_x, shell_y))
    //And this is the smallest shell value
    if (shell_z < shell)
        shell = shell_z;
    
    return shell;
}

//This function finds the shell to which one coordinate belongs to
int Background_vectors::Find_shell(const double &x_in, const double &x_0, const double &len_x, const double &dx, const double &win_min_x, const double &win_max_x, const int &num_shells)
{
    //Calculate the middle point of the sample
    double x_m = x_0 + len_x/2;
    //Calculate the minimum coordinate of the largest observation window
    //below this corrdinate, x is outside the largest observation window,
    //or in the boundary layer if the largest observation window equals the sample size
    double x_layer = x_0 + (len_x - win_max_x)/2;
    //The minimum coordinate of the smallesr observation window
    //any point above this coordinate is inside the core shell
    double x_core = x_0 + (len_x - win_min_x)/2;
    //An observation window grows a total of dx. However, as it is centered, it grows dx/2 on each direction
    double dx_half = dx/2;
    
    //Effective coordinate
    double x;
    //If the point is greater than the middle point, map it to the mirrored range
    if (x_in > x_m) x = 2*x_m - x_in; 
    else x = x_in;
    
    //Check if x is in the outer shell (it is in the outer shell when the point is below the x_min)
    //hout<<x_in<<" x="<<x<<" shell";
    if (x < x_layer) {
        //hout<<"m=0"<<endl;
        return 0;
    } else if( x > x_core ) {
        //Check if x is in the inner shell
        //hout<<"M="<<shells_cnt.size()-1<<endl;
        return num_shells-1;
    } else {
        //I need the integer part of (x - x_min)/dx_half + 1
        //The Zero is for foating point errors
        int shell = (int) ((x - x_layer + Zero)/dx_half) + 1;
        //hout<<"=f("<<((x - x_layer + Zero)/dx_half)+1<<")="<<shell<<endl;
        return shell;
    }
}

//This function finds the corresponding shell sub-region where the point is located. Then it adds the CNT number to the
//corresponding vector in the 2D vector shells_cnt
int Background_vectors::Add_to_shells(const struct Geom_RVE &sample, const GCH &hybrid, vector<vector<int> > &shells_gnp)
{
    //Number of shells
    int num_shells = (int)shells_gnp.size();
    
    //Range of shells in which the GNP will be stored
    int shell_max = -1; //Initialized in a negative value to find a larger value
    int shell_min = num_shells+1; //Initialized in a large value to find a smaller value
    
    //Base point to find the eight corner points
    Point_3D corner_base( -hybrid.gnp.len_x/2, -hybrid.gnp.wid_y/2, -hybrid.gnp.hei_z/2);
    
    //Loop over the eigt possible corners of the cube
    for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
            for(int k=0; k<2; k++) {
                
                //If i (j,k) is zero, then add nothing to the x (y,z) coordinate
                //If i (j,k) is one, then add the length on direction x (y,z) to the x (y,z) coordinate
                Point_3D corner(((double)i)*hybrid.gnp.len_x, ((double)j)*hybrid.gnp.wid_y, ((double)k)*hybrid.gnp.hei_z);
                
                //Add the center and the "adjustment" so that the loop calculates the eight coordinates are calculated
                corner = corner_base + corner;
                
                //Map to global coordinates
                corner = corner.rotation(hybrid.rotation, hybrid.center);
                
                //Find the shell that corresponds to the point
                int shell = Find_minimum_shell(sample, corner, num_shells);
                
                //Update maximum and minimum shells
                if (shell < shell_min)
                    shell_min = shell;
                if (shell > shell_max)
                    shell_max = shell;
            }
    
    //Add to all shells in range
    for (int i = shell_min; i <= shell_max; i++) {
        //Add the GNP number to a shell sub-region
        shells_gnp[i].push_back(hybrid.flag);
    }    
    
    return 1;
}