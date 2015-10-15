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

int Background_vectors::Generate_shells_and_structure(const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, const vector<Point_3D> &points_in, vector<vector<int> > &shells_cnt)
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
    
    if (!Fill_structure_and_shell(sample, points_in, shells_cnt)) {
        hout << "Error in Generate_Background_vectorss." << endl;
        return 0;
    }
    //hout <<"shells_cnt[0].size()="<<shells_cnt[0].size()<<" shells_cnt[1].size()="<<shells_cnt[2].size()<<endl;
    
    return 1;
}

//This function generates the vectors shells_cnt and structure
int Background_vectors::Fill_structure_and_shell(const struct Geom_RVE &sample, const vector<Point_3D> &points_in, vector<vector<int> > &shells_cnt)
{
    //Vector to increase the size of structure
    vector<long int> empty;
    //Scan all points to determine in which sub-regions the CNTs are located and construct the structure vector
    for (long int i = 0; i < (long int)points_in.size(); i++) {
        //Add to the corresponding shell
        if (!Add_to_shell(sample, points_in[i], shells_cnt)) {
            hout << "Error in Fill_structure_and_shell"<< endl;
            return 0;
        }
    }
    
    return 1;
}

//This function finds the corresponding shell sub-region where the point is located. Then it adds the CNT number to the
//corresponding vector in the 2D vector shells_cnt
int Background_vectors::Add_to_shell(const struct Geom_RVE &sample, const Point_3D &point, vector<vector<int> > &shells_cnt)
{
    //Find the shell based on the x coordinate
    int shell_x = Find_shell(point.x, sample.origin.x, sample.len_x, sample.win_delt_x, sample.win_min_x, sample.win_max_x, shells_cnt);
    //Find the shell based on the y coordinate
    int shell_y = Find_shell(point.y, sample.origin.y, sample.wid_y, sample.win_delt_y, sample.win_min_y, sample.win_max_y, shells_cnt);
    //Find the shell based on the z coordinate
    int shell_z = Find_shell(point.z, sample.origin.z, sample.hei_z, sample.win_delt_z, sample.win_min_z, sample.win_max_z, shells_cnt);
    
    //The shell to which the CNT of point is the outer-most shell from the three coordinates, that is, the minimum shell number overall
    int shell = shell_x;
    //With this if-statement I have shell = min(shell_x, shell_y)
    if (shell_y < shell)
        shell = shell_y;
    //With this if-statement I have shell = min( shell_z, min(shell_x, shell_y))
    //And this is the smallest shell value
    if (shell_z < shell)
        shell = shell_z;
            
    //Finally add the CNT on the corresponding sectioned domain
    //hout << "shell="<<shell<<endl;
    //Add the CNT number to the shell sub-region, only if it is empty or if the last CNT is not the current CNT
    if ( (!shells_cnt[shell].size()) || (shells_cnt[shell].back() != point.flag) ){
        	shells_cnt[shell].push_back(point.flag);
    }
    
    return 1;
}

//This function finds the shell to which one coordinate belongs to
int Background_vectors::Find_shell(double x_in, double x_min, double len_x, double dx, double win_min_x, double win_max_x, vector<vector<int> > &shells_cnt)
{
    //Calculate the middle point of the observation window
    double x_m = x_min + len_x/2;
    //The maximum coordinate
    double x_max = x_min + len_x;
    //An observation window grows a total of dx. However, as it is centered, it gros dx/2 on each direction
    dx = dx/2;
    
    //Effective coordinate
    double x;
    //If the point is greater than the middle point, map it to the mirrored range
    if (x_in > x_m) {
        double m = (x_m - x_min)/(x_m - x_max);
        x = m*(x_in - x_m) + x_m;
    } else {
        x = x_in;
    }
    //Check if x is in the outer shell (it is in the outer shell when the point is below the x)
    //hout << "x_in="<<x_in<<" x="<<x<<" shell=";
    if (x < x_m - win_max_x/2) {
        //hout<<'0'<<endl;
        return 0;
    } else if( x > x_m - (win_min_x/2) ) {
        //Check if x is in the inner shell
        //hout<<	shells_cnt.size()-1<<endl;
        return (int)	shells_cnt.size()-1;
    } else {
        //Shell 0 is the boundary layer, so I need to increase in 1 the shell number
        int shell = (int) ((x - x_min + Zero)/dx);
        //hout<<shell+1<<endl;
        return shell+1;
    }
    
    hout << "Should not reach this part" << endl;
    return 0;
}