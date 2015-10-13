//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Background_grid.cpp
//OBJECTIVE:	Create a background grid to make the CNT trimming faster
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Background_grid.h"

/*
 
 This classs creates a background grid that is used to determine which CNTs will be trimmed for each observation window. This is done
 in the class defined in Cutoff_Wins.cpp
 The grid will be created as follows: The smallest observation will be one sub-region. The rest of the subregions will have the size of the step size on each direction. So if the step size for increasing the observation window are Dx, Dy and Dz, then each subregion will be a paralelepiped with dimensions Dx, Dy and Dz.
 Then, for every point the code will work as follows:
 First check if the point is inside the central sub-region, i.e. the smallest observation window. If not then check if it belongs to another sub-region in the grid.
 
 Input:
    struct Geom_RVE sample
        Geometry of the generated sample
    struct Nanotube_Geo cnts
        Geometry of the CNTs
    vector<Point_3D> points_in
        List of points
 
 Output (These three are class variables):
    vector<vector<int> > sectioned_domain_cnt
    vector<vector<long int> > structure
 
 Modified inputs:
 
 */

int Background_grid::Generate_background_grids(const struct Geom_RVE &sample, const struct Nanotube_Geo &cnts, vector<Point_3D> &points_in)
{
    //Initialize the sectioned_domain_cnt vector
    //Empty vector to initialize the sectioned_domain_cnt vector
    vector<int> empty;
    //sample.cut_num is the number of increments from the smallest to the largest observation widows
    //Then, there are sample.cut_num+1 observation windows
    //Then, there are sample.cut_num+2 shells or sub-regions if we include the boundary layer.
    //Hence the sectioned_domain_cnt vector will have sample.cut_num+2 elements
    sectioned_domain_cnt.assign(sample.cut_num+2, empty);
    
    if (!Fill_structure_and_shell(sample, points_in)) {
        hout << "Error in Generate_background_grids." << endl;
        return 0;
    }
    
    return 1;
}

int Background_grid::Fill_structure_and_shell(const struct Geom_RVE &sample, vector<Point_3D> &points_in)
{
    //Vector to increase the size of structure
    vector<long int> empty;
    //Variable to count the number of CNTs
    int CNT = -1;
    //Scan all points to determine in which sub-regions the CNTs are located and construct the structure vector
    for (long int i = 0; i < points_in.size(); i++) {
        if (!points_in[i].flag) {
            //Increase the size of the structure
            structure.push_back(empty);
            CNT++;
        }
        //Add current point number
        structure.back().push_back(i);
        
        //Update CNT number
        points_in[i].flag = CNT;
        
        //Add to the corresponding shell
        if (!Add_to_shell(sample, points_in[i], sectioned_domain_cnt)) {
            hout << "Error in Fill_structure_and_shell"<< endl;
            return 0;
        }
        
        //Next point number
        CNT++;
    }
    
    return 1;
}

//This function finds the corresponding sub-region or shell where the point is located. Then it adds the CNT number to the
//corresponding vector in the 2D vector sectioned_domain_cnt
int Background_grid::Add_to_shell(struct Geom_RVE sample, Point_3D point, vector<vector<int> > &sectioned_domain_cnt)
{
    //Find the shell based on the x coordinate
    int shell_x = Find_shell(point.x, sample.origin.x, sample.len_x, sample.win_delt_x);
    //Find the shell based on the y coordinate
    int shell_y = Find_shell(point.y, sample.origin.y, sample.wid_y, sample.win_delt_y);
    //Find the shell based on the z coordinate
    int shell_z = Find_shell(point.z, sample.origin.z, sample.hei_z, sample.win_delt_z);
    
    //The shell to which the CNT of point is the largest shell from the three coordinates
    int shell;
    //With this if-statement I have shell = max(shell_x, shell_y)
    if (shell_x > shell_y)
        shell = shell_x;
    else
        shell = shell_y;
    //With this if-statement I have shell = max( shell_z, max(shell_x, shell_y))
    //And this is the largest shell value
    if (shell_z > shell)
        shell = shell_z;
            
    //Finally add the CNT on the corresponding sectioned domain
    //Add the CNT number to the sub-region, only if it is empty or if the last CNT is not the current CNT
    if ( (!sectioned_domain_cnt[shell].size()) || (sectioned_domain_cnt[shell].back() != point.flag) ){
        sectioned_domain_cnt[shell].push_back(point.flag);
    }
    
    return 1;
}

//This function finds the shell to which one coordinate belongs to
int Background_grid::Find_shell(double x_in, double x_min, double len_x, double dx)
{
    //Calculate the middle point of the observation window
    double x;
    double x_m = x_min + len_x/2;
    double x_max = x_min + len_x;
    
    //If the point is greater than the middle point, map it to the mirrored range
    if (x_in > x_m) {
        double m = (x_m - x_min)/(x_m - x_max);
        x = m*(x_in - x_m) + x_m;
    } else {
        x = x_in;
    }
    
    //Check if x is the outer shell
    if (x < x_min) {
        return 0;
    } else if(x > (x_max - dx)/2) {
        //Check if x is in the inner shell
        return (int)sectioned_domain_cnt.size()-1;
    } else {
        //Shell 0 is the boundary layer, so I need to increase in 1 the shell number
        int shell = (int) ((x + Zero)/dx);
        return shell+1;
    }
    
    hout << "Should not reach this part" << endl;
    return 0;
}