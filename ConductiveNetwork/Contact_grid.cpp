//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Contact_grid.cpp
//OBJECTIVE:	Create a background grid to find the points in contact faster
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Contact_grid.h"

/*
 
 This function groups the points inside the observation window into sub-regions.
 
 The task of finding overlapping points can be computationally expensive. To reduce computational cost, the sample was divided into smaller cubic sub-regions. The cubic sample of size $[a \times a \times a ]$ was divided into $m$ segments in each direction resulting in $m^3$ cubic sub-regions. Then, overlapping points are searched for only on each smaller sub-region rather than in the whole sample. It may happen that two overlapping points belong to different sub-regions. In order to take into account these overlapping points, each sub-region is ``extended". If a point lies exactly at a boundary, then an overlapping point can be at a maximum distance equal to $r_{max}+r_{max}+d_W = 2r_{max}+d_W$ from the boundary. Here, $r_{max}$ is the maximum radii of the CNTs. Then, each boundary plane of each cubic sub-region is translated externally to a parallel plane at a distance $2r_{max}+d_W$.
 
 Input:
    struct Geom_RVE sample
        Goemetry of the simulated sample
    struct Cutoff_dist cutoffs
        Structure that contains the cutoff for tunneling and overlapping
    struct Nanotube_Geo cnts
        Structure that contains the geometry of the CNTs
    vector<int> cnts_inside
        List of CNTs that are inside the observation window
    vector<vector<long int> > structure
        Vector with the structure
    vector<Point_3D> points_in
        List of points
    int window
        Current observation window. window=0 means the largest observation window
 
 Output (These three are class variables):
    vector<vector< long int> > sectioned_domain
        Vector with overlaping sub-regions. This vector is used to look for contact points in CNTs faster
    vector<vector< long int> > sectioned_domain_gnps
        Vector with overlaping sub-regions. This vector is used to look for contact points in GNPs faster
 
 */

int Contact_grid::Generate_contact_grid(const int &window, const struct Geom_RVE &sample, const struct Cutoff_dist &cutoffs, const struct Nanotube_Geo &cnts, const vector<int> &cnts_inside, vector<Point_3D> &points_in, const vector<vector<long int> > &structure, const vector<int> &gnps_inside, const vector<Point_3D> &points_gnp, const vector<vector<long int> > &structure_gnp)
{
    //Generate the window geometry
    struct Geom_RVE window_geom;
    if (!Generate_window_geometry(window, sample, window_geom)) {
        hout << "Error in Generate_contact_grid when calling Generate_window_geometry" << endl;
        return 0;
    }
    
    //Number of regions on each direction
    int sx = (int)(window_geom.len_x/window_geom.gs_minx);
    int sy = (int)(window_geom.wid_y/window_geom.gs_miny);
    int sz = (int)(window_geom.hei_z/window_geom.gs_minz);
    
    //Maximum distance between two points in contact inside the sample
    double cutoff = cutoffs.tunneling_dist + 2*cnts.rad_max;
    
    //Check if the number of regions needs to be adjusted
    if (!Adjust_regions_if_needed(cutoff, window_geom, sx, sy, sz)) {
        hout << "Error in Generate_contact_grid when calling Adjust_regions_if_needed" << endl;
        return 0;
    }
    //hout<<"There are "<<sx*sy*sz<<" overlapping sub-regions."<<endl;
    
    //Fill the vector for the sectioned domain for CNTs only when there are CNTs in the structure
    //Equivalently, this is done when the structure is not made of only GNPs
    if (sample.particle_type != "GNP_cuboids") {
        
        if (!Fill_sectioned_domain(window_geom, cnts_inside, structure, points_in, cutoff, sx, sy, sz, sectioned_domain)) {
            hout << "Error in Generate_contact_grid when calling Generate_sectioned_domain_cnts" << endl;
            return 0;
        }
    }
    
    //Fill the vector for the sectioned domain for discrete GNPs only when there are GNPs in the structure
    //Equivalently, this is done when the structure is not made of only CNTs
    if (sample.particle_type != "CNT_wires") {
        
        //Fill the sectioned domain corresponding to GNP points
        if (!Fill_sectioned_domain(window_geom, gnps_inside, structure_gnp, points_gnp, cutoff, sx, sy, sz, sectioned_domain_gnps)) {
            hout << "Error in Generate_contact_grid when calling Generate_sectioned_domain_cnts" << endl;
            return 0;
        }
        
        //Fill the sectioned domain corresponding to GNP numbers
        if (!Fill_sectioned_domain(window_geom, gnps_inside, structure_gnp, points_gnp, cutoff, sx, sy, sz, sectioned_domain_hyb)) {
            hout << "Error in Generate_contact_grid when calling Generate_sectioned_domain_cnts" << endl;
            return 0;
        }
    }
    
    return 1;
}
//
int Contact_grid::Generate_window_geometry(const int &window, const struct Geom_RVE &sample, struct Geom_RVE &window_geom)
{
    //Dimensions of the current observation window
    window_geom.len_x = sample.win_max_x - ((double)window)*sample.win_delt_x;
    window_geom.wid_y = sample.win_max_y - ((double)window)*sample.win_delt_y;
    window_geom.hei_z = sample.win_max_z - ((double)window)*sample.win_delt_z;
    
    //Check that the current observation window is not smaller than the minimum observation window
    if (window_geom.len_x < sample.win_min_x) {
        //If the current observation window is smaller than the minimum observation window, set the observation window equal to the minimum
        window_geom.len_x = sample.win_min_x;
    }
    if (window_geom.wid_y < sample.win_min_y) {
        //If the current observation window is smaller than the minimum observation window, set the observation window equal to the minimum
        window_geom.wid_y = sample.win_min_y;
    }
    if (window_geom.hei_z < sample.win_min_z) {
        //If the current observation window is smaller than the minimum observation window, set the observation window equal to the minimum
        window_geom.hei_z = sample.win_min_z;
    }
    
    //These variables are the coordinates of the lower corner of the observation window
    window_geom.origin.x = sample.origin.x + (sample.len_x - window_geom.len_x)/2;
    window_geom.origin.y = sample.origin.y + (sample.wid_y - window_geom.wid_y)/2;
    window_geom.origin.z = sample.origin.z + (sample.hei_z - window_geom.hei_z)/2;
    
    //Sizes of each region
    window_geom.gs_minx = sample.gs_minx;
    window_geom.gs_miny = sample.gs_miny;
    window_geom.gs_minz = sample.gs_minz;
    
    //hout<<"Observation window geometry:"<<endl;
    //hout<<"window_geom.origin.x="<<window_geom.origin.x<<" window_geom.origin.y="<<window_geom.origin.y<<" window_geom.origin.z="<<window_geom.origin.z<<endl;
    //hout<<"window_geom.len_x="<<window_geom.len_x<<" window_geom.wid_y="<<window_geom.wid_y<<" window_geom.hei_z="<<window_geom.hei_z<<endl;
    
    return 1;
}
//
int Contact_grid::Adjust_regions_if_needed(const double &cutoff, struct Geom_RVE &window_geom, int &sx, int &sy, int &sz)
{
    //hout << "sx = " << sx << '\t' << "dx = " << window_geom.gs_minx << "\n";
    //hout << "sy = " << sy << '\t' << "dy = " << window_geom.gs_miny << "\n";
    //hout << "sz = " << sz << '\t' << "dz = " << window_geom.gs_minz  << "\n";
    
    //Check that the regions are not too small for the maximum cutoff distance 2r_max+tunnel
    //If they are, then change the number of sections to the maximum possible
    if (window_geom.gs_minx < 2*cutoff) {
        sx = (int)(window_geom.len_x/(2*(cutoff+Zero)));
        window_geom.gs_minx = window_geom.len_x/(double)sx;
        hout << "Modified the number of sections along x. " << "sx = " << sx << '\t' << "dx = " << window_geom.gs_minx << endl;
    }
    if (window_geom.gs_miny < 2*cutoff) {
        sy = (int)(window_geom.wid_y/(2*(cutoff+Zero)));
        window_geom.gs_miny = window_geom.wid_y/(double)sy;
        hout << "Modified the number of sections along y. " << "sy = " << sy << '\t' << "dy = " << window_geom.gs_miny << endl;
    }
    if (window_geom.gs_minz < 2*cutoff) {
        sz = (int)(window_geom.hei_z/(2*(cutoff+Zero)));
        window_geom.gs_minz = window_geom.hei_z/(double)sz;
        hout << "Modified the number of sections along z. " << "sz = " << sz << '\t' << "dz = " << window_geom.gs_minz  << endl;
    }
    return 1;
}
//
int Contact_grid::Fill_sectioned_domain(const struct Geom_RVE &window_geom, const vector<int> &particles_inside, const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const double &cutoff, const int &sx, const int &sy, const int &sz, vector<vector< long int> > &sectioned_domain)
{
    //There will be sx*sy*sz different regions
    sectioned_domain.clear();
    vector<long int> empty_long;
    sectioned_domain.assign(sx*sy*sz, empty_long);
    
    //First loop over the particles inside the box, then loop over the points inside each particle
    for (int i = 0; i < (int)particles_inside.size(); i++) {
        //hout<<"particles_inside["<<i<<"]="<<particles_inside[i]<<' ';
        
        //Current particle number (either CNT or GNP)
        int particle = particles_inside[i];
        //hout<<"structure["<<particle<<"].size()="<<structure[particle].size()<<endl;
        
        //Scan each point in the particle
        for (int j = 0; j < (int)structure[particle].size(); j++) {
            
            //Current point number
            long int P = structure[particle][j];
            //hout<<"P=structure["<<particle<<"]["<<j<<"]="<<structure[particle][j]<< endl;
            
            //Calculate the region-coordinates
            int a, b, c;
            //hout << "points_in[P]=("<<points_in[P].x<<", "<<points_in[P].y<<", "<<points_in[P].z<<")"<< endl;
            if (!Calculate_region_coordinates(window_geom, points_in[P], sx, sy, sz, a, b, c)) {
                hout << "Error in Generate_sectioned_domain_cnts when calling Calculate_region_coordinates" << endl;
                return 0;
            }
            
            //Initialize flags for overlaping regions
            int fx = 0;
            int fy = 0;
            int fz = 0;
            
            //Assign value of flag according to position of point
            if (!Calculate_postion_flags(window_geom, points_in[P], cutoff, a, b, c, sx, sy, sz, fx, fy, fz)) {
                hout << "Error in Generate_sectioned_domain_cnts when calling Calculate_postion_flags" << endl;
                return 0;
            }
            
            //Assign the point P to the correspoding region or regions
            if (!Assign_point_to_region(a, b, c, fx, fy, fz, sx, sy, P, sectioned_domain)) {
                hout << "Error in Generate_sectioned_domain_cnts when calling Assign_point_to_region" << endl;
                return 0;
            }
            
        }
    }
    
    return 1;
}
//
int Contact_grid::Calculate_region_coordinates(const struct Geom_RVE &window_geom, const Point_3D &point, const int &sx, const int &sy, const int &sz, int &a, int &b, int &c)
{
    
    //Calculate the region-coordinates
    a = (int)((point.x-window_geom.origin.x)/window_geom.gs_minx);
    //Limit the value of a as it has to go from 0 to sx-1
    if (a == sx) {
        a--;
        //hout << " a-- ";
    }
    b = (int)((point.y-window_geom.origin.y)/window_geom.gs_miny);
    //Limit the value of b as it has to go from 0 to sy-1
    if (b == sy) {
        b--;
        //hout << " b-- ";
    }
    c = (int)((point.z-window_geom.origin.z)/window_geom.gs_minz);
    //Limit the value of c as it has to go from 0 to sz-1
    if (c == sz){
        c--;
        //hout << " c-- ";
    }
    
    return 1;
}
//
int Contact_grid::Calculate_postion_flags(const struct Geom_RVE &window_geom, const Point_3D &point, const double &cutoff, const int &a, const int &b, const int &c, const int &sx, const int &sy, const int &sz, int &fx, int &fy, int &fz)
{
    //Coordinates of non-overlaping region the point belongs to
    double x1 = a*window_geom.gs_minx +  window_geom.origin.x;
    double x2 = x1 + window_geom.gs_minx;
    double y1 = b*window_geom.gs_miny +  window_geom.origin.y;
    double y2 = y1 + window_geom.gs_miny;
    double z1 = c*window_geom.gs_minz +  window_geom.origin.z;
    double z2 = z1 + window_geom.gs_minz;
    
    //Assign value of flag according to position of point
    //The first operand eliminates the periodicity on the boundary
    if ((a > 0) && (point.x >= x1) && (point.x <= x1+cutoff))
        fx = -1;
    else if ((a < sx-1) && (point.x >= x2-cutoff) && (point.x <= x2 ))
        fx = 1;
    if ((b > 0) && (point.y >= y1) && (point.y <= y1+cutoff))
        fy = -1;
    else if ((b < sy-1) && (point.y >= y2-cutoff) && (point.y <= y2 ))
        fy = 1;
    if ((c > 0) && (point.z >= z1) && (point.z <= z1+cutoff))
        fz = -1;
    else if ((c < sz-1) && (point.z >= z2-cutoff) && (point.z <= z2 ))
        fz = 1;
    
    return 1;
}
//
int Contact_grid::Assign_point_to_region(const int &a, const int &b, const int &c, const int &fx, const int &fy, const int &fz, const int &sx, const int &sy, const long int &P, vector<vector< long int> > &sectioned_domain)
{
    //Create array for loop over overlaping regions
    int temp[2][3] = { {a+fx, b+fy, c+fz}, {a, b, c}};
    
    //In this loop I check all regions a point can belong to when it is in an overlaping zone
    for (int ii = 0; ii < 2; ii++) {
        if (!fx) ii++; //if flag is zero, do this loop only once
        for (int jj = 0; jj < 2; jj++) {
            if (!fy) jj++; //if flag is zero, do this loop only once
            for (int kk = 0; kk < 2; kk++) {
                if (!fz) kk++; //if flag is zero, do this loop only once
                //hout <<"a="<<a<<" fx="<<fx<<" b="<<b<<" fy="<<fy<<" c="<<c<<" fz="<<fz;
                int t = Calculate_t(temp[ii][0],temp[jj][1],temp[kk][2],sx,sy);
                //hout<<" t="<<t<<" sectioned_domain["<<t<<"].size()="<<sectioned_domain[t].size();
                sectioned_domain[t].push_back(P);
                //hout<<'.'<<endl;
            }
        }
    }
    return 1;
}
//Calculates the region to which a point corresponds
int Contact_grid::Calculate_t(const int &a, const int &b, const int &c, const int &sx, const int &sy)
{
    return a + b*sx + c*sx*sy;
}
//
int Contact_grid::Fill_sectioned_domain(const struct Geom_RVE &window_geom, const vector<int> &hybs_inside, const vector<vector<long int> > &structure, const vector<Point_3D> &points_in, const double &cutoff, const int &sx, const int &sy, const int &sz, vector<vector<int> > &sectioned_domain_hyb)
{
    //There will be sx*sy*sz different regions
    sectioned_domain_hyb.clear();
    vector<int> empty_int;
    sectioned_domain_hyb.assign(sx*sy*sz, empty_int);
    
    //First loop over the hybrids inside the box, then loop over the points inside each GNP
    for (int i = 0; i < (int)hybs_inside.size(); i++) {
        //hout<<"hybs_inside["<<i<<"]="<<cnts_inside[i]<<' ';
        int GNP = hybs_inside[i];
        //hout<<"structure["<<GNP<<"].size()="<<structure[GNP].size()<<endl;
        for (int j = 0; j < (int)structure[GNP].size(); j++) {
            long int P = structure[GNP][j];
            //hout<<"P=structure["<<GNP<<"]["<<j<<"]="<<structure[GNP][j]<< endl;
            
            //Calculate the region-coordinates
            int a, b, c;
            //hout << "points_in[P]=("<<points_in[P].x<<", "<<points_in[P].y<<", "<<points_in[P].z<<")"<< endl;
            if (!Calculate_region_coordinates(window_geom, points_in[P], sx, sy, sz, a, b, c)) {
                hout << "Error in Generate_sectioned_domain_cnts when calling Calculate_region_coordinates" << endl;
                return 0;
            }
            
            //Initialize flags for overlaping regions
            int fx = 0;
            int fy = 0;
            int fz = 0;
            
            //Assign value of flag according to position of point
            if (!Calculate_postion_flags(window_geom, points_in[P], cutoff, a, b, c, sx, sy, sz, fx, fy, fz)) {
                hout << "Error in Generate_sectioned_domain_cnts when calling Calculate_postion_flags" << endl;
                return 0;
            }
            
            //Create array for loop over overlaping regions
            int temp[2][3] = { {a+fx, b+fy, c+fz}, {a, b, c}};
            
            //In this loop I check all regions a point can belong to when it is in an overlaping zone
            for (int ii = 0; ii < 2; ii++) {
                if (!fx) ii++; //if flag is zero, do this loop only once
                for (int jj = 0; jj < 2; jj++) {
                    if (!fy) jj++; //if flag is zero, do this loop only once
                    for (int kk = 0; kk < 2; kk++) {
                        if (!fz) kk++; //if flag is zero, do this loop only once
                        //hout <<"a="<<a<<" fx="<<fx<<" b="<<b<<" fy="<<fy<<" c="<<c<<" fz="<<fz;
                        int t = Calculate_t(temp[ii][0],temp[jj][1],temp[kk][2],sx,sy);
                        //hout<<" t="<<t<<" sectioned_domain["<<t<<"].size()="<<sectioned_domain[t].size();
                        //Add the GNP number if the sectioned domain is empty or the last added GNP is not the current GNP
                        if (!sectioned_domain_hyb[t].size() || sectioned_domain_hyb[t].back() != GNP) {
                            sectioned_domain_hyb[t].push_back(GNP);
                        }
                        //hout<<'.'<<endl;
                    }
                }
            }
            
        }
    }
    
    return 1;
}