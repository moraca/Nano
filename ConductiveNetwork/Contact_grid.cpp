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
        Vector with overlaping sub-regions. This vector is used to look for contact points faster
 
 */

int Contact_grid::Generate_contact_grid(const struct Geom_RVE &sample, struct Cutoff_dist &cutoffs, const struct Nanotube_Geo &cnts, const vector<int> &cnts_inside, const vector<vector<long int> > &structure, vector<Point_3D> &points_in, int window)
{
    //Dimensions of the current observation window
    double w_x = sample.win_max_x - window*sample.win_delt_x;
    double w_y = sample.win_max_y - window*sample.win_delt_y;
    double w_z = sample.win_max_z - window*sample.win_delt_z;

    //These variables are the coordinates of the lower corner of the observation window
    double xmin = sample.origin.x + (sample.len_x - w_x)/2;
    double ymin = sample.origin.y + (sample.wid_y - w_y)/2;
    double zmin = sample.origin.z + (sample.hei_z - w_z)/2;
    
    //Sizes of each region
    double dx = sample.gs_minx;
    double dy = sample.gs_miny;
    double dz = sample.gs_minz;
    
    //Number of regions on each direction
    int sx = (int)(w_x/dx);
    int sy = (int)(w_y/dy);
    int sz = (int)(w_z/dz);
    
    
    //hout << "sx = " << sx << '\t' << "dx = " << dx << "\n";
    //hout << "sy = " << sy << '\t' << "dy = " << dy << "\n";
    //hout << "sz = " << sz << '\t' << "dz = " << dz  << "\n";
    
    //Maximum distance between two points in contact inside the sample
    double cutoff = cutoffs.tunneling_dist + 2*cnts.rad_max;
    
    //Check that the regions are not too small for the maximum cutoff distance 2r_max+tunnel
    //If they are, then change the number of sections to the maximum possible
    if (dx < 2*cutoff) {
        sx = (int)(w_x/(2*(cutoff+Zero)));
        dx = w_x/(double)sx;
        hout << "Modified the number of sections along x. " << "sx = " << sx << '\t' << "dx = " << dx << endl;
    }
    if (dy < 2*cutoff) {
        sy = (int)(w_y/(2*(cutoff+Zero)));
        dy = w_y/(double)sy;
        hout << "Modified the number of sections along y. " << "sy = " << sy << '\t' << "dy = " << dy << endl;
    }
    if (dz < 2*cutoff) {
        sz = (int)(w_z/(2*(cutoff+Zero)));
        dz = w_z/(double)sz;
        hout << "Modified the number of sections along z. " << "sz = " << sz << '\t' << "dz = " << dz  << endl;
    }
    
    
    //These variables will give me the region cordinates of the region that a point belongs to
    int a, b, c;
    int t;
    
    //These variables are to reduce operations when accesing the coordinates of a point and it's CNT number
    double x, y, z;
    int CNT;
    long int P;
    
    //There will be sx*sy*sz different regions
    sectioned_domain.clear();
    vector<long int> empty_long;
    sectioned_domain.assign(sx*sy*sz, empty_long);
    hout<<"There are "<<sx*sy*sz<<" overlapping sub-regions."<<endl;
    
    //First loop over the CNTs inside the box, then loop over the points inside each CNT
    for (int i = 0; i < (int)cnts_inside.size(); i++) {
        //hout<<"cnts_inside["<<i<<"]="<<cnts_inside[i]<<' ';
        CNT = cnts_inside[i];
        //hout<<"structure["<<CNT<<"].size()="<<structure[CNT].size()<<endl;
        for (int j = 0; j < (int)structure[CNT].size(); j++) {
            P = structure[CNT][j];
            //hout<<"P=structure["<<CNT<<"]["<<j<<"]="<<structure[CNT][j]<<' ';
            //Save coordinates of the point
            x = points_in[P].x;
            y = points_in[P].y;
            z = points_in[P].z;
            //hout<<"x="<<x<<" y="<<y<<" z="<<z<<endl;
            
            //Calculate the region-coordinates
            a = (int)((x-xmin)/dx);
            //Limit the value of a as it has to go from 0 to sx-1
            if (a == sx) {
                a--;
                //hout << " a-- x="<<x<<" xmin="<<xmin<<" dx="<<dx<<" a="<<((x-xmin)/dx)<<" sx="<<sx<<endl;
            }
            b = (int)((y-ymin)/dy);
            //Limit the value of b as it has to go from 0 to sy-1
            if (b == sy) {
                b--;
                //hout << " b-- y="<<y<<" ymin="<<ymin<<" dy="<<dy<<" b="<<((y-ymin)/dy)<<" sy="<<sy<<endl;
            }
            c = (int)((z-zmin)/dz);
            //Limit the value of c as it has to go from 0 to sz-1
            if (c == sz){
                c--;
                //hout << " c-- z="<<z<<" zmin="<<zmin<<" dz="<<dz<<" c="<<((z-zmin)/dz)<<" sz="<<sz<<endl;
            }
            
            //Coordinates of non-overlaping region the point belongs to
            double x1 = a*dx +  xmin;
            double x2 = x1 + dx;
            double y1 = b*dy +  ymin;
            double y2 = y1 + dy;
            double z1 = c*dz +  zmin;
            double z2 = z1 + dz;
            
            //Initialize flags for overlaping regions
            int fx = 0;
            int fy = 0;
            int fz = 0;
            
            //Assign value of flag according to position of point
            //The first operand eliminates the periodicity on the boundary
            if ((x > cutoff + xmin) && (x >= x1) && (x <= x1+cutoff))
                fx = -1;
            else if ((x < w_x+xmin-cutoff) && (x >= x2-cutoff) && (x <= x2 ))
                fx = 1;
            if ((y > cutoff + ymin) && (y >= y1) && (y <= y1+cutoff))
                fy = -1;
            else if ((y < w_y+ymin-cutoff) && (y >= y2-cutoff) && (y <= y2 ))
                fy = 1;
            if ((z > cutoff + zmin) && (z >= z1) && (z <= z1+cutoff))
                fz = -1;
            else if ((z < w_z+zmin-cutoff) && (z >= z2-cutoff) && (z <= z2 ))
                fz = 1;
            
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
                        t = calculate_t(temp[ii][0],temp[jj][1],temp[kk][2],sx,sy);
                        //hout<<" t="<<t<<" sectioned_domain["<<t<<"].size()="<<sectioned_domain[t].size();
                        sectioned_domain[t].push_back(P);
                        //hout<<'.'<<endl;
                    }
                }
            }
        }
    }
    
    return 1;
}

//Calculates the region to which a point corresponds
int Contact_grid::calculate_t(int a, int b, int c, int sx, int sy)
{
    return a + b*sx + c*sx*sy;
}