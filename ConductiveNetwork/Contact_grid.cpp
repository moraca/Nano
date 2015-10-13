//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Contact_grid.cpp
//OBJECTIVE:	Create a background grid to find the points in contact faster
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Contact_grid.h"

int Contact_grid::Generate_contact_grid(const vector<vector<long int> > &structure, const vector<int> &cnts_inside, vector<Point_3D> &points_in, const struct Geom_RVE &sample, struct Cutoff_dist &cutoffs, const struct Nanotube_Geo &cnts, double lx, double ly, double lz)
{
    //These variables are the coordinates of the lower corner of the RVE that defines its geometry
    double xmin = sample.origin.x;
    double ymin = sample.origin.y;
    double zmin = sample.origin.z;
    
    //Sizes of each region
    double dx = sample.gs_minx;
    double dy = sample.gs_miny;
    double dz = sample.gs_minz;
    
    //Number of regions on each direction
    int sx = (int)lx/dx;
    int sy = (int)ly/dy;
    int sz = (int)lz/dz;
    
    
    hout << "sx = " << sx << '\t' << "dx = " << dx << "\n";
    hout << "sy = " << sy << '\t' << "dy = " << dy << "\n";
    hout << "sz = " << sz << '\t' << "dz = " << dz  << "\n";
    
    //Maximum distance between two points in contact inside the sample
    double cutoff = cutoffs.tunneling_dist + 2*cnts.rad_max;
    
    //Check that the regions are not too small for the maximum cutoff distance 2r_max+tunnel
    //If they are, then change the number of sections to the maximum possible
    if (dx < 2*cutoff) {
        sx = (int)(lx/(2*(cutoff+Zero)));
        dx = lx/(double)sx;
        hout << "Modified the number of sections along x. " << "sx = " << sx << '\t' << "dx = " << dx << endl;
    }
    if (dy < 2*cutoff) {
        sy = (int)(ly/(2*(cutoff+Zero)));
        dy = ly/(double)sy;
        hout << "Modified the number of sections along y. " << "sy = " << sy << '\t' << "dy = " << dy << endl;
    }
    if (dz < 2*cutoff) {
        sz = (int)(lz/(2*(cutoff+Zero)));
        dz = lz/(double)sz;
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
    
    //First loop over the CNTs inside the box, then loop over the points inside each CNT
    for (int i = 0; i < (int)cnts_inside.size(); i++) {
        CNT = cnts_inside[i];
        for (int j = 0; j < (int)structure[CNT].size(); j++) {
            P = structure[CNT][j];
            //Save coordinates of the point
            x = points_in[P].x;
            y = points_in[P].y;
            z = points_in[P].z;
            
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
            else if ((x < lx+xmin-cutoff) && (x >= x2-cutoff) && (x <= x2 ))
                fx = 1;
            if ((y > cutoff + ymin) && (y >= y1) && (y <= y1+cutoff))
                fy = -1;
            else if ((y < ly+ymin-cutoff) && (y >= y2-cutoff) && (y <= y2 ))
                fy = 1;
            if ((z > cutoff + zmin) && (z >= z1) && (z <= z1+cutoff))
                fz = -1;
            else if ((z < lz+zmin-cutoff) && (z >= z2-cutoff) && (z <= z2 ))
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
                        t = calculate_t(temp[ii][0],temp[jj][1],temp[kk][2],sx,sy);
                        sectioned_domain[t].push_back(P);
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