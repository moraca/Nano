//
//  RNetwork.cpp
//  Nanocode
//
//  Created by Angel Mora Cordova on 3/24/14.
//  Copyright (c) 2014 Angel Mora. All rights reserved.
//

/* v16
 Improved the new way to create clusters. The improvemente in time is around 5%. Not so ggod but is something.
 
 //*/
#include "RNetwork.h"


int RNetwork::Construct(ifstream &infile, string structure_type, const int &samples_count, struct CNT_Geo cnts_geo, vector<Point_3D> points_in, const struct RVE_Geo &cell_geo, struct Region_Geo cnt_regions,vector<double> cnts_radius, vector<vector<long int> > structure, vector<vector<int> > sectioned_domain_cnt)
{
    //Initial number of points
    points0 = points_in.size();
    
    //To check the time of each function
    clock_t ct_begin, ct_end;
    
    //Read the parameters for the function
    ct_begin = clock();
    if(!Read_parameters(infile, cell_geo, cnts_geo)) return 0;
    ct_end = clock();
    hout << "Read_parameters " << (double)(ct_end-ct_begin)/CLOCKS_PER_SEC << " secs " << endl;
    
    if (structure_type == "structure"){
        //Write input files so they can be used when using the option "file"
        ct_begin = clock();
        if(!Print_input_data_files(sectioned_domain_cnt, points_in, cnts_radius, cell_geo, cnts_geo, cnt_regions)) return 0;
        ct_end = clock();
        hout << "Print_input_data_files  " << (double)(ct_end-ct_begin)/CLOCKS_PER_SEC << " secs " << endl;
    } else if(structure_type != "file") {
        hout << "Invalid input. Only 'structure' and 'file' are accepted values. Inpur string: " << structure_type << endl;
        return 0;
    }
    
    //initialize the size of the box with the largest possible value
    double window = window0;
    
    //Variables to use the command line
    int s, it = 0;
    char command[100];
    
    //loop over the size of the box
    while (window <= window_max) {
        hout << "iteration " << it << " ==========================================" << endl;
        
        ct_begin = clock();
        if (!Analysis(window, cell_geo, cnts_geo, cnt_regions, points_in, cnts_radius, sectioned_domain_cnt, structure)) return 0;
        ct_end = clock();
        hout << "Analysis  " << (double)(ct_end-ct_begin)/CLOCKS_PER_SEC << " secs " << endl;
        
        //Move the visualization files to a new folder
        s = sprintf(command, "mkdir iter_%.4d", it);
        system(command);
        s = sprintf(command, "mv Single*.dat iter_%.4d", it);
        system(command);//*/
        it++;
        
        //update the size of the box
        window = window + dwindow;
        //update box geomtry
        Update_box_geometry();
        
        //Send to output file the data for the restart size of the box:
        hout << "Restart from: " << window_max << ' ' << dwindow << ' ' << window <<endl;
    }
    
    
    
    return 1;
}

//Read necesary input parameters from input file
int RNetwork::Read_parameters(ifstream &infile, const struct RVE_Geo cell_geo, CNT_Geo cnts_geo)
{
    
    // Read the cutoff for tunneling
    istringstream istr_tunnel(Get_Line(infile));
    istr_tunnel >> tunnel;
    
    //Calculate the maximum cutoff distance for tunneling
    cutoff = tunnel + 2*cnts_geo.rad_max;
    
    //First read the sizes of the regions
    istringstream istr_regions(Get_Line(infile));
    istr_regions >> secx >> secy >> secz;
    
    //Read the corner point and the size of the box inside the sample
    //istringstream istr_box(Get_Line(infile));
    //istr_box >> box_geometry.poi_min.x >> box_geometry.poi_min.y >> box_geometry.poi_min.z;
    //istr_box >> box_geometry.len_x >> box_geometry.wid_y >> box_geometry.hei_z;
    
    //Get the value of the initial observation window
    istringstream istr_window(Get_Line(infile));
    istr_window >> window_max >> dwindow >> window0;
    
    //Assign the size of the box. The maximum size is used, as the iterative process will start from the largest sample size
    //to the smallest.
    //window_max is the size of the sample is the simulation starts from the generation of the structure or it is the
    //restart size if data is read from a file
    box_geometry.len_x = window0;
    box_geometry.wid_y = window0;
    box_geometry.hei_z = window0;
    //The reference point of the box will be calculated using the reference point of the cell. As the observation window will
    //always be centered there is no need to specify a new point.
    box_geometry.poi_min.x = cell_geo.poi_min.x + (cell_geo.len_x - window0)/2;
    box_geometry.poi_min.y = cell_geo.poi_min.y + (cell_geo.wid_y - window0)/2;
    box_geometry.poi_min.z = cell_geo.poi_min.z + (cell_geo.hei_z - window0)/2;
    //These variables are to store the size of the RVE and reduces operations when accessing them
    lx = box_geometry.len_x;
    ly = box_geometry.wid_y;
    lz = box_geometry.hei_z;
    //These variables are the coordinates of the lower corner of the RVE that defines its geometry
    xmin = box_geometry.poi_min.x;
    ymin = box_geometry.poi_min.y;
    zmin = box_geometry.poi_min.z;
    //These variables are the coordinates of the upper corner of the RVE that defines its geometry
    xmax = lx + xmin;
    ymax = ly + ymin;
    zmax = lz + zmin;
    
    //Check the geometry of the box is correct
    if ((lx>cell_geo.len_x)||(ly>cell_geo.wid_y)||(lz>cell_geo.hei_z)) {
        hout << "Box size is bigger than RVE"<<endl;
        return 0;
    } else if ((lx < 0)||(ly < 0)||(lz < 0)) {
        hout << "Negative size found. Box cannot have negative size"<<endl;
        return 0;
    } else if ((xmin<cell_geo.poi_min.x)||(ymin<cell_geo.poi_min.y)||(zmin<cell_geo.poi_min.z)) {
        hout << "Minimum point of box is outside RVE"<<endl;
        return 0;
    } else if ((lx+xmin>cell_geo.len_x)||(ly+ymin>cell_geo.wid_y)||(lz+zmin>cell_geo.hei_z)) {
        hout << "Given the minimum point and box dimensions, box goes outside the RVE"<<endl;
        return 0;
    }
    
    /*/ Check that the regions are not too small for the maximum cutoff distance 2r_max+tunnel
     if ((lx/(double)secx) < 2*cutoff) {
     hout << "The regions along x are too many for the cutoff distance for tunneling." <<  endl;
     hout << "The length of the region along each direction has to be at least twice the cutoff distance" << endl;
     hout << "2cutoff = 2(2r_max+tunnel) = " << 2*cutoff << ", size of region along x = " << (lx/(double)secx) << endl;
     return 0;
     }
     if ((ly/(double)secy) < 2*cutoff) {
     hout << "The regions along y are too many for the cutoff distance for tunneling." <<  endl;
     hout << "The length of the region along each direction has to be at least twice the cutoff distance" << endl;
     hout << "2cutoff = 2(2r_max+tunnel) = " << 2*cutoff << ", size of region along y = " << (ly/(double)secy) << endl;
     return 0;
     }
     if ((lz/(double)secz) < 2*cutoff) {
     hout << "The regions along z are too many for the cutoff distance for tunneling." <<  endl;
     hout << "The length of the region along each direction has to be at least twice the cutoff distance" << endl;
     hout << "2cutoff = 2(2r_max+tunnel) = " << 2*cutoff << ", size of region along z = " << (lz/(double)secz) << endl;
     return 0;
     }//*/
    
    //Get the magnitude of voltage between these two directions
    istringstream istr_magnitude(Get_Line(infile));
    istr_magnitude >> magnitude;
    
    //Get the value of resistivity of carbon fiber/CNT
    istringstream istr_resistivity(Get_Line(infile));
    istr_resistivity >> resistivity;
    
    return 1;
}

//Save (print) the input data into a file
int RNetwork::Print_input_data_files(vector<vector<int> > sectioned_domain_cnt, vector<Point_3D> point_list, vector<double> cnt_radius, RVE_Geo cell_geo, CNT_Geo cnts_geo, Region_Geo cnt_regions)
{
    //Write the list of points to a file
    ofstream points("point_list_out.txt");
    //hout << "the number of points generated is " << (int)point_list.size() <<endl;
    for (long int i=0; i < point_list.size(); i++) {
        points << point_list[i].x << '\t' << point_list[i].y << '\t' ;
        points << point_list[i].z << '\t' << point_list[i].flag << '\n';
    }
    points.close();
    
    //Compress the points file and delete the original one (to save space on disk)
    system("tar -zcvf point_list_out.txt.tar.gz point_list_out.txt; rm point_list_out.txt");
    
    //Write the list of CNT radii to a file
    ofstream radii("radii_out.txt");
    //hout << "the number of radius " << (int)cnt_radius.size() <<endl;
    for (int i = 0; i < (int)cnt_radius.size(); i++) {
        radii << cnt_radius[i] << '\n' ;
    }
    radii.close();
    
    //Write the data of the geometry of the RVE to a file
    ofstream rve("rve_out.txt");
    rve << cell_geo.poi_min.x << '\t' ;
    rve << cell_geo.poi_min.y << '\t' ;
    rve << cell_geo.poi_min.z << '\t' ;
    rve << cell_geo.poi_min.flag << '\n' ;
    rve << cell_geo.len_x << '\t' ;
    rve << cell_geo.wid_y << '\t' ;
    rve << cell_geo.hei_z << '\n' ;
    rve << cell_geo.volume << '\n' ;
    //These if-statements are placed because in tests a small number is stored in the variables of O(10^-310)
    //which is way below e_machine. It caused some errors when using the option to read the data
    //from files. Somehow after reading these small values, reading the following values gave a wrong number
    //For instance, the type (CNT/CF) did not appear and more values were O(10^-310).
    //With this the cutoff for tunneling was practically zero, so all CNTs were maked as isolated.
    //Probably these are just unused variables.
    if (cell_geo.density < Zero)
        rve << 0 << '\n' ;
    else
        rve << cell_geo.density << '\n' ;
    if (cell_geo.delt_x < Zero)
        rve << 0 << '\t' ;
    else
        rve << cell_geo.delt_x << '\t' ;
    if (cell_geo.delt_y < Zero)
        rve << 0 << '\t' ;
    else
        rve << cell_geo.delt_y << '\t' ;
    if (cell_geo.delt_z < Zero)
        rve << 0;
    else
        rve << cell_geo.delt_z ;
    rve.close();
    
    
    //Write the data of the geometry of the CNTs to a file
    //The tunnel is written along with the CNT geometry just so I don't write a file with only one number
    ofstream cnt("cnt_out.txt");
    cnt << cnts_geo.criterion << '\n' ;
    cnt << cnts_geo.step_length << '\n' ;
    cnt << cnts_geo.len_dist_type << '\t' ;
    cnt << cnts_geo.len_min << '\t' ;
    cnt << cnts_geo.len_max << '\n' ;
    cnt << cnts_geo.rad_dist_type << '\t' ;
    cnt << cnts_geo.rad_min << '\t' ;
    cnt << cnts_geo.rad_max << '\n' ;
    cnt << cnts_geo.dir_dist_type << '\t' ;
    //These if-statements are placed because in tests a small number is stored in the variables of O(10^-310)
    //which is way below e_machine. It caused some errors when using the option to read the data
    //from files. Somehow after reading these small values, reading the following values gave a wrong number
    //For instance, the type (CNT/CF) did not appear and more values we O(10^-310).
    //With this the cutoff for tunneling was practically zero, so all CNTs were maked as isolated.
    //Probably these are just unused variables.
    if (cnts_geo.ini_sita < Zero)
        cnt << 0 << '\t' ;
    else
        cnt << cnts_geo.ini_sita << '\t' ;
    if (cnts_geo.ini_pha < Zero)
        cnt << 0 << '\t' ;
    else
        cnt << cnts_geo.ini_pha << '\t' ;
    cnt << cnts_geo.angle_max << '\n' ;
    cnt << cnts_geo.volume_fraction << '\t' ;
    cnt << cnts_geo.real_volume << '\t' ;
    cnt << cnts_geo.weight_fraction << '\t' ;
    cnt << cnts_geo.real_weight << '\t' ;
    cnt << cnts_geo.linear_density << '\n' ;
    cnt << cnts_geo.type << '\n' ;
    cnt.close();
    
    //Write the data of the geometry of the CNT regions
    ofstream cnt_region("cnt_regions_out.txt");
    cnt_region << cnt_regions.secx << '\t' << cnt_regions.secy << '\t' << cnt_regions.secz << '\n';
    cnt_region << cnt_regions.lx << '\t' << cnt_regions.ly << '\t' << cnt_regions.lz ;
    cnt_region.close();
    
    //Write the list of CNT regions to a file
    ofstream sectioned_domain("sectioned_domain_out.txt");
    for (long int i = 0; i < sectioned_domain_cnt.size(); i++){
        sectioned_domain << sectioned_domain_cnt[i].size() << '\t';
        for (long int j = 0; j < sectioned_domain_cnt[i].size(); j++) {
            sectioned_domain << sectioned_domain_cnt[i][j] << '\t' ;
        }
        sectioned_domain << '\n';
    }
    sectioned_domain << -1;
    sectioned_domain.close();
    
    return 1;
}

//This function groups some major functions so that the proces of finding the RVE can be done iteratively
int RNetwork::Analysis(double window, const struct RVE_Geo cell_geo, struct CNT_Geo cnts_geo, struct Region_Geo cnt_regions, vector<Point_3D> points_in, vector<double> cnts_radius, vector<vector<int> > sectioned_domain_cnt, vector<vector<long int> > structure)
{
    //To check the time of each function
    clock_t ct_begin, ct_end;
    
    ct_begin = clock();
    if (!Region_of_interest(points_in, cnt_regions, structure, sectioned_domain_cnt, cnts_radius)) return 0;
    ct_end = clock();
    hout << "Region_of_interest " << (double)(ct_end-ct_begin)/CLOCKS_PER_SEC << " secs " << endl;
    
    ct_begin = clock();
    if(!Assign_region(box_geometry, points_in, structure)) return 0;
    ct_end = clock();
    hout << "Assign_region " << (double)(ct_end-ct_begin)/CLOCKS_PER_SEC << " secs " << endl;
    
    ct_begin = clock();
    if(!Check_contacts(points_in, cnts_radius, tunnel, structure.size())) return 0;
    ct_end = clock();
    hout << "Check_contacts " << (double)(ct_end-ct_begin)/CLOCKS_PER_SEC << " secs " << endl;
    
    ct_begin = clock();
    if(!Make_CNT_clusters(cnts_geo, points_in, structure)) return 0;
    ct_end = clock();
    hout << "Make_CNTC_lusters " << (double)(ct_end-ct_begin)/CLOCKS_PER_SEC << " secs " << endl;
    
    ct_begin = clock();
    if(Check_clusters_percolation(points_in, structure)) {
        /*string filename;
         vector<int> clust(1, 0);
         for (int i = 0; i < clusters_cnt[0].size(); i++) {
         clust.front() = clusters_cnt[0][i];
         ostringstream number;
         number << clusters_cnt[0][i];
         filename = filename.append("CNT");
         filename = filename.append(number.str());
         filename = filename.append(".dat");
         Export_cnt_networks_meshes(cell_geo, structure, clust, points_in, cnts_radius, filename);
         filename.clear();
         }*/
        ct_end = clock();
        hout << "Check_clusters_percolation " << (double)(ct_end-ct_begin)/CLOCKS_PER_SEC << " secs " << endl;
        
        ct_begin = clock();
        if(!Split_cnts(points_in, cnts_radius, structure)) return 0;
        ct_end = clock();
        hout << "Split_cnts " << (double)(ct_end-ct_begin)/CLOCKS_PER_SEC << " secs " << endl;
    } else{
        ct_end = clock();
        hout << "Check_clusters_percolation " << (double)(ct_end-ct_begin)/CLOCKS_PER_SEC << " secs " << endl;
    }
    
    ct_begin = clock();
    if (!Clusters_length(points_in, structure)) return 0;
    ct_end = clock();
    hout << "Clusters_length " << (double)(ct_end-ct_begin)/CLOCKS_PER_SEC << " secs " << endl;
    
    ct_begin = clock();
    //if (!Export_visualization_files(cell_geo, points_in, cnts_radius, structure)) return 0;
    if (!Export_visualization_files(box_geometry, points_in, cnts_radius, structure)) return 0;
    ct_end = clock();
    hout << "Export_visualization_files " << (double)(ct_end-ct_begin)/CLOCKS_PER_SEC << " secs " << endl;//*/
    
    Clear_vectors();
    
    return 1;
}

//Make the clusters according to the box size
int RNetwork::Region_of_interest(vector<Point_3D> &points_in, struct Region_Geo cnt_regions, vector<vector<long int> > &structure, vector<vector<int> > sectioned_domain_cnt, vector<double> &radii)
{
    //Check if the maximum size of the box is equal to the whole sample
    cnts_inside_flags.assign(structure.size(), '0');
    //Print1DVec(cnts_inside_flags, "cnts_inside_flags.txt");
    //Identify the regions where the boundaries of the box are located
    int sx0, sy0, sz0, sx1, sy1, sz1;
    
    //hout << "5 ";
    //Gather all CNTs in the boundary
    Get_boundary_cnts(sx0, sx1, sy0, sy1, sz0, sz1, cnt_regions, sectioned_domain_cnt);
    
    //hout << "6 ";
    //Scan every Nanotube that is the boundary region. Delete and trim CNTs when needed.
    if (!Locate_and_trim_boundary_cnts(points_in, structure, radii)){
        hout << "Error in Locate_and_trim_boundary_cnts (initial)" << endl;
        return 0;
    }
    
    //hout << "7 ";
    //Remove CNTs that are empty
    long int t;
    int CNT;
    for (int i = (int)cnts_inside.size()-1; i >= 0; i--) {
        CNT = cnts_inside[i];
        if (!structure[CNT].size()) {
            cnts_inside.erase(cnts_inside.begin()+i);
        }
    }
    
    //hout << "8 ";
    //Add the rest of CNTs inside the box
    for (int i = sx0+1; i < sx1 ; i++) {
        for (int j = sy0+1; j < sy1; j++) {
            for (int k = sz0+1; k < sz1; k++) {
                t = calculate_t(i, j, k, cnt_regions.secx, cnt_regions.secy);
                //Add all CNTs from the current region
                Update_cnts_inside_flags(sectioned_domain_cnt, t);
            }
        }
    }//*/
    
    //hout << "9 ";
    
    return 1;
}

void RNetwork::Get_boundary_cnts(int &sx0, int &sx1, int &sy0, int &sy1, int &sz0, int &sz1, struct Region_Geo cnt_regions, vector<vector<int> > sectioned_domain_cnt)
{
    //hout << "1 ";
    sx0 = (int)(box_geometry.poi_min.x/cnt_regions.lx);
    if (sx0 == cnt_regions.secx) sx0--;
    sx1 = (int)((box_geometry.poi_min.x + box_geometry.len_x)/cnt_regions.lx);
    if (sx1 == cnt_regions.secx) sx1--;
    
    sy0 = (int)(box_geometry.poi_min.y/cnt_regions.ly);
    if (sy0 == cnt_regions.secy) sy0--;
    sy1 = (int)((box_geometry.poi_min.y + box_geometry.wid_y)/cnt_regions.ly);
    if (sy1 == cnt_regions.secy) sy1--;
    
    sz0 = (int)(box_geometry.poi_min.z/cnt_regions.lz);
    if (sz0 == cnt_regions.secz) sz0--;
    sz1 = (int)((box_geometry.poi_min.z + box_geometry.hei_z)/cnt_regions.lz);
    if (sz1 == cnt_regions.secz) sz1--;
    
    //Gather all CNTs on the boundary regions
    long int t;
    /*hout << "2 ";
     hout << "s..." << sx0 << ' ' << sx1 << ' ' << sy0 << ' ' << sy1 << ' ' << sz0 << ' ' << sz1 << ' ';
     hout << "box..." << box_geometry.poi_min.x << ' ' << box_geometry.poi_min.y << ' ' << box_geometry.poi_min.z << ' ';
     hout << "cnts..." << cnt_regions.lx << ' ' << cnt_regions.ly << ' ' << cnt_regions.lz << ' ';//*/
    //Fixed sx:
    for (int j = sy0; j <= sy1; j++) {
        for (int k = sz0; k <= sz1; k++) {
            t = calculate_t((long int)sx0, (long int)j, (long int)k, cnt_regions.secx, cnt_regions.secy);
            //Add all CNTs from the current region
            Update_cnts_inside_flags(sectioned_domain_cnt, t);
            t = calculate_t((long int)sx1, (long int)j, (long int)k, cnt_regions.secx, cnt_regions.secy);
            //Add all CNTs from the current region
            Update_cnts_inside_flags(sectioned_domain_cnt, t);
        }
    }
    //hout << "3 ";
    //Fixed sy
    for (int i = sx0; i <=sx1 ; i++) {
        for (int k = sz0; k <= sz1; k++) {
            t = calculate_t((long int)i, (long int)sy0, (long int)k, cnt_regions.secx, cnt_regions.secy);
            //Add all CNTs from the current region
            Update_cnts_inside_flags(sectioned_domain_cnt, t);
            t = calculate_t((long int)i, (long int)sy1, (long int)k, cnt_regions.secx, cnt_regions.secy);
            //Add all CNTs from the current region
            Update_cnts_inside_flags(sectioned_domain_cnt, t);
        }
    }
    //hout << "4 ";
    //Fixed sz
    for (int i = sx0; i <=sx1 ; i++) {
        for (int j = sy0; j <= sy1; j++) {
            t = calculate_t((long int)i, (long int)j, (long int)sz0, cnt_regions.secx, cnt_regions.secy);
            //Add all CNTs from the current region
            Update_cnts_inside_flags(sectioned_domain_cnt, t);
            t = calculate_t((long int)i, (long int)j, (long int)sz1, cnt_regions.secx, cnt_regions.secy);
            //Add all CNTs from the current region
            Update_cnts_inside_flags(sectioned_domain_cnt, t);
        }
    }//*/
}

//
void RNetwork::Update_cnts_inside_flags(vector<vector<int> > sectioned_domain_cnt, long int t)
{
    //Add all CNTs from the current region
    int CNT;
    for (int i = 0; i < (int)sectioned_domain_cnt[t].size(); i++) {
        CNT = sectioned_domain_cnt[t][i];
        if (cnts_inside_flags[CNT] == '0') {
            cnts_inside.push_back(CNT);
            cnts_inside_flags[CNT] = '1';
        }
    }
    
}

//This Functions deletes the repeated items added to the vector vec
void RNetwork::Discard_repeated(vector<int> &vec)
{
    for (int i = 0; i <  (int) vec.size()-1; i++) {
        for (int k = i+1; k < (int) vec.size(); k++) {
            if (vec[i] == vec[k]){
                vec.erase(vec.begin() + k);
                k--;
            }
        }
    }
}//*/

//This Functions deletes the repeated items added to the vector vec
void RNetwork::Discard_repeated(vector<long int> &vec)
{
    for (long int i = 0; i <  vec.size()-1; i++) {
        for (long int k = i+1; k < vec.size(); k++) {
            if (vec[i] == vec[k]){
                vec.erase(vec.begin() + k);
                k--;
            }
        }
    }
}//*/

//
int RNetwork::Locate_and_trim_boundary_cnts(vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<double> &radii)
{
    //These variables will help me locate the point with respect with the box
    string currentPoint, nextPoint;
    //Variables for current and next points and current CNT
    long int P1, P2;
    int CNT;
    //Empty vector to increase size of other vectors
    vector<int> empty;
    //Resize the vector of boundary_flags with empty vector
    boundary_flags.assign(points_in.size(), empty);
    //Empty vector to increase size of other vectors
    vector<short int> empty_short;
    //Resize the vector of boundary_flags_cnt with empty vector
    boundary_flags_cnt.assign(structure.size(), empty_short);
    //hout << "cnts_inside.size() = "<<cnts_inside.size()<<endl;
    for (long int i = 0; i < cnts_inside.size(); i++) {
        CNT = cnts_inside[i];
        //hout << "Check0 " <<  CNT << ' ' << structure[CNT].size() << ' ' << endl;
        //hout << " all=" << structure.size() << " contacts_cnt_point=" << contacts_cnt_point.size() << ' ' ;
        P1 = structure[CNT][0];
        currentPoint = Where_is(points_in[P1].x, points_in[P1].y, points_in[P1].z);
        //hout << "Check0.1 " ;
        for (long int j = 1; j <= structure[CNT].size(); j++) {
            //hout << "Check1 " << i << ' ' << j << " CNT=" << CNT << ' ' << "P1=" << P1 << " s[CNT].s=" << structure[CNT].size() << ' ' << currentPoint << ' ';
            //Handle the last point
            if (j == structure[CNT].size()) {
                nextPoint = "nothing";
                //hout << " nothing ";
            } else {
                P2 = structure[CNT][j];
                nextPoint = Where_is(points_in[P2].x, points_in[P2].y, points_in[P2].z);
                //hout << nextPoint << ' ';
            }
            //hout << " (" << points_in[P1].x << ", " << points_in[P1].y << ", " << points_in[P1].z << ") ";
            //hout << "Check2 " << endl;//*/
            
            //Check where is the current point and where the next point
            //First I check if the point is inside the box. If it is, what follows makes that assumption
            if (currentPoint == "inside") {
                if (nextPoint == "outside") {
                    //hout << "Check3 ";
                    //Make the projection and add the new point to the point_list and the the structure vector
                    //j-1 is the current point
                    //j is the next point
                    if (!New_boundary_point(points_in, structure, j-1, j, CNT, currentPoint)) {
                        hout  << "Could not make projection on the boundary" << endl;
                        return 0;
                    }
                    //Trim the CNT from the projected boundary point, which now is in position j
                    Trim_CNT(points_in, structure, cnts_inside, radii, j, CNT);
                    //Now the position of nextPoint is for the boundary point, so I need to update nextPoint and P2
                    nextPoint = "boundary";
                    P2 = points_in.size()-1;
                    //hout << "Check5 ";
                }
            } else if (currentPoint == "outside") {
                if (nextPoint == "inside") {
                    //hout << "Check6 ";
                    //When the next point is inside, take the intersection and add the new point
                    //to the point_list and to the structure vector
                    if (!New_boundary_point(points_in, structure, j, j-1, CNT, currentPoint)) {
                        hout  << "Could not make projection on the boundary" << endl;
                        return 0;
                    }
                    //Now the position of nextPoint is for the boundary point, so I need to update nextPoint
                    nextPoint = "boundary";
                    P2 = points_in.size()-1;
                    //hout << "Check7 ";
                }
                //hout << "Check11 ";
                //Remove the current point from the structure vector
                structure[CNT].erase(structure[CNT].begin()+j-1);
                //Update the iterator j. This is needed beacuase one element was deleted, so all positions shift one place
                j--;
                //hout << "Check12 ";
            } else if (currentPoint == "boundary") {
                //hout << "Check13 ";
                //If the next point is outside, then the CNT might need to be trimmed
                if (nextPoint == "outside") {
                    //hout << "Check14 ";
                    //If the boundary point is the first point, then this point needs to be deleted and
                    //treated as if it was outside
                    if (j==1){
                        structure[CNT].erase(structure[CNT].begin());
                        j--;
                    } else {
                        //If the boundary is not the first point then trim the CNT
                        Trim_CNT(points_in, structure, cnts_inside, radii, j-1, CNT);
                        //Add the current point to the corresponding boundary vector
                        long int point_number = structure[CNT][j-1];
                        Add_to_boundary_vectors(points_in[point_number], point_number);
                    }
                } else if ((nextPoint == "nothing")&&(j==1)) {
                    //If this is the last and only point, then delete it
                    structure[CNT].erase(structure[CNT].begin());
                    j--;
                } else {
                    //Add the current point to the corresponding boundary vector
                    long int point_number = structure[CNT][j-1];
                    Add_to_boundary_vectors(points_in[point_number], point_number);
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

string RNetwork::Where_is(double x, double y, double z)
{
    //If the point is outside the box then any these conditions needs to be true
    if ((x < xmin)||(x > xmin+lx)||(y < ymin)||(y > ymin+ly)||(z < zmin)||(z > zmin+lz))
        return "outside";
    //If the point it's at a boundary of the box, then only one of these conditions needs to be true,
    //provided that it is not outside
    else if ((x == xmin)||(x == xmin+lx)||(y == ymin)||(y == ymin+ly)||(z == zmin)||(z == zmin+lz))
        return "boundary";
    //If the point is not outside nor at the boundary, then it's inside
    else
        return "inside";//*/
    
    /*/============================
     //If the point is outside the box then any these conditions needs to be true
     if ((x < xmin)||(x > xmin+lx)||(y < ymin)||(y > ymin+ly)||(z < zmin)||(z > zmin+lz))
     return "outside";
     //If the point it's at a boundary of the box, then only one of these conditions needs to be true,
     //provided that it is not outside
     else if (((x - xmin) < Zero)||((x - (xmin+lx)) < Zero)||((y - ymin) < Zero)||((y - (ymin+ly)) < Zero)||((z - zmin) < Zero)||((z - (zmin+lz)) < Zero))
     return "boundary";
     //If the point is not outside nor at the boundary, then it's inside
     else
     return "inside";//*/
}

//This function created the projection of the pint on the boundary of the box, and adds it to the
//point vectors and updates the total number of points
//The new point is added at the end on the point_list. Otherwise, if it's added between the two points
//of the CNT it belongs to, then structure and contacst should also be modified. Changing structure and update
//all other variables will be particularly challenging
int RNetwork::New_boundary_point(vector<Point_3D> &points_in, vector<vector<long int> > &structure, long int insidePoint, long int outsidePoint, int CNT, string currentLocation)
{
    //Variables for the parametrization of the line btween two consecutive points
    double x, y, z, x1, y1, z1;
    
    //Proyected variables
    double xp, yp ,zp;
    
    //To make things easier, the inside point will be x,y,z and the outside point x1,y1,z1
    //Global coordinate
    //Coordinates of the inside point
    long int global_i = structure[CNT][insidePoint];
    x = points_in[global_i].x;
    y = points_in[global_i].y;
    z = points_in[global_i].z;
    //Coordinates of the outside point
    long int global_o = structure[CNT][outsidePoint];
    x1 = points_in[global_o].x;
    y1 = points_in[global_o].y;
    z1 = points_in[global_o].z;
    
    string point_loc=Where_is(x1,y1,z1);
    
    while(point_loc!="boundary"){
        if(!Projected_Point(x,y,z,x1,y1,z1,xp,yp,zp)) {
            hout<< "error in getting projected point" <<endl;
            return 0;
        }
        point_loc = Where_is(xp,yp,zp);
        x1=xp; y1=yp; z1=zp;
    }
    
    
    //Projected point
    Point_3D projection(xp ,yp, zp); //Coordinates
    projection.flag = CNT; //CNT number
    
    //Add the projected point to the end of the point_list
    points_in.push_back(projection);
    
    //Update total number of points
    long int new_P = points_in.size();
    
    //Increase the size of the boundary_flags vector
    vector<int> empty;
    boundary_flags.push_back(empty);
    
    //Add the new point to the CNT (in structure vector)
    if (currentLocation == "inside") {
        //If the current point is inside, add the new point after the inside point
        structure[CNT].insert(structure[CNT].begin()+insidePoint+1,new_P-1);
    } else if (currentLocation == "outside") {
        //If the corrent point is outside, then the new point is located after the outside point
        structure[CNT].insert(structure[CNT].begin()+outsidePoint+1,new_P-1);
    } else {
        hout << "Invalid position of point to create projection on a box boundary. Valid options are: ";
        hout << "inside or outside. Input was: " << currentLocation << endl;
        return 0;
    }
    
    return 1;
}

int RNetwork::Projected_Point(double x, double y, double z, double x1, double y1, double z1, double &xp, double &yp, double &zp )
{
    //Check in which boundary the point is
    if ((x > xmin) && (x1 < xmin) ) {
        //In this case, the projection is at the plane x = xmin
        xp = xmin;
        Projection(xmin, x, x1, y, y1, z, z1, yp, zp);
    } else if ((x < xmax) && (x1 > xmax) ) {
        //In this case, the projection is at the plane x = xmax
        xp = xmax;
        Projection(xmax, x, x1, y, y1, z, z1, yp, zp);
    } else if ((y > ymin) && (y1 < ymin) ) {
        //In this case, the projection is at the plane y = ymin
        yp = ymin;
        Projection(ymin, y, y1, x, x1, z, z1, xp, zp);
    } else if ((y < ymax) && (y1 > ymax) ) {
        //In this case, the projection is at the plane y = ymax
        yp = ymax;
        Projection(ymax, y, y1, x, x1, z, z1, xp, zp);
    } else if ((z > zmin) && (z1 < zmin) ) {
        //In this case, the projection is at the plane z = zmin
        zp = zmin;
        Projection(zmin, z, z1, x, x1, y, y1, xp, yp);
    } else if ((z < zmax) && (z1 > zmax) ) {
        //In this case, the projection is at the plane z = zmax
        zp = zmax;
        Projection(zmax, z, z1, x, x1, y, y1, xp, yp);
    } else {
        hout << "Error when calculating projection. currentPoint and nextPoint are not intersecting a ";
        hout << "box boundary." << endl;
        return 0;
    }
    
    return 1;
}

//This function calculates the projection of a point into a plane between x00 and x01
//x1 and x2 are the projected coordinates. It is assumed x0 = a, as that is the intersecting plane between
//the points (x0_0,x1_0,x2_0) and (x0_1,x1_1,x2_1)
void RNetwork::Projection(double a, double x0_0, double x0_1, double x1_0, double x1_1, double x2_0, double x2_1, double &x1_p, double &x2_p)
{
    
    //The directional components
    double d0, d1, d2;
    d0 = x0_1 - x0_0;
    d1 = x1_1 - x1_0;
    d2 = x2_1 - x2_0;
    
    //Calculate parameter t
    double t = (a - x0_0)/d0;
    
    //Calculate the other two coordinates
    x1_p = x1_0 + t*d1;
    x2_p = x2_0 + t*d2;
}

void RNetwork::Trim_CNT(vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<int> &cnts_inside, vector<double> &radii, long int boundary, int CNT)
{
    //Here check which point is the current point
    vector<long int> empty;
    structure.push_back(empty);
    vector<short int> empty_short;
    boundary_flags_cnt.push_back(empty_short);
    //The CNT will bre trimed from the point after the boundary point and until the last point
    structure.back().insert(structure.back().begin(),structure[CNT].begin()+boundary+1, structure[CNT].end());
    //Erase form CNT
    structure[CNT].erase(structure[CNT].begin()+boundary+1, structure[CNT].end());
    //Add the new CNT at the end of the CNTs inside the box
    cnts_inside.push_back(structure.size()-1);
    cnts_inside_flags.push_back('1');
    
    //Update the CNT number of the points that were moved
    long int P;
    for (long int i = 0; i < structure.back().size(); i++) {
        P = structure.back()[i];
        points_in[P].flag = ((int)structure.size()) - 1;
    }
    
    //Update the radii vector
    //The new CNT is just a segment of the old one so they should have the same radius
    radii.push_back(radii[CNT]);
}

void RNetwork::Add_to_boundary_vectors(Point_3D point3d, long int point)
{
    //Add point and CNT to the boundary vector
    double x = point3d.x;
    double y = point3d.y;
    double z = point3d.z;
    int CNT = point3d.flag;
    if (x == xmin){
        //Add_CNT_to_boundary(bbdyx1_cnt, CNT);
        boundary_flags[point].push_back(0);
        boundary_flags[point].push_back(0);
        boundary_flags_cnt[CNT].push_back(0);
    } else if (x == xmin+lx){
        //Add_CNT_to_boundary(bbdyx2_cnt, CNT);
        boundary_flags[point].push_back(0);
        boundary_flags[point].push_back(1);
        boundary_flags_cnt[CNT].push_back(1);
    } else if (y == ymin){
        //Add_CNT_to_boundary(bbdyy1_cnt, CNT);
        boundary_flags[point].push_back(1);
        boundary_flags[point].push_back(0);
        boundary_flags_cnt[CNT].push_back(2);
    } else if (y == ymin+ly){
        //Add_CNT_to_boundary(bbdyy2_cnt, CNT);
        boundary_flags[point].push_back(1);
        boundary_flags[point].push_back(1);
        boundary_flags_cnt[CNT].push_back(3);
    } else if (z == zmin) {
        //Add_CNT_to_boundary(bbdyz1_cnt, CNT);
        boundary_flags[point].push_back(2);
        boundary_flags[point].push_back(0);
        boundary_flags_cnt[CNT].push_back(4);
    } else if (z == zmin+lz) {
        //Add_CNT_to_boundary(bbdyz2_cnt, CNT);
        boundary_flags[point].push_back(2);
        boundary_flags[point].push_back(1);
        boundary_flags_cnt[CNT].push_back(5);
    }
}

void RNetwork::Add_CNT_to_boundary(vector<int> &boundary, int CNT){
    if (!boundary.size()) {
        boundary.push_back(CNT);
    } else if(boundary.back() != CNT){
        boundary.push_back(CNT);
    }
}

//This function deletes an specified number from a vector
//If the number is not there, nothing happens
void RNetwork::Remove_from_vector(int num, vector<int> &vec)
{
    for (long int i = 0; i < vec.size(); i++)
        if (vec[i] == num) {
            vec.erase(vec.begin()+i);
            break;
        }
}

//This function deletes an specified number from a vector
//If the number is not there, nothing happens
void RNetwork::Remove_from_vector(long int num, vector<long int> &vec)
{
    for (long int i = 0; i < vec.size(); i++)
        if (vec[i] == num) {
            vec.erase(vec.begin()+i);
            break;
        }
}

//This function puts a point on a region where contacts will be looked for.
int RNetwork::Assign_region(const struct RVE_Geo cell_geo, vector<Point_3D> points_in, vector<vector<long int> > structure)
{
    //These variables are to store the size of the RVE and reduces operations when accessing them
    double lx = cell_geo.len_x;
    double ly = cell_geo.wid_y;
    double lz = cell_geo.hei_z;
    //These variables are the coordinates of the lower corner of the RVE that defines its geometry
    double xmin = cell_geo.poi_min.x;
    double ymin = cell_geo.poi_min.y;
    double zmin = cell_geo.poi_min.z;
    
    int sx = secx;
    int sy = secy;
    int sz = secz;
    
    //Sizes of each region
    double dx = lx/(double)sx;
    double dy = ly/(double)sy;
    double dz = lz/(double)sz;
    
    //hout << "sx = " << sx << '\t' << "dx = " << dx << "\n";
    //hout << "sy = " << sy << '\t' << "dy = " << dy << "\n";
    //hout << "sz = " << sz << '\t' << "dz = " << dz  << "\n";
    
    //Check that the regions are not too small for the maximum cutoff distance 2r_max+tunnel
    //If they are, then change the number of sections to the maximum possible
    if (dx < 2*cutoff) {
        dx = 2*cutoff;
        sx = (int)(lx/dx);
        //hout << "Modified the number of sections along x. " << "sx = " << sx << '\t' << "dx = " << dx << endl;
    }
    if (dy < 2*cutoff) {
        dy = 2*cutoff;
        sy = (int)(ly/dy);
        //hout << "Modified the number of sections along y. " << "sy = " << sy << '\t' << "dy = " << dy << endl;
    }
    if (dz < 2*cutoff) {
        dz = 2*cutoff;
        sz = (int)(lz/dz);
        //hout << "Modified the number of sections along z. " << "sz = " << sz << '\t' << "dz = " << dz  << endl;
    }
    
    
    //These variables will give me the region cordinates of the region that a point belongs to
    int a, b, c;
    long int t;
    
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
            if (a == sx) a--;
            b = (int)((y-ymin)/dy);
            //Limit the value of b as it has to go from 0 to sy-1
            if (b == sy) b--;
            c = (int)((z-zmin)/dz);
            //Limit the value of c as it has to go from 0 to sz-1
            if (c == sz) c--;
            
            //Coordinates of non-overlaping region the point belongs to
            double x1 = a*dx +  xmin;
            double y1 = b*dy +  ymin;
            double z1 = c*dz +  zmin;
            double x2, y2, z2;
            if (a == sx-1) x2 = lx +  xmin;
            else x2 = (a+1)*dx +  xmin;
            if (b == sy-1) y2 = ly +  ymin;
            else y2 = (b+1)*dy +  ymin;
            if (c == sz-1) z2 = lz +  zmin;
            else z2 = (c+1)*dz +  zmin;
            
            //Initialize flags for overlaping regions
            int fx = 0;
            int fy = 0;
            int fz = 0;
            
            //Assign value of flag according to position of point
            //The third operand eliminates the periodicity on the boundary
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
            long int temp[2][3] = { {(long int)a, (long int)b, (long int)c}, {(long int)a+fx, (long int)b+fy, (long int)c+fz}};
            
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
long int RNetwork::calculate_t(long int a, long int b, long int c, int sx, int sy)
{
    return a + b*((long int)sx) + c*((long int)sx*sy);
}

//This function makes the lists of contats:
//contacts : vector of vectors that stores point to point contact.
//contacts_cnt : vector of vectors that stores CNT to CNT contact.
//contacts_cnt_point : vector of vectors that stores the points that have a contact of each CNT.
int RNetwork::Check_contacts(vector<Point_3D> points_in, vector<double> radii, double tunnel, long int cnts_t)
{
    //size of the vector for the outer loop
    long int outer = sectioned_domain.size();
    //hout << "outer=" << outer;
    //just need to declare the variable as will change in each outer loop
    //int inner;
    long int inner;
    //These ints are just to store the global point number and the CNTs they belong to.
    //They are just intermediate variables and I only use them to make the code more readable
    long int P1, P2;
    int CNT1, CNT2;
    //Temporary vector of int's to store the contact pair
    vector<long int> empty;
    vector<int> empty_int;
    //the list of contacts has to be the same size as the list of points
    contacts.assign(points_in.size(), empty);
    //This list has to have a size equal to the number of CNTs
    contacts_cnt.assign(cnts_t, empty_int);
    //This list has to have a size equal to the number of CNTs
    contacts_cnt_point.assign(cnts_t, empty);
    //Varable to calculate the cutoff for tunneling
    double cutoff_t;
    //double separation;
    //hout << "tunnel=" << tunnel << endl;
    for (long int i = 0; i < outer; i++) {
        inner = sectioned_domain[i].size();
        for (long int j = 0; j < inner-1; j++) {
            P1 = sectioned_domain[i][j];
            CNT1 = points_in[P1].flag;
            for (long int k = j+1; k<inner; k++) {
                P2 = sectioned_domain[i][k];
                CNT2 = points_in[P2].flag;
                //If distance below the cutoff and points belong to different CNT
                cutoff_t = radii[CNT1] + radii[CNT2] + tunnel;
                //separation = points_in[P1].distance_to(points_in[P2]);
                //hout <<"P1="<<P1<<" CNT1="<<CNT1<<" P2="<<P2<<" CNT2="<<CNT2;
                //hout <<" separation="<<separation<<" cutoff_t="<<cutoff_t;
                //hout <<" CNT1="<<CNT1<<" CNT2="<<CNT2;
                //hout<<" r1="<<radii[CNT1]<<" r2="<<radii[CNT2]<<endl;//
                if ((CNT1!=CNT2)&&(points_in[P1].distance_to(points_in[P2]) <= cutoff_t)) {
                    if (!Check_repeated(contacts[P1], P2)) {
                        //Add point contacts
                        contacts[P2].push_back(P1);
                        contacts[P1].push_back(P2);
                        //Add point contacts in CNT1
                        if (!Check_repeated(contacts_cnt_point[CNT1], P1)) {
                            contacts_cnt_point[CNT1].push_back(P1);
                        }
                        //Add point contacts in CNT2
                        if (!Check_repeated(contacts_cnt_point[CNT2], P2)) {
                            contacts_cnt_point[CNT2].push_back(P2);
                        }
                        //Add CNT contacts
                        if (!Check_repeated(contacts_cnt[CNT1], CNT2)) {
                            contacts_cnt[CNT2].push_back(CNT1);
                            contacts_cnt[CNT1].push_back(CNT2);
                        }
                        
                    }
                    
                }
            }
        }
    }
    //Print2DVec(contacts_cnt, "contacts_cnt.txt");
    //Print2DVec(contacts, "contacts.txt");
    //Print2DVec(contacts_cnt_point, "contacts_cnt_point.txt");
    return 1;
}

//This function checks if the point Point is in the vector region
int RNetwork::Check_repeated(vector<long int> region, long int Point)
{
    for (long int i = 0; i < region.size(); i++) {
        if (Point == region[i]) {
            return 1;
        }
    }
    return 0;
}

//This function checks if the point Point is in the vector region
int RNetwork::Check_repeated(vector<int> region, int Point)
{
    for (long int i = 0; i < region.size(); i++) {
        if (Point == region[i]) {
            return 1;
        }
    }
    return 0;
}

//This function deletes contact from the contact list given by contacts_vector
//Written in this way, this function can be used to delete a point or a cnt contact.
int RNetwork::Delete_contact(int contact, vector<vector<int> > &contacts_vector)
{
    int other_contact, flag = 0;
    //First delete all contacts that others have with "contact"
    for (long int i = 0; i < contacts_vector[contact].size(); i++) {
        //Finde the other contact of contact
        other_contact = contacts_vector[contact][i];
        //Go to that other contact, find contact and delete it
        for (long int j = 0; j < contacts_vector[other_contact].size(); j++) {
            if (contacts_vector[other_contact][j] == contact) {
                contacts_vector[other_contact].erase(contacts_vector[other_contact].begin()+j);
                flag = 1;
                break;
            }
        }
        if (!flag) {
            hout << "ERROR. Could not find contact of point " << contact << " with point " << other_contact;
            hout << ". No contact was deleted. " << endl;
            return 0;
        }
    }
    
    //Second, delete all contacts that "contact" has. That is empty contacts_vector[contact]
    contacts_vector[contact].clear();
    
    return 1;
}

//This function deletes contact from the contact list given by contacts_vector
//Written in this way, this function can be used to delete a point or a cnt contact.
int RNetwork::Delete_contact(long int contact, vector<vector<long int> > &contacts_vector)
{
    long int other_contact, flag = 0;
    //First delete all contacts that others have with "contact"
    for (long int i = 0; i < contacts_vector[contact].size(); i++) {
        //Finde the other contact of contact
        other_contact = contacts_vector[contact][i];
        //Go to that other contact, find contact and delete it
        for (long int j = 0; j < contacts_vector[other_contact].size(); j++) {
            if (contacts_vector[other_contact][j] == contact) {
                contacts_vector[other_contact].erase(contacts_vector[other_contact].begin()+j);
                flag = 1;
                break;
            }
        }
        if (!flag) {
            hout << "ERROR. Could not find contact of point " << contact << " with point " << other_contact;
            hout << ". No contact was deleted. " << endl;
            return 0;
        }
    }
    
    //Second, delete all contacts that "contact" has. That is empty contacts_vector[contact]
    contacts_vector[contact].clear();
    
    return 1;
}

void RNetwork::Clear_vectors()
{
    //Cluster vectors
    clusters_cnt.clear();
    directional_clusters.clear();
    isolated.clear();
    dead_branches.clear();
    //Family number
    family.clear();
    //Boundary vectors
    bbdyx1_cnt.clear();
    bbdyx2_cnt.clear();
    bbdyy1_cnt.clear();
    bbdyy2_cnt.clear();
    bbdyz1_cnt.clear();
    bbdyz2_cnt.clear();
    //Percolation flags
    boundary_flags_cnt.clear();
    percolation_flags.clear();
    //contacts vectors
    contacts.clear();
    contacts_cnt.clear();
    contacts_cnt_point.clear();
    //Length vectors
    clusters_lengths.clear();
    dead_branches_lengths.clear();
    clusters_fractions.clear();
    //Sphere vectors
    sphere_c.clear();
    sphere_r.clear();
    closest_distance.clear();
    //Boundary flags
    boundary_flags.clear();
    //Inside flags
    cnts_inside.clear();
    cnts_inside_flags.clear();
    //Sub-domains to check contacts
    sectioned_domain.clear();
}

//This function will generate the clusters of CNT numbers
int RNetwork::Make_CNT_clusters(struct CNT_Geo cnts_geo, vector<Point_3D> points_in, vector<vector<long int> > structure)
{
    int CNT;
    vector <int> empty;
    //Scan the CNTs inside and find single isloated CNTs
    
    //If the CNT has no contacts, there is still the chance that it conducts.
    //This conduction without contacts can only happen if the CNT spans from one boundary to the other and
    //this happens only if the (maximum) length of the CNT is greater than the size of the box
    //So I check if this is the case, otherwise the CNT is isolated
    if ( (cnts_geo.len_max > lx) || (cnts_geo.len_max > ly) || (cnts_geo.len_max > lz) ){
        //hout << "Check1 ";
        for (int i = 0; i < (int)cnts_inside.size(); i++) {
            CNT = cnts_inside[i];
            //hout << "Check2 "<< "CNT=" << CNT <<" contacts_cnt.size="<<contacts_cnt.size()<<' ';
            if (!contacts_cnt[CNT].size()) {
                //If the CNT has no contacts, check if it spans the whole sample, i.e. the CNT is
                //in itself one cluster
                //hout << "Check3 " <<" structure[CNT].size="<< structure[CNT].size();
                long int front = structure[CNT].front();
                long int back = structure[CNT].back();
                //hout << "Check3.1 " << endl;
                if ( (Where_is(points_in[front].x, points_in[front].y, points_in[front].z)=="boundary") && (Where_is(points_in[back].x, points_in[back].y, points_in[back].z)=="boundary")) {
                    //Add the CNT as a new cluster
                    vector<int> empty;
                    clusters_cnt.push_back(empty);
                    clusters_cnt.back().push_back(CNT);
                    //Update the boundary_flags_cnt vector
                    //Initialize the vector of flags the cluster
                    vector<short int> cluster_flag(7,0);
                    //Add corresponding flags
                    for (int a = 0; a < (int)boundary_flags_cnt[CNT].size(); a++) {
                        short int flag = boundary_flags_cnt[CNT][a];
                        cluster_flag[flag] = 1;
                    }
                    //Add to vector of cluster flags
                    percolation_flags.push_back(cluster_flag);
                } else {
                    //hout << "Check3b " << endl;
                    isolated.push_back(empty);
                    isolated.back().push_back(CNT);
                }
            }
        }
    } else {
        //hout << "Check4 ";
        for (int i = 0; i < (int)cnts_inside.size(); i++) {
            //hout << "Check5 ";
            CNT = cnts_inside[i];
            if (!contacts_cnt[CNT].size()) {
                //hout << "Checkt6 ";
                isolated.push_back(empty);
                isolated.back().push_back(CNT);
            }
        }
    }
    //hout << "Check7 " << "clusters_cnt.size() " << clusters_cnt.size() << ' ';
    
    //Then make clusters of CNTs
    int count;
    vector<int> vec;
    
    //I need to modify contacts_cnt to make the clusters, so I create a temporary vector and do all changes there
    vector<vector<int> > contacts_cnt_tmp = contacts_cnt;
    
    for (int i = 0; i < (int)cnts_inside.size(); i++) {
        //hout << "cnts_inside[i]=" << cnts_inside[i] << ", i=" << i << endl;
        CNT = cnts_inside[i];
        //hout << "Check7.1 ";
        if (contacts_cnt_tmp[CNT].size()) {
            vec = contacts_cnt_tmp[CNT];
            vector<short int> cluster_flag(7,0);
            //hout << "Check7.2 ";
            count = -1;
            //This will help save time when creating the clusters. Instead of scanning all CNTs from the beginning at every loop
            //in the function Single_cluster, it wil scan only the newyle added CNTs
            int start = 0;
            do {
                //hout << "Check7.3 ";
                //count = Single_cluster(CNT, start, vec, contacts_cnt_tmp);
                count = Single_cluster(CNT, start, vec, contacts_cnt_tmp,cluster_flag);
                
                Discard_repeated(vec);
                //At the following loop, the search of CNT contacts has to start from the first newly added CNT
                //This CNT will start at index contacts_cnt_tmp[CNT].size() since contacts_cnt_tmp[CNT] is the initial cluster
                //and vec is contacts_cnt_tmp[CNT] with more CNTs
                start = (int)contacts_cnt_tmp[CNT].size();
                //hout << "Check7.5 ";
                //Update the vector that contains the new cluster
                contacts_cnt_tmp[CNT] = vec;
                /*/
                 for (unsigned j = 0; j < contacts_cnt_tmp[CNT].size(); j++) {
                 hout << contacts_cnt_tmp[CNT][j] << ' ';
                 }
                 hout << endl;//*/
            } while (count != 0);
            //hout << "Check7.6 ";
            //Now in contacts_cnt_tmp[CNT] we have a cluster. So save it in the vector of CNT clusters
            clusters_cnt.push_back(contacts_cnt_tmp[CNT]);
            percolation_flags.push_back(cluster_flag);
            //Clear the vector variable to find the next cluster
            vec.clear();
            //hout << "Check7.7 ";
        }
    }
    return 1;
}

//This function creates one cluster at a time
int RNetwork::Single_cluster(int cnt_seed, int start, vector<int> &vec, vector<vector<int> > &contacts_vector)
{
    int CNT;
    int count = 0;
    int limit = (int)vec.size();
    //hout << "\nCNT base = " << index << endl;
    for (int i = start; i < limit; i++) {
        CNT = vec[i];
        //hout << "CNT=" << CNT << " v2[CNT].size()=" << v2[CNT].size() << ' ';
        if ( (contacts_vector[CNT].size()) && (CNT != cnt_seed) ) {
            for (int k = 0; k < (int) contacts_vector[CNT].size(); k++) {
                //hout << ' ' << v2[CNT][k];
                vec.push_back(contacts_vector[CNT][k]);
                count++;
            }
            //hout << endl;
            contacts_vector[CNT].clear();
        }
    }
    return count;
}

//This function creates one cluster at a time
int RNetwork::Single_cluster(int cnt_seed, int start, vector<int> &vec, vector<vector<int> > &contacts_vector, vector<short int> &cluster_flag)
{
    int CNT;
    int count = 0;
    int limit = (int)vec.size();
    //hout << "\nCNT base = " << index << endl;
    for (int i = start; i < limit; i++) {
        CNT = vec[i];
        //hout << "CNT=" << CNT << " v2[CNT].size()=" << v2[CNT].size() << ' ';
        if ( (contacts_vector[CNT].size()) && (CNT != cnt_seed) ) {
            for (int k = 0; k < (int) contacts_vector[CNT].size(); k++) {
                //hout << ' ' << v2[CNT][k];
                vec.push_back(contacts_vector[CNT][k]);
                for (int a = 0; a < (int)boundary_flags_cnt[CNT].size(); a++) {
                    short int flag = boundary_flags_cnt[CNT][a];
                    cluster_flag[flag] = 1;
                }
                count++;
            }
            //hout << endl;
            contacts_vector[CNT].clear();
        }
    }
    return count;
}

//This function checks percolation in any direction in all clusters and removes the non-percolating clusters
int RNetwork::Check_clusters_percolation(vector<Point_3D> points_in, vector<vector<long int> > structure)
{
    int count = 0;
    //Check if there is any cluster at all. If the size of clusters_cnt is non-zero, then there are clusters
    if (clusters_cnt.size()) {
        //Move all non-percolating clusters to the corresponding vector<vector>
        //I start from the end of the vector to avoid issues with the index when removing an element
        int size = (int)clusters_cnt.size();
        //This variable will store the family number given by CheckPercolationSingleCluster
        int fam;
        for (int i = size - 1; i >= 0 ; i--) {
            //hout <<"Check_clusters_percolation " << i << " size " << clusters_cnt.size() << endl;
            //hout << "percolation_flags.size()=" << percolation_flags.size() << endl;
            //if (!Check_percolation_single_cluster(clusters_cnt[i], fam)){
            if (!Check_percolation_single_cluster(percolation_flags[i], fam)){
                //percolation_flags and clusters_cnt have the same size
                //hout << "NO\n";
                isolated.push_back(clusters_cnt[i]);
                //remove the non-percolating cluster
                clusters_cnt.erase(clusters_cnt.begin()+i);
                count++;
            } else {
                //Add the family number to the vector of families
                family.insert(family.begin(), fam);
            }
        }
        //Find the sphere that encloses each isolated cluster
        //Create the spheres that enclose the isolated clusters
        if (!Find_spheres(points_in, structure)) {
            hout << "Error in Find_spheres." << endl;
            return 0;
        }
        
        //hout << "fam size = " << family.size()<< "\t cluster size = "<<  clusters_cnt.size() << endl;
        //Just a check. family and clusters_cnt MUST have the same size
        if (family.size() != clusters_cnt.size()) {
            hout <<  "ERROR the vector of family number and clusters_cnt do not have the same size" << endl;
            hout << "fam size = " << family.size()<< "\t cluster size = "<<  clusters_cnt.size() << endl;
            return 0;
        }
        //If the count variable is equal to the initial number of clusters, that means that none of the clusters percolated
        if (count == size) {
            hout <<  "None of the clusters percolated" << endl;
            //Print2DVec(clusters_cnt_iso, "isolated.txt");
            return 0;
        }
    } else {
        //Find the sphere that encloses each isolated cluster
        //Create the spheres that enclose the isolated clusters
        if (!Find_spheres(points_in, structure)) {
            hout << "Error in Find_spheres." << endl;
            return 0;
        }
        
        hout << "The size of clusters_cnt is found to be zero i.e there are no clusters ";
        hout << "but only isolated CNTs." << endl;
        //Print2DVec(clusters_cnt_iso, "isolated.txt");
        return 0;
    }
    //Print2DVec(clusters_cnt_iso, "isolated.txt");
    return 1;
}

//This function will check if there is percolation for a single cluster in x, y and/or z directions
int RNetwork::Check_percolation_single_cluster(vector<int> cluster, int &family)
{
    ///*
    //Flags that will let me know if the cluster is in contact with a boundary
    int Boundary[] = {0, 0, 0, 0, 0, 0};
    //Falgs that will tell me the family the cluster belogns to
    int percolating[] = {0, 0, 0, 0, 0, 0, 0};
    //Check contact with x-boundary
    Boundary[0] = Intersection(cluster, bbdyx1_cnt);
    Boundary[1] = Intersection(cluster, bbdyx2_cnt);
    //Check contact with y-boundary
    Boundary[2] = Intersection(cluster, bbdyy1_cnt);
    Boundary[3] = Intersection(cluster, bbdyy2_cnt);
    //Check contact with z-boundary
    Boundary[4] = Intersection(cluster, bbdyz1_cnt);
    Boundary[5] = Intersection(cluster, bbdyz2_cnt);
    
    percolating[0] = Boundary[0] && Boundary[1]; //x-x only
    percolating[1] = Boundary[2] && Boundary[3]; //y-y only
    percolating[2] = Boundary[4] && Boundary[5]; //z-z only
    percolating[3] = percolating[0] && percolating[1]; //x-x and y-y only
    percolating[4] = percolating[0] && percolating[2]; //x-x and z-z only
    percolating[5] = percolating[1] && percolating[2]; //y-y and z-z only
    percolating[6] = percolating[0] && percolating[1] && percolating[2]; //x-x, y-y and z-z
    
    //Scan the percolating array backwards
    for (int i = 6; i >=0 ; i--)
        if (percolating[i]){
            //If there is a non-zero percolating[i], then this cluster percolates and belongs to family i
            family = i;
            return 1;
        }
    //If all percolating[i] were 0, then there is no percolation in the cluster
    return 0;
    
}

//This function will check if there is percolation for a single cluster in x, y and/or z directions
int RNetwork::Check_percolation_single_cluster(vector<short int> cluster_flag, int &family)
{
    ///*
    //Falgs that will tell me the family the cluster belogns to
    int percolating[] = {0, 0, 0, 0, 0, 0, 0};
    
    //Boolean operations to find the different percolation directions
    percolating[0] = cluster_flag[0] && cluster_flag[1]; //x-x only
    percolating[1] = cluster_flag[2] && cluster_flag[3]; //y-y only
    percolating[2] = cluster_flag[4] && cluster_flag[5]; //z-z only
    percolating[3] = percolating[0] && percolating[1]; //x-x and y-y only
    percolating[4] = percolating[0] && percolating[2]; //x-x and z-z only
    percolating[5] = percolating[1] && percolating[2]; //y-y and z-z only
    percolating[6] = percolating[0] && percolating[1] && percolating[2]; //x-x, y-y and z-z
    
    //Scan the percolating array backwards
    for (int i = 6; i >=0 ; i--)
        if (percolating[i]){
            //If there is a non-zero percolating[i], then this cluster percolates and belongs to family i
            family = i;
            return 1;
        }
    //If all percolating[i] were 0, then there is no percolation in the cluster
    return 0;
    
}

//This function finds if a number in vec1 is in vec2
//Returns 1 when a number in vec1 is in vec 2, i.e. there is non empty intersection
//Returns 0 when no number in vec1 is in vec 2, i.e. the intersection is empty
int RNetwork::Intersection(vector<int> vec1, vector<int> vec2)
{
    for (int i = 0; i < (int)vec1.size(); i++)
        for (int j = 0; j < (int)vec2.size(); j++)
            if (vec1[i] == vec2[j])
                return 1;
    return 0;
}

int RNetwork::Find_spheres(vector<Point_3D> points_in, vector<vector<long int> > structure) {
    //Variables for the extreme coordinates.
    double minx, maxx, miny, maxy, minz, maxz;
    //Sphere geometry
    double sx, sy, sz, sr;
    //Point coordinates
    double x, y, z, r;
    //
    int CNT;
    long int P;
    //Scan an isolated cluster and find the extreme corrdinates
    for (int i = 0; i < (int)isolated.size(); i++) {
        //hout << "Check1 ";
        //hout << "Check2 cluster=" << i << endl;
        minx = xmin + lx;
        maxx = xmin;
        miny = ymin + ly;
        maxy = ymin;
        minz = zmin + lz;
        maxz = zmin; //*/
        for (int j = 0; j < (int)isolated[i].size(); j++) {
            //hout << "Check3 ";
            //Set the current CNT
            CNT = isolated[i][j];
            //Look for extreme points
            //hout << "Check4 CNT="<< CNT << " size="<< structure[CNT].size() << ' ' << endl;
            for (int k = 0; k < (int)structure[CNT].size(); k++) {
                //hout << "Check5 ";
                //
                P = structure[CNT][k];
                x = points_in[P].x;
                y = points_in[P].y;
                z = points_in[P].z;
                //hout << "Check6 " << x << '\t' << y << '\t' << z << endl;
                if (x < minx)
                    minx = x;
                if (x > maxx)
                    maxx = x;
                if (y < miny)
                    miny = y;
                if (y > maxy)
                    maxy = y;
                if (z <  minz)
                    minz = z;
                if (z > maxz)
                    maxz = z;
            }
        }
        //After the whole cluster is scanned I have the extrme variables, so I calculate the center of the sphere
        sx = (maxx + minx)/2;
        sy = (maxy + miny)/2;
        sz = (maxz + minz)/2;
        //hout << "Extremes " << minx << ' ' << miny << ' ' << minz << ' ' << maxx << ' ' << maxy << ' ' << maxz << endl;
        //Add the point that represents the sphere center
        sphere_c.push_back(Point_3D(sx,sy,sz));
        //hout << "Check7 ";
    }
    
    //Scan again the cluster and find the farthest point. This will be the radius
    for (int i = 0; i < (int)isolated.size(); i++) {
        //hout << "Check8 ";
        sr = 0;
        for (int j = 0; j < (int)isolated[i].size(); j++) {
            //hout << "Check9 ";
            //Set the current CNT
            CNT = isolated[i][j];
            for (int k = 0; k < (int)structure[CNT].size(); k++) {
                //hout << "Check10 ";
                P = structure[CNT][k];
                r = sphere_c[i].distance_to(points_in[P]);
                //hout << "Check11 ";
                if (r > sr) {
                    sr = r;
                }
            }
        }
        //hout << "Check12 ";
        //After a whole cluster is scanned I have the radius of the sphere so I just add it to the vector of radii
        sphere_r.push_back(sr);
        //hout << "Check13 ";
    }
    
    if (sphere_c.size() != sphere_r.size()) {
        hout << "ERROR. The vector of sphere centers has a different size than the vector of sphere radii: ";
        hout << "sphere_c.size() = " << sphere_c.size() << " sphere_c.size() = "<< sphere_r.size() << endl;
        return 0;
    }
    
    //--------------------------------------------------------------------------------------------------------
    //Once I have all centers and radii I scan the sphere_c vector to find the closest distance to each sphere
    
    //Resize the vectors that store the closest sphere to the necessary size
    //Initialize each element of the closest_distance vector with a large value
    closest_distance.assign(sphere_c.size(), lx + ly + lz);
    for (int i = 0; i < (int)sphere_c.size()-1; i++) { //This for-loop has to end one element before the last
        //hout << "Check14 ";
        for (int j = i+1; j < (int)sphere_c.size(); j++) {
            //hout << "Check15 ";
            r = sphere_c[i].distance_to(sphere_c[j]);
            if (r < closest_distance[i]) {
                //hout << "Check16 ";
                closest_distance[i] = r;
                closest_distance[j] = r;
                sphere_c[i].flag = j;
                sphere_c[j].flag = i;
                //hout << "Check17 ";
            }
        }
    }
    //hout << "Check18 " << endl;
    
    /*/struct elliparam
     {
     double x, y, z;		//Õ÷«Úµƒ÷––ƒµ„
     double a, b, c;		//Õ÷«Úµƒ≥§°¢÷–∫Õ∂Ã÷·
     double alpha1, alpha2, alpha3; //Õ÷«Úµƒæ≈∏ˆº–Ω«≤Œ ˝
     double beta1, beta2, beta3;
     double gamma1, gamma2, gamma3;
     };//*/
    
    //--------------------------------------------------------------------------------------------------------
    //Save the data to a file
    ofstream sphere("spheres.txt");
    /*/Visualization files
     elliparam ellipse;
     ellipse.alpha1 = 0;
     ellipse.alpha2 = 0;
     ellipse.alpha3 = 0;
     ellipse.beta1 = 0;
     ellipse.beta2 = 0;
     ellipse.beta3 = 0;
     ellipse.gamma1 = 0;
     ellipse.gamma2 = 0;
     ellipse.gamma3 = 0;
     vector<struct elliparam> ellipses;
     for (int i=0; i < (int)sphere_c.size(); i++) {
     sphere << sphere_c[i].x << "\t" << sphere_c[i].y << "\t" << sphere_c[i].z << "\t" << sphere_r[i] << "\t" ;
     sphere << closest_distance[i] << "\t" << sphere_c[i].flag << "\n";
     ellipse.x = sphere_c[i].x;
     ellipse.y = sphere_c[i].y;
     ellipse.z = sphere_c[i].z;
     ellipse.a = sphere_r[i];
     ellipse.b = sphere_r[i];
     ellipse.c = sphere_r[i];
     ellipses.push_back(ellipse);
     }
     sphere.close();
     Export_cluster_ellipsoids_mesh(cell_geometry, ellipses, 0);//*/
    
    return 1;
}


//This function applies the electirfying algorithm to each cluster of CNTs
int RNetwork::Split_cnts(vector<Point_3D> &points_in, vector<double> &cnts_radius, vector<vector<long int> > &structure)
{
    //Initialize matrices
    vector<int> LM_matrix, dead, empty;
    //Set the sizes for the clusters vectors
    directional_clusters.resize(7, empty);
    dead_branches.resize(7, empty);
    
    //Matrix to store the solution
    MathMatrix voltages;
    
    //Variables for SSS algorithm
    vector<long int> col_ind, row_ptr;
    vector<double> values, diagonal;
    MathMatrix R;
    
    //Variable for the number of nodes
    int global_nodes;
    
    hout << "Percolating clusters = "<< clusters_cnt.size() << endl;
    //Print2DVec(clusters_cnt, "clusters_cnt.txt");
    
    //The split CNTs has to be done for each CNT cluster
    for (int i = 0; i < (int)clusters_cnt.size(); i++) {
        //hout << "Check1 ";
        //Print1DVec(clusters_cnt[i], "cluster.txt");
        //If there is only one CNT in the cluster, that means that it percolates in a given direction
        //and is in itself a percolating cluster. Only  when there is more than one CNT, i.e. when the
        //size of the cluster is greater than one, the DEA is applied. Otherwise The CNT number is added
        //to the corresponding directional_cluster
        if (clusters_cnt[i].size() > 1) {
            LM_matrix.assign(points_in.size(), -1);
            global_nodes = Get_LM_matrix(family[i], clusters_cnt[i], structure, contacts_cnt_point, LM_matrix);
            //hout << "Check2 global_nodes=" << global_nodes << ' ';
            //Check that the number of nodes is not zero
            if (!global_nodes)
                return 0;
            
            //hout << "Check3 ";
            //Fill the sparse stiffness matrix
            Fill_sparse_stiffness_matrix(global_nodes, clusters_cnt[i], structure, LM_matrix, col_ind, row_ptr, values, diagonal, R);
            //Extract the electric backbone using the direct electrifying algorithm. To solve the system of
            //equations, the conjugate gradient (CG) is used. To reduce the computational cost of vector-matrix
            //multiplication, and given that the stiffness matrix is symmetric and sparse, the Symmetric Sparse
            //Skyline (SSS) format is used.
            //hout << "Check4 " << endl;
            voltages = Solve_DEA_equations_CG_SSS(global_nodes, col_ind, row_ptr, values, diagonal, R);
            //hout << "Check5 " << endl;
            
            //Remove CNTs that do not carry any current
            if (!Remove_currentless_CNTs(voltages, clusters_cnt[i], LM_matrix, dead, cnts_radius, structure)){
                hout << "Error when removing the non percolating branches." << endl;
                return 0;
            }
            //hout << "Check6 " << endl;
            
            //Add the current cluster to the family it corresponds to in directional_clusters
            directional_clusters[family[i]].insert(directional_clusters[family[i]].begin(), clusters_cnt[i].begin(),clusters_cnt[i].end());
            //hout << "Check7 " << endl;
            //If there are dead branches then add the dead_branches vector in the index corresponding to its family
            dead_branches[family[i]].insert(dead_branches[family[i]].begin(), dead.begin(),dead.end());
            //Clear the dead vector for the next iteration
            //hout << "Check8 " << endl;
            dead.clear();
        } else {
            //Add the only CNT of the current cluster to the family it corresponds to in directional_clusters
            directional_clusters[family[i]].push_back(clusters_cnt[i][0]);
        }
    }
    
    //After all clusters are
    return 1;
}

//Build the LM matrix
int RNetwork::Get_LM_matrix(int family, vector<int> cluster, vector<vector<long int> > structure, vector<vector<long int> > &contacts_cnt_point, vector<int> &LM_matrix)
{
    //Variables
    int CNT;
    long int P;
    
    //Node numbers start in 2, as 0 and 1 are reserved for the boundaries where the voltage is applied
    int node = 2;
    //Assign a node number to each contact point and CNT endpoints
    vector<int> empty;
    vector<vector<int> > LM, points;
    for (int i = 0; i < (int)cluster.size(); i++) {
        CNT = cluster[i];
        
        //Sort the contatcs vector
        if (!Sort_vector(contacts_cnt_point[CNT], structure[CNT])) {
            hout << "ERROR in sorting the contacts vector contacts_cnt_point["<<CNT<<']'<<endl;
            return 0;
        }
        
        //Scan the CNT contacts
        LM.push_back(empty);
        points.push_back(empty);
        for (int j = 0; j < (int)contacts_cnt_point[CNT].size(); j++) {
            //Point number
            P = contacts_cnt_point[CNT][j];
            //check if the point is in a relevant boudary
            if ((boundary_flags[P].size()==2) && Is_in_relevant_boundary(family, boundary_flags[P][0])) {
                //If the point is in a relevant boundary add the reserved node number
                LM_matrix[P] = boundary_flags[P][1];
                LM.back().push_back(boundary_flags[P][1]);
                points.back().push_back(P);
            } else {
                //If the point is not in a boundary, then add a new node number to the CNT
                LM_matrix[P] = node;
                LM.back().push_back(node);
                points.back().push_back(P);
                //Increase the number of nodes
                node++;
            }
        }
    }
    //Print2DVec(LM, "LM.txt");
    //Print2DVec(points, "points.txt");
    //Print2DVec(contacts_cnt_point, "contacts_cnt_point_sorted.txt");
    
    return node;
}

//This function sorts the contacts_cnt_point vector, including the case when some of the contact points are
//new points (added due to the CNT crossing the inner box of the sample).
//The back and front of the CNT are also added
int RNetwork::Sort_vector(vector<long int> &contacts, vector<long int> cnt){
    
    //Sort the vector using the built-in function
    sort(contacts.begin(), contacts.end());
    
    //I will use the size of the vector to locate the previous to last element
    int size = (int)contacts.size();
    
    //Check the last two points for added boundary points
    if ( (size >=2) && (contacts[size-2] >= points0) ){ //Size has to be at least 2
        //hout << "Sort10 ";
        //If the previous to last point is a new point, so is the last point.
        //Then I remove these two points and add the front of the CNT to the front of the contacts
        //and the back of the CNT to the back of the contacts
        contacts.pop_back();
        contacts.pop_back();
        contacts.push_back(cnt.back());
        contacts.insert(contacts.begin(), cnt.front());
    } else if ( (size >= 1) && (contacts.back() >= points0) ){ //Size has to be at least 1
        //hout << "Sort20 ";
        //If only the last point is a new point, I remove it and add the endpoints to the contacts vector
        contacts.pop_back();
        //Because only one is for sure an endpoint, I need to check if the endpoint is already there
        if (contacts.back() != cnt.back())
            contacts.push_back(cnt.back());
        if (contacts.front() != cnt.front())
            contacts.insert(contacts.begin(), cnt.front());
    } else {
        //hout << "Sort30 ";
        if (size) {
            //If there are no new points in the contact vector, then check if the back and front of the CNT
            //need to be added to the contacts vector
            if (contacts.back() != cnt.back())
                contacts.push_back(cnt.back());
            if (contacts.front() != cnt.front())
                contacts.insert(contacts.begin(), cnt.front());
        } else {
            //If the contacts vector is empty, just add both endpoints. This might mean that the CNT
            //spans the whole sample
            contacts.push_back(cnt.front());
            contacts.push_back(cnt.back());
        }
    }
    return 1;
}

//This function checks if a boundary point is in a relevant boundary depending on the family the cluster belongs to.
int RNetwork::Is_in_relevant_boundary(int family, int boundary_node)
{
    if ( (boundary_node==0) && ((family==0)||(family==3)||(family==4)||(family==6)) )
        return 1;
    if ( (boundary_node==1) && ((family==1)||(family==3)||(family==5)||(family==6)) )
        return 1;
    if ( (boundary_node==2) && ((family==2)||(family>=4)) )
        return 1;
    else
        return 0;
}

void RNetwork::Fill_sparse_stiffness_matrix(long int nodes, vector<int> cluster, vector<vector<long int> > structure, vector<int> LM_matrix, vector<long int> &col_ind, vector<long int> &row_ptr, vector<double> &values, vector<double> &diagonal, MathMatrix &R)
{//col_ind, row_ptr, values, diagonal
    //Variables
    int CNT;
    long int P1, P2, node1, node2;
    vector<vector<long int> > col_ind_2d;
    vector<long int> empty_long;
    col_ind_2d.assign(nodes, empty_long);
    vector<vector<double> > values_2d;
    vector<double> empty_double;
    values_2d.assign(nodes, empty_double);
    diagonal.clear();
    diagonal.assign(nodes, 0);
    
    //Scan every CNT in the cluster
    for (long int i = 0; i < cluster.size(); i++) {
        CNT = cluster[i];
        for (long int j = 0; j < contacts_cnt_point[CNT].size()-1; j++) {
            //Find node numbers of the first two elements
            P1 = contacts_cnt_point[CNT][j];
            node1 = LM_matrix[P1];
            P2 = contacts_cnt_point[CNT][j+1];
            node2 = LM_matrix[P2];
            //hout << " P1="<<P1<<" node1="<<node1<<" P2="<<P2<<" node2="<<node2<<endl;
            
            //Fill the diagonal elements of the stiffness matrix
            diagonal[node1] += 1;
            diagonal[node2] += 1;
            //Fill the off diagonal elements of the stiffness matrix
            Add_to_sparse_stiffness(node1, node2, col_ind_2d, values_2d);
            
            //Check if the current node1 has any contacts and add the corresponding contributions to the
            //stiffness matrix
            if (contacts[P1].size()) {
                for (int k = 0; k < (int)contacts[P1].size(); k++) {
                    P2 = contacts[P1][k];
                    node2 = LM_matrix[P2];
                    //hout << " P1="<<P1<<" node1="<<node1<<" contact P2="<<P2<<" node2="<<node2<<endl;
                    //Fill the diagonal elements of the stiffness matrix
                    diagonal[node1] += 1;
                    diagonal[node2] += 1;
                    //Fill the stiffness matrix
                    Add_to_sparse_stiffness(node1, node2, col_ind_2d, values_2d);
                    //hout << "Added ";
                    Remove_from_vector(P1, contacts[P2]);
                    //hout << "Removed ";
                }
                contacts[P1].clear();
            }
        }
        //Check if the last node has any contacts and add the corresponding contributions to the
        //stiffness matrix
        P1 = contacts_cnt_point[CNT].back();
        node1 = LM_matrix[P1];
        if (contacts[P1].size()) {
            for (int k = 0; k < (int)contacts[P1].size(); k++) {
                P2 = contacts[P1][k];
                node2 = LM_matrix[P2];
                //hout << " P1="<<P1<<" node1="<<node1<<" last contact P2="<<P2<<" node2="<<node2<<endl;
                //Fill the diagonal elements of the stiffness matrix
                diagonal[node1] += 1;
                diagonal[node2] += 1;
                //Fill the stiffness matrix
                Add_to_sparse_stiffness(node1, node2, col_ind_2d, values_2d);
                //hout << "Added ";
                Remove_from_vector(P1, contacts[P2]);
                //hout << "Removed ";
            }
            contacts[P1].clear();
        }
    }
    
    //------------------------------------------------------------------------
    //Convert from 2D vectors to 1D vectors
    
    //Clear vectors
    values.clear();
    row_ptr.clear();
    col_ind.clear();
    R.element.clear();
    
    //The two first elements of row_ptr are zero
    //row_ptr.push_back(0);
    row_ptr.push_back(0);
    empty_double.push_back(0);
    R.element.assign(nodes-2, empty_double);
    
    //Fill the rest of the elements
    //Nodes 0 and 1 are to be ignored
    for (long int i = 2; i < col_ind_2d.size(); i++) {
        for (long int j = 0; j < col_ind_2d[i].size(); j++) {
            if (col_ind_2d[i][j] > 1) {
                values.push_back(values_2d[i][j]);
                //The column numbers that I need are the numbers in the full matrix-2
                col_ind.push_back(col_ind_2d[i][j]-2);
            } else if (!col_ind_2d[i][j]) {
                //When the column index is zero, I need to save the value of the resistance on the vector R
                //that will be used in the SSS algorithm for the CG algorithm
                R.element[i-2][0] = values_2d[i][j];
            }
        }
        row_ptr.push_back(values.size());
    }
    //Rmove the first two elements of the diagonal as they are not used
    diagonal.erase(diagonal.begin());
    diagonal.erase(diagonal.begin());
    
    //Print2DVec(col_ind_2d, "col_ind_2d.txt");
}

//This functions assumes node1 > node2
//Only store the value of the resistance between two nodes on the node with the largest node number
//In this way I will store only the lower triangular part of the stiffness matrix
void RNetwork::Add_to_sparse_stiffness(long int node1, long int node2, vector<vector<long int> > &col_ind, vector<vector<double> > &values)
{
    if (node1 > node2) {
        col_ind[node1].push_back(node2);
        //This is the resistance between the two nodes
        values[node1].push_back(-1);
    } else {
        col_ind[node2].push_back(node1);
        //This is the resistance between the two nodes
        values[node2].push_back(-1);
    }
}

//This function solves the equation of the electric circuit as done in the Direct Electrifing Algorithm (DEA)
MathMatrix RNetwork::Solve_DEA_equations_CG_SSS(long int nodes, vector<long int> col_ind, vector<long int> row_ptr, vector<double> values, vector<double> diagonal, MathMatrix R)
{
    //Vectors for the CG algorithm
    MathMatrix X(nodes-2,1), P(nodes-2,1);
    
    //Voltage applied to the sample
    double voltage = (double)(nodes);
    //double voltage = 1;
    
    //Use the adecuate voltage for the vector R
    for (long int i = 0; i < R.element.size(); i++) {
        //Only if non-zero perform the multiplication
        if (R.element[i][0]) {
            R.element[i][0] = -voltage*R.element[i][0];
            P.element[i][0] = R.element[i][0];
        }
    }
    
    //Known vector
    MathMatrix d1(2,1);
    
    //Boundary conditions
    d1.element[0][0] = voltage;
    d1.element[1][0] = 0;
    
    //Variables of the algorithm
    MathMatrix AP(nodes-2,1);
    double alpha, beta, rr0, rr;
    
    //=========================================
    // Conjugate Gradient Algorithm
    
    //Iteration
    long int k;
    long int max_iter = 10*nodes;
    int test = 500;
    
    //Initial residual
    double R0 = 1.0E-10*sqrt(fabs(V_dot_v(R, R)));
    //double R0 = Zero*sqrt(fabs(V_dot_v(R, R)));
    
    for (k = 1; k <= max_iter; k++) {
        //Calculate Ap
        AP = spM_V_SSS(P, row_ptr, col_ind, diagonal, values);
        //Calculate norm or residual of step k-1. Will be used later as convergence criteria and to calculate beta
        rr0 = V_dot_v(R, R);
        //Step length
        alpha = rr0/(V_dot_v(P, AP));
        //Approximate solution
        X = X + P*alpha;
        //Residual
        R = R - AP*alpha;
        //Calculate norm or residual of step k. Used as convergence criteria and to calculate beta
        rr = V_dot_v(R, R);
        //Status update: print every hundred iterations
        if ( k == test){
            hout << "CG iteration " << k << endl;
            test = test + 100;
        }
        //Convergence criteria
        if (sqrt(fabs(rr)) <= R0)
            break;
        //Improvement of step
        beta = rr/rr0;
        //Search direction
        P = R + P*beta;
    }
    
    if (k >= max_iter)
        hout << "CG reached maximum number of iterations" << endl;
    hout << "CG iterations: " << k << endl;
    
    //Need to put together d1 and solution X
    X.element.insert(X.element.begin(), d1.element[1]);
    X.element.insert(X.element.begin(), d1.element[0]);
    /*Print2DVec(X.element, "voltages.txt");
     Print1DVec(col_ind, "col_ind.txt");
     Print1DVec(row_ptr, "row_ptr.txt");
     Print1DVec(values, "values.txt");
     Print1DVec(diagonal, "diagonal.txt");//*/
    
    return X;
}

//Multiplication of a symmetric sparse matrix and a vector using SSS
MathMatrix RNetwork::spM_V_SSS(MathMatrix V, vector<long int> rowptr, vector<long int> colind, vector<double> diagonal, vector<double> values)
{
    //Size of the system
    long int N = V.element.size();
    //Initialize result vector
    MathMatrix R(N,1);
    //SSS
    long int c;
    for (long int r = 0; r < N; r++) {
        R.element[r][0] = diagonal[r]*V.element[r][0];
        for (long int j = rowptr[r]; j < rowptr[r+1]; j++) {
            c = colind[j];
            R.element[r][0] = R.element[r][0] + values[j]*V.element[c][0];
            R.element[c][0] = R.element[c][0] + values[j]*V.element[r][0];
        }
    }
    
    return R;
}

//This function assumes that the two vectores are of size (n,1)
double RNetwork::V_dot_v(MathMatrix v1, MathMatrix v2)
{
    //initialize sum
    double sum = 0;
    
    //calculate dot product
    for (long int i = 0; i < v1.element.size(); i++)
        sum += v1.element[i][0]*v2.element[i][0];
    
    return sum;
}

int RNetwork::Remove_currentless_CNTs(MathMatrix voltages, vector<int> &cluster, vector<int> LM_matrix, vector<int> &dead, vector<double> &cnts_radius, vector<vector<long int> > &structure)
{
    //Variable for the current
    double I;
    //Variable to compare to zero
    double zero = 1.0E-8;
    //Variables
    int CNT;
    long int P1, node1, P2, node2;
    
    //Scan all the CNTs that belong to the current cluster
    for (long int i = 0; i < cluster.size(); i++) {
        CNT = cluster[i];
        //The dead branches are in th extremes of a CNT, they cannot be in the middle.
        //Hence I just need to check the current of the first and last elements of the CNT
        
        //hout << "CheckFront " <<endl;
        long int j, skip_back = 0;
        //Check all elements from the begining until one that carries current is found (if any)
        for (j = 0; j < contacts_cnt_point[CNT].size()-1; j++) {
            P1 = contacts_cnt_point[CNT][j];
            node1 = LM_matrix[P1];
            P2 = contacts_cnt_point[CNT][j+1];
            node2 = LM_matrix[P2];
            I = abs(voltages.element[node2].front() - voltages.element[node1].front());
            //hout << "CNT="<<CNT<<" P1="<<P1<<" node1="<<node1<<" P2="<<P2<<" node2="<<node2<<" I="<<I<<endl;
            if (I > zero)
                break;
        }
        //The CNT will be trimmed depending on the value of j and the last calculated I
        if (j == contacts_cnt_point[CNT].size()-1) {
            //If the loop reached the last element I need to check if it carries any current or not
            //If the last element carries current then it stays, otherwise all the CNT is trimmed
            if (I > zero) {
                P1 = contacts_cnt_point[CNT].front();
                P2 = contacts_cnt_point[CNT][j];
                node2 = LM_matrix[P2];
                Trim_currenles_CNT(structure, cnts_radius, cluster, dead, CNT, P1, P2);
                //hout << "TrimF1 j="<<j<<" size-1="<<(int)contacts_cnt_point[CNT].size()-1<<endl;
            } else {
                P1 = contacts_cnt_point[CNT].front();
                P2 = contacts_cnt_point[CNT].back();
                Trim_currenles_CNT(structure, cnts_radius, cluster, dead, CNT, P1, P2);
                //Since all the CNT has to be removed, the size of the cluster changes, so I need to
                //adjust the counter i because of that
                i--;
                //hout << "TrimF2 j="<<j<<" size-1="<<(int)contacts_cnt_point[CNT].size()-1<<endl;
            }
            //If this part is reached then:
            //1) there is onle element left in the CNT and it carries current
            //2) the whole CNT does not carry any current
            //So in any case there is no need to look for a dead branch from teh back so the
            //skip_back flag is set to 1
            skip_back = 1;
        } else if (j > 0) {
            //If j i zero that means that the first element carries current so it stays and there is
            //nothing to be done. only when j is non zero and did not reached the last element there
            //is only a segment of the CNT to be trimmed
            P1 = contacts_cnt_point[CNT].front();
            P2 = contacts_cnt_point[CNT][j];
            node2 = LM_matrix[P2];
            Trim_currenles_CNT(structure, cnts_radius, cluster, dead, CNT, P1, P2);
            //hout << "TrimF3 j="<<j<<" size-1="<<(int)contacts_cnt_point[CNT].size()-1<<endl;
        }
        
        //hout << "CheckBack " << endl;
        //Get the nodes of the back element
        //If the last element carries no current, then check how many of the previous elements
        //are not carrying any current either
        if (!skip_back){
            for (j = contacts_cnt_point[CNT].size()-1; j >= 1; j--) {
                //Note: although I check again the last element, I do so in case there is only one element
                //so to avoid a segmentation fault
                P1 = contacts_cnt_point[CNT][j-1];
                node1 = LM_matrix[P1];
                P2 = contacts_cnt_point[CNT][j];
                node2 = LM_matrix[P2];
                I = abs(voltages.element[node2].front() - voltages.element[node1].front());
                //hout << "CNT="<<CNT<<" P1="<<P1<<" node1="<<node1<<" P2="<<P2<<" node2="<<node2<<" I="<<I<<endl;
                if (I > zero)
                    break;
            }
            //There is a dead branch at the back of the CNT if j < contacts_cnt_point[CNT].size()-1,
            //or equivalently if j != contacts_cnt_point[CNT].size()-1, as it was the value it
            //was initialized with
            if (j != contacts_cnt_point[CNT].size()-1) {
                //There is at least one conducting element at the front of the CNT, so when
                //checking from the back there could never be the case when the whole CNT has
                //to be deleted. If the loop reaches the frontal element for sure it carries current
                P1 = P2;
                P2 = contacts_cnt_point[CNT].back();
                Trim_currenles_CNT(structure, cnts_radius, cluster, dead, CNT, P1, P2);
            }
        }
    }
    
    //hout << "EndCheck ";
    return 1;
}

//Once a dead branch is found, this function removes it
int RNetwork::Trim_currenles_CNT(vector<vector<long int> > &structure, vector<double> &cnts_radius, vector<int> &cluster, vector<int> &dead, int CNT, long int P1, long int P2)
{
    //Empty vector
    vector<long int> empty;
    
    //The dead branches can only be on extremes of a CNT, so I check if P1 is the front of the CNT
    //or if P2 is its back.
    if (P1 == structure[CNT].front()) {
        if (P2 == structure[CNT].back()) {
            //hout << "TrimAll ";
            //In this case all the CNT is to be trimmed
            //Remove CNT from cluster
            Remove_from_vector(CNT, cluster);
            //Add the CNT to the dead vector
            dead.push_back(CNT);
        } else {
            //hout << "TrimFront ";
            //hout << " CNT=" << CNT << " P1=" << P1 << " P2=" << P2 << " cnt.front()=" << structure[CNT].front();
            //hout << " cnt.back()=" << structure[CNT].back() << ' ';
            //If P1 is the front, then I calculate the index of P2 inside the CNT. For this I use the
            //second point of the CNT, which for sure has a consecutive numbering. Besides, because of the
            //if-statement above, for sure the second point is not the last point of the CNT so it follows
            //the consecuteve numbering of the structure
            long int i = P2 - structure[CNT][1] + 1;
            //Add the new CNT number to the dead_braches vector
            dead.push_back(structure.size());
            //A new CNT is created so it has to be added to the cnts_inside vector
            cnts_inside.push_back(structure.size());
            //Cut the CNT and add it to the back of the structure as a new CNT
            structure.push_back(empty);
            structure.back().insert(structure.back().begin(), structure[CNT].begin(), structure[CNT].begin()+i+1);
            //Here I need to keep the last point of the CNT segment in both the new and old CNTs
            structure[CNT].erase(structure[CNT].begin(), structure[CNT].begin()+i);
            //Add the radius of the new CNT, which is the same as the old one
            cnts_radius.push_back(cnts_radius[CNT]);
        }
        //hout << "Trimmed " << endl;
        return 1;
    } else if (P2 == structure[CNT].back()){
        //hout << "TrimBack ";
        //hout << " CNT=" << CNT << " P1=" << P1 << " P2=" << P2 << " cnt.front()=" << structure[CNT].front();
        //hout << " cnt.back()=" << structure[CNT].back() << ' ';
        //If P2 is the back, then I calculate the index of P1 inside the CNT. For this I use the
        //second point of the CNT, which for sure has a consecutive numbering. Besides, because when
        //chacking from the back the whole CNT is never trimmed, for sure the first point is not the
        //first point of the CNT. So it follows the consecuteve numbering of the structure
        long int i = P1 - structure[CNT][1] + 1;
        //Add the new CNT number to the dead_braches vector
        dead.push_back(structure.size());
        //A new CNT is created so it has to be added to the cnts_inside vector
        cnts_inside.push_back(structure.size());
        //Cut the CNT and add it to the back of the structure as a new CNT
        structure.push_back(empty);
        structure.back().insert(structure.back().begin(), structure[CNT].begin()+i, structure[CNT].end());
        //Here I need to keep the first point of the CNT segment in both the new and old CNTs
        structure[CNT].erase(structure[CNT].begin()+i+1, structure[CNT].end());
        //Add the radius of the new CNT, which is the same as the old one
        cnts_radius.push_back(cnts_radius[CNT]);
        //hout << "Trimmed " << endl;
        return 1;
    } else {
        //If it reaches here there is an error because a middle fragment of a CNT cannot be currentless
        hout << "Currenless segment does not start at the front nor at the back of CNT.";
        hout << " CNT=" << CNT << " P1=" << P1 << " P2=" << P2 << " cnt.front()=" << structure[CNT].front();
        hout << " cnt.back()=" << structure[CNT].back() << endl;
        return 0;
    }
}

//This function calculates the lengths of the CNTs in each cluster
int RNetwork::Clusters_length(vector<Point_3D> points, vector<vector<long int> > structure)
{
    //Set the correct size of the vector that will store the lengths of the CNTs in each cluster
    clusters_lengths.resize(7,0);
    
    //These variables are just to make the code easier to read
    int CNT;
    long int point1, point2;
    //Here I will save the total length of percolating clusters
    double total_percolating = 0;
    //If directional clusters has size zero then, total_percolating stays in zero
    for (int i = 0; i < (int)directional_clusters.size(); i++) {
        for (int j = 0; j < (int)directional_clusters[i].size(); j++) {
            CNT = directional_clusters[i][j];
            for (int k = 0; k < (int)structure[CNT].size()-1; k++) {
                point1 = structure[CNT][k];
                point2 = structure[CNT][k+1];
                clusters_lengths[i]+= points[point1].distance_to(points[point2]);
            }
        }
        total_percolating += clusters_lengths[i];
    }
    
    //Calculate the lengths of the dead branches
    double dead_length = 0;
    dead_branches_lengths.resize(dead_branches.size(),0);
    for (int i = 0; i < (int)dead_branches.size(); i++) {
        for (int j = 0; j < (int)dead_branches[i].size(); j++) {
            CNT = dead_branches[i][j];
            for (int k = 0; k < (int)structure[CNT].size()-1; k++) {
                point1 = structure[CNT][k];
                point2 = structure[CNT][k+1];
                dead_branches_lengths[i]+= points[point1].distance_to(points[point2]);
            }
        }
        dead_length += dead_branches_lengths[i];
    }
    //Add the length of dead branches to the length vector
    clusters_lengths.push_back(dead_length);
    
    //Calculate the lengths of the non percolating clusters
    double total_iso_length = 0;
    for (int i = 0; i < (int)isolated.size(); i++) {
        for (int j = 0; j < (int) isolated[i].size(); j++) {
            CNT = isolated[i][j];
            for (int k = 0; k < (int) structure[CNT].size()-1; k++) {
                point1 = structure[CNT][k];
                point2 = structure[CNT][k+1];
                total_iso_length += points[point1].distance_to(points[point2]);
            }
        }
    }
    //Add the length of isolated clusters to the length vector
    clusters_lengths.push_back(total_iso_length);
    
    //Add the sum of lengths of isolated CNTS and dead branches to the length vector
    total_iso_length += dead_length;
    clusters_lengths.push_back(total_iso_length);
    
    //calculate the total length of CNTs inside the box. This is done as a check of the code.
    double total_length = 0;
    for (int i = 0; i < (int)cnts_inside.size(); i++) {
        CNT = cnts_inside[i];
        for (int j = 0; j < (int)structure[CNT].size()-1; j++) {
            point1 = structure[CNT][j];
            point2 = structure[CNT][j+1];
            total_length += points[point1].distance_to(points[point2]);
        }
    }
    //Add the total length of CNTs inside the box to the length vector
    clusters_lengths.push_back(total_length);
    
    //Calculate the fractions of each cluster
    clusters_fractions.resize(clusters_lengths.size(), 0);
    for (int i = 0; i < (int)clusters_lengths.size(); i++) {
        clusters_fractions[i] = clusters_lengths[i]/total_length;
    }
    
    //Print1DVec(clusters_lengths, "clusters_lengths_.txt");
    //Print1DVec(dead_branches_lengths, "dead_branches_lengths_.txt");
    //Print1DVec(clusters_fractions, "clusters_fractions_.txt");
    Append1DVec(clusters_lengths, "clusters_lengths.txt");
    Append1DVec(dead_branches_lengths, "dead_branches_lengths.txt");
    Append1DVec(clusters_fractions, "clusters_fractions.txt");
    
    return 1;
}

//This fuctions exports all tecplot files corresponding to each cluster, dead branches and non percolating clusters
int RNetwork::Export_visualization_files(const struct RVE_Geo &cell_geo, vector<Point_3D> points_in, vector<double> cnts_radius, vector<vector<long int> > structure)
{
    //Save a separate file for each percolating direction
    vector<string> filenames;
    filenames.push_back("SingleZone_00_xx.dat");
    filenames.push_back("SingleZone_01_yy.dat");
    filenames.push_back("SingleZone_02_zz.dat");
    filenames.push_back("SingleZone_03_xx_yy.dat");
    filenames.push_back("SingleZone_04_xx_zz.dat");
    filenames.push_back("SingleZone_05_yy_zz.dat");
    filenames.push_back("SingleZone_06_xx_yy_zz.dat");
    //Save a separate file for the dead branches of each percolating direction
    filenames.push_back("SingleZoneDead_00_xx.dat");
    filenames.push_back("SingleZoneDead_01_yy.dat");
    filenames.push_back("SingleZoneDead_02_zz.dat");
    filenames.push_back("SingleZoneDead_03_xx_yy.dat");
    filenames.push_back("SingleZoneDead_04_xx_zz.dat");
    filenames.push_back("SingleZoneDead_05_yy_zz.dat");
    filenames.push_back("SingleZoneDead_06_xx_yy_zz.dat");
    
    for (int i = 0; i < (int)directional_clusters.size(); i++){
        //Check if the family is non empty. If it is non empty then a visualization file can be created
        if (directional_clusters[i].size())
            if(!Export_cnt_networks_meshes(cell_geo, structure, directional_clusters[i], points_in, cnts_radius, filenames[i])) {
                hout << "Error while translating and exporting directional clusters" <<endl;
                return 0;
            }
        //hout << " i+7=" << i+7;
        //It is possible that a CNT spans from one boundary to the other and has no contacts. In this case
        //The CNT is a cluster itself and has no dead branches, so I need to check separately if dead_branches[i]
        //is non empty. That is, the fact that a cluster percolates does not mean that there will be dead branches
        if (dead_branches[i].size())
            if(!Export_cnt_networks_meshes(cell_geo, structure, dead_branches[i], points_in, cnts_radius, filenames[i+7])) {
                hout << "Error while translating and exporting directional clusters" <<endl;
                return 0;
            }
    }
    hout << "Tecplot files exported for percolating clusters." << endl;
    
    //Gather all isolated clusters into one so it goes into the exporting functions
    vector<int> all_iso;
    for (int i = 0; i < (int)isolated.size(); i++)
        all_iso.insert(all_iso.begin(), isolated[i].begin(), isolated[i].end());
    //Once all isolated CNTs are in one single cluster, export it as a single file
    if(!Export_cnt_networks_meshes(cell_geo, structure, all_iso, points_in, cnts_radius, "SingleZone_isolated.dat")) {
        hout << "Error while translating and exporting directional clusters" <<endl;
        return 0;
    }
    
    
    return 1;
}//

//This fuction updates the box geometry and related variables
void RNetwork::Update_box_geometry()
{
    //size
    box_geometry.len_x = box_geometry.len_x + dwindow;
    box_geometry.wid_y = box_geometry.wid_y + dwindow;
    box_geometry.hei_z = box_geometry.hei_z + dwindow;
    //reference point
    box_geometry.poi_min.x = box_geometry.poi_min.x - dwindow/2;
    box_geometry.poi_min.y = box_geometry.poi_min.y - dwindow/2;
    box_geometry.poi_min.z = box_geometry.poi_min.z - dwindow/2;
    //These variables are to store the size of the RVE and reduces operations when accessing them
    lx = box_geometry.len_x;
    ly = box_geometry.wid_y;
    lz = box_geometry.hei_z;
    //These variables are the coordinates of the lower corner of the RVE that defines its geometry
    xmin = box_geometry.poi_min.x;
    ymin = box_geometry.poi_min.y;
    zmin = box_geometry.poi_min.z;
    //These variables are the coordinates of the upper corner of the RVE that defines its geometry
    xmax = lx + xmin;
    ymax = ly + ymin;
    zmax = lz + zmin;
}


//===================================================================================================
//===================================================================================================
//===================================================================================================
//Funtions that print

//Print a vector of 3D points
void RNetwork::Print1DVec(const vector<Point_3D> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < list.size(); i++) {
        otec << list[i].x << "\t" << list[i].y << "\t" << list[i].z << "\t" << list[i].flag << "\n";
    }
    otec.close();
}

//Print a vector of chars with the specified filename
void RNetwork::Print1DVec(const vector<char> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of integers with the specified filename
void RNetwork::Print1DVec(const vector<int> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    //otec << "Positions \n";
    for (long int i=0; i < list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of doubles with the specified filename
void RNetwork::Print1DVec(const vector<double> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of doubles with the specified filename
void RNetwork::Append1DVec(const vector<double> &list, const string &filename)
{
    ofstream otec(filename.c_str(), std::ios_base::app);
    hout << "Appending to file: " << filename << "\n";
    for (long int i=0; i < list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\t";
    }
    otec << "\n";
    otec.close();
}

//Print a vector of long integers with the specified filename
void RNetwork::Print1DVec(const vector<long int> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of vectors of integers with the specified filename
void RNetwork::Print2DVec(const vector<vector<int> > &num_mat, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < num_mat.size(); i++) {
        for (long int j = 0; j < num_mat[i].size(); j++) {
            otec << num_mat[i][j] << '\t' ;
        }
        otec << '\n' ;
    }
    otec.close();
}

//Print a vector of vectors of integers with the specified filename
void RNetwork::Print2DVec(const vector<vector<long int> > &num_mat, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < num_mat.size(); i++) {
        for (long int j = 0; j < num_mat[i].size(); j++) {
            otec << num_mat[i][j] << '\t' ;
        }
        otec << '\n' ;
    }
    otec.close();
}

//Print a vector of vectors of doubles with the specified filename
void RNetwork::Print2DVec(const vector<vector<double> > &num_mat, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < num_mat.size(); i++) {
        for (long int j = 0; j < num_mat[i].size(); j++) {
            otec << num_mat[i][j] << '\t' ;
        }
        otec << '\n' ;
    }
    otec.close();
}

//This function has to be declared every time in every class
string RNetwork::Get_Line(ifstream &infile)const
{
	string s;
	//∂¡»Î–≈œ¢“ª––
	getline(infile,s);
	//Ã¯π˝◊¢ Õ––
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}

//===================================================================================================
//===================================================================================================
//===================================================================================================
//Just copy paste of some functions that I need but are defined as private in GeoNano
//I also made some modifications so the functions fit the variables I have in this class
int RNetwork::Export_cnt_networks_meshes(const struct RVE_Geo &cell, vector<vector<long int> > structure, vector<int> cluster, vector<Point_3D> points_in, vector<double> cnts_radius, const string &filename)
{
	//…˙≥…”√”⁄±Ì æƒ…√◊π‹µƒΩ⁄µ„º∞Àƒ√ÊÃÂÕ¯∏Ò
	vector<vector<Node> > cnts_nodes;
	vector<vector<Element> > cnts_eles;
    
    //hout << "Export1 ";
	if(Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, structure, cluster, points_in, cnts_radius)==0) return 0;
    
	// ‰≥ˆƒ…√◊π‹œﬂÕ¯∏Òµ•Zone in Tecplot
    //hout << "Export2 ";
    if(Export_cnts_meshes_singlezone(cell, cnts_nodes, cnts_eles,filename)==0) return 0;
    
	return 1;
}


int RNetwork::Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, vector<vector<long int> > structure, vector<int> cluster, vector<Point_3D> points_in, vector<double> cnts_radius)
{
    //Variables for the points
    long int Pij, Pij_m, Pij_p; //Pij = [i][j], Pij_m = [i][j-1], Pij_p = [i][j+1]
    int CNT;
	//—≠ª∑ƒ…√◊π‹œﬂ
	for(int i=0; i<(int)cluster.size(); i++){
		vector<Node> nod_temp;
		vector<Element> ele_temp;
        
        //Current CNT of the cluster
        CNT = cluster[i];
        
		const int cps = (int)structure[CNT].size();
        //hout << "CNT="<<CNT<<" cps="<<cps<<endl;
		for(int j=0; j<cps; j++){
            //Set the current point
            Pij = structure[CNT][j];
			//µ√µΩ∆Ω√Êµƒ∑®œÚœÚ¡ø
			Point_3D plane_normal;
			if(j==0){
                Pij_p = structure[CNT][j+1];
                plane_normal.x = points_in[Pij].x - points_in[Pij_p].x;
				plane_normal.y = points_in[Pij].y - points_in[Pij_p].y;
				plane_normal.z = points_in[Pij].z - points_in[Pij_p].z;
                
			} else if(j==cps-1){
                Pij_m = structure[CNT][j-1];
				plane_normal.x = points_in[Pij_m].x - points_in[Pij].x;
				plane_normal.y = points_in[Pij_m].y - points_in[Pij].y;
				plane_normal.z = points_in[Pij_m].z - points_in[Pij].z;
			} else{
                Pij_p = structure[CNT][j+1];
                Pij_m = structure[CNT][j-1];
				const Point_3D vect[3] = { points_in[Pij_m], points_in[Pij], points_in[Pij_p] };
				const double A = pow(vect[0].x-vect[1].x,2)+pow(vect[0].y-vect[1].y,2)+pow(vect[0].z-vect[1].z,2);
				const double B = pow(vect[2].x-vect[1].x,2)+pow(vect[2].y-vect[1].y,2)+pow(vect[2].z-vect[1].z,2);
				const double tt = sqrt(A/B);
                
				//◊¯±Í±‰ªªªπ‘≠Ωªµ„Œª÷√
				double x, y, z;
				x=vect[1].x+tt*(vect[2].x-vect[1].x);
				y=vect[1].y+tt*(vect[2].y-vect[1].y);
				z=vect[1].z+tt*(vect[2].z-vect[1].z);
                
				plane_normal.x = vect[0].x - x;
				plane_normal.y = vect[0].y - y;
				plane_normal.z = vect[0].z - z;
			}
			//µ√µΩ∆Ω√Ê…œµƒµ„
			Point_3D plane_center = points_in[Pij];
            
			//∂®“Â‘≤…œ∑÷µƒ∑› ˝
			const int num_sec = 36;
			if(j==0){
				double normal_sita, normal_pha;  //º–Ω«
				//µ√µΩœÚ¡øplane_normal‘⁄«Ú◊¯±Íœµ÷–µƒº–Ω«
				if(Get_angles_vector_in_spherial_coordinates(plane_normal, normal_sita, normal_pha)==0) return 0;
                
				//µ√µΩµ„∑® ∏∆Ω√Ê…œ(“‘µ„Œ™÷––ƒ“‘“ª∂®÷µŒ™∞Îæ∂)“ª∏ˆ‘≤ª∑…œµƒ“ª◊Èµ„, ∑® ∏Ω«∂»(«Ú√Ê◊¯±Í)“—÷™
				if(Get_points_circle_in_plane(plane_center, normal_sita, normal_pha, cnts_radius[i], num_sec, nod_temp)==0) return 0;
			} else{
				//º∆À„«∞“ª∏ˆ‘≤…œµƒµ„—ÿ«∞“ª∏ˆœﬂ∂Œ∑ΩœÚ(line_vec)‘⁄π˝µ„plane_center∑® ∏plane_normalµƒ∆Ω√Ê…œµƒÕ∂”∞µ„
				Point_3D line_vec;
				line_vec.x = points_in[Pij_m].x - points_in[Pij].x;
				line_vec.y = points_in[Pij_m].y - points_in[Pij].y;
				line_vec.z = points_in[Pij_m].z - points_in[Pij].z;
				if(Get_projected_points_in_plane(plane_center, plane_normal, line_vec, num_sec, nod_temp)==0) return 0;
			}
            
			//…˙≥…“ª◊Èµ•‘™œÚ¡ø
			if(j!=0){
				int nodes_num[6];
				nodes_num[0] = (j-1)*(num_sec+1);   //÷––ƒµ„±‡∫≈
				nodes_num[3] = j*(num_sec+1);
				for(int k=1; k<=num_sec; k++){
					nodes_num[1] = (j-1)*(num_sec+1) + k;
					nodes_num[2] = (j-1)*(num_sec+1) + 1 + k%num_sec;
					nodes_num[4] = j*(num_sec+1) + k;
					nodes_num[5] = j*(num_sec+1) + 1 + k%num_sec;
                    
					Element eles_num[3];
					//----------------------------------------------------------------
					//≤Â»ÎΩ⁄µ„
					eles_num[0].nodes_id.push_back(nodes_num[0]);
					eles_num[0].nodes_id.push_back(nodes_num[1]);
					eles_num[0].nodes_id.push_back(nodes_num[2]);
					eles_num[0].nodes_id.push_back(nodes_num[3]);
                    
					eles_num[1].nodes_id.push_back(nodes_num[1]);
					eles_num[1].nodes_id.push_back(nodes_num[2]);
					eles_num[1].nodes_id.push_back(nodes_num[3]);
					eles_num[1].nodes_id.push_back(nodes_num[5]);
                    
					eles_num[2].nodes_id.push_back(nodes_num[1]);
					eles_num[2].nodes_id.push_back(nodes_num[3]);
					eles_num[2].nodes_id.push_back(nodes_num[4]);
					eles_num[2].nodes_id.push_back(nodes_num[5]);
					//----------------------------------------------------------------
					//≤Â»Îµ•‘™
					ele_temp.push_back(eles_num[0]);
					ele_temp.push_back(eles_num[1]);
					ele_temp.push_back(eles_num[2]);
				}
			}
		}
        
		//≤Â»Î◊‹œÚ¡ø
		nodes.push_back(nod_temp);
		eles.push_back(ele_temp);
	}
    //hout << "nodes " << nodes.size() << " eles " << eles.size() << endl;
	return 1;
}


int RNetwork::Export_cnts_meshes_singlezone(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const string &filename)const
{
	//ƒ…√◊π‹œﬂµƒ ˝¡ø
	int cnts_account = (int)nodes.size();
	if(cnts_account!=(int)eles.size()) { hout << "Ω⁄µ„∫Õµ•‘™œ‘ æµƒƒ…√◊π‹œﬂ ˝¡ø≤ª“ª÷¬£° «ÎºÏ≤È£°" << endl; return 0; }
	
	ofstream otec(filename.c_str());
	otec << "TITLE = CNT_Meshes_Singlezone" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;
	
	//---------------------------------------------------------------------------
	// ‰≥ˆµ•∞˚øÚº‹
	otec << "ZONE N=" << 8 << ", E=" << 1 << ", F=FEPOINT, ET=BRICK" << endl;
	double cell_x[2] = {cell.poi_min.x, cell.poi_min.x+cell.len_x};
	double cell_y[2] = {cell.poi_min.y, cell.poi_min.y+cell.wid_y};
	double cell_z[2] = {cell.poi_min.z, cell.poi_min.z+cell.hei_z};
	for(int i=0; i<2; i++)
		for(int j=0; j<2; j++)
			for(int k=0; k<2; k++)
			{
				otec << cell_x[i] << "  " << cell_y[j] << "  " << cell_z[k] << endl;
			}
	
	otec << "1 2 4 3 5 6 8 7" << endl;
	otec << endl << endl;
    
	//---------------------------------------------------------------------------
	// ‰≥ˆƒ…√◊π‹Õ¯∏Ò
    
	int nodes_num = 0;
	int eles_num = 0;
    
	for(int i=0; i<cnts_account; i++)
	{
		nodes_num +=  (int)nodes[i].size();
		eles_num += (int)eles[i].size();
	}
    
	otec << "ZONE N=" << nodes_num << ", E=" << eles_num << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
    //hout << "\nNodes:" << endl;
	for(int i=0; i<cnts_account; i++)
	{
        //hout << "\ni= " << i <<"\nj: \t";
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
            //hout << j << ' ';
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
	}
	otec << endl;
    
    //hout << "\nElements" << endl;
	nodes_num = 0;
	for(int i=0; i<cnts_account; i++)
	{
        //hout << "\ni= " << i <<"\nj: \t";
		if(i!=0)  nodes_num +=  (int)nodes[i-1].size();
		for(int j=0; j<(int)eles[i].size(); j++)
		{
            //hout << j << ' ';
			otec	<< eles[i][j].nodes_id[0]+1+nodes_num << "  " << eles[i][j].nodes_id[1]+1+nodes_num << "  "
            << eles[i][j].nodes_id[2]+1+nodes_num << "  " << eles[i][j].nodes_id[3]+1+nodes_num << endl;
		}
	}
    
	otec.close();
    
	return 1;
}

int RNetwork::Get_angles_vector_in_spherial_coordinates(const Point_3D &normal, double &sita, double &pha)const
{
	if(normal.x==0&&normal.y==0&&normal.z==0) { hout << "œÚ¡øµƒ»˝∏ˆ∑÷¡ø∂ºŒ™¡„£° «ÎºÏ≤È£°" << endl; return 0; }
	sita =  acos(normal.z/sqrt(normal.x*normal.x+normal.y*normal.y+normal.z*normal.z));
	if(normal.x==0&&normal.y==0) pha = 0;
	else if(normal.y>=0) pha = acos(normal.x/sqrt(normal.x*normal.x+normal.y*normal.y));
	else if(normal.y<0) pha = 2*PI - acos(normal.x/sqrt(normal.x*normal.x+normal.y*normal.y));
    
	return 1;
}
//---------------------------------------------------------------------------
int RNetwork::Get_points_circle_in_plane(const Point_3D &center, const double &trans_sita, const double &trans_pha, const double &radius, const int &num_sec, vector<Node> &nod_temp)const
{
	//≤Â»Î÷––ƒµ„
	Node new_node(center.x, center.y, center.z);
	nod_temp.push_back(new_node);
    
	//∂®“Â±‰ªªæÿ’Û
	MathMatrix trans_mat(3,3);
	trans_mat = Get_transformation_matrix(trans_sita, trans_pha);
    
	//“ªŒ¨œÚ¡ø
	MathMatrix Rvec(3,1);
	Rvec.element[0][0] = 0;
	Rvec.element[1][0] = 0;
	Rvec.element[2][0] = radius;
    
	//“ªŒ¨œÚ¡ø
	MathMatrix Res(3,1);
    
	double sita, pha;
	sita = 0.5*PI; //‘⁄xy∆Ω√Ê…œ‘À∂Ø
	for(int i=0; i<num_sec; i++)
	{
		pha = i*2*PI/num_sec;
		MathMatrix matrix_temp = trans_mat*Get_transformation_matrix(sita, pha);
		Res = matrix_temp*Rvec;
        
		new_node.x = center.x + Res.element[0][0];
		new_node.y = center.y + Res.element[1][0];
		new_node.z = center.z + Res.element[2][0];
        
		//≤Â»Î‘≤ª∑…œµƒµ„
		nod_temp.push_back(new_node);
	}
	
	return 1;
}
//---------------------------------------------------------------------------
int RNetwork::Get_projected_points_in_plane(const Point_3D &center, const Point_3D &normal, const Point_3D &line, const int &num_sec, vector<Node> &nod_temp)const
{
	//º«¬º«∞“ª¥Œ…˙≥…µƒ◊‹Ω⁄µ„∏ˆ ˝
	const int nod_size = (int)nod_temp.size();
    
	//≤Â»Î÷––ƒµ„
	Node new_node(center.x, center.y, center.z);
	nod_temp.push_back(new_node);
    
	const double vectors_dot_product = normal.x*line.x+normal.y*line.y+normal.z*line.z;
	if(vectors_dot_product==0.0)
	{
		hout << "¡Ω∏ˆ∑®œÚ¡ø≥…÷±Ω«(∂‘”¶»˝∏ˆµ„0°¢1°¢2, ÷ª”–µ±Ω«012–°”⁄90∂»£¨«“œﬂ∂Œ12>œﬂ∂Œ01 ±, ≤≈ø…ƒ‹≥ˆœ÷’‚÷÷«Èøˆ£¨∂¯’‚ª˘±æ…œ «≤ªø…ƒ‹µƒ)£° «ÎºÏ≤È£°" << endl;
		return 0;
	}
    
	for(int i=num_sec; i>0; i--)
	{
		Point_3D point(center.x-nod_temp[nod_size-i].x, center.y-nod_temp[nod_size-i].y, center.z-nod_temp[nod_size-i].z);
		new_node.x = nod_temp[nod_size-i].x + (normal.x*point.x+normal.y*point.y+normal.z*point.z)*line.x/ vectors_dot_product;
		new_node.y = nod_temp[nod_size-i].y + (normal.x*point.x+normal.y*point.y+normal.z*point.z)*line.y/ vectors_dot_product;
		new_node.z = nod_temp[nod_size-i].z + (normal.x*point.x+normal.y*point.y+normal.z*point.z)*line.z/ vectors_dot_product;
        
		//≤Â»Î‘≤ª∑…œµƒµ„
		nod_temp.push_back(new_node);
	}
    
	return 1;
}

MathMatrix RNetwork::Get_transformation_matrix(const double &sita, const double &pha)const
{
	//cnt_sitaΩ«±‰ªªæÿ’Û
	MathMatrix Msita(3,3);
	Msita.element[0][0] = cos(sita);
	Msita.element[0][2] = sin(sita);
	Msita.element[1][1] = 1;
	Msita.element[2][0] = -sin(sita);
	Msita.element[2][2] = cos(sita);
    
	//cnt_phaΩ«±‰ªªæÿ’Û
	MathMatrix Mpha(3,3);
	Mpha.element[0][0] = cos(pha);
	Mpha.element[0][1] = -sin(pha);
	Mpha.element[1][0] = sin(pha);
	Mpha.element[1][1] = cos(pha);
	Mpha.element[2][2] = 1;
    
	return Mpha*Msita;
}
