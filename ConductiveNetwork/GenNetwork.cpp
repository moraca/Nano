//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	GenNetwork.cpp
//OBJECTIVE:	To generate networks with overlapping
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "GenNetwork.h"

//Generate 3D nantube networks with ovelapping
int GenNetwork::Generate_nanotube_networks(const struct Geom_RVE &geom_rve, const struct Cluster_Geo &clust_geo, const struct Nanotube_Geo &nanotube_geo, const struct Cutoff_dist &cutoffs, vector<Point_3D> &cpoints, vector<double> &cnts_radius, vector<vector<long int> > &cstructures)const
{
    //Define a two-dimensional vector of three-dimensional points for storing the CNT threads
    vector<vector<Point_3D> > cnts_points;
    
    //Generate a network defined by points and connections
    //Use rand() for the random number generation
    //if(Generate_network_threads(geom_rve, clust_geo, nanotube_geo, cutoffs, cnts_points, cnts_radius)==0) return 0;
    //Use the Mersenne Twister for the random number generation
    if (Generate_network_threads_mt(geom_rve, clust_geo, nanotube_geo, cutoffs, cnts_points, cnts_radius)==0) return 0;
    
    //Checking the angle between two segments in one nanotube (if less than PI/2, provide an alarm)
    if(CNTs_quality_testing(cnts_points)==0) return 0;
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Transform the 2D cnts_points into 1D cpoints and 2D cstructuers
    if(Transform_cnts_points(cnts_points, cpoints, cstructures)==0) return 0;
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    /*/A new class of Tecplot_Export
     Tecplot_Export *Tecexpt = new Tecplot_Export;
     
     struct cuboid cub;														//Generate a cuboid for RVE
     cub.poi_min = geom_rve.ex_origin;
     cub.len_x = geom_rve.ex_len;
     cub.wid_y = geom_rve.ey_wid;
     cub.hei_z = geom_rve.ez_hei;
     
     //The geometric structure of CNT network (by threads in Tecplot)
     if(Tecexpt->Export_network_threads(cub, cnts_points)==0) return 0;
     
     //The geometric structure of CNT network (by tetrahedron meshes in Tecplot) //Attention: little parts of nanotube volumes out of the cuboid
     if(Tecexpt->Export_cnt_network_meshes(cub, cnts_points, cnts_radius)==0) return 0;//*/
    
    return 1;
}

//Generate a network defined by points and connections
//Use rand() for the random number generation
int GenNetwork::Generate_network_threads(const struct Geom_RVE &geom_rve, const struct Cluster_Geo &clust_geo, const struct Nanotube_Geo &nanotube_geo, const struct Cutoff_dist &cutoffs,
                                         vector<vector<Point_3D> > &cnts_points,  vector<double> &cnts_radius)const
{
    //Generate random seed in terms of local time
    //unsigned int time_seed = 1453384844;
    unsigned int time_seed = (unsigned int)time(NULL);
    hout << "Time seed "<<time_seed<<endl;
    srand(time_seed);
    //srand((unsigned int)time(NULL));
    
    //---------------------------------------------------------------------------
    //Generate the data for nanotube clusters limited in ellipsoid surfaces (each ellipsoid is within the RVE and is separated with each other)
    if(clust_geo.vol_fra_criterion>0.0)
    {
        struct cuboid cub;										//Generate a cuboid for RVE
        cub.poi_min = geom_rve.origin;
        cub.len_x = geom_rve.len_x;
        cub.wid_y = geom_rve.wid_y;
        cub.hei_z = geom_rve.hei_z;
        cub.volume = cub.len_x*cub.wid_y*cub.hei_z;
        
        if(Get_ellip_clusters(cub, clust_geo)==0) return 0;
        //Generate a number of sperical clusters in regular arrangement
        //if(Get_spherical_clusters_regular_arrangement(cub, clust_geo)==0) return 0;
    }
    
    //---------------------------------------------------------------------------
    double vol_sum = 0;  //the sum of volume of generated CNTs
    double wt_sum = 0;   //the sum of weight of generated CNTs
    int cnt_seed_count = 0; //to record the number of generated seeds of a CNT (If the growth of a CNT fails, but this seed will be used for the next growth of a CNT)
    int cnt_reject_count = 0; //to record the number of CNTs that were deleted due to penetration
    int point_overlap_count = 0; //to record the number of times that a point had to be relocated
    int point_overlap_count_unique = 0; //to record the number of points that were overlapping other points
    
    //---------------------------------------------------------------------------
    //Define the variable MAX_ATTEMPTS
    const int MAX_ATTEMPTS = 5;
    //Vectors for handling CNT penetration
    //global_coordinates[i][0] stores the CNT number of global point i
    //global_coordinates[i][1] stores the local point number of global point i
    vector<vector<int> > global_coordinates;
    //sectioned_domain[i] contains all the points in sub-region i.
    //Sub-region i is an overlapping subregion to check for penetrations
    vector<vector<long int> > sectioned_domain;
    //n_subregions[0] is the number of subregions along x
    //n_subregions[1] is the number of subregions along y
    //n_subregions[2] is the number of subregions along z
    vector<int> n_subregions;
    //Define cutoff for overlapping
    double overlap_max_cutoff = 2*nanotube_geo.rad_max + cutoffs.van_der_Waals_dist;
    //Initialize the vector sub-regions
    Initialize_subregions(geom_rve, n_subregions, sectioned_domain);
    //This flag will be used to skip overlapping functions
    //1 = non-penetrating model
    //0 = penetrating model
    int penetrating_model_flag = 1;
    
    //Generate cuboids that represent the extended domain and the composite domain
    //To calculate the effective portion (length) which falls into the given region (RVE)
    struct cuboid gvcub;					//generate a cuboid to represent the composite domain
    gvcub.poi_min = geom_rve.origin;
    gvcub.len_x = geom_rve.len_x;
    gvcub.wid_y = geom_rve.wid_y;
    gvcub.hei_z = geom_rve.hei_z;
    gvcub.volume = geom_rve.volume;
    struct cuboid excub;					//generate a cuboid to represent the extended domain
    excub.poi_min = geom_rve.ex_origin;
    excub.len_x = geom_rve.ex_len;
    excub.wid_y = geom_rve.ey_wid;
    excub.hei_z = geom_rve.ez_hei;
    excub.volume = excub.len_x*excub.wid_y*excub.hei_z;

    
    //---------------------------------------------------------------------------
    while((nanotube_geo.criterion == "vol"&&vol_sum < nanotube_geo.real_volume)||
          (nanotube_geo.criterion == "wt"&&wt_sum < nanotube_geo.real_weight))
    {
        //---------------------------------------------------------------------------
        //Define a vector for a new nanotube
        vector<Point_3D> new_cnt;
        int new_cnt_size = (int)new_cnt.size();
        
        //---------------------------------------------------------------------------
        //Randomly generate a length of a CNT
        double cnt_length;
        if(Get_random_value(nanotube_geo.len_distrib_type, nanotube_geo.len_min, nanotube_geo.len_max, cnt_length)==0) return 0;
        //Calculate the total number of growth step for a CNT
        int step_num = (int)(cnt_length/nanotube_geo.step_length) + 1;
        
        //---------------------------------------------------------------------------
        //Randomly generate a radius of a CNT
        double cnt_rad;
        if(Get_random_value(nanotube_geo.rad_distrib_type, nanotube_geo.rad_min, nanotube_geo.rad_max, cnt_rad)==0) return 0;
        
        //---------------------------------------------------------------------------
        //Randomly generate a direction in the spherical coordinates as the intital direction of CNT segments
        double cnt_sita, cnt_pha;
        if(Get_uniform_direction(nanotube_geo, cnt_sita, cnt_pha)==0) return 0;

        MathMatrix multiplier(3,3);
        multiplier = Get_transformation_matrix(cnt_sita, cnt_pha);
        
        //---------------------------------------------------------------------------
        //The increased volume of each segement (growth step) of nanotube (Here the overlapping volume is ignored)
        const double step_vol_para = PI*cnt_rad*cnt_rad;
        //---------------------------------------------------------------------------
        //The increased weight of each segement (growth step) of nanotube (If the different radii of nanotube are considered, the linear_density may be different in every nanotube)
        const double step_wei_para = nanotube_geo.linear_density;
        
        //---------------------------------------------------------------------------
        //Randomly generate a seed (initial point) of a CNT in the extended RVE	(Comments: the seed generation is after the radius generation for the non-overlapping nanotubes generation)
        
        Point_3D cnt_poi;
        if(Get_seed_point(excub, cnt_poi)==0) return 0;
        
        //Check overlapping of the intial point
        int counter = 1;
        while (penetrating_model_flag && !Check_penetration(geom_rve, nanotube_geo, cnts_points, global_coordinates, sectioned_domain, cnts_radius, new_cnt, n_subregions, cnt_rad, cutoffs.van_der_Waals_dist, MAX_ATTEMPTS, point_overlap_count, point_overlap_count_unique, cnt_poi)) {
            if(Get_seed_point(excub, cnt_poi)==0) return 0;
            cnt_seed_count++;					//record the number of seed generations
            //hout << "Seed deleted" << endl;
            if (counter == MAX_ATTEMPTS) {
                hout << "Too many attempts to resolve overlapping of an intial CNT point (" << counter << " attempts). ";
                hout << cnt_poi.x << ' ' << cnt_poi.y << ' ' << cnt_poi.z << endl;
                return 0;
            }
            counter ++;
        }//*/
        
        new_cnt.push_back(cnt_poi);	//store this seed point in the vector for a new nanotube
        
        //---------------------------------------------------------------------------
        cnt_seed_count++;					//record the number of seed generations
        int max_seed = 1E9;
        if(cnt_seed_count>max_seed)
        {
            hout << "The number of seed genrations is lager than "<<max_seed<<", but the nanotube generation still fails to acheive the demanded volume fraction." << endl;
            return 0;
        }
        //---------------------------------------------------------------------------
        //
        
        //---------------------------------------------------------------------------
        //The growth process of nanotube
        int ellip_num = -1; //For recording the serial number of ellipsoid cluster which a nanotube penetrates out. It is no use if the cluster generation is not considered.
        for(int i=0; i<step_num; i++)
        {
            //Randomly generate a direction in the spherical coordinates
            //To have the positive Z-axis to be a central axis
            //Then, the radial angle, sita, obeys a normal distribution (sita \in fabs[(-omega,+omega)]) and the zonal angle, pha, obeys a uniform distribution (pha \in (0,2PI))
            if(Get_normal_direction(nanotube_geo.angle_max, cnt_sita, cnt_pha)==0) return 0;
            
            //To calculate the new multiplier for transformation of coordinates
            multiplier = multiplier*Get_transformation_matrix(cnt_sita, cnt_pha);
            
            //To calculate the coordinates of the new CNT point (transformation of coordinates)
            cnt_poi = cnt_poi + Get_new_point(multiplier, nanotube_geo.step_length);
            cnt_poi.flag = 1;							//1 means that point is not the intial point
            
            //---------------------------------------------------------------------------
            //If a CNT penetrates the ellipsoidal surface of a cluster from inside, the growth of this CNT will be headed back in the cluster in a probability p, (0<=p<=1).
            if(clust_geo.ellips.size()>0)
            {
                return 0;
            }
            
            //---------------------------------------------------------------------------
            //If the new CNT point grows out of the RVE, the intersecting point at the surfaces of RVE will be calculated.
            //The new segment will be cut and (the outside part will be translated into the RVE for the periodical boundary condition)
            bool touch_end = false;
            if(Judge_RVE_including_point(excub, cnt_poi)==0)
            {
                //Calculate all intersection points between the new segment and surfaces of RVE
                //(using a parametric equatio:  the parameter 0<t<1, and sort all intersection points from the smaller t to the greater t)
                vector<Point_3D> ipoi_vec;  //a vector for intersection points
                if(Get_intersecting_point_RVE_surface(excub, new_cnt.back(), cnt_poi, ipoi_vec)==0) {
                    hout << "Error in Generate_network_threads"<<endl;
                    return 0;
                }
                cnt_poi = ipoi_vec[0];
                //Break the for-loop
                touch_end = true;
            }
            
            //---------------------------------------------------------------------------
            //Check for overlapping
            //hout << "Check for overlapping ";
            if (!penetrating_model_flag || Check_penetration(geom_rve, nanotube_geo, cnts_points, global_coordinates, sectioned_domain, cnts_radius, new_cnt, n_subregions, cnt_rad, cutoffs.van_der_Waals_dist, MAX_ATTEMPTS, point_overlap_count, point_overlap_count_unique, cnt_poi)) {

                //---------------------------------------------------------------------------
                //If the overlapping model is used or if overlapping was successfully solved, then proceed as usual
                //---------------------------------------------------------------------------
                
                double temp_length = Effective_length_given_region(gvcub, new_cnt.back(),cnt_poi);
                //double temp_length = new_cnt.back().distance_to(cnt_poi);
                if (temp_length > 0.0)
                {
                    vol_sum += temp_length*step_vol_para;		//a accumulation on the volume
                    wt_sum += temp_length*step_wei_para;		//a accumulation on the weight
                }
                
                new_cnt.push_back(cnt_poi);							//store a new point
                new_cnt_size = (int)new_cnt.size()-1;				//calculate the size of new point
            } else {
                //---------------------------------------------------------------------------
                //If the overlapping was not solved, then delete the current CNT and generate a new one
                //---------------------------------------------------------------------------
                //hout << "Penetrating point could not be accommodated" <<endl;
                //Remove the volume and weight fractions corresponding to the new_cnt
                double temp_length;
                for(int j=0; j<(int)new_cnt.size()-1; j++) {
                    if(new_cnt[j+1].flag!=0) {
                        //Calculate segment length
                        temp_length = Effective_length_given_region(gvcub, new_cnt.back(),cnt_poi);
                        //new_cnt[j].distance_to(new_cnt[j+1]);
                        //Subtract corresponding volume
                        vol_sum -= temp_length*step_vol_para;
                        //Subtract corresponding weight
                        wt_sum -= temp_length*step_wei_para;
                    }
                }
                
                //Clear the new_cnt vector so that it is not added to the rest of CNTs
                new_cnt.clear();
                
                //Increase the count of rejected cnts
                cnt_reject_count++;
                
                //Break the for-loop
                touch_end = true;
            }
            //hout << "done" << endl;
            
            
            //---------------------------------------------------------------------------
            //Judge the new volume or weight
            if(nanotube_geo.criterion == "vol"&&vol_sum >= nanotube_geo.real_volume) break;		//Break out when the volume reaches the critical value
            else if(nanotube_geo.criterion == "wt"&&wt_sum >= nanotube_geo.real_weight) break;		//Break out when the weight reaches the critical value
            else if (touch_end) break; //Enforce break
        }
        
        //---------------------------------------------------------------------------
        //Store the CNT points
        if(new_cnt.size() >= 2)
        {
            //If the new_cnt vector has at least two points, then it can be added to the rest of the points
            vector<Point_3D> cnt_temp;
            //Variables needed for updating global_coordinates
            vector<int> empty;
            for(int i=0; i<(int)new_cnt.size(); i++)
            {
                cnt_temp.push_back(new_cnt[i]); //insert the first CNT point
                
                //Perform these operations when the non-overlapping model is used
                if (penetrating_model_flag) {
                    //Add global coordinate
                    global_coordinates.push_back(empty);
                    global_coordinates.back().push_back((int)cnts_points.size());
                    global_coordinates.back().push_back((int)cnt_temp.size()-1);
                    //Add point to an overlapping region in the vector sectioned_domain
                    Add_to_overlapping_regions(geom_rve, overlap_max_cutoff, new_cnt[i], (long int)global_coordinates.size()-1, n_subregions, sectioned_domain);
                }
                
                //Two cases (when using periodic boundary conditions)
                //a) Check if reached the end of the new CNT (new_cnt.size()-1)
                //b) Check if reached the end of a CNT but there are still points in new_cnt (flag=0)
                if(i==(int)new_cnt.size()-1||new_cnt[i+1].flag==0)
                {
                    if((int)cnt_temp.size()>1)									//if the size is equal to 1, that means the whole CNT only include one point, it doesn't creat a CNT segment
                    {
                        cnts_points.push_back(cnt_temp);				//to store the points of a CNT
                        cnts_radius.push_back(cnt_rad);					//to store the radius of a CNT
                    }
                    cnt_temp.clear();												//to clean the temporary cnt vector
                }
            }
        }
        
    }
    
    if(nanotube_geo.criterion == "wt") hout << "    The volume fraction of generated CNTs is about : " << vol_sum/geom_rve.volume << endl;
    
    hout << "There were " << point_overlap_count_unique << " overlapping points and ";
    hout << point_overlap_count << " overlaps, " << cnt_reject_count << " CNTs were rejected." << endl;
    
    return 1;
}
//Generate a network defined by points and connections
//Use the Mersenne Twister for the random number generation
int GenNetwork::Generate_network_threads_mt(const struct Geom_RVE &geom_rve, const struct Cluster_Geo &clust_geo, const struct Nanotube_Geo &nanotube_geo, const struct Cutoff_dist &cutoffs, vector<vector<Point_3D> > &cnts_points,  vector<double> &cnts_radius)const
{
    //Generate random seed in terms of local time
    //unsigned int time_seed = 1453384844;
    //unsigned int time_seed = (unsigned int)time(NULL);
    //hout << "Time seed "<<time_seed<<endl;
    //srand(time_seed);
    //srand((unsigned int)time(NULL));
    
    //---------------------------------------------------------------------------
    //Set up the Mersenne Twisters used for the different variables
    // Use random_device to generate a seed for Mersenne twister engine.
    std::random_device rd;
    // Use Mersenne twister engine to generate pseudo-random numbers.
    //Generate differnet engines for different variables
    std::mt19937 engine_x(rd());
    std::mt19937 engine_y(rd());
    std::mt19937 engine_z(rd());
    std::mt19937 engine_pha(rd());
    std::mt19937 engine_sita(rd());
    std::mt19937 engine_rand(rd());
    
    // "Filter" MT's output to generate double values, uniformly distributed on the closed interval [0, 1].
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    
    //---------------------------------------------------------------------------
    //Generate the data for nanotube clusters limited in ellipsoid surfaces (each ellipsoid is within the RVE and is separated with each other)
    if(clust_geo.vol_fra_criterion>0.0)
    {
        struct cuboid cub;										//Generate a cuboid for RVE
        cub.poi_min = geom_rve.origin;
        cub.len_x = geom_rve.len_x;
        cub.wid_y = geom_rve.wid_y;
        cub.hei_z = geom_rve.hei_z;
        cub.volume = cub.len_x*cub.wid_y*cub.hei_z;
        
        if(Get_ellip_clusters(cub, clust_geo)==0) return 0;
        //Generate a number of sperical clusters in regular arrangement
        //		if(Get_spherical_clusters_regular_arrangement(cub, clust_geo)==0) return 0;
    }
    
    //---------------------------------------------------------------------------
    double vol_sum = 0;  //the sum of volume of generated CNTs
    double wt_sum = 0;   //the sum of weight of generated CNTs
    int cnt_seed_count = 0; //to record the number of generated seeds of a CNT (If the growth of a CNT fails, but this seed will be used for the next growth of a CNT)
    int cnt_reject_count = 0; //to record the number of CNTs that were deleted due to penetration
    int point_overlap_count = 0; //to record the number of times that a point had to be relocated
    int point_overlap_count_unique = 0; //to record the number of points that were overlapping other points
    
    //---------------------------------------------------------------------------
    //Define the variable MAX_ATTEMPTS
    const int MAX_ATTEMPTS = 5;
    //Vectors for handling CNT penetration
    //global_coordinates[i][0] stores the CNT number of global point i
    //global_coordinates[i][1] stores the local point number of global point i
    vector<vector<int> > global_coordinates;
    //sectioned_domain[i] contains all the points in sub-region i.
    //Sub-region i is an overlapping subregion to check for penetrations
    vector<vector<long int> > sectioned_domain;
    //n_subregions[0] is the number of subregions along x
    //n_subregions[1] is the number of subregions along y
    //n_subregions[2] is the number of subregions along z
    vector<int> n_subregions;
    //Define cutoff for overlapping
    double overlap_max_cutoff = 2*nanotube_geo.rad_max + cutoffs.van_der_Waals_dist;
    //Initialize the vector sub-regions
    Initialize_subregions(geom_rve, n_subregions, sectioned_domain);
    //This flag will be used to skip overlapping functions
    //1 = non-penetrating model
    //0 = penetrating model
    int penetrating_model_flag = 1;
    
    //Generate cuboids that represent the extended domain and the composite domain
    //To calculate the effective portion (length) which falls into the given region (RVE)
    struct cuboid gvcub;					//generate a cuboid to represent the composite domain
    gvcub.poi_min = geom_rve.origin;
    gvcub.len_x = geom_rve.len_x;
    gvcub.wid_y = geom_rve.wid_y;
    gvcub.hei_z = geom_rve.hei_z;
    gvcub.volume = geom_rve.volume;
    struct cuboid excub;					//generate a cuboid to represent the extended domain
    excub.poi_min = geom_rve.ex_origin;
    excub.len_x = geom_rve.ex_len;
    excub.wid_y = geom_rve.ey_wid;
    excub.hei_z = geom_rve.ez_hei;
    excub.volume = excub.len_x*excub.wid_y*excub.hei_z;
    
    
    //---------------------------------------------------------------------------
    while((nanotube_geo.criterion == "vol"&&vol_sum < nanotube_geo.real_volume)||
          (nanotube_geo.criterion == "wt"&&wt_sum < nanotube_geo.real_weight))
    {
        //---------------------------------------------------------------------------
        //Define a vector for a new nanotube
        vector<Point_3D> new_cnt;
        int new_cnt_size = (int)new_cnt.size();
        
        //---------------------------------------------------------------------------
        //Randomly generate a length of a CNT
        double cnt_length;
        if(Get_random_value_mt(nanotube_geo.len_distrib_type, engine_rand, dist, nanotube_geo.len_min, nanotube_geo.len_max, cnt_length)==0) return 0;
        //Calculate the total number of growth step for a CNT
        int step_num = (int)(cnt_length/nanotube_geo.step_length) + 1;
        
        //---------------------------------------------------------------------------
        //Randomly generate a radius of a CNT
        double cnt_rad;
        if(Get_random_value_mt(nanotube_geo.rad_distrib_type, engine_rand, dist, nanotube_geo.rad_min, nanotube_geo.rad_max, cnt_rad)==0) return 0;
        
        //---------------------------------------------------------------------------
        //Randomly generate a direction in the spherical coordinates as the intital direction of CNT segments
        double cnt_sita, cnt_pha;
        if(Get_uniform_direction_mt(nanotube_geo, cnt_sita, cnt_pha, engine_sita, engine_pha, dist)==0) return 0;
        
        MathMatrix multiplier(3,3);
        multiplier = Get_transformation_matrix(cnt_sita, cnt_pha);
        
        //---------------------------------------------------------------------------
        //The increased volume of each segement (growth step) of nanotube (Here the overlapping volume is ignored)
        const double step_vol_para = PI*cnt_rad*cnt_rad;
        //---------------------------------------------------------------------------
        //The increased weight of each segement (growth step) of nanotube (If the different radii of nanotube are considered, the linear_density may be different in every nanotube)
        const double step_wei_para = nanotube_geo.linear_density;
        
        //---------------------------------------------------------------------------
        //Randomly generate a seed (initial point) of a CNT in the extended RVE	(Comments: the seed generation is after the radius generation for the non-overlapping nanotubes generation)
        
        Point_3D cnt_poi;
        if(Get_seed_point_mt(excub, cnt_poi, engine_x, engine_y, engine_z, dist)==0) return 0;
        
        //Check overlapping of the intial point
        int counter = 1;
        while (penetrating_model_flag && !Check_penetration(geom_rve, nanotube_geo, cnts_points, global_coordinates, sectioned_domain, cnts_radius, new_cnt, n_subregions, cnt_rad, cutoffs.van_der_Waals_dist, MAX_ATTEMPTS, point_overlap_count, point_overlap_count_unique, cnt_poi)) {
            if(Get_seed_point_mt(excub, cnt_poi, engine_x, engine_y, engine_z, dist)==0) return 0;
            cnt_seed_count++;					//record the number of seed generations
            //hout << "Seed deleted" << endl;
            if (counter == MAX_ATTEMPTS) {
                hout << "Too many attempts to resolve overlapping of an intial CNT point (" << counter << " attempts). ";
                hout << cnt_poi.x << ' ' << cnt_poi.y << ' ' << cnt_poi.z << endl;
                return 0;
            }
            counter ++;
        }//*/
        
        new_cnt.push_back(cnt_poi);	//store this seed point in the vector for a new nanotube
        
        //---------------------------------------------------------------------------
        cnt_seed_count++;					//record the number of seed generations
        int max_seed = 1E9;
        if(cnt_seed_count>max_seed)
        {
            hout << "The number of seed genrations is lager than "<<max_seed<<", but the nanotube generation still fails to acheive the demanded volume fraction." << endl;
            return 0;
        }
        //---------------------------------------------------------------------------
        //
        
        //---------------------------------------------------------------------------
        //The growth process of nanotube
        int ellip_num = -1; //For recording the serial number of ellipsoid cluster which a nanotube penetrates out. It is no use if the cluster generation is not considered.
        for(int i=0; i<step_num; i++)
        {
            //Randomly generate a direction in the spherical coordinates
            //To have the positive Z-axis to be a central axis
            //Then, the radial angle, sita, obeys a normal distribution (sita \in fabs[(-omega,+omega)]) and the zonal angle, pha, obeys a uniform distribution (pha \in (0,2PI))
            if(Get_normal_direction_mt(nanotube_geo.angle_max, cnt_sita, cnt_pha, engine_sita, engine_pha, dist)==0) return 0;
            
            //To calculate the new multiplier for transformation of coordinates
            multiplier = multiplier*Get_transformation_matrix(cnt_sita, cnt_pha);
            
            //To calculate the coordinates of the new CNT point (transformation of coordinates)
            cnt_poi = cnt_poi + Get_new_point(multiplier, nanotube_geo.step_length);
            cnt_poi.flag = 1;							//1 means that point is not the intial point
            
            //---------------------------------------------------------------------------
            //If a CNT penetrates the ellipsoidal surface of a cluster from inside, the growth of this CNT will be headed back in the cluster in a probability p, (0<=p<=1).
            if(clust_geo.ellips.size()>0)
            {
                return 0;
            }
            
            //---------------------------------------------------------------------------
            //If the new CNT point grows out of the RVE, the intersecting point at the surfaces of RVE will be calculated.
            //The new segment will be cut and (the outside part will be translated into the RVE for the periodical boundary condition)
            bool touch_end = false;
            if(Judge_RVE_including_point(excub, cnt_poi)==0)
            {
                //Calculate all intersection points between the new segment and surfaces of RVE
                //(using a parametric equatio:  the parameter 0<t<1, and sort all intersection points from the smaller t to the greater t)
                vector<Point_3D> ipoi_vec;  //a vector for intersection points
                if(Get_intersecting_point_RVE_surface(excub, new_cnt.back(), cnt_poi, ipoi_vec)==0) {
                    hout << "Error in Generate_network_threads"<<endl;
                    return 0;
                }
                cnt_poi = ipoi_vec[0];
                //Break the for-loop
                touch_end = true;
            }
            
            //---------------------------------------------------------------------------
            //Check for overlapping
            //hout << "Check for overlapping ";
            if (!penetrating_model_flag || Check_penetration(geom_rve, nanotube_geo, cnts_points, global_coordinates, sectioned_domain, cnts_radius, new_cnt, n_subregions, cnt_rad, cutoffs.van_der_Waals_dist, MAX_ATTEMPTS, point_overlap_count, point_overlap_count_unique, cnt_poi)) {
                
                //---------------------------------------------------------------------------
                //If the overlapping model is used or if overlapping was successfully solved, then proceed as usual
                //---------------------------------------------------------------------------
                
                double temp_length = Effective_length_given_region(gvcub, new_cnt.back(),cnt_poi);
                //double temp_length = new_cnt.back().distance_to(cnt_poi);
                if (temp_length > 0.0)
                {
                    vol_sum += temp_length*step_vol_para;		//a accumulation on the volume
                    wt_sum += temp_length*step_wei_para;		//a accumulation on the weight
                }
                
                new_cnt.push_back(cnt_poi);							//store a new point
                new_cnt_size = (int)new_cnt.size()-1;				//calculate the size of new point
            } else {
                //---------------------------------------------------------------------------
                //If the overlapping was not solved, then delete the current CNT and generate a new one
                //---------------------------------------------------------------------------
                //hout << "Penetrating point could not be accommodated" <<endl;
                //Remove the volume and weight fractions corresponding to the new_cnt
                double temp_length;
                for(int j=0; j<(int)new_cnt.size()-1; j++) {
                    if(new_cnt[j+1].flag!=0) {
                        //Calculate segment length
                        temp_length = Effective_length_given_region(gvcub, new_cnt.back(),cnt_poi);
                        //new_cnt[j].distance_to(new_cnt[j+1]);
                        //Subtract corresponding volume
                        vol_sum -= temp_length*step_vol_para;
                        //Subtract corresponding weight
                        wt_sum -= temp_length*step_wei_para;
                    }
                }
                
                //Clear the new_cnt vector so that it is not added to the rest of CNTs
                new_cnt.clear();
                
                //Increase the count of rejected cnts
                cnt_reject_count++;
                
                //Break the for-loop
                touch_end = true;
            }
            //hout << "done" << endl;
            
            
            //---------------------------------------------------------------------------
            //Judge the new volume or weight
            if(nanotube_geo.criterion == "vol"&&vol_sum >= nanotube_geo.real_volume) break;		//Break out when the volume reaches the critical value
            else if(nanotube_geo.criterion == "wt"&&wt_sum >= nanotube_geo.real_weight) break;		//Break out when the weight reaches the critical value
            else if (touch_end) break; //Enforce break
        }
        
        //---------------------------------------------------------------------------
        //Store the CNT points
        if(new_cnt.size() >= 2)
        {
            //If the new_cnt vector has at least two points, then it can be added to the rest of the points
            vector<Point_3D> cnt_temp;
            //Variables needed for updating global_coordinates
            vector<int> empty;
            for(int i=0; i<(int)new_cnt.size(); i++)
            {
                cnt_temp.push_back(new_cnt[i]); //insert the first CNT point
                
                //Perform these operations when the non-overlapping model is used
                if (penetrating_model_flag) {
                    //Add global coordinate
                    global_coordinates.push_back(empty);
                    global_coordinates.back().push_back((int)cnts_points.size());
                    global_coordinates.back().push_back((int)cnt_temp.size()-1);
                    //Add point to an overlapping region in the vector sectioned_domain
                    Add_to_overlapping_regions(geom_rve, overlap_max_cutoff, new_cnt[i], (long int)global_coordinates.size()-1, n_subregions, sectioned_domain);
                }
                
                //Two cases (when using periodic boundary conditions)
                //a) Check if reached the end of the new CNT (new_cnt.size()-1)
                //b) Check if reached the end of a CNT but there are still points in new_cnt (flag=0)
                if(i==(int)new_cnt.size()-1||new_cnt[i+1].flag==0)
                {
                    if((int)cnt_temp.size()>1)									//if the size is equal to 1, that means the whole CNT only include one point, it doesn't creat a CNT segment
                    {
                        cnts_points.push_back(cnt_temp);				//to store the points of a CNT
                        cnts_radius.push_back(cnt_rad);					//to store the radius of a CNT
                    }
                    cnt_temp.clear();												//to clean the temporary cnt vector
                }
            }
        }
        
    }
    
    if(nanotube_geo.criterion == "wt") hout << "    The volume fraction of generated CNTs is about : " << vol_sum/geom_rve.volume << endl;
    
    hout << "There were " << point_overlap_count_unique << " overlapping points and ";
    hout << point_overlap_count << " overlaps, " << cnt_reject_count << " CNTs were rejected." << endl;
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Generate a random value through a probability distribution function
int GenNetwork::Get_random_value_mt(const string &dist_type, mt19937 &engine, uniform_real_distribution<double> &dist, const double &min, const double &max, double &value)const
{
    //Check if limits are correctly defined
    if(min>max) { hout << "Error, the minimum value is larger than the maximum value (Get_random_value)!" << endl; return 0; }
    
    //Check if the interval has 0 length
    if (max == min) {
        //In this case, value is either of the limits
        //To be consistent with the formulation below, value is set equal to min
        value = min;
        return 1;
    }
    
    if(dist_type=="uniform")	//uniform distribution
    {
        value = (max-min)*dist(engine) + min;
    }
    else if(dist_type=="normal")	//normal distribution
    {
        double sum=0;
        for(int i=0; i<12; i++)
        {
            sum = sum + dist(engine);
        }
        value = (max-min)*sum/12.0 + min;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Randomly generate a seed (intial point) of a CNT in the RVE
int GenNetwork::Get_seed_point_mt(const struct cuboid &cub, Point_3D &point, mt19937 &engine_x, mt19937 &engine_y, mt19937 &engine_z, uniform_real_distribution<double> &dist)const
{
    
    point.x = cub.poi_min.x + cub.len_x*dist(engine_x);
    
    point.y = cub.poi_min.y + cub.wid_y*dist(engine_y);
    
    point.z = cub.poi_min.z + cub.hei_z*dist(engine_z);
    
    
    point.flag = 0; //0 denotes this point is the initial point of a CNT
    
    return 1;
}
//---------------------------------------------------------------------------
//Randomly generate a direction in the spherical coordinates as the initial direction of CNT segments
int GenNetwork::Get_uniform_direction_mt(const struct Nanotube_Geo &nanotube_geo, double &cnt_sita, double &cnt_pha, mt19937 &engine_sita, mt19937 &engine_pha, uniform_real_distribution<double> &dist)const
{
    if(nanotube_geo.dir_distrib_type=="random")
    {
        //sita is chosen in [0, PI] with uniform distribution
        cnt_sita = PI*dist(engine_sita);
        
        //pha is chosen in [0, 2PI] with uniform distribution
        cnt_pha = 2.0*PI*dist(engine_pha);//*/
        
    }
    else if(nanotube_geo.dir_distrib_type=="specific")
    {
        //Use the probability of a random number to be even
        if( (engine_sita()+engine_pha())%2==0)
        {
            cnt_sita = nanotube_geo.ini_sita;		//"positive" direction
            cnt_pha = nanotube_geo.ini_pha;		//"positive" direction
        }
        else
        {
            cnt_sita = PI - nanotube_geo.ini_sita;		//"negative" (opposite) direction
            cnt_pha = PI + nanotube_geo.ini_pha;	//"negative" (opposite) direction
        }
    }
    
    return 1;
}
//---------------------------------------------------------------------------
int GenNetwork::Get_normal_direction_mt(const double &omega, double &cnt_sita, double &cnt_pha, mt19937 &engine_sita, mt19937 &engine_pha, uniform_real_distribution<double> &dist)const
{
    
    //sita centers around 0 and obeys a normal distribution in (-omega, +omega)
    double sum=0;
    for(int i=0; i<12; i++)
    {
        sum = sum + dist(engine_sita);
    }
    cnt_sita = fabs(omega*(sum/6 - 1));
    
    //pha satisfies a uniform distribution in (0, 2PI)
    cnt_pha = 2.0*PI*dist(engine_pha);//*/
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//This functions initializes the vectors n_subregions and sectioned_domain
//
//The n_subregions vector is defined to avoid calculating the number of sub-regions for every point when the functions
//Default_region and Add_to_subregion are called. Thus, saving computational time.
//n_subregions[0] is the number of subregions along x
//n_subregions[1] is the number of subregions along y
//n_subregions[2] is the number of subregions along z
//
//The vector sectioned_domain contains the sub-regions to look for overlapping
//It is initialized with the number of sub-regions in the sample
void GenNetwork::Initialize_subregions(const struct Geom_RVE &geom_rve, vector<int> &nsubregions, vector<vector<long int> > &sectioned_domain)const
{
    //Initialize nsubregions
    
    //variable to store the number of subregions
    int s;
    //Number of subregions along x
    s = (int)(geom_rve.len_x/geom_rve.gs_minx);
    //Add the number of sub-regions to the vector
    nsubregions.push_back(s);
    //Number of subregions along y
    s = (int)(geom_rve.wid_y/geom_rve.gs_miny);
    //Add the number of sub-regions to the vector
    nsubregions.push_back(s);
    //Number of subregions along z
    s = (int)(geom_rve.hei_z/geom_rve.gs_minz);
    //Add the number of sub-regions to the vector
    nsubregions.push_back(s);
    
    //Initialize sectioned_domain
    vector<long int> empty;
    sectioned_domain.assign(nsubregions[0]*nsubregions[1]*nsubregions[2], empty);
}
//Check if the current CNT is penetrating another CNT, i.e. is the new point is overlapping other point
//1: a) No penetration
//   b) No need to check for penetration (point is in boundary layer or there are no other points in the same sub-region)
//   c) There was penetration but it was succesfully resolved
//0: There was penetration but could not be resolved
int GenNetwork::Check_penetration(const struct Geom_RVE &geom_rve, const struct Nanotube_Geo &nanotube_geo, const vector<vector<Point_3D> > &cnts, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain, const vector<double> &radii, const vector<Point_3D> &cnt_new, const vector<int> &n_subregions, const double &cnt_rad, const double &d_vdw, const int &MAX_ATTEMPTS, int &point_overlap_count, int &point_overlap_count_unique, Point_3D &point)const
{
    //Get the sub-region the point belongs to
    int subregion = Get_subregion(geom_rve, n_subregions, point);
    
    //If the sub-region is -1, then the point is in te boundary layer, so there is no need to check penetration
    if (subregion == -1) {
        return 1;
    }
    //This vector will store the coordintes of the points that the input "point" is penetrating
    vector<vector<int> > affected_points;
    //This vector stores the distance at which the two points should be
    vector<double> cutoffs_p;
    //This vector stores the distance at which the two points actually are
    vector<double> distances;
    //Need to keep the count of the number of attempts outside the for-loop
    int attempts;
    //I move the point up to max_attempts times. If there is still penetration then I delete it
    for (attempts = 0; attempts <= MAX_ATTEMPTS; attempts++) {
        
        //Check if there are any penetrations in the corresponding sub-region
        Get_penetrating_points(cnts, global_coordinates, sectioned_domain[subregion], radii, cnt_rad, d_vdw, point, affected_points, cutoffs_p, distances);
        
        //--------------------------------------------------------------------------------------------
        //Check if there are any penetrating points
        if (affected_points.size()) {
            //Update the counter of overlaps
            point_overlap_count++;
            
            //Update the counter of overlapping points only when an overlapping point was found the first time
            //i.e. attempts = 0
            if (!attempts) {
                point_overlap_count_unique++;
            }
            
            //If this is the last iteration and there are still affected points then the point could not be accommodated
            if (attempts == MAX_ATTEMPTS) {
                //hout << "Deleted CNT number " << cnts.size() << " of size " << cnt_new.size();
                //hout << " (reached maximum number of attempts for relocation)" << endl;//*/
                return 0;
            }
            
            //Find the new point
            //hout << "Point " << global_coordinates.size()-1+cnt_new.size() << " in CNT " << cnts.size() << " is overlapping." <<endl;
            //hout << "Moved a point from initial position (" << point.x << ", " << point.y << ", " << point.z << ")." << endl;
            Move_point(geom_rve, nanotube_geo, cnts, cnt_new, point, cutoffs_p, distances, affected_points);
            //hout << "Moved a point to final position (" << point.x << ", " << point.y << ", " << point.z << ")." << endl;
            
            //Check that the new point is within the permited orientation repect to the previous segment
            if (!Check_segment_orientation(point, cnt_new)) {
                //hout << "Deleted CNT number " << cnts.size() << " of size " << cnt_new.size();
                //hout << " (the point is not in a valid orientation)" << endl;//*/
                //When not in a valid position it cannot be moved again so a new CNT is needed
                return 0;
            }
            
            //Need to update point sub-region as it could be relocated to a new sub-region
            subregion = Get_subregion(geom_rve, n_subregions, point);
            //Check if after moving the point it is now in the boundary layer
            if (subregion == -1) {
                //If the point is now in the boundary layer, terminate the function
                //there is no need to continue checking
                return 1;
            }
            
            //Need to clear the vectors affected_points, contact_coordinates and temporal_contacts so they are used again with the new point
            affected_points.clear();
            cutoffs_p.clear();
            distances.clear();
            
        } else {
            //if the size of affected_points is zero, then terminate the function
            return 1;
        }
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//This function returns the subregion a point belongs to
int GenNetwork::Get_subregion(const struct Geom_RVE &geom_rve, const vector<int> &n_subregions, const Point_3D &point)const
{
    if (Judge_RVE_including_point(geom_rve, point)) {
        //These variables will give me the region cordinates of the region that a point belongs to
        int a, b, c;
        //Calculate the region-coordinates
        a = (int)((point.x-geom_rve.origin.x)/geom_rve.gs_minx);
        //Limit the value of a as it has to go from 0 to n_subregions[0]-1
        if (a == n_subregions[0]) a--;
        b = (int)((point.y-geom_rve.origin.y)/geom_rve.gs_miny);
        //Limit the value of b as it has to go from 0 to n_subregions[1]-1
        if (b == n_subregions[1]) b--;
        c = (int)((point.z-geom_rve.origin.z)/geom_rve.gs_minz);
        //Limit the value of c as it has to go from 0 to n_subregions[2]-1
        if (c == n_subregions[2]) c--;
        
        return (a + (b*n_subregions[0]) + (c*n_subregions[0]*n_subregions[1]));
    } else {
        //If the point is in the boundary layer, then there is no need to calculate is sub-region
        return -1;
    }
}
//This functions iterates over a sub-region and determines if there are any penetrating points
//If there are penetrating points, they are stored in the vector affected_points
void GenNetwork::Get_penetrating_points(const vector<vector<Point_3D> > &cnts, const vector<vector<int> > &global_coordinates, const vector<long int> &subregion_vec, const vector<double> &radii, const double &cnt_rad, const double &d_vdw, Point_3D &point, vector<vector<int> > &affected_points, vector<double> &cutoffs_p, vector<double> &distances)const
{
    //They are just intermediate variables and I only use them to make the code more readable
    long int coord2;
    int P2, CNT2;
    //cutoff_p is used for the cutoff between two points (of different CNTs), distance is the actual distance
    //between those two points
    double cutoff_p, distance;
    
    //hout << "Check0 " << subregion_vec.size() << ' ';
    for (long int i = 0; i < (long int)subregion_vec.size(); i++) {
        //hout << "Check1 " ;
        coord2 = subregion_vec[i];
        //hout << "Check2 ";
        CNT2 = global_coordinates[coord2][0];
        P2 = global_coordinates[coord2][1];
        //hout << "Check4 ";
        cutoff_p = cnt_rad + radii[CNT2] + d_vdw;
        //Check is the second point is in the cube of size 2cutoff_p and centered in P1
        //This is easier and faster to check than calculating the distance from poin to point every time
        if ( (cnts[CNT2][P2].x<point.x+cutoff_p)&&(cnts[CNT2][P2].x>point.x-cutoff_p)&&(cnts[CNT2][P2].y<point.y+cutoff_p)&&(cnts[CNT2][P2].y>point.y-cutoff_p)&&(cnts[CNT2][P2].z<point.z+cutoff_p)&&(cnts[CNT2][P2].z>point.z-cutoff_p) ) {
            distance = point.distance_to(cnts[CNT2][P2]);
            //If it is inside the cube, then it is worth to take the time to calculate the distance from point ot point
            if (distance < cutoff_p) {
                affected_points.push_back(global_coordinates[coord2]);
                cutoffs_p.push_back(cutoff_p);
                distances.push_back(distance);
                /*/hout << "CNT=" << CNT1 << " Point=" << P1 << " r1=" << cnt_rad;
                hout << " Penetrating points="<< affected_points.size();
                hout << " CNT2=" << CNT2 << " P2=" << P2 << " r2=" << radii[CNT2] << " (" << cnts[CNT2][P2].x << ", " << cnts[CNT2][P2].y << ", " << cnts[CNT2][P2].z << ") ";
                hout << endl;//*/
            }
        }
        //hout << "Check5 ";
    }
    //hout << "Check6 " << endl;
}
//This function moves a point according to the number of points it is overlapping
void GenNetwork::Move_point(const struct Geom_RVE &geom_rve, const struct Nanotube_Geo &nanotube_geo, const vector<vector<Point_3D> > &cnts, const vector<Point_3D> &cnt_new, Point_3D &point, vector<double> &cutoffs, vector<double> &distances, vector<vector<int> > &affected_points)const
{
    //The number of overlapings will determine how the new point is moved
    //However, first I need to eliminate invalid-points, which are the points of perfect overlapping, i.e.
    //points in the exact same location
    int initial_overlappings = (int)affected_points.size();
    int overlappings = Check_points_in_same_position(cutoffs, distances, affected_points);
    if (overlappings == 1){
        //Handle the case with only one overlapping point
        //hout << "ovelappings == 1"<<endl;
        One_overlapping_point(cnts, cutoffs, distances, affected_points, point);
    } else if (overlappings == 2){
        //Handle the case with two overlapping points
        //hout << "ovelappings == 2"<<endl;
        Two_overlapping_points(cnts, cutoffs, affected_points, point);
    } else if (overlappings >= 3) {
        //Handle the case with three overlapping points
        //This actually finds the two closest point and calls the function
        //that handles the case with two overlapping points
        //hout << "ovelappings == 3"<<endl;
        Three_or_more_overlapping_points(cnts, cutoffs, distances, affected_points, point);
    } else {//if (!ovelappings) {
        //If after cheking for points in the same position there are no ovelappings,
        //then all points are overlapping are in the same location
        //hout << "ovelappings == 0"<<endl;
        Overlapping_points_same_position(geom_rve, nanotube_geo, cnt_new, point);
        hout << "There were " << initial_overlappings - overlappings << " points overlapping in exactly the same location: ";
        hout << "P("<<point.x<<", "<<point.y<<", "<<point.z<<")."<<endl;
    }
}
//This point checks if two points are actually in the same position
//This used to happen because of a bug in the code. It seems now like a remote possibility
//so I'll keep it just in case
int GenNetwork::Check_points_in_same_position(vector<double> &cutoffs, vector<double> &distances, vector<vector<int> > &affected_points)const
{
    //hout << "Initial overlaps = " << distances.size();
    //check if any distance is less than the variable Zero
    for (int i = (int)distances.size()-1; i >= 0 ; i--) {
        if (distances[i] < Zero){
            //Delete points that are in the exact same location
            distances.erase(distances.begin()+i);
            cutoffs.erase(cutoffs.begin()+i);
            affected_points.erase(affected_points.begin()+i);
        }
    }
    //hout << ", valid overlaps = " << distances.size() <<endl;
    //This is the number of valid overlapping points
    return (int)distances.size();
}
//Move the point when all points are in the same location
//After solving the bug that caused multiple points to be in the same location,
//probably this function is not needed. I leave it here just in case
void GenNetwork::Overlapping_points_same_position(const struct Geom_RVE &geom_rve, const struct Nanotube_Geo &nanotube_geo, const vector<Point_3D> &cnt_new, Point_3D &point)const
{
    //The new point will be moved depending on wheter it is the first point, second point or other point after the second
    if (!cnt_new.size()) {
        //Point is the first point
        //Move the point to a random location
        point.x = ((double)rand()/RAND_MAX)*geom_rve.ex_len + geom_rve.ex_origin.x;
        point.y = ((double)rand()/RAND_MAX)*geom_rve.ey_wid + geom_rve.ex_origin.y;
        point.z = ((double)rand()/RAND_MAX)*geom_rve.ez_hei + geom_rve.ex_origin.z;
    } else if (cnt_new.size() == 1) {
        //point is the second point
        //If there is only one point, just move the new point to a new direction
        
        //Generate a rotation matrix
        double phi = ((double)rand()/RAND_MAX)*2*PI;
        //For the second point, the limitation on the angles is not important since there are no other
        //segments to compare with
        double theta = ((double)rand()/RAND_MAX)*2*PI - PI;
        MathMatrix rotation(3,3);
        rotation = Get_transformation_matrix(theta, phi);
        
        //Calculate new point
        //The operation point = cnt_new.front() + Get_new_point(rotation, nanotube_geo.len_max)
        //has to be done in two steps because I get an error. The types of Point_3D are different
        //one is const and the other is not
        point = cnt_new.front();
        point = point + Get_new_point(rotation, nanotube_geo.len_max);
    } else {
        //point is the third point of higher
        //If there are 2 or more points, calculate the z_i unit vector and then move the new point
        //to a random direction
        
        //Generate a rotation matrix
        double phi = ((double)rand()/RAND_MAX)*2*PI;
        double theta = ((double)rand()/RAND_MAX)*PI - PI/2;
        MathMatrix rotation(3,3);
        rotation = Get_transformation_matrix(theta, phi);
        
        //Calculate the z_i unit vector
        //The operation Point_3D z_i = cnt_new.back() - cnt_new[cnt_new.size()-2];
        //has to be done in two steps because I get an error.
        Point_3D z_i = cnt_new.back();
        z_i = z_i - cnt_new[cnt_new.size()-2];
        z_i = z_i/(z_i.distance_to(z_i)); //unit vector
        
        //Calculate new point
        //temporary matrix to store a matrix vector multiplication
        MathMatrix vec(3,1);
        vec.element[0][0] = z_i.x;
        vec.element[1][0] = z_i.y;
        vec.element[2][0] = z_i.z;
        //Rotate unit vector, and multiply by the magnitude of the segment length
        vec = (rotation*vec)*nanotube_geo.step_length;
        //Add coordintes to last point in cnt_new to create new point
        point = cnt_new.back();
        point.x = point.x + vec.element[0][0];
        point.y = point.y + vec.element[1][0];
        point.z = point.z + vec.element[2][0];
        
    }
}
//This function finds the new location for an overlapping point when it overlaps only one point
void GenNetwork::One_overlapping_point(const vector<vector<Point_3D> > &cnts, const vector<double> &cutoffs, const vector<double> &distances, const vector<vector<int> > &affected_points, Point_3D &point)const
{
    //When there is overlapping with one point only, then this is the simplest and easiest case
    //Just move the point in the direction form P_old to P_new a distance cutoff from P_old
    
    //Get the penetrating point. Its coordinates are in the first (and only) element of affected_points
    Point_3D P = cnts[affected_points[0][0]][affected_points[0][1]];
    //Calculate direction unit vector. distances[0] already has the lenght of the vector form point to P
    Point_3D direction = (point - P)/distances[0];
    //So the new point is P_old + d*n. d is the cutoff and n the unit vector
    point = P + direction*(cutoffs[0]+Zero); //The Zero is to avoid machine precision errors. Without it, when comparing
    //the new point with the other points in the same region, the program was judging them to be below the cutoff
    //for the van der Waals distance. Even though they were in the limit. After adding this Zero that issue
    //was eliminated
}
//This function finds the new location for an overlapping point when it overlaps two points
void GenNetwork::Two_overlapping_points(const vector<vector<Point_3D> > &cnts, const vector<double> &cutoffs, const vector<vector<int> > &affected_points, Point_3D &point)const
{
    //Point variables
    Point_3D P, Q, R, P1, P2;
    //distance variables
    double a, b, c, d;
    
    //Get the penetrating points.
    P1 = cnts[affected_points[0][0]][affected_points[0][1]];
    P2 = cnts[affected_points[1][0]][affected_points[1][1]];
    //Calculate P vector
    P = P2 - P1;
    //Calculate Q vector
    Q = point - P1;
    //Calculate normal vector PxQ
    R = (P.cross(Q)).cross(P);
    //Sides of the triangle
    a = cutoffs[0];
    b = cutoffs[1];
    c = P1.distance_to(P2);
    //Distance from P1 to M
    d = (b*b - a*a - c*c)/(-2*c);
    //Make P a unit vector
    P = P/sqrt(P.dot(P));
    //Make R a unit vector
    R = R/sqrt(R.dot(R));
    //Calculate new position
    point = P1 + P*(d + Zero) + R*(sqrt(a*a - d*d)+Zero);//The Zero is to avoid machine precision errors. Without it, when comparing
    //the new point with the other points in the same region, the program was judging them to be below the cutoff
    //for the van der Waals distance. Even though they were in the limit. After adding this Zero that issue
    //was eliminated
}
//This function finds the two closest points and calls the function that moves a point that overlaps two other points
//When a point overlaps three or more points, it becomes too difficult to find the new location
//Hence, this function that finds the two closest points
//The two closest point are chosen since those would be the more critical ones
void GenNetwork::Three_or_more_overlapping_points(const vector<vector<Point_3D> > &cnts, const vector<double> &cutoffs, const vector<double> &distances, const vector<vector<int> > &affected_points, Point_3D &point)const
{
    //Use the distances vector to find the two closest points
    //i1 and i2 will be the indices of the closest affected_points, they will be initialized with the first two
    //i1 will have the index of the closest point, while i2 will have the index of the second closest point
    int i1, i2;
    //d1 and d2 will be the distances to the closest affected_points, they will be initialized with the first two
    //d1 will have the distance to the closest point, while d2 will have the distance to the second closest point
    double d1, d2;
    //Sort the first two distances
    if (distances[0] < distances[1]) {
        i1 = 0;
        d1 = distances[0];
        i2 = 1;
        d2 = distances[1];
    } else {
        i1 = 1;
        d1 = distances[1];
        i2 = 0;
        d2 = distances[0];
    }
    
    //Once the closest points have been initialized, then scan the rest of the affected points
    for (int i = 2; i < (int)distances.size(); i++) {
        //First check if the distances[i] is smaller than d2
        if (distances[i] < d2) {
            //In this case, distances[i] is one of the two closest point
            //Now I need to check against d1 in case distances[i] is the closest point now
            if (distances[i] < d1) {
                //distances[i] is the closest point
                //update the distances
                d2 = d1;
                d1 = distances[i];
                //update the indices
                i2 = i1;
                i1 = i;
            } else {
                //distances[i] is the second closest point
                //update the distances
                d2 = distances[i];
                i2 = i;
            }
        }
    }
    
    //Generate the necessary vectors for the case of two overlapping points
    vector<vector<int> > two_affected_points;
    two_affected_points.push_back(affected_points[i1]);
    two_affected_points.push_back(affected_points[i2]);
    vector<double> two_cutoffs;
    two_cutoffs.push_back(cutoffs[i1]);
    two_cutoffs.push_back(cutoffs[i2]);
    
    //Now, call the function that moves a point that overlaps two points
    Two_overlapping_points(cnts, two_cutoffs, two_affected_points, point);
}
//This function checks that "point" is within the bounds of the segment orientation.
//The criterion is just checking the point is not more than pi/2 respect with the previous
//segment. In the limiting case we have a straight triangle. So I calculate the hypotenuse.
//I also measure the distance between "point" and the second before that.
//If the distance  between points is less than the hypotenuse, then it has an
//incorrect orientation
int GenNetwork::Check_segment_orientation(const Point_3D &point, const vector<Point_3D> &cnt_new)const
{
    //If at least two points have already been generated, then check if the new point has a valid orientation
    if (cnt_new.size()>=2) {
        int last = (int)cnt_new.size()-1;
        //I don't need the square root of the function distance_to, so to save time I perform these operations
        //distance_to(cnt_new[last]);
        double d1 = (point.x - cnt_new[last].x)*(point.x - cnt_new[last].x) + (point.y - cnt_new[last].y)*(point.y - cnt_new[last].y) + (point.z - cnt_new[last].z)*(point.z - cnt_new[last].z);
        //cnt_new[last].distance_to(cnt_new[last-1]);
        double d2 = (cnt_new[last].x - cnt_new[last-1].x)*(cnt_new[last].x - cnt_new[last-1].x) + (cnt_new[last].y - cnt_new[last-1].y)*(cnt_new[last].y - cnt_new[last-1].y) + (cnt_new[last].z - cnt_new[last-1].z)*(cnt_new[last].z - cnt_new[last-1].z);
        //hypotenuse squared
        double hypotenuse = d1 + d2;
        
        //When hypotenuse is equal to segment, there could be a machine precision error.
        //In this case, they are equal but judged to be different. Then I use the Zero to reduce the
        //incidence of that issue. (Just as I did in the function Move_point)
        //I subtract the zero so the limiting case is included when segment => hypothenuse
        //
        //This actually calculates the squared of cnt_new[last-1].distance_to(point)
        double segment = (point.x - cnt_new[last-1].x)*(point.x - cnt_new[last-1].x) + (point.y - cnt_new[last-1].y)*(point.y - cnt_new[last-1].y) + (point.z - cnt_new[last-1].z)*(point.z - cnt_new[last-1].z) - Zero;
        
        //Since both segment and hypotenuse are lentghs both are greater or equal to zero
        //Hence, comparing them is the same as comparing their squares
        //So instead of calculating a square root, I calculate segement*segment and compare with
        //hypotenuse squared
        if (segment<hypotenuse) {
            //The point is not in a valid position
            return 0;
        } else {
            //The point is in a valid position
            return 1;
        }
    } else {
        //If point is the first or second point, its orientation does not matter
        return 1;
    }
}
//This function adds a point to a region so penetration can be checked
void GenNetwork::Add_to_overlapping_regions(const struct Geom_RVE &geom_rve, double overlap_max_cutoff, Point_3D point, long int global_num, const vector<int> &n_subregions, vector<vector<long int> > &sectioned_domain)const
{
    //A point is added only if it is in the composite domain
    //If the point is in the boundary layer, overlapping is not important
    if (Judge_RVE_including_point(geom_rve, point)) {
        //Save coordinates of the point
        double x = point.x;
        double y = point.y;
        double z = point.z;
        
        //These variables will give me the region cordinates of the region that a point belongs to
        int a, b, c;
        //Calculate the region-coordinates
        a = (int)((x-geom_rve.origin.x)/geom_rve.gs_minx);
        //Limit the value of a as it has to go from 0 to n_subregions[0]-1
        if (a == n_subregions[0]) a--;
        b = (int)((y-geom_rve.origin.y)/geom_rve.gs_miny);
        //Limit the value of b as it has to go from 0 to n_subregions[1]-1
        if (b == n_subregions[1]) b--;
        c = (int)((z-geom_rve.origin.z)/geom_rve.gs_minz);
        //Limit the value of c as it has to go from 0 to n_subregions[2]-1
        if (c == n_subregions[2]) c--;
        
        //These variables are the coordinates of the lower corner of the RVE that defines its geometry
        double xmin = geom_rve.origin.x;
        double ymin = geom_rve.origin.y;
        double zmin = geom_rve.origin.z;
        
        //Coordinates of non-overlaping region the point belongs to
        double x1 = a*geom_rve.gs_minx +  xmin;
        double x2 = x1 + geom_rve.gs_minx;
        double y1 = b*geom_rve.gs_miny +  ymin;
        double y2 = y1 + geom_rve.gs_miny;
        double z1 = c*geom_rve.gs_minz +  zmin;
        double z2 = z1 + geom_rve.gs_minz;
        
        //Initialize flags for overlaping regions
        int fx = 0;
        int fy = 0;
        int fz = 0;
        
        //Assign value of flag according to position of point
        //The first operand eliminates the periodicity on the boundary
        if ((x > overlap_max_cutoff + xmin) && (x >= x1) && (x <= x1+overlap_max_cutoff))
            fx = -1;
        else if ((x < geom_rve.len_x+xmin-overlap_max_cutoff) && (x >= x2-overlap_max_cutoff) && (x <= x2 ))
            fx = 1;
        if ((y > overlap_max_cutoff + ymin) && (y >= y1) && (y <= y1+overlap_max_cutoff))
            fy = -1;
        else if ((y < geom_rve.wid_y+ymin-overlap_max_cutoff) && (y >= y2-overlap_max_cutoff) && (y <= y2 ))
            fy = 1;
        if ((z > overlap_max_cutoff + zmin) && (z >= z1) && (z <= z1+overlap_max_cutoff))
            fz = -1;
        else if ((z < geom_rve.hei_z+zmin-overlap_max_cutoff) && (z >= z2-overlap_max_cutoff) && (z <= z2 ))
            fz = 1;
        
        //Create array for loop over overlaping regions
        int temp[2][3] = { {a+fx, b+fy, c+fz}, {a, b, c}};
        int t;
        
        //In this loop I check all regions a point can belong to when it is in an overlaping zone
        for (int ii = 0; ii < 2; ii++) {
            if (!fx) ii++; //if flag is zero, do this loop only once
            for (int jj = 0; jj < 2; jj++) {
                if (!fy) jj++; //if flag is zero, do this loop only once
                for (int kk = 0; kk < 2; kk++) {
                    if (!fz) kk++; //if flag is zero, do this loop only once
                    //hout <<"a="<<a<<" fx="<<fx<<" b="<<b<<" fy="<<fy<<" c="<<c<<" fz="<<fz;
                    t = Calculate_t(temp[ii][0],temp[jj][1],temp[kk][2],n_subregions[0],n_subregions[1]);
                    //hout<<" t="<<t<<" sectioned_domain["<<t<<"].size()="<<sectioned_domain[t].size();
                    sectioned_domain[t].push_back(global_num);
                    //hout<<'.'<<endl;
                }
            }
        }


    }
}
//Calculates the region to which a point corresponds
int GenNetwork::Calculate_t(int a, int b, int c, int sx, int sy)const
{
    return a + b*sx + c*sx*sy;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Generate a number of ellipsoids
int GenNetwork::Get_ellip_clusters(const struct cuboid &cub, const struct Cluster_Geo &clust_geo)const
{
    double epsilon = 0.01;						//A ratio for extending ellipsoids
    double ellip_volume = 0.0;
    vector<struct elliparam> ellips;			//Define the temporary vector of ellipsoids for nanotube cluster zones
    double real_volume_fraction;				//Define the real volume fraction of ellips in the RVE
    
    const int N_times=1000;					//A maximum number for generation
    int times = 0;										//Count the number of generation
    do
    {
        //-------------------------------------------------------------
        //Ready to generate an ellipsoid
        struct elliparam ell_temp;
        //Generate the center point of an ellipsoid
        ell_temp.x=cub.poi_min.x + ((double)rand()/RAND_MAX)*cub.len_x;
        
        ell_temp.y=cub.poi_min.y + ((double)rand()/RAND_MAX)*cub.wid_y;
        
        ell_temp.z=cub.poi_min.z + ((double)rand()/RAND_MAX)*cub.hei_z;
        
        //Generate the lengths of half-axes of an ellipsoid
        ell_temp.a=clust_geo.amin + ((double)rand()/RAND_MAX)*(clust_geo.amax - clust_geo.amin);
        if(!(clust_geo.bmin==0&&clust_geo.cmin==0))
        {
            ell_temp.b = clust_geo.bmin + ((double)rand()/RAND_MAX)*(ell_temp.a - clust_geo.bmin);
            
            ell_temp.c = clust_geo.cmin + ((double)rand()/RAND_MAX)*(ell_temp.b - clust_geo.cmin);
        }
        else
        {
            ell_temp.b = ell_temp.a;
            ell_temp.c = ell_temp.a;
        }
        
        //Generate 9 angles: [(alpha1,beta1,gamma1),(alpha2,beta2,gamma2),(alpha3,beta3,gamma3)]
        //between three axes of ellipsoid (a,b,c) with three coordinate axes (ox,oy,oz)
        double alpha1 = ((double)rand()/RAND_MAX)*PI;
        double beta1 = 0;
        if(alpha1>PI/2.0)
        {
            beta1 = (alpha1-PI/2.0) + ((double)rand()/RAND_MAX)*2*(PI-alpha1);
        }
        else
        {
            beta1 = (PI/2.0-alpha1) + ((double)rand()/RAND_MAX)*2*alpha1;
        }
        
        ell_temp.alpha1	=	cos(alpha1);																//alpha1 is chosen from (0, PI)
        ell_temp.beta1	=	cos(beta1);																//beta1 is chosen from (pi/2-r1) to (pi/2+r1)
        ell_temp.gamma1 = pow(-1.0, fmod(rand(), 2.0)+1.0)*sqrt(1.0-pow(ell_temp.alpha1,2)-pow(ell_temp.beta1,2));	 //Calculate the value of gamma but randomly choose "positive" or "negative"
        double alpha2 = 0;																					//alpha2 is chosen from (pi/2-r1) to (pi/2+r1)
        if(alpha1>PI/2.0)
        {
            alpha2  = (alpha1-PI/2.0) + ((double)rand()/RAND_MAX)*2*(PI-alpha1);
        }
        else
        {
            alpha2  = (PI/2.0-alpha1) + ((double)rand()/RAND_MAX)*2*alpha1;
        }
        ell_temp.alpha2 = cos(alpha2);
        
        double A, B, C;
        A = 1+pow(ell_temp.beta1/ell_temp.gamma1,2);
        B = 2*(ell_temp.alpha1*ell_temp.alpha2*ell_temp.beta1)/pow(ell_temp.gamma1,2);
        C = pow(ell_temp.alpha1*ell_temp.alpha2/ell_temp.gamma1,2)+pow(ell_temp.alpha2,2)-1.0;
        
        ell_temp.beta2 = (-B+pow(-1.0, fmod(rand(),2.0)+1)*sqrt(pow(B,2)-4*A*C))/(2*A);
        ell_temp.gamma2 = -(ell_temp.beta1/ell_temp.gamma1)*ell_temp.beta2-(ell_temp.alpha1*ell_temp.alpha2/ell_temp.gamma1);
        
        double sign;
        sign = (ell_temp.alpha1*ell_temp.beta2)/fabs(ell_temp.alpha1*ell_temp.beta2);
        ell_temp.alpha3 = sign*sqrt(1-pow(ell_temp.alpha1,2)-pow(ell_temp.alpha2,2));
        ell_temp.beta3 = -(ell_temp.alpha1*ell_temp.beta1+ell_temp.alpha2*ell_temp.beta2)/ell_temp.alpha3;
        ell_temp.gamma3 = -(ell_temp.alpha1*ell_temp.gamma1+ell_temp.alpha2*ell_temp.gamma2)/ell_temp.alpha3;
        
        ell_temp.a = (1+epsilon)*ell_temp.a;          //Extend axes of ellipsoid a little bit for checking intersection or not
        ell_temp.b = (1+epsilon)*ell_temp.b;
        ell_temp.c = (1+epsilon)*ell_temp.c;
        
        //-------------------------------------------------------------
        //To check if an intersection happens between this ellipsoid with the surfaces of RVE or other generated ellipsoids
        double delt_h = ell_temp.c/50;						//Attention: if devided too much, it will spend too much computing time
        int k1 = (int)(sqrt(pow(ell_temp.a,2)+pow(ell_temp.b,2))/delt_h);
        int K = 4*(k1+1);
        double sita = 2*PI/K;
        
        //To check if an intersection happens between this ellipsoid with the surfaces of RVE
        for(int i=0; i<=K/2; i++)
        {
            int l1 = (int)(sqrt(pow(ell_temp.a*sin(i*sita),2)+pow(ell_temp.b*sin(i*sita),2))/delt_h);
            int L = 4*(l1+1);
            double phi = 2*PI/L;
            
            for(int j=1; j<=L; j++)
            {
                double x=ell_temp.a*sin(i*sita)*cos(j*phi);
                double y=ell_temp.b*sin(i*sita)*sin(j*phi);
                double z=ell_temp.c*cos(i*sita);
                
                double x1=ell_temp.x+x*ell_temp.alpha1+y*ell_temp.alpha2+z*ell_temp.alpha3;
                double y1=ell_temp.y+x*ell_temp.beta1+y*ell_temp.beta2+z*ell_temp.beta3;
                double z1=ell_temp.z+x*ell_temp.gamma1+y*ell_temp.gamma2+z*ell_temp.gamma3;
                
                if(x1-cub.poi_min.x<Zero||x1-cub.poi_min.x>cub.len_x-Zero||
                   y1-cub.poi_min.y<Zero||y1-cub.poi_min.y>=cub.wid_y-Zero||
                   z1-cub.poi_min.z<Zero||z1-cub.poi_min.z>=cub.hei_z-Zero)
                {
                    times=times+1;
                    goto gen_again;
                }
            }
        }
        //To check if an intersection happens between this ellipsoid with other generated ellipsoids
        for(int i=0; i<(int)ellips.size(); i++)
        {
            //Rough estimate
            double dist = sqrt(pow(ell_temp.x-ellips[i].x, 2) + pow(ell_temp.y-ellips[i].y, 2) + pow(ell_temp.z-ellips[i].z, 2));
            
            if(dist>ell_temp.a+ellips[i].a+Zero)
            {
                goto gene;
            }
            else if((dist<ell_temp.c+ellips[i].c+Zero))
            {
                times=times+1;
                goto gen_again;
            }
            else
            {
                //accurate estimate
                for(int j=1; j<=K/2; j++)
                {
                    int l1=(int)(sqrt(pow(ell_temp.a*sin(j*sita),2)+pow(ell_temp.b*sin(j*sita),2))/delt_h);
                    int L=4*(l1+1);
                    double phi=2*PI/L;
                    
                    for(int m=1;m<=L;m++)
                    {
                        double x=ell_temp.a*sin(j*sita)*cos(m*phi);
                        double y=ell_temp.b*sin(j*sita)*sin(m*phi);
                        double z=ell_temp.c*cos(j*sita);
                        
                        double x1=ell_temp.x+x*ell_temp.alpha1+y*ell_temp.alpha2+z*ell_temp.alpha3;
                        double y1=ell_temp.y+x*ell_temp.beta1+y*ell_temp.beta2+z*ell_temp.beta3;
                        double z1=ell_temp.z+x*ell_temp.gamma1+y*ell_temp.gamma2+z*ell_temp.gamma3;
                        
                        x=x1-ellips[i].x;
                        y=y1-ellips[i].y;
                        z=z1-ellips[i].z;
                        
                        x1=x*ellips[i].alpha1+y*ellips[i].beta1+z*ellips[i].gamma1;
                        y1=x*ellips[i].alpha2+y*ellips[i].beta2+z*ellips[i].gamma2;
                        z1=x*ellips[i].alpha3+y*ellips[i].beta3+z*ellips[i].gamma3;
                        
                        double f=pow(x1,2)/pow(ellips[i].a, 2)+pow(y1,2)/pow(ellips[i].b, 2)+pow(z1,2)/pow(ellips[i].c, 2)-1.0;
                        
                        if(f<0.0)
                        {
                            times=times+1;
                            goto gen_again;
                        }
                    }
                }
            }
        gene: ;
        }
        //---------------------------------------------------------------------
        //To clear the number of times and to insert an ellipsoid to the vector
        times=0;
        ellips.push_back(ell_temp);
        //---------------------------------------------------------------------
        //Calculate the sum of ellipsoid volume
        ellip_volume += 4*PI*ell_temp.a*ell_temp.b*ell_temp.c/(3*pow(1+epsilon, 3.0));
        real_volume_fraction = ellip_volume/cub.volume;
    gen_again:	;
    }while(times<=N_times&&real_volume_fraction<clust_geo.vol_fra_criterion);
    
    //---------------------------------------------------------------------
    //Shrink back to original ellipsoids
    for(int i=0; i<(int)ellips.size(); i++)
    {
        ellips[i].a=ellips[i].a/(1+epsilon);
        ellips[i].b=ellips[i].b/(1+epsilon);
        ellips[i].c=ellips[i].c/(1+epsilon);
    }
    //---------------------------------------------------------------------
    //Print the ellipsoid surfaces by grids
    if(clust_geo.print_key ==2)	Export_cluster_ellipsoids_mesh(cub, ellips);
    
    //---------------------------------------------------------------------
    //Export the data of ellipsoid surfaces
    if(clust_geo.print_key==1||clust_geo.print_key==2)	Export_cluster_ellipsoids_data(ellips, real_volume_fraction);
    
    //To print the number of ellipsoids and volume fraction
    hout << "    The number of clusters and the sum of their volume fraction:" << (int)ellips.size() << "  " << real_volume_fraction << endl;
    
    return 1;
}
//---------------------------------------------------------------------------
//Print the ellipsoid surfaces by grids
void GenNetwork::Export_cluster_ellipsoids_mesh(const struct cuboid &cub, const vector<struct elliparam> &ellips)const
{
    ofstream otec("Cluster_Ellipsoid_Mesh.dat");
    otec << "TITLE = Cluster_Ellipsoid_Mesh" << endl;
    otec << "VARIABLES = X, Y, Z" << endl;
    
    //---------------------------------------------------------------------------
    //Print the frame of RVE
    otec << "ZONE N=" << 8 << ", E=" << 1 << ", F=FEPOINT, ET=BRICK" << endl;
    double cell_x[2] = {cub.poi_min.x, cub.poi_min.x+cub.len_x};
    double cell_y[2] = {cub.poi_min.y, cub.poi_min.y+cub.wid_y};
    double cell_z[2] = {cub.poi_min.z, cub.poi_min.z+cub.hei_z};
    for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
            for(int k=0; k<2; k++)
            {
                otec << cell_x[i] << "  " << cell_y[j] << "  " << cell_z[k] << endl;
            }
    
    otec << "1 2 4 3 5 6 8 7" << endl;
    otec << endl;
    
    for(int i=0; i<(int)ellips.size(); i++)
    {
        const int num_sita = 20;
        const int num_phi = int(2*num_sita*ellips[i].a/ellips[i].c+0.5);		//Rounded to the nearest whole number
        otec << "ZONE I=" << num_phi+1 << ", J=" << num_sita+1 << ", K=1, F=POINT" << endl;
        double x, y, z;
        double x1, y1, z1;
        double sita = PI/num_sita;
        double phi=2*PI/num_phi;
        for(int j=0; j<=num_sita; j++)
        {
            for(int m=1; m<=num_phi; m++)
            {
                x=ellips[i].a*sin(j*sita)*cos(m*phi);
                y=ellips[i].b*sin(j*sita)*sin(m*phi);
                z=ellips[i].c*cos(j*sita);
                
                x1=ellips[i].x+x*ellips[i].alpha1+y*ellips[i].alpha2+z*ellips[i].alpha3;
                y1=ellips[i].y+x*ellips[i].beta1+y*ellips[i].beta2+z*ellips[i].beta3;
                z1=ellips[i].z+x*ellips[i].gamma1+y*ellips[i].gamma2+z*ellips[i].gamma3;
                
                otec << x1 << "  " << y1 << "  " << z1 << endl;
            }
            
            x=ellips[i].a*sin(j*sita)*cos(phi);
            y=ellips[i].b*sin(j*sita)*sin(phi);
            z=ellips[i].c*cos(j*sita);
            
            x1=ellips[i].x+x*ellips[i].alpha1+y*ellips[i].alpha2+z*ellips[i].alpha3;
            y1=ellips[i].y+x*ellips[i].beta1+y*ellips[i].beta2+z*ellips[i].beta3;
            z1=ellips[i].z+x*ellips[i].gamma1+y*ellips[i].gamma2+z*ellips[i].gamma3;
            
            otec << x1 << "  " << y1 << "  " << z1 << endl;
        }
        otec << endl;
    }
    otec.close();
}
//---------------------------------------------------------------------
//Export the data of ellipsoid surfaces
void GenNetwork::Export_cluster_ellipsoids_data(const vector<struct elliparam> &ellips, const double &ellip_ratio)const
{
    ofstream out("Cluster_Ellipsoids_Data.dat");
    out <<"%The number of clusters and the sum of their volume fraction" << endl;
    out << (int)ellips.size() << " " << ellip_ratio <<endl;
    out <<"%Print 15 parameters: center point (x,y,z), size of axis (a,b,c), angles (alpha1,beta1,gamma1; alpha2,beta2,gamma2; alpha3,beta3,gamma3)" << endl;
    for(int i=0; i<(int)ellips.size (); i++)
    {
        out	<< i << "  "
        << ellips[i].x << " " << ellips[i].y << " " << ellips[i].z << " "
        << ellips[i].a << " " << ellips[i].b << " " << ellips[i].c << " "
        << ellips[i].alpha1 << " " << ellips[i].beta1 << " " << ellips[i].gamma1 << " "
        << ellips[i].alpha2 << " " << ellips[i].beta2 << " " << ellips[i].gamma2 << " "
        << ellips[i].alpha3 << " " << ellips[i].beta3 << " " << ellips[i].gamma3 << " "
        << endl;
    }
    out.close();
}
//---------------------------------------------------------------------------
//Generate a number of sperical clusters in regular arrangement (Increase the number of clusters which cannot be achieved by random distribution)
int GenNetwork::Get_spherical_clusters_regular_arrangement(const struct cuboid &cub, struct Cluster_Geo &clust_geo)const
{
    int snum = 2;			//The number of spheres on each side of RVE
    double sd_x = 0.5*cub.len_x/snum;
    double sd_y = 0.5*cub.wid_y/snum;
    double sd_z = 0.5*cub.hei_z/snum;
    if(sd_x<=clust_geo.amin||sd_y<=clust_geo.amin||sd_z<=clust_geo.amin) { hout << "Error: the number of spheres on each side of RVE is too many, please check again." << endl; return 0; }
    
    double real_volume_fraction;				//Define the real volume fraction of ellips in the RVE
    double ellip_volume = 0.0;
    vector<struct elliparam> ellips;			//Define the temporary vector of ellipsoids for nanotube cluster zones
    for(int i=0; i<snum; i++)
        for(int j=0; j<snum; j++)
            for(int k=0; k<snum; k++)
            {
                //-------------------------------------------------------------
                //Generate a sphere
                struct elliparam ell_temp;
                ell_temp.x	=	(2*k+1)*sd_x;
                ell_temp.y	=	(2*j+1)*sd_y;
                ell_temp.z	=	(2*i+1)*sd_z;
                
                ell_temp.a=clust_geo.amin;
                ell_temp.b = ell_temp.a;
                ell_temp.c = ell_temp.a;
                
                //Generate 9 angles: [(alpha1,beta1,gamma1),(alpha2,beta2,gamma2),(alpha3,beta3,gamma3)]
                //between three axes of ellipsoid (a,b,c) with three coordinate axes (ox,oy,oz)
                double alpha1 = ((double)rand()/RAND_MAX)*PI;
                double beta1 = 0;
                if(alpha1>PI/2.0)
                {
                    beta1 = (alpha1-PI/2.0) + ((double)rand()/RAND_MAX)*2*(PI-alpha1);
                }
                else
                {
                    beta1 = (PI/2.0-alpha1) + ((double)rand()/RAND_MAX)*2*alpha1;
                }
                
                ell_temp.alpha1	=	cos(alpha1);																//alpha1 is chosen from (0, PI)
                ell_temp.beta1	=	cos(beta1);																//beta1 is chosen from (pi/2-r1) to (pi/2+r1)
                ell_temp.gamma1 = pow(-1.0, fmod(rand(), 2.0)+1.0)*sqrt(1.0-pow(ell_temp.alpha1,2)-pow(ell_temp.beta1,2));	  //Calculate the value of gamma but randomly choose "positive" or "negative"
                double alpha2 = 0;																					//alpha2 is chosen from (pi/2-r1) to (pi/2+r1)
                if(alpha1>PI/2.0)
                {
                    alpha2  = (alpha1-PI/2.0) + ((double)rand()/RAND_MAX)*2*(PI-alpha1);
                }
                else
                {
                    alpha2  = (PI/2.0-alpha1) + ((double)rand()/RAND_MAX)*2*alpha1;
                }
                ell_temp.alpha2 = cos(alpha2);
                
                double A, B, C;
                A = 1+pow(ell_temp.beta1/ell_temp.gamma1,2);
                B = 2*(ell_temp.alpha1*ell_temp.alpha2*ell_temp.beta1)/pow(ell_temp.gamma1,2);
                C = pow(ell_temp.alpha1*ell_temp.alpha2/ell_temp.gamma1,2)+pow(ell_temp.alpha2,2)-1.0;
                
                ell_temp.beta2 = (-B+pow(-1.0, fmod(rand(),2.0)+1)*sqrt(pow(B,2)-4*A*C))/(2*A);
                ell_temp.gamma2 = -(ell_temp.beta1/ell_temp.gamma1)*ell_temp.beta2-(ell_temp.alpha1*ell_temp.alpha2/ell_temp.gamma1);
                
                double sign;
                sign = (ell_temp.alpha1*ell_temp.beta2)/fabs(ell_temp.alpha1*ell_temp.beta2);
                ell_temp.alpha3 = sign*sqrt(1-pow(ell_temp.alpha1,2)-pow(ell_temp.alpha2,2));
                ell_temp.beta3 = -(ell_temp.alpha1*ell_temp.beta1+ell_temp.alpha2*ell_temp.beta2)/ell_temp.alpha3;
                ell_temp.gamma3 = -(ell_temp.alpha1*ell_temp.gamma1+ell_temp.alpha2*ell_temp.gamma2)/ell_temp.alpha3;
                
                //---------------------------------------------------------------------
                //To insert an ellipsoid to the vector
                ellips.push_back(ell_temp);
                
                //---------------------------------------------------------------------
                //Calculate the sum of ellipsoid volume
                ellip_volume += 4*PI*ell_temp.a*ell_temp.b*ell_temp.c/3;
                real_volume_fraction = ellip_volume/cub.volume;
            }
    
    //To check if the volume fraction is less than the criterion value
    if(real_volume_fraction<clust_geo.vol_fra_criterion)
    {
        hout << "The sum of volume fraction of spheres: " << real_volume_fraction;
        hout << " which is less than the criterion value: " << clust_geo.vol_fra_criterion << " , please check it again!" << endl;
        return 0;
    }
    
    //---------------------------------------------------------------------
    //Print the ellipsoid surfaces by grids
    if(clust_geo.print_key==2)	Export_cluster_ellipsoids_mesh(cub, ellips);
    
    //---------------------------------------------------------------------
    //Export the data of ellipsoid surfaces
    if(clust_geo.print_key==1||clust_geo.print_key==2)	Export_cluster_ellipsoids_data(ellips, real_volume_fraction);
    
    //To print the number of ellipsoids and volume fraction
    hout << "    The number of clusters and the sum of their volume fraction:" << (int)ellips.size() << "  " << real_volume_fraction << endl;
    
    return 1;
}
//---------------------------------------------------------------------------
//Generate a random value through a probability distribution function
int GenNetwork::Get_random_value(const string &dist_type, const double &min, const double &max, double &value)const
{
    if(min>max) { hout << "Error, the minimum value is larger than the maximum value (Get_random_value)!" << endl; return 0; }
    
    if(dist_type=="uniform")	//uniform distribution
    {
        value = (max-min)*((double)rand()/RAND_MAX) + min;
    }
    else if(dist_type=="normal")	//normal distribution
    {
        long long int sum=0;
        for(int i=0; i<12; i++)
        {
            sum = sum + rand();
        }
        value = ((double)sum/RAND_MAX-6.0)*(max-min)/12.0 + 0.5*(max+min);
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Randomly generate a seed (intial point) of a CNT in the RVE
int GenNetwork::Get_seed_point(const struct cuboid &cub, Point_3D &point)const
{
    
    point.x = cub.poi_min.x + cub.len_x*((double)rand()/RAND_MAX);
    
    point.y = cub.poi_min.y + cub.wid_y*((double)rand()/RAND_MAX);
    
    point.z = cub.poi_min.z + cub.hei_z*((double)rand()/RAND_MAX);//*/

    
    point.flag = 0; //0 denotes this point is the initial point of a CNT
    
    return 1;
}
//---------------------------------------------------------------------------
//Randomly generate a direction in the spherical coordinates as the initial direction of CNT segments
int GenNetwork::Get_uniform_direction(const struct Nanotube_Geo &nanotube_geo, double &cnt_sita, double &cnt_pha)const
{
    if(nanotube_geo.dir_distrib_type=="random")
    {
        //Using only rand()
        //sita is chosen in [0, PI] with uniform distribution
        cnt_sita = PI*((double)rand()/RAND_MAX);
        
        //pha is chosen in [0, 2PI] with uniform distribution
        cnt_pha = 2.0*PI*((double)rand()/RAND_MAX);//*/

    }
    else if(nanotube_geo.dir_distrib_type=="specific")
    {
        //Use the probability of a random number to be even
        if(rand()%2==0)
        {
            cnt_sita = nanotube_geo.ini_sita;		//"positive" direction
            cnt_pha = nanotube_geo.ini_pha;		//"positive" direction
        }
        else
        {
            cnt_sita = PI - nanotube_geo.ini_sita;		//"negative" (opposite) direction
            cnt_pha = PI + nanotube_geo.ini_pha;	//"negative" (opposite) direction
        }
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Transform angles into matrix
MathMatrix GenNetwork::Get_transformation_matrix(const double &sita, const double &pha)const
{
    //M = M_pha*M_sita
    //          |cos(pha) -sin(pha) 0|
    // M_pha  = |sin(pha)  cos(pha) 0|
    //          |   0         0     1|
    //
    //          | cos(sita)  0  sin(sita)|
    // M_sita = |     0      1      0    |
    //          |-sin(sita)  0  cos(sita)|
    //Calculate the matrix elements directly, instead of multiplying two matrices
    MathMatrix M(3,3);
    M.element[0][0] = cos(pha)*cos(sita);
    M.element[0][1] = -sin(pha);
    M.element[0][2] = cos(pha)*sin(sita);
    
    M.element[1][0] = sin(pha)*cos(sita);
    M.element[1][1] = cos(pha);
    M.element[1][2] = sin(pha)*sin(sita);
    
    M.element[2][0] = -sin(sita);
    M.element[2][2] = cos(sita);
    
    return M;
}
//---------------------------------------------------------------------------
//Randomly generate a direction in the spherical coordinates
//To have the positive Z-axis to be a central axis
//Then, the radial angle, sita, obeys a normal distribution (sita \in fabs[(-omega,+omega)]) and the zonal angle, pha, obeys a uniform distribution (pha \in (0,2PI))
int GenNetwork::Get_normal_direction(const double &omega, double &cnt_sita, double &cnt_pha)const
{

    //Using only rand()
    //sita centres 0 and obeys a normal distribution in (-omega, +omega)
    long long int sum=0;
    for(int i=0; i<12; i++)
    {
        sum = sum + rand();
    }
    double sum_d = (double)sum;
    cnt_sita = fabs((sum_d*omega)/(6.0*((double)RAND_MAX))-omega);
    
    //pha satisfies a uniform distribution in (0, 2PI)
    cnt_pha = 2.0*PI*((double)rand()/RAND_MAX);//*/
    
    return 1;
}
//---------------------------------------------------------------------------
//To calculate the coordinates of the new CNT point (transformation of coordinates)
Point_3D GenNetwork::Get_new_point(MathMatrix &Matrix, const double &Rad)const
{
    //Point = Matrix*v
    //v = [0; 0; Rad]
    //Calculate the new point directly
    Point_3D Point(Matrix.element[0][2]*Rad, Matrix.element[1][2]*Rad, Matrix.element[2][2]*Rad);
    
    return Point;
}
//---------------------------------------------------------------------------
//To judge if a point is included in a cuboid
int GenNetwork::Judge_RVE_including_point(const struct cuboid &cub, const Point_3D &point)const
{
    if(point.x<cub.poi_min.x||point.x>cub.poi_min.x+cub.len_x||
       point.y<cub.poi_min.y||point.y>cub.poi_min.y+cub.wid_y||
       point.z<cub.poi_min.z||point.z>cub.poi_min.z+cub.hei_z) return 0;
    
    return 1;
}
//---------------------------------------------------------------------------
//To judge if a point is included in a RVE
int GenNetwork::Judge_RVE_including_point(const struct Geom_RVE &geom_rve, const Point_3D &point)const
{
    if(point.x<geom_rve.origin.x||point.x>geom_rve.origin.x+geom_rve.len_x||
       point.y<geom_rve.origin.y||point.y>geom_rve.origin.y+geom_rve.wid_y||
       point.z<geom_rve.origin.z||point.z>geom_rve.origin.z+geom_rve.hei_z) return 0;
    
    return 1;
}
//---------------------------------------------------------------------------
//Calculate all intersection points between the new segment and surfaces of RVE
//(using a parametric equatio:  the parameter 0<t<1, and sort all intersection points from the smaller t to the greater t)
int GenNetwork::Get_intersecting_point_RVE_surface(const struct cuboid &cub, const Point_3D &point0, const Point_3D &point1, vector<Point_3D> &ipoi_vec)const
{
    double t_temp[6];
    //The planes (surfaces of RVE) perpendicular to X axis
    t_temp[0] = (cub.poi_min.x - point0.x)/(point1.x - point0.x);
    t_temp[1] = (cub.poi_min.x + cub.len_x - point0.x)/(point1.x - point0.x);
    //The planes (surfaces of RVE) perpendicular to Y axis
    t_temp[2] = (cub.poi_min.y - point0.y)/(point1.y - point0.y);
    t_temp[3] = (cub.poi_min.y + cub.wid_y - point0.y)/(point1.y - point0.y);
    //The planes (surfaces of RVE) perpendicular to Z axis
    t_temp[4] = (cub.poi_min.z - point0.z)/(point1.z - point0.z);
    t_temp[5] = (cub.poi_min.z + cub.hei_z - point0.z)/(point1.z - point0.z);
    
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
        hout << "Cuboid P_min="<<cub.poi_min.x<<' '<<cub.poi_min.y<<' '<<cub.poi_min.z<<endl;
        hout << "Cuboid size="<<cub.len_x<<' '<<cub.wid_y<<' '<<cub.hei_z<<endl;
        hout << "P0= "<<point0.x<<' '<<point0.y<<' '<<point0.z<<' '<<endl;
        hout << "Judge_RVE_including_point="<<Judge_RVE_including_point(cub, point0)<<endl;
        hout << "P1= "<<point1.x<<' '<<point1.y<<' '<<point1.z<<' '<<endl;
        hout << "Judge_RVE_including_point="<<Judge_RVE_including_point(cub, point1)<<endl;
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
        if(fabs(point_temp.x-cub.poi_min.x)<Zero) point_temp.x = cub.poi_min.x;
        else if(fabs(point_temp.x-cub.poi_min.x-cub.len_x)<Zero) point_temp.x = cub.poi_min.x + cub.len_x;
        
        if(fabs(point_temp.y-cub.poi_min.y)<Zero) point_temp.y = cub.poi_min.y;
        else if(fabs(point_temp.y-cub.poi_min.y-cub.wid_y)<Zero) point_temp.y = cub.poi_min.y + cub.wid_y;
        
        if(fabs(point_temp.z-cub.poi_min.z)<Zero) point_temp.z = cub.poi_min.z;
        else if(fabs(point_temp.z-cub.poi_min.z-cub.hei_z)<Zero) point_temp.z = cub.poi_min.z + cub.hei_z;
        
        //---------------------------------------------------------------------------
        //Insert a new point
        ipoi_vec.push_back(point_temp);
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//To calculate the effective portion (length) which falls into the given region (RVE)
double GenNetwork::Effective_length_given_region(const struct cuboid &cub, const Point_3D last_point, const Point_3D new_point)const
{
    //Check if the last point is inside the given region
    int last_bool = Judge_RVE_including_point(cub, last_point);
    //Check if the new point is inside the given region
    int new_bool = Judge_RVE_including_point(cub, new_point);
    
    //Vector to store the intersecting point
    vector<Point_3D> ipoi_vec;
    
    //Decide the corresponding case and calculate volume fraction
    if (last_bool&&new_bool)
        return last_point.distance_to(new_point); //both points are inside so add the total length
    else if (last_bool&&(!new_bool))  //if the last point is inside and the new point is outside
    {
        if(Get_intersecting_point_RVE_surface(cub, last_point, new_point, ipoi_vec)==0){
            hout << "Error in Effective_length_given_region, case last_bool&&(!new_bool) "<<endl;
            return 0;
        }
        return last_point.distance_to(ipoi_vec[0]);
    }
    else if ((!last_bool)&&new_bool)  //if the last point is outside and the new point is inside
    {
        if(Get_intersecting_point_RVE_surface(cub, new_point, last_point, ipoi_vec)==0) {
            hout << "Error in Effective_length_given_region, case (!last_bool)&&new_bool"<<endl;
            return 0;
        }
        return new_point.distance_to(ipoi_vec[0]);
    }
    else
        return 0.0; //if both points are outside
}
//---------------------------------------------------------------------------
//Checking the angle between two segments in one nanotube (if less than PI/2, provide an alarm)
int GenNetwork::CNTs_quality_testing(const vector<vector<Point_3D> > &cnts_points)const
{
    //---------------------------------------------------------------------------
    //Checking the angle between two segments
    for(int i=0; i<(int)cnts_points.size(); i++)
    {
        for (int j=1; j<(int)cnts_points[i].size()-1; j++)
        {
            double x01 = cnts_points[i][j-1].x - cnts_points[i][j].x;
            double y01 = cnts_points[i][j-1].y - cnts_points[i][j].y;
            double z01 = cnts_points[i][j-1].z - cnts_points[i][j].z;
            
            double x21 = cnts_points[i][j+1].x - cnts_points[i][j].x;
            double y21 = cnts_points[i][j+1].y - cnts_points[i][j].y;
            double z21 = cnts_points[i][j+1].z - cnts_points[i][j].z;
            
            //Checking the angle between two segments (nod21ºÍnod01)
            double cos_ang210 = (x21*x01+y21*y01+z21*z01)/(sqrt(x21*x21+y21*y21+z21*z21)*sqrt(x01*x01+y01*y01+z01*z01));
            
            //Check if cos_ang210 is small enough to be considered zero
            //if (abs(cos_ang210) < Zero)
                //cos_ang210 = 0;
            
            //Judge the cos value(<= 90degree and >=0degree, cos value<=1.0 and >=0.0)
            if( (cos_ang210<=1.0-Zero )&& cos_ang210>=Zero )
            {
                hout << "Error: there exists at least one angle which is larger than PI/2!" << endl;
                hout << setwp(1,20)<<"cos_ang210="<<cos_ang210<<endl;
                hout << "CNT# = "<<i<<endl;
                hout << "P1 ("<<j-1<<") = "<<cnts_points[i][j-1].x<<' '<<cnts_points[i][j-1].y<<' '<<cnts_points[i][j-1].z<<endl;
                hout << "P2 ("<<j<<") = "<<cnts_points[i][j].x<<' '<<cnts_points[i][j].y<<' '<<cnts_points[i][j].z<<endl;
                hout << "P3 ("<<j+1<<") = "<<cnts_points[i][j+1].x<<' '<<cnts_points[i][j+1].y<<' '<<cnts_points[i][j+1].z<<endl;
                return 0;
            }
        }
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Generate the nodes and tetrahedron elements of nanotubes (No const following this function because a sum operation on two Point_3D points inside)
int GenNetwork::Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius)
{
    //Looping the generated nanotubes
    for(int i=0; i<(int)cnts_points.size(); i++)
    {
        vector<Node> nod_temp;
        vector<Element> ele_temp;
        
        const int cps = (int)cnts_points[i].size();
        for(int j=0; j<cps; j++)
        {
            //Calculate the normal vector of the plane
            Point_3D plane_normal;
            if(j==0)
            {
                plane_normal.x = cnts_points[i][j].x - cnts_points[i][j+1].x;
                plane_normal.y = cnts_points[i][j].y - cnts_points[i][j+1].y;
                plane_normal.z = cnts_points[i][j].z - cnts_points[i][j+1].z;
            }
            else if(j==cps-1)
            {
                plane_normal.x = cnts_points[i][j-1].x - cnts_points[i][j].x;
                plane_normal.y = cnts_points[i][j-1].y - cnts_points[i][j].y;
                plane_normal.z = cnts_points[i][j-1].z - cnts_points[i][j].z;
            }
            else
            {
                const Point_3D vect[3] = { cnts_points[i][j-1], cnts_points[i][j], cnts_points[i][j+1] };
                const double A = pow(vect[0].x-vect[1].x,2)+pow(vect[0].y-vect[1].y,2)+pow(vect[0].z-vect[1].z,2);
                const double B = pow(vect[2].x-vect[1].x,2)+pow(vect[2].y-vect[1].y,2)+pow(vect[2].z-vect[1].z,2);
                const double tt = sqrt(A/B);
                
                //Coordinate transformation to find the intersection points
                double x, y, z;
                x=vect[1].x+tt*(vect[2].x-vect[1].x);
                y=vect[1].y+tt*(vect[2].y-vect[1].y);
                z=vect[1].z+tt*(vect[2].z-vect[1].z);
                
                plane_normal.x = vect[0].x - x;
                plane_normal.y = vect[0].y - y;
                plane_normal.z = vect[0].z - z;
            }
            //The center point of the circle on the plane
            Point_3D plane_center = cnts_points[i][j];
            
            //Define the number of sections along the circumference
            const int num_sec = 36;
            if(j==0)
            {
                double normal_sita, normal_pha;  //Direction angles
                //Calculate the angles of the normal verctor of the plane in the spherical coordinate
                if(Get_angles_vector_in_spherial_coordinates(plane_normal, normal_sita, normal_pha)==0) return 0;
                
                //Calculate a group of equidistant points along the circumference which is on the plane defined by the center point of the circle and the normal vector
                if(Get_points_circle_in_plane(plane_center, normal_sita, normal_pha, cnts_radius[i], num_sec, nod_temp)==0) return 0;
            }
            else
            {
                //Calculate a group of projected points (which are on the plane with the center point of the circle and the normal vector)
                //which are projected from a group of points on the previous circumference and projected along the direction of line_vec
                Point_3D line_vec;
                line_vec.x = cnts_points[i][j-1].x - cnts_points[i][j].x;
                line_vec.y = cnts_points[i][j-1].y - cnts_points[i][j].y;
                line_vec.z = cnts_points[i][j-1].z - cnts_points[i][j].z;
                if(Get_projected_points_in_plane(plane_center, plane_normal, line_vec, num_sec, nod_temp)==0) return 0;
            }
            
            //Generate a vector of elements
            if(j!=0)
            {
                int nodes_num[6];
                nodes_num[0] = (j-1)*(num_sec+1);   //The number of the center
                nodes_num[3] = j*(num_sec+1);
                for(int k=1; k<=num_sec; k++)
                {
                    nodes_num[1] = (j-1)*(num_sec+1) + k;
                    nodes_num[2] = (j-1)*(num_sec+1) + 1 + k%num_sec;
                    nodes_num[4] = j*(num_sec+1) + k;
                    nodes_num[5] = j*(num_sec+1) + 1 + k%num_sec;
                    
                    Element eles_num[3];
                    //----------------------------------------------------------------
                    //Insert the numbers of nodes to the elements
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
                    //Insert the number of elements to element vector
                    ele_temp.push_back(eles_num[0]);
                    ele_temp.push_back(eles_num[1]);
                    ele_temp.push_back(eles_num[2]);
                }
            }
        }
        
        nodes.push_back(nod_temp);
        eles.push_back(ele_temp);
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Generate the nodes and tetrahedron elements of nanotubes (No const following this function because a sum operation on two Point_3D points inside). This function uses a 1D point vector and a 2D structure vector that references the point vector
int GenNetwork::Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<Point_3D> &cnts_points, const vector<double> &cnts_radius, const vector<vector<long int> > &structure)
{
    //Looping the generated nanotubes
    for(int i=0; i<(int)structure.size(); i++)
    {
        vector<Node> nod_temp;
        vector<Element> ele_temp;
        
        const int cps = (int)structure[i].size();
        for(int j=0; j<cps; j++)
        {
            //Calculate the normal vector of the plane
            Point_3D plane_normal;
            if(j==0)
            {
                long int P1 = structure[i][j];
                long int P2 = structure[i][j+1];
                plane_normal.x = cnts_points[P1].x - cnts_points[P2].x;
                plane_normal.y = cnts_points[P1].y - cnts_points[P2].y;
                plane_normal.z = cnts_points[P1].z - cnts_points[P2].z;
            }
            else if(j==cps-1)
            {
                long int P1 = structure[i][j-1];
                long int P2 = structure[i][j];
                plane_normal.x = cnts_points[P1].x - cnts_points[P2].x;
                plane_normal.y = cnts_points[P1].y - cnts_points[P2].y;
                plane_normal.z = cnts_points[P1].z - cnts_points[P2].z;
            }
            else
            {
                long int P1 = structure[i][j-1];
                long int P2 = structure[i][j];
                long int P3 = structure[i][j+1];
                const Point_3D vect[3] = { cnts_points[P1], cnts_points[P2], cnts_points[P3] };
                const double A = pow(vect[0].x-vect[1].x,2)+pow(vect[0].y-vect[1].y,2)+pow(vect[0].z-vect[1].z,2);
                const double B = pow(vect[2].x-vect[1].x,2)+pow(vect[2].y-vect[1].y,2)+pow(vect[2].z-vect[1].z,2);
                const double tt = sqrt(A/B);
                
                //Coordinate transformation to find the intersection points
                double x, y, z;
                x=vect[1].x+tt*(vect[2].x-vect[1].x);
                y=vect[1].y+tt*(vect[2].y-vect[1].y);
                z=vect[1].z+tt*(vect[2].z-vect[1].z);
                
                plane_normal.x = vect[0].x - x;
                plane_normal.y = vect[0].y - y;
                plane_normal.z = vect[0].z - z;
            }
            //The center point of the circle on the plane
            Point_3D plane_center = cnts_points[structure[i][j]];
            
            //Define the number of sections along the circumference
            const int num_sec = 36;
            if(j==0)
            {
                double normal_sita, normal_pha;  //Direction angles
                //Calculate the angles of the normal verctor of the plane in the spherical coordinate
                if(Get_angles_vector_in_spherial_coordinates(plane_normal, normal_sita, normal_pha)==0) return 0;
                
                //Calculate a group of equidistant points along the circumference which is on the plane defined by the center point of the circle and the normal vector
                if(Get_points_circle_in_plane(plane_center, normal_sita, normal_pha, cnts_radius[i], num_sec, nod_temp)==0) return 0;
            }
            else
            {
                //Calculate a group of projected points (which are on the plane with the center point of the circle and the normal vector)
                //which are projected from a group of points on the previous circumference and projected along the direction of line_vec
                Point_3D line_vec;
                long int P1 = structure[i][j-1];
                long int P2 = structure[i][j];
                line_vec.x = cnts_points[P1].x - cnts_points[P2].x;
                line_vec.y = cnts_points[P1].y - cnts_points[P2].y;
                line_vec.z = cnts_points[P1].z - cnts_points[P2].z;
                if(Get_projected_points_in_plane(plane_center, plane_normal, line_vec, num_sec, nod_temp)==0) return 0;
            }
            
            //Generate a vector of elements
            if(j!=0)
            {
                int nodes_num[6];
                nodes_num[0] = (j-1)*(num_sec+1);   //The number of the center
                nodes_num[3] = j*(num_sec+1);
                for(int k=1; k<=num_sec; k++)
                {
                    nodes_num[1] = (j-1)*(num_sec+1) + k;
                    nodes_num[2] = (j-1)*(num_sec+1) + 1 + k%num_sec;
                    nodes_num[4] = j*(num_sec+1) + k;
                    nodes_num[5] = j*(num_sec+1) + 1 + k%num_sec;
                    
                    Element eles_num[3];
                    //----------------------------------------------------------------
                    //Insert the numbers of nodes to the elements
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
                    //Insert the number of elements to element vector
                    ele_temp.push_back(eles_num[0]);
                    ele_temp.push_back(eles_num[1]);
                    ele_temp.push_back(eles_num[2]);
                }
            }
        }
        
        nodes.push_back(nod_temp);
        eles.push_back(ele_temp);
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Calculate the angles of a verctor in the spherical coordinate
int GenNetwork::Get_angles_vector_in_spherial_coordinates(const Point_3D &normal, double &sita, double &pha)const
{
    if(normal.x==0&&normal.y==0&&normal.z==0) { hout << "Error, three elements of the vector are all zero!" << endl; return 0; }
    sita =  acos(normal.z/sqrt(normal.x*normal.x+normal.y*normal.y+normal.z*normal.z));
    if(normal.x==0&&normal.y==0) pha = 0;
    else if(normal.y>=0) pha = acos(normal.x/sqrt(normal.x*normal.x+normal.y*normal.y));
    else if(normal.y<0) pha = 2*PI - acos(normal.x/sqrt(normal.x*normal.x+normal.y*normal.y));
    
    return 1;
}
//---------------------------------------------------------------------------
//Calculate a group of equidistant points along the circumference which is on the plane defined by the center point of the circle and the normal vector
int GenNetwork::Get_points_circle_in_plane(const Point_3D &center, const double &trans_sita, const double &trans_pha, const double &radius, const int &num_sec, vector<Node> &nod_temp)const
{
    //Insert the center point firstly
    Node new_node(center.x, center.y, center.z);
    nod_temp.push_back(new_node);
    
    //Define the transformation matrix
    MathMatrix trans_mat(3,3);
    trans_mat = Get_transformation_matrix(trans_sita, trans_pha);
    
    //1D vector defined by a matrix
    MathMatrix Rvec(3,1);
    Rvec.element[0][0] = 0;
    Rvec.element[1][0] = 0;
    Rvec.element[2][0] = radius;
    
    //1D vector defined by a matrix
    MathMatrix Res(3,1);
    
    double sita, pha;
    sita = 0.5*PI;	//Defined on the XOY plane
    for(int i=0; i<num_sec; i++)
    {
        pha = i*2*PI/num_sec;
        MathMatrix matrix_temp = trans_mat*Get_transformation_matrix(sita, pha);
        Res = matrix_temp*Rvec;
        
        new_node.x = center.x + Res.element[0][0];
        new_node.y = center.y + Res.element[1][0];
        new_node.z = center.z + Res.element[2][0]; 
        
        //Insert the points on the circumference
        nod_temp.push_back(new_node);
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Calculate a group of projected points (which are on the plane with the center point of the circle and the normal vector) 
//which are projected from a group of points on the previous circumference and projected along the direction of line_vec
int GenNetwork::Get_projected_points_in_plane(const Point_3D &center, const Point_3D &normal, const Point_3D &line, const int &num_sec, vector<Node> &nod_temp)const
{
    //Record the total number of nodes after the previous generation
    const int nod_size = (int)nod_temp.size();  
    
    //Insert the center point
    Node new_node(center.x, center.y, center.z);
    nod_temp.push_back(new_node);
    
    const double vectors_dot_product = normal.x*line.x+normal.y*line.y+normal.z*line.z;
    
    if(vectors_dot_product==0.0) 
    {
        //Corresponding to three points: number 0, 1 and 2, the peak of this angle is at the point number 1. 
        hout << "Error: these two normal vectors are perpendicular to each other!" << endl;
        return 0; 
    }
    
    for(int i=num_sec; i>0; i--)
    {
        Point_3D point(center.x-nod_temp[nod_size-i].x, center.y-nod_temp[nod_size-i].y, center.z-nod_temp[nod_size-i].z);
        new_node.x = nod_temp[nod_size-i].x + (normal.x*point.x+normal.y*point.y+normal.z*point.z)*line.x/ vectors_dot_product;
        new_node.y = nod_temp[nod_size-i].y + (normal.x*point.x+normal.y*point.y+normal.z*point.z)*line.y/ vectors_dot_product;
        new_node.z = nod_temp[nod_size-i].z + (normal.x*point.x+normal.y*point.y+normal.z*point.z)*line.z/ vectors_dot_product;
        
        //Insert the points on the circumference
        nod_temp.push_back(new_node);
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Transform the 2D cnts_points into 1D cpoints and 2D cstructuers
int GenNetwork::Transform_cnts_points(const vector<vector<Point_3D> > &cnts_points, vector<Point_3D> &cpoints, vector<vector<long int> > &cstructures)const
{
    hout << "There are "<<cnts_points.size()<<" CNTs."<<endl;
    long int count = 0;
    for(int i=0; i<(int)cnts_points.size(); i++)
    {
        vector<long int> struct_temp;
        for(int j=0; j<(int)cnts_points[i].size(); j++)
        {
            cpoints.push_back(cnts_points[i][j]);
            cpoints.back().flag = i;
            struct_temp.push_back(count);
            count++;
        }
        cstructures.push_back(struct_temp);
    }
    
    hout << "There are "<<cpoints.size() << " points."<<endl;
    
    return 1;
}
//===========================================================================
