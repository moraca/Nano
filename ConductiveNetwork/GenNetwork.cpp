//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	GenNetwork.cpp
//OBJECTIVE:	To generate networks with overlapping
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "GenNetwork.h"

//Generate 3D nantube networks with ovelapping
int GenNetwork::Generate_nanotube_networks(const struct Geom_RVE &geom_rve, const struct Cluster_Geo &clust_geo, const struct Nanotube_Geo &nanotube_geo, 
																		vector<Point_3D> &cpoints, vector<double> &cnts_radius, vector<vector<long int> > &cstructures)const
{
	//Define a two-dimensional vector of three-dimensional points for storing the CNT threads
	vector<vector<Point_3D> > cnts_points;

	//Generate a network defined by points and connections 
	if(Generate_network_threads(geom_rve, clust_geo, nanotube_geo, cnts_points, cnts_radius)==0) return 0;
    
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
int GenNetwork::Generate_network_threads(const struct Geom_RVE &geom_rve, const struct Cluster_Geo &clust_geo, const struct Nanotube_Geo &nanotube_geo, 
																		vector<vector<Point_3D> > &cnts_points,  vector<double> &cnts_radius)const
{
	//Generate random seed in terms of local time
    srand((unsigned int)time(NULL));

	//---------------------------------------------------------------------------
	//Generate the data for nanotube clusters limited in ellipsoid surfaces (each ellipsoid is within the RVE and is separated with each other)
	if(clust_geo.vol_fra_criterion>0.0)
	{
		int seed_ellip_poi = rand()%MAX_INT;
		int seed_ellip_axis = rand()%MAX_INT;
		int seed_ellip_angle = rand()%MAX_INT;

		struct cuboid cub;										//Generate a cuboid for RVE
		cub.poi_min = geom_rve.origin;
		cub.len_x = geom_rve.len_x;
		cub.wid_y = geom_rve.wid_y;
		cub.hei_z = geom_rve.hei_z;
		cub.volume = cub.len_x*cub.wid_y*cub.hei_z;

		if(Get_ellip_clusters(cub, clust_geo, seed_ellip_poi, seed_ellip_axis, seed_ellip_angle)==0) return 0;
		//Generate a number of sperical clusters in regular arrangement
//		if(Get_spherical_clusters_regular_arrangement(cub, clust_geo, seed_ellip_poi, seed_ellip_axis, seed_ellip_angle)==0) return 0;
	}

	//---------------------------------------------------------------------------
	//Generate random seeds for cnts
	int seed_cnt_origin = rand()%MAX_INT;     //MAX_INT==2^15-1 in a 16-bit computer, the range of rand() is [0, MAX_INT]
	int seed_cnt_length = rand()%MAX_INT;		//MAX_INT==2^31-1 in a 32-bit computer, the range of rand() is [0, MAX_INT]
	int seed_cnt_radius = rand()%MAX_INT;
	int seed_cnt_sita = rand()%MAX_INT;
	int seed_cnt_pha = rand()%MAX_INT;
	int seed_growth_probability = rand()%MAX_INT;

	double vol_sum = 0;  //the sum of volume of generated CNTs 
	double wt_sum = 0;   //the sum of weight of generated CNTs
	int cnt_seed_count =0; //to record the number of generated seeds of a CNT (If the growth of a CNT fails, but this seed will be used for the next growth of a CNT)
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
		if(Get_random_value(nanotube_geo.len_distrib_type, nanotube_geo.len_min, nanotube_geo.len_max, seed_cnt_length, cnt_length)==0) return 0;
		//Calculate the total number of growth step for a CNT 
		int step_num = (int)(cnt_length/nanotube_geo.step_length) + 1;

		//---------------------------------------------------------------------------
		//Randomly generate a radius of a CNT
		double cnt_rad;
		if(Get_random_value(nanotube_geo.rad_distrib_type, nanotube_geo.rad_min, nanotube_geo.rad_max, seed_cnt_radius, cnt_rad)==0) return 0;

		//---------------------------------------------------------------------------
		//Randomly generate a direction in the spherical coordinates as the original direction of CNT segments
		double cnt_sita, cnt_pha;
		if(Get_uniform_direction(nanotube_geo, seed_cnt_sita, seed_cnt_pha, cnt_sita, cnt_pha)==0) return 0;
		MathMatrix multiplier(3,3);
		multiplier = Get_transformation_matrix(cnt_sita, cnt_pha);

		//---------------------------------------------------------------------------
		//The increased volume of each segement (growth step) of nanotube (Here the overlapping volume is ignored)
		const double step_vol_para = PI*cnt_rad*cnt_rad;
		//---------------------------------------------------------------------------
		//The increased weight of each segement (growth step) of nanotube (If the different radii of nanotube are considered, the linear_density may be different in every nanotube)
		const double step_wei_para = nanotube_geo.linear_density;

		//---------------------------------------------------------------------------
        //Randomly generate a seed (original point) of a CNT in the extended RVE	(Comments: the seed generation is after the radius generation for the non-overlapping nanotubes generation)
		struct cuboid excub;										//generate a cuboid to represent the RVE
		excub.poi_min = geom_rve.ex_origin;
		excub.len_x = geom_rve.ex_len;
		excub.wid_y = geom_rve.ey_wid;
		excub.hei_z = geom_rve.ez_hei;
		excub.volume = excub.len_x*excub.wid_y*excub.hei_z;

		Point_3D cnt_poi;
		if(Get_seed_point(excub, seed_cnt_origin, cnt_poi)==0) return 0;

		new_cnt.push_back(cnt_poi);	//store this seed point in the vector for a new nanotube
		cnt_seed_count++;					//record the number of seed generations
		if(cnt_seed_count>1E6) 
		{ 
			hout << "The number of seed genrations is lager than 1E6, but the nanotube generation still fails to acheive the demanded volume fraction." << endl; 
			return 0; 
		}

		//---------------------------------------------------------------------------
		//The growth process of nanotube
		int ellip_num = -1; //For recording the serial number of ellipsoid cluster which a nanotube penetrates out. It is no use if the cluster generation is not considered.
		for(int i=0; i<step_num; i++)
		{
			//Randomly generate a direction in the spherical coordinates
			//To have the positive Z-axis to be a central axis
			//Then, the radial angle, sita, obeys a normal distribution (sita \in fabs[(-omega,+omega)]) and the zonal angle, pha, obeys a uniform distribution (pha \in (0,2PI))
			if(Get_normal_direction(nanotube_geo.angle_max, seed_cnt_sita, seed_cnt_pha, cnt_sita, cnt_pha)==0) return 0;

			//To calculate the new multiplier for transformation of coordinates
			multiplier = multiplier*Get_transformation_matrix(cnt_sita, cnt_pha);
			
			//To calculate the coordinates of the new CNT point (transformation of coordinates)
			cnt_poi = cnt_poi + Get_new_point(multiplier, nanotube_geo.step_length);
			cnt_poi.flag = 1;							//1 means that point is not the original point

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
				if(Get_intersecting_point_RVE_surface(excub, new_cnt.back(), cnt_poi, ipoi_vec)==0) return 0;
				cnt_poi = ipoi_vec[0];
				touch_end = true;
			}

			//---------------------------------------------------------------------------
			//To calculate the effective portion (length) which falls into the given region (RVE)
			struct cuboid gvcub;										//generate a cuboid to represent the give RVE
			gvcub.poi_min = geom_rve.origin;
			gvcub.len_x = geom_rve.len_x;
			gvcub.wid_y = geom_rve.wid_y;
			gvcub.hei_z = geom_rve.hei_z;
			gvcub.volume = geom_rve.volume;

            double temp_length = Effective_length_given_region(gvcub, new_cnt.back(),cnt_poi);
            //double temp_length = new_cnt.back().distance_to(cnt_poi);
            if (temp_length > 0.0) 
			{
                vol_sum += temp_length*step_vol_para;		//a accumulation on the volume
                wt_sum += temp_length*step_wei_para;		//a accumulation on the weight
            }

			new_cnt.push_back(cnt_poi);							//store a new point
            new_cnt_size = (int)new_cnt.size()-1;				//calculate the size of new point

			//---------------------------------------------------------------------------
			//Judge the new volume or weight 
			if(nanotube_geo.criterion == "vol"&&vol_sum >= nanotube_geo.real_volume) break;		//Break out when the volume reaches the critical value
			else if(nanotube_geo.criterion == "wt"&&wt_sum >= nanotube_geo.real_weight) break;		//Break out when the weight reaches the critical value
			else if (touch_end) break; //Enforce break
		}

		//---------------------------------------------------------------------------
		//To store the CNT points
		if((int)new_cnt.size()<2) 
		{
			hout << "Error, very few CNT points are generated (n<2)." << endl; 
			return 0;
		}

		vector<Point_3D> cnt_temp;
		for(int i=0; i<(int)new_cnt.size(); i++)
		{
			cnt_temp.push_back(new_cnt[i]); //insert the first CNT point
			//Judge if to the end of a CNT
			if(i==(int)new_cnt.size()-1||new_cnt[i+1].flag==0)	//two case: to the end of all CNTs or to the end of a CNT
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

	if(nanotube_geo.criterion == "wt") hout << "    The volume fraction of generated CNTs is about : " << vol_sum/geom_rve.volume << endl;

	return 1;
}
//---------------------------------------------------------------------------
//Generate a number of ellipsoids
int GenNetwork::Get_ellip_clusters(const struct cuboid &cub, const struct Cluster_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle)const
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
		seed_poi = (2053*seed_poi + 13849)%MAX_INT;
		ell_temp.x=cub.poi_min.x + seed_poi*cub.len_x/MAX_INT;
        
		seed_poi = (2053*seed_poi + 13849)%MAX_INT;
		ell_temp.y=cub.poi_min.y + seed_poi*cub.wid_y/MAX_INT;
        
		seed_poi = (2053*seed_poi + 13849)%MAX_INT;
		ell_temp.z=cub.poi_min.z + seed_poi*cub.hei_z/MAX_INT;
        
		//Generate the lengths of half-axes of an ellipsoid
		seed_axis = (2053*seed_axis + 13849)%MAX_INT;
		ell_temp.a=clust_geo.amin + seed_axis*(clust_geo.amax - clust_geo.amin)/MAX_INT;
		if(!(clust_geo.bmin==0&&clust_geo.cmin==0))
		{
			seed_axis = (2053*seed_axis + 13849)%MAX_INT;
			ell_temp.b = clust_geo.bmin + seed_axis*(ell_temp.a - clust_geo.bmin)/MAX_INT;
            
			seed_axis = (2053*seed_axis + 13849)%MAX_INT;
			ell_temp.c = clust_geo.cmin + seed_axis*(ell_temp.b - clust_geo.cmin)/MAX_INT;
		}
		else
		{
			ell_temp.b = ell_temp.a;
			ell_temp.c = ell_temp.a;
		}
        
		//Generate 9 angles: [(alpha1,beta1,gamma1),(alpha2,beta2,gamma2),(alpha3,beta3,gamma3)]  
		//between three axes of ellipsoid (a,b,c) with three coordinate axes (ox,oy,oz)
		seed_angle = (2053*seed_angle + 13849)%MAX_INT;
		double alpha1 = seed_angle*PI/MAX_INT;
		double beta1 = 0;
		if(alpha1>PI/2.0)
		{
			seed_angle = (2053*seed_angle + 13849)%MAX_INT;
			beta1 = (alpha1-PI/2.0) + seed_angle*2*(PI-alpha1)/MAX_INT;
		}
		else
		{
			seed_angle	= (2053*seed_angle + 13849)%MAX_INT;
			beta1 = (PI/2.0-alpha1) + seed_angle*2*alpha1/MAX_INT;
		}
        
		ell_temp.alpha1	=	cos(alpha1);																//alpha1 is chosen from (0, PI)
		ell_temp.beta1	=	cos(beta1);																//beta1 is chosen from (pi/2-r1) to (pi/2+r1)
		seed_angle	= (2053*seed_angle + 13849)%MAX_INT;
		ell_temp.gamma1 = pow(-1.0, fmod(seed_angle, 2.0)+1.0)*sqrt(1.0-pow(ell_temp.alpha1,2)-pow(ell_temp.beta1,2));	 //Calculate the value of gamma but randomly choose "positive" or "negative"
		double alpha2 = 0;																					//alpha2 is chosen from (pi/2-r1) to (pi/2+r1)
		if(alpha1>PI/2.0)
		{
			seed_angle	= (2053*seed_angle + 13849)%MAX_INT;
			alpha2  = (alpha1-PI/2.0) + seed_angle*2*(PI-alpha1)/MAX_INT;
		}
		else
		{
			seed_angle	= (2053*seed_angle + 13849)%MAX_INT;
			alpha2  = (PI/2.0-alpha1) + seed_angle*2*alpha1/MAX_INT;
		}
		ell_temp.alpha2 = cos(alpha2);
        
		double A, B, C;
		A = 1+pow(ell_temp.beta1/ell_temp.gamma1,2);
		B = 2*(ell_temp.alpha1*ell_temp.alpha2*ell_temp.beta1)/pow(ell_temp.gamma1,2);
		C = pow(ell_temp.alpha1*ell_temp.alpha2/ell_temp.gamma1,2)+pow(ell_temp.alpha2,2)-1.0;
        
		seed_angle	= (2053*seed_angle + 13849)%MAX_INT;
		ell_temp.beta2 = (-B+pow(-1.0, fmod(seed_angle,2.0)+1)*sqrt(pow(B,2)-4*A*C))/(2*A);
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
int GenNetwork::Get_spherical_clusters_regular_arrangement(const struct cuboid &cub, struct Cluster_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle)const
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
				seed_angle = (2053*seed_angle + 13849)%MAX_INT;
				double alpha1 = seed_angle*PI/MAX_INT;
				double beta1 = 0;
				if(alpha1>PI/2.0)
				{
					seed_angle = (2053*seed_angle + 13849)%MAX_INT;
					beta1 = (alpha1-PI/2.0) + seed_angle*2*(PI-alpha1)/MAX_INT;
				}
				else
				{
					seed_angle	= (2053*seed_angle + 13849)%MAX_INT;
					beta1 = (PI/2.0-alpha1) + seed_angle*2*alpha1/MAX_INT;
				}
                
				ell_temp.alpha1	=	cos(alpha1);																//alpha1 is chosen from (0, PI)
				ell_temp.beta1	=	cos(beta1);																//beta1 is chosen from (pi/2-r1) to (pi/2+r1)
				seed_angle	= (2053*seed_angle + 13849)%MAX_INT;
				ell_temp.gamma1 = pow(-1.0, fmod(seed_angle, 2.0)+1.0)*sqrt(1.0-pow(ell_temp.alpha1,2)-pow(ell_temp.beta1,2));	  //Calculate the value of gamma but randomly choose "positive" or "negative"
				double alpha2 = 0;																					//alpha2 is chosen from (pi/2-r1) to (pi/2+r1)
				if(alpha1>PI/2.0)
				{
					seed_angle	= (2053*seed_angle + 13849)%MAX_INT;
					alpha2  = (alpha1-PI/2.0) + seed_angle*2*(PI-alpha1)/MAX_INT;
				}
				else
				{
					seed_angle	= (2053*seed_angle + 13849)%MAX_INT;
					alpha2  = (PI/2.0-alpha1) + seed_angle*2*alpha1/MAX_INT;
				}
				ell_temp.alpha2 = cos(alpha2);
                
				double A, B, C;
				A = 1+pow(ell_temp.beta1/ell_temp.gamma1,2);
				B = 2*(ell_temp.alpha1*ell_temp.alpha2*ell_temp.beta1)/pow(ell_temp.gamma1,2);
				C = pow(ell_temp.alpha1*ell_temp.alpha2/ell_temp.gamma1,2)+pow(ell_temp.alpha2,2)-1.0;
                
				seed_angle	= (2053*seed_angle + 13849)%MAX_INT;
				ell_temp.beta2 = (-B+pow(-1.0, fmod(seed_angle,2.0)+1)*sqrt(pow(B,2)-4*A*C))/(2*A);
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
int GenNetwork::Get_random_value(const string &dist_type, const double &min, const double &max, int &seed, double &value)const
{
	if(min>max) { hout << "Error, the minimum value is larger than the maximum value (Get_random_value)!" << endl; return 0; }
    
	if(dist_type=="uniform")	//uniform distribution
	{
		seed = (2053*seed + 13849)%MAX_INT;
		value = seed*(max-min)/MAX_INT + min;
	}
	else if(dist_type=="normal")	//normal distribution
	{
		int sum=0;
		for(int i=0; i<12; i++)
		{
			seed = (2053*seed + 13849)%MAX_INT;
			sum += seed;
		}
		value = ((double)sum/MAX_INT-6.0)*(max-min)/12.0 + 0.5*(max+min);
	}

	return 1;
}
//---------------------------------------------------------------------------
//Randomly generate a seed (original point) of a CNT in the RVE
int GenNetwork::Get_seed_point(const struct cuboid &cub, int &seed, Point_3D &point)const
{
	seed = (2053*seed + 13849)%MAX_INT;
	point.x = cub.poi_min.x + seed*cub.len_x/MAX_INT;
    
	seed = (2053*seed + 13849)%MAX_INT;
	point.y = cub.poi_min.y + seed*cub.wid_y/MAX_INT;
    
	seed = (2053*seed + 13849)%MAX_INT;
	point.z = cub.poi_min.z + seed*cub.hei_z/MAX_INT;
    
	point.flag = 0; //0 denotes this point is the original point of a CNT
    
	return 1;
}
//---------------------------------------------------------------------------
//Randomly generate a direction in the spherical coordinates as the original direction of CNT segments
int GenNetwork::Get_uniform_direction(const struct Nanotube_Geo &nanotube_geo, int &seed_sita, int &seed_pha, double &cnt_sita, double &cnt_pha)const
{
	if(nanotube_geo.dir_distrib_type=="random")
	{
		//sita is chosen in [0, PI] with uniform distribution
		seed_sita = (2053*seed_sita + 13849)%MAX_INT;
		cnt_sita = seed_sita*PI/MAX_INT;
        
		//pha is chosen in [0, 2PI] with uniform distribution
		seed_pha = (2053*seed_pha + 13849)%MAX_INT;
		cnt_pha = 2.0*seed_pha*PI/MAX_INT;
	}
	else if(nanotube_geo.dir_distrib_type=="specific")
	{
		//A specific original-direction
		seed_sita = (2053*seed_sita + 13849)%MAX_INT;
		seed_pha = (2053*seed_pha + 13849)%MAX_INT;
        
		if((seed_sita+seed_pha)%2==0)
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
	//cnt_sita�Ǳ任����
	MathMatrix Msita(3,3);
	Msita.element[0][0] = cos(sita);
	Msita.element[0][2] = sin(sita);
	Msita.element[1][1] = 1;
	Msita.element[2][0] = -sin(sita);
	Msita.element[2][2] = cos(sita);
    
	//cnt_pha�Ǳ任����
	MathMatrix Mpha(3,3);
	Mpha.element[0][0] = cos(pha);
	Mpha.element[0][1] = -sin(pha);
	Mpha.element[1][0] = sin(pha);
	Mpha.element[1][1] = cos(pha);
	Mpha.element[2][2] = 1;
    
	return Mpha*Msita;
}
//---------------------------------------------------------------------------
//Randomly generate a direction in the spherical coordinates
//To have the positive Z-axis to be a central axis
//Then, the radial angle, sita, obeys a normal distribution (sita \in fabs[(-omega,+omega)]) and the zonal angle, pha, obeys a uniform distribution (pha \in (0,2PI))
int GenNetwork::Get_normal_direction(const double &omega, int &seed_sita, int &seed_pha, double &cnt_sita, double &cnt_pha)const
{
	//sita centres 0 and obeys a normal distribution in (-omega, +omega)
	int sum=0;
	for(int i=0; i<12; i++)
	{
		seed_sita = (2053*seed_sita + 13849)%MAX_INT;
		sum += seed_sita;
	}
	cnt_sita = fabs((sum*omega)/(6.0*MAX_INT)-omega);
    
	//pha satisfies a uniform distribution in (0, 2PI)
	seed_pha = (2053*seed_pha + 13849)%MAX_INT;
	cnt_pha = 2.0*seed_pha*PI/MAX_INT;

	return 1;
}
//---------------------------------------------------------------------------
//To calculate the coordinates of the new CNT point (transformation of coordinates)
Point_3D GenNetwork::Get_new_point(MathMatrix &Matrix, const double &Rad)const
{
	//1D vector
	MathMatrix Rvec(3,1);
	Rvec.element[0][0] = 0;
	Rvec.element[1][0] = 0;
	Rvec.element[2][0] = Rad;
    
	//1D vector
	MathMatrix Res(3,1);
	Res = Matrix*Rvec;
    
	Point_3D Point(Res.element[0][0], Res.element[1][0], Res.element[2][0]);
    
	return Point;
}
//---------------------------------------------------------------------------
//To judge if a point is included in a RVE
int GenNetwork::Judge_RVE_including_point(const struct cuboid &cub, const Point_3D &point)const
{
	if(point.x<cub.poi_min.x||point.x>cub.poi_min.x+cub.len_x||
		point.y<cub.poi_min.y||point.y>cub.poi_min.y+cub.wid_y||
		point.z<cub.poi_min.z||point.z>cub.poi_min.z+cub.hei_z) return 0;

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
		hout << "Error, the number of intersection points between the segement and the surfaces of RVE is" << (int)t_ratio.size() << ", less than one or more than three!" << endl;
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
    if (last_bool&&new_bool) return last_point.distance_to(new_point); //both points are inside so add the total length
    else if (last_bool&&(!new_bool))  //if the new point is outside
	{
        if(Get_intersecting_point_RVE_surface(cub, last_point, new_point, ipoi_vec)==0) return 0;
	    return last_point.distance_to(ipoi_vec[0]);
    } 
	else if (new_bool)  //if the new point is inside
	{
        if(Get_intersecting_point_RVE_surface(cub, new_point, last_point, ipoi_vec)==0) return 0;
	    return last_point.distance_to(ipoi_vec[0]);
    }
	else return 0.0;
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

			//Checking the angle between two segments (nod21��nod01)
			double cos_ang210 = (x21*x01+y21*y01+z21*z01)/(sqrt(x21*x21+y21*y21+z21*z21)*sqrt(x01*x01+y01*y01+z01*z01));

			//Judge the cos value(<= 90degree and >=0degree��cos value<=1.0 and >=0.0)
			if(cos_ang210<=1.0&&cos_ang210>=0.0)
			{
				hout << "Error: there exists at least one angle which is larger than PI/2!" << endl;
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

	return 1;
}
//===========================================================================
