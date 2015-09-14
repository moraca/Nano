//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	GenNetwork.cpp
//OBJECTIVE:	To generate networks with overlapping
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "GenNetwork.h"

//Generate 3D networks with ovelapping
int GenNetwork::Generate_geometric_networks(const struct Geom_RVE &geom_rve, struct Cluster_Geo &clust_geo, struct Nanotube_Geo &nanotube_geo)const
{
	//Generate random seed in terms of local time
    srand((unsigned int)time(NULL));

	//---------------------------------------------------------------------------
	//Generate the data for nanotube clusters limited in ellipsoid surfaces (each ellipsoid is within the RVE and is separated with each other)
	if(clust_geo.vol_fra_criterion>0.0)
	{
		int seed_ellip_poi = rand()%RAND_MAX;
		int seed_ellip_axis = rand()%RAND_MAX;
		int seed_ellip_angle = rand()%RAND_MAX;

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
	int seed_cnt_origin = rand()%RAND_MAX;     //RAND_MAX==2^15-1 in a 16-bit computer, the range of rand() is [0, RAND_MAX]
	int seed_cnt_length = rand()%RAND_MAX;		//RAND_MAX==2^31-1 in a 32-bit computer, the range of rand() is [0, RAND_MAX]
	int seed_cnt_radius = rand()%RAND_MAX;
	int seed_cnt_sita = rand()%RAND_MAX;
	int seed_cnt_pha = rand()%RAND_MAX;
	int seed_growth_probability = rand()%RAND_MAX;

	double vol_sum = 0;  //the sum of volume of generated CNTs 
	double wt_sum = 0;   //the sum of weight of generated CNTs
	int cnt_seed_count =0; //to record the number of generated seeds of a CNT (If the growth of a CNT fails, but this seed will be used for the next growth of a CNT)
	//---------------------------------------------------------------------------
    while((nanotube_geo.criterion == "vol"&&vol_sum < nanotube_geo.real_volume)||
			 (nanotube_geo.criterion == "wt"&&wt_sum < nanotube_geo.real_weight))
	{
		//---------------------------------------------------------------------------
        //Randomly generate a seed (original point) of a CNT in the extended RVE	
		struct cuboid excub;										//generate a cuboid to represent the RVE
		excub.poi_min = geom_rve.ex_origin;
		excub.len_x = geom_rve.ex_len;
		excub.wid_y = geom_rve.ey_wid;
		excub.hei_z = geom_rve.ez_hei;
		excub.volume = excub.len_x*excub.wid_y*excub.hei_z;

		Point_3D cnt_poi;
		if(Get_seed_point(excub, seed_cnt_origin, cnt_poi)==0) return 0;

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

return 0;
	}

	return 1;
}
//---------------------------------------------------------------------------
//Generate a number of ellipsoids
int GenNetwork::Get_ellip_clusters(const struct cuboid &cub, struct Cluster_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle)const
{
	double epsilon = 0.01;						//A ratio for extending ellipsoids
	double ellip_volume = 0.0;
	vector<struct elliparam> ellips;			//Define the temporary vector of ellipsoids for nanotube cluster zones

	const int N_times=1000;					//A maximum number for generation
	int times = 0;										//Count the number of generation
	do
	{
		//-------------------------------------------------------------
		//Ready to generate an ellipsoid
		struct elliparam ell_temp;
		//Generate the center point of an ellipsoid
		seed_poi = (2053*seed_poi + 13849)%RAND_MAX;
		ell_temp.x=cub.poi_min.x + seed_poi*cub.len_x/RAND_MAX;
        
		seed_poi = (2053*seed_poi + 13849)%RAND_MAX;
		ell_temp.y=cub.poi_min.y + seed_poi*cub.wid_y/RAND_MAX;
        
		seed_poi = (2053*seed_poi + 13849)%RAND_MAX;
		ell_temp.z=cub.poi_min.z + seed_poi*cub.hei_z/RAND_MAX;
        
		//Generate the lengths of half-axes of an ellipsoid
		seed_axis = (2053*seed_axis + 13849)%RAND_MAX;
		ell_temp.a=clust_geo.amin + seed_axis*(clust_geo.amax - clust_geo.amin)/RAND_MAX;
		if(!(clust_geo.bmin==0&&clust_geo.cmin==0))
		{
			seed_axis = (2053*seed_axis + 13849)%RAND_MAX;
			ell_temp.b = clust_geo.bmin + seed_axis*(ell_temp.a - clust_geo.bmin)/RAND_MAX;
            
			seed_axis = (2053*seed_axis + 13849)%RAND_MAX;
			ell_temp.c = clust_geo.cmin + seed_axis*(ell_temp.b - clust_geo.cmin)/RAND_MAX;
		}
		else
		{
			ell_temp.b = ell_temp.a;
			ell_temp.c = ell_temp.a;
		}
        
		//Generate 9 angles: [(alpha1,beta1,gamma1),(alpha2,beta2,gamma2),(alpha3,beta3,gamma3)]  
		//between three axes of ellipsoid (a,b,c) with three coordinate axes (ox,oy,oz)
		seed_angle = (2053*seed_angle + 13849)%RAND_MAX;
		double alpha1 = seed_angle*PI/RAND_MAX;
		double beta1 = 0;
		if(alpha1>PI/2.0)
		{
			seed_angle = (2053*seed_angle + 13849)%RAND_MAX;
			beta1 = (alpha1-PI/2.0) + seed_angle*2*(PI-alpha1)/RAND_MAX;
		}
		else
		{
			seed_angle	= (2053*seed_angle + 13849)%RAND_MAX;
			beta1 = (PI/2.0-alpha1) + seed_angle*2*alpha1/RAND_MAX;
		}
        
		ell_temp.alpha1	=	cos(alpha1);																//alpha1 is chosen from (0, PI)
		ell_temp.beta1	=	cos(beta1);																//beta1 is chosen from (pi/2-r1) to (pi/2+r1)
		seed_angle	= (2053*seed_angle + 13849)%RAND_MAX;
		ell_temp.gamma1 = pow(-1.0, fmod(seed_angle, 2.0)+1.0)*sqrt(1.0-pow(ell_temp.alpha1,2)-pow(ell_temp.beta1,2));	 //Calculate the value of gamma but randomly choose "positive" or "negative"
		double alpha2 = 0;																					//alpha2 is chosen from (pi/2-r1) to (pi/2+r1)
		if(alpha1>PI/2.0)
		{
			seed_angle	= (2053*seed_angle + 13849)%RAND_MAX;
			alpha2  = (alpha1-PI/2.0) + seed_angle*2*(PI-alpha1)/RAND_MAX;
		}
		else
		{
			seed_angle	= (2053*seed_angle + 13849)%RAND_MAX;
			alpha2  = (PI/2.0-alpha1) + seed_angle*2*alpha1/RAND_MAX;
		}
		ell_temp.alpha2 = cos(alpha2);
        
		double A, B, C;
		A = 1+pow(ell_temp.beta1/ell_temp.gamma1,2);
		B = 2*(ell_temp.alpha1*ell_temp.alpha2*ell_temp.beta1)/pow(ell_temp.gamma1,2);
		C = pow(ell_temp.alpha1*ell_temp.alpha2/ell_temp.gamma1,2)+pow(ell_temp.alpha2,2)-1.0;
        
		seed_angle	= (2053*seed_angle + 13849)%RAND_MAX;
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
		ellip_volume = ellip_volume + 4*PI*ell_temp.a*ell_temp.b*ell_temp.c/(3*pow(1+epsilon, 3.0));
		clust_geo.real_volume_fraction = ellip_volume/cub.volume;
    gen_again:	;
	}while(times<=N_times&&clust_geo.real_volume_fraction<clust_geo.vol_fra_criterion);
	
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
	if(clust_geo.print_key==1||clust_geo.print_key==2)	Export_cluster_ellipsoids_data(ellips, clust_geo.real_volume_fraction);

	//To print the number of ellipsoids and volume fraction
	hout << "    The number of clusters and the sum of their volume fraction:" << (int)ellips.size() << "  " << clust_geo.real_volume_fraction << endl;

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
				seed_angle = (2053*seed_angle + 13849)%RAND_MAX;
				double alpha1 = seed_angle*PI/RAND_MAX;
				double beta1 = 0;
				if(alpha1>PI/2.0)
				{
					seed_angle = (2053*seed_angle + 13849)%RAND_MAX;
					beta1 = (alpha1-PI/2.0) + seed_angle*2*(PI-alpha1)/RAND_MAX;
				}
				else
				{
					seed_angle	= (2053*seed_angle + 13849)%RAND_MAX;
					beta1 = (PI/2.0-alpha1) + seed_angle*2*alpha1/RAND_MAX;
				}
                
				ell_temp.alpha1	=	cos(alpha1);																//alpha1 is chosen from (0, PI)
				ell_temp.beta1	=	cos(beta1);																//beta1 is chosen from (pi/2-r1) to (pi/2+r1)
				seed_angle	= (2053*seed_angle + 13849)%RAND_MAX;
				ell_temp.gamma1 = pow(-1.0, fmod(seed_angle, 2.0)+1.0)*sqrt(1.0-pow(ell_temp.alpha1,2)-pow(ell_temp.beta1,2));	  //Calculate the value of gamma but randomly choose "positive" or "negative"
				double alpha2 = 0;																					//alpha2 is chosen from (pi/2-r1) to (pi/2+r1)
				if(alpha1>PI/2.0)
				{
					seed_angle	= (2053*seed_angle + 13849)%RAND_MAX;
					alpha2  = (alpha1-PI/2.0) + seed_angle*2*(PI-alpha1)/RAND_MAX;
				}
				else
				{
					seed_angle	= (2053*seed_angle + 13849)%RAND_MAX;
					alpha2  = (PI/2.0-alpha1) + seed_angle*2*alpha1/RAND_MAX;
				}
				ell_temp.alpha2 = cos(alpha2);
                
				double A, B, C;
				A = 1+pow(ell_temp.beta1/ell_temp.gamma1,2);
				B = 2*(ell_temp.alpha1*ell_temp.alpha2*ell_temp.beta1)/pow(ell_temp.gamma1,2);
				C = pow(ell_temp.alpha1*ell_temp.alpha2/ell_temp.gamma1,2)+pow(ell_temp.alpha2,2)-1.0;
                
				seed_angle	= (2053*seed_angle + 13849)%RAND_MAX;
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
				ellip_volume = ellip_volume + 4*PI*ell_temp.a*ell_temp.b*ell_temp.c/3;
				clust_geo.real_volume_fraction = ellip_volume/cub.volume;
			}
    
	//To check if the volume fraction is less than the criterion value
	if(clust_geo.real_volume_fraction<clust_geo.vol_fra_criterion)
	{
		hout << "The sum of volume fraction of spheres: " << clust_geo.real_volume_fraction;
		hout << " which is less than the criterion value: " << clust_geo.vol_fra_criterion << " , please check it again!" << endl;
		return 0;
	}
    
	//---------------------------------------------------------------------
	//Print the ellipsoid surfaces by grids
	if(clust_geo.print_key==2)	Export_cluster_ellipsoids_mesh(cub, ellips);
    
	//---------------------------------------------------------------------
	//Export the data of ellipsoid surfaces
	if(clust_geo.print_key==1||clust_geo.print_key==2)	Export_cluster_ellipsoids_data(ellips, clust_geo.real_volume_fraction);
    
	//To print the number of ellipsoids and volume fraction
	hout << "    The number of clusters and the sum of their volume fraction:" << (int)ellips.size() << "  " << clust_geo.real_volume_fraction << endl;
    
	return 1;
}
//---------------------------------------------------------------------------
//Generate a random value through a probability distribution function
int GenNetwork::Get_random_value(const string &dist_type, const double &min, const double &max, int &seed, double &value)const
{
	if(min>max) { hout << "Error, the minimum value is larger than the maximum value (Get_random_value)!" << endl; return 0; }
    
	if(dist_type=="uniform")	//uniform distribution
	{
		seed = (2053*seed + 13849)%RAND_MAX;
		value = seed*(max-min)/RAND_MAX + min;
	}
	else if(dist_type=="normal")	//normal distribution
	{
		int sum=0;
		for(int i=0; i<12; i++)
		{
			seed = (2053*seed + 13849)%RAND_MAX;
			sum += seed;
		}
		value = ((double)sum/RAND_MAX-6.0)*(max-min)/12.0 + 0.5*(max+min);
	}

	return 1;
}
//---------------------------------------------------------------------------
//Randomly generate a seed (original point) of a CNT in the RVE
int GenNetwork::Get_seed_point(const struct cuboid &cub, int &seed, Point_3D &point)const
{
	seed = (2053*seed + 13849)%RAND_MAX;
	point.x = cub.poi_min.x + seed*cub.len_x/RAND_MAX;
    
	seed = (2053*seed + 13849)%RAND_MAX;
	point.y = cub.poi_min.y + seed*cub.wid_y/RAND_MAX;
    
	seed = (2053*seed + 13849)%RAND_MAX;
	point.z = cub.poi_min.z + seed*cub.hei_z/RAND_MAX;
    
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
		seed_sita = (2053*seed_sita + 13849)%RAND_MAX;
		cnt_sita = seed_sita*PI/RAND_MAX;
        
		//pha is chosen in [0, 2PI] with uniform distribution
		seed_pha = (2053*seed_pha + 13849)%RAND_MAX;
		cnt_pha = 2.0*seed_pha*PI/RAND_MAX;
	}
	else if(nanotube_geo.dir_distrib_type=="specific")
	{
		//A specific original-direction
		seed_sita = (2053*seed_sita + 13849)%RAND_MAX;
		seed_pha = (2053*seed_pha + 13849)%RAND_MAX;
        
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
	//cnt_sita½Ç±ä»»¾ØÕó
	MathMatrix Msita(3,3);
	Msita.element[0][0] = cos(sita);
	Msita.element[0][2] = sin(sita);
	Msita.element[1][1] = 1;
	Msita.element[2][0] = -sin(sita);
	Msita.element[2][2] = cos(sita);
    
	//cnt_pha½Ç±ä»»¾ØÕó
	MathMatrix Mpha(3,3);
	Mpha.element[0][0] = cos(pha);
	Mpha.element[0][1] = -sin(pha);
	Mpha.element[1][0] = sin(pha);
	Mpha.element[1][1] = cos(pha);
	Mpha.element[2][2] = 1;
    
	return Mpha*Msita;
}
//===========================================================================
