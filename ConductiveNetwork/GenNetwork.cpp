//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	GenNetwork.cpp
//OBJECTIVE:	To generate networks with overlapping
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "GenNetwork.h"

//Generate 3D networks with ovelapping
int GenNetwork::Generate_geometric_networks(const struct Geom_RVE &geom_rve, struct Cluster_Geo &clust_geo)const
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
//===========================================================================
