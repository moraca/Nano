//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	GenNetwork.cpp
//OBJECTIVE:	To generate networks with overlapping
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "GenNetwork.h"

//Generate 3D networks with ovelapping
int GenNetwork::Generate_geometric_networks(const struct Geom_RVE &geom_rve, struct Cluster_Geo &cluster_geo)const
{
	//Generate random seed in terms of local time
    srand((unsigned int)time(NULL));

	//---------------------------------------------------------------------------
	//Generate the data for nanotube clusters limited in ellipsoid surfaces (each ellipsoid is within the RVE and is separated with each other)
	if(cluster_geo.vol_fra_criterion>0.0)
	{
		int seed_ellip_poi = rand()%RAND_MAX;
		int seed_ellip_axis = rand()%RAND_MAX;
		int seed_ellip_angle = rand()%RAND_MAX;
		if(Get_ellip_clusters(geom_rve, cluster_geo, seed_ellip_poi, seed_ellip_axis, seed_ellip_angle)==0) return 0;
//		if(Get_specific_sphere_clusters(cell_geo, clust_geo, seed_ellip_poi, seed_ellip_axis, seed_ellip_angle, ellips, 1)==0) return 0;  //生成特定位置的圆球团簇
	}

	return 1;
}
//---------------------------------------------------------------------------
//Generate a number of ellipsoids
int GenNetwork::Get_ellip_clusters(const struct Geom_RVE &cell, struct Cluster_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle)const
{
	double epsilon = 0.01;
	int times = 0;
	double ellip_volume = 0.0;
	vector<struct elliparam> ellips;  //Define the temporary vector of ellipsoids for nanotube cluster zones
	Point_3D poi_min = cell.ex_origin;  //Here we choose the extended RVE
	do
	{
		//-------------------------------------------------------------
		//Ready to generate an ellipsoid
		struct elliparam ell_temp;
		//Generate the center point of an ellipsoid
		seed_poi = (2053*seed_poi + 13849)%RAND_MAX;
		ell_temp.x=poi_min.x + seed_poi*cell.len_x/RAND_MAX;
        
		seed_poi = (2053*seed_poi + 13849)%RAND_MAX;
		ell_temp.y=poi_min.y + seed_poi*cell.wid_y/RAND_MAX;
        
		seed_poi = (2053*seed_poi + 13849)%RAND_MAX;
		ell_temp.z=poi_min.z + seed_poi*cell.hei_z/RAND_MAX;
        
		//Generate the lengths of axes of an ellipsoid
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
        
		ell_temp.a = (1+epsilon)*ell_temp.a;
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
                
				if(x1-poi_min.x<Zero||x1-poi_min.x>cell.len_x-Zero||
                   y1-poi_min.y<Zero||y1-poi_min.y>=cell.wid_y-Zero||
                   z1-poi_min.z<Zero||z1-poi_min.z>=cell.hei_z-Zero)
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
		//Calculate the volume
		ellip_volume = ellip_volume + 4*PI*ell_temp.a*ell_temp.b*ell_temp.c/(3*pow(1+epsilon, 3.0));
		clust_geo.real_volume_fraction = ellip_volume/cell.volume;
    gen_again:	;
	}while(times<=N_times&&clust_geo.real_volume_fraction<clust_geo.vol_fra_criterion);
	
	//---------------------------------------------------------------------
	//椭球收缩
	for(int i=0; i<(int)ellips.size(); i++)
	{
		ellips[i].a=ellips[i].a/(1+epsilon);
		ellips[i].b=ellips[i].b/(1+epsilon);
		ellips[i].c=ellips[i].c/(1+epsilon);
	}
	//---------------------------------------------------------------------
	//输出纳米管团簇椭球网格
	if(export_mod==2)	Export_cluster_ellipsoids_mesh(cell, ellips);
    
	//---------------------------------------------------------------------
	//输出纳米管团簇椭球数据
	if(export_mod==1||export_mod==2)	Export_cluster_ellipsoids_data(ellips, clust_geo.real_volume_fraction);
    
	//检查输出模式
	if(export_mod>2||export_mod<0) { hout <<  "生成纳米管团簇椭球程序中输出模式错误！ 请检查！" << endl; return 0; }
	else {	hout << "    纳米管团簇椭球的个数和实际体积：" << (int)ellips.size() << "  " << clust_geo.real_volume_fraction << endl; }
    
	return 1;
}