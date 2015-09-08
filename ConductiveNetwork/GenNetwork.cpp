//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	GenNetwork.cpp
//OBJECTIVE:	To generate networks with overlapping
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "GenNetwork.h"

//Generate 3D networks with ovelapping
int GenNetwork::Generate_geometric_networks(const Input *Init)const
{
	//Generate random seed in terms of local time
    srand((unsigned int)time(NULL));

	cout << "RAND_MAX: " << RAND_MAX << endl;
	//---------------------------------------------------------------------------
	//�������׹��Ŵ�������Ϣ�������ڵ����ڲ����໥���룩
	if(Init->cluster_geo.vol_fra_criterion>0.0)
	{
		int seed_ellip_poi = rand()%RAND_MAX;
		int seed_ellip_axis = rand()%RAND_MAX;
		int seed_ellip_angle = rand()%RAND_MAX;
		//���һλʵ�Σ�0��ʾ�������1��ʾֻ������׹��Ŵ��������ݣ�2��ʾ����Ŵ��������ݼ�����������
//		if(Get_ellip_clusters(cell_geo, clust_geo, seed_ellip_poi, seed_ellip_axis, seed_ellip_angle, ellips, 1)==0) return 0;
        //		if(Get_specific_sphere_clusters(cell_geo, clust_geo, seed_ellip_poi, seed_ellip_axis, seed_ellip_angle, ellips, 1)==0) return 0;  //�����ض�λ�õ�Բ���Ŵ�
	}

	return 1;
}
