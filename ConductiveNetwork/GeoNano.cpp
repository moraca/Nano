//===========================================================================
// GeoNano.cpp
// NanoComposites΢�ṹ���ν�ģ��Ա����
// Member functions in a class of cell of geometric modelling
//===========================================================================
#include "GeoNano.h"

//---------------------------------------------------------------------------
//���ɼ������ݿ�
int GeoNano::Geometric_modeling(ifstream &infile, const int &samples_count)
{
    
	//����RVE��CNT�ļ�������
	if(Import_geometric_data(infile, cell_geo, overlap_regions, cnt_regions, cnts_geo, clust_geo, samples_count, sectioned_domain, sectioned_domain_cnt, d_vdw, cutoff)==0) return 0;
    hout <<  "Import_geometric_data" << endl;
    
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�������׹�����
	vector<struct elliparam> ellips;					//�������׹��Ŵ�������Ϣ
	vector<vector<Point_3D> > cnts_points; //������ά������(���ڼ�¼���׹�����)
	if(cnts_geo.criterion=="nwt")
	{
		if(New_generate_nanotube_networks(cell_geo, cnts_geo, cnts_points, clust_geo, ellips)==0) return 0;
	}
	else if(cnts_geo.criterion=="wt"||cnts_geo.criterion=="vol")
	{
        if (cnts_geo.type =="CF") {
            hout << "Carbor Fibers not yet implemented" << endl;
            return 0;
        }
        else if (cnts_geo.type =="CNT") {
            if(Generate_nanotube_networks(cell_geo, cnts_geo, cnts_points, clust_geo, ellips)==0) return 0;
        } else {
            hout << "ERROR. Invalid material: " << cnts_geo.type << ". Only options are CF or CNT."  << endl;
            return 0;
        }
	}
	else	{ hout << "���׹�����׼��Ȳ������(vol)Ҳ��������(wt)�� ���飡" << endl; return 0; }
    
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�������׹ܵ��������
    //	if(CNTs_quality_testing(cnts_points)==0) return 0;
    
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//������׹����缸�ι���(��Tecplot)
	//if(Export_cnt_networks_threads(cell_geo, cnts_points)==0) return 0;
    
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//������׹����缸��(Template for Sofiane)
    //	if(Export_cnt_fibers(cnts_points)==0) return 0;
    
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//������׹����缸�ι���(����������Tecplot) //ע����RVE���洦��û�п��Ǳ���Ҫ�е����׹ܵ�һ����
    //if(Export_cnt_networks_meshes(cell_geo, cnts_points, ellips)==0) return 0;
    
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//����������Լ�Ϳ�㼸�ι���(����������(Ϳ��)��������Tecplot)
    //	if(Export_cnt_threads_coating_mesh(cell_geo, cnts_points)==0) return 0;
    
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//������׹��Լ�Ϳ�㼸�ι���(�����嵥Ԫ(���׹�)�����嵥Ԫ(Ϳ��)Tecplot)(Ϳ�����ڴ˺���������)
    //	if(Export_cnt_coating_meshes(cell_geo, cnts_points)==0) return 0;
    
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�������׹��Ŵ����������׹ܵı��غ��������
	//if(Estimate_cnt_cluster_parameters(cell_geo, cnts_points, ellips, cnts_geo, clust_geo)==0) return 0;
    
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//��¼���׹ܵ���Ϣ(Point��flag==0��ʾ���׹ܵ���ʼ��; flag>0��ʾ���׹��ڵĵ�, �������α��)
    //================================================================
    //Modified the flag to have the CNT number
	if(Record_cnt_points_information(cnts_points)==0) return 0;
    
	return 1;
}
//---------------------------------------------------------------------------
//����RVE��CNT�ļ�������
int GeoNano::Import_geometric_data(ifstream &infile, struct RVE_Geo &cell_geo, struct Region_Geo &overlap_regions, struct Region_Geo &cnt_regions, struct CNT_Geo &cnts_geo, struct Clust_Geo &clust_geo, const int &samples_count, vector<vector<long int> > &sectioned_domain, vector<vector<int> > &sectioned_domain_cnt, double &d_vdw, double &cutoff)const
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//����RVE���������½ǵ㣨���������������Сֵ�㣩�ͳ����ֵ
	istringstream istr_rve(Get_Line(infile));
	istr_rve >> cell_geo.poi_min.x >> cell_geo.poi_min.y >> cell_geo.poi_min.z;
	istr_rve >> cell_geo.len_x >> cell_geo.wid_y >> cell_geo.hei_z;
	if(cell_geo.len_x<0||cell_geo.wid_y<0||cell_geo.hei_z<0){	hout << "RVE dimensions are negative" << endl;	return 0; }
	cell_geo.volume = cell_geo.len_x*cell_geo.wid_y*cell_geo.hei_z;		//RVE���
    
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Get the type of the network (CNT or CF)
    istringstream istr_type(Get_Line(infile));
    istr_type >> cnts_geo.type;
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�������׹ܳ�ʼ��������
	istringstream istr_initial_direction(Get_Line(infile));
	istr_initial_direction >> cnts_geo.dir_dist_type;
	if(cnts_geo.dir_dist_type!="random"&&cnts_geo.dir_dist_type!="specific"){ hout << "The direction distribution type must be either random or specific" << endl;	return 0; }
	if(cnts_geo.dir_dist_type=="specific")
	{
		istr_initial_direction >> cnts_geo.ini_sita >> cnts_geo.ini_pha;
		if(cnts_geo.ini_sita<0||cnts_geo.ini_sita>PI||cnts_geo.ini_pha<0||cnts_geo.ini_pha>=2*PI)
		{
			hout << "The specified angle is not in the acceptable range of 0-2PI" << endl;
			return 0;
		}
	}
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�������׹�ÿ������ʱ�Ƕȷ�����̬�ֲ������ȡֵ��Χ
	istringstream istr_angle_range(Get_Line(infile));
	istr_angle_range >> cnts_geo.angle_max;
	if(cnts_geo.angle_max>0.5*PI){ hout << "The specified angle is not in the acceptable range" << endl;	return 0; }
    
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�������׹���������
	istringstream istr_step_len(Get_Line(infile));
	istr_step_len >> cnts_geo.step_length;
	if(cnts_geo.step_length<=0||cnts_geo.step_length>=0.25*cell_geo.len_x||cnts_geo.step_length>=0.25*cell_geo.wid_y||
       cnts_geo.step_length>=0.25*cell_geo.hei_z){ hout << "The step length must be positive and 0.25 times lesser than the dimension of the box" << endl;	return 0; }
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Get the distribution type of the length of the whole CNT/CF and the minimum and maximum range
    istringstream istr_cnt_len(Get_Line(infile));
    istr_cnt_len >> cnts_geo.len_dist_type;
    if(cnts_geo.len_dist_type!="uniform"&&cnts_geo.len_dist_type!="normal"){ hout << "The distribution of the length should be either normal or uniform" << endl;	return 0; }
    istr_cnt_len >> cnts_geo.len_min >> cnts_geo.len_max;
    if(cnts_geo.len_min<0||cnts_geo.len_max<0||cnts_geo.len_max<cnts_geo.len_min){ hout << "The length must be non-negative and min must be smaller than max" << endl; return 0; }
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //Get the distribution type for the radius of the CNT/CF and the minimum and maximum range
    istringstream istr_cnt_rad(Get_Line(infile));
    istr_cnt_rad >> cnts_geo.rad_dist_type;
    if(cnts_geo.rad_dist_type!="uniform"&&cnts_geo.rad_dist_type!="normal"){ hout << "The distribution of the radius should be either normal or uniform" << endl;	return 0; }
    istr_cnt_rad >> cnts_geo.rad_min >> cnts_geo.rad_max;
    if(cnts_geo.rad_min<0||cnts_geo.rad_max<0||cnts_geo.rad_max<cnts_geo.rad_min){ hout << "The radius must be non-negative and min must be smaller than max" << endl; return 0; }//||cnts_geo.rad_max>0.05*cnts_geo.len_min
    //-----------------------------------------------------------------------------------------------------------------------------------------
	//The final criterion for the geometry generation (wt or vol)
	istringstream istr_cnt_vol(Get_Line(infile));
	istr_cnt_vol >> cnts_geo.criterion;
	if(cnts_geo.criterion=="vol")
	{
		istr_cnt_vol >> cnts_geo.volume_fraction;
		if(cnts_geo.volume_fraction>1||cnts_geo.volume_fraction<0){ hout << "The volume fraction must be between 0 and 1" << endl; return 0; }
		hout << "    The volume fraction is "<< cnts_geo.volume_fraction << endl;
        
		int mode_accum;
		istr_cnt_vol >> mode_accum;
		if(mode_accum<0&&mode_accum>2){ hout <<"The mode of accumulation should be between 0 and 2" << endl; return 0; }
		if(mode_accum==1)	 cnts_geo.volume_fraction = cnts_geo.volume_fraction*samples_count;
		else if(mode_accum==2)	 cnts_geo.volume_fraction = cnts_geo.volume_fraction*pow(2.0, samples_count-1);
		hout <<"    The volume fraction is "  << cnts_geo.volume_fraction << endl;
        
		//The total volume of the CNT/CF network
		cnts_geo.real_volume = cnts_geo.volume_fraction*cell_geo.volume;
	}
	else if(cnts_geo.criterion=="wt"||cnts_geo.criterion=="nwt")
	{
		istr_cnt_vol >> cnts_geo.weight_fraction;
		if(cnts_geo.weight_fraction>1||cnts_geo.weight_fraction<0){ hout << "���׹���RVE����ռ���ط���������� ���飡" << endl; return 0; }
		hout << "    ��������ֵ�ǣ� " << cnts_geo.weight_fraction << endl;
        
		int mode_accum;
		istr_cnt_vol >> mode_accum;
		if(mode_accum<0&&mode_accum>2){ hout << "���������ۼƷ�ʽ���� ���飡" << endl; return 0; }
		if(mode_accum==1)	 cnts_geo.weight_fraction = cnts_geo.weight_fraction*samples_count;
		else if(mode_accum==2)	 cnts_geo.weight_fraction = cnts_geo.weight_fraction*pow(2.0, samples_count-1);
		hout << "    ʵ�ʱ���ֵ�ǣ� " << cnts_geo.weight_fraction << endl;
        
		istr_cnt_vol >> cnts_geo.linear_density;  //�������׹ܵ����ܶ�
		if(cnts_geo.linear_density<0){ hout << "���׹ܵ����ܶ�������� ���飡" << endl; return 0; }
		istr_cnt_vol >> cell_geo.density;				//�������(polymer)���ܶȣ��������ں������׹ܵ���������Ի�����ܶȽ���Ҳ�ǵ������ܶ�
		if(cell_geo.density<0){ hout << "��������ܶ�������� ���飡" << endl; return 0; }
		if(cnts_geo.linear_density>=cell_geo.density){ hout << "���׹ܵ����ܶȻ��߻�������ܶ�������� ���飡" << endl; return 0; }
        
		//ʵ�����׹ܱ���(ʵ�������׹ܵ����ܶ��ǿ�����һ����Χ�ڱ仯�ģ���ͬ�뾶)
		cnts_geo.real_weight = cnts_geo.weight_fraction*cell_geo.volume*cell_geo.density;
	}
	else { hout << "���׹�����׼��Ȳ������(vol)Ҳ��������(wt)��������� ���飡" << endl; return 0; }
	
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�������׹��Ŵ��������
	istringstream istr_clust_para(Get_Line(infile));
	istr_clust_para >> clust_geo.vol_fra_criterion;
	if(clust_geo.vol_fra_criterion>1||clust_geo.vol_fra_criterion<0){ hout << "���׹��Ŵ�������RVE����ռ�������������� ���飡" << endl; return 0; }
    clust_geo.vol_fra_criterion = clust_geo.vol_fra_criterion*pow(2.0, samples_count-1);
	if(clust_geo.vol_fra_criterion!=0)  //���������׹��Ŵ�ʱ�����������Ŵ�������Ϣ
	{
		istr_clust_para >> clust_geo.amin >> clust_geo.amax;
		if(clust_geo.amin<0||clust_geo.amax<0||clust_geo.amin>clust_geo.amax){ hout << "���׹��Ŵ�����ĳ���ȡֵ��Χ������� ���飡" << clust_geo.amin << "  " << clust_geo.amax << endl; return 0; }
		istr_clust_para >> clust_geo.bmin >> clust_geo.cmin;
		if(clust_geo.bmin<0||clust_geo.bmin>clust_geo.amin||clust_geo.cmin<0||clust_geo.cmin>clust_geo.bmin||
           clust_geo.cmin>clust_geo.amin){ hout << "���׹��Ŵ��������������ȡֵ��Χ������� ���飡" << endl; return 0; }
		istr_clust_para >> clust_geo.growth_probability;
		if(clust_geo.growth_probability<0||clust_geo.growth_probability>1){ hout << "���׹����Ŵ������������ĸ����������(<0����>1)�� ���飡" << endl; return 0; }
		istr_clust_para >> clust_geo.wt_fra_cluster;
		if(clust_geo.wt_fra_cluster<0||clust_geo.wt_fra_cluster>1) { hout << "���׹����Ŵ����������׹ܵ������������(<0����>1)�� ���飡" << endl; return 0; }
        
	}
	//-----------------------------------------------------------------------------------------------------------------------------------------
    //------------- AMC
    // Read the parameters for the number of regions in which the sample will be divided.
    // This will be used to determine penetration of CNTs.
    istringstream istr_regions(Get_Line(infile));
    istr_regions >> overlap_regions.secx >> overlap_regions.secy >> overlap_regions.secz;
    
    //There will be secx*secy*secz different regions
    //hout << (size_t)overlap_regions.secx*overlap_regions.secy*overlap_regions.secz << endl;
    //sectioned_domain.resize((size_t)overlap_regions.secx*overlap_regions.secy*overlap_regions.secz);
    vector<long int> empty;
    //hout << "empty" << endl;
    sectioned_domain.assign((size_t)overlap_regions.secx*overlap_regions.secy*overlap_regions.secz, empty);
    hout << "resize1 ";
    
    //Calculate the sizes of each region and save them in the variable for the region geometry
    overlap_regions.lx = cell_geo.len_x/overlap_regions.secx;
    overlap_regions.ly = cell_geo.wid_y/overlap_regions.secy;
    overlap_regions.lz = cell_geo.hei_z/overlap_regions.secz;
    
    //Now calculate the size of the regions for locating the CNTs
    if (cnts_geo.len_max < cell_geo.len_x) {
        cnt_regions.secx = (int)(cell_geo.len_x/cnts_geo.len_max)+1;
        cnt_regions.lx = cell_geo.len_x/cnt_regions.secx;
    } else {
        cnt_regions.secx = 5;
        cnt_regions.lx = cell_geo.len_x/cnt_regions.secx;
    }
    
    if (cnts_geo.len_max < cell_geo.wid_y) {
        cnt_regions.secy = (int)(cell_geo.wid_y/cnts_geo.len_max)+1;
        cnt_regions.ly = cell_geo.wid_y/cnt_regions.secy;
    } else {
        cnt_regions.secy = 5;
        cnt_regions.ly = cell_geo.wid_y/cnt_regions.secy;
    }
    
    if (cnts_geo.len_max < cell_geo.hei_z) {
        cnt_regions.secz = (int)(cell_geo.hei_z/cnts_geo.len_max)+1;
        cnt_regions.lz = cell_geo.hei_z/cnt_regions.secz;
    } else {
        cnt_regions.secz = 5;
        cnt_regions.lz = cell_geo.hei_z/cnt_regions.secz;
    }
    
    //Resize the vector for locating the CNTs
    //sectioned_domain_cnt.resize((long int)cnt_regions.secx*cnt_regions.secy*cnt_regions.secz);
    vector<int> empty_int;
    sectioned_domain_cnt.assign((size_t)cnt_regions.secx*cnt_regions.secy*cnt_regions.secz, empty_int);
    hout << "resize2" << endl;
    
    // Read the van der Waals distance
    istringstream istr_vdw(Get_Line(infile));
    istr_vdw >> d_vdw;

    
    // Check that the regions are not too small for the maximum cutoff distance 2r_max+d_vdw
    cutoff = 2*cnts_geo.rad_max + d_vdw;
    if (overlap_regions.lx < 2*cutoff) {
        hout << "The regions along x are too many for the cutoff distance for tunneling." <<  endl;
        hout << "The length of the region along each direction has to be at least twice the cutoff distance" << endl;
        hout << "2cutoff = 2(2r_max+d_vdw) = " << 2*cutoff << ", size of region along x = " << overlap_regions.lx << endl;
        return 0;
    }
    if (overlap_regions.ly < 2*cutoff) {
        hout << "The regions along y are too many for the cutoff distance for tunneling." <<  endl;
        hout << "The length of the region along each direction has to be at least twice the cutoff distance" << endl;
        hout << "2cutoff = 2(2r_max+d_vdw) = " << 2*cutoff << ", size of region along y = " << overlap_regions.ly << endl;
        return 0;
    }
    if (overlap_regions.lz < 2*cutoff) {
        hout << "The regions along z are too many for the cutoff distance for tunneling." <<  endl;
        hout << "The length of the region along each direction has to be at least twice the cutoff distance" << endl;
        hout << "2cutoff = 2(2r_max+d_vdw) = " << 2*cutoff << ", size of region along z = " << overlap_regions.lz << endl;
        return 0;
    }
    hout << "Import data end" << endl;

    //------------- AMC
	
	return 1;
}
//---------------------------------------------------------------------------
//���㷨�������׹�����
int GeoNano::New_generate_nanotube_networks(const struct RVE_Geo &cell_geo, const struct CNT_Geo &cnts_geo, vector<vector<Point_3D> > &cnts_points, struct Clust_Geo &clust_geo, vector<struct elliparam> &ellips)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�������׹�����
	//����������ɵ���ʼʱ��
    //	srand((unsigned int)time(NULL));
	
	//---------------------------------------------------------------------------
	//�������׹��Ŵ�������Ϣ�������ڵ����ڲ����໥���룩
	if(clust_geo.vol_fra_criterion>0.0)
	{
		int seed_ellip_poi = rand()%MAX_INT;
		int seed_ellip_axis = rand()%MAX_INT;
		int seed_ellip_angle = rand()%MAX_INT;
		//���һλʵ�Σ�0��ʾ�������1��ʾֻ������׹��Ŵ��������ݣ�2��ʾ����Ŵ��������ݼ�����������
		if(Get_ellip_clusters(cell_geo, clust_geo, seed_ellip_poi, seed_ellip_axis, seed_ellip_angle, ellips, 1)==0) return 0;
        //		if(Get_specific_sphere_clusters(cell_geo, clust_geo, seed_ellip_poi, seed_ellip_axis, seed_ellip_angle, ellips, 1)==0) return 0;  //�����ض�λ�õ�Բ���Ŵ�
	}
    
	//---------------------------------------------------------------------------
	//�������׹���ÿ���Ŵؿ����еĳ����Լ���������ֲ������׹ܳ���
	double cnt_weight_random = cnts_geo.real_weight;
	vector<double> cnt_weight_cluster(ellips.size());
	for(int i=0; i<(int)cnt_weight_cluster.size(); i++)
	{
		cnt_weight_cluster[i] = clust_geo.wt_fra_cluster*cell_geo.density*4*PI*ellips[i].a*ellips[i].b*ellips[i].c/3.0;
		cnt_weight_random -= cnt_weight_cluster[i];
	}
	if(cnt_weight_random<=0) { hout << "����! ���׹��ڵ����е�������С�����׹����Ŵؿ������������ܺͣ� ���飡" << endl;
        hout << "if(cnt_weight_random<=0) \n" << endl;return 0; }
    
	//---------------------------------------------------------------------------
	//������������
	int seed_cnt_origin = rand()%MAX_INT;     //���RAND_MAX==2^15, ��ôrand()ȡ[0, RAND_MAX]��ע��16λ����
	int seed_cnt_length = rand()%MAX_INT;		//���RAND_MAX==2^31, ��ôrand()ȡ[0, MAX_INT]��ע��32λ����, ��ȡ�������˵ [0, MAX_INT]��Χ�㹻��
	int seed_cnt_radius = rand()%MAX_INT;
	int seed_cnt_sita = rand()%MAX_INT;
	int seed_cnt_pha = rand()%MAX_INT;
	int seed_growth_probability = rand()%MAX_INT;
    
	double wt_total = 0;
	double vol_total = 0;
	//---------------------------------------------------------------------------
	//��������ֲ��������׹�
	double wt_sum = 0;
	double vol_sum = 0;
	const double one_step_weight = cnts_geo.step_length*cnts_geo.linear_density;
	while(wt_sum<cnt_weight_random&&one_step_weight<cnt_weight_random)
	{
		//---------------------------------------------------------------------------
		//���嵥�����׹�����
		vector<Point_3D> new_cnt;
		int new_cnt_size = (int)new_cnt.size();
        
		//---------------------------------------------------------------------------
		//�ֱ���RVE�ĳ�����ϰ����ȷֲ�ѡȡ���׹ܵ���ʼ��
		Point_3D cnt_poi;
		if(Get_seed_point_outside_clusters(cell_geo, seed_cnt_origin, cnt_poi, ellips)==0) return 0;
		new_cnt.push_back(cnt_poi);	//�洢�ڵ�
        
		//---------------------------------------------------------------------------
		//���ֲ��������ѡȡ���׹ܳ���
		double cnt_length;
		if(Get_random_value(cnts_geo.len_dist_type, cnts_geo.len_min, cnts_geo.len_max, seed_cnt_length, cnt_length)==0) return 0;
		//��������������
		int step_num = (int)(cnt_length/cnts_geo.step_length) + 1;
        
		//---------------------------------------------------------------------------
		//���ֲ��������ѡȡ���׹ܰ뾶
		double cnt_rad;
		if(Get_random_value(cnts_geo.rad_dist_type, cnts_geo.rad_min, cnts_geo.rad_max, seed_cnt_radius, cnt_rad)==0) return 0;
        
		//---------------------------------------------------------------------------
		//�����������ϰ����ȷֲ����ѡȡһ��������Ϊ���׹��ߵ��׸�����
		double cnt_sita, cnt_pha;
		if(Get_uniform_direction(cnts_geo, seed_cnt_sita, seed_cnt_pha, cnt_sita, cnt_pha)==0) return 0;
		MathMatrix multiplier(3,3);
		multiplier = Get_transformation_matrix(cnt_sita, cnt_pha);
        
		//---------------------------------------------------------------------------
		//���׹�ÿ����һ�������ӵ����ϵ��(����û�п������������������������)
		const double step_vol_para = PI*cnt_rad*cnt_rad;
		//---------------------------------------------------------------------------
		//���׹�ÿ����һ�������ӵ�����ϵ��(���������״�ϸ��ͬʱ�����ܶ�Ҳ����һ����Χ�ڿɱ�ģ�ÿ�����׹ܿ��ܲ�ͬ)
		const double step_wei_para = cnts_geo.linear_density;
        
		//---------------------------------------------------------------------------
		//���׹�����
		for(int i=0; i<step_num; i++)
		{
			//����Z������һ�µķ���Ϊ����, ������̬�ֲ�sita = fabs[(-omega,+omega)], ͨ��omega<=PI/2, γ����ȷֲ�pha=(0,2PI)
			//���ѡȡ����
			if(Get_normal_direction(cnts_geo.angle_max, seed_cnt_sita, seed_cnt_pha, cnt_sita, cnt_pha)==0) return 0;
			
			//�޸�����任����
			multiplier = multiplier*Get_transformation_matrix(cnt_sita, cnt_pha);
			
			//�����½ڵ�����(����任)
			cnt_poi = cnt_poi + Get_new_point(multiplier, cnts_geo.step_length);
			cnt_poi.flag = 1;							//����ʼ�˵�
            
			//-------------------------------------------------------------
			//�ж��Ƿ���������
			int ei0 = Judge_ellipses_including_point(ellips, cnt_poi);
			if(ei0>=0)
			{
				//�������׹��߶���������Ľ���(�Ҳ������㱨��)
				Point_3D intersect_point(0, 0 ,0);
				if(Get_intersection_cnt_ellipsoid(ellips[ei0], new_cnt.back(), cnt_poi, intersect_point)==0) { hout << "����û���ҵ����׹��ߺ�������Ľ��㣬���飡" << endl; return 0; }
				if(intersect_point.flag==-1)  { hout << new_cnt_size << "  " << "�������׹��ϵĵ����Ŵ��������ϵĽ�����������" << endl; return 0; }
				
				new_cnt.push_back(intersect_point);		//�洢�ڵ�
                
				//---------------------------------------------------------------------------
				//�����������������
				for(int j=new_cnt_size; j<(int)new_cnt.size()-1; j++)
				{
					if(new_cnt[j+1].flag!=0)
					{
						double temp_length = new_cnt[j].distance_to(new_cnt[j+1]);
						vol_sum += temp_length*step_vol_para;		//�������
						wt_sum += temp_length*step_wei_para;		//��������
					}
				}
				new_cnt_size = (int)new_cnt.size()-1;			//�����²�����λ��
                
				break;
			}
            
			//---------------------------------------------------------------------------
			//����������ڵ������⣬ƽ�����׹�����
			if(Judge_cell_including_point(cell_geo, cnt_poi)==0)
			{
				//�����������˵��뵥����������н���(�������˵������֮�䣬����0<t<1)����������t��С����������Щ����
				vector<Point_3D> ipoi_vec;  //��������
				if(Get_intersecting_point_RVE_surface(cell_geo, new_cnt.back(), cnt_poi, ipoi_vec)==0) return 0;
				if(new_cnt.back()!=ipoi_vec[0]) new_cnt.push_back(ipoi_vec[0]);
				for(int j=0; j<(int)ipoi_vec.size(); j++)
				{
					Point_3D inter_point[2];
                    
					inter_point[0] = ipoi_vec[j];
					if(j==(int)ipoi_vec.size()-1) inter_point[1] = cnt_poi;
					else inter_point[1] = ipoi_vec[j+1];
                    
					//����ƽ������任
					if(Periodical_coordinate_transformation(cell_geo, inter_point)==0) return 0;
                    
					if(j==(int)ipoi_vec.size()-1) cnt_poi = inter_point[1];  //����λ��
                    
					new_cnt.push_back(inter_point[0]);
					new_cnt.push_back(inter_point[1]);
				}
                
				//-------------------------------------------------------------
				//�ж��Ƿ���������
				int ei1 = Judge_ellipses_including_point(ellips, new_cnt.back());
				if(ei1>=0)
				{
					//�������׹��߶���������Ľ���(�Ҳ������㱨��)
					Point_3D intersect_point(0, 0 ,0);
					if(Get_intersection_cnt_ellipsoid(ellips[ei1], new_cnt[(int)new_cnt.size()-2], new_cnt.back(), intersect_point)==0) { hout << "����û���ҵ����׹��ߺ�������Ľ��㣬���飡" << endl; return 0; }
					if(intersect_point.flag==-1)  { hout << new_cnt_size << "  " << "�������׹��ϵĵ����Ŵ��������ϵĽ�����������" << endl; return 0; }
                    
					new_cnt.push_back(intersect_point);		//�洢�ڵ�
                    
					//---------------------------------------------------------------------------
					//�����������������
					for(int j=new_cnt_size; j<(int)new_cnt.size()-1; j++)
					{
						if(new_cnt[j+1].flag!=0)
						{
							double temp_length = new_cnt[j].distance_to(new_cnt[j+1]);
							vol_sum += temp_length*step_vol_para;		//�������
							wt_sum += temp_length*step_wei_para;		//��������
						}
					}
					new_cnt_size = (int)new_cnt.size()-1;			//�����²�����λ��
                    
					break;
				}
			}
			else
			{
				new_cnt.push_back(cnt_poi);		//�洢�ڵ�
			}
            
			//---------------------------------------------------------------------------
			//�����������������
			for(int j=new_cnt_size; j<(int)new_cnt.size()-1; j++)
			{
				if(new_cnt[j+1].flag!=0)
				{
					double temp_length = new_cnt[j].distance_to(new_cnt[j+1]);
					vol_sum += temp_length*step_vol_para;		//�������
					wt_sum += temp_length*step_wei_para;		//��������
				}
			}
			new_cnt_size = (int)new_cnt.size()-1;			//�����²�����λ��
            
			//---------------------------------------------------------------------------
			//�ж����������������
			if(wt_sum >= cnt_weight_random) break;			//��������, ����
		}
        
		//---------------------------------------------------------------------------
		//�洢���׹ܽڵ�����
		if((int)new_cnt.size()==1)
		{
			hout << "new_cnt_size: " <<  new_cnt.size() << endl;
			hout << "��������ֲ����ֵ����׹��������� ���飡" << endl;
			return 0;
		}
        
		vector<Point_3D> cnt_temp;
		for(int i=0; i<(int)new_cnt.size(); i++)
		{
			cnt_temp.push_back(new_cnt[i]); //���Ȳ���ڵ�
			//�ж��Ƿ��ǽ���
			if(i==(int)new_cnt.size()-1||new_cnt[i+1].flag==0) //ĩβ���������ת������һ�����׹�
			{
				if((int)cnt_temp.size()>1)						//����1�����˵��ֻ��һ����ʼ�㣬����������������׹��߶�
				{
					cnts_points.push_back(cnt_temp);		//��¼���׹���Ϣ
					cnts_radius.push_back(cnt_rad);			//��¼���׹ܰ뾶
				}
				cnt_temp.clear();									//�����ʱ���׹�
			}
		}
	}
	//��¼���׹����������
	vol_total += vol_sum;
	wt_total += wt_sum;
    
	//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	//ѭ�����׹��Ŵ�����, ���������Ŵ����������׹�
	for(int i=0; i<(int)ellips.size(); i++)
	{
		wt_sum = 0;  //�������
		vol_sum = 0; //������
		while(wt_sum < cnt_weight_cluster[i])
		{
			//---------------------------------------------------------------------------
			//���嵥�����׹�����
			vector<Point_3D> new_cnt;
			int new_cnt_size = (int)new_cnt.size();
            
			//---------------------------------------------------------------------------
			//�ֱ���������ĳ��ж���2��Ϊ����ߵ��������ϰ����ȷֲ�ѡȡ���׹ܵ���ʼ��
			//���жϴ˵��Ƿ��������ڣ�����������ڣ��ٱ任����������ϵ��
			Point_3D cnt_poi;
			if(Get_seed_point_inside_clusters(ellips[i], seed_cnt_origin, cnt_poi)==0) return 0;
			new_cnt.push_back(cnt_poi);	//�洢�ڵ�
            
			//---------------------------------------------------------------------------
			//���ֲ��������ѡȡ���׹ܳ���
			double cnt_length;
			if(Get_random_value(cnts_geo.len_dist_type, cnts_geo.len_min, cnts_geo.len_max, seed_cnt_length, cnt_length)==0) return 0;
			//��������������
			int step_num = (int)(cnt_length/cnts_geo.step_length) + 1;
            
			//---------------------------------------------------------------------------
			//���ֲ��������ѡȡ���׹ܰ뾶
			double cnt_rad;
			if(Get_random_value(cnts_geo.rad_dist_type, cnts_geo.rad_min, cnts_geo.rad_max, seed_cnt_radius, cnt_rad)==0) return 0;
            
			//---------------------------------------------------------------------------
			//�����������ϰ����ȷֲ����ѡȡһ��������Ϊ���׹��ߵ��׸�����
			double cnt_sita, cnt_pha;
			if(Get_uniform_direction(cnts_geo, seed_cnt_sita, seed_cnt_pha, cnt_sita, cnt_pha)==0) return 0;
			MathMatrix multiplier(3,3);
			multiplier = Get_transformation_matrix(cnt_sita, cnt_pha);
            
			//---------------------------------------------------------------------------
			//���׹�ÿ����һ�������ӵ����(����û�п������������������������)
			const double step_vol_para = PI*cnt_rad*cnt_rad;
			//---------------------------------------------------------------------------
			//���׹�ÿ����һ�������ӵ�����(���������״�ϸ��ͬʱ�����ܶ�Ҳ����һ����Χ�ڿɱ�ģ�ÿ�����׹ܿ��ܲ�ͬ)
			const double step_wei_para = cnts_geo.linear_density;
            
			//---------------------------------------------------------------------------
			//���׹�����
			for(int j=0; j<step_num; j++)
			{
				//����Z������һ�µķ���Ϊ����, ������̬�ֲ�sita = fabs[(-omega,+omega)], ͨ��omega<=PI/2, γ����ȷֲ�pha=(0,2PI)
				//���ѡȡ����
				if(Get_normal_direction(cnts_geo.angle_max, seed_cnt_sita, seed_cnt_pha, cnt_sita, cnt_pha)==0) return 0;
                
				//�޸�����任����
				multiplier = multiplier*Get_transformation_matrix(cnt_sita, cnt_pha);
                
				//�����½ڵ�����(����任)
				cnt_poi = cnt_poi + Get_new_point(multiplier, cnts_geo.step_length);
				cnt_poi.flag = 1;							//����ʼ�˵�
                
				//---------------------------------------------------------------------------
				//������׹�������Խ�Ŵؽ��棬�������׹ܰ�һ���ĸ����ۻ�
				while(Cnt_go_through_a_cluster(clust_geo.growth_probability, ellips[i], seed_growth_probability, new_cnt.back(), cnt_poi)) //��Ҫ�ۻ�
				{
					//�������׹��߶���������Ľ���(�Ҳ������㱨��)
					Point_3D intersect_point(0, 0 ,0);
					if(Get_intersection_cnt_ellipsoid(ellips[i], new_cnt.back(), cnt_poi, intersect_point)==0) { hout << "����û���ҵ����׹��ߺ�������Ľ��㣬���飡" << endl; return 0; }
					if(intersect_point.flag==-1)  break;  //�ж���������������߶ε�һ�������������ϣ�����ʵ�����������ⲿ
                    
					//�ƶ����׹��߶�����������˵㵽�������ڵĶԳƵ㣨������Ϊ�Գ��棩ע�⣺������һ��cnt_poi��λ���Ѿ������˱仯
					if(Move_cnt_point_symmetry_point(ellips[i], intersect_point, cnt_poi)==0) { hout << "����û���ҵ����׹��ߺ�������Ľ��㣬���飡" << endl; return 0; }
                    
					//��������н�
					if(Cos_angle_in_three_points(new_cnt.back(), intersect_point, cnt_poi)<cos(PI/6)) //�нǴ���30��
					{
						//���뽻�㵽���׹����У������flagֵ�ڼ���ʱ�Ѿ����裩
						new_cnt.push_back(intersect_point);
						
						//�޸�����任����
						multiplier = Get_vector_transformation_matrix(intersect_point, cnt_poi);
					}
					else //ȷ��û��̫��Ľ�
					{
						//�޸�����任����
						multiplier = Get_vector_transformation_matrix(new_cnt.back(), cnt_poi);
					}
				}
                
				new_cnt.push_back(cnt_poi);		//�洢�ڵ�
                
				//---------------------------------------------------------------------------
				//�����������������
				for(int k=new_cnt_size; k<(int)new_cnt.size()-1; k++)
				{
					if(new_cnt[k+1].flag!=0)
					{
						double temp_length = new_cnt[k].distance_to(new_cnt[k+1]);
						vol_sum += temp_length*step_vol_para;		//�������
						wt_sum += temp_length*step_wei_para;		//��������
					}
				}
				new_cnt_size = (int)new_cnt.size()-1;			//�����²�����λ��
                
				//---------------------------------------------------------------------------
				//�ж����������������
				if(wt_sum >= cnt_weight_cluster[i]) break;			//��������, ����
			}
            
			//---------------------------------------------------------------------------
			//�洢���׹ܽڵ�����
			if((int)new_cnt.size()==1)
			{
				hout << "�����Ŵ�����" << i << "�ڲ������׹��������� ���飡" << endl;
				return 0;
			}
            
			vector<Point_3D> cnt_temp;
			for(int j=0; j<(int)new_cnt.size(); j++)
			{
				cnt_temp.push_back(new_cnt[j]); //���Ȳ���ڵ�
				//�ж��Ƿ��ǽ���
				if(j==(int)new_cnt.size()-1||new_cnt[j+1].flag==0) //ĩβ���������ת������һ�����׹�
				{
					if((int)cnt_temp.size()>1)						//����1�����˵��ֻ��һ����ʼ�㣬����������������׹��߶�
					{
						cnts_points.push_back(cnt_temp);		//��¼���׹���Ϣ
						cnts_radius.push_back(cnt_rad);			//��¼���׹ܰ뾶
					}
					cnt_temp.clear();									//�����ʱ���׹�
				}
			}
		}
		//��¼���׹����������
		vol_total += vol_sum;
		wt_total += wt_sum;
	}
    
	hout << "    ���������׹ܵ��������" << vol_total << endl;
	hout << "    ���������׹ܵ���������" << wt_sum << endl;
    
	return 1;
}
//---------------------------------------------------------------------------
//Code for generating the nanotube network when the criterion is "wt" or "vol". now modified to generate carbon fiber strands
//with waviness and misalignment.
int GeoNano::Generate_nanotube_networks(const struct RVE_Geo &cell_geo, const struct CNT_Geo &cnts_geo, vector<vector<Point_3D> > &cnts_points, struct Clust_Geo &clust_geo, vector<struct elliparam> &ellips)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�������׹�����
	//����������ɵ���ʼʱ��
	srand((unsigned int)time(NULL));
	//---------------------------------------------------------------------------
	//�������׹��Ŵ�������Ϣ�������ڵ����ڲ����໥���룩
	if(clust_geo.vol_fra_criterion>0.0)
	{
		int seed_ellip_poi = rand()%MAX_INT;
		int seed_ellip_axis = rand()%MAX_INT;
		int seed_ellip_angle = rand()%MAX_INT;
		//���һλʵ�Σ�0��ʾ�������1��ʾֻ������׹��Ŵ��������ݣ�2��ʾ����Ŵ��������ݼ�����������
		if(Get_ellip_clusters(cell_geo, clust_geo, seed_ellip_poi, seed_ellip_axis, seed_ellip_angle, ellips, 1)==0) return 0;
	}
    
	//---------------------------------------------------------------------------
	int seed_cnt_origin = rand()%MAX_INT;     //���RAND_MAX==2^15, ��ôrand()ȡ[0, RAND_MAX]��ע��16λ����
	int seed_cnt_length = rand()%MAX_INT;		//���RAND_MAX==2^31, ��ôrand()ȡ[0, MAX_INT]��ע��32λ����, ��ȡ�������˵ [0, MAX_INT]��Χ�㹻��
	int seed_cnt_radius = rand()%MAX_INT;
	int seed_cnt_sita = rand()%MAX_INT;
	int seed_cnt_pha = rand()%MAX_INT;
	int seed_growth_probability = rand()%MAX_INT;
	
	double vol_sum = 0;  //��¼���׹�����ռ�����
	double wt_sum = 0;   //��¼���׹�����ռ������
	int cnt_seed_count =0; //��¼���˶��ٴ����׹����ӣ����л�������ʧ�ܵģ��������������ڣ�
    
    //----------- AMC
    //Global variable for the regions.
    //Each element of the vector saves the CNT number and point number within the CNT
    vector<long int> empty;
    //current_region = stores the region of the point
    long int current_region;
    //counter is the limit for attempts for relocating the seed
    int counter;
    //boundary_flag is a flag to determine if a point is now in the boundary after it was moved
    // if 1, then the point is at the boundary and the CNT does not grow more
    int boundary_flag;
    //delete_flag is a flag to determine if a CNT can be incorporated to the rest of CNTs
    // if 1 then the CNT penetration could not be avoided so it has to be deleted
    int delete_flag = 0;
    //This is just to keep the count of rejected CNTs
    int rejects = 0;
    //To count the points that penetrate other points
    penetrating_points = 0;

    while((cnts_geo.criterion == "vol"&&vol_sum < cnts_geo.real_volume)||
          (cnts_geo.criterion == "wt"&&wt_sum < cnts_geo.real_weight))
	{
		//---------------------------------------------------------------------------
		//���ֲ��������ѡȡ���׹ܳ���
		double cnt_length;
		if(Get_random_value(cnts_geo.len_dist_type, cnts_geo.len_min, cnts_geo.len_max, seed_cnt_length, cnt_length)==0) return 0;
		//��������������
		int step_num = (int)(cnt_length/cnts_geo.step_length) + 1;
        
		//---------------------------------------------------------------------------
		//���ֲ��������ѡȡ���׹ܰ뾶
		double cnt_rad;
		if(Get_random_value(cnts_geo.rad_dist_type, cnts_geo.rad_min, cnts_geo.rad_max, seed_cnt_radius, cnt_rad)==0) return 0;
        
		//---------------------------------------------------------------------------
		//�����������ϰ����ȷֲ����ѡȡһ��������Ϊ���׹��ߵ��׸�����
		double cnt_sita, cnt_pha;
		if(Get_uniform_direction(cnts_geo, seed_cnt_sita, seed_cnt_pha, cnt_sita, cnt_pha)==0) return 0;
		MathMatrix multiplier(3,3);
		multiplier = Get_transformation_matrix(cnt_sita, cnt_pha);
        
		//---------------------------------------------------------------------------
		//���׹�ÿ����һ�������ӵ����(����û�п������������������������)
		const double step_vol_para = PI*cnt_rad*cnt_rad;
		//---------------------------------------------------------------------------
		//���׹�ÿ����һ�������ӵ�����(���������״�ϸ��ͬʱ�����ܶ�Ҳ����һ����Χ�ڿɱ�ģ�ÿ�����׹ܿ��ܲ�ͬ)
		const double step_wei_para = cnts_geo.linear_density;
        
		//---------------------------------------------------------------------------
		//Changed the oreder because I neded the radius
		vector<Point_3D> new_cnt;
		int new_cnt_size = (int)new_cnt.size();
        
		//---------------------------------------------------------------------------
		//�ֱ���RVE�ĳ�����ϰ����ȷֲ�ѡȡ���׹ܵ���ʼ��
        //Generates the seed for the CNT
		Point_3D cnt_poi;
		if(Get_seed_point(cell_geo, seed_cnt_origin, cnt_poi)==0) return 0;
        //hout << "random seed crated, ";
        
        //----------- AMC
        //Obtain the region to which the point belongs to wihtout considering overlapping
        current_region = Default_region(cnt_poi);
        //Chek if the new seed is penetrating any CNT. When 0 I need to make a new seed
        counter = 0;
        //Increase the size of the structure vector
        structure.push_back(empty);
        //COMMENT HERE TO ALLOW OVERLAPPING. JUST USE A * AT THE BEGINNIN OF THIS LINE
        /*/hout << "sp ";
        while(!Check_penetration(cnt_poi, current_region, cnt_rad, new_cnt, cnts_points, cell_geo, boundary_flag)){
            if(Get_seed_point(cell_geo, seed_cnt_origin, cnt_poi)==0) return 0;
            counter ++;
            //hout << "Seed deleted" << endl;
            if (counter == 10) {
                hout << "Too many attempts to create a new seed (" << counter << " attempts). ";
                hout << cnt_poi.x << ' ' << cnt_poi.y << ' ' << cnt_poi.z << endl;
                return 0;
            }
        }//*/
        
		new_cnt.push_back(cnt_poi);	//�洢�ڵ�
		cnt_seed_count++;					//��¼ѡ�ֳɹ�һ��
		//if(cnt_seed_count>100000) { hout << "���׹����ֳ���ʮ��Σ���û�дﵽҪ������������ ���飡" << endl; return 0; }

		//---------------------------------------------------------------------------
		//���׹�����
		int ellip_num = -1; //�����������Ŵ�ʱ�����ڱ�����׹ܴ����������ţ������׹�û�п�ʼ����ʱ���ᴩ���κ�������棻�������������Ŵ�ʱ���˱���û���õ�
		for(int i=0; i<step_num; i++)
		{
			//����Z������һ�µķ���Ϊ����, ������̬�ֲ�sita = fabs[(-omega,+omega)], ͨ��omega<=PI/2, γ����ȷֲ�pha=(0,2PI)
			//���ѡȡ����
			if(Get_normal_direction(cnts_geo.angle_max, seed_cnt_sita, seed_cnt_pha, cnt_sita, cnt_pha)==0) return 0;
			
			//�޸�����任����
			multiplier = multiplier*Get_transformation_matrix(cnt_sita, cnt_pha);
			
			//�����½ڵ�����(����任)
			cnt_poi = cnt_poi + Get_new_point(multiplier, cnts_geo.step_length);
			cnt_poi.flag = 1;							//����ʼ�˵�
            
			//---------------------------------------------------------------------------
			//������׹�������Խ�Ŵؽ��棬�������׹ܰ�һ���ĸ����ۻ�
			if(ellips.size()>0)
			{
				while(Cnt_cross_cluster_surface(clust_geo.growth_probability, ellips, ellip_num, seed_growth_probability, new_cnt.back(), cnt_poi)) //��Ҫ�ۻ�
				{
					if(ellip_num<0||ellip_num>=(int)ellips.size()) { hout << "���󣡼�¼���Ŵ������Ų�������������ŷ�Χ�ڣ����飡" << endl; return 0; }
					//�������׹��߶���������Ľ���(�Ҳ������㱨��)
					Point_3D intersect_point(0, 0 ,0);
					if(Get_intersection_cnt_ellipsoid(ellips[ellip_num], new_cnt.back(), cnt_poi, intersect_point)==0) { hout << "����û���ҵ����׹��ߺ�������Ľ��㣬���飡" << endl; return 0; }
					if(intersect_point.flag==-1)  break;  //�ж���������������߶ε�һ�������������ϣ�����ʵ�����������ⲿ
                    
					//�ƶ����׹��߶�����������˵㵽�������ڵĶԳƵ㣨������Ϊ�Գ��棩ע�⣺������һ��cnt_poi��λ���Ѿ������˱仯
					if(Move_cnt_point_symmetry_point(ellips[ellip_num], intersect_point, cnt_poi)==0) { hout << "����û���ҵ����׹��ߺ�������Ľ��㣬���飡" << endl; return 0; }
                    
					//��������н�
					if(Cos_angle_in_three_points(new_cnt.back(), intersect_point, cnt_poi)<cos(PI/6)) //�нǴ���30��
					{
						//���뽻�㵽���׹����У������flagֵ�ڼ���ʱ�Ѿ����裩
						new_cnt.push_back(intersect_point);
						
						//�޸�����任����
						multiplier = Get_vector_transformation_matrix(intersect_point, cnt_poi);
					}
					else //ȷ��û��̫��Ľ�
					{
						//�޸�����任����
						multiplier = Get_vector_transformation_matrix(new_cnt.back(), cnt_poi);
					}
				}
			}
            
			//---------------------------------------------------------------------------
			//����������ڵ������⣬ƽ�����׹�����
			if(Judge_cell_including_point(cell_geo, cnt_poi)==0)
			{
				//�����������˵��뵥����������н���(�������˵������֮�䣬����0<t<1)����������t��С����������Щ����
				vector<Point_3D> ipoi_vec;  //��������
				if(Get_intersecting_point_RVE_surface(cell_geo, new_cnt.back(), cnt_poi, ipoi_vec)==0) return 0;
                cnt_poi = ipoi_vec[0]; //The condition below now is evaluated when updating the volume and weigth fractions 
				//if(new_cnt.back()!=ipoi_vec[0]) cnt_poi = ipoi_vec[0];
                    //new_cnt.push_back(ipoi_vec[0]);
                //If it reaches the boundary and we are not considering periodicity of the boundary,
                //then the CNT is done, so we set the counter "i" to step_num
                i = step_num;
                //This for-loop handles the periodical boundary. Comment it when not needed.
				/*/for(int j=0; j<(int)ipoi_vec.size(); j++)
				{
					Point_3D inter_point[2];
                    
					inter_point[0] = ipoi_vec[j];
					if(j==(int)ipoi_vec.size()-1) inter_point[1] = cnt_poi;
					else inter_point[1] = ipoi_vec[j+1];
                    
					//����ƽ������任
					if(Periodical_coordinate_transformation(cell_geo, inter_point)==0) return 0;
                    
					if(j==(int)ipoi_vec.size()-1) cnt_poi = inter_point[1];  //����λ��
                    
					new_cnt.push_back(inter_point[0]);
					new_cnt.push_back(inter_point[1]);
				}//*/
			}
            //COMMENT HERE TO ALLOW OVERLAPPING. JUST USE A * IN THE FOLLOWING LINE
            /*/Obtain the region to which the point belongs to wihtout considering overlapping
            current_region = Default_region(cnt_poi);
            //Chek if the new point is penetrating any CNT
            //hout << "pp ";
            if(!Check_penetration(cnt_poi, current_region, cnt_rad, new_cnt, cnts_points, cell_geo, boundary_flag)){
                //When it reaches this part is because the point that was overlapping could not be accommodated
                //hout << "Penetrating point could not be accommodated" <<endl;
                //Remove the volume and weight fractions corresponding to the cnt
                for(int j=0; j<(int)new_cnt.size()-1; j++)
                {
                    if(new_cnt[j+1].flag!=0)
                    {
                        //Calculate segment length
                        double temp_length = new_cnt[j].distance_to(new_cnt[j+1]);
                        //Subtract corresponding volume
                        vol_sum -= temp_length*step_vol_para;		//�������
                        //Subtract corresponding weight
                        wt_sum -= temp_length*step_wei_para;		//��������
                    }
                }
                //Set the delete flag to 1
                delete_flag = 1;
                //Then break the for loop so a new seed is generated
                break;
            }
            //It could happen that, if the point had to be moved, the new one lies on the boundary.
            //In that case I also need to set the counter "i" to step_num as the CNT is done
            if (boundary_flag) {
                i = step_num;
            }//*/
            
			/*/---------------------------------------------------------------------------
			//�����������������
			for(int j=new_cnt_size; j<(int)new_cnt.size()-1; j++)
			{
				if(new_cnt[j+1].flag!=0)
				{
					double temp_length = new_cnt[j].distance_to(new_cnt[j+1]);
					vol_sum += temp_length*step_vol_para;		//�������
					wt_sum += temp_length*step_wei_para;		//��������
				}
			}//*/
			//---------------------------------------------------------------------------
            //Below is a version of the commented section above. It is made as a follow up to the changes to
            //the handling of the intersection with the sample's boundary. So because of the removal of the
            //periodicity of the boundary then only one point is added at a time.
            double temp_length = new_cnt.back().distance_to(cnt_poi);
            if (temp_length > 0) {
                vol_sum += temp_length*step_vol_para;		//�������
                wt_sum += temp_length*step_wei_para;		//��������
                new_cnt.push_back(cnt_poi);		//�洢�ڵ�
                new_cnt_size = (int)new_cnt.size()-1;			//�����²�����λ��
            }
            
			//---------------------------------------------------------------------------
			//�ж����������������
			if(cnts_geo.criterion == "vol"&&vol_sum >= cnts_geo.real_volume) break;			//�������, ����
			else if(cnts_geo.criterion == "wt"&&wt_sum >= cnts_geo.real_weight) break;		//��������, ����
		}
        
		//---------------------------------------------------------------------------
        if (delete_flag) {
            //If the delete flag is on, reset it and go to the next while loop to create a new seed
            delete_flag = 0;
            //Delete the last empty vector added to structure
            structure.pop_back();
            //Increase the count of rejects
            rejects++;
        } else {
            //�洢���׹ܽڵ�����
            if(new_cnt.size()<2)
            {
                //If the CNT has size less than 2, then is only one point or it's empty.
                //Delete the CNT and create a new one
                structure.pop_back();
                //Increase the count of rejects
                rejects++;
                hout << "���ɵ����׹��������� ���飡" << endl;
                //return 0;
            } else {
                //This is done only when the delete_flag is zero
                vector<Point_3D> cnt_temp;
                long int coord;
                //hout << new_cnt.size() << ": ";
                //hout << "CNTs#";
                for(int i=0; i<(int)new_cnt.size(); i++)
                {
                    //�ж��Ƿ��ǽ���
                    //Add the globlal coordinates of the points of the new CNT
                    global_point_coord.push_back(empty);
                    //The coordinate of the new point is (cnts_points.size(), new_cnt.size()-1)
                    global_point_coord.back().push_back(cnts_points.size());
                    global_point_coord.back().push_back(i);
                    //Now add to its corresponding region
                    coord = global_point_coord.size()-1;
                    structure.back().push_back(coord);
                    Add_to_regions(new_cnt[i],  coord);
                    //Change the flag to the global number
                    new_cnt[i].flag = (int)coord;
                    
                    /*/The code below is only necessary when the boundary is periodical.
                     //Otherwise I can comment it and place after this for loop the lines that are nested inside the double if
                     cnt_temp.push_back(new_cnt[i]); //���Ȳ���ڵ�
                     if(i==(int)new_cnt.size()-1||new_cnt[i+1].flag==0) //ĩβ���������ת������һ�����׹�
                     {
                     if((int)cnt_temp.size()>1)						//����1�����˵��ֻ��һ����ʼ�㣬����������������׹��߶�
                     {
                     cnts_points.push_back(cnt_temp);		//��¼���׹���Ϣ
                     cnts_radius.push_back(cnt_rad);			//��¼���׹ܰ뾶
                     }
                     cnt_temp.clear();									//�����ʱ���׹�
                     }//*/
                }
                cnts_points.push_back(new_cnt);		//��¼���׹���Ϣ
                //hout << cnts_points.size() << " P#"<<global_point_coord.size()<< endl;
                //hout << "CNT#=" << cnts_points.size()-1 << " of size " << cnts_points.back().size();
                //hout << " Struc#=" << structure.size()-1 << " of size " << structure.back().size() << endl;
                cnts_radius.push_back(cnt_rad);			//��¼���׹ܰ뾶
                if (!(cnts_points.size() % 5000))
                    hout << "CNTs generated: " << cnts_points.size() << " points generated: " << global_point_coord.size() << endl;
            }
        }
	}
    
	if(cnts_geo.criterion == "wt")
	{
		hout << "    ���ɵ����׹���ռ����ٷֱ�ԼΪ��" << vol_sum/cell_geo.volume << endl;
	}
    hout << "Total number of CNTs generated: " << cnts_points.size() << endl;
    hout << "Total number of points generated: " << global_point_coord.size() << endl;
    hout << "Total number of penetrating points: " << penetrating_points << endl;
    hout << "Total number of CNTs rejected: " << rejects << endl;
	return 1;
}

//This function adds a point to a region so penetration can be checked
long int GeoNano::Default_region(Point_3D point)
{
    //These variables are the coordinates of the lower corner of the RVE that defines its geometry
    double xmin = cell_geo.poi_min.x;
    double ymin = cell_geo.poi_min.y;
    double zmin = cell_geo.poi_min.z;
    
    //Save coordinates of the point
    double x = point.x;
    double y = point.y;
    double z = point.z;
    
    //These variables will give me the region cordinates of the region that a point belongs to
    int a, b, c;
    //Calculate the region-coordinates
    a = (int)((x-xmin)/overlap_regions.lx);
    //Limit the value of a as it has to go from 0 to secx-1
    if (a == overlap_regions.secx) a--;
    b = (int)((y-ymin)/overlap_regions.ly);
    //Limit the value of b as it has to go from 0 to secy-1
    if (b == overlap_regions.secy) b--;
    c = (int)((z-zmin)/overlap_regions.lz);
    //Limit the value of c as it has to go from 0 to secz-1
    if (c == overlap_regions.secz) c--;
    
    //I also need the "default" zone, so instead of adding one more operation to the nested for-loops above
    //I just do it again
    return (a + b*((long int)overlap_regions.secx) + c*((long int)overlap_regions.secx*overlap_regions.secy));
}

//
int GeoNano::Check_penetration(Point_3D &point, long int point_region, double cnt_rad, vector<Point_3D> cnt_new, vector<vector<Point_3D> > cnts, RVE_Geo cell_geo, int &boundary_flag)
{
    //The boundary_flag is always 0 unless in the specific case below
    boundary_flag = 0;
    //Temporary vector of int's to store the contact pair
    vector<int> empty;
    //This vector will store the points that "point" is penetrating
    vector<vector<long int> > affected_points;
    //This is the point I will compare the input "point" to.
    Point_3D region_point;
    //These ints are just to store the global point number and the CNTs they belong to.
    //They are just intermediate variables and I only use them to make the code more readable
    //int CNT1 = (int)cnts.size(), P1 = (int)cnt_new.size();
    long int coord2, P2, CNT2;
    // Calculate (maximum) cutoff distance for the region overlap.
    double cutoff_p;// = 2*cnts_geo.rad_max + d_vdw;
    vector<double> cutoffs_p;
    //This variable is to store the distance between points so it won't be calculated every time
    double distance;
    //This will be used to count the maximum number of attempts to move a point.
    int max_attempts = 5;
    int attempts;
    //I move the point up to max_attempts times. If there is still penetration then I delete it
    for (attempts = 0; attempts <= max_attempts; attempts++) {
        //Check the point agains all points in the same region
        for (long int i = 0; i < sectioned_domain[point_region].size(); i++) {
            //hout << "Check1 " << " region size = " << sectioned_domain[point_region].size() << ' ';
            coord2 = sectioned_domain[point_region][i];
            //hout << "Check2 ";
            CNT2 = global_point_coord[coord2][0];
            P2 = global_point_coord[coord2][1];
            //The points of the current CNT are in cnt_new.
            //hout << "Check3 ";
            region_point = cnts[CNT2][P2];
            //hout << "Check4 ";
            cutoff_p = cnt_rad + cnts_radius[CNT2] + d_vdw;
            //Check is the second point is in the cube of size 2cutoff_p and centered in P1
            //This is easier and faster to check than calculating the distance from poin to point every time
            if ( (region_point.x<point.x+cutoff_p)&&(region_point.x>point.x-cutoff_p)&&(region_point.y<point.y+cutoff_p)&&(region_point.y>point.y-cutoff_p)&&(region_point.z<point.z+cutoff_p)&&(region_point.z>point.z-cutoff_p) ) {
                distance = point.distance_to(region_point);
                //If it is inside the cube, then it is worth to take the time to calculate the distance from point ot point
                if (distance < cutoff_p) {
                    affected_points.push_back(global_point_coord[coord2]);
                    cutoffs_p.push_back(cutoff_p);
                    /*/hout << "CNT=" << CNT1 << " Point=" << P1 << " r1=" << cnt_rad;
                    hout << " Penetrating points="<< affected_points.size();
                    hout << " CNT2=" << CNT2 << " P2=" << P2 << " r2=" << cnts_radius[CNT2] << " (" << region_point.x << ", " << region_point.y << ", " << region_point.z << ") ";
                     hout << endl;//*/
                }
            }
            //hout << "Check5 " << endl;
        }
        //-----------------------
        //Find new point according to the number of penetrations.
        if (affected_points.size()) {
            //If this is the last iteration and there are still affected points then the point could not be accommodated
            if (attempts == max_attempts) {
                /*hout << "Delete CNT number " << CNT1 << " of size " << P1 ;
                 hout << " (reached maximum number of attempts for relocation)" << endl;//*/
                return 0;
            }
            
            //Find the new point
            //hout << "Moved a point from initial position (" << point.x << ", " << point.y << ", " << point.z << ")." << endl;
            Move_point(point, cutoffs_p, affected_points, cnts);
            //hout << "Moved a point to final position (" << point.x << ", " << point.y << ", " << point.z << ")." << endl;
            //hout << "Check6 ";
            
            //Check that the new point is within the permitted orientation repect to the previous segment
            if (!Check_segment_orientation(point, cnt_new)) {
                /*hout << "Delete CNT number " << CNT1 << " of size " << P1 ;
                 hout << " (the point is not in a valid orientation)" << endl;//*/
                //When not in a valid position it cannot be moved again so a new CNT is needed
                return 0;
            }
            //hout << "Check7 ";
            
            //If the new point is outside, find the intersection wih the boundary and use that point
            if (!Judge_cell_including_point(cell_geo, point)) {
                vector<Point_3D> ipoi_vec;
                //hout << "Check7.1 " << cnt_new.size();
                //When the size of cnt_new is zero then we have a seed. If the seed is outside the boundary
                //it seems easier to just create a new one
                if (!cnt_new.size()) {
                    return 0;
                }
				if(Get_intersecting_point_RVE_surface(cell_geo, cnt_new.back(), point, ipoi_vec)==0) {
                    //If the point could not be projected into the boundary, then delete it
                    /*hout << "Delete CNT number " << CNT1 << " of size " << P1 ;
                     hout << " (point coud not be projected into the boundary)" << endl;//*/
                    return 0;
                }
                //hout << "Check7.2 ";
                point = ipoi_vec[0];
                boundary_flag = 1;
                //hout << "Moved to intersection with boundary (" << point.x << ", " << point.y << ", " << point.z << ')' << endl;
            }
            //hout << "Check8 ";
            
            //Need to update point region just in case
            point_region = Default_region(point);
            
            //Need to clear the vectors affected_points, contact_coordinates and temporal_contacts so they are used again with the new point
            affected_points.clear();
            cutoffs_p.clear();
            //hout << "Check9 ";
            
        } else {
            break; //if the size of affected_points is zero, then break the loop
        }
    }
    
    //This is just to know in how many attempts a point was sucessfully moved
     if (attempts) {
         penetrating_points++;
     //hout << "Point relocated in " << attempts << " attempts." << endl;
     }//*/
    //hout << "Check10 " << endl;
    return 1;
}

//This function moves a point according to the number of points it is penetrating
void GeoNano::Move_point(Point_3D &point, vector<double> cutoff, vector<vector<long int> > affected_points, vector<vector<Point_3D> > cnts)
{
    Point_3D direction;
    //Point_3D original = point;
    Point_3D P, Q, R, P1, P2, P3;
    double length, a, b, c, d;
    //When more than four, in this variable the closset three are stored
    vector<vector<long int> > three;
    int indices[] = {0, 1, 2};
    //This will determine how the new point is moved
    int intersections = (int) affected_points.size();
    //It is rare but it might happen that a penetrating point is too close to "point" that their distance is too close to Zero or it is
    //actually 0. In that case we make an adjustment to the intersection number
    for (int i = 0; i < (int) affected_points.size(); i++) {
        P = cnts[affected_points[i][0]][affected_points[i][1]];
        if (point.distance_to(P) <= Zero) {
            intersections = 0;
            break;
        }
    }
    if (!intersections) {
        //In the case that intersections is zero then artificially move each coordinate of point 0.1steplength in each direction
        //Hope his might not affect the result
        length = 0.01*cnts_geo.step_length;
        point.x += length;
        point.y += length;
        //point.z += length;
        hout << "Had to move the point 0.01steplength in each direction to avoid errors. Hopefully this will not affect. Reduction in case from " << affected_points.size();
        hout << " to 0." << endl;
    } else {
        switch (intersections) {
            case 1:
                //When there is penetration with one point only, then this is the simplest and easiest case
                //Just translate in the direction form P_old to P_new a distance cutoff from P_old
                
                //Get the penetrating point. Its coordinates are in the first (and only) element of affected_points
                P = cnts[affected_points[0][0]][affected_points[0][1]];
                //Calculate direction unit vector
                direction = (point - P)/(point.distance_to(P));
                //So the new point is P_old + d*n. d is the cutoff and n the unit vector
                point = P + direction*(cutoff[0]+Zero); //The Zero is to avoid machine precision errors. Without it, when comparing
                //the new point with the other points in the same region, the program was judging them to be below the cutoff
                //for the van der Waals distance. Even though they were in the limit. After adding this Zero that issue
                //was eliminated
                break;
            case 2:
                //
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
                a = cutoff[0];
                b = cutoff[1];
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
                break;
            default:
                //When there are 3 or more penetrating points, then we look for the three smallest ones and reduce
                //the case to three by recursively calling this function
                
                //Get the first three penetrationg points.
                P1 = cnts[affected_points[0][0]][affected_points[0][1]];
                P2 = cnts[affected_points[1][0]][affected_points[1][1]];
                
                //Find the distances to each point
                a = point.distance_to(P1);
                b = point.distance_to(P2);
                
                if (a < b) {
                    //Change the order of the points
                    indices[0] = 1;
                    indices[1] = 0;
                    length = a;
                    a = b;
                    b = length;
                }
                
                //Check the rest of points and if a closer point is found then add it
                for (int i = 2; i < (int)affected_points.size(); i++) {
                    length = point.distance_to(cnts[affected_points[i][0]][affected_points[i][1]]);
                    if (length < b) {
                        if (length < a) {
                            //The new point is the smallest
                            b = a;
                            a = length;
                            indices[1] = indices[0];
                            indices[0] = i;
                        } else {
                            //The new point is the second smallest
                            b = length;
                            indices[1] = i;
                        }
                    }
                }
                
                //Use the indices in indices[] to reduce the case
                three.push_back(affected_points[indices[0]]);
                three.push_back(affected_points[indices[1]]);
                a = cutoff[indices[0]];
                b = cutoff[indices[1]];
                cutoff.clear();
                cutoff.push_back(a);
                cutoff.push_back(b);
                //Now that I have the three closest points I call this function recursively with case 3
                //hout << "Reducing case from " << intersections << " to 3. ";
                Move_point(point, cutoff, three, cnts);
                break;
        }
        //hout << "Moved a point from (" << original.x << ", " << original.y << ", " << original.z;
        //hout << ") to (" << point.x << ", " << point.y << ", " << point.z << ')' << endl;
    }
}

//This function checks that "point" is within the bounds of the segment orientation.
//The criterion is just checking the point is not more than pi/2 respect with the previous
//segment. In the limiting case we have a straight triangle. So I calculate the hypotenuse.
//I also measure the distance between "point" and the second before that.
//If the distance  between points is less than the hypotenuse, then it has an
//incorrect orientation
int GeoNano::Check_segment_orientation(Point_3D point, vector<Point_3D> cnt)
{
    //Check the size of the cnt. If has to be at least 2, otherwise it does not matter
    //where "point" is. It will always be in a valid position
    if (cnt.size()>=2) {
        int last = (int)cnt.size()-1;
        double hypotenuse = point.distance_to(cnt[last])*point.distance_to(cnt[last]);
        hypotenuse  = hypotenuse + cnt[last].distance_to(cnt[last-1])*cnt[last].distance_to(cnt[last-1]);
        hypotenuse = sqrt(hypotenuse);
        //When hypotenuse is equal to segment, there could be a machine precision error.
        //In this case, they are equal but judged to be different. Then I use the Zero to reduce the
        //incidence of that issue. (Just as I did in the function Move_point)
        //I subtract the zero so the limiting case is included when segment => hypothenuse
        double segment = cnt[last-1].distance_to(point) - Zero;
        if (segment<hypotenuse) {
            //The point is not in a valid position
            return 0;
        } else {
            //The point is in a valid position
            return 1;
        }
    } else
        return 1;
}

//This function adds a point to a region so penetration can be checked
void GeoNano::Add_to_regions(Point_3D point, long int point_coordinate)
{
    //Save coordinates of the point
    double x = point.x;
    double y = point.y;
    double z = point.z;

    //These variables are the coordinates of the lower corner of the RVE that defines its geometry
    double xmin = cell_geo.poi_min.x;
    double ymin = cell_geo.poi_min.y;
    double zmin = cell_geo.poi_min.z;
    
    //These variables will give me the region cordinates of the region that a point belongs to
    int a, b, c;
    long unsigned t;

    //COMMENT HERE STARTING ON LINE BELOW TO ALLOW OVERLAPPING. JUST ADD '*'
    /*/These variables are to store the size of the RVE and reduces operations when accessing them
    double lx = cell_geo.len_x;
    double ly = cell_geo.wid_y;
    double lz = cell_geo.hei_z;
    
    //Sizes of each region
    double dx = overlap_regions.lx;
    double dy = overlap_regions.ly;
    double dz = overlap_regions.lz;
    
    
    //Calculate the region-coordinates
    a = (int)((x-xmin)/overlap_regions.lx);
    //Limit the value of a as it has to go from 0 to secx-1
    if (a == overlap_regions.secx) a--;
    b = (int)((y-ymin)/overlap_regions.ly);
    //Limit the value of b as it has to go from 0 to secy-1
    if (b == overlap_regions.secy) b--;
    c = (int)((z-zmin)/overlap_regions.lz);
    //Limit the value of c as it has to go from 0 to secz-1
    if (c == overlap_regions.secz) c--;
    
    //-------------------------------
    //Manage the overlaping of the regions
    
    //Coordinates of non-overlaping region the point belongs to
    double x1 = a*dx +  xmin;
    double y1 = b*dy +  ymin;
    double z1 = c*dz +  zmin;
    double x2, y2, z2;
    if (a == overlap_regions.secx-1) x2 = lx +  xmin;
    else x2 = (a+1)*dx +  xmin;
    if (b == overlap_regions.secy-1) y2 = ly +  ymin;
    else y2 = (b+1)*dy +  ymin;
    if (c == overlap_regions.secz-1) z2 = lz +  zmin;
    else z2 = (c+1)*dz +  zmin;
    
    //Initialize flags for overlaping regions
    int fx = 0;
    int fy = 0;
    int fz = 0;
    
    //Assign value of flag according to position of point
    //The first operand eliminates the periodicity on the boundary and as a consequence invalid coordinates a,b,c
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
    long int region_coord[2][3] = { {(long int)a, (long int)b, (long int)c}, {(long int)a+fx, (long int)b+fy, (long int)c+fz}};
    
    //In this loop I check all 8 regions a point can belong to when that point is in an overlaping zone
    for (int ii = 0; ii < 2; ii++) {
        if (!fx) ii++; //if flag is zero, do this loop only once
        for (int jj = 0; jj < 2; jj++) {
            if (!fy) jj++; //if flag is zero, do this loop only once
            for (int kk = 0; kk < 2; kk++) {
                if (!fz) kk++; //if flag is zero, do this loop only once
                //t = a + b*sx + c*sx*sy;
                //a = temp[ii][0], b = temp[jj][1], c = temp[kk][2]
                t = region_coord[ii][0] + region_coord[jj][1]*((long int)overlap_regions.secx) + region_coord[kk][2]*((long int)overlap_regions.secx*overlap_regions.secy);
                //Store the global point number on its corresponding region or regions
                sectioned_domain[t].push_back(point_coordinate);
                //Store the CNT number on its corresponding region or regions
            }
        }
    }//*/
    
    //---------------------------------
    //Now assign the region for the CNT
    
    //Calculate the region-coordinates
    a = (int)((x-xmin)/cnt_regions.lx);
    //Limit the value of a as it has to go from 0 to secx-1
    if (a == cnt_regions.secx) a--;
    b = (int)((y-ymin)/cnt_regions.ly);
    //Limit the value of b as it has to go from 0 to secy-1
    if (b == cnt_regions.secy) b--;
    c = (int)((z-zmin)/cnt_regions.lz);
    //Limit the value of c as it has to go from 0 to secz-1
    if (c == cnt_regions.secz) c--;
    
    //Calculate the region the point belongs to
    t = a + b*cnt_regions.secx + c*cnt_regions.secx*cnt_regions.secy;
    
    //Add the CNT if the corresponding region is empty or the last element is not the current cnt
    int CNT = (int)global_point_coord[point_coordinate][0];
    //hout << " " << global_point_coord.size() << " ";
    if ((!sectioned_domain_cnt[t].size())||sectioned_domain_cnt[t].back() != CNT) {
        sectioned_domain_cnt[t].push_back(CNT);
    }
    //hout << "sectioned_domain_cnt.size() " << sectioned_domain_cnt.size() << endl;
}

//---------------------------------------------------------------------------
//����ռ�һ���������λ�ù�ϵ��������ֵ
double GeoNano::Estimate_position_point_ellipsoid(const struct elliparam &temp_ellips, const Point_3D &point)const
{
	double x, y, z;
	x=point.x-temp_ellips.x;
	y=point.y-temp_ellips.y;
	z=point.z-temp_ellips.z;
	
	double x1, y1, z1;
	x1=x*temp_ellips.alpha1+y*temp_ellips.beta1+z*temp_ellips.gamma1;
	y1=x*temp_ellips.alpha2+y*temp_ellips.beta2+z*temp_ellips.gamma2;
	z1=x*temp_ellips.alpha3+y*temp_ellips.beta3+z*temp_ellips.gamma3;
    
	return pow(x1,2)/pow(temp_ellips.a, 2)+pow(y1,2)/pow(temp_ellips.b, 2)+pow(z1,2)/pow(temp_ellips.c, 2)-1.0;
}
//---------------------------------------------------------------------------
//�������׹��Ŵ���������
int GeoNano::Get_ellip_clusters(const struct RVE_Geo &cell, struct Clust_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle, vector<struct elliparam> &ellips, const int &export_mod)const
{
	double epsilon = 0.01;
	int times = 0;
	double ellip_volume = 0.0;
	//������������
	do
	{
		//-------------------------------------------------------------
		//��������
		struct elliparam ell_temp;
		//�������ĵ�λ��
		seed_poi = (2053*seed_poi + 13849)%MAX_INT;
		ell_temp.x=cell.poi_min.x + seed_poi*cell.len_x/MAX_INT;
        
		seed_poi = (2053*seed_poi + 13849)%MAX_INT;
		ell_temp.y=cell.poi_min.y + seed_poi*cell.wid_y/MAX_INT;
        
		seed_poi = (2053*seed_poi + 13849)%MAX_INT;
		ell_temp.z=cell.poi_min.z + seed_poi*cell.hei_z/MAX_INT;
        
		//�����ж���
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
        
		//����Ÿ������[(alpha1,beta1,gamma1),(alpha2,beta2,gamma2),(alpha3,beta3,gamma3)]�ֱ���������(a,b,c)��(x,y,z)����ļн�
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
        
		ell_temp.alpha1	=	cos(alpha1);																//alpha1��(0, PI)��ѡ
		ell_temp.beta1	=	cos(beta1);																	//beta1��(pi/2-r1)��(pi/2+r1)ѡȡ
		seed_angle	= (2053*seed_angle + 13849)%MAX_INT;
		ell_temp.gamma1 = pow(-1.0, fmod(seed_angle, 2.0)+1.0)*sqrt(1.0-pow(ell_temp.alpha1,2)-pow(ell_temp.beta1,2));	 //����gamma��ֵ�����ȡ��������
		double alpha2 = 0;																						//alpha2��(pi/2-r1)��(pi/2+r1)ѡȡ
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
        
		ell_temp.a = (1+epsilon)*ell_temp.a;
		ell_temp.b = (1+epsilon)*ell_temp.b;
		ell_temp.c = (1+epsilon)*ell_temp.c;
        
		//-------------------------------------------------------------
		//��������Ƿ��뵥���������������ཻ
		double delt_h = ell_temp.c/50;						//�ָ�Ĵ���̫��Ӱ�����Ч?
		int k1 = (int)(sqrt(pow(ell_temp.a,2)+pow(ell_temp.b,2))/delt_h);
		int K = 4*(k1+1);
		double sita = 2*PI/K;
        
		//����Ƿ��뵥���߽��ཻ
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
                
				if(x1-cell.poi_min.x<Zero||x1-cell.poi_min.x>cell.len_x-Zero||
                   y1-cell.poi_min.y<Zero||y1-cell.poi_min.y>=cell.wid_y-Zero||
                   z1-cell.poi_min.z<Zero||z1-cell.poi_min.z>=cell.hei_z-Zero)
				{
					times=times+1;
					goto gen_again;
				}
			}
		}
		//����Ƿ������������ཻ
		for(int i=0; i<(int)ellips.size(); i++)
		{
			//���Թ���
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
				//��ȷ����
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
		//�������㲢������������
		times=0;
		ellips.push_back(ell_temp);
		//---------------------------------------------------------------------
		//���������
		ellip_volume = ellip_volume + 4*PI*ell_temp.a*ell_temp.b*ell_temp.c/(3*pow(1+epsilon, 3.0));
		clust_geo.real_volume_fraction = ellip_volume/cell.volume;
    gen_again:	;
	}while(times<=N_times&&clust_geo.real_volume_fraction<clust_geo.vol_fra_criterion);
	
	//---------------------------------------------------------------------
	//��������
	for(int i=0; i<(int)ellips.size(); i++)
	{
		ellips[i].a=ellips[i].a/(1+epsilon);
		ellips[i].b=ellips[i].b/(1+epsilon);
		ellips[i].c=ellips[i].c/(1+epsilon);
	}
	//---------------------------------------------------------------------
	//������׹��Ŵ���������
	if(export_mod==2)	Export_cluster_ellipsoids_mesh(cell, ellips);
    
	//---------------------------------------------------------------------
	//������׹��Ŵ���������
	if(export_mod==1||export_mod==2)	Export_cluster_ellipsoids_data(ellips, clust_geo.real_volume_fraction);
    
	//������ģʽ
	if(export_mod>2||export_mod<0) { hout <<  "�������׹��Ŵ�������������ģʽ���� ���飡" << endl; return 0; }
	else {	hout << "    ���׹��Ŵ�����ĸ�����ʵ�������" << (int)ellips.size() << "  " << clust_geo.real_volume_fraction << endl; }
    
	return 1;
}
//---------------------------------------------------------------------------
//�����ض�λ�õ�Բ���Ŵ�����(����������ʱ�����������޷�����)
int GeoNano::Get_specific_sphere_clusters(const struct RVE_Geo &cell, struct Clust_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle, vector<struct elliparam> &ellips, const int &export_mod)const
{
	int snum = 2;	//����ÿһ��������Բ�������
	double sd_x = 0.5*cell.len_x/snum;
	double sd_y = 0.5*cell.wid_y/snum;
	double sd_z = 0.5*cell.hei_z/snum;
	
	if(clust_geo.amin!=clust_geo.amax||clust_geo.bmin!=0||clust_geo.cmin!=0) { hout << "ע�⣡�����������ض�λ��Բ���Ŵ����е����������飡" << endl; return 0; }
    
	if(sd_x<=clust_geo.amin||sd_y<=clust_geo.amin||sd_z<=clust_geo.amin) { hout << "ע�⣡������һ��������Բ����������࣬����" << endl; return 0; }
    
	double ellip_volume = 0.0;
	for(int i=0; i<snum; i++)
		for(int j=0; j<snum; j++)
			for(int k=0; k<snum; k++)
			{
				//-------------------------------------------------------------
				//��������
				struct elliparam ell_temp;
				ell_temp.x	=	(2*k+1)*sd_x;
				ell_temp.y	=	(2*j+1)*sd_y;
				ell_temp.z	=	(2*i+1)*sd_z;
                
				ell_temp.a=clust_geo.amin;
				ell_temp.b = ell_temp.a;
				ell_temp.c = ell_temp.a;
                
				//����Ÿ������[(alpha1,beta1,gamma1),(alpha2,beta2,gamma2),(alpha3,beta3,gamma3)]�ֱ���������(a,b,c)��(x,y,z)����ļн�
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
                
				ell_temp.alpha1	=	cos(alpha1);																//alpha1��(0, PI)��ѡ
				ell_temp.beta1	=	cos(beta1);																	//beta1��(pi/2-r1)��(pi/2+r1)ѡȡ
				seed_angle	= (2053*seed_angle + 13849)%MAX_INT;
				ell_temp.gamma1 = pow(-1.0, fmod(seed_angle, 2.0)+1.0)*sqrt(1.0-pow(ell_temp.alpha1,2)-pow(ell_temp.beta1,2));	 //����gamma��ֵ�����ȡ��������
				double alpha2 = 0;																						//alpha2��(pi/2-r1)��(pi/2+r1)ѡȡ
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
				//������������
				ellips.push_back(ell_temp);
				//---------------------------------------------------------------------
				//���������
				ellip_volume = ellip_volume + 4*PI*ell_temp.a*ell_temp.b*ell_temp.c/3;
				clust_geo.real_volume_fraction = ellip_volume/cell.volume;
			}
    
	//������������Ƿ���
	if(clust_geo.real_volume_fraction<clust_geo.vol_fra_criterion)
	{
		hout << "ע�⣡����Բ���Ŵ����е��������" << clust_geo.real_volume_fraction;
		hout << "С�����������׼" << clust_geo.vol_fra_criterion << "������" << endl;
		return 0;
	}
    
	//---------------------------------------------------------------------
	//������׹��Ŵ���������
	if(export_mod==2)	Export_cluster_ellipsoids_mesh(cell, ellips);
    
	//---------------------------------------------------------------------
	//������׹��Ŵ���������
	if(export_mod==1||export_mod==2)	Export_cluster_ellipsoids_data(ellips, clust_geo.real_volume_fraction);
    
	//������ģʽ
	if(export_mod>2||export_mod<0) { hout <<  "�������׹��Ŵ�������������ģʽ���� ���飡" << endl; return 0; }
	else {	hout << "    ���׹��Ŵ�����ĸ�����ʵ�������" << (int)ellips.size() << "  " << clust_geo.real_volume_fraction << endl; }
    
	return 1;
}
//---------------------------------------------------------------------------
//������׹��Ŵ���������
void GeoNano::Export_cluster_ellipsoids_mesh(const struct RVE_Geo &cell, const vector<struct elliparam> &ellips)const
{
	ofstream otec("Cluster_Ellipsoid_Mesh.dat");
	otec << "TITLE = Cluster_Ellipsoid_Mesh" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;
    
	//---------------------------------------------------------------------------
	//����������
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
	otec << endl;
    
	for(int i=0; i<(int)ellips.size(); i++)
	{
		const int num_sita = 20;
		const int num_phi = int(2*num_sita*ellips[i].a/ellips[i].c+0.5);  //��������ȡ��
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
//������׹��Ŵ���������
void GeoNano::Export_cluster_ellipsoids_data(const vector<struct elliparam> &ellips, const double &ellip_ratio)const
{
	ofstream out("Cluster_Ellipsoids_Data.dat");
	out <<"%�������Ŵ�����ĸ�����ʵ����ռ�������" << endl;
	out << (int)ellips.size() << " " << ellip_ratio <<endl;
	out <<"%�������15������:������������(x,y,z)������ߴ�(a,b,c)����λ��(alpha1,beta1,gamma1; alpha2,beta2,gamma2; alpha3,beta3,gamma3)" << endl;
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
//�ж����׹��Ƿ�Խ�Ŵ��������
int GeoNano::Cnt_cross_cluster_surface(const double &growth_probability, const vector<struct elliparam> &ellips, int &ellip_num, int &seed, const Point_3D &point0, const Point_3D &point1)const
{
	double x, y, z;
	double x1, y1, z1;
	double f0, f1;
	if(ellip_num<-1||ellip_num>=(int)ellips.size()) { hout << "ע�⣡��¼�����Ŵ��󣬲������ڷ�Χ֮�ڣ����飡" << endl; return 0; }
	if(ellip_num==-1)  //ѭ������
	{
		for(int i=0; i<(int)ellips.size(); i++)
		{
			x=point0.x-ellips[i].x;
			y=point0.y-ellips[i].y;
			z=point0.z-ellips[i].z;
			
			x1=x*ellips[i].alpha1+y*ellips[i].beta1+z*ellips[i].gamma1;
			y1=x*ellips[i].alpha2+y*ellips[i].beta2+z*ellips[i].gamma2;
			z1=x*ellips[i].alpha3+y*ellips[i].beta3+z*ellips[i].gamma3;
			
			f0=pow(x1,2)/pow(ellips[i].a, 2)+pow(y1,2)/pow(ellips[i].b, 2)+pow(z1,2)/pow(ellips[i].c, 2)-1.0;
            
			if(f0<=Zero) //�ڴ������ڲ�
			{
				x=point1.x-ellips[i].x;
				y=point1.y-ellips[i].y;
				z=point1.z-ellips[i].z;
                
				x1=x*ellips[i].alpha1+y*ellips[i].beta1+z*ellips[i].gamma1;
				y1=x*ellips[i].alpha2+y*ellips[i].beta2+z*ellips[i].gamma2;
				z1=x*ellips[i].alpha3+y*ellips[i].beta3+z*ellips[i].gamma3;
                
				f1=pow(x1,2)/pow(ellips[i].a, 2)+pow(y1,2)/pow(ellips[i].b, 2)+pow(z1,2)/pow(ellips[i].c, 2)-1.0;
                
				if(f1>Zero)
				{
					seed = (2053*seed + 13849)%MAX_INT;  //�����Ŵ��ڲ�������?
					if((double)seed/MAX_INT<=growth_probability)
					{
						ellip_num = i;
						return 1;	//���׹�Ҫ���Ŵ����ۻ�����
					}
					else return 0; //��������
				}
				else return 0;	//��������
			}
		}
	}
	else
	{
		//��γ���Ϊ�˼�����������ཻ��ļ��������������ʱ���˶γ������ע�͵�
		x=point0.x-ellips[ellip_num].x;
		y=point0.y-ellips[ellip_num].y;
		z=point0.z-ellips[ellip_num].z;
        
		x1=x*ellips[ellip_num].alpha1+y*ellips[ellip_num].beta1+z*ellips[ellip_num].gamma1;
		y1=x*ellips[ellip_num].alpha2+y*ellips[ellip_num].beta2+z*ellips[ellip_num].gamma2;
		z1=x*ellips[ellip_num].alpha3+y*ellips[ellip_num].beta3+z*ellips[ellip_num].gamma3;
        
		f0=pow(x1,2)/pow(ellips[ellip_num].a, 2)+pow(y1,2)/pow(ellips[ellip_num].b, 2)+pow(z1,2)/pow(ellips[ellip_num].c, 2)-1.0;
        
		if(f0>Zero) { hout << "����! �˵㱾��Ӧ��������" << ellip_num << "�ڲ�, �����ڼ���ȴ�������ⲿ(f0>Zero), ���飡 " << f0 << endl; return 0; } //�ж����
        
		x=point1.x-ellips[ellip_num].x;
		y=point1.y-ellips[ellip_num].y;
		z=point1.z-ellips[ellip_num].z;
        
		x1=x*ellips[ellip_num].alpha1+y*ellips[ellip_num].beta1+z*ellips[ellip_num].gamma1;
		y1=x*ellips[ellip_num].alpha2+y*ellips[ellip_num].beta2+z*ellips[ellip_num].gamma2;
		z1=x*ellips[ellip_num].alpha3+y*ellips[ellip_num].beta3+z*ellips[ellip_num].gamma3;
        
		f1=pow(x1,2)/pow(ellips[ellip_num].a, 2)+pow(y1,2)/pow(ellips[ellip_num].b, 2)+pow(z1,2)/pow(ellips[ellip_num].c, 2)-1.0;
        
		if(f1>Zero)
		{
			seed = (2053*seed + 13849)%MAX_INT;  //�����Ŵ��ڲ�������?
			if((double)seed/MAX_INT<=growth_probability)	return 1;	//���׹�Ҫ���Ŵ����ۻ�����
			else { ellip_num = -1; return 0; } //��������
		}
		else return 0;	//��������
	}
    
	return 0;  //û�д�Խ�Ŵر��棬��������
}
//---------------------------------------------------------------------------
//�ж����׹��Ƿ�Խһ���Ŵ��������
int GeoNano::Cnt_go_through_a_cluster(const double &growth_probability, const struct elliparam &ellip, int &seed, const Point_3D &point0, const Point_3D &point1)const
{
	//--------------------------------------------------------------------------------------------------------
	//��γ���Ϊ�˼�����������ཻ��ļ�����
	double x=point0.x-ellip.x;
	double y=point0.y-ellip.y;
	double z=point0.z-ellip.z;
    
	double x1=x*ellip.alpha1+y*ellip.beta1+z*ellip.gamma1;
	double y1=x*ellip.alpha2+y*ellip.beta2+z*ellip.gamma2;
	double z1=x*ellip.alpha3+y*ellip.beta3+z*ellip.gamma3;
    
	double ef=pow(x1,2)/pow(ellip.a, 2)+pow(y1,2)/pow(ellip.b, 2)+pow(z1,2)/pow(ellip.c, 2)-1.0;
    
	if(ef>Zero) { hout << "����! �˵㱾��Ӧ���������ڲ�, �����ڼ���ȴ�������ⲿ(ef>Zero), ���飡 " << ef << endl; return 0; } //�ж����
    
	x=point1.x-ellip.x;
	y=point1.y-ellip.y;
	z=point1.z-ellip.z;
    
	x1=x*ellip.alpha1+y*ellip.beta1+z*ellip.gamma1;
	y1=x*ellip.alpha2+y*ellip.beta2+z*ellip.gamma2;
	z1=x*ellip.alpha3+y*ellip.beta3+z*ellip.gamma3;
    
	ef=pow(x1,2)/pow(ellip.a, 2)+pow(y1,2)/pow(ellip.b, 2)+pow(z1,2)/pow(ellip.c, 2)-1.0;
    
	if(ef>Zero)
	{
		seed = (2053*seed + 13849)%MAX_INT;  //�����Ŵ��ڲ�������?
		if((double)seed/MAX_INT<=growth_probability)	return 1;	//���׹�Ҫ���Ŵ����ۻ�����
		else  return 0;
	}
	
	return 0;  //û�д�Խ�Ŵر��棬��������
}

//---------------------------------------------------------------------------
//�����׹�����������Ľ���
int GeoNano::Get_intersection_cnt_ellipsoid(struct elliparam &temp_ellips, const Point_3D &point0, const Point_3D &point1, Point_3D &intersect_point)const
{
	double x, y, z;
	Point_3D poi[2];
	//���׹������˵�����任
	x=point0.x-temp_ellips.x;
	y=point0.y-temp_ellips.y;
	z=point0.z-temp_ellips.z;
    
	poi[0].x=x*temp_ellips.alpha1+y*temp_ellips.beta1+z*temp_ellips.gamma1;
	poi[0].y=x*temp_ellips.alpha2+y*temp_ellips.beta2+z*temp_ellips.gamma2;
	poi[0].z=x*temp_ellips.alpha3+y*temp_ellips.beta3+z*temp_ellips.gamma3;
    
	x=point1.x-temp_ellips.x;
	y=point1.y-temp_ellips.y;
	z=point1.z-temp_ellips.z;
    
	poi[1].x=x*temp_ellips.alpha1+y*temp_ellips.beta1+z*temp_ellips.gamma1;
	poi[1].y=x*temp_ellips.alpha2+y*temp_ellips.beta2+z*temp_ellips.gamma2;
	poi[1].z=x*temp_ellips.alpha3+y*temp_ellips.beta3+z*temp_ellips.gamma3;
    
	//���㽻��Ĳ���ֵ
	double A=pow(poi[1].x-poi[0].x,2)/pow(temp_ellips.a,2)+pow(poi[1].y-poi[0].y,2)/pow(temp_ellips.b,2)+pow(poi[1].z-poi[0].z,2)/pow(temp_ellips.c,2);
	double B=2.0*((poi[1].x-poi[0].x)*poi[0].x/pow(temp_ellips.a,2)+(poi[1].y-poi[0].y)*poi[0].y/pow(temp_ellips.b,2)+(poi[1].z-poi[0].z)*poi[0].z/pow(temp_ellips.c,2));
	double C=pow(poi[0].x,2)/pow(temp_ellips.a,2)+pow(poi[0].y,2)/pow(temp_ellips.b,2)+pow(poi[0].z,2)/pow(temp_ellips.c,2)-1;
	double t[2]={0, 0};
	t[0]=(-B+sqrt(B*B-4*A*C))/(2*A);
	t[1]=(-B-sqrt(B*B-4*A*C))/(2*A);
    
	int count = 0;
	double tt = 0;
	for(int i=0; i<2; i++)
	{
		if(t[i]>Zero&&t[i]<1-Zero)
		{
			tt=t[i];
			count++;
		}
	}
	if(count!=1)
	{
		if(count==0&&t[0]>-1.0E-6)	{ intersect_point.flag = -1; return 1; }   //��������ʵ�������߶����⣬ԭ������Ϊ0�ŵ�������������򣨽���-1.0E-6����һ�γ������Ϣ�����ֵ��
		else { hout << "�����߶���������û�н�������߶���������ͬʱ����������, ���飡" << " t[0]=" << t[0] << " t[1]=" << t[1] << endl; return 0; }
	}
    
	//����任��ԭ����λ��
	x=poi[0].x+tt*(poi[1].x-poi[0].x);
	y=poi[0].y+tt*(poi[1].y-poi[0].y);
	z=poi[0].z+tt*(poi[1].z-poi[0].z);
    
	intersect_point.x=temp_ellips.x+x*temp_ellips.alpha1+y*temp_ellips.alpha2+z*temp_ellips.alpha3;
	intersect_point.y=temp_ellips.y+x*temp_ellips.beta1+y*temp_ellips.beta2+z*temp_ellips.beta3;
	intersect_point.z=temp_ellips.z+x*temp_ellips.gamma1+y*temp_ellips.gamma2+z*temp_ellips.gamma3;
	intersect_point.flag = 1; //�������������ʵ�ҵ��ĵ�
	
	return 1;
}
//---------------------------------------------------------------------------
//�ƶ����׹��߶�����������˵㵽�������ڵĶԳƵ㣨������Ϊ�Գ��棩
int GeoNano::Move_cnt_point_symmetry_point(struct elliparam &temp_ellips, const Point_3D &intersect_point, Point_3D &cnt_poi)const
{
	//����任
	double x, y, z;
	Point_3D poi[2];
	x=intersect_point.x-temp_ellips.x;
	y=intersect_point.y-temp_ellips.y;
	z=intersect_point.z-temp_ellips.z;
    
	poi[0].x=x*temp_ellips.alpha1+y*temp_ellips.beta1+z*temp_ellips.gamma1;
	poi[0].y=x*temp_ellips.alpha2+y*temp_ellips.beta2+z*temp_ellips.gamma2;
	poi[0].z=x*temp_ellips.alpha3+y*temp_ellips.beta3+z*temp_ellips.gamma3;
    
	x=cnt_poi.x-temp_ellips.x;
	y=cnt_poi.y-temp_ellips.y;
	z=cnt_poi.z-temp_ellips.z;
    
	poi[1].x=x*temp_ellips.alpha1+y*temp_ellips.beta1+z*temp_ellips.gamma1;
	poi[1].y=x*temp_ellips.alpha2+y*temp_ellips.beta2+z*temp_ellips.gamma2;
	poi[1].z=x*temp_ellips.alpha3+y*temp_ellips.beta3+z*temp_ellips.gamma3;
    
	double A = pow(poi[0].x,2)/pow(temp_ellips.a,4)+pow(poi[0].y,2)/pow(temp_ellips.b,4)+pow(poi[0].z,2)/pow(temp_ellips.c,4);
	double B = 1-poi[0].x*poi[1].x/pow(temp_ellips.a,2)-poi[0].y*poi[1].y/pow(temp_ellips.b,2)-poi[0].z*poi[1].z/pow(temp_ellips.c,2);
    
	x=poi[1].x+2.0*poi[0].x*B/(pow(temp_ellips.a,2)*A);
	y=poi[1].y+2.0*poi[0].y*B/(pow(temp_ellips.b,2)*A);
	z=poi[1].z+2.0*poi[0].z*B/(pow(temp_ellips.c,2)*A);
    
	//����任
	cnt_poi.x=temp_ellips.x+x*temp_ellips.alpha1+y*temp_ellips.alpha2+z*temp_ellips.alpha3;
	cnt_poi.y=temp_ellips.y+x*temp_ellips.beta1+y*temp_ellips.beta2+z*temp_ellips.beta3;
	cnt_poi.z=temp_ellips.z+x*temp_ellips.gamma1+y*temp_ellips.gamma2+z*temp_ellips.gamma3;
    
	return 1;
}
//---------------------------------------------------------------------------
//�������������м�һ�㴦���еļнǵ�Cosֵ
double GeoNano::Cos_angle_in_three_points(const Point_3D &point0, const Point_3D &point1, const Point_3D &point2)const
{
	double A2 = pow(point2.x-point0.x,2)+pow(point2.y-point0.y,2)+pow(point2.z-point0.z,2);
	double B2 = pow(point0.x-point1.x,2)+pow(point0.y-point1.y,2)+pow(point0.z-point1.z,2);
	double C2 = pow(point2.x-point1.x,2)+pow(point2.y-point1.y,2)+pow(point2.z-point1.z,2);
	
	return (B2+C2-A2)/(2*sqrt(B2)*sqrt(C2));
}
//---------------------------------------------------------------------------
//�ֱ���RVE�ĳ�����ϰ����ȷֲ�ѡȡ���׹ܵ���ʼ��
int GeoNano::Get_seed_point(const struct RVE_Geo &cell, int &seed, Point_3D &point)const
{
	seed = (2053*seed + 13849)%MAX_INT;
	point.x = cell.poi_min.x + seed*cell.len_x/MAX_INT;
    
	seed = (2053*seed + 13849)%MAX_INT;
	point.y = cell.poi_min.y + seed*cell.wid_y/MAX_INT;
    
	seed = (2053*seed + 13849)%MAX_INT;
    if (cnts_geo.type=="CF")
    point.z=0;
    else
	point.z = cell.poi_min.z + seed*cell.hei_z/MAX_INT;
    
	point.flag = 0; //���������׹ܵ���ʼ��
    
	return 1;
}
//---------------------------------------------------------------------------
//�ֱ���RVE�ĳ�����ϰ����ȷֲ�ѡȡ���׹ܵ���ʼ��
int GeoNano::Get_seed_point_outside_clusters(const struct RVE_Geo &cell, int &seed, Point_3D &point, const vector<struct elliparam> &ellips)const
{
	int cnt_seed_count =0; //��¼���˶��ٴ����׹�����
	do
	{
		//-------------------------------------------------------------
		//��������
		seed = (2053*seed + 13849)%MAX_INT;
		point.x = cell.poi_min.x + seed*cell.len_x/MAX_INT;
        
		seed = (2053*seed + 13849)%MAX_INT;
		point.y = cell.poi_min.y + seed*cell.wid_y/MAX_INT;
        
		seed = (2053*seed + 13849)%MAX_INT;
		point.z = cell.poi_min.z + seed*cell.hei_z/MAX_INT;
        
		point.flag = 0; //���������׹ܵ���ʼ��
        
		cnt_seed_count++;					//��¼ѡ��һ��
		if(cnt_seed_count>100000) { hout << "���׹������ڵ����ڵ������Ŵ�������ʮ��Σ���û�ҵ����ʵ����ӵ㣡 ���飡" << endl; return 0; }
        
	}while(Judge_ellipses_including_point(ellips, point)>=0); //�ж��Ƿ���������,-1��ʾ��������������
    
	return 1;
}
//---------------------------------------------------------------------------
//�ֱ���������ĳ��ж���2��Ϊ����ߵ��������ϰ����ȷֲ�ѡȡ���׹ܵ���ʼ��
//���жϴ˵��Ƿ��������ڣ�����������ڣ��ٱ任����������ϵ��
int GeoNano::Get_seed_point_inside_clusters(const struct elliparam ellip, int &seed, Point_3D &point)const
{
	do
	{
		//-------------------------------------------------------------
		//��������
		seed = (2053*seed + 13849)%MAX_INT;
		point.x = -ellip.a + 2*seed*ellip.a/MAX_INT;
        
		seed = (2053*seed + 13849)%MAX_INT;
		point.y = -ellip.b + 2*seed*ellip.b/MAX_INT;
        
		seed = (2053*seed + 13849)%MAX_INT;
		point.z = -ellip.c + 2*seed*ellip.c/MAX_INT;
        
	}while(pow(point.x,2)/pow(ellip.a, 2)+pow(point.y,2)/pow(ellip.b, 2)+pow(point.z,2)/pow(ellip.c, 2)-1.0>=-Zero);  //�������ⲿҪ��������
    
	//-------------------------------------------------------------
	//ת�������ڵ���������ϵ��
	const double x = point.x;
	const double y = point.y;
	const double z = point.z;
    
	point.x = ellip.x+x*ellip.alpha1+y*ellip.alpha2+z*ellip.alpha3;
	point.y = ellip.y+x*ellip.beta1+y*ellip.beta2+z*ellip.beta3;
	point.z = ellip.z+x*ellip.gamma1+y*ellip.gamma2+z*ellip.gamma3;
    
	point.flag = 0; //���������׹ܵ���ʼ��
    
	return 1;
}
//---------------------------------------------------------------------------
//���ֲ��������ѡȡ���ֵ
int GeoNano::Get_random_value(const string &dist_type, const double &min, const double &max, int &seed, double &value)const
{
	if(min>max) { hout << "�������Сֵ�������ֵ(Get_random_value)�� ���飡" << endl; return 0; }
    
	if(dist_type=="uniform")
	{
		seed = (2053*seed + 13849)%MAX_INT;
		value = seed*(max-min)/MAX_INT + min;
	}
	else if(dist_type=="normal")
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
//�����������ϰ����ȷֲ����ѡȡһ��������Ϊ���׹��ߵ��׸�����
int GeoNano::Get_uniform_direction(const struct CNT_Geo &cnts_geo, int &seed_sita, int &seed_pha, double &cnt_sita, double &cnt_pha)const
{
	if(cnts_geo.dir_dist_type=="random")
	{
		//sita��(0, PI)�Ͼ��ȷֲ�
		seed_sita = (2053*seed_sita + 13849)%MAX_INT;
		cnt_sita = seed_sita*PI/MAX_INT;
        
		//pha��(0, 2PI)�Ͼ��ȷֲ�
		seed_pha = (2053*seed_pha + 13849)%MAX_INT;
		cnt_pha = 2.0*seed_pha*PI/MAX_INT;
	}
	else if(cnts_geo.dir_dist_type=="specific")
	{
		//�ض��ĳ�ʼ�����
		seed_sita = (2053*seed_sita + 13849)%MAX_INT;
		seed_pha = (2053*seed_pha + 13849)%MAX_INT;
        
		if((seed_sita+seed_pha)%2==0)
		{
			cnt_sita = cnts_geo.ini_sita; //������
			cnt_pha = cnts_geo.ini_pha; //������
		}
		else
		{
			cnt_sita = PI - cnts_geo.ini_sita; //������
			cnt_pha = PI + cnts_geo.ini_pha; //������
		}
	}
	
	return 1;
}
//---------------------------------------------------------------------------
//����Z������һ�µķ���Ϊ����, ������̬�ֲ�, γ����ȷֲ�, ���ѡȡ����
int GeoNano::Get_normal_direction(const double &omega, int &seed_sita, int &seed_pha, double &cnt_sita, double &cnt_pha)const
{
	//sita��0Ϊ������(-omega, +omega)����̬�ֲ�
	int sum=0;
	for(int i=0; i<12; i++)
	{
		seed_sita = (2053*seed_sita + 13849)%MAX_INT;
		sum += seed_sita;
	}
	cnt_sita = fabs((sum*omega)/(6.0*MAX_INT)-omega);
    
	//pha��(0, 2PI)�Ͼ��ȷֲ�
	seed_pha = (2053*seed_pha + 13849)%MAX_INT;
	cnt_pha = 2.0*seed_pha*PI/MAX_INT;
	
	return 1;
}
//---------------------------------------------------------------------------
//��¼�����
MathMatrix GeoNano::Get_transformation_matrix(const double &sita, const double &pha)const
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
//��¼��������Ǿ���
MathMatrix GeoNano::Get_vector_transformation_matrix(const Point_3D &point0, const Point_3D &point1)const
{
	Point_3D temp_poi(point1.x-point0.x, point1.y-point0.y, point1.z-point0.z); //������
	const double Sx2y2 = sqrt(temp_poi.x*temp_poi.x+temp_poi.y*temp_poi.y);
	const double Sx2y2z2 = sqrt(temp_poi.x*temp_poi.x+temp_poi.y*temp_poi.y+temp_poi.z*temp_poi.z);
	const double cos_sita = temp_poi.x/Sx2y2z2;
	const double sin_sita = Sx2y2/Sx2y2z2;
    
	//cnt_sita�Ǳ任����
	MathMatrix Msita(3,3);
	Msita.element[0][0] = cos_sita;
	Msita.element[0][2] = sin_sita;
	Msita.element[1][1] = 1;
	Msita.element[2][0] = -sin_sita;
	Msita.element[2][2] = cos_sita;
    
	const double cos_pha = temp_poi.x/Sx2y2;
	const double sin_pha = temp_poi.y/Sx2y2;
	//cnt_pha�Ǳ任����
	MathMatrix Mpha(3,3);
	Mpha.element[0][0] = cos_pha;
	Mpha.element[0][1] = -sin_pha;
	Mpha.element[1][0] = sin_pha;
	Mpha.element[1][1] = cos_pha;
	Mpha.element[2][2] = 1;
    
	return Mpha*Msita;
}
//---------------------------------------------------------------------------
//�����½ڵ����꣨����任��
Point_3D GeoNano::Get_new_point(MathMatrix &Matrix, const double &Rad)const
{
	//һά����
	MathMatrix Rvec(3,1);
	Rvec.element[0][0] = 0;
	Rvec.element[1][0] = 0;
	Rvec.element[2][0] = Rad;
    
	//һά����
	MathMatrix Res(3,1);
	Res = Matrix*Rvec;
    
	Point_3D Point(Res.element[0][0], Res.element[1][0], Res.element[2][0]);
    
	return Point;
}
//---------------------------------------------------------------------------
//�жϵ��Ƿ񱻰����ڵ�����
int GeoNano::Judge_cell_including_point(const struct RVE_Geo &cell, const Point_3D &point)const
{
	if(point.x<cell.poi_min.x||point.x>cell.poi_min.x+cell.len_x||
       point.y<cell.poi_min.y||point.y>cell.poi_min.y+cell.wid_y||
       point.z<cell.poi_min.z||point.z>cell.poi_min.z+cell.hei_z) return 0;
    
	return 1;
}
//---------------------------------------------------------------------------
//�жϵ��Ƿ񱻰�����������
int GeoNano::Judge_ellip_including_point(const struct elliparam &ellip, const Point_3D &point)const
{
	const double x=point.x-ellip.x;
	const double y=point.y-ellip.y;
	const double z=point.z-ellip.z;
    
	const double ex = x*ellip.alpha1+y*ellip.beta1+z*ellip.gamma1;
	const double ey = x*ellip.alpha2+y*ellip.beta2+z*ellip.gamma2;
	const double ez = x*ellip.alpha3+y*ellip.beta3+z*ellip.gamma3;
    
	const double ef=pow(ex,2)/pow(ellip.a, 2)+pow(ey,2)/pow(ellip.b, 2)+pow(ez,2)/pow(ellip.c, 2)-1.0;
    
	if(ef>Zero) return 0;
    
	return 1;
}
//---------------------------------------------------------------------------
//�жϵ��Ƿ񱻰�����������(����������ı��, ǰ�����������򲢲��ཻ)
int GeoNano::Judge_ellipses_including_point(const vector<struct elliparam> &ellips, const Point_3D &point)const
{
	for(int i=0; i<(int)ellips.size(); i++)
	{
		if(Judge_ellip_including_point(ellips[i], point)) return i;
	}
    
	return -1;
}
//---------------------------------------------------------------------------
//�����������˵��뵥����������н���(�������˵������֮�䣬����0<t<1)����������t��С����������Щ����
int GeoNano::Get_intersecting_point_RVE_surface(const struct RVE_Geo &cell, Point_3D &point0, Point_3D &point1, vector<Point_3D> &ipoi_vec)const
{
	double t_temp[6];
	//��X�ᴹֱƽ��
	t_temp[0] = (cell.poi_min.x - point0.x)/(point1.x - point0.x);
	t_temp[1] = (cell.poi_min.x + cell.len_x - point0.x)/(point1.x - point0.x);
	//��Y�ᴹֱƽ��
	t_temp[2] = (cell.poi_min.y - point0.y)/(point1.y - point0.y);
	t_temp[3] = (cell.poi_min.y + cell.wid_y - point0.y)/(point1.y - point0.y);
	//��Z�ᴹֱƽ��
	t_temp[4] = (cell.poi_min.z - point0.z)/(point1.z - point0.z);
	t_temp[5] = (cell.poi_min.z + cell.hei_z - point0.z)/(point1.z - point0.z);
    
	vector<double> t_ratio;
	for(int i=0; i<6; i++)
	{
		if(t_temp[i]>=0&&t_temp[i]<1)
		{
			//���ַ�����
			int left = 0;
			int right = (int)t_ratio.size()-1;
			while(right>=left)
			{
				int middle = (left + right)/2;
				if(fabs(t_ratio[middle] - t_temp[i])<Zero) goto T_Value_Same; //ֵ��ͬ�����
				else if(t_ratio[middle] > t_temp[i]) right = middle - 1;
				else left = middle + 1;
			}
			t_ratio.insert(t_ratio.begin()+left, t_temp[i]);	//����
        T_Value_Same: ;
		}
	}
    
	if((int)t_ratio.size()<1||(int)t_ratio.size()>3)
	{
		hout << "�������׹������߶���RVE����Ľ������Ϊ" << (int)t_ratio.size() << ", ����һ�����߶��������� ���飡 " << endl;
		return 0;
	}
    
	Point_3D point_temp;
	for(int i=0; i<(int)t_ratio.size(); i++)
	{
		point_temp = point0+(point1-point0)*t_ratio[i];
		point_temp.flag = 1;		//�м��
		
		//---------------------------------------------------------------------------
		//�������
		if(fabs(point_temp.x-cell.poi_min.x)<Zero) point_temp.x = cell.poi_min.x;
		else if(fabs(point_temp.x-cell.poi_min.x-cell.len_x)<Zero) point_temp.x = cell.poi_min.x + cell.len_x;
        
		if(fabs(point_temp.y-cell.poi_min.y)<Zero) point_temp.y = cell.poi_min.y;
		else if(fabs(point_temp.y-cell.poi_min.y-cell.wid_y)<Zero) point_temp.y = cell.poi_min.y + cell.wid_y;
        
		if(fabs(point_temp.z-cell.poi_min.z)<Zero) point_temp.z = cell.poi_min.z;
		else if(fabs(point_temp.z-cell.poi_min.z-cell.hei_z)<Zero) point_temp.z = cell.poi_min.z + cell.hei_z;
        
		//---------------------------------------------------------------------------
		//�����µ�
		ipoi_vec.push_back(point_temp);
	}
    
	return 1;
}
//---------------------------------------------------------------------------
//����ƽ������任
int GeoNano::Periodical_coordinate_transformation(const struct RVE_Geo &cell_geo, Point_3D point[])const
{
	Point_3D mid_poi = (point[0] + point[1])*0.5;
	if(Judge_cell_including_point(cell_geo, mid_poi)==1)
	{
		hout << "���� RVE��һ����RVE������һ�����ߵ��е���RVE�ڲ������飡 " << endl;
		return 0;
	}
    
	Point_3D opoit(cell_geo.poi_min);
	if(mid_poi.x<cell_geo.poi_min.x) opoit.x -= cell_geo.len_x;
	else if(mid_poi.x>cell_geo.poi_min.x+cell_geo.len_x) opoit.x += cell_geo.len_x;
    
	if(mid_poi.y<cell_geo.poi_min.y) opoit.y -= cell_geo.wid_y;
	else if(mid_poi.y>cell_geo.poi_min.y+cell_geo.wid_y) opoit.y += cell_geo.wid_y;
    
	if(mid_poi.z<cell_geo.poi_min.z) opoit.z -= cell_geo.hei_z;
	else if(mid_poi.z>cell_geo.poi_min.z+cell_geo.hei_z) opoit.z += cell_geo.hei_z;
    
	//����任
	point[0] = point[0] - opoit;
	point[0].flag = 0;
    
	point[1] = point[1] - opoit;
	point[1].flag = 1;
    
	if(Judge_cell_including_point(cell_geo, point[0])==0||Judge_cell_including_point(cell_geo, point[1])==0)
	{
		hout << "���� ����任��ĵ���Ȼ����RVE�ڲ������飡 " << endl;
		return 0;
	}
    
	return 1;
}
//---------------------------------------------------------------------------
//������׹����缸�ι���(��Tecplot)
int GeoNano::Export_cnt_networks_threads(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points)const
{
	ofstream otec("CNT_Wires.dat");
	otec << "TITLE = CNT_Wires" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;
    
	//---------------------------------------------------------------------------
	//����������
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
	//������׹�����ά����
	for(int i=0; i<(int)cnts_points.size(); i++)
	{
		otec << "ZONE T=\"Line\"" << endl;
		otec << "i=1," << "j=" << (int)cnts_points[i].size() << ", f=point" << endl;
		for (int j=0; j<(int)cnts_points[i].size(); j++)
		{
			otec << cnts_points[i][j].x << "  " << cnts_points[i][j].y << "  " << cnts_points[i][j].z << endl;
		}
		otec << endl << endl;
	}
    
	otec.close();
    
	return 1;
}
//---------------------------------------------------------------------------
//����������Լ�Ϳ�㼸�ι���(����������(Ϳ��)��������Tecplot)
int GeoNano::Export_cnt_threads_coating_mesh(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points)
{
	//�������ڱ�ʾ���׹ܵĽڵ㼰����������
	vector<vector<Node> > cnts_nodes;
	vector<vector<Element> > cnts_eles;
    
	if(Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, cnts_points)==0) return 0;
    
	//��������߼�(������)Ϳ����������(�ں���(Export_cnt_threads_coating_mesh))�е���
    //	if(Export_cnts_threads_coating_singlezone(cell, cnts_nodes, cnts_eles, cnts_points)==0) return 0;
    
	//��������߼�(������)Ϳ����������(�ں���(Export_cnt_threads_coating_mesh))�е���
	if(Export_cnts_threads_coating_multizones(cell, cnts_nodes, cnts_eles, cnts_points)==0) return 0;
    
	return 1;
}
//---------------------------------------------------------------------------
//��������߼�(������)Ϳ����������(�ں���(Export_cnt_threads_coating_mesh))�е���
int GeoNano::Export_cnts_threads_coating_singlezone(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points)const
{
	//���׹��ߵ�����
	int cnts_account = (int)nodes.size();
	if(cnts_account!=(int)eles.size()) { hout << "�ڵ�͵�Ԫ��ʾ�����׹���������һ�£� ���飡" << endl; return 0; }
	
	ofstream otec("CNT_Threads_Coating_Meshes_Singlezone.dat");
	otec << "TITLE = CNT_Threads_Coating_Meshes_Singlezone" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;
	
	//---------------------------------------------------------------------------
	//����������
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
	//������׹�����
    
	int nodes_num = 0;
	int eles_num = 0;
    
	for(int i=0; i<cnts_account; i++)
	{
		nodes_num +=  (int)nodes[i].size();
		eles_num += (int)eles[i].size();
	}
    
	otec << "ZONE N=" << nodes_num << ", E=" << eles_num << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for(int i=0; i<cnts_account; i++)
	{
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
	}
	otec << endl;
    
	nodes_num = 0;
	for(int i=0; i<cnts_account; i++)
	{
		if(i!=0)  nodes_num +=  (int)nodes[i-1].size();
		for(int j=0; j<(int)eles[i].size(); j++)
		{
			otec	<< eles[i][j].nodes_id[0]+1+nodes_num << "  " << eles[i][j].nodes_id[1]+1+nodes_num << "  "
            << eles[i][j].nodes_id[2]+1+nodes_num << "  " << eles[i][j].nodes_id[3]+1+nodes_num << endl;
		}
	}
    
	//---------------------------------------------------------------------------
	//������׹�����ά����
	for(int i=0; i<(int)cnts_points.size(); i++)
	{
		otec << "ZONE T=\"Line\"" << endl;
		otec << "i=1," << "j=" << (int)cnts_points[i].size() << ", f=point" << endl;
		for (int j=0; j<(int)cnts_points[i].size(); j++)
		{
			otec << cnts_points[i][j].x << "  " << cnts_points[i][j].y << "  " << cnts_points[i][j].z << endl;
		}
		otec << endl << endl;
	}
    
	otec.close();
    
	return 1;
}
//---------------------------------------------------------------------------
//��������߼�(������)Ϳ����������(�ں���(Export_cnt_threads_coating_mesh))�е���
int GeoNano::Export_cnts_threads_coating_multizones(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points)const
{
	//���׹��ߵ�����
	int cnts_account = (int)nodes.size();
	if(cnts_account!=(int)eles.size()) { hout << "�ڵ�͵�Ԫ��ʾ�����׹���������һ�£� ���飡" << endl; return 0; }
	
	ofstream otec("CNT_Threads_Coating_Meshes_Multizones.dat");
	otec << "TITLE = CNT_Threads_Coating_Meshes_Multizones" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;
	
	//---------------------------------------------------------------------------
	//����������
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
	//������׹�����
	for(int i=0; i<cnts_account; i++)
	{
		otec << "ZONE N=" << (int)nodes[i].size() << ", E=" << (int)eles[i].size() << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
		otec << endl;
		for(int j=0; j<(int)eles[i].size(); j++)
		{
			otec	<< eles[i][j].nodes_id[0]+1 << "  " << eles[i][j].nodes_id[1]+1 << "  "
            << eles[i][j].nodes_id[2]+1 << "  " << eles[i][j].nodes_id[3]+1 << endl;
		}
		otec << endl << endl;
	}
    
	//---------------------------------------------------------------------------
	//������׹�����ά����
	for(int i=0; i<(int)cnts_points.size(); i++)
	{
		otec << "ZONE T=\"Line\"" << endl;
		otec << "i=1," << "j=" << (int)cnts_points[i].size() << ", f=point" << endl;
		for (int j=0; j<(int)cnts_points[i].size(); j++)
		{
			otec << cnts_points[i][j].x << "  " << cnts_points[i][j].y << "  " << cnts_points[i][j].z << endl;
		}
		otec << endl << endl;
	}
    
	otec.close();
    
	return 1;
}
//---------------------------------------------------------------------------
//������׹����缸��(Template for Sofiane)
int GeoNano::Export_cnt_fibers(const vector<vector<Point_3D> > &cnts_points)const
{
	ofstream otec("CNT_Fibers.dat");
    
	//---------------------------------------------------------------------------
	//������׹�����ά����
	for(int i=0; i<(int)cnts_points.size(); i++)
		for (int j=0; j<(int)cnts_points[i].size(); j++)
		{
			otec << i+1 << "  " << cnts_points[i][j].x << "  " << cnts_points[i][j].y << "  " << cnts_points[i][j].z << endl;
		}
    
	otec.close();
    
	return 1;
}
//---------------------------------------------------------------------------
//������׹����缸�ι���(����������Tecplot)
int GeoNano::Export_cnt_networks_meshes(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points, const vector<struct elliparam> &ellips)
{
	//�������ڱ�ʾ���׹ܵĽڵ㼰����������
	vector<vector<Node> > cnts_nodes;
	vector<vector<Element> > cnts_eles;
    
	if(Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, cnts_points)==0) return 0;
    
	//������׹��������Zones in Tecplot
	if(Export_cnts_meshes_multizones(cell, cnts_nodes, cnts_eles)==0) return 0;
    
	//������׹�������Zone in Tecplot
    //	if(Export_cnts_meshes_singlezone(cell, cnts_nodes, cnts_eles)==0) return 0;
    
	//������׹�������Zone in Tecplot�Լ���������
    //	if(Export_cnts_ellipsoids_meshes(cell, cnts_nodes, cnts_eles, ellips)==0) return 0;
    
	//�ֱ�������Ŵ�����Ŵ��ڵ����׹���
    //	if(Export_cnts_discerned_clusters(cell, cnts_nodes, cnts_eles, ellips)==0) return 0;
    
	return 1;
}
//---------------------------------------------------------------------------
//������׹��Լ�Ϳ�㼸�ι���(�����嵥Ԫ(���׹�)�����嵥Ԫ(Ϳ��))
int GeoNano::Export_cnt_coating_meshes(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points)
{
	//�������ڱ�ʾ���׹ܵĽڵ㼰����������
	vector<vector<Node> > cnts_nodes;
	vector<vector<Element> > cnts_eles;
	//�������ڱ�ʾͿ��Ľڵ㼰����������
	vector<vector<Node> > coating_nodes;
	vector<vector<Element> > coating_eles;
    
	//���׹�Ϳ��ĺ��
	const double coating_thickness = 0.005;
	if(Generate_cnts_coating_nodes_elements(cnts_nodes, cnts_eles, coating_nodes, coating_eles, cnts_points, coating_thickness)==0) return 0;
    
	//������׹��Լ�Ϳ������
	if(Export_cnts_coating_zones(cell, cnts_nodes, cnts_eles, coating_nodes, coating_eles)==0) return 0;
    
	return 1;
}
//---------------------------------------------------------------------------
//�������ڱ�ʾ���׹ܵĽڵ㼰����������
int GeoNano::Generate_cnts_coating_nodes_elements(vector<vector<Node> > &cnts_nodes, vector<vector<Element> > &cnts_eles, vector<vector<Node> > &coating_nodes,
                                                  vector<vector<Element> > &coating_eles, const vector<vector<Point_3D> > &cnts_points, const double &coating_thickness)
{
	//ѭ�����׹���
	for(int i=0; i<(int)cnts_points.size(); i++)
	{
		vector<Node> cnts_nod_temp;
		vector<Element> cnts_ele_temp;
		vector<Node> coating_total_nodes;
		vector<Node> coating_nod_temp;
		vector<Element> coating_ele_temp;
        
		const int cps = (int)cnts_points[i].size();
		for(int j=0; j<cps; j++)
		{
			//�õ�ƽ��ķ�������
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
                
				//����任��ԭ����λ��
				double x, y, z;
				x=vect[1].x+tt*(vect[2].x-vect[1].x);
				y=vect[1].y+tt*(vect[2].y-vect[1].y);
				z=vect[1].z+tt*(vect[2].z-vect[1].z);
                
				plane_normal.x = vect[0].x - x;
				plane_normal.y = vect[0].y - y;
				plane_normal.z = vect[0].z - z;
			}
			//�õ�ƽ���ϵĵ�
			Point_3D plane_center = cnts_points[i][j];
            
			//����Բ�Ϸֵķ���
			const int num_sec = 36;
			if(j==0)
			{
				double normal_sita, normal_pha;  //�н�
				//�õ�����plane_normal��������ϵ�еļн�
				if(Get_angles_vector_in_spherial_coordinates(plane_normal, normal_sita, normal_pha)==0) return 0;
                
				//�õ��㷨ʸƽ����(�Ե�Ϊ������һ��ֵΪ�뾶)һ��Բ���ϵ�һ���, ��ʸ�Ƕ�(��������)��֪
				if(Get_points_circle_in_plane(plane_center, normal_sita, normal_pha, cnts_radius[i], num_sec, cnts_nod_temp)==0) return 0;	//��¼cnt�߽����ϵĵ�
				if(Get_points_circle_in_plane(plane_center, normal_sita, normal_pha, cnts_radius[i]+coating_thickness, num_sec, coating_nod_temp)==0) return 0;	//��¼coating�߽����ϵĵ�
			}
			else
			{
				//����ǰһ��Բ�ϵĵ���ǰһ���߶η���(line_vec)�ڹ���plane_center��ʸplane_normal��ƽ���ϵ�ͶӰ��
				Point_3D line_vec;
				line_vec.x = cnts_points[i][j-1].x - cnts_points[i][j].x;
				line_vec.y = cnts_points[i][j-1].y - cnts_points[i][j].y;
				line_vec.z = cnts_points[i][j-1].z - cnts_points[i][j].z;
				if(Get_projected_points_in_plane(plane_center, plane_normal, line_vec, num_sec, cnts_nod_temp)==0) return 0;			//��¼cnt�߽����ϵĵ�
				if(Get_projected_points_in_plane(plane_center, plane_normal, line_vec, num_sec, coating_nod_temp)==0) return 0;	//��¼coating�߽����ϵĵ�
			}
			
			//���뵽Ϳ��ڵ�����(ȥ�����ĵ�)�������еĵ�����͸�ֵ
			int nodsize = (int)cnts_nod_temp.size()-1;
			for(int k=num_sec; k>=0; k--)
			{
				if(k==num_sec) cnts_nod_temp[nodsize-k].type = 0;
				else
				{
					cnts_nod_temp[nodsize-k].type = 1;
					coating_total_nodes.push_back(cnts_nod_temp[nodsize-k]);
				}
			}
			nodsize = (int)coating_nod_temp.size()-1;
			for(int k=num_sec-1; k>=0; k--)
			{
				coating_nod_temp[nodsize-k].type = 2;
				coating_total_nodes.push_back(coating_nod_temp[nodsize-k]);
			}
            
			//����һ�鵥Ԫ����
			if(j!=0)
			{
				//---------------------------------------------------------------
				//���׹ܵĵ�Ԫ
				int nodes_cnt[6];
				nodes_cnt[0] = (j-1)*(num_sec+1);   //���ĵ���
				nodes_cnt[3] = j*(num_sec+1);
				for(int k=1; k<=num_sec; k++)
				{
					nodes_cnt[1] = (j-1)*(num_sec+1) + k;
					nodes_cnt[2] = (j-1)*(num_sec+1) + 1 + k%num_sec;
					nodes_cnt[4] = j*(num_sec+1) + k;
					nodes_cnt[5] = j*(num_sec+1) + 1 + k%num_sec;
                    
					Element eles_cnt[3];
					//----------------------------------------------------------------
					//����ڵ�
					eles_cnt[0].nodes_id.push_back(nodes_cnt[0]);
					eles_cnt[0].nodes_id.push_back(nodes_cnt[1]);
					eles_cnt[0].nodes_id.push_back(nodes_cnt[2]);
					eles_cnt[0].nodes_id.push_back(nodes_cnt[3]);
                    
					eles_cnt[1].nodes_id.push_back(nodes_cnt[1]);
					eles_cnt[1].nodes_id.push_back(nodes_cnt[2]);
					eles_cnt[1].nodes_id.push_back(nodes_cnt[3]);
					eles_cnt[1].nodes_id.push_back(nodes_cnt[5]);
                    
					eles_cnt[2].nodes_id.push_back(nodes_cnt[1]);
					eles_cnt[2].nodes_id.push_back(nodes_cnt[3]);
					eles_cnt[2].nodes_id.push_back(nodes_cnt[4]);
					eles_cnt[2].nodes_id.push_back(nodes_cnt[5]);
					//----------------------------------------------------------------
					//���뵥Ԫ
					cnts_ele_temp.push_back(eles_cnt[0]);
					cnts_ele_temp.push_back(eles_cnt[1]);
					cnts_ele_temp.push_back(eles_cnt[2]);
				}
                
				//---------------------------------------------------------------
				//Ϳ��ĵ�Ԫ
				int nodes_coat[8];
				for(int k=0; k<num_sec; k++)
				{
					nodes_coat[0] = (j-1)*2*num_sec + k;
					nodes_coat[1] = (j-1)*2*num_sec + (k+1)%num_sec;
					nodes_coat[2] = (j-1)*2*num_sec + num_sec + (k+1)%num_sec;
					nodes_coat[3] = (j-1)*2*num_sec + num_sec + k;
                    
					nodes_coat[4] = j*2*num_sec + k;
					nodes_coat[5] = j*2*num_sec + (k+1)%num_sec;
					nodes_coat[6] = j*2*num_sec + num_sec + (1+k)%num_sec;
					nodes_coat[7] = j*2*num_sec + num_sec + k;
                    
					Element eles_coat;
					//----------------------------------------------------------------
					//����ڵ�
					for(int m=0; m<8; m++) eles_coat.nodes_id.push_back(nodes_coat[m]);
                    
					//----------------------------------------------------------------
					//���뵥Ԫ
					coating_ele_temp.push_back(eles_coat);
				}
			}
		}
        
		//����������
		cnts_nodes.push_back(cnts_nod_temp);
		cnts_eles.push_back(cnts_ele_temp);
		coating_nodes.push_back(coating_total_nodes);
		coating_eles.push_back(coating_ele_temp);
	}
    
	return 1;
}
//---------------------------------------------------------------------------
//�������ڱ�ʾ���׹ܵĽڵ㼰����������
int GeoNano::Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points)
{
	//ѭ�����׹���
	for(int i=0; i<(int)cnts_points.size(); i++)
	{
		vector<Node> nod_temp;
		vector<Element> ele_temp;
        
		const int cps = (int)cnts_points[i].size();
		for(int j=0; j<cps; j++)
		{
			//�õ�ƽ��ķ�������
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
                
				//����任��ԭ����λ��
				double x, y, z;
				x=vect[1].x+tt*(vect[2].x-vect[1].x);
				y=vect[1].y+tt*(vect[2].y-vect[1].y);
				z=vect[1].z+tt*(vect[2].z-vect[1].z);
                
				plane_normal.x = vect[0].x - x;
				plane_normal.y = vect[0].y - y;
				plane_normal.z = vect[0].z - z;
			}
			//�õ�ƽ���ϵĵ�
			Point_3D plane_center = cnts_points[i][j];
            
			//����Բ�Ϸֵķ���
			const int num_sec = 36;
			if(j==0)
			{
				double normal_sita, normal_pha;  //�н�
				//�õ�����plane_normal��������ϵ�еļн�
				if(Get_angles_vector_in_spherial_coordinates(plane_normal, normal_sita, normal_pha)==0) return 0;
                
				//�õ��㷨ʸƽ����(�Ե�Ϊ������һ��ֵΪ�뾶)һ��Բ���ϵ�һ���, ��ʸ�Ƕ�(��������)��֪
				if(Get_points_circle_in_plane(plane_center, normal_sita, normal_pha, cnts_radius[i], num_sec, nod_temp)==0) return 0;
			}
			else
			{
				//����ǰһ��Բ�ϵĵ���ǰһ���߶η���(line_vec)�ڹ���plane_center��ʸplane_normal��ƽ���ϵ�ͶӰ��
				Point_3D line_vec;
				line_vec.x = cnts_points[i][j-1].x - cnts_points[i][j].x;
				line_vec.y = cnts_points[i][j-1].y - cnts_points[i][j].y;
				line_vec.z = cnts_points[i][j-1].z - cnts_points[i][j].z;
				if(Get_projected_points_in_plane(plane_center, plane_normal, line_vec, num_sec, nod_temp)==0) return 0;
			}
            
			//����һ�鵥Ԫ����
			if(j!=0)
			{
				int nodes_num[6];
				nodes_num[0] = (j-1)*(num_sec+1);   //���ĵ���
				nodes_num[3] = j*(num_sec+1);
				for(int k=1; k<=num_sec; k++)
				{
					nodes_num[1] = (j-1)*(num_sec+1) + k;
					nodes_num[2] = (j-1)*(num_sec+1) + 1 + k%num_sec;
					nodes_num[4] = j*(num_sec+1) + k;
					nodes_num[5] = j*(num_sec+1) + 1 + k%num_sec;
                    
					Element eles_num[3];
					//----------------------------------------------------------------
					//����ڵ�
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
					//���뵥Ԫ
					ele_temp.push_back(eles_num[0]);
					ele_temp.push_back(eles_num[1]);
					ele_temp.push_back(eles_num[2]);
				}
			}
		}
        
		//����������
		nodes.push_back(nod_temp);
		eles.push_back(ele_temp);
	}
    
	return 1;
}
//---------------------------------------------------------------------------
//�õ�������������ϵ�еļн�
int GeoNano::Get_angles_vector_in_spherial_coordinates(const Point_3D &normal, double &sita, double &pha)const
{
	if(normal.x==0&&normal.y==0&&normal.z==0) { hout << "����������������Ϊ�㣡 ���飡" << endl; return 0; }
	sita =  acos(normal.z/sqrt(normal.x*normal.x+normal.y*normal.y+normal.z*normal.z));
	if(normal.x==0&&normal.y==0) pha = 0;
	else if(normal.y>=0) pha = acos(normal.x/sqrt(normal.x*normal.x+normal.y*normal.y));
	else if(normal.y<0) pha = 2*PI - acos(normal.x/sqrt(normal.x*normal.x+normal.y*normal.y));
    
	return 1;
}
//---------------------------------------------------------------------------
//�õ��㷨ʸƽ����(�Ե�Ϊ������һ��ֵΪ�뾶)һ��Բ���ϵ�һ���, ��ʸ�Ƕ�(��������)��֪
int GeoNano::Get_points_circle_in_plane(const Point_3D &center, const double &trans_sita, const double &trans_pha, const double &radius, const int &num_sec, vector<Node> &nod_temp)const
{
	//�������ĵ�
	Node new_node(center.x, center.y, center.z);
	nod_temp.push_back(new_node);
    
	//����任����
	MathMatrix trans_mat(3,3);
	trans_mat = Get_transformation_matrix(trans_sita, trans_pha);
    
	//һά����
	MathMatrix Rvec(3,1);
	Rvec.element[0][0] = 0;
	Rvec.element[1][0] = 0;
	Rvec.element[2][0] = radius;
    
	//һά����
	MathMatrix Res(3,1);
    
	double sita, pha;
	sita = 0.5*PI; //��xyƽ�����˶�
	for(int i=0; i<num_sec; i++)
	{
		pha = i*2*PI/num_sec;
		MathMatrix matrix_temp = trans_mat*Get_transformation_matrix(sita, pha);
		Res = matrix_temp*Rvec;
        
		new_node.x = center.x + Res.element[0][0];
		new_node.y = center.y + Res.element[1][0];
		new_node.z = center.z + Res.element[2][0];
        
		//����Բ���ϵĵ�
		nod_temp.push_back(new_node);
	}
	
	return 1;
}
//---------------------------------------------------------------------------
//����ǰһ��Բ�ϵĵ���ǰһ���߶η���(line_vec)�ڹ���plane_center��ʸplane_normal��ƽ���ϵ�ͶӰ��
int GeoNano::Get_projected_points_in_plane(const Point_3D &center, const Point_3D &normal, const Point_3D &line, const int &num_sec, vector<Node> &nod_temp)const
{
	//��¼ǰһ�����ɵ��ܽڵ����
	const int nod_size = (int)nod_temp.size();
    
	//�������ĵ�
	Node new_node(center.x, center.y, center.z);
	nod_temp.push_back(new_node);
    
	const double vectors_dot_product = normal.x*line.x+normal.y*line.y+normal.z*line.z;
	if(vectors_dot_product==0.0)
	{
		hout << "������������ֱ��(��Ӧ������0��1��2, ֻ�е���012С��90�ȣ����߶�12>�߶�01ʱ, �ſ��ܳ����������������������ǲ����ܵ�)�� ���飡" << endl;
		return 0;
	}
    
	for(int i=num_sec; i>0; i--)
	{
		Point_3D point(center.x-nod_temp[nod_size-i].x, center.y-nod_temp[nod_size-i].y, center.z-nod_temp[nod_size-i].z);
		new_node.x = nod_temp[nod_size-i].x + (normal.x*point.x+normal.y*point.y+normal.z*point.z)*line.x/ vectors_dot_product;
		new_node.y = nod_temp[nod_size-i].y + (normal.x*point.x+normal.y*point.y+normal.z*point.z)*line.y/ vectors_dot_product;
		new_node.z = nod_temp[nod_size-i].z + (normal.x*point.x+normal.y*point.y+normal.z*point.z)*line.z/ vectors_dot_product;
        
		//����Բ���ϵĵ�
		nod_temp.push_back(new_node);
	}
    
	return 1;
}
//---------------------------------------------------------------------------
//�������׹ܵ��������
int GeoNano:: CNTs_quality_testing(const vector<vector<Point_3D> > &cnts_points)const
{
	//---------------------------------------------------------------------------
	//���׹ܼнǼ��
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
            
			//���������߶�(nod21��nod01)�ļн�
			double cos_ang210 = (x21*x01+y21*y01+z21*z01)/(sqrt(x21*x21+y21*y21+z21*z21)*sqrt(x01*x01+y01*y01+z01*z01));
            
			//���Ƕ�(<=90�Ȳ�>=0�ȣ�cosֵ<=1.0��>=0.0)
			if(cos_ang210<=1.0&&cos_ang210>=0.0)
			{
				hout << "�����׹ܵ��۽ǽǶ�С��PI/2�� ���飡" << endl;
				return 0;
			}
		}
	}
    
	return 1;
}
//---------------------------------------------------------------------------
//������׹��������Zones in Tecplot
int GeoNano::Export_cnts_meshes_multizones(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const
{
	//���׹��ߵ�����
	int cnts_account = (int)nodes.size();
	if(cnts_account!=(int)eles.size()) { hout << "�ڵ�͵�Ԫ��ʾ�����׹���������һ�£� ���飡" << endl; return 0; }
	
	ofstream otec("CNT_Meshes_Multizones.dat");
	otec << "TITLE = CNT_Meshes_Multizones" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;
	
	//---------------------------------------------------------------------------
	//����������
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
	//������׹�����
	for(int i=0; i<cnts_account; i++)
	{
		otec << "ZONE N=" << (int)nodes[i].size() << ", E=" << (int)eles[i].size() << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
		otec << endl;
		for(int j=0; j<(int)eles[i].size(); j++)
		{
			otec	<< eles[i][j].nodes_id[0]+1 << "  " << eles[i][j].nodes_id[1]+1 << "  "
            << eles[i][j].nodes_id[2]+1 << "  " << eles[i][j].nodes_id[3]+1 << endl;
		}
		otec << endl << endl;
	}
    
	otec.close();
    
	return 1;
}
//---------------------------------------------------------------------------
//������׹�������Zone in Tecplot
int GeoNano::Export_cnts_meshes_singlezone(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const
{
	//���׹��ߵ�����
	int cnts_account = (int)nodes.size();
	if(cnts_account!=(int)eles.size()) { hout << "�ڵ�͵�Ԫ��ʾ�����׹���������һ�£� ���飡" << endl; return 0; }
	
	ofstream otec("CNT_Meshes_Singlezone.dat");
	otec << "TITLE = CNT_Meshes_Singlezone" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;
	
	//---------------------------------------------------------------------------
	//����������
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
	//������׹�����
    
	int nodes_num = 0;
	int eles_num = 0;
    
	for(int i=0; i<cnts_account; i++)
	{
		nodes_num +=  (int)nodes[i].size();
		eles_num += (int)eles[i].size();
	}
    
	otec << "ZONE N=" << nodes_num << ", E=" << eles_num << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for(int i=0; i<cnts_account; i++)
	{
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
	}
	otec << endl;
    
	nodes_num = 0;
	for(int i=0; i<cnts_account; i++)
	{
		if(i!=0)  nodes_num +=  (int)nodes[i-1].size();
		for(int j=0; j<(int)eles[i].size(); j++)
		{
			otec	<< eles[i][j].nodes_id[0]+1+nodes_num << "  " << eles[i][j].nodes_id[1]+1+nodes_num << "  "
            << eles[i][j].nodes_id[2]+1+nodes_num << "  " << eles[i][j].nodes_id[3]+1+nodes_num << endl;
		}
	}
    
	otec.close();
    
	return 1;
}
//---------------------------------------------------------------------------
//������׹��Լ�Ϳ������
int GeoNano::Export_cnts_coating_zones(const struct RVE_Geo &cell, const vector<vector<Node> > &cnts_nodes, const vector<vector<Element> > &cnts_eles,
                                       const vector<vector<Node> > &coating_nodes, const vector<vector<Element> > &coating_eles)const
{
	//���׹��ߵ�����
	const int cnts_account = (int)cnts_nodes.size();
	if(cnts_account!=(int)cnts_eles.size()) { hout << "�ڵ�͵�Ԫ��ʾ�����׹���������һ�£� ���飡" << endl; return 0; }
	const int coating_account = (int)coating_nodes.size();
	if(coating_account!=(int)coating_eles.size()) { hout << "�ڵ�͵�Ԫ��ʾ��Ϳ��������һ�£� ���飡" << endl; return 0; }
    
	ofstream otec("CNT_Coating_Meshes.dat");
	otec << "TITLE = CNT_Coating_Meshes" << endl;
	otec << "VARIABLES = X, Y, Z, U" << endl;
	
	//---------------------------------------------------------------------------
	//����������
	otec << "ZONE N=" << 8 << ", E=" << 1 << ", F=FEPOINT, ET=BRICK" << endl;
	double cell_x[2] = {cell.poi_min.x, cell.poi_min.x+cell.len_x};
	double cell_y[2] = {cell.poi_min.y, cell.poi_min.y+cell.wid_y};
	double cell_z[2] = {cell.poi_min.z, cell.poi_min.z+cell.hei_z};
	for(int i=0; i<2; i++)
		for(int j=0; j<2; j++)
			for(int k=0; k<2; k++)
			{
				otec << cell_x[i] << "  " << cell_y[j] << "  " << cell_z[k] << "  " << 3 << endl;
			}
	
	otec << "1 2 4 3 5 6 8 7" << endl;
	otec << endl << endl;
    
	//---------------------------------------------------------------------------
	//������׹�����
	int nodes_num = 0;
	int eles_num = 0;
    
	for(int i=0; i<cnts_account; i++)
	{
		nodes_num +=  (int)cnts_nodes[i].size();
		eles_num += (int)cnts_eles[i].size();
	}
    
	otec << "ZONE N=" << nodes_num << ", E=" << eles_num << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for(int i=0; i<cnts_account; i++)
	{
		for(int j=0; j<(int)cnts_nodes[i].size(); j++)
		{
			otec << cnts_nodes[i][j].x << "  " << cnts_nodes[i][j].y << "  " << cnts_nodes[i][j].z << "  " << cnts_nodes[i][j].type<< endl;
		}
	}
	otec << endl;
    
	nodes_num = 0;
	for(int i=0; i<cnts_account; i++)
	{
		if(i!=0)  nodes_num +=  (int)cnts_nodes[i-1].size();
		for(int j=0; j<(int)cnts_eles[i].size(); j++)
		{
			otec	<< cnts_eles[i][j].nodes_id[0]+1+nodes_num << "  " << cnts_eles[i][j].nodes_id[1]+1+nodes_num << "  "
            << cnts_eles[i][j].nodes_id[2]+1+nodes_num << "  " << cnts_eles[i][j].nodes_id[3]+1+nodes_num << endl;
		}
	}
    
	//---------------------------------------------------------------------------
	//���Ϳ������
	nodes_num = 0;
	eles_num = 0;
    
	for(int i=0; i<coating_account; i++)
	{
		nodes_num +=  (int)coating_nodes[i].size();
		eles_num += (int)coating_eles[i].size();
	}
    
	otec << "ZONE N=" << nodes_num << ", E=" << eles_num << ", F=FEPOINT, ET=BRICK" << endl;
	for(int i=0; i<coating_account; i++)
	{
		for(int j=0; j<(int)coating_nodes[i].size(); j++)
		{
			otec << coating_nodes[i][j].x << "  " << coating_nodes[i][j].y << "  " << coating_nodes[i][j].z << "  " << coating_nodes[i][j].type << endl;
		}
	}
	otec << endl;
    
	nodes_num = 0;
	for(int i=0; i<coating_account; i++)
	{
		if(i!=0)  nodes_num +=  (int)coating_nodes[i-1].size();
		for(int j=0; j<(int)coating_eles[i].size(); j++)
		{
			otec	<< coating_eles[i][j].nodes_id[0]+1+nodes_num << "  " << coating_eles[i][j].nodes_id[1]+1+nodes_num << "  "
            << coating_eles[i][j].nodes_id[2]+1+nodes_num << "  " << coating_eles[i][j].nodes_id[3]+1+nodes_num << "  "
            << coating_eles[i][j].nodes_id[4]+1+nodes_num << "  " << coating_eles[i][j].nodes_id[5]+1+nodes_num << "  "
            << coating_eles[i][j].nodes_id[6]+1+nodes_num << "  " << coating_eles[i][j].nodes_id[7]+1+nodes_num << endl;
		}
	}
    
	otec.close();
    
	return 1;
}
//---------------------------------------------------------------------------
//������׹��߼����Ŵ���������
int GeoNano::Export_cnts_ellipsoids_meshes(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const vector<struct elliparam> &ellips)const
{
	//���׹��ߵ�����
	int cnts_account = (int)nodes.size();
	if(cnts_account!=(int)eles.size()) { hout << "�ڵ�͵�Ԫ��ʾ�����׹���������һ�£� ���飡" << endl; return 0; }
	
	ofstream otec("CNTs_Cluster_Ellipsoids_Meshes.dat");
	otec << "TITLE = CNT_Meshes_Singlezone" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;
	
	//---------------------------------------------------------------------------
	//����������
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
	//������׹�����
    
	int nodes_num = 0;
	int eles_num = 0;
    
	for(int i=0; i<cnts_account; i++)
	{
		nodes_num +=  (int)nodes[i].size();
		eles_num += (int)eles[i].size();
	}
    
	otec << "ZONE N=" << nodes_num << ", E=" << eles_num << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for(int i=0; i<cnts_account; i++)
	{
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
	}
	otec << endl;
    
	nodes_num = 0;
	for(int i=0; i<cnts_account; i++)
	{
		if(i!=0)  nodes_num +=  (int)nodes[i-1].size();
		for(int j=0; j<(int)eles[i].size(); j++)
		{
			otec	<< eles[i][j].nodes_id[0]+1+nodes_num << "  " << eles[i][j].nodes_id[1]+1+nodes_num << "  "
            << eles[i][j].nodes_id[2]+1+nodes_num << "  " << eles[i][j].nodes_id[3]+1+nodes_num << endl;
		}
	}
    
	//---------------------------------------------------------------------------
	//������׹��Ŵ�����������
	for(int i=0; i<(int)ellips.size(); i++)
	{
		const int num_sita = 20;
		const int num_phi = int(2*num_sita*ellips[i].a/ellips[i].c+0.5);  //��������ȡ��
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
    
	return 1;
}
//---------------------------------------------------------------------------
//�ֱ�������Ŵ�����Ŵ��ڵ����׹���
int GeoNano::Export_cnts_discerned_clusters(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const vector<struct elliparam> &ellips)const
{
	//���׹��ߵ�����
	int cnts_account = (int)nodes.size();
	if(cnts_account!=(int)eles.size()) { hout << "�ڵ�͵�Ԫ��ʾ�����׹���������һ�£� ���飡" << endl; return 0; }
	
	ofstream otec("CNTs_Discerned_Clusters.dat");
	otec << "TITLE = CNT_Meshes_Doublezone" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;
	
	//---------------------------------------------------------------------------
	//����������
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
	//���㲢�ֱ�洢�������ڻ���������׹�����Ԫ
	vector<Element> ele_clus[2];  //�ֱ�洢�������ڻ�������׹�����Ԫ
	Point_3D poi_clus;					 //����Ԫ�����ĵ�λ��
	int nodes_num = 0;					 //��¼��Ԫ�ڵ��ʵ��λ��
	for(int i=0; i<cnts_account; i++)
	{
		if(i!=0) nodes_num += (int)nodes[i-1].size();
		for(int j=0; j<(int)eles[i].size(); j++)
		{
			poi_clus.x = 0.0;
			poi_clus.y = 0.0;
			poi_clus.z = 0.0;
            
			Element ele_temp;
			for(int k=0; k<4; k++)
			{
				ele_temp.nodes_id.push_back(eles[i][j].nodes_id[k] + nodes_num);  //��¼�˵�Ԫ
                
				poi_clus.x += nodes[i][eles[i][j].nodes_id[k]].x;
				poi_clus.y += nodes[i][eles[i][j].nodes_id[k]].y;
				poi_clus.z += nodes[i][eles[i][j].nodes_id[k]].z;
			}
            
			poi_clus.x = poi_clus.x/4.0;
			poi_clus.y = poi_clus.y/4.0;
			poi_clus.z = poi_clus.z/4.0;
            
			//�жϴ˵��Ƿ���ĳ��������
			double x1, y1, z1;
			double x, y, z;
			double position_value;
			int key = 0;
            
			for(int k=0; k<(int)ellips.size(); k++)
			{
				x1 = poi_clus.x - ellips[k].x;
				y1 = poi_clus.y - ellips[k].y;
				z1 = poi_clus.z - ellips[k].z;
                
				x = x1*ellips[k].alpha1 + y1*ellips[k].beta1 + z1*ellips[k].gamma1;
				y = x1*ellips[k].alpha2 + y1*ellips[k].beta2 + z1*ellips[k].gamma2;
				z = x1*ellips[k].alpha3 + y1*ellips[k].beta3 + z1*ellips[k].gamma3;
                
				position_value = x*x/(ellips[k].a*ellips[k].a) + y*y/(ellips[k].b*ellips[k].b) + z*z/(ellips[k].c*ellips[k].c) - 1.0;
				if(position_value<0.0) key++;
			}
            
			if(key==0)	ele_clus[0].push_back(ele_temp);
			else if(key==1)	ele_clus[1].push_back(ele_temp);
			else { hout << "һ����ͬʱ��" << key << "�������ڲ��� ���飡" << endl; return 0; }
		}
	}
    
	//---------------------------------------------------------------------------
	//����Ŵ����������׹�
	nodes_num +=  (int)nodes[cnts_account-1].size();
    
	otec << "ZONE N=" << nodes_num << ", E=" << (int)ele_clus[0].size() << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for(int i=0; i<cnts_account; i++)
	{
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
	}
	otec << endl;
    
	for(int i=0; i<(int)ele_clus[0].size(); i++)
	{
		otec	<< ele_clus[0][i].nodes_id[0]+1 << "  " << ele_clus[0][i].nodes_id[1]+1 << "  "
        << ele_clus[0][i].nodes_id[2]+1 << "  " << ele_clus[0][i].nodes_id[3]+1 << endl;
	}
	otec << endl;
    
	//---------------------------------------------------------------------------
	//����Ŵ����������׹�
	otec << "ZONE N=" << nodes_num << ", E=" << (int)ele_clus[1].size() << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for(int i=0; i<cnts_account; i++)
	{
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
	}
	otec << endl;
    
	for(int i=0; i<(int)ele_clus[1].size(); i++)
	{
		otec	<< ele_clus[1][i].nodes_id[0]+1 << "  " << ele_clus[1][i].nodes_id[1]+1 << "  "
        << ele_clus[1][i].nodes_id[2]+1 << "  " << ele_clus[1][i].nodes_id[3]+1 << endl;
	}
    
	otec.close();
    
	return 1;
}
//---------------------------------------------------------------------------
//�������׹��Ŵ����������׹ܵı��غ��������
int GeoNano::Estimate_cnt_cluster_parameters(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points, const vector<struct elliparam> &ellips, const struct CNT_Geo &cnts_geo, struct Clust_Geo &clust_geo)const
{
	double cluster_cnt_length = 0;
	double total_cnt_length = 0;
	double total_cnt_volume = 0;
	double total_cnt_weight = 0;
	//---------------------------------------------------------------------------
	//����ֵ
	clust_geo.cnt_real_volume = 0;
	clust_geo.cnt_real_weight = 0;
	//---------------------------------------------------------------------------
	//�������׹��Ŵ�
	if(clust_geo.vol_fra_criterion!=0&&(int)ellips.size()>0)
	{
		//ѭ���������׹��߶�
		for(int i=0; i<(int)cnts_points.size(); i++)
			for(int j=0; j<(int)cnts_points[i].size()-1; j++)
			{
				double temp_len = cnts_points[i][j].distance_to(cnts_points[i][j+1]);
				total_cnt_length += temp_len;		//���׹ܵ��ܳ���
                
				double temp_vol = temp_len*cnts_radius[i]*cnts_radius[i]*PI;
				total_cnt_volume += temp_vol;		//���׹ܵ������
                
				Point_3D mid_poi;
				mid_poi.x = 0.5*(cnts_points[i][j].x+cnts_points[i][j+1].x);
				mid_poi.y = 0.5*(cnts_points[i][j].y+cnts_points[i][j+1].y);
				mid_poi.z = 0.5*(cnts_points[i][j].z+cnts_points[i][j+1].z);
				
				double x, y, z, x1, y1, z1, f0;
				for(int k=0; k<(int)ellips.size(); k++)
				{
					x=mid_poi.x-ellips[k].x;
					y=mid_poi.y-ellips[k].y;
					z=mid_poi.z-ellips[k].z;
                    
					x1=x*ellips[k].alpha1+y*ellips[k].beta1+z*ellips[k].gamma1;
					y1=x*ellips[k].alpha2+y*ellips[k].beta2+z*ellips[k].gamma2;
					z1=x*ellips[k].alpha3+y*ellips[k].beta3+z*ellips[k].gamma3;
                    
					f0=pow(x1,2)/pow(ellips[k].a, 2)+pow(y1,2)/pow(ellips[k].b, 2)+pow(z1,2)/pow(ellips[k].c, 2)-1.0;
                    
					//�ڴ������ڲ�
					if(f0<=Zero) 
					{ 
						cluster_cnt_length += temp_len;				 //�Ŵ����������׹ܵĳ���
						clust_geo.cnt_real_volume += temp_vol;   //�Ŵ����������׹ܵ����
						break; 
					}
				}
			}
		
		//�Ŵ����������׹ܵ�����
		clust_geo.cnt_real_weight = cluster_cnt_length*cnts_geo.linear_density;
		total_cnt_weight = total_cnt_length*cnts_geo.linear_density;
	}
    
	if(cnts_geo.criterion=="vol")
	{
		hout << "    ���׹��Ŵ����������׹���ռ���������" << clust_geo.cnt_real_volume << endl;
		hout << "    �������������׹���ռ���������" << total_cnt_volume << endl;
	}
	else if(cnts_geo.criterion=="wt"||cnts_geo.criterion=="nwt")
	{
		hout << "    ���׹��Ŵ����������׹���ռ���أ�" << clust_geo.cnt_real_weight/(cell_geo.volume*cell_geo.density) << endl;
		hout << "    �������������׹���ռ���أ�" << total_cnt_weight/(cell_geo.volume*cell_geo.density) << endl;
		hout << "    �Ŵ����������׹�ռ�������������׹ܱ�����" << clust_geo.cnt_real_weight/total_cnt_weight << endl;
	}
    
	return 1;
}
//---------------------------------------------------------------------------
//��¼���׹ܵ���Ϣ(CNTs points��flag==0��ʾ���׹ܵ���ʼ��; flag>0��ʾ���׹��ڵĵ�, �������α��)
int GeoNano::Record_cnt_points_information(const vector<vector<Point_3D> > &cnts_points)
{
    clock_t ct_begin,ct_end;
	ct_begin = clock();
    for(int i=0; i<(int)cnts_points.size(); i++)
    {
        for(int j=0; j<(int)cnts_points[i].size(); j++)
        {
            cnps.push_back(cnts_points[i][j]);
            cnps.back().flag = i;
        }
    }
    ct_end = clock();
    hout << "Record_cnt_points_information time: "<<(double)(ct_end-ct_begin)/CLOCKS_PER_SEC<<"�롣"<<endl;
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//����һ����Ϣ��������ע���У���"%"��ͷ����
string GeoNano::Get_Line(ifstream &infile)const
{
	string s;
	//������Ϣһ��
	getline(infile,s);
	//����ע����     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===========================================================================
