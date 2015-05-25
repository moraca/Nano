//===========================================================================
// GeoNano.h
// NanoComposites΢�ṹ���ν�ģͷ�ļ�
// A class of nanocomposites cell of geometric modelling
//===========================================================================

#ifndef GEONANO_H
#define GEONANO_H

#include<iomanip>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include "Geometry_3D.h"
#include "MathMatrix.h"
#include "Fem_3D.h"
#include "Gauss.h"
#include "Hns.h"
using namespace hns;

const int MAX_INT = 65536; //2^16 ����ȡ�����
const int N_times=1000;    //�������������������

//---------------------------------------------------------------------------
//��¼RVE����������Ϣ
struct RVE_Geo
{	
	Point_3D poi_min;	//RVE������������Сֵ�㣨�������������궼����Сֵ��
	double len_x, wid_y, hei_z;	//������͸�
	double volume;	//���
	double density;  //������������polymer�����ܶ�
	double delt_x, delt_y, delt_z; //���ڼ�¼����������ϸ�ֵĳ߶ȣ���Ϊ���㴫�ݣ����������
};
//---------------------------------------------------------------------------
//��¼CNT�ļ���������Ϣ
struct CNT_Geo
{	
	string criterion;						//���۱�׼
	double step_length;				//���׹���������
	string	len_dist_type;			//���׹ܳ��ȷֲ�����
	string	rad_dist_type;			//���׹ܰ뾶�ֲ�����
	string	dir_dist_type;
	double ini_sita, ini_pha;		//��ʼ�����Ƕ�(����specific����)����������
	double angle_max;					//ÿ�������Ƕ���̬�ֲ������ȡֵ��Χ
	double len_min, len_max;	//���׹ܳ��ȷ�Χ
	double rad_min, rad_max;	//���׹ܰ뾶��Χ
	double volume_fraction;		//���׹���ռ�������
	double real_volume;				//���׹�ʵ�����
	double weight_fraction;        //���׹���ռ��������
	double real_weight;				//���׹�ʵ������
	double linear_density;          //���׹ܵ����ܶ�
    //-----Lakshmi
    string type;
    //-----Lakshmi
};
//---------------------------------------------------------------------------
//��¼���׹��Ŵ����������Ϣ
struct Clust_Geo
{
	double wt_fra_cluster;				//���׹��Ŵ������׹ܵ���������
	double vol_fra_criterion;			//���׹��Ŵ�����������������ޣ�
	double amin;								//���׹��Ŵ�������ȡֵ��Χ��Сֵ
	double amax;								//���׹��Ŵ�������ȡֵ��Χ���ֵ
	double bmin;								//���׹��Ŵ���������ȡֵ��Χ��Сֵ
	double cmin;								//���׹��Ŵ��������ȡֵ��Χ��Сֵ
	double growth_probability;		//���׹����Ŵ������������ĸ���
	double real_volume_fraction;	//���׹��Ŵ������ʵ�����
	double cnt_real_weight;			//���׹��Ŵ������������������׹�ʵ�ʱ���
	double cnt_real_volume;			//���׹��Ŵ������������������׹�ʵ�����
};
//---------------------------------------------------------------------------
//��¼���������Ϣ
struct elliparam		
{	
	double x, y, z;		//��������ĵ�
	double a, b, c;		//����ĳ����кͶ���
	double alpha1, alpha2, alpha3; //����ľŸ��нǲ���
	double beta1, beta2, beta3;
	double gamma1, gamma2, gamma3;
};
//-------------------------------------------------------
//Defining this type of structure will be more helpful
struct Region_Geo {
    double lx, ly, lz;  //Length of each region
    int secx, secy, secz; //Number of regions on each direction
};
//-------------------------------------------------------
class GeoNano
{
public:
    
    //���ݳ�Ա
    struct RVE_Geo cell_geo;		//RVE����������Ϣ
    struct CNT_Geo cnts_geo;		//CNT�ļ�����Ϣ
    struct Clust_Geo clust_geo;	//���׹��Ŵ����������Ϣ
    vector<Point_3D> cnps;			//������ά������(����������׹���Ϣ)
    vector<double> cnts_radius;	//�������׹�������ÿ�����׹ܵİ뾶
    vector<Point_3D> seed_points;
    //----------AMC
    struct Region_Geo overlap_regions, cnt_regions; //Here I will store the sizes of the regions used for finding contacts and locating CNTs
    double d_vdw, cutoff; //van der Waals distance and cutoff for overlapping
    int penetrating_points;
    vector<vector<long int> > sectioned_domain;
    vector<vector<int> > sectioned_domain_cnt; //variable to store the diffeent points according to their region number
    vector<vector<long int> > global_point_coord, structure;
    //----------AMC
    
    //���캯��
    GeoNano(){};
    
    //��Ա����
    //��������ģ��
    int Geometric_modeling(ifstream &infile, const int &samples_count);
    
private:
    
    //��Ա����
    
    //����RVE��CNT�ļ�������
    int Import_geometric_data(ifstream &infile, struct RVE_Geo &cell_geo, struct Region_Geo &overlap_regions, struct Region_Geo &cnt_regions, struct CNT_Geo &cnts_geo, struct Clust_Geo &clust_geo, const int &samples_count, vector<vector<long int> > &sectioned_domain, vector<vector<int> > &sectioned_domain_cnt, double &d_vdw, double &cutoff)const;
    //�������׹�����
    int Generate_nanotube_networks(const struct RVE_Geo &cell_geo, const struct CNT_Geo &cnts_geo, vector<vector<Point_3D> > &cnts_points, struct Clust_Geo &clust_geo, vector<struct elliparam> &ellips);
    long int Default_region(Point_3D point);
    void Add_to_regions(Point_3D point, long int point_coordinate);
    int Check_penetration(Point_3D &point, long int point_region, double cnt_rad, vector<Point_3D> cnt_new, vector<vector<Point_3D> > cnts, RVE_Geo cell_geo, int &boundary_flag);
    void Move_point(Point_3D &point, vector<double> cutoff, vector<vector<long int> > affected_points, vector<vector<Point_3D> > cnts);
    void Move_point(Point_3D &point, double cutoff, vector<vector<int> > affected_points, vector<vector<Point_3D> > cnts);
    int Check_segment_orientation(Point_3D point, vector<Point_3D> cnt);
    //���㷨�������׹�����
    int New_generate_nanotube_networks(const struct RVE_Geo &cell_geo, const struct CNT_Geo &cnts_geo, vector<vector<Point_3D> > &cnts_points, struct Clust_Geo &clust_geo, vector<struct elliparam> &ellips);
    //�������׹��Ŵ���������
    int Get_ellip_clusters(const struct RVE_Geo &cell, struct Clust_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle, vector<struct elliparam> &ellips, const int &export_mod)const;
    //�����ض�λ�õ�Բ���Ŵ�����(����������ʱ�����������޷�����)
    int Get_specific_sphere_clusters(const struct RVE_Geo &cell, struct Clust_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle, vector<struct elliparam> &ellips, const int &export_mod)const;
    //������׹��Ŵ���������
    void Export_cluster_ellipsoids_mesh(const struct RVE_Geo &cell, const vector<struct elliparam> &ellips)const;
    //������׹��Ŵ���������
    void Export_cluster_ellipsoids_data(const vector<struct elliparam> &ellips, const double &ellip_ratio)const;
    //����ռ�һ���������λ�ù�ϵ��������ֵ
    double Estimate_position_point_ellipsoid(const struct elliparam &temp_ellips, const Point_3D &point)const;
    //�ж����׹��Ƿ�Խ�Ŵ��������
    int Cnt_cross_cluster_surface(const double &growth_probability, const vector<struct elliparam> &ellips, int &ellip_num, int &seed, const Point_3D &point0, const Point_3D &point1)const;
    //�ж����׹��Ƿ�Խһ���Ŵ��������
    int Cnt_go_through_a_cluster(const double &growth_probability, const struct elliparam &ellip, int &seed, const Point_3D &point0, const Point_3D &point1)const;
    //�����׹�����������Ľ���
    int Get_intersection_cnt_ellipsoid(struct elliparam &temp_ellips, const Point_3D &point0, const Point_3D &point1, Point_3D &intersect_point)const;
    //�ƶ����׹��߶�����������˵㵽�������ڵĶԳƵ㣨������Ϊ�Գ��棩
    int Move_cnt_point_symmetry_point(struct elliparam &temp_ellips, const Point_3D &intersect_point, Point_3D &cnt_poi)const;
    //��¼��������Ǿ���
    MathMatrix Get_vector_transformation_matrix(const Point_3D &point0, const Point_3D &point1)const;
    //�������������м�һ�㴦���еļнǵ�Cosֵ
    double Cos_angle_in_three_points(const Point_3D &point0, const Point_3D &point1, const Point_3D &point2)const;
    //�ֱ���RVE�ĳ�����ϰ����ȷֲ�ѡȡ���׹ܵ���ʼ��
    int Get_seed_point(const struct RVE_Geo &cell, int &seed, Point_3D &point)const;
    //�ֱ���RVE�ĳ�����ϰ����ȷֲ�ѡȡ���׹ܵ���ʼ��(����Ҫ�������Ŵؿ�����������)
    int Get_seed_point_outside_clusters(const struct RVE_Geo &cell, int &seed, Point_3D &point, const vector<struct elliparam> &ellips)const;
    //�ֱ���������ĳ��ж���2��Ϊ����ߵ��������ϰ����ȷֲ�ѡȡ���׹ܵ���ʼ��
    //���жϴ˵��Ƿ��������ڣ�����������ڣ��ٱ任����������ϵ��
    int Get_seed_point_inside_clusters(const struct elliparam ellip, int &seed, Point_3D &point)const;
    //���ֲ��������ѡȡ���ֵ
    int Get_random_value(const string &dist_type, const double &min, const double &max, int &seed, double &value)const;
    //�����������ϰ����ȷֲ����ѡȡһ��������Ϊ���׹��ߵ��׸�����
    int Get_uniform_direction(const struct CNT_Geo &cnts_geo, int &seed_sita, int &seed_pha, double &cnt_sita, double &cnt_pha)const;
    //����Z������һ�µķ���Ϊ����, ������̬�ֲ�, γ����ȷֲ�, ���ѡȡ����
    int Get_normal_direction(const double &omega, int &seed_sita, int &seed_pha, double &cnt_sita, double &cnt_pha)const;
    //��¼�����
    MathMatrix Get_transformation_matrix(const double &sita, const double &pha)const;
    //�����½ڵ����꣨����任��
    Point_3D Get_new_point(MathMatrix &Matrix, const double &Rad)const;
    //�жϵ��Ƿ񱻰����ڵ�����
    int Judge_cell_including_point(const struct RVE_Geo &cell, const Point_3D &point)const;
    //�жϵ��Ƿ񱻰�����������
    int Judge_ellip_including_point(const struct elliparam &ellip, const Point_3D &point)const;
    //�жϵ��Ƿ񱻰�����һ������������(����������ı��, ǰ�����������򲢲��ཻ)
    int Judge_ellipses_including_point(const vector<struct elliparam> &ellips, const Point_3D &point)const;
    //�����������˵��뵥����������н���(�������˵������֮�䣬����0<t<1)����������t��С����������Щ����
    int Get_intersecting_point_RVE_surface(const struct RVE_Geo &cell, Point_3D &point0, Point_3D &point1, vector<Point_3D> &ipoi_vec)const;
    //����ƽ������任
    int Periodical_coordinate_transformation(const struct RVE_Geo &cell_geo, Point_3D point[])const;
    //������׹����缸�ι���(��Tecplot)
    int Export_cnt_networks_threads(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points)const;
    //������׹����缸��(Template for Sofiane)
    int Export_cnt_fibers(const vector<vector<Point_3D> > &cnts_points)const;
    //������׹����缸�ι���(����������Tecplot)//����������const���������������ڲ����õĺ�������һ����������Ĳ���
    int Export_cnt_networks_meshes(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points, const vector<struct elliparam> &ellips);
    //����������Լ�Ϳ�㼸�ι���(����������(Ϳ��)��������Tecplot)
    int Export_cnt_threads_coating_mesh(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points);
    //��������߼�(������)Ϳ����������(�ں���(Export_cnt_threads_coating_mesh))�е���
    int Export_cnts_threads_coating_singlezone(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points)const;
    //��������߼�(������)Ϳ����������(�ں���(Export_cnt_threads_coating_mesh))�е���
    int Export_cnts_threads_coating_multizones(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points)const;
    //������׹��Լ�Ϳ�㼸�ι���(�����嵥Ԫ(���׹�)�����嵥Ԫ(Ϳ��))
    int Export_cnt_coating_meshes(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points);
    //�������ڱ�ʾ���׹ܵĽڵ㼰����������
    int Generate_cnts_coating_nodes_elements(vector<vector<Node> > &cnts_nodes, vector<vector<Element> > &cnts_eles, vector<vector<Node> > &coating_nodes,
                                             vector<vector<Element> > &coating_eles, const vector<vector<Point_3D> > &cnts_points, const double &coating_thickness);
    //�������ڱ�ʾ���׹ܵĽڵ㼰����������//����������const���������������ڲ���һ����������Ĳ���
    int Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points);
    //�õ�������������ϵ�еļн�
    int Get_angles_vector_in_spherial_coordinates(const Point_3D &normal, double &sita, double &pha)const;
    //�õ��㷨ʸƽ����(�Ե�Ϊ������һ��ֵΪ�뾶)һ��Բ���ϵ�һ���, ��ʸ�Ƕ�(��������)��֪
    int Get_points_circle_in_plane(const Point_3D &center, const double &sita, const double &pha, const double &radius, const int &num_sec, vector<Node> &nod_temp)const;
    //����ǰһ��Բ�ϵĵ���ǰһ���߶η���(line_vec)�ڹ���plane_center��ʸplane_normal��ƽ���ϵ�ͶӰ��
    int Get_projected_points_in_plane(const Point_3D &center, const Point_3D &normal, const Point_3D &line, const int &num_sec, vector<Node> &nod_temp)const;
    //�������׹ܵ��������
    int CNTs_quality_testing(const vector<vector<Point_3D> > &cnts_points)const;
    //������׹��������Zones in Tecplot
    int Export_cnts_meshes_multizones(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const;
    //������׹�������Zone in Tecplot
    int Export_cnts_meshes_singlezone(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const;
    //������׹��Լ�Ϳ������
    int Export_cnts_coating_zones(const struct RVE_Geo &cell, const vector<vector<Node> > &cnts_nodes, const vector<vector<Element> > &cnts_eles,
                                  const vector<vector<Node> > &coating_nodes, const vector<vector<Element> > &coating_eles)const;
    //������׹��߼����Ŵ���������
    int Export_cnts_ellipsoids_meshes(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const vector<struct elliparam> &ellips)const;
    //��¼���׹ܵ���Ϣ(CNTs points��flag==0��ʾ���׹ܵ���ʼ��; flag>0��ʾ���׹��ڵĵ�, �������α��)
    int Record_cnt_points_information(const vector<vector<Point_3D> > &cnts_points);
    //�ֱ�������Ŵ�����Ŵ��ڵ����׹���
    int Export_cnts_discerned_clusters(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const vector<struct elliparam> &ellips)const;
    //�������׹��Ŵ����������׹ܵı��غ��������
    int Estimate_cnt_cluster_parameters(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points, const vector<struct elliparam> &ellips, const struct CNT_Geo &cnts_geo, struct Clust_Geo &clust_geo)const;
    
    //������Ϣһ�У�����ע���У���%��ͷ��
    string Get_Line(ifstream &infile)const;
};
//-------------------------------------------------------
#endif
//===========================================================================
