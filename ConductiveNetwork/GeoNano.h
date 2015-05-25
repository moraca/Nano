//===========================================================================
// GeoNano.h
// NanoComposites微结构几何建模头文件
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

const int MAX_INT = 65536; //2^16 用于取随机数
const int N_times=1000;    //用于生成椭球的最大次数

//---------------------------------------------------------------------------
//记录RVE单胞几何信息
struct RVE_Geo
{	
	Point_3D poi_min;	//RVE单胞的坐标最小值点（其各个方向的坐标都是最小值）
	double len_x, wid_y, hei_z;	//长、宽和高
	double volume;	//体积
	double density;  //单胞（纯基体polymer）的密度
	double delt_x, delt_y, delt_z; //用于记录单胞中网格细分的尺度（因为方便传递，定义在这里）
};
//---------------------------------------------------------------------------
//记录CNT的几何物理信息
struct CNT_Geo
{	
	string criterion;						//评价标准
	double step_length;				//纳米管生长步长
	string	len_dist_type;			//纳米管长度分布类型
	string	rad_dist_type;			//纳米管半径分布类型
	string	dir_dist_type;
	double ini_sita, ini_pha;		//初始生长角度(用于specific类型)，球面坐标
	double angle_max;					//每步生长角度正态分布的最大取值范围
	double len_min, len_max;	//纳米管长度范围
	double rad_min, rad_max;	//纳米管半径范围
	double volume_fraction;		//纳米管所占体积分数
	double real_volume;				//纳米管实际体积
	double weight_fraction;        //纳米管所占重量分数
	double real_weight;				//纳米管实际重量
	double linear_density;          //纳米管的线密度
    //-----Lakshmi
    string type;
    //-----Lakshmi
};
//---------------------------------------------------------------------------
//记录纳米管团簇椭球参数信息
struct Clust_Geo
{
	double wt_fra_cluster;				//纳米管团簇中纳米管的重量分数
	double vol_fra_criterion;			//纳米管团簇椭球体积分数（界限）
	double amin;								//纳米管团簇椭球长轴取值范围最小值
	double amax;								//纳米管团簇椭球长轴取值范围最大值
	double bmin;								//纳米管团簇椭球中轴取值范围最小值
	double cmin;								//纳米管团簇椭球短轴取值范围最小值
	double growth_probability;		//纳米管在团簇椭球中生长的概率
	double real_volume_fraction;	//纳米管团簇椭球的实际体积
	double cnt_real_weight;			//纳米管团簇椭球中所包含的纳米管实际比重
	double cnt_real_volume;			//纳米管团簇椭球中所包含的纳米管实际体积
};
//---------------------------------------------------------------------------
//记录椭球参数信息
struct elliparam		
{	
	double x, y, z;		//椭球的中心点
	double a, b, c;		//椭球的长、中和短轴
	double alpha1, alpha2, alpha3; //椭球的九个夹角参数
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
    
    //数据成员
    struct RVE_Geo cell_geo;		//RVE单胞几何信息
    struct CNT_Geo cnts_geo;		//CNT的几何信息
    struct Clust_Geo clust_geo;	//纳米管团簇椭球参数信息
    vector<Point_3D> cnps;			//定义三维点向量(用于输出纳米管信息)
    vector<double> cnts_radius;	//定义纳米管网络中每根纳米管的半径
    vector<Point_3D> seed_points;
    //----------AMC
    struct Region_Geo overlap_regions, cnt_regions; //Here I will store the sizes of the regions used for finding contacts and locating CNTs
    double d_vdw, cutoff; //van der Waals distance and cutoff for overlapping
    int penetrating_points;
    vector<vector<long int> > sectioned_domain;
    vector<vector<int> > sectioned_domain_cnt; //variable to store the diffeent points according to their region number
    vector<vector<long int> > global_point_coord, structure;
    //----------AMC
    
    //构造函数
    GeoNano(){};
    
    //成员函数
    //建立几何模型
    int Geometric_modeling(ifstream &infile, const int &samples_count);
    
private:
    
    //成员函数
    
    //读入RVE和CNT的几何数据
    int Import_geometric_data(ifstream &infile, struct RVE_Geo &cell_geo, struct Region_Geo &overlap_regions, struct Region_Geo &cnt_regions, struct CNT_Geo &cnts_geo, struct Clust_Geo &clust_geo, const int &samples_count, vector<vector<long int> > &sectioned_domain, vector<vector<int> > &sectioned_domain_cnt, double &d_vdw, double &cutoff)const;
    //生成纳米管网络
    int Generate_nanotube_networks(const struct RVE_Geo &cell_geo, const struct CNT_Geo &cnts_geo, vector<vector<Point_3D> > &cnts_points, struct Clust_Geo &clust_geo, vector<struct elliparam> &ellips);
    long int Default_region(Point_3D point);
    void Add_to_regions(Point_3D point, long int point_coordinate);
    int Check_penetration(Point_3D &point, long int point_region, double cnt_rad, vector<Point_3D> cnt_new, vector<vector<Point_3D> > cnts, RVE_Geo cell_geo, int &boundary_flag);
    void Move_point(Point_3D &point, vector<double> cutoff, vector<vector<long int> > affected_points, vector<vector<Point_3D> > cnts);
    void Move_point(Point_3D &point, double cutoff, vector<vector<int> > affected_points, vector<vector<Point_3D> > cnts);
    int Check_segment_orientation(Point_3D point, vector<Point_3D> cnt);
    //新算法生成纳米管网络
    int New_generate_nanotube_networks(const struct RVE_Geo &cell_geo, const struct CNT_Geo &cnts_geo, vector<vector<Point_3D> > &cnts_points, struct Clust_Geo &clust_geo, vector<struct elliparam> &ellips);
    //生成纳米管团簇椭球序列
    int Get_ellip_clusters(const struct RVE_Geo &cell, struct Clust_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle, vector<struct elliparam> &ellips, const int &export_mod)const;
    //生成特定位置的圆球团簇序列(避免随机情况时大数量椭球无法生成)
    int Get_specific_sphere_clusters(const struct RVE_Geo &cell, struct Clust_Geo &clust_geo, int &seed_poi, int &seed_axis, int &seed_angle, vector<struct elliparam> &ellips, const int &export_mod)const;
    //输出纳米管团簇椭球网格
    void Export_cluster_ellipsoids_mesh(const struct RVE_Geo &cell, const vector<struct elliparam> &ellips)const;
    //输出纳米管团簇椭球数据
    void Export_cluster_ellipsoids_data(const vector<struct elliparam> &ellips, const double &ellip_ratio)const;
    //计算空间一点与椭球的位置关系，并返回值
    double Estimate_position_point_ellipsoid(const struct elliparam &temp_ellips, const Point_3D &point)const;
    //判断纳米管是否穿越团簇椭球表面
    int Cnt_cross_cluster_surface(const double &growth_probability, const vector<struct elliparam> &ellips, int &ellip_num, int &seed, const Point_3D &point0, const Point_3D &point1)const;
    //判断纳米管是否穿越一个团簇椭球表面
    int Cnt_go_through_a_cluster(const double &growth_probability, const struct elliparam &ellip, int &seed, const Point_3D &point0, const Point_3D &point1)const;
    //求纳米管线与椭球面的交点
    int Get_intersection_cnt_ellipsoid(struct elliparam &temp_ellips, const Point_3D &point0, const Point_3D &point1, Point_3D &intersect_point)const;
    //移动纳米管线段在椭球面外端点到椭球面内的对称点（以切面为对称面）
    int Move_cnt_point_symmetry_point(struct elliparam &temp_ellips, const Point_3D &intersect_point, Point_3D &cnt_poi)const;
    //记录向量方向角矩阵
    MathMatrix Get_vector_transformation_matrix(const Point_3D &point0, const Point_3D &point1)const;
    //计算三个点中中间一点处所夹的夹角的Cos值
    double Cos_angle_in_three_points(const Point_3D &point0, const Point_3D &point1, const Point_3D &point2)const;
    //分别在RVE的长宽高上按均匀分布选取纳米管的起始点
    int Get_seed_point(const struct RVE_Geo &cell, int &seed, Point_3D &point)const;
    //分别在RVE的长宽高上按均匀分布选取纳米管的起始点(但是要在所有团簇颗粒区域以外)
    int Get_seed_point_outside_clusters(const struct RVE_Geo &cell, int &seed, Point_3D &point, const vector<struct elliparam> &ellips)const;
    //分别在以椭球的长中短轴2倍为长宽高的立方体上按均匀分布选取纳米管的起始点
    //并判断此点是否在椭球内，如果在椭球内，再变换到椭球坐标系下
    int Get_seed_point_inside_clusters(const struct elliparam ellip, int &seed, Point_3D &point)const;
    //按分布函数随机选取相关值
    int Get_random_value(const string &dist_type, const double &min, const double &max, int &seed, double &value)const;
    //在球面坐标上按均匀分布随机选取一个方向作为纳米管线的首个方向
    int Get_uniform_direction(const struct CNT_Geo &cnts_geo, int &seed_sita, int &seed_pha, double &cnt_sita, double &cnt_pha)const;
    //按与Z轴正向一致的方向为中心, 径向正态分布, 纬向均匀分布, 随机选取方向
    int Get_normal_direction(const double &omega, int &seed_sita, int &seed_pha, double &cnt_sita, double &cnt_pha)const;
    //记录方向角
    MathMatrix Get_transformation_matrix(const double &sita, const double &pha)const;
    //计算新节点坐标（坐标变换）
    Point_3D Get_new_point(MathMatrix &Matrix, const double &Rad)const;
    //判断点是否被包含在单胞中
    int Judge_cell_including_point(const struct RVE_Geo &cell, const Point_3D &point)const;
    //判断点是否被包含在椭球中
    int Judge_ellip_including_point(const struct elliparam &ellip, const Point_3D &point)const;
    //判断点是否被包含在一组椭球向量中(并返回椭球的编号, 前提是这组椭球并不相交)
    int Judge_ellipses_including_point(const vector<struct elliparam> &ellips, const Point_3D &point)const;
    //计算这两个端点与单胞界面的所有交点(在两个端点的连线之间，参数0<t<1)，并按参数t从小到大排列这些交点
    int Get_intersecting_point_RVE_surface(const struct RVE_Geo &cell, Point_3D &point0, Point_3D &point1, vector<Point_3D> &ipoi_vec)const;
    //周期平移坐标变换
    int Periodical_coordinate_transformation(const struct RVE_Geo &cell_geo, Point_3D point[])const;
    //输出纳米管网络几何构造(线Tecplot)
    int Export_cnt_networks_threads(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points)const;
    //输出纳米管网络几何(Template for Sofiane)
    int Export_cnt_fibers(const vector<vector<Point_3D> > &cnts_points)const;
    //输出纳米管网络几何构造(四面体网格Tecplot)//函数本身不加const都是由于是由于内部调用的函数中有一处两点相减的操作
    int Export_cnt_networks_meshes(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points, const vector<struct elliparam> &ellips);
    //输出纳米线以及涂层几何构造(四面体网格(涂层)和纳米线Tecplot)
    int Export_cnt_threads_coating_mesh(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points);
    //输出纳米线及(单区域)涂层区域网格(在函数(Export_cnt_threads_coating_mesh))中调用
    int Export_cnts_threads_coating_singlezone(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points)const;
    //输出纳米线及(多区域)涂层区域网格(在函数(Export_cnt_threads_coating_mesh))中调用
    int Export_cnts_threads_coating_multizones(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points)const;
    //输出纳米管以及涂层几何构造(四面体单元(纳米管)六面体单元(涂层))
    int Export_cnt_coating_meshes(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points);
    //生成用于表示纳米管的节点及四面体网格
    int Generate_cnts_coating_nodes_elements(vector<vector<Node> > &cnts_nodes, vector<vector<Element> > &cnts_eles, vector<vector<Node> > &coating_nodes,
                                             vector<vector<Element> > &coating_eles, const vector<vector<Point_3D> > &cnts_points, const double &coating_thickness);
    //生成用于表示纳米管的节点及四面体网格//函数本身不加const都是由于是由于内部有一处两点相减的操作
    int Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points);
    //得到向量在球坐标系中的夹角
    int Get_angles_vector_in_spherial_coordinates(const Point_3D &normal, double &sita, double &pha)const;
    //得到点法矢平面上(以点为中心以一定值为半径)一个圆环上的一组点, 法矢角度(球面坐标)已知
    int Get_points_circle_in_plane(const Point_3D &center, const double &sita, const double &pha, const double &radius, const int &num_sec, vector<Node> &nod_temp)const;
    //计算前一个圆上的点沿前一个线段方向(line_vec)在过点plane_center法矢plane_normal的平面上的投影点
    int Get_projected_points_in_plane(const Point_3D &center, const Point_3D &normal, const Point_3D &line, const int &num_sec, vector<Node> &nod_temp)const;
    //生成纳米管的质量检测
    int CNTs_quality_testing(const vector<vector<Point_3D> > &cnts_points)const;
    //输出纳米管线网格多Zones in Tecplot
    int Export_cnts_meshes_multizones(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const;
    //输出纳米管线网格单Zone in Tecplot
    int Export_cnts_meshes_singlezone(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const;
    //输出纳米管以及涂层区域
    int Export_cnts_coating_zones(const struct RVE_Geo &cell, const vector<vector<Node> > &cnts_nodes, const vector<vector<Element> > &cnts_eles,
                                  const vector<vector<Node> > &coating_nodes, const vector<vector<Element> > &coating_eles)const;
    //输出纳米管线及其团簇椭球网格
    int Export_cnts_ellipsoids_meshes(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const vector<struct elliparam> &ellips)const;
    //记录纳米管点信息(CNTs points的flag==0表示纳米管的起始点; flag>0表示纳米管内的点, 按序依次编号)
    int Record_cnt_points_information(const vector<vector<Point_3D> > &cnts_points);
    //分别输出在团簇外和团簇内的纳米管线
    int Export_cnts_discerned_clusters(const struct RVE_Geo &cell, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, const vector<struct elliparam> &ellips)const;
    //计算纳米管团簇椭球中纳米管的比重和体积分数
    int Estimate_cnt_cluster_parameters(const struct RVE_Geo &cell, const vector<vector<Point_3D> > &cnts_points, const vector<struct elliparam> &ellips, const struct CNT_Geo &cnts_geo, struct Clust_Geo &clust_geo)const;
    
    //读入信息一行，跳过注释行（以%开头）
    string Get_Line(ifstream &infile)const;
};
//-------------------------------------------------------
#endif
//===========================================================================
