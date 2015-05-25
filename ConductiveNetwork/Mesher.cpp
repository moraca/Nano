//===========================================================================
// Mesher.cpp
// 三维网格剖分成员函数
// Member functions in a class of 3D mesh generation
//===========================================================================
#include "Mesher.h"

//---------------------------------------------------------------------------
//网格生成函数；
int Mesher::Generate_mesh(ifstream &infile, struct RVE_Geo &cell_geo, vector<Point_3D> &cnps)
{
	//生成背景六面体网格
	if(Brick_background_mesh(infile, cell_geo)== 0) return 0;
	//一维存储周期边界节点编号及起始点信息
	if(Record_periodic_boundary_nodes_number(cell_geo)==0) return 0;
	//赋予网格单元材料属性
	if(Element_material_property(infile)==0) return 0;
	//记录单元的相关纳米管线段
	if(Element_relative_cnts(cnps, cell_geo)==0) return 0;
	//输出检测单元是否匹配了正确的相关纳米管线段信息
//	if(Test_element_relative_cnts(cnps)== 0) return 0;
	//输出Tecplot可视化网格数据
//	if(Export_mesh_data_tecplot("Mesh_in_tecplot.dat")== 0) return 0;

	return 1;
}
//---------------------------------------------------------------------------
//背景六面体网格生成
int Mesher::Brick_background_mesh(ifstream &infile, struct RVE_Geo &cell_geo)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//读入单胞的长、宽和高的分割细度
	istringstream istr2(Get_Line(infile));	
	istr2 >> dx >> dy >> dz;
	if(dx<0||dy<0||dz<0)
	{
		hout << "RVE长宽高的分割细度输入错误！ 请检查！ ......... Brick_background_mesh" << endl;
		return 0;
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//背景网格(节点和单元)

	//生成各个方向的点数
	int x_sec_num = int(cell_geo.len_x/dx) +1;
	int y_sec_num = int(cell_geo.wid_y/dy)+1;
	int z_sec_num = int(cell_geo.hei_z/dz) +1;

	//微调整剖分步长，保证单胞的尺寸不变
	dx = cell_geo.len_x/(x_sec_num-1);
	dy = cell_geo.wid_y/(y_sec_num-1);
	dz = cell_geo.hei_z/(z_sec_num-1);

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//记录分割细度到单胞信息中
	cell_geo.delt_x = dx;
	cell_geo.delt_y = dy;
	cell_geo.delt_z = dz;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//生成节点
	for( int i=0; i<z_sec_num; i++ )
	{
		double z = cell_geo.poi_min.z + i * dz ;
		for( int j=0; j<y_sec_num; j++ )
		{
			double y = cell_geo.poi_min.y + j * dy ;
			for( int k=0; k<x_sec_num; k++ )
			{
				double x = cell_geo.poi_min.x + k * dx ;

				Node nd(x,y,z);
				nd.type = Deter_node_type(i, j, k, z_sec_num-1, y_sec_num-1, x_sec_num-1);	//标注节点的位置(角点、边界线点、边界面点、内点)
				nodes.push_back(nd);
			}
		}
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//生成六面体单元

	//定义六面体向量
	for( int i=0; i<z_sec_num-1; i++ )
	{
		for( int j=0; j<y_sec_num-1; j++ )
		{
			for( int k=0; k<x_sec_num-1; k++ )
			{
				//六面体的八个顶点
				Element ele_temp;
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + j * x_sec_num + k);
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + j * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k);
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k );
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k);

				elements.push_back(ele_temp);
			}
		}
	}

	return 1;
}
//---------------------------------------------------------------------------
//赋予网格单元材料属性
int Mesher::Element_material_property(ifstream &infile)
{
	istringstream istr(Get_Line(infile));		//跳过注释行
	int define_case;
	double fiber_radius;
	istr >> define_case;

	if(define_case!=0&&define_case!=1) { hout <<"注意：定义的材料单元属性情况错误！请检查....Element_material_property 1" << endl;	return 0; }
	if(define_case==1) istr >> fiber_radius;
	if(fiber_radius<0.0) { hout <<"注意：纤维材料的半径输入错误！请检查....Element_material_property 2" << endl;	return 0; }

	if(define_case==0)	//纯材料
	{
		for(int i=0; i<(int)elements.size(); i++)
		{
			elements[i].type = 381;		//对应单元的形状类型(三个数字xyz： x 维度；y单元的节点个数；z单元的形函数幂次)；
														//例如，121: 一维两节点(线段)线性形函数；231: 二维三节点(三角形)线性形函数；241: 二维四节点(四边形)线性形函数；
														//341: 三维四节点(四面体)线性形函数；361: 三维六节点(三棱柱)线性形函数；381：三维八节点(六面体)线性形函数
			elements[i].mat = 0;			//对应材料属性编号
		}
	}
	else if(define_case==1) //Z轴单向纤维
	{
		for(int i=0; i<(int)elements.size(); i++)
		{
			elements[i].type = 381;	
			Point_3D centp(0,0,0);
			for(int j=0; j<(int)elements[i].nodes_id.size(); j++)
			{
				centp.x += nodes[elements[i].nodes_id[j]].x;
				centp.y += nodes[elements[i].nodes_id[j]].y;
				centp.z += nodes[elements[i].nodes_id[j]].z;
			}
			centp = centp/elements[i].nodes_id.size();
			if((centp.x-0.5)*(centp.x-0.5)+(centp.y-0.5)*(centp.y-0.5)<fiber_radius*fiber_radius) elements[i].mat = 1;
			else elements[i].mat = 0;
		}
	}

	return 1;
}
//---------------------------------------------------------------------------
//记录单元的相关纳米管线段
int Mesher::Element_relative_cnts(vector<Point_3D> &cnps, const struct RVE_Geo &cell_geo)
{
	//生成各个方向的段数(与背景网格的生成部分一致)
	int x_sec_num = int(cell_geo.len_x/dx)+1;
	int y_sec_num = int(cell_geo.wid_y/dy)+1;
	int z_sec_num = int(cell_geo.hei_z/dz)+1;

	for(int i=0; i<(int)cnps.size()-1; i++)
	{
		if(cnps[i+1].flag==0) continue;
		
		int n_temp;
		Point_3D point[2] = {cnps[i], cnps[i+1]};
		
		int nx[2] = { int((point[0].x-cell_geo.poi_min.x)/dx), int((point[1].x-cell_geo.poi_min.x)/dx) };		
		if(nx[0]>nx[1]) { n_temp = nx[0]; nx[0] = nx[1]; nx[1] = n_temp; } //按从小到大
		if(nx[1]==x_sec_num-1) nx[1] -= 1; //最大边界上的点属于前一段内

		int ny[2] = { int((point[0].y-cell_geo.poi_min.y)/dy), int((point[1].y-cell_geo.poi_min.y)/dy) };
		if(ny[0]>ny[1]) { n_temp = ny[0]; ny[0] = ny[1]; ny[1] = n_temp; } //按从小到大
		if(ny[1]==y_sec_num-1) ny[1] -= 1; //最大边界上的点属于前一段内

		int nz[2] = { int((point[0].z-cell_geo.poi_min.z)/dz), int((point[1].z-cell_geo.poi_min.z)/dz) };
		if(nz[0]>nz[1]) { n_temp = nz[0]; nz[0] = nz[1]; nz[1] = n_temp; } //按从小到大	  
		if(nz[1]==z_sec_num-1) nz[1] -= 1; //最大边界上的点属于前一段内

		for(int jz=nz[0]; jz<=nz[1]; jz++)
			for(int jy=ny[0]; jy<=ny[1]; jy++)
				for(int jx=nx[0]; jx<=nx[1]; jx++)
				{
					int elen = jz*(x_sec_num-1)*(y_sec_num-1) + jy*(x_sec_num-1) + jx;

					Point_3D Hexc[2];  //单元的最小值点和最大值点
					Hexc[0].x = nodes[elements[elen].nodes_id[0]].x;
					Hexc[0].y = nodes[elements[elen].nodes_id[0]].y;
					Hexc[0].z = nodes[elements[elen].nodes_id[0]].z;
					Hexc[1].x = nodes[elements[elen].nodes_id[6]].x;
					Hexc[1].y = nodes[elements[elen].nodes_id[6]].y;
					Hexc[1].z = nodes[elements[elen].nodes_id[6]].z;

					//---------------------------------------------------------------------------
					//线段的某个端点在单元内部
					if(point[0].x>=Hexc[0].x&&point[0].x<=Hexc[1].x&&
						point[0].y>=Hexc[0].y&&point[0].y<=Hexc[1].y&&
						point[0].z>=Hexc[0].z&&point[0].z<=Hexc[1].z)
					{
						elements[elen].relative_cnts.push_back(i);
						continue;
					}
					if(point[1].x>=Hexc[0].x&&point[1].x<=Hexc[1].x&&
						point[1].y>=Hexc[0].y&&point[1].y<=Hexc[1].y&&
						point[1].z>=Hexc[0].z&&point[1].z<=Hexc[1].z)
					{
						elements[elen].relative_cnts.push_back(i);
						continue;
					}

					//---------------------------------------------------------------------------
					//线段的端点单元外部但是线段与单元相交
					double t_ratio[6];
					//线段所在直线与X轴垂直平面相交
					t_ratio[0] = (Hexc[0].x - point[0].x)/(point[1].x - point[0].x);
					t_ratio[1] = (Hexc[1].x - point[0].x)/(point[1].x - point[0].x);
					//线段所在直线与Y轴垂直平面相交
					t_ratio[2] = (Hexc[0].y - point[0].y)/(point[1].y - point[0].y);
					t_ratio[3] = (Hexc[1].y - point[0].y)/(point[1].y - point[0].y);
					//线段所在直线与Z轴垂直平面相交
					t_ratio[4] = (Hexc[0].z - point[0].z)/(point[1].z - point[0].z);
					t_ratio[5] = (Hexc[1].z - point[0].z)/(point[1].z - point[0].z);
		
					int count=0;
					for(int k=0; k<6; k++)
					{
						if(t_ratio[k]>0&&t_ratio[k]<1)
						{
							Point_3D ipoi = point[0]+(point[1]-point[0])*t_ratio[k];
							if((k==0||k==1)&&
								ipoi.y>=Hexc[0].y&&ipoi.y<=Hexc[1].y&&
								ipoi.z>=Hexc[0].z&&ipoi.z<=Hexc[1].z) count++;
							if((k==2||k==3)&&
								ipoi.x>=Hexc[0].x&&ipoi.x<=Hexc[1].x&&
								ipoi.z>=Hexc[0].z&&ipoi.z<=Hexc[1].z) count++;
							if((k==4||k==5)&&
								ipoi.x>=Hexc[0].x&&ipoi.x<=Hexc[1].x&&
								ipoi.y>=Hexc[0].y&&ipoi.y<=Hexc[1].y) count++;
						}
					}
					if(count==2) 	elements[elen].relative_cnts.push_back(i);
					else if(count>2) { hout <<"错误！直线与单元表面的交点个数大于2个，请检查！ ...Element_relative_cnts" << endl;  return 0; }
				}
	}

	return 1;
}
//---------------------------------------------------------------------------
//输出检测单元是否匹配了正确的相关纳米管线段信息
int Mesher::Test_element_relative_cnts(vector<Point_3D> &cnps)const
{
	ofstream otec("Test_element_cnts.dat");
	otec << "TITLE = Test_element_cnts" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	const int count = 10;
	//随机输出N个四面体单元
	srand((unsigned int)time(NULL));	//用于随机生成的起始时间
	for(int i=0; i<count; i++)
	{
		int ele_num = rand()%(int)elements.size();
		//---------------------------------------------------------------------------
		//输出单胞框架
		otec << "ZONE N=" << 8 << ", E=" << 1 << ", F=FEPOINT, ET=BRICK" << endl;
		for (int j=0; j<(int)elements[ele_num].nodes_id.size(); j++)
		{
			otec << nodes[elements[ele_num].nodes_id[j]].x << "  "
					<< nodes[elements[ele_num].nodes_id[j]].y << "  "
					<< nodes[elements[ele_num].nodes_id[j]].z << endl;
		}
	
		otec << "1 2 3 4 5 6 7 8" << endl;
		otec << endl << endl;
	
		//---------------------------------------------------------------------------
		//输出纳米管线三维构造
		vector<int> veci_temp;
		vector<vector<int> > veci;
		for(int j=0; j<(int)elements[ele_num].relative_cnts.size(); j++)
		{
			veci_temp.push_back(elements[ele_num].relative_cnts[j]);  //插入节点

			if(j==(int)elements[ele_num].relative_cnts.size()-1||
				elements[ele_num].relative_cnts[j]+1!=elements[ele_num].relative_cnts[j+1])
			{
				veci.push_back(veci_temp);
				veci_temp.clear();
			}
		}

		for(int j=0; j<(int)veci.size(); j++)
		{
			otec << "ZONE T=\"Line\"" << endl;
			otec << "i=1," << "j=" << (int)veci[j].size() << ", f=point" << endl;
			for (int k=0; k<(int)veci[j].size(); k++)
			{
				otec << cnps[veci[j][k]].x << "  " << cnps[veci[j][k]].y << "  " << cnps[veci[j][k]].z << endl;
			}
			otec << endl << endl;
		}
	}

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//输出Tecplot可视化网格数据
int Mesher::Export_mesh_data_tecplot(string output_file_name)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE = Tecplot_Meshes" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size() << ", F=FEPOINT, ET=BRICK" << endl;
	for (int i=0; i < (int)nodes.size(); i++)
	{
		otec << nodes[i].x << "  " << nodes[i].y << "  " << nodes[i].z << endl;
	}
	otec << endl;
	for (int i=0; i < (int)elements.size(); i++)
	{
		otec	 << elements[i].nodes_id[0]+1 << "  " << elements[i].nodes_id[1]+1 << "  " 
				 << elements[i].nodes_id[2]+1 << "  " << elements[i].nodes_id[3]+1 << "  " 
				 << elements[i].nodes_id[4]+1 << "  " << elements[i].nodes_id[5]+1 << "  " 
				 << elements[i].nodes_id[6]+1 << "  " << elements[i].nodes_id[7]+1 << endl;
	}
	otec.close();
	return 1;
}
//---------------------------------------------------------------------------
//根据给定的i j k以及最大i_max j_max k_max，决定给定节点的位置（角点、边界线点、边界面点、内点）
//返回值：0内点，1边界面点，2边界线点，3角点
//i--z坐标； j--y坐标；k--x坐标
int Mesher::Deter_node_type(const int i, const int j, const int k, const int i_max, const int j_max, const int k_max)const
{

	int type=0;
	if( i == 0 )	type++;
	else if( i == i_max ) type++;

	if( j == 0 )	type++;
	else if( j == j_max ) type++;

	if( k == 0 )	type++;
	else if( k == k_max ) type++;

	return type;
}
//---------------------------------------------------------------------------
//一维存储周期边界节点编号及起始点信息
int Mesher::Record_periodic_boundary_nodes_number(const struct RVE_Geo &cell_geo)
{
	//用于记录那些节点已经被记录过(0:没有记录过，1:记录过)
	vector<int> nflag(nodes.size(),0);

	//生成各个方向点的最大编号
	int k_max = int(cell_geo.len_x/dx);
	int j_max = int(cell_geo.wid_y/dy);
	int i_max = int(cell_geo.hei_z/dz);

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//周期边界节点编号及相关信息
	int count = 0;  //记录整体节点编号
	for( int i=0; i<=i_max; i++ )
		for( int j=0; j<=j_max; j++ )
			for( int k=0; k<=k_max; k++ )
			{
				if(nodes[count].type!=0&&nflag[count]==0)  //不是内点并且没有记录过
				{
					int key = 0;
					for(int m=i; m<=i_max; m+=i_max)
						for(int n=j; n<=j_max; n+=j_max)
							for(int p=k; p<=k_max; p+=k_max)
							{
								int num = m*(j_max+1)*(k_max+1) + n*(k_max+1) + p;
								peri_bnods[0].push_back(num);
								nflag[num]=1;	//此点已经记录过
								peri_bnods[1].push_back(key);  //标记位置
								key++;
							}
				}
				count++;  //整体节点编号加1
			}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//对一维存储按相关点从大到小重新排序
	int left = 0;
	int right = 0;
	for(int i=0; i<(int)peri_bnods[1].size(); i++)
	{
		if(i+1==(int)peri_bnods[1].size()||peri_bnods[1][i+1]==0) 
		{
			right = i;
			//---------------------------------------------------------------
			//将从小到大排序变为从大到小排序
			while(left<right)
			{
				int temp_peri = peri_bnods[0][right];
				peri_bnods[0][right] = peri_bnods[0][left];
				peri_bnods[0][left] = temp_peri;
				left++;
				right--;
			}
			//---------------------------------------------------------------
			left = i+1;
		}
	}

	//输出检查
	//hout << "peri_bnods:" << endl << endl;
	//for(int i=0; i<(int)peri_bnods[0].size(); i++)
	//{
	//	hout << peri_bnods[0][i] << "  " << peri_bnods[1][i] << endl;
	//}

	return 1;
}
//---------------------------------------------------------------------------
//生成格子用于预估局部模型等效模量
int Mesher::Generate_grids_for_effective_stiffness(ifstream &infile, const double &decayR, double &grid_vol)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//读入格子尺寸(长、宽和高)
	istringstream istr(Get_Line(infile));	
	istr >> dx >> dy >> dz;
	if(dx<0||dy<0||dz<0)
	{
		hout << "预估刚度矩阵时，生成格子长宽高的分割细度输入错误！ 请检查！ ...Generate_grids_for_effective_stiffness" << endl;
		return 0;
	}

	//格子的体积
	grid_vol = dx*dy*dz;

	//生成各个方向的点数
	int x_sec_num = 2*(int(decayR/dx-Zero)+2);  //加一个小值为了避免取整的错误
	int y_sec_num = 2*(int(decayR/dy-Zero)+2);  //+2：一个+1是点的个数的一半比段数多1；一个+1是保证decayR<=n*dx(y,z)时，网格数要多1
	int z_sec_num = 2*(int(decayR/dz-Zero)+2);

	Point_3D origin(0,0,0);
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//生成节点
	for( int i=0; i<z_sec_num; i++ )
	{
		double z = origin.z + i * dz ;
		for( int j=0; j<y_sec_num; j++ )
		{
			double y = origin.y + j * dy ;
			for( int k=0; k<x_sec_num; k++ )
			{
				double x = origin.x + k * dx ;

				Node nd(x,y,z);
				nodes.push_back(nd);
			}
		}
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//生成六面体单元

	//定义六面体向量
	for( int i=0; i<z_sec_num-1; i++ )
		for( int j=0; j<y_sec_num-1; j++ )
			for( int k=0; k<x_sec_num-1; k++ )
			{
				//六面体的八个顶点
				Element ele_temp;
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + j * x_sec_num + k);
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + j * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k);
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k );
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1);
				ele_temp.nodes_id.push_back(( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k);

				elements.push_back(ele_temp);
			}

	return 1;
}
//---------------------------------------------------------------------------
//读入一行信息，并跳过注释行（以"%"开头）；
string Mesher::Get_Line(ifstream &infile)const
{
	string s;
	//读入信息一行
	getline(infile,s);
	//跳过注释行     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===========================================================================
