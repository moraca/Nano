//===========================================================================
// MatBase.cpp
// 材料库类的成员函数
// Member functions in a class of material database
//===========================================================================
#include "MatBase.h"

//---------------------------------------------------------------------------
//生成材料数据库
int MatBase::Generate_matbase(ifstream &infile)
{
	//---------------------------------------------------------------------------
	//读取材料的刚度参数
	MatPro Smat;
	//读取材料刚度信息
	int mat_type;
	istringstream imat0(Get_Line(infile));
	imat0 >> mat_type;

	if(mat_type==0)//各向同性材料
	{
		double E, Nu;
		imat0 >> E >> Nu;
		Smat.set_ela_para(E, Nu);
	}
	else if(mat_type==1)//横观各向同性
	{
		double E1, E2, Nu1, Nu2, G2;
		imat0 >> E1 >> E2 >> Nu1 >> Nu2 >> G2;
		Smat.set_ela_para(E1, E2, Nu1, Nu2, G2);
	}
	else if(mat_type==2)//正交各向异性
	{
		double E1, E2, E3, Nu12, Nu23, Nu13, G12, G23, G13;
		imat0 >> E1 >> E2 >> E3 >> Nu12 >> Nu23 >> Nu13 >> G12 >> G23 >> G13;
		Smat.set_ela_para(E1, E2, E3, Nu12, Nu23, Nu13, G12, G23, G13);
	}
	else
	{
		hout << "读取基本材料类型" << mat_type << "不在0和2之间。" << endl;
		return 0;
	}
	if(Smat.Generate_local_elas_matrix()==0) return 0; //生成材料的刚度矩阵
	mats_vec.push_back(Smat);

	//---------------------------------------------------------------------------
	//读取并计算长程力的等效刚度参数
	MatPro known_mat;
	//读取材料刚度信息
	istringstream imat(Get_Line(infile));
	imat >> mat_type;

	if(mat_type==0)//各向同性材料
	{
		double E;
		imat >> E;
		known_mat.set_ela_para(E);
	}
	else if(mat_type==1)//横观各向同性(输入泊松比)
	{
		double E1,E2,Nu2;
		imat >> E1 >> E2 >> Nu2;
		known_mat.set_ela_para(E1, E2, Nu2);
	}
	else if(mat_type==3)//横观各向同性(输入剪切模量)
	{
		double E1,E2,G2;
		imat >> E1 >> E2 >> G2;
		known_mat.set_ela_para_transverse(E1, E2, G2);
	}
	else if(mat_type==2)//正交各向异性
	{
		double E1,E2,E3,Nu12,Nu23,Nu13;
		imat >> E1 >> E2 >> E3 >> Nu12 >> Nu23 >> Nu13;
		known_mat.set_ela_para(E1, E2, E3, Nu12, Nu23, Nu13);
	}
	else
	{
		hout << "读取长程力等效刚度材料类型" << mat_type << "不在0和2之间。" << endl;
		return 0;
	}

	if(known_mat.Generate_nonlocal_elas_matrix()==0) return 0; //生成材料的刚度矩阵

	//---------------------------------------------------------------------------
	//读取衰减函数参数
	istringstream in2(Get_Line(infile));
	in2 >> decay.R >> decay.radius;
	if(decay.R<0||decay.radius<=0) { hout << "注意！衰减函数参数小于0或者小于等于0，请重新输入！" << endl; return 0; }

	//修改材料系数，根据非局部长程力
	if(decay.R>0.0) 
	{
		MatPro mat;
		//---------------------------------------------------------------------------
		//生成格子用于预估局部模型等效刚度阵
		Mesher grid;
		double grid_vol;
		if(grid.Generate_grids_for_effective_stiffness(infile, decay.R, grid_vol)==0) return 0;

		//取高斯点及加权系数
		Gauss gau;
		if(gau.Generate_gauss(infile)==0) return 0;
		
		//预估局部模型等效刚度阵
		MathMatrix Km(6, 6);
		for(int i=0; i<6; i++)
		{
			double a[6] = {0};
			a[i] = 1.0;
			if(Estimate_stiffness_long_range_interaction(grid.elements, grid.nodes, grid_vol, gau.gauss, gau.weight, a, mat, 0)==0) return 0;
			Km.element[0][i] = mat.elas_matrix[0][0];
			Km.element[1][i] = mat.elas_matrix[1][1];
			Km.element[2][i] = mat.elas_matrix[2][2];
			Km.element[3][i] = mat.elas_matrix[0][1];
			Km.element[4][i] = mat.elas_matrix[1][2];
			Km.element[5][i] = mat.elas_matrix[0][2];
		}

		//重新计算
		double Cm[6] = {	known_mat.elas_matrix[0][0],	known_mat.elas_matrix[1][1],	known_mat.elas_matrix[2][2], 
										known_mat.elas_matrix[0][1], known_mat.elas_matrix[1][2], known_mat.elas_matrix[0][2] };

		MathMatrix Dm(6,6);
		Dm = (Km.Transpose()*Km).Inverse()*Km.Transpose();
		for(int i=0; i<6; i++)
		{
			decay.acoe[i] = 0.0;
			for(int j=0; j<6; j++) decay.acoe[i] += Dm.element[i][j]*Cm[j];
		}

		if(Estimate_stiffness_long_range_interaction(grid.elements, grid.nodes, grid_vol, gau.gauss, gau.weight, decay.acoe, mat, 0)==0) return 0; //最后一位形参如果输入1，表示输出error矩阵

		//---------------------------------------------------------------------------
		//根据理论公式反求系数，并验证一致性（由于要除以Horizon半径，所以Horizon半径应该大于1，避免误差过大）
//		mat.Compare_coef_by_analysis_formula(mat_type, decay);
		//---------------------------------------------------------------------------
		//根据理论公式反求刚度矩阵，并验证一致性
//		mat.Compare_matrix_by_analysis_formula(decay);

		//---------------------------------------------------------------------------
		//根据弹性矩阵，反求材料参数
		mat.Get_ele_para_by_ela_matrix();

		mats_vec.push_back(mat); //插入材料向量
	}
	else  istringstream istr_precision(Get_Line(infile));   //跳过读入

	//---------------------------------------------------------------------------
	//输出所有材料的刚度矩阵及对应的杨氏模量、泊松比
	Print_stiffness_matrix("Material_stiffness_matrix.dat");
	//根据弹性矩阵，反求材料参数
	known_mat.Get_ele_para_by_ela_matrix();
	Print_stiffness_matrix("Material_stiffness_matrix.dat", known_mat, 1);

	return 1;
}
//------------------------------------------------------------------------------
//输出所有材料的刚度矩阵到文件中
void MatBase::Print_stiffness_matrix(string print_name)const
{
	ofstream opri(print_name.c_str());
	//材料刚度矩阵输出
	for(int num=0; num<(int)mats_vec.size(); num++)
	{
		opri << "Material " << num+1 << "  stiffness_matrix:" << endl;
		for(int i=0; i<6; i++)
		{
			for(int j=0; j<6; j++) opri << setprecision(18) << setw(24) << setiosflags(ios::right)  << mats_vec[num].elas_matrix[i][j] << "  ";
			opri << endl;
		}
		opri << endl;

		//材料的杨氏模量、泊松比和剪切模量输出
		opri << "%Young's Modulus：E11，E22，E33 ; Poisson's Ratio：Nu12，Nu23，Nu13 ; Shear Modulus：G12，G23，G13 " << endl;
		opri << mats_vec[num].E11 << " " << mats_vec[num].E22 << " " << mats_vec[num].E33 << "   ";
		opri << mats_vec[num].Nu12 << " " << mats_vec[num].Nu23 << " " << mats_vec[num].Nu13 << "   ";
		opri << mats_vec[num].G12 << " " << mats_vec[num].G23 << " " << mats_vec[num].G13 << "   ";

		opri << endl << endl;
	}

	opri.close();
}
//------------------------------------------------------------------------------
//输出材料的刚度矩阵到文件中
void MatBase::Print_stiffness_matrix(string print_name, MatPro &mat, int key)const
{
	ofstream opri;
	if(key==0) opri.open(print_name.c_str());
	else opri.open(print_name.c_str(), ios_base::app);
	//材料刚度矩阵输出
	opri << "Compared Material  stiffness_matrix:" << endl;
	for(int i=0; i<6; i++)
	{
		for(int j=0; j<6; j++) opri << setprecision(18) << setw(24) << setiosflags(ios::right)  << mat.elas_matrix[i][j] << "  ";
		opri << endl;
	}
	opri << endl;

	//材料的杨氏模量、泊松比和剪切模量输出
	opri << "%Young's Modulus：E11，E22，E33 ; Poisson's Ratio：Nu12，Nu23，Nu13 ; Shear Modulus：G12，G23，G13 " << endl;
	opri << mat.E11 << " " << mat.E22 << " " << mat.E33 << "   ";
	opri << mat.Nu12 << " " << mat.Nu23 << " " << mat.Nu13 << "   ";
	opri << mat.G12 << " " << mat.G23 << " " << mat.G13 << "   ";

		opri << endl << endl;

	opri.close();
}
//------------------------------------------------------------------------------
//预估局部模型等效刚度阵
int MatBase::Estimate_stiffness_long_range_interaction(const vector<Element> &elements, const vector<Node> &nodes, const double &grid_vol,
						const vector<Node> &gauss, const vector<double> &weight, const double a[], MatPro &mat, const int &output_key)const
{
	//--------------------------------------------------	
	//计算主单元编号及高斯点坐标
	if((int)elements.size()==0) { hout << "预估刚度阵所用网格单元数为0, 请检查！" << endl; return 0; }
	int nume = ((int)elements.size()-1)/2;

	//--------------------------------------------
	//高斯点个数
	const int gau_num = (int)weight.size();

	//--------------------------------------------
	//计算主单元高斯点坐标
	Point_3D gd_temp(0,0,0);
	vector<Point_3D> gd_main(gau_num, gd_temp);

	for(int count=0; count<gau_num; count++)
	{
		double Nshape[8] = {0};
		
		//计算高斯积分点整体坐标
		Nshape[0]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
		Nshape[1]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
		Nshape[2]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
		Nshape[3]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
		Nshape[4]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
		Nshape[5]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
		Nshape[6]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);
		Nshape[7]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);

		//计算高斯积分点坐标
		for(int j=0; j<8; j++) 
		{
			gd_main[count].x += Nshape[j]*nodes[elements[nume].nodes_id[j]].x;
			gd_main[count].y += Nshape[j]*nodes[elements[nume].nodes_id[j]].y;
			gd_main[count].z += Nshape[j]*nodes[elements[nume].nodes_id[j]].z;
		}
	}

	//--------------------------------------------------	
	//长程力等效刚度阵
	double (*Long_range_stiffness)[6][6] = new double [gau_num][6][6];
	//长程力等效刚度阵初始化为零
	for(int count =0; count<gau_num; count++)		
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				Long_range_stiffness[count][i][j] = 0.0;

	//--------------------------------------------------	
	//循环所有单元
	for(int i=0; i<(int)elements.size(); i++)
	{
		//循环主单元高斯点计算积分
		for(int count1=0; count1<gau_num; count1++)
		{
			//循环外单元高斯点计算积分
			for(int count2=0; count2<gau_num; count2++)
			{
				double Nshape[8] = {0};
				//--------------------------------------------
				//计算高斯积分点整体坐标
				Nshape[0]=0.125*(1.0-gauss[count2].x)*(1.0-gauss[count2].y)*(1.0-gauss[count2].z);
				Nshape[1]=0.125*(1.0+gauss[count2].x)*(1.0-gauss[count2].y)*(1.0-gauss[count2].z);
				Nshape[2]=0.125*(1.0+gauss[count2].x)*(1.0+gauss[count2].y)*(1.0-gauss[count2].z);
				Nshape[3]=0.125*(1.0-gauss[count2].x)*(1.0+gauss[count2].y)*(1.0-gauss[count2].z);
				Nshape[4]=0.125*(1.0-gauss[count2].x)*(1.0-gauss[count2].y)*(1.0+gauss[count2].z);
				Nshape[5]=0.125*(1.0+gauss[count2].x)*(1.0-gauss[count2].y)*(1.0+gauss[count2].z);
				Nshape[6]=0.125*(1.0+gauss[count2].x)*(1.0+gauss[count2].y)*(1.0+gauss[count2].z);
				Nshape[7]=0.125*(1.0-gauss[count2].x)*(1.0+gauss[count2].y)*(1.0+gauss[count2].z);

				//--------------------------------------------
				//计算高斯积分点坐标
				Point_3D gd_assist(0,0,0);
				for(int j=0; j<8; j++) 
				{
					gd_assist.x += Nshape[j]*nodes[elements[i].nodes_id[j]].x;
					gd_assist.y += Nshape[j]*nodes[elements[i].nodes_id[j]].y;
					gd_assist.z += Nshape[j]*nodes[elements[i].nodes_id[j]].z;
				}				
			
				//--------------------------------------------
				//计算长程作用衰减函数值
				const double x = gd_assist.x-gd_main[count1].x;
				const double y = gd_assist.y-gd_main[count1].y;
				const double z = gd_assist.z-gd_main[count1].z;	

				const double dis_squr = x*x + y*y + z*z;
				const double poi_dis = sqrt(dis_squr);		//两点之间的距离
			
				if(poi_dis>decay.R+Zero||poi_dis<Zero) continue;		//圆形积分域
			
				const double cos2sita = z*z/dis_squr;		//cos(sita)^2
				double cos2pha;
				if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
				else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

				double sum = a[0] + a[1]*0.5*(3*cos2sita-1) + a[2]*(2*cos2pha-1)*3*(1-cos2sita) + a[3]*0.125*(35*cos2sita*cos2sita - 30*cos2sita + 3)
										+ a[4]*(2*cos2pha-1)*7.5*(7*cos2sita-1)*(1-cos2sita) + a[5]*(8*cos2pha*cos2pha-8*cos2pha+1)*105*(1-cos2sita)*(1-cos2sita);

				const double gv = exp(-poi_dis/decay.radius)*sum*weight[count2]; //衰减以及权重函数值
//				const double gv = sum*weight[count2]/dis_squr; //衰减以及权重函数值  (用于对比分析公式和数值公式)

				const double r_comp[6] = {x*x, y*y, z*z, x*y, y*z, z*x}; //等效长程力衰减计算
				const double temp_gmat[6] = {gv*r_comp[0], gv*r_comp[1], gv*r_comp[2], gv*r_comp[3], gv*r_comp[4], gv*r_comp[5]}; //衰减函数矩阵的对称项

				//--------------------------------------------
				//计算长程力等效刚度矩阵
				for(int j=0; j<6; j++)
					for(int k=0; k<6; k++)
						Long_range_stiffness[count1][j][k] += temp_gmat[j]*r_comp[k]; //相当于gv*r_comp[j]*r_comp[k]
			}
		}
	}

	//计算雅可比行列式的值
	//--------------------------------------------
	//形函数N对中心点的偏导矩阵
	double diff[3][8];
	diff[0][0]=-0.125;
	diff[0][1]=-diff[0][0];                         
	diff[0][2]=0.125;
	diff[0][3]=-diff[0][2];
	diff[0][4]=-0.125;
	diff[0][5]=-diff[0][4];
	diff[0][6]=0.125;
	diff[0][7]=-diff[0][6];

	diff[1][0]=-0.125;
	diff[1][1]=-0.125;
	diff[1][2]=-diff[1][1];
	diff[1][3]=-diff[1][0];
	diff[1][4]=-0.125;
	diff[1][5]=-0.125;
	diff[1][6]=-diff[1][5];
	diff[1][7]=-diff[1][4];

	diff[2][0]=-0.125;
	diff[2][1]=-0.125;
	diff[2][2]=-0.125;
	diff[2][3]=-0.125;
	diff[2][4]=-diff[2][0];
	diff[2][5]=-diff[2][1];
	diff[2][6]=-diff[2][2];
	diff[2][7]=-diff[2][3];
	//--------------------------------------------------
	//单元节点坐标矩阵
	double elenode[8][3];
	for(int j=0; j<8; j++)
	{
		elenode[j][0]=nodes[elements[nume].nodes_id[j]].x;
		elenode[j][1]=nodes[elements[nume].nodes_id[j]].y;
		elenode[j][2]=nodes[elements[nume].nodes_id[j]].z;
	}
	//--------------------------------------------------
	//J矩阵
	double Jmatrix[3][3];
	//以上两个矩阵的积
	for(int j=0; j<3; j++)
		for(int k=0; k<3; k++)
		{
			Jmatrix[j][k]=0;
			for(int m=0; m<8; m++)
			Jmatrix[j][k] += diff[j][m]*elenode[m][k];
		}
	//--------------------------------------------------
	//求出J矩阵的行列式
	double J_val = Jmatrix[0][0]*(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])
								-Jmatrix[0][1]*(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])
								+Jmatrix[0][2]*(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0]);

	//--------------------------------------------------
	//误差的最大矩阵和最小矩阵
	double error_max[6][6], error_min[6][6];	

	//--------------------------------------------------
	//循环主单元高斯点计算各个高斯点上的等效长程力矩阵
	for(int i=0; i<6; i++)
		for(int j=0; j<6; j++)
		{
			mat.elas_matrix[i][j] = 0;
			for(int count=0; count<gau_num; count++)
			{
				if(count==0) { error_max[i][j] = Long_range_stiffness[count][i][j];  error_min[i][j] = Long_range_stiffness[count][i][j]; }
				else 
				{
					if(error_max[i][j] < Long_range_stiffness[count][i][j]) error_max[i][j] = Long_range_stiffness[count][i][j];
					if(error_min[i][j] > Long_range_stiffness[count][i][j]) error_min[i][j] = Long_range_stiffness[count][i][j];
				}
				
				mat.elas_matrix[i][j] += Long_range_stiffness[count][i][j]*weight[count];
			}
			mat.elas_matrix[i][j] = mat.elas_matrix[i][j]*0.5*J_val*J_val/grid_vol;
		}

	if(output_key==1)
	{
		//用于输出检测
		//for(int count=0; count<gau_num; count++)
		//{
		//	hout << endl << "The " << count << " effective matrix:" << endl;
		//	for(int i=0; i<6; i++)
		//	{
		//		for(int j=0; j<6; j++)
		//			hout << Long_range_stiffness[count][i][j] << "  ";
		//		hout << endl;
		//	}
		//}

		//输出误差用于检测
		hout << endl << "The final error effective matrix:" << endl;
		for(int i=0; i<6; i++)
		{
			for(int j=0; j<6; j++)
				hout << (error_max[i][j]-error_min[i][j])*(0.5*J_val*J_val/grid_vol)/mat.elas_matrix[i][j] << "  ";
			hout << endl;
		}
	}

	return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//读入一行信息，并跳过注释行（以"%"开头）；
string MatBase::Get_Line(ifstream &infile)const
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
