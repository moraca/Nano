//===========================================================================
// Algorithm_FEM.cpp
// 算法及有限元实现模块类成员函数
// Member functions in classes to implement algorithm and finite element method
//===========================================================================
#include "Algorithm_FEM.h"

//--------------------------------------------------------------------------
//算法的有限元求解过程
int Algorithm_FEM::Solve(ifstream &infile, vector<Node> &nodes, const vector<int>* peri_bnods, vector<Element> &elements,  const vector<MatPro> &mats, const struct Decay_Para &decay, const struct RVE_Geo &cell, 
											   const struct CNT_Geo &cnts, const vector<Point_3D> &cnps)
{
	clock_t ct0,ct1;  //定义变量记录执行的开始和结束时间

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//读入计算类型数据
	if(Import_computational_data(infile)==0) return 0;
	backup_file_name = "Displacement_Solution_Backup";  //备份位移解及等效应变能数据文件名
	if(wr_mod=="read_all")  return 1;   //读入全部数据包括位移解和等效应变能
	
	Glosmat = new Global_Stiff_Matrix;		//定义总体刚度矩阵类，每个点三个自由度，由于wr_mod读写模式的缘由要提前在此定义

	ct0 = clock();
	hout << "-_- 开始统计单元的相关单元和节点的相关节点" << endl;
	if(Deter_relative_nodes_elements(nodes, elements, decay.R, cell, com_mod)==0) return 0; //决定单元的相关单元和节点的相关节点, 单元的相关单元信息要在估计等效应变能时用到

	ct1 = clock();
	hout << "    统计单元的相关单元和节点的相关节点耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
	hout << "^_^ 统计单元的相关单元和节点的相关节点完毕!" << endl << endl;

	//---------------------------------------------------------------------------
	//定义存储向量, 用于等效应变能的计算
	vector<long> backup_iz;
	vector<int> backup_ig;
	vector<double> backup_matrix;
	vector<double> backup_equright[9];
	double backup_ege[9] = { 0 };

	if(wr_mod=="write")
	{
		//-----------------------------------------------------------------------------------------------------------------------------------------
		//生成紧缩存储节点关系向量
		ct0 = clock();
		hout << "-_- 开始紧缩存储刚度矩阵" << endl;

		Solv = new SolveEqu;	//定义解方程组类
		vector<long> Iz;		//Iz,Ig是变带宽存储的相关信息
		vector<int> Ig;			//Iz,Ig是变带宽存储的相关信息
		if(Solv->izig(nodes, Iz, Ig)==0) return 0;
		//注释：节点的相关节点信息和单元向量的相关单元信息要在后面用到, 在此不清零

		ct1 = clock();
		hout << "    紧缩存储刚度矩阵耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
		hout << "^_^ 紧缩存储刚度矩阵完毕!" << endl << endl;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//生成总体刚度矩阵
		ct0 = clock();
		hout << "-_- 开始生成总纲阵" << endl;

		//申请总刚度矩阵存储空间(单刚叠加总刚的方式量身申请尺度)
		const long Total_Matrix_Size = 6*(long)nodes.size() + 9*Iz.back();
		vector<double> total_matrix(Total_Matrix_Size, 0);
		if(Total_Matrix_Size>(long)total_matrix.max_size()) { hout << "注意！申请的刚度矩阵尺寸大于向量vector所能容许的最大值！" << endl; return 0; }
		hout << "total_matrix.size: " << 	total_matrix.size() << "  " << "Iz.back: " << Iz.back() << endl;

		vector<double> nl_equright[9];		//用于记录长程非局部力在边界层处由于周期边界条件所得的右端项
		for(int i=0; i<9; i++) nl_equright[i].assign(3*nodes.size(), 0.0);	//初始化

		if(Glosmat->Gen_global_stiff_matrix(infile, cnps, nodes, elements, mats, decay, cell, cnts, com_mod, Iz, Ig, total_matrix, nl_equright, backup_ege)==0) return 0; //生成总体刚度矩阵
		//后面在计算等效应变能的时候仍要用到这个类，所以这里没有删除

		//------------------------------------------------------------------------
		//数据读入或者写入二进制文件，以备以后单独求解线性方程组时用到
		string maeq_mod = "write";
//		string maeq_mod = "read";
		if(maeq_mod=="read") { Get_Line(infile); Get_Line(infile); }
		if(Write_read_matrix_equright_data(Total_Matrix_Size, total_matrix, nl_equright, backup_ege, maeq_mod)==0) return 0;

		//---------------------------------------------------------------------------
		backup_iz = Iz;									//记录原始Iz
		backup_ig = Ig;									//记录原始Ig
		backup_matrix = total_matrix;			//记录原始刚度矩阵
		for(int i=0; i<9; i++) backup_equright[i] = nl_equright[i];  //记录原始右端向量
		//---------------------------------------------------------------------------

		ct1 = clock();
		hout << "    生成总纲阵耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
		hout << "^_^ 生成总纲阵完毕!" << endl << endl;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//处理位移约束
		ct0 = clock();
		hout << "-_- 开始处理位移约束" << endl;

		int bnod_num = 0;		//用于记录边界节点的个数，注意要初始化零值
		vector<int> ip;			//用于记录边界节点类型信息
		vector<double> vp;	//用于记录边界节点上的值
		if(Solv->Fixed_displacement_constraints(infile, nodes, bnod_num, ip, vp)==0) return 0;		//设置固定位移约束
		ct1 = clock();
		hout << "    位移约束处理耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
		hout << "^_^ 位移约束处理完毕!" << endl << endl;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//处理周期边界条件
		ct0 = clock();
		hout << "-_- 开始处理周期边界条件" << endl;

	//	Solv->Complete_matrix_equright_testing_periodical_constraints(peri_bnods, total_matrix, nl_equright, nodes, Iz, Ig);   //用完全的刚度矩阵及右端项检测周期边界条件的处理
		Solv->Periodical_boundary_constraints(nodes, peri_bnods, Iz, Ig, total_matrix, nl_equright);
		for(int i=0; i<(int)nodes.size(); i++) nodes[i].relative_nods.clear();  //节点的相关节点信息在此清零
	//	Solv->Export_complete_matrix_equright(total_matrix, nl_equright, nodes, Iz, Ig);		//输出完全的刚度矩阵及右端项用于检测
		ct1 = clock();
		hout << "    周期边界条件处理耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
		hout << "^_^ 周期边界条件处理完毕!" << endl << endl;

		//---------------------------------------------------------------------------
		//处理三维固定位移约束点的零位移值条件
		ct0 = clock();
		hout << "-_- 开始三维固定位移约束点的零位移值条件" << endl;

		Solv->Deal_with_displacement_zero_value(bnod_num, (int)nodes.size(), Iz, Ig, ip, vp, nl_equright, total_matrix);
		ct1 = clock();
		hout << "    三维固定位移约束点的零位移值条件处理耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
		hout << "^_^ 三维固定位移约束点的零位移值条件处理完毕!" << endl << endl;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//循环并行求解线性方程组
		ct0 = clock();
		hout << "-_- 开始循环并行求解线性方程组" << endl;

		vector<double> temp_solu;			//临时定义
		U_Solution.assign(9, temp_solu);	//初始化
		//执行openmp
		#pragma omp parallel
		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<9; i++)
		{
			clock_t ct2, ct3;	//计时器
			ct2 = clock();
			vector<double> solution(3*nodes.size(), 0.0); //定义
			Solv->Solve_linear_equations(bnod_num, (int)nodes.size(), Iz, Ig, ip, vp, total_matrix, nl_equright[i], solution);	

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//用于输出周期边界上节点的结果用于比对
//			Solv->Compared_periodical_bounday_solutions(i, peri_bnods, solution);
			//---------------------------------------------------------------------------
			//记录这组解
			#pragma omp critical
			{
				U_Solution[i] = solution;
			}
			ct3 = clock();
			hout << "    ^_^ 第" << i+1 <<"组求解线性方程组完毕! 共耗时：" << (double)(ct3 - ct2)/CLOCKS_PER_SEC << "秒"  << endl;
		}

		//-----------------------------------------------------------------------------------------------------------------------------------------		
		//根据周期边界条件附加其他节点的解
		Solv->Add_periodical_solution(peri_bnods, nodes, U_Solution);

		//-----------------------------------------------------------------------------------------------------------------------------------------		
		//输出位移解用于检测
		//for(int i=0; i<9; i++)
		//{
		//	hout << endl << "//==============================================================" << endl;	
		//	hout << "    ^_^ 第" << i+1 <<"组位移解："<< endl;		
		//	for(int j=0; j<(int)nodes.size(); j++)
		//	{
		//		hout << setwp(6) << j << "  ";
		//		hout << setwp(15,6) << U_Solution[i][3*j] << "  " << setwp(15,6) << U_Solution[i][3*j+1] << "  " << setwp(15,6) << U_Solution[i][3*j+2]<< endl;
		//	}
		//}
		ct1 = clock();
		hout << "    循环并行求解线性方程组处理耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
		hout << "^_^ 循环并行求解线性方程组处理完毕!" << endl << endl;
		delete Solv;		//删除解方程组类
	}
	else if(wr_mod!="read_solution")	{ hout << "注意！读写指令即不是读read_all或read_solution也不是写指令，请检查！" << endl;	return 0; }

	//-----------------------------------------------------------------------------------------------------------------------------------------		
	//读取位移解
	if(wr_mod=="read_solution") read_u_solution(backup_file_name);

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//计算等效应变能, 注意本函数中还需要GaussDats相关数据
	if(Glosmat->Calculate_equivalent_energy_simple(backup_iz, backup_ig, backup_matrix, backup_equright, backup_ege, U_Solution, equivalent_energy)==0) return 0; //生成总体刚度矩阵
//	if(Glosmat->Calculate_equivalent_energy(elements, mats, decay, U_Solution, cell, com_mod, equivalent_energy)==0) return 0; //生成总体刚度矩阵
	delete Glosmat;		//删除生成总刚度阵类

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//输出位移解及等效应变能数据用于记录
	write_u_solution_equivalent_energy(backup_file_name);
	
	return 1;
}
//---------------------------------------------------------------------------
//读入计算类型数据
int Algorithm_FEM::Import_computational_data(ifstream &infile)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//读入计算模型类型参数
	istringstream istr_com_mod(Get_Line(infile));
	istr_com_mod >> com_mod; 
	if(com_mod!="fem"&&com_mod!="nonlocal"&&com_mod!="hybrid")
	{ 
		hout <<"注意！计算模型选择参数输入错误，请重新输入计算模型参数！" << endl; 
		return 0; 
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//位移解的读写模式参数
	istringstream istr_wr_mod(Get_Line(infile));
	istr_wr_mod >> wr_mod;
	if(wr_mod=="read_all")	hout <<"注意！读入位移解及等效应变能备份数据！" << endl; 
	else if(wr_mod=="read_solution")	hout <<"注意！读入位移解备份数据！" << endl; 
	else if(wr_mod!="write")
	{
		hout <<"注意！位移解的读写模式参数输入错误，请重新输入读写模式参数！" << endl; 
		return 0;
	}

	return 1;
}
//---------------------------------------------------------------------------
//判断节点的相关节点和单元的相关单元
int Algorithm_FEM::Deter_relative_nodes_elements(vector<Node> &nodes, vector<Element> &elements, const double &dist, const struct RVE_Geo &cell, string mod)
{
	const int ES = (int)elements.size();
	const int NS = (int)nodes.size();
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//清零(避免出错)
	for(int i=0; i<ES; i++) elements[i].relative_eles.clear();
	for(int i=0; i<NS; i++) nodes[i].relative_nods.clear();

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//单元的相关单元
	if(mod=="fem")
	{
		for(int i=0; i<ES; i++) 
			elements[i].relative_eles.push_back(i);
	}
	else if(mod=="hybrid"||mod=="nonlocal") 
	{
		//生成单元的形心点向量
		vector<Point_3D> centrele(ES);
		for(int i=0; i<ES; i++)
		{
			Point_3D poi_tem(0,0,0);
			for(int j=0; j<8; j++)  //六面体单元节点个数8
			{
				poi_tem.x += nodes[elements[i].nodes_id[j]].x;
				poi_tem.y += nodes[elements[i].nodes_id[j]].y;
				poi_tem.z += nodes[elements[i].nodes_id[j]].z;
			}
			centrele[i] = poi_tem/8.0;
		}

		//保证了dist <= n*cell.delt_x(y,z)的情况能对单元中每个高斯点都取到完全的dist邻域
		double distx, disty, distz;
		distx = (int(dist/cell.delt_x-Zero)+1)*cell.delt_x;
		disty = (int(dist/cell.delt_y-Zero)+1)*cell.delt_y;
		distz = (int(dist/cell.delt_z-Zero)+1)*cell.delt_z;

		for(int i=0; i<ES; i++)
			for(int j=0; j<ES; j++)
			{
				if((fabs(centrele[j].x-centrele[i].x)<distx+Zero||fabs(centrele[j].x-centrele[i].x)>cell.len_x-distx-Zero)&&
					(fabs(centrele[j].y-centrele[i].y)<disty+Zero||fabs(centrele[j].y-centrele[i].y)>cell.wid_y-disty-Zero)&&
					(fabs(centrele[j].z-centrele[i].z)<distz+Zero||fabs(centrele[j].z-centrele[i].z)>cell.hei_z-distz-Zero))
				elements[i].relative_eles.push_back(j);
			}	
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//节点的相关节点(按从小到大排序)
	//经过反复验证，按如下方式写OpenMP代码，执行最快
	#pragma omp parallel
	{	
		vector<int> temp_num;
		vector<vector<int> > temp_rnod;

		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++)
			for(int j=0; j<8; j++)
			{
				vector<int> rela_nods;
				const int nodn = elements[i].nodes_id[j];
				for(int k=0; k<(int)elements[i].relative_eles.size(); k++)
				{
					const int eiek = elements[i].relative_eles[k];
					for(int m=0; m<8; m++)
					{
						const int renod = elements[eiek].nodes_id[m];
						if(nodn!=renod)
						{
							//二分法插入
							int left = 0;
							int right = (int)rela_nods.size()-1;
							while(right>=left)
							{
								int middle = (left + right)/2;
								if(rela_nods[middle] == renod) goto Node_Num_Same_0; //节点编号相同的情况
								else if(rela_nods[middle] > renod) right = middle - 1;
								else left = middle + 1;
							}

							//用这种方式比用vector.insert函数稍快一些
							rela_nods.push_back(0);
							for(int n=(int)rela_nods.size()-1; n>left; n--) rela_nods[n] = rela_nods[n-1];
							rela_nods[left] = renod;

Node_Num_Same_0: ;					
						}
					}
				}
				temp_num.push_back(nodn);   
				temp_rnod.push_back(rela_nods);
			}

		#pragma omp critical
		{
			for(int i=0; i<(int)temp_num.size(); i++)
			{
				if(nodes[temp_num[i]].relative_nods.size()==0)
				{
					nodes[temp_num[i]].relative_nods = temp_rnod[i];
				}
				else
				{
					int left = 0;
					for(int j=0; j<(int)temp_rnod[i].size(); j++)
					{
						//二分法插入
						int right = (int)nodes[temp_num[i]].relative_nods.size()-1;
						while(right>=left)
						{
							int middle = (left + right)/2;
							if(nodes[temp_num[i]].relative_nods[middle] == temp_rnod[i][j]) goto Node_Num_Same_1; //节点编号相同的情况
							else if(nodes[temp_num[i]].relative_nods[middle] > temp_rnod[i][j]) right = middle - 1;
							else left = middle + 1;
						}

						//用这种方式比用vector.insert函数稍快一些
						nodes[temp_num[i]].relative_nods.push_back(0);
						for(int n=(int)nodes[temp_num[i]].relative_nods.size()-1; n>left; n--) nodes[temp_num[i]].relative_nods[n] = nodes[temp_num[i]].relative_nods[n-1];
						nodes[temp_num[i]].relative_nods[left] = temp_rnod[i][j];

Node_Num_Same_1: ;					
					}
				}
			}
		}
	}

	//以下程序是没有执行OpenMP的程序
//	for(int i=0; i<ES; i++)
//		for(int j=0; j<8; j++)
//		{
//			const int nodn = elements[i].nodes_id[j];
//			for(int k=0; k<(int)elements[i].relative_eles.size(); k++)
//			{
//				const int eiek = elements[i].relative_eles[k];
//				for(int m=0; m<8; m++)
//				{
//					const int renod = elements[eiek].nodes_id[m];
//					if(nodn!=renod)
//					{
//						//二分法插入
//						int left = 0;
//						int right = (int)nodes[nodn].relative_nods.size()-1;
//						while(right>=left)
//						{
//							int middle = (left + right)/2;
//							if(nodes[nodn].relative_nods[middle] == renod) goto Node_Num_Same; //节点编号相同的情况
//							else if(nodes[nodn].relative_nods[middle] > renod) right = middle - 1;
//							else left = middle + 1;
//						}
//
//						nodes[nodn].relative_nods.insert(nodes[nodn].relative_nods.begin()+left, renod);
//Node_Num_Same: ;					
//					}
//				}
//			}
//		}

	//用于输出比较和检测
	//for(int i=0; i<(int)nodes.size(); i++)
	//{
	//	hout << i << "  " << (int)nodes[i].relative_nods.size() << "  ";
	//	for(int j=0; j<(int)nodes[i].relative_nods.size(); j++)
	//	{
	//		hout << nodes[i].relative_nods[j] << "  ";
	//	}
	//	hout << endl;
	//}

	return 1;
}
//---------------------------------------------------------------------------
//数据读入或者写入二进制文件，以备以后单独求解线性方程组时用到
int Algorithm_FEM::Write_read_matrix_equright_data(const long &Total_Matrix_Size, vector<double> &total_matrix, vector<double>* nl_equright, double backup_ege[], string &maeq_mod)
{
	if(maeq_mod=="write")
	{
		//------------------------------------------------------------------------
		//数据写入二进制文件
		fstream fout( "Matrix_Equright_Data", ios::out|ios::binary);
		for(long i=0; i<(long)total_matrix.size(); i++)
			fout.write((char *)&total_matrix[i], sizeof(double));
		for(int i=0; i<9; i++)
			for(int j=0; j<(int)nl_equright[i].size(); j++)
				fout.write((char *)&nl_equright[i][j], sizeof(double));
		for(int i=0; i<9; i++)
			fout.write((char *)&backup_ege[i], sizeof(double));
		fout.close();
	}
	else if(maeq_mod=="read")
	{
		//------------------------------------------------------------------------
		//数据读入二进制文件, 用于计算等效应变能时用
		fstream fin( "Matrix_Equright_Data", ios::in|ios::binary);
		for(long i=0; i<(long)total_matrix.size(); i++)
			fin.read((char *)&total_matrix[i], sizeof(double));
		for(int i=0; i<9; i++)
			for(int j=0; j<(int)nl_equright[i].size(); j++)
				fin.read((char *)&nl_equright[i][j], sizeof(double));
		for(int i=0; i<9; i++)
				fin.read((char *)&backup_ege[i], sizeof(double));
		fin.close();
	}
	else return 0;

	return 1;
}
//---------------------------------------------------------------------------
//输出或读取位移解
void Algorithm_FEM::read_u_solution(const string &output_file_name)
{
	//-------------------------------------------------------------------------------------------------
	//二进制读取数据
	int us, uis;
	fstream idata(output_file_name.c_str(), ios::in|ios::binary);
	idata.read((char *)&us, sizeof(int));
	idata.read((char *)&uis, sizeof(int));

	vector<double> temp_u(uis);
	U_Solution.assign(us, temp_u);
	for(int i=0; i<us; i++)
		for(int j=0; j<uis; j++)
			idata.read((char *)&U_Solution[i][j], sizeof(double));
	idata.close();
}
//---------------------------------------------------------------------------
//输出位移解及等效应变能
void Algorithm_FEM::write_u_solution_equivalent_energy(const string &output_file_name)const
{
	//-------------------------------------------------------------------------------------------------
	//二进制存储数据
	int us = (int)U_Solution.size();
	int uis = (int)U_Solution[0].size();
	fstream odata(output_file_name.c_str(), ios::out|ios::binary);
	odata.write((char *)&us, sizeof(int));
	odata.write((char *)&uis, sizeof(int));
	for(int i=0; i<us; i++)
		for(int j=0; j<uis; j++)
			odata.write((char *)&U_Solution[i][j], sizeof(double));

	odata.write((char *)equivalent_energy, sizeof(double)*9);

	odata.close();
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//读入一行信息，并跳过注释行（以"%"开头）；
string Algorithm_FEM::Get_Line(ifstream &infile)const
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
