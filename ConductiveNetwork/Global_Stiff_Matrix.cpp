//===========================================================================
// Global_Stiff_Matrix.cpp
// 求解总刚阵类成员函数
// Member Functions in a Class of the Global Stiff Matrix
//===========================================================================

#include "Global_Stiff_Matrix.h"

//---------------------------------------------------------------------------
//生成总体刚度矩阵
int Global_Stiff_Matrix::Gen_global_stiff_matrix(ifstream &infile, const vector<Point_3D> &cnps, const vector<Node> &nodes, const vector<Element> &elements, const vector<MatPro> &mats, const struct Decay_Para &decay, 
																						const struct RVE_Geo &cell, const struct CNT_Geo &cnts, const string &com_mod, const vector<long> &Iz, const vector<int> &Ig, vector<double> &total_matrix, vector<double>* nl_equright, double backup_ege[])
{
	//---------------------------------------------------------------------------
	//取高斯点
	Gauss gau;		//高斯点向量
	if(gau.Generate_gauss(infile)==0) return 0;

	clock_t ct0,ct1;
	//------------------------------------------------------------------------
	//计算单元每个高斯点的数据和相关纳米管权重
	ct0 = clock();
	hout << "-_- 开始计算高斯点的相关数据和权重" << endl;

	const int GS = (int)gau.gauss.size();
	const int ES = (int)elements.size();
	double (*gauss_ns)[8] = new double [GS][8];				//记录单元高斯点的形函数
	double (*gauss_nw)[8] = new double [GS][8];				//记录单元高斯点的形函数(带权重系数)
	double (*gauss_po)[3] = new double [GS][3];				//标准立方体单元的高斯点坐标
	double (*ele_cent)[3] = new double [ES][3];					//记录单元中心点位置向量(x,y,z)分别放于[0],[1]和[2]中（以标准立方体单元中心点为原点）
	double Jacobi = 0;																//标准立方体单元的雅可比值
	double (*gauss_dfx)[3][8] = new double [GS][3][8];		//记录单元高斯点形函数的导数
	Generate_element_gauss_data(elements, nodes, com_mod, gau.gauss, gau.weight, gauss_ns, gauss_nw, gauss_dfx, Jacobi, gauss_po, ele_cent);
	double (*gauss_wcnt)[10] = new double[ES*GS][10];	//纳米管的权重（[0]记录距离权重; [1]-[9]分别记录坐标轴旋转夹角余弦; 依次顺序为x'轴与x, y, z轴的夹角 [cos(a1), cos(b1), cos(c1)];
																								//y'轴与x, y, z轴的夹角[cos(a2), cos(b2), cos(c2)]; z'轴与x, y, z轴的夹角[cos(a3), cos(b3), cos(c3)]; z'轴为纳米管方向, x'轴为高斯点到纳米管垂直方向, y'=z'叉乘x'）
	if(Generate_element_gauss_weight_cnt(infile, decay.R, cell, cnts, cnps, elements, nodes, gau.gauss, gau.weight, gauss_po, ele_cent, gauss_wcnt, Jacobi)==0) return 0;

	//------------------------------------------------------------------------
	//数据写入二进制文件，以备后面计算等效应变能时用到
	//fstream fout( "GaussDats", ios::out|ios::binary);
	//fout.write((char *)&GS, sizeof(int));
	//fout.write((char *)&ES, sizeof(int));
	//for(int i=0; i<(int)gau.weight.size(); i++)
	//	fout.write((char *)&gau.weight[i], sizeof(double));

	//if(com_mod=="hybrid"||com_mod=="nonlocal") 
	//{
	//	fout.write((char *)gauss_po, sizeof(double)*GS*3);
	//	fout.write((char *)gauss_ns, sizeof(double)*GS*8);
	//	fout.write((char *)ele_cent, sizeof(double)*ES*3);
	//	fout.write((char *)gauss_wcnt, sizeof(double)*ES*GS*10);
	//}
	//if(com_mod=="hybrid"||com_mod=="fem") 
	//{
	//	fout.write((char *)&Jacobi, sizeof(double));
	//	fout.write((char *)gauss_dfx, sizeof(double)*GS*3*8);
	//}
	//fout.close();

	//------------------------------------------------------------------------
	ct1 = clock();
	hout << "    生成高斯点的相关数据和权重耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
	hout << "^_^ 计算高斯点的相关数据和权重操作完毕！" << endl << endl;

	//------------------------------------------------------------------------
	//循环所有单元
	ct0 = clock();
	hout << "-_- 开始计算单刚并加到总刚" << endl;
	//执行openmp
	#pragma omp parallel
	{
		//定义单刚
		double element_stiff_matrix1[24][24];
		double element_stiff_matrix2[24][24];
		double element_equright[9][24];
		double (*gwc_left)[10] = new double [GS][10];
		double (*gwc_right)[10] = new double [GS][10];

		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++)
		{
			//---------------------------------------------------------------------------
			//生成RVE背景网格单元(长程力作用)的单刚阵并添加到总刚
			if(com_mod=="hybrid"||com_mod=="nonlocal")
			{
				//开始计时
				clock_t ctn1,ctn2;
				ctn1 = clock();

				//--------------------------------------------------
				//初始化变量
				double elec_left[3] = { ele_cent[i][0], ele_cent[i][1], ele_cent[i][2] };

				//--------------------------------------------------
				//单元的每个高斯点关于最近距离纳米管的权重值
				const int IGS = i*GS;
				for(int j=0; j<GS; j++)
					for(int k=0; k<10; k++)
						gwc_left[j][k] = gauss_wcnt[IGS+j][k];

				//------------------------------------------------------------------------
				//长程力作用
				for(int j=0; j<(int)elements[i].relative_eles.size(); j++)
				{
					//--------------------------------------------------
					//初始化变量
					const int ere = elements[i].relative_eles[j];
					double elec_right[3] = { ele_cent[ere][0], ele_cent[ere][1], ele_cent[ere][2] };

					//判断那些由于周期边界条件产生的相关单元的情况
					bool mark = false;
					double peri_dist[3] = {0};
					if(elec_right[0]-elec_left[0] > cell.len_x-decay.R-Zero) { elec_right[0] -= cell.len_x; peri_dist[0] = -cell.len_x; mark = true; }
					else if(elec_right[0]-elec_left[0] < -cell.len_x+decay.R+Zero) { elec_right[0] += cell.len_x; peri_dist[0] = cell.len_x; mark = true; }

					if(elec_right[1]-elec_left[1] > cell.wid_y-decay.R-Zero) { elec_right[1] -= cell.wid_y; peri_dist[1] = -cell.wid_y; mark = true; }
					else if(elec_right[1]-elec_left[1] < -cell.wid_y+decay.R+Zero) { elec_right[1] += cell.wid_y; peri_dist[1] = cell.wid_y; mark = true; }

					if(elec_right[2]-elec_left[2] > cell.hei_z-decay.R-Zero) { elec_right[2] -= cell.hei_z; peri_dist[2] = -cell.hei_z; mark = true; }
					else if(elec_right[2]-elec_left[2] < -cell.hei_z+decay.R+Zero) { elec_right[2] += cell.hei_z; peri_dist[2] = cell.hei_z; mark = true; }

					//--------------------------------------------------
					//单元的每个高斯点关于最近距离纳米管的权重值
					const int EREGS = ere*GS;
					for(int k=0; k<GS; k++)
						for(int m=0; m<10; m++)
							gwc_right[k][m] = gauss_wcnt[EREGS+k][m];

					//------------------------------------------------------------------------
					//边界单元所以长程力相关单元是周期平移的体外单元
					if(mark)
					{
						double ege[9] = { 0 };
						Generate_Longforce_Equivalent_Equright(element_equright, Jacobi, gauss_nw, gau.weight, decay, gauss_po, elec_left, elec_right, gwc_left, gwc_right, peri_dist, ege);
						//单刚添加到右端项
						#pragma omp critical
						{
							for(int k=0; k<9; k++)
								for(int m=0; m<8; m++)
									for(int p=0; p<3; p++)
										nl_equright[k][3*elements[i].nodes_id[m]+p] += element_equright[k][3*m+p];

							for(int k=0; k<9; k++)	backup_ege[k] += ege[k];
						}
					}

					//------------------------------------------------------------------------
					//生成RVE背景网格单元(长程力作用)单刚阵
					Generate_Longforce_Elestiff(element_stiff_matrix1, element_stiff_matrix2, gauss_ns, gauss_nw, gau.weight, decay, gauss_po, elec_left, elec_right, gwc_left, gwc_right);

					//单刚添加到总刚
					#pragma omp critical
					{
						Add_to_gsmatrix(element_stiff_matrix1, Iz, Ig, total_matrix, elements[i]);	//由于在openmp中, 所以不能返回0值, 终止程序
						Add_to_gsmatrix(element_stiff_matrix2, Iz, Ig, total_matrix, elements[i], elements[elements[i].relative_eles[j]]);	//由于在openmp中, 所以不能返回0值, 终止程序
					}
				}

				//------------------------------------------------------------------------
				ctn2 = clock();
				hout << "Total num of elements: " << (int)elements.size() << "; Element " << i << " took time: " << (double)(ctn2-ctn1)/CLOCKS_PER_SEC << "sec; " << endl;
			}

			//---------------------------------------------------------------------------
			//生成RVE背景网格单元(接触力作用)的单刚阵并添加到总刚
			//生成RVE网格单元(接触力作用)单刚阵（经典有限元单刚阵）
			if(com_mod=="hybrid"||com_mod=="fem")
			{
				if(com_mod=="fem")  //局部连续, 纳米管增强部分
				{
					//--------------------------------------------------
					//单元的每个高斯点关于最近距离纳米管的权重值
					const int IGS = i*GS;
					for(int j=0; j<GS; j++)
						for(int k=0; k<10; k++)
							gwc_left[j][k] = gauss_wcnt[IGS+j][k];

					CNT_Reinforcement_Local_Continuum(element_stiff_matrix1, mats, Jacobi, gauss_dfx, gau.weight, gwc_left);
					
					//单刚添加到总刚
					#pragma omp critical
					{					
						Add_to_gsmatrix(element_stiff_matrix1, Iz, Ig, total_matrix, elements[i]);  //由于在openmp中, 所以不能返回0值, 终止程序
					}
				}

				Generate_Contactforce_Elestiff(element_stiff_matrix1, mats, Jacobi, gauss_dfx, gau.weight);
				//单刚添加到总刚
				#pragma omp critical
				{
					Add_to_gsmatrix(element_stiff_matrix1, Iz, Ig, total_matrix, elements[i]);  //由于在openmp中, 所以不能返回0值, 终止程序
				}
			}
		}
		//---------------------------------------------------------------------------
		//删除指针
		delete[] gwc_left;
		delete[] gwc_right;
	}

	//---------------------------------------------------------------------------
	//删除指针
	delete[] gauss_ns;
	delete[] gauss_nw;
	delete[] gauss_dfx;
	delete[] gauss_po;
	delete[] ele_cent;
	delete[] gauss_wcnt;

	//---------------------------------------------------------------------------
	ct1 = clock();
	hout << "    计算单刚并加到总刚耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
	hout << "^_^ 计算单刚并加到总刚操作完毕！" << endl << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//生成RVE背景网格单元(接触力作用)单刚阵(六面体)（CNT增强界面层局部连续模型部分）
void Global_Stiff_Matrix::CNT_Reinforcement_Local_Continuum(double (*element_stiff_matrix)[24], const vector<MatPro> &mats, const double &Jacobi, const double (*gauss_dfx)[3][8], 
																													 const vector<double> &weight, const double (*gwc_left)[10])const
{
	//--------------------------------------------------
	//初始化
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++) 
			element_stiff_matrix[i][j] = 0;
	//--------------------------------------------------
	//循环高斯点计算积分
	const int ws = (int)weight.size();
	for(int count=0; count<ws; count++)
	{
		if(gwc_left[count][0]==0) continue;
		//-----------------------------------------------------------------
		//刚度矩阵的旋转变换矩阵
		double Ro[6][6];
		Ro[0][0] = gwc_left[count][1]*gwc_left[count][1];
		Ro[0][1] = gwc_left[count][2]*gwc_left[count][2];
		Ro[0][2] = gwc_left[count][3]*gwc_left[count][3];
		Ro[0][3] = 2*gwc_left[count][1]*gwc_left[count][2];
		Ro[0][4] = 2*gwc_left[count][2]*gwc_left[count][3];
		Ro[0][5] = 2*gwc_left[count][1]*gwc_left[count][3];

		Ro[1][0] = gwc_left[count][4]*gwc_left[count][4];
		Ro[1][1] = gwc_left[count][5]*gwc_left[count][5];
		Ro[1][2] = gwc_left[count][6]*gwc_left[count][6];
		Ro[1][3] = 2*gwc_left[count][4]*gwc_left[count][5];
		Ro[1][4] = 2*gwc_left[count][5]*gwc_left[count][6];
		Ro[1][5] = 2*gwc_left[count][4]*gwc_left[count][6];

		Ro[2][0] = gwc_left[count][7]*gwc_left[count][7];
		Ro[2][1] = gwc_left[count][8]*gwc_left[count][8];
		Ro[2][2] = gwc_left[count][9]*gwc_left[count][9];
		Ro[2][3] = 2*gwc_left[count][7]*gwc_left[count][8];
		Ro[2][4] = 2*gwc_left[count][8]*gwc_left[count][9];
		Ro[2][5] = 2*gwc_left[count][7]*gwc_left[count][9];

		Ro[3][0] = gwc_left[count][1]*gwc_left[count][4];
		Ro[3][1] = gwc_left[count][2]*gwc_left[count][5];
		Ro[3][2] = gwc_left[count][3]*gwc_left[count][6];
		Ro[3][3] = gwc_left[count][1]*gwc_left[count][5] + gwc_left[count][2]*gwc_left[count][4];
		Ro[3][4] = gwc_left[count][2]*gwc_left[count][6] + gwc_left[count][3]*gwc_left[count][5];
		Ro[3][5] = gwc_left[count][1]*gwc_left[count][6] + gwc_left[count][3]*gwc_left[count][4];

		Ro[4][0] = gwc_left[count][4]*gwc_left[count][7];
		Ro[4][1] = gwc_left[count][5]*gwc_left[count][8];
		Ro[4][2] = gwc_left[count][6]*gwc_left[count][9];
		Ro[4][3] = gwc_left[count][4]*gwc_left[count][8] + gwc_left[count][5]*gwc_left[count][7];
		Ro[4][4] = gwc_left[count][5]*gwc_left[count][9] + gwc_left[count][6]*gwc_left[count][8];
		Ro[4][5] = gwc_left[count][4]*gwc_left[count][9] + gwc_left[count][6]*gwc_left[count][7];

		Ro[5][0] = gwc_left[count][1]*gwc_left[count][7];
		Ro[5][1] = gwc_left[count][2]*gwc_left[count][8];
		Ro[5][2] = gwc_left[count][3]*gwc_left[count][9];
		Ro[5][3] = gwc_left[count][1]*gwc_left[count][8] + gwc_left[count][2]*gwc_left[count][7];
		Ro[5][4] = gwc_left[count][2]*gwc_left[count][9] + gwc_left[count][3]*gwc_left[count][8];
		Ro[5][5] = gwc_left[count][1]*gwc_left[count][9] + gwc_left[count][3]*gwc_left[count][7];

		//-----------------------------------------------------------------
		//计算此单元所对应的材料弹性矩阵
		double temp_ele[6][6];
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
			{
				temp_ele[i][j] = 0.0;
				for(int k=0; k<6; k++)
				temp_ele[i][j] += Ro[i][k]*mats[1].elas_matrix[k][j];
			}

		//-----------------------------------------------------------------
		//乘以Ro的转置
		double ele_elas[6][6];
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
			{
				ele_elas[i][j] = 0.0;
				for(int k=0; k<6; k++)
				ele_elas[i][j] += temp_ele[i][k]*Ro[j][k];
			}

		//--------------------------------------------------------
		//B矩阵
		double B[6][24] = {{0}, {0}, {0}, {0}, {0}, {0}};

		for(int i=0; i<8; i++)
		{
			B[0][i*3+0]=gauss_dfx[count][0][i];
			B[1][i*3+1]=gauss_dfx[count][1][i];
			B[2][i*3+2]=gauss_dfx[count][2][i];
			B[3][i*3+0]=gauss_dfx[count][1][i];
			B[3][i*3+1]=gauss_dfx[count][0][i];
			B[4][i*3+1]=gauss_dfx[count][2][i];
			B[4][i*3+2]=gauss_dfx[count][1][i];
			B[5][i*3+0]=gauss_dfx[count][2][i];
			B[5][i*3+2]=gauss_dfx[count][0][i];
		}

		//--------------------------------------------------------
		//计算B_trans
		double B_trans[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
				B_trans[i][j] = B[j][i];
		//--------------------------------------------------------------------------------------------
		//求出B_trans矩阵与ele_elas矩阵的乘积array1
		double array1[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
			{
				array1[i][j]=0; 
				for(int k=0; k<6; k++)
					array1[i][j] += B_trans[i][k]*ele_elas[k][j];
			}
		//求出array1矩阵与B矩阵的乘积array2
		double array2[24][24];
		for(int i=0; i<24; i++)
			for(int j=0; j<24; j++)
			{
				array2[i][j]=0;
				for(int k=0; k<6; k++)
					array2[i][j] += array1[i][k]*B[k][j];
				element_stiff_matrix[i][j] += array2[i][j]*weight[count]*gwc_left[count][0];
			}
	}

	//在高斯点循环外乘雅可比值减少乘法次数
	for(int j=0; j<24; j++) 
		for (int k=0; k<24; k++) 
			element_stiff_matrix[j][k] = element_stiff_matrix[j][k]*Jacobi;

	//输出单刚阵，用于检查
	//hout << "element stiffness matrix: " << endl;
	//for(int i=0; i<24; i++){
	//	for(int j=0; j<24; j++)
	//		hout << element_stiff_matrix[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl << endl;
}

//-----------------------------------------------------------------------------------------------
//生成RVE背景网格单元(接触力作用)单刚阵(六面体)
void Global_Stiff_Matrix::Generate_Contactforce_Elestiff(double (*element_stiff_matrix)[24], const vector<MatPro> &mats, const double &Jacobi, const double (*gauss_dfx)[3][8], const vector<double> &weight)const
{
	//--------------------------------------------------
	//初始化
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++) 
			element_stiff_matrix[i][j] = 0;
	//--------------------------------------------------
	//循环高斯点计算积分
	const int ws = (int)weight.size();
	for(int count=0; count<ws; count++)
	{
		//-----------------------------------------------------------------
		//计算此单元所对应的材料弹性矩阵
		double ele_elas[6][6];		
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				ele_elas[i][j] = mats[0].elas_matrix[i][j];
		
		//--------------------------------------------------------
		//B矩阵
		double B[6][24] = {{0}, {0}, {0}, {0}, {0}, {0}};

		for(int i=0; i<8; i++)
		{
			B[0][i*3+0]=gauss_dfx[count][0][i];
			B[1][i*3+1]=gauss_dfx[count][1][i];
			B[2][i*3+2]=gauss_dfx[count][2][i];
			B[3][i*3+0]=gauss_dfx[count][1][i];
			B[3][i*3+1]=gauss_dfx[count][0][i];
			B[4][i*3+1]=gauss_dfx[count][2][i];
			B[4][i*3+2]=gauss_dfx[count][1][i];
			B[5][i*3+0]=gauss_dfx[count][2][i];
			B[5][i*3+2]=gauss_dfx[count][0][i];
		}

		//--------------------------------------------------------
		//计算B_trans
		double B_trans[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
				B_trans[i][j] = B[j][i];
		//--------------------------------------------------------------------------------------------
		//求出B_trans矩阵与ele_elas矩阵的乘积array1
		double array1[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
			{
				array1[i][j]=0; 
				for(int k=0; k<6; k++)
					array1[i][j] += B_trans[i][k]*ele_elas[k][j];
			}
		//求出array1矩阵与B矩阵的乘积array2
		double array2[24][24];
		for(int i=0; i<24; i++)
			for(int j=0; j<24; j++)
			{
				array2[i][j]=0;
				for(int k=0; k<6; k++)
					array2[i][j] += array1[i][k]*B[k][j];
				element_stiff_matrix[i][j] += array2[i][j]*weight[count];
			}
	}

	//在高斯点循环外乘雅可比值减少乘法次数
	for(int j=0; j<24; j++) 
		for (int k=0; k<24; k++) 
			element_stiff_matrix[j][k] = element_stiff_matrix[j][k]*Jacobi;

	//输出单刚阵，用于检查
	//hout << "element stiffness matrix: " << endl;
	//for(int i=0; i<24; i++){
	//	for(int j=0; j<24; j++)
	//		hout << element_stiff_matrix[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl << endl;
}
//-----------------------------------------------------------------------------------------------
//生成RVE背景网格单元(长程力作用)单刚阵(六面体)
void Global_Stiff_Matrix::Generate_Longforce_Elestiff(double (*element_stiff_matrix1)[24], double (*element_stiff_matrix2)[24], const double (*gauss_ns)[8], const double (*gauss_nw)[8], const vector<double> &weight, 
																								   const struct Decay_Para &decay, const double (*gauss_po)[3], const double elec_left[], const double elec_right[], const double (*gwc_left)[10], const double (*gwc_right)[10])const
{
	//--------------------------------------------------
	//初始化
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++)
		{
			element_stiff_matrix1[i][j] = 0;
			element_stiff_matrix2[i][j] = 0;
		}

	//--------------------------------------------------	
	//循环本单元高斯点计算积分
	const int gau_num = (int)weight.size();
	for(int count1=0; count1<gau_num; count1++)
	{
		//--------------------------------------------------	
		double Gmat[6] = {0};
		double GNmat[6][8] = { {0}, {0}, {0}, {0}, {0}, {0} };
		//--------------------------------------------------	
		//左端高斯点坐标
		Point_3D gaupoi_left(0, 0, 0);
		gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
		gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
		gaupoi_left.z = gauss_po[count1][2] + elec_left[2];

		//------------------------------------------------------------------------------------------------------------------------
		//循环外单元高斯点计算积分
		for(int count2=0; count2<gau_num; count2++)
		{
			//--------------------------------------------
			//右端高斯点坐标
			Point_3D gaupoi_right(0, 0, 0);
			gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
			gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
			gaupoi_right.z = gauss_po[count2][2] + elec_right[2];
	
			//--------------------------------------------
			//权重值判断
			if(gwc_left[count1][0]<Zero&&gwc_right[count2][0]<Zero) continue;  //等于零值，没有长程效果

			//--------------------------------------------
			//计算长程作用衰减函数值
			const double xleft = gaupoi_right.x-gaupoi_left.x;
			const double yleft = gaupoi_right.y-gaupoi_left.y;
			const double zleft = gaupoi_right.z-gaupoi_left.z;

			//--------------------------------------------
			//计算bond的长度
			const double dis_squr = xleft*xleft + yleft*yleft + zleft*zleft;
			const double poi_dis = sqrt(dis_squr);		//两点之间的距离
			
			if(poi_dis>decay.R+Zero||poi_dis<Zero) continue;		//圆形积分域

			//--------------------------------------------
			//计算left点权重
			double sum_left = 0;
			if(gwc_left[count1][0]>Zero)
			{
				const double x = gwc_left[count1][1]*xleft + gwc_left[count1][2]*yleft + gwc_left[count1][3]*zleft;
				const double y = gwc_left[count1][4]*xleft + gwc_left[count1][5]*yleft + gwc_left[count1][6]*zleft;
				const double z = gwc_left[count1][7]*xleft + gwc_left[count1][8]*yleft + gwc_left[count1][9]*zleft;

				const double cos2sita = z*z/dis_squr;		//cos(sita)^2
				double cos2pha;
				if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
				else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

				sum_left = decay.acoe[0] + decay.acoe[1]*0.5*(3*cos2sita-1) + decay.acoe[2]*(2*cos2pha-1)*3*(1-cos2sita) + decay.acoe[3]*0.125*(35*cos2sita*cos2sita - 30*cos2sita + 3)
									+ decay.acoe[4]*(2*cos2pha-1)*7.5*(7*cos2sita-1)*(1-cos2sita) + decay.acoe[5]*(8*cos2pha*cos2pha-8*cos2pha+1)*105*(1-cos2sita)*(1-cos2sita);
			}

			//--------------------------------------------
			//计算right点权重
			double sum_right = 0;
			if(gwc_right[count2][0]>Zero)
			{
				const double xright = gaupoi_left.x-gaupoi_right.x;
				const double yright = gaupoi_left.y-gaupoi_right.y;
				const double zright = gaupoi_left.z-gaupoi_right.z;

				const double x = gwc_right[count2][1]*xright + gwc_right[count2][2]*yright + gwc_right[count2][3]*zright;
				const double y = gwc_right[count2][4]*xright + gwc_right[count2][5]*yright + gwc_right[count2][6]*zright;
				const double z = gwc_right[count2][7]*xright + gwc_right[count2][8]*yright + gwc_right[count2][9]*zright;

				const double cos2sita = z*z/dis_squr;		//cos(sita)^2
				double cos2pha;
				if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
				else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

				sum_right = decay.acoe[0] + decay.acoe[1]*0.5*(3*cos2sita-1) + decay.acoe[2]*(2*cos2pha-1)*3*(1-cos2sita) + decay.acoe[3]*0.125*(35*cos2sita*cos2sita - 30*cos2sita + 3)
									+ decay.acoe[4]*(2*cos2pha-1)*7.5*(7*cos2sita-1)*(1-cos2sita) + decay.acoe[5]*(8*cos2pha*cos2pha-8*cos2pha+1)*105*(1-cos2sita)*(1-cos2sita);
			}

			//--------------------------------------------
			//计算权重函数值
			const double sum = 0.5*(gwc_left[count1][0]*sum_left+gwc_right[count2][0]*sum_right);
			if(fabs(sum)<Zero) continue;  //等于零值，没有长程效果

			//--------------------------------------------
			//计算长程力矩阵
			const double gv = exp(-poi_dis/decay.radius)*sum*weight[count2]; //衰减以及权重函数值

			const double temp_gmat[6] = {gv*xleft*xleft, gv*yleft*yleft, gv*zleft*zleft, gv*xleft*yleft, gv*yleft*zleft, gv*zleft*xleft}; //衰减函数矩阵的对称项(用xleft或者xright的乘积都是一样的)

			for(int i=0; i<6; i++) Gmat[i] += temp_gmat[i];  //关于count2叠加

			for(int i=0; i<6; i++)
				for(int j=0; j<8; j++) GNmat[i][j] +=  temp_gmat[i]*gauss_ns[count2][j]; //关于count2并与形函数乘积叠加
		}

		//------------------------------------------------------------------------------------------------------------------------
		int ni = 0, nj = 0;
		double temxy = 0, nxy[6] = {0};
		for(int i=0; i<8; i++)
		{
			nj = 0;
			for(int j=0; j<8; j++)
			{
				temxy = gauss_nw[count1][i]*gauss_ns[count1][j];
				for(int k=0; k<6; k++)
				nxy[k] = temxy*Gmat[k];

				element_stiff_matrix1[ni][nj] += nxy[0];
				element_stiff_matrix1[ni+1][nj+1] += nxy[1];
				element_stiff_matrix1[ni+2][nj+2] += nxy[2];
				element_stiff_matrix1[ni][nj+1] += nxy[3];
				element_stiff_matrix1[ni+1][nj] += nxy[3];
				element_stiff_matrix1[ni+1][nj+2] += nxy[4];
				element_stiff_matrix1[ni+2][nj+1] += nxy[4];
				element_stiff_matrix1[ni][nj+2] += nxy[5];
				element_stiff_matrix1[ni+2][nj] += nxy[5];

				for(int k=0; k<6; k++)
				nxy[k] = gauss_nw[count1][i]*GNmat[k][j];

				element_stiff_matrix2[ni][nj] -= nxy[0];
				element_stiff_matrix2[ni+1][nj+1] -= nxy[1];
				element_stiff_matrix2[ni+2][nj+2] -= nxy[2];
				element_stiff_matrix2[ni][nj+1] -= nxy[3];
				element_stiff_matrix2[ni+1][nj] -= nxy[3];
				element_stiff_matrix2[ni+1][nj+2] -= nxy[4];
				element_stiff_matrix2[ni+2][nj+1] -= nxy[4];
				element_stiff_matrix2[ni][nj+2] -= nxy[5];
				element_stiff_matrix2[ni+2][nj] -= nxy[5];

				nj += 3;
			}
			ni += 3;
		}
	}

	////输出单刚阵，用于检查
	//hout << "element stiffness matrix1: " << endl;
	//for(int i=0; i<24; i++){
	//	for(int j=0; j<24;j++)
	//		hout << element_stiff_matrix1[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl;

	//hout << "element stiffness matrix2: " << endl;
	//for(int i=0; i<24; i++){
	//	for(int j=0; j<24;j++)
	//		hout << element_stiff_matrix2[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl;
}
//-----------------------------------------------------------------------------------------------
//生成边界单元(长程力作用)单刚阵(六面体)等效方程右端向量
void Global_Stiff_Matrix::Generate_Longforce_Equivalent_Equright(double (*element_equright)[24], const double &Jacobi, const double (*gauss_nw)[8], const vector<double> &weight, const struct Decay_Para &decay,
																														const double (*gauss_po)[3], const double elec_left[], const double elec_right[], const double (*gwc_left)[10], const double (*gwc_right)[10], const double peri_dist[], double ege[])const
{
	//--------------------------------------------------
	//初始化
	for(int i=0; i<9; i++)
		for(int j=0; j<24; j++)
			element_equright[i][j] = 0;

	//--------------------------------------------------	
	//循环本单元高斯点计算积分
	double TGmat[6] = { 0 };
	double NGmat[8][6] = { {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0} };
	const int gau_num = (int)weight.size();
	for(int count1=0; count1<gau_num; count1++)
	{
		//--------------------------------------------------	
		double Gmat[6] = {0};

		//--------------------------------------------------	
		//左端高斯点坐标
		Point_3D gaupoi_left(0, 0, 0);
		gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
		gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
		gaupoi_left.z = gauss_po[count1][2] + elec_left[2];

		//------------------------------------------------------------------------------------------------------------------------
		//循环外单元高斯点计算积分
		for(int count2=0; count2<gau_num; count2++)
		{
			//--------------------------------------------
			//右端高斯点坐标
			Point_3D gaupoi_right(0, 0, 0);
			gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
			gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
			gaupoi_right.z = gauss_po[count2][2] + elec_right[2];
	
			//--------------------------------------------
			//权重值判断
			if(gwc_left[count1][0]<Zero&&gwc_right[count2][0]<Zero) continue;  //等于零值，没有长程效果

			//--------------------------------------------
			//计算长程作用衰减函数值
			const double xleft = gaupoi_right.x-gaupoi_left.x;
			const double yleft = gaupoi_right.y-gaupoi_left.y;
			const double zleft = gaupoi_right.z-gaupoi_left.z;

			//--------------------------------------------
			//计算bond的长度
			const double dis_squr = xleft*xleft + yleft*yleft + zleft*zleft;
			const double poi_dis = sqrt(dis_squr);		//两点之间的距离
			
			if(poi_dis>decay.R+Zero||poi_dis<Zero) continue;		//圆形积分域

			//--------------------------------------------
			//计算left点权重
			double sum_left = 0;
			if(gwc_left[count1][0]>Zero)
			{
				const double x = gwc_left[count1][1]*xleft + gwc_left[count1][2]*yleft + gwc_left[count1][3]*zleft;
				const double y = gwc_left[count1][4]*xleft + gwc_left[count1][5]*yleft + gwc_left[count1][6]*zleft;
				const double z = gwc_left[count1][7]*xleft + gwc_left[count1][8]*yleft + gwc_left[count1][9]*zleft;

				const double cos2sita = z*z/dis_squr;		//cos(sita)^2
				double cos2pha;
				if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
				else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

				sum_left = decay.acoe[0] + decay.acoe[1]*0.5*(3*cos2sita-1) + decay.acoe[2]*(2*cos2pha-1)*3*(1-cos2sita) + decay.acoe[3]*0.125*(35*cos2sita*cos2sita - 30*cos2sita + 3)
									+ decay.acoe[4]*(2*cos2pha-1)*7.5*(7*cos2sita-1)*(1-cos2sita) + decay.acoe[5]*(8*cos2pha*cos2pha-8*cos2pha+1)*105*(1-cos2sita)*(1-cos2sita);
			}

			//--------------------------------------------
			//计算right点权重
			double sum_right = 0;
			if(gwc_right[count2][0]>Zero)
			{
				const double xright = gaupoi_left.x-gaupoi_right.x;
				const double yright = gaupoi_left.y-gaupoi_right.y;
				const double zright = gaupoi_left.z-gaupoi_right.z;

				const double x = gwc_right[count2][1]*xright + gwc_right[count2][2]*yright + gwc_right[count2][3]*zright;
				const double y = gwc_right[count2][4]*xright + gwc_right[count2][5]*yright + gwc_right[count2][6]*zright;
				const double z = gwc_right[count2][7]*xright + gwc_right[count2][8]*yright + gwc_right[count2][9]*zright;

				const double cos2sita = z*z/dis_squr;		//cos(sita)^2
				double cos2pha;
				if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
				else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

				sum_right = decay.acoe[0] + decay.acoe[1]*0.5*(3*cos2sita-1) + decay.acoe[2]*(2*cos2pha-1)*3*(1-cos2sita) + decay.acoe[3]*0.125*(35*cos2sita*cos2sita - 30*cos2sita + 3)
									+ decay.acoe[4]*(2*cos2pha-1)*7.5*(7*cos2sita-1)*(1-cos2sita) + decay.acoe[5]*(8*cos2pha*cos2pha-8*cos2pha+1)*105*(1-cos2sita)*(1-cos2sita);
			}

			//--------------------------------------------
			//计算权重函数值
			const double sum = 0.5*(gwc_left[count1][0]*sum_left+gwc_right[count2][0]*sum_right);
			if(fabs(sum)<Zero) continue;  //等于零值，没有长程效果

			//--------------------------------------------
			//计算长程力矩阵			
			const double gv = exp(-poi_dis/decay.radius)*sum*weight[count2]; //衰减以及权重函数值

			//------------------------------------------------------------------------------------------------
			//注意这里 temp_gmat[6]的赋值顺序不同于以往长程力计算时temp_gmat[6]的顺序
			const double temp_gmat[6] = { gv*xleft*xleft, gv*xleft*yleft, gv*yleft*yleft, gv*xleft*zleft, gv*yleft*zleft, gv*zleft*zleft }; //衰减函数矩阵的对称项

			for(int i=0; i<6; i++) Gmat[i] += temp_gmat[i];  //关于count2叠加
		}

		for(int i=0; i<8; i++)
			for(int j=0; j<6; j++)
				NGmat[i][j] += gauss_nw[count1][i]*Gmat[j];

		for(int i=0; i<6; i++) TGmat[i] += Gmat[i]*weight[count1];
	}

	for(int i=0; i<9; i++)
	{
		double E[3][3] = {{0}, {0}, {0}}; //用于周期边界条件处理时的变形矩阵条件
		switch(i)	//设置均匀化应变做为周期边界条件
		{
		case 0: E[0][0]=0.1; break;
		case 1: E[1][1]=0.1; break;
		case 2: E[2][2]=0.1; break;
		case 3: E[0][1]=0.05; E[1][0]=0.05; break;
		case 4: E[1][2]=0.05; E[2][1]=0.05; break;
		case 5: E[0][2]=0.05; E[2][0]=0.05; break;
		case 6: E[0][0]=0.1; E[1][1]=0.1; break;
		case 7: E[1][1]=0.1; E[2][2]=0.1; break;
		case 8: E[0][0]=0.1; E[2][2]=0.1; break;
		default: hout << "错误！ 处理周期边界条件约束时，循环次数值等于" << i << "小于0或者大于8！请检查！" << endl;
		}

		//均匀化位移差
		double uni_disp[3] = {0};
		for(int j=0; j<3; j++)
			for(int k=0; k<3; k++) 
				uni_disp[j] += E[j][k]*peri_dist[k];

		for(int j=0; j<3; j++) uni_disp[j] = uni_disp[j]*Jacobi; //为了减少乘法的次数，这里先乘上雅可比常数

		for(int j=0; j<8; j++)
		{
			for(int k=0; k<=2; k++)
				for(int m=0; m<=2; m++)
				{
					if(k==2||m==2)	
						element_equright[i][3*j+k] += NGmat[j][k+m+1]*uni_disp[m];
					else 
						element_equright[i][3*j+k] += NGmat[j][k+m]*uni_disp[m];
				}
		}

		double Tu[3] = { 0 };
		for(int j=0; j<=2; j++)
			for(int k=0; k<=2; k++)
			{
				if(j==2||k==2)	
					Tu[j] += TGmat[j+k+1]*uni_disp[k];
				else 
					Tu[j] += TGmat[j+k]*uni_disp[k];
			}

		for(int j=0; j<=2; j++) ege[i] += uni_disp[j]*Tu[j];
	}
}
//-----------------------------------------------------------------------------------------------
//将单刚添加到总刚
void Global_Stiff_Matrix::Add_to_gsmatrix(const double (*element_stiff_matrix)[24], const vector<long> &Iz, const vector<int> &Ig, vector<double> &total_matrix, const Element &element)const
{
	for(int i=0; i<(int)element.nodes_id.size(); i++)
		for(int j=0; j<(int)element.nodes_id.size(); j++)
		{
			const int row_node = element.nodes_id[i];
			const int col_node = element.nodes_id[j];
			if(row_node < col_node) continue;
			if(row_node > col_node )
			{
				//二分法查找位置
				bool mark = false;
				long left = Iz[row_node-1];
				long middle = 0;
				long right = Iz[row_node];
				while(right>=left)
				{
					middle = (left + right)/2;
					if(Ig[middle] == col_node) { mark = true; break; } //找到对应节点
					else if(Ig[middle] > col_node) right = middle - 1;
					else left = middle + 1;
				}					
				if(!mark) { hout << "错误！在节点" << row_node << "的Iz, Ig中没有找到节点" << col_node << "请检查！" << endl; }
				
				const long Mt = 6*(long)row_node + 9*middle;	//找到对应一维存储的起始位置
if( Mt>=(long)total_matrix.size() ) hout << "$$$1" <<  Mt << " row_now: " << row_node << " middle: " << middle << "  col_node: " << col_node << " Iz[row_node-1]: "  <<  Iz[row_node-1] << " Iz[row_node]: "  <<  Iz[row_node]<< endl;
				for(int k=0; k<=2; k++)
					for(int m=0; m<=2; m++)
						total_matrix[Mt+3*k+m] += element_stiff_matrix[3*i+k][3*j+m];
			}
			else
			{
				const long Mt = 6*(long)row_node + 9*Iz[row_node]; //找到对应一维存储的起始位置
if( Mt>=(long)total_matrix.size() ) hout << "$$$2" <<  Mt << " row_now: " << row_node << "  col_node: " << col_node << " Iz[row_node-1]: "  <<  Iz[row_node-1] << " Iz[row_node]: "  <<  Iz[row_node] << endl;							
				for(int k=0; k<=2; k++)
					for(int m=0; m<=k; m++)	//注意此处是小于等于k, 因为是下三角阵
					{
						if(k==2)
							total_matrix[Mt+k+m+1] += element_stiff_matrix[3*i+k][3*j+m];
						else
							total_matrix[Mt+k+m] += element_stiff_matrix[3*i+k][3*j+m];
					}
			}			             
		}
}
//-----------------------------------------------------------------------------------------------
//将单刚添加到总刚(行列单元节点不同)
void Global_Stiff_Matrix::Add_to_gsmatrix(const double (*element_stiff_matrix)[24], const vector<long> &Iz, const vector<int> &Ig, vector<double> &total_matrix, const Element &ele_row, const Element &ele_col)const
{
	for(int i=0; i<(int)ele_row.nodes_id.size(); i++)
		for(int j=0; j<(int)ele_col.nodes_id.size(); j++)
		{
			const int row_node = ele_row.nodes_id[i];
			const int col_node = ele_col.nodes_id[j];
			if(row_node < col_node) continue;
			if(row_node > col_node )
			{
				//二分法查找位置
				bool mark = false;
				long left = Iz[row_node-1];
				long middle = 0;
				long right = Iz[row_node];
				while(right>=left)
				{
					middle = (left + right)/2;
					if(Ig[middle] == col_node) { mark = true; break; } //找到对应节点
					else if(Ig[middle] > col_node) right = middle - 1;
					else left = middle + 1;
				}					
				if(!mark) { hout << "错误！在节点" << row_node << "的Iz, Ig中没有找到节点" << col_node << "请检查！" << endl; }
				
				const long Mt = 6*(long)row_node + 9*middle;	//找到对应一维存储的起始位置

				for(int k=0; k<=2; k++)
					for(int m=0; m<=2; m++)
						total_matrix[Mt+3*k+m] += element_stiff_matrix[3*i+k][3*j+m];
			}
			else
			{
				const long Mt = 6*(long)row_node + 9*Iz[row_node]; //找到对应一维存储的起始位置
				
				for(int k=0; k<=2; k++)
					for(int m=0; m<=k; m++)	//注意此处是小于等于k, 因为是下三角阵
					{
						if(k==2)
							total_matrix[Mt+k+m+1] += element_stiff_matrix[3*i+k][3*j+m];
						else
							total_matrix[Mt+k+m] += element_stiff_matrix[3*i+k][3*j+m];
					}
			}			             
		}
}
//-----------------------------------------------------------------------------------------------
//计算单元每个高斯点关于最近距离纳米管的权重值
int Global_Stiff_Matrix::Generate_element_gauss_weight_cnt(ifstream &infile, const double &dist, const struct RVE_Geo &cell, const struct CNT_Geo &cnts, const vector<Point_3D> &cnps, const vector<Element> &elements,
																											 const vector<Node> &nodes, const vector<Node> &gauss, const vector<double> &weight, const double (*gauss_po)[3], const double (*ele_cent)[3], double (*gauss_wcnt)[10], double &Jacobi)const
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//读入权重函数的作用距离以及权重函数的阶数
	double w_dist, w_constant;
	int w_order;
	istringstream istrw(Get_Line(infile));
	istrw >> w_dist;
	if(w_dist<0&&w_dist!=-1) { hout << "注意！权重函数的最大作用距离等于" << w_dist << "，小于0且不等于-1（特殊情况），请重新输入！" << endl; return 0; }
	if(w_dist>=0)
	{
		istrw >> w_order;
		if(w_order>3||w_order<-2) { hout << "注意！权重函数的幂次等于" << w_order << "，高于3次幂或者低于-2次幂，请重新输入！" << endl; return 0; }
		if(w_order==0)
		{
			istrw >> w_constant;		// 如果是0次幂表示常数，这里还需要输入常数值
			if(w_constant>1.0||w_constant<0.0) { hout << "注意！权重函数的幂次等于0，但是常数系数输入有误，大于1或者小于0或者根本没有输入！" << endl; return 0; }
		}
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//决定单元的相关单元, 用于计算高斯点距纳米管的最小距离
	vector<vector<int> > relemd;
	Deter_relative_elements_min_dist(nodes, elements, cell, w_dist, relemd);

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//计算高斯点到最近纳米管的距离和方向
	double total_cnt_volume = 0;
	double total_interface_volume = 0;
	double total_value_gps = 0;
	#pragma omp parallel
	{
		//定义临时权重变量
		double temp_cnt_volume =0;
		double temp_interface_volume = 0;
		vector<double> temp_gwcnt[10];
		vector<int> temp_num;

		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<(int)elements.size(); i++)
		{
			const int GS = (int)gauss.size();
			for(int count=0; count<GS; count++)
			{	
				temp_num.push_back(i*GS+count);  //记录高斯节点在一维存储中的整体编号，在omp后存储是用到
				//特殊情况, 每一个高斯点都是权重都是1，且旋转阵为单位阵
				if(w_dist==-1) 
				{ 
					temp_gwcnt[0].push_back(1.0);
					//保证旋转矩阵时单位矩阵
					for(int j=0; j<3; j++)
						for(int k=0; k<3; k++)
						{
							int tgn = j*3+k+1;
							if(j==k) temp_gwcnt[tgn].push_back(1.0);
							else temp_gwcnt[tgn].push_back(0.0);
						}
					continue; 
				}

				//--------------------------------------------
				Point_3D gaupoi(0, 0, 0);		//高斯积分点
				gaupoi.x = gauss_po[count][0] + ele_cent[i][0];
				gaupoi.y = gauss_po[count][1] + ele_cent[i][1];
				gaupoi.z = gauss_po[count][2] + ele_cent[i][2];

				//--------------------------------------------
				//循环单元的相关单元
				double min_dist = w_dist;
				Point_3D min_ep[2];		//最小距离对应线段的两个端点
				Point_3D min_fp;				//最小距离对应垂足
				for(int j=0; j<(int)relemd[i].size(); j++) //包括i单元自身
				{
					int ere = relemd[i][j];
					//保证了dist <= n*cell.delt_x(y,z)的情况能对单元中每个高斯点都取到完全的w_dist邻域
					double distx, disty, distz;
					distx = (int(w_dist/cell.delt_x-Zero)+1)*cell.delt_x;
					disty = (int(w_dist/cell.delt_y-Zero)+1)*cell.delt_y;
					distz = (int(w_dist/cell.delt_z-Zero)+1)*cell.delt_z;
					//判断那些由于周期边界条件产生的相关单元的情况
					Point_3D delt_dist(0,0,0);
					if(ele_cent[ere][0]-ele_cent[i][0] > cell.len_x-distx-Zero) delt_dist.x = ele_cent[ere][0] - cell.len_x;
					else if(ele_cent[ere][0]-ele_cent[i][0] < distx+Zero-cell.len_x) delt_dist.x = ele_cent[ere][0] + cell.len_x;

					if(ele_cent[ere][1]-ele_cent[i][1] > cell.wid_y-disty-Zero) delt_dist.y = ele_cent[ere][1] - cell.wid_y;
					else if(ele_cent[ere][1]-ele_cent[i][1] < disty+Zero-cell.wid_y) delt_dist.y = ele_cent[ere][1] + cell.wid_y;

					if(ele_cent[ere][2]-ele_cent[i][2] > cell.hei_z-distz-Zero) delt_dist.z = ele_cent[ere][2] - cell.hei_z;
					else if(ele_cent[ere][2]-ele_cent[i][2] < distz+Zero-cell.hei_z) delt_dist.z = ele_cent[ere][2] + cell.hei_z;

					//循环单元的相关纳米管线段
					for(int k=0; k<(int)elements[ere].relative_cnts.size(); k++)
					{
						int num = elements[ere].relative_cnts[k];
						//线段的两个端点
						Point_3D point[2] = { delt_dist+cnps[num], delt_dist+cnps[num+1] };

						//计算空间直线上一点的参数
						double ratio_t = ((gaupoi.x-point[0].x)*(point[1].x-point[0].x)+(gaupoi.y-point[0].y)*(point[1].y-point[0].y)+(gaupoi.z-point[0].z)*(point[1].z-point[0].z))				
													/((point[1].x-point[0].x)*(point[1].x-point[0].x)+(point[1].y-point[0].y)*(point[1].y-point[0].y)+(point[1].z-point[0].z)*(point[1].z-point[0].z));

						//计算直线上最近距离点（垂足）
						Point_3D footpoi(point[0].x+ratio_t*(point[1].x-point[0].x), point[0].y+ratio_t*(point[1].y-point[0].y), point[0].z+ratio_t*(point[1].z-point[0].z));

						if(ratio_t>0.0&&ratio_t<1.0)
						{
							//计算距离
							double temp_dist = sqrt((gaupoi.x-footpoi.x)*(gaupoi.x-footpoi.x)+(gaupoi.y-footpoi.y)*(gaupoi.y-footpoi.y)+(gaupoi.z-footpoi.z)*(gaupoi.z-footpoi.z));
							if(temp_dist<min_dist)
							{
								min_dist = temp_dist;
								min_ep[0] = point[0];
								min_ep[1] = point[1];
								min_fp = footpoi;
							}
						}
						else
						{
							for(int m=0; m<2; m++)
							{
								//计算距离
								double temp_dist = sqrt((gaupoi.x-point[m].x)*(gaupoi.x-point[m].x)+(gaupoi.y-point[m].y)*(gaupoi.y-point[m].y)+(gaupoi.z-point[m].z)*(gaupoi.z-point[m].z));
								if(temp_dist<min_dist)
								{
									min_dist = temp_dist;
									min_ep[0] = point[0];
									min_ep[1] = point[1];
									min_fp = footpoi;
								}
							}
						}

						//特殊情况高斯点和垂足是同一点, 说明高斯点在纳米管线上, 这是需要修改min_fp
						if(gaupoi.distance_to(footpoi)<Zero)
						{
							//在以point[1]-point[0]为法向量，以gaupoi为通过点的平面上随机选取一点，计算垂足点（即x'轴）
							//随机产生一个种子
							int seed = rand()%MAX_INT;
							seed = (2053*seed + 13849)%MAX_INT;

							Point_3D gpoi(0,0,0);
							gpoi.x = 2.0*seed/MAX_INT-1.0;  //产生[-1, 1]之间的一个数;
							seed = (2053*seed + 13849)%MAX_INT;
							gpoi.y = 2.0*seed/MAX_INT-1.0;  //产生[-1, 1]之间的一个数;
							if(fabs(point[1].z-point[0].z)>Zero) gpoi.z =  ((point[1].x-point[0].x)*(gpoi.x-gaupoi.x) + (point[1].y-point[0].y)*(gpoi.y-gaupoi.y))/(point[0].z-point[1].z) + gaupoi.z;
							else gpoi.z = 0.0;  //任意值都可以，因为系数(point[1].z-point[0].z)等于0

							min_fp = gaupoi + gaupoi - gpoi;
						}
					}
				}
				//计算权重值(前面只记录比w_dist小的, 所以不可能大于w_dist)
				if(min_dist<cnts.rad_max)	temp_cnt_volume += weight[count];  //记录高斯点的权重之和, 实际纳米管所占体积分数
				if(min_dist==w_dist)
				{
					temp_gwcnt[0].push_back(0.0);
					//保证旋转矩阵时单位矩阵
					for(int j=0; j<3; j++)
						for(int k=0; k<3; k++)
						{
							int tgn = j*3+k+1;
							if(j==k) temp_gwcnt[tgn].push_back(1.0);
							else temp_gwcnt[tgn].push_back(0.0);
						}
				}
				else if(min_dist<w_dist)
				{
					temp_interface_volume += weight[count];  //记录高斯点的权重之和, 实际界面层所占体积分数

					if(w_order==0) temp_gwcnt[0].push_back(w_constant);
					else if(w_order==1) temp_gwcnt[0].push_back(1.0-min_dist/w_dist);
					else if(w_order==3) temp_gwcnt[0].push_back(1.0+min_dist*min_dist*(2*min_dist-3*w_dist)/(w_dist*w_dist*w_dist));
					else if(w_order==-1) temp_gwcnt[0].push_back(2.0/(1.0+min_dist/w_dist)-1.0);
					else if(w_order==-2) temp_gwcnt[0].push_back(4.0/((1.0+min_dist/w_dist)*(1.0+min_dist/w_dist))-4.0/(1.0+min_dist/w_dist)+1.0);
					else	{	hout << "w_order=" << w_order << " 没有对应次幂的权重函数计算公式，请重新输入！" << endl; }

					//计算旋转矩阵值
					double tgve[3][4] = {{0}, {0}, {0}};
					
					//x'轴
					tgve[0][0] = gaupoi.x - min_fp.x;
					tgve[0][1] = gaupoi.y - min_fp.y;
					tgve[0][2] = gaupoi.z - min_fp.z;
					tgve[0][3] = sqrt(tgve[0][0]*tgve[0][0]+tgve[0][1]*tgve[0][1]+tgve[0][2]*tgve[0][2]);
					
					//z'轴
					tgve[2][0] = min_ep[1].x - min_ep[0].x;
					tgve[2][1] = min_ep[1].y - min_ep[0].y;
					tgve[2][2] = min_ep[1].z - min_ep[0].z;
					tgve[2][3] = sqrt(tgve[2][0]*tgve[2][0]+tgve[2][1]*tgve[2][1]+tgve[2][2]*tgve[2][2]);
					
					//y'轴 (z'叉乘x')
					tgve[1][0] = tgve[2][1]*tgve[0][2] - tgve[2][2]*tgve[0][1];
					tgve[1][1] = tgve[2][2]*tgve[0][0] - tgve[2][0]*tgve[0][2];
					tgve[1][2] = tgve[2][0]*tgve[0][1] - tgve[2][1]*tgve[0][0];
					tgve[1][3] = sqrt(tgve[1][0]*tgve[1][0]+tgve[1][1]*tgve[1][1]+tgve[1][2]*tgve[1][2]);
					
					//计算旋转矩阵
					for(int j=0; j<3; j++)
						for(int k=0; k<3; k++)
							temp_gwcnt[j*3+k+1].push_back(tgve[j][k]/tgve[j][3]);
				}
				else { hout << "注意！所得最小距离大于设定距离，请检查！" << endl; }
			}
		}

		#pragma omp critical
		{
			total_cnt_volume += temp_cnt_volume;
			total_interface_volume += temp_interface_volume;
			for(int i=0; i<(int)temp_num.size(); i++)
			{
				total_value_gps += temp_gwcnt[0][i];
				for(int j=0; j<10; j++)
				{
					gauss_wcnt[temp_num[i]][j] = temp_gwcnt[j][i];
				}
				//为检测10个权重值正确与否
				//if(gauss_wcnt[temp_num[i]][0]<-Zero||gauss_wcnt[temp_num[i]][0]>1.0+Zero) hout << "注意错误！gauss_wcnt[" << temp_num[i] << "][0]=" << gauss_wcnt[temp_num[i]][0] << "，请检查！" << endl;
				//for(int j=0; j<3; j++)
				//{
				//	double sum[2] = { 0 };
				//	for(int k=0; k<3; k++)	sum[0] +=  gauss_wcnt[temp_num[i]][3*j+k+1]*gauss_wcnt[temp_num[i]][3*j+k+1];
				//	for(int k=0; k<3; k++)	sum[1] +=  gauss_wcnt[temp_num[i]][j+3*k+1]*gauss_wcnt[temp_num[i]][j+3*k+1];
				//	if(fabs(sum[0]-1.0)>Zero||fabs(sum[1]-1.0)>Zero) hout << "注意错误！sum[0]=" << sum[0] << "sum[1]=" << sum[1] << "，请检查！" << endl;
				//}
			}
		}
	}

	hout << "    实际纳米管所占的体积分数为：" << total_cnt_volume*Jacobi << endl;
	hout << "    实际界面层所占的体积分数为：" << total_interface_volume*Jacobi << endl;
	hout << "    纳米管界面层的平均权重值：" << total_value_gps/((int)elements.size()*(int)gauss.size()) << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//决定单元的相关单元, 用于计算高斯点距纳米管的最小距离
void Global_Stiff_Matrix::Deter_relative_elements_min_dist(const vector<Node> &nodes, const vector<Element> &elements, const struct RVE_Geo &cell, const double &dist, vector<vector<int> > &relemd)const
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//单元的相关单元
	const int ES = (int)elements.size();
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
	{
		vector<int> temp_rele;
		for(int j=0; j<ES; j++)
		{
			if((fabs(centrele[j].x-centrele[i].x)<distx+Zero||fabs(centrele[j].x-centrele[i].x)>cell.len_x-distx-Zero)&&
				(fabs(centrele[j].y-centrele[i].y)<disty+Zero||fabs(centrele[j].y-centrele[i].y)>cell.wid_y-disty-Zero)&&
				(fabs(centrele[j].z-centrele[i].z)<distz+Zero||fabs(centrele[j].z-centrele[i].z)>cell.hei_z-distz-Zero))
			temp_rele.push_back(j);
		}
		relemd.push_back(temp_rele);
	}
}
//-----------------------------------------------------------------------------------------------
//计算单元每个高斯点的数据
void Global_Stiff_Matrix::Generate_element_gauss_data(const vector<Element> &elements, const vector<Node> &nodes, const string &com_mod, const vector<Node> &gauss, const vector<double> &weight,
																									double (*gauss_ns)[8], double (*gauss_nw)[8], double (*gauss_dfx)[3][8], double &Jacobi,  double (*gauss_po)[3], double (*ele_cent)[3])const
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//循环单元
	Point_3D scenter(0,0,0);	//标准立方体单元的中心点坐标
	int ie = 0; //设0号单元为标准立方体单元
	for(int count=0; count<(int)gauss.size(); count++)
	{	
		double Nshape[8] = {0};
		//--------------------------------------------
		//计算高斯积分点整体坐标
		Nshape[0]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
		Nshape[1]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
		Nshape[2]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
		Nshape[3]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
		Nshape[4]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
		Nshape[5]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
		Nshape[6]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);
		Nshape[7]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);

		//--------------------------------------------
		//计算高斯积分点坐标
		Point_3D gaupoi(0, 0, 0);		
		for(int j=0; j<8; j++) 
		{
			gaupoi.x += Nshape[j]*nodes[elements[ie].nodes_id[j]].x;
			gaupoi.y += Nshape[j]*nodes[elements[ie].nodes_id[j]].y;
			gaupoi.z += Nshape[j]*nodes[elements[ie].nodes_id[j]].z;
		}

		//记录高斯积分点坐标
		gauss_po[count][0] = gaupoi.x;
		gauss_po[count][1] = gaupoi.y;
		gauss_po[count][2] = gaupoi.z;

		//--------------------------------------------
		//计算Ｊ矩阵
		//--------------------------------------------
		//形函数N对gauss[count].x, gauss[count].y, gauss[count].z的偏导矩阵
		double diff[3][8];
		diff[0][0]=-0.125*(1.0-gauss[count].y)*(1.0-gauss[count].z);
		diff[0][1]=-diff[0][0];                         
		diff[0][2]=0.125*(1.0+gauss[count].y)*(1.0-gauss[count].z);
		diff[0][3]=-diff[0][2];
		diff[0][4]=-0.125*(1.0-gauss[count].y)*(1.0+gauss[count].z);
		diff[0][5]=-diff[0][4];
		diff[0][6]=0.125*(1.0+gauss[count].y)*(1.0+gauss[count].z);
		diff[0][7]=-diff[0][6];

		diff[1][0]=-0.125*(1.0-gauss[count].x)*(1.0-gauss[count].z);
		diff[1][1]=-0.125*(1.0+gauss[count].x)*(1.0-gauss[count].z);
		diff[1][2]=-diff[1][1];
		diff[1][3]=-diff[1][0];
		diff[1][4]=-0.125*(1.0-gauss[count].x)*(1.0+gauss[count].z);
		diff[1][5]=-0.125*(1.0+gauss[count].x)*(1.0+gauss[count].z);
		diff[1][6]=-diff[1][5];
		diff[1][7]=-diff[1][4];

		diff[2][0]=-0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y);
		diff[2][1]=-0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y);
		diff[2][2]=-0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y);
		diff[2][3]=-0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y);
		diff[2][4]=-diff[2][0];
		diff[2][5]=-diff[2][1];
		diff[2][6]=-diff[2][2];
		diff[2][7]=-diff[2][3];

		//--------------------------------------------------
		//单元节点坐标矩阵
		double elenode[8][3];
		for(int j=0; j<8; j++)
		{
			elenode[j][0]=nodes[elements[ie].nodes_id[j]].x;
			elenode[j][1]=nodes[elements[ie].nodes_id[j]].y;
			elenode[j][2]=nodes[elements[ie].nodes_id[j]].z;
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
		Jacobi = Jmatrix[0][0]*(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])
						-Jmatrix[0][1]*(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])
						+Jmatrix[0][2]*(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0]);

		if(com_mod=="hybrid"||com_mod=="nonlocal") 
		{
			for(int j=0; j<8; j++) 
			{
				gauss_ns[count][j] = Nshape[j]*Jacobi; 
				gauss_nw[count][j] = gauss_ns[count][j]*weight[count];
			} 	
		}

		if(com_mod=="hybrid"||com_mod=="fem") 
		{
			//----------------------------------------------------
			//求出J矩阵的逆矩阵
			double Jinverse[3][3];
			
			Jinverse[0][0]=(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])/Jacobi;
			Jinverse[1][1]=(Jmatrix[0][0]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][0])/Jacobi;
			Jinverse[2][2]=(Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0])/Jacobi;

			Jinverse[0][1]=-(Jmatrix[0][1]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][1])/Jacobi;
			Jinverse[0][2]=(Jmatrix[0][1]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][1])/Jacobi;

			Jinverse[1][0]=-(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])/Jacobi;
			Jinverse[1][2]=-(Jmatrix[0][0]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][0])/Jacobi;

			Jinverse[2][0]=(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0])/Jacobi;
			Jinverse[2][1]=-(Jmatrix[0][0]*Jmatrix[2][1]-Jmatrix[0][1]*Jmatrix[2][0])/Jacobi;

			//-------------------------------------------------------
			//求出N对x,y,z的偏导
			double diffxy[3][8];
			for(int j=0; j<3; j++)
				for(int k=0; k<8; k++)
				{
					diffxy[j][k]=0;
					for(int m=0; m<3; m++)
						diffxy[j][k] += Jinverse[j][m]*diff[m][k];
				
					//记录J矩阵的行列式值
					gauss_dfx[count][j][k] = diffxy[j][k];
				}
		}
	}

	//计算该单元的中心点坐标
	for(int j=0; j<8; j++) 
	{
		scenter.x += nodes[elements[ie].nodes_id[j]].x;
		scenter.y += nodes[elements[ie].nodes_id[j]].y;
		scenter.z += nodes[elements[ie].nodes_id[j]].z;
	}
	scenter = scenter/8;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//循环单元, 计算各个单元中心点的相对坐标（相对于0号单元）
	for(int i=0; i<(int)elements.size(); i++)
	{
		for(int j=0; j<3; j++) ele_cent[i][j] = 0.0;
		//计算该单元的中心点坐标
		for(int j=0; j<8; j++) 
		{
			ele_cent[i][0] += nodes[elements[i].nodes_id[j]].x;
			ele_cent[i][1] += nodes[elements[i].nodes_id[j]].y;
			ele_cent[i][2] += nodes[elements[i].nodes_id[j]].z;
		}
		for(int j=0; j<3; j++) ele_cent[i][j] = ele_cent[i][j]/8;
		ele_cent[i][0] = ele_cent[i][0] - scenter.x;
		ele_cent[i][1] = ele_cent[i][1] - scenter.y;
		ele_cent[i][2] = ele_cent[i][2] - scenter.z;
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------
//计算等效应变能
int Global_Stiff_Matrix::Calculate_equivalent_energy_simple(const vector<long> &Iz, const vector<int> &Ig, const vector<double> &total_matrix, const vector<double>* nl_equright, const double backup_ege[],
																											const vector<vector<double> > &disp_solu,  double equivalent_energy[])const
{
	//------------------------------------------------------------------------
	//向量初始化
	for(int i=0; i<9; i++) equivalent_energy[i] = 0;

	clock_t ct0,ct1;
	//------------------------------------------------------------------------
	//循环所有节点的相关节点
	ct0 = clock();
	hout << "-_- 开始计算等效应变能" << endl;

	for(int i=0; i<(int)Iz.size(); i++)
	{
		//-----------------------------------------------------------
		//相关点
		if(i!=0)
		{
			for(long j=Iz[i-1]; j<Iz[i]; j++)
			{
				const long Ij = 6*(long)i + 9*j;
				for(int k=0; k<(int)disp_solu.size(); k++)
					for(int m=0; m<=2; m++)
						for(int n=0; n<=2; n++)
						{
							equivalent_energy[k] += disp_solu[k][3*i+m]*total_matrix[Ij+3*m+n]*disp_solu[k][3*Ig[j]+n];   //本身
							equivalent_energy[k] += disp_solu[k][3*Ig[j]+m]*total_matrix[Ij+m+3*n]*disp_solu[k][3*i+n];	 //转置
						}
			}
		}
		//-----------------------------------------------------------
		//对角阵
		const long II = 6*(long)i + 9*Iz[i];
		for(int j=0; j<(int)disp_solu.size(); j++)
			for(int k=0; k<=2; k++)
				for(int m=0; m<=2; m++)
				{
					if(k==2||m==2)
						equivalent_energy[j] += disp_solu[j][3*i+k]*total_matrix[II+k+m+1]*disp_solu[j][3*i+m];
					else
						equivalent_energy[j] += disp_solu[j][3*i+k]*total_matrix[II+k+m]*disp_solu[j][3*i+m];
				}
	}

	//------------------------------------------------------------------------
	//右端项周期边界条件部分
	for(int i=0; i<(int)disp_solu.size(); i++)
	{
		double sum = 0.0;
		for(int j=0; j<(int)disp_solu[i].size(); j++) 
			sum += disp_solu[i][j]*nl_equright[i][j];
		equivalent_energy[i] -= 2*sum;
		equivalent_energy[i] += 0.5*backup_ege[i];
	}

	//---------------------------------------------------------------------------
	//hout << "Equivalent Energy(simple) Vector:" << endl;
	//for(int i=0; i<9; i++)	hout << equivalent_energy[i] << " ";
	//hout << endl;

	//---------------------------------------------------------------------------
	ct1 = clock();
	hout << "    计算等效应变能耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
	hout << "^_^ 计算等效应变能操作完毕！" << endl << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------
//计算等效应变能
int Global_Stiff_Matrix::Calculate_equivalent_energy(const vector<Element> &elements, const vector<MatPro> &mats, const struct Decay_Para &decay, const vector<vector<double> > &disp_solu, const struct RVE_Geo &cell,
																							 const string &com_mod,  double equivalent_energy[])const
{
	//------------------------------------------------------------------------
	//向量初始化
	for(int i=0; i<9; i++) equivalent_energy[i] = 0;

	clock_t ct0,ct1;
	//------------------------------------------------------------------------
	//数据读入二进制文件, 用于计算等效应变能时用
	fstream fin( "GaussDats", ios::in|ios::binary);
	int GS, ES; 	//变量定义	
	fin.read((char *)&GS, sizeof(int));	
	fin.read((char *)&ES, sizeof(int));
	double *weight = new double[GS];									//记录单元高斯点的权重值
	fin.read((char *)weight, sizeof(double)*GS);

	double (*gauss_po)[3] = new double [GS][3];				//标准立方体单元的高斯点坐标
	double (*gauss_ns)[8] = new double [GS][8];				//记录单元高斯点的形函数
	double (*ele_cent)[3] = new double [ES][3];					//记录单元中心点位置向量(x,y,z)分别放于[0],[1]和[2]中（以标准立方体单元中心点为原点）
	double (*gauss_wcnt)[10] = new double[ES*GS][10];	//纳米管的权重
	double Jacobi;																	//标准立方体单元的雅可比值
	double (*gauss_dfx)[3][8] = new double [GS][3][8];		//记录单元高斯点形函数的导数

	if(com_mod=="hybrid"||com_mod=="nonlocal")
	{
		fin.read((char *)gauss_po, sizeof(double)*GS*3);
		fin.read((char *)gauss_ns, sizeof(double)*GS*8);	
		fin.read((char *)ele_cent, sizeof(double)*ES*3);
		fin.read((char *)gauss_wcnt, sizeof(double)*ES*GS*10);
	}
	if(com_mod=="hybrid"||com_mod=="fem")
	{
		fin.read((char *)&Jacobi, sizeof(double));
		fin.read((char *)gauss_dfx, sizeof(double)*GS*3*8);
	}
	fin.close();

	//------------------------------------------------------------------------
	//循环所有单元
	ct0 = clock();
	hout << "-_- 开始计算等效应变能" << endl;
	//执行openmp
	#pragma omp parallel
	{
		//定义单刚
		double elements_energy[9] = {0};
		double (*gwc_left)[10] = new double [GS][10];
		double (*gwc_right)[10] = new double [GS][10];

		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++)
		{
			//--------------------------------------------------
			//单元顶点的位移值
			double disp_left[9][24];
			for(int j=0; j<9; j++)
				for(int k=0; k<8; k++)
				{
					disp_left[j][3*k] = disp_solu[j][3*elements[i].nodes_id[k]];
					disp_left[j][3*k+1] = disp_solu[j][3*elements[i].nodes_id[k]+1];
					disp_left[j][3*k+2] = disp_solu[j][3*elements[i].nodes_id[k]+2];
				}

			//---------------------------------------------------------------------------
			//计算RVE背景网格单元(长程力作用)的等效应变能
			if(com_mod=="hybrid"||com_mod=="nonlocal")
			{
				//开始计时
				clock_t ctn1,ctn2;
				ctn1 = clock();

				//--------------------------------------------------
				//初始化变量
				double elec_left[3] = { ele_cent[i][0], ele_cent[i][1], ele_cent[i][2] };

				//--------------------------------------------------
				//单元的每个高斯点关于最近距离纳米管的权重值
				const int IGS = i*GS;
				for(int j=0; j<GS; j++)
					for(int k=0; k<10; k++)
						gwc_left[j][k] = gauss_wcnt[IGS+j][k];

				//------------------------------------------------------------------------
				//长程力作用
				for(int j=0; j<(int)elements[i].relative_eles.size(); j++)
				{
					//--------------------------------------------------
					//初始化变量
					const int ere = elements[i].relative_eles[j];
					double elec_right[3] = { ele_cent[ere][0], ele_cent[ere][1], ele_cent[ere][2] };

					//--------------------------------------------------
					//单元顶点的位移值
					double disp_right[9][24];
					for(int k=0; k<9; k++)
						for(int m=0; m<8; m++)
						{
							disp_right[k][3*m] = disp_solu[k][3*elements[ere].nodes_id[m]];
							disp_right[k][3*m+1] = disp_solu[k][3*elements[ere].nodes_id[m]+1];
							disp_right[k][3*m+2] = disp_solu[k][3*elements[ere].nodes_id[m]+2];
						}
					//--------------------------------------------------
					//判断那些由于周期边界条件产生的相关单元的情况，并修改相关值
					Verify_Periodical_Condition_Modify_Values(elec_left, elec_right, disp_right, cell, decay);

					//--------------------------------------------------
					//单元的每个高斯点关于最近距离纳米管的权重值
					const int EREGS = ere*GS;
					for(int k=0; k<GS; k++)
						for(int m=0; m<10; m++)
							gwc_right[k][m] = gauss_wcnt[EREGS+k][m];

					//------------------------------------------------------------------------
					//计算RVE背景网格单元(长程力作用)等效应变能
					Calculate_Longforce_Energy(elements_energy, weight, gauss_ns, decay, gauss_po, elec_left, elec_right, gwc_left, gwc_right, disp_left, disp_right, GS);
				}
				//------------------------------------------------------------------------
				ctn2 = clock();
				hout << "Equivalent Energy: total num of elements: " << (int)elements.size() << "; Element " << i << " took time: " << (double)(ctn2-ctn1)/CLOCKS_PER_SEC << "sec;" << endl;
			}

			//---------------------------------------------------------------------------
			//计算RVE背景网格单元接触力等效应变能
			if(com_mod=="hybrid"||com_mod=="fem")
			{
				Calculate_Contactforce_Energy(elements_energy, disp_left, mats, Jacobi, gauss_dfx, weight, GS);
			}
		}

		//添加到总应变能
		#pragma omp critical
		{
			for(int i=0; i<9; i++)	equivalent_energy[i] += elements_energy[i];
		}

		//---------------------------------------------------------------------------
		//删除指针
		delete[] gwc_left;
		delete[] gwc_right;
	}

	//---------------------------------------------------------------------------
	//删除指针
	delete[] weight;
	delete[] gauss_ns;
	delete[] gauss_po;
	delete[] ele_cent;
	delete[] gauss_wcnt;
	delete[] gauss_dfx;

	//---------------------------------------------------------------------------
	hout << "Equivalent Energy Vector:" << endl;
	for(int i=0; i<9; i++)	hout << equivalent_energy[i] << " ";
	hout << endl;

	//---------------------------------------------------------------------------
	ct1 = clock();
	hout << "    计算等效应变能耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
	hout << "^_^ 计算等效应变能操作完毕！" << endl << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------
//判断那些由于周期边界条件产生的相关单元的情况，并修改相关值
void Global_Stiff_Matrix::Verify_Periodical_Condition_Modify_Values(const double elec_left[], double elec_right[], double (*disp_right)[24], const struct RVE_Geo &cell, const struct Decay_Para &decay)const
{
	//---------------------------------------------------------------------------
	//判断那些由于周期边界条件产生的相关单元
	bool mark = false;
	double peri_dist[3] = {0};
	if(elec_right[0]-elec_left[0] > cell.len_x-decay.R-Zero) { elec_right[0] -= cell.len_x; peri_dist[0] = -cell.len_x; mark = true; }
	else if(elec_right[0]-elec_left[0] < -cell.len_x+decay.R+Zero) { elec_right[0] += cell.len_x; peri_dist[0] = cell.len_x; mark = true; }

	if(elec_right[1]-elec_left[1] > cell.wid_y-decay.R-Zero) { elec_right[1] -= cell.wid_y; peri_dist[1] = -cell.wid_y; mark = true; }
	else if(elec_right[1]-elec_left[1] < -cell.wid_y+decay.R+Zero) { elec_right[1] += cell.wid_y; peri_dist[1] = cell.wid_y; mark = true; }

	if(elec_right[2]-elec_left[2] > cell.hei_z-decay.R-Zero) { elec_right[2] -= cell.hei_z; peri_dist[2] = -cell.hei_z; mark = true; }
	else if(elec_right[2]-elec_left[2] < -cell.hei_z+decay.R+Zero) { elec_right[2] += cell.hei_z; peri_dist[2] = cell.hei_z; mark = true; }
	
	//---------------------------------------------------------------------------
	if(mark)
	{
		for(int i=0; i<9; i++)
		{
			//---------------------------------
			double E[3][3] = {{0}, {0}, {0}}; //用于周期边界条件处理时的变形矩阵条件
			switch(i)	//设置均匀化应变做为周期边界条件
			{
			case 0: E[0][0]=0.1; break;
			case 1: E[1][1]=0.1; break;
			case 2: E[2][2]=0.1; break;
			case 3: E[0][1]=0.05; E[1][0]=0.05; break;
			case 4: E[1][2]=0.05; E[2][1]=0.05; break;
			case 5: E[0][2]=0.05; E[2][0]=0.05; break;
			case 6: E[0][0]=0.1; E[1][1]=0.1; break;
			case 7: E[1][1]=0.1; E[2][2]=0.1; break;
			case 8: E[0][0]=0.1; E[2][2]=0.1; break;
			default: hout << "错误！ 处理周期边界条件约束时，循环次数值等于" << i << "小于0或者大于8！请检查！" << endl;
			}

			//---------------------------------
			//均匀化位移差
			double uni_disp[3] = {0};
			for(int j=0; j<3; j++)
				for(int k=0; k<3; k++) 
					uni_disp[j] += E[j][k]*peri_dist[k];

			//---------------------------------
			//计算相关位移值
			for(int j=0; j<8; j++)
				for(int k=0; k<3; k++)
					disp_right[i][3*j+k] += uni_disp[k];
		}
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------
//计算RVE背景网格单元(长程力作用)等效应变能
void Global_Stiff_Matrix::Calculate_Longforce_Energy(double elements_energy[], const double weight[], const double (*gauss_ns)[8], const struct Decay_Para &decay, const double (*gauss_po)[3], const double elec_left[],
																								  const double elec_right[], const double (*gwc_left)[10], const double (*gwc_right)[10], const double (*disp_left)[24], const double (*disp_right)[24], const int &GS)const
{
	//--------------------------------------------------	
	//循环本单元高斯点计算积分
	for(int count1=0; count1<GS; count1++)
	{
		//记录每个本单元高斯点上的等效应变能相关值
		double ele_energy[9] = { 0 };

		//--------------------------------------------------	
		double u_left[9][3] = { {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0} };  //初始化
		for(int i=0; i<9; i++) 
			for(int j=0; j<3; j++)
				for(int k=0; k<8; k++)
					u_left[i][j] += gauss_ns[count1][k]*disp_left[i][3*k+j];

		//--------------------------------------------------	
		//左端高斯点坐标
		Point_3D gaupoi_left(0, 0, 0);
		gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
		gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
		gaupoi_left.z = gauss_po[count1][2] + elec_left[2];

		//------------------------------------------------------------------------------------------------------------------------
		//循环外单元高斯点计算积分
		for(int count2=0; count2<GS; count2++)
		{
			//--------------------------------------------------	
			double u_right[9][3] = { {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0} };  //初始化
			for(int i=0; i<9; i++) 
				for(int j=0; j<3; j++)
					for(int k=0; k<8; k++)
						u_right[i][j] += gauss_ns[count2][k]*disp_right[i][3*k+j];

			double u_diff[9][3];
			for(int i=0; i<9; i++)
				for(int j=0; j<3; j++)
					u_diff[i][j] = u_right[i][j] - u_left[i][j];

			//--------------------------------------------
			//右端高斯点坐标
			Point_3D gaupoi_right(0, 0, 0);
			gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
			gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
			gaupoi_right.z = gauss_po[count2][2] + elec_right[2];

			//--------------------------------------------
			//权重值判断
			if(gwc_left[count1][0]<Zero&&gwc_right[count2][0]<Zero) continue;  //等于零值，没有长程效果

			//--------------------------------------------
			//计算长程作用衰减函数值
			const double xleft = gaupoi_right.x-gaupoi_left.x;
			const double yleft = gaupoi_right.y-gaupoi_left.y;
			const double zleft = gaupoi_right.z-gaupoi_left.z;

			//--------------------------------------------
			//计算bond的长度
			const double dis_squr = xleft*xleft + yleft*yleft + zleft*zleft;
			const double poi_dis = sqrt(dis_squr);		//两点之间的距离
			
			if(poi_dis>decay.R+Zero||poi_dis<Zero) continue;		//圆形积分域

			//--------------------------------------------
			//计算left点权重
			double sum_left = 0;
			if(gwc_left[count1][0]>Zero)
			{
				const double x = gwc_left[count1][1]*xleft+ gwc_left[count1][2]*yleft + gwc_left[count1][3]*zleft;
				const double y = gwc_left[count1][4]*xleft + gwc_left[count1][5]*yleft + gwc_left[count1][6]*zleft;
				const double z = gwc_left[count1][7]*xleft + gwc_left[count1][8]*yleft + gwc_left[count1][9]*zleft;

				const double cos2sita = z*z/dis_squr;		//cos(sita)^2
				double cos2pha;
				if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
				else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

				sum_left = decay.acoe[0] + decay.acoe[1]*0.5*(3*cos2sita-1) + decay.acoe[2]*(2*cos2pha-1)*3*(1-cos2sita) + decay.acoe[3]*0.125*(35*cos2sita*cos2sita - 30*cos2sita + 3)
									+ decay.acoe[4]*(2*cos2pha-1)*7.5*(7*cos2sita-1)*(1-cos2sita) + decay.acoe[5]*(8*cos2pha*cos2pha-8*cos2pha+1)*105*(1-cos2sita)*(1-cos2sita);
			}

			//--------------------------------------------
			//计算right点权重
			double sum_right = 0;
			if(gwc_right[count2][0]>Zero)
			{
				const double xright = gaupoi_left.x-gaupoi_right.x;
				const double yright = gaupoi_left.y-gaupoi_right.y;
				const double zright = gaupoi_left.z-gaupoi_right.z;

				const double x = gwc_right[count2][1]*xright + gwc_right[count2][2]*yright + gwc_right[count2][3]*zright;
				const double y = gwc_right[count2][4]*xright + gwc_right[count2][5]*yright + gwc_right[count2][6]*zright;
				const double z = gwc_right[count2][7]*xright + gwc_right[count2][8]*yright + gwc_right[count2][9]*zright;

				const double cos2sita = z*z/dis_squr;		//cos(sita)^2
				double cos2pha;
				if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
				else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

				sum_right = decay.acoe[0] + decay.acoe[1]*0.5*(3*cos2sita-1) + decay.acoe[2]*(2*cos2pha-1)*3*(1-cos2sita) + decay.acoe[3]*0.125*(35*cos2sita*cos2sita - 30*cos2sita + 3)
									+ decay.acoe[4]*(2*cos2pha-1)*7.5*(7*cos2sita-1)*(1-cos2sita) + decay.acoe[5]*(8*cos2pha*cos2pha-8*cos2pha+1)*105*(1-cos2sita)*(1-cos2sita);
			}

			//--------------------------------------------
			//计算权重函数值
			const double sum = 0.5*(gwc_left[count1][0]*sum_left+gwc_right[count2][0]*sum_right);
			if(fabs(sum)<Zero) continue;  //等于零值，没有长程效果

			//--------------------------------------------
			//计算长程力矩阵					
			const double gv = exp(-poi_dis/decay.radius)*sum*weight[count2]; //衰减以及权重函数值

			double Gmatrix[3][3];
			Gmatrix[0][0] = gv*xleft*xleft;
			Gmatrix[1][1] = gv*yleft*yleft;
			Gmatrix[2][2] = gv*zleft*zleft;
			Gmatrix[0][1] = gv*xleft*yleft;
			Gmatrix[0][2] = gv*xleft*zleft;
			Gmatrix[1][2] = gv*yleft*zleft;
			Gmatrix[1][0] = Gmatrix[0][1];
			Gmatrix[2][0] = Gmatrix[0][2];
			Gmatrix[2][1] = Gmatrix[1][2];

			for(int i=0; i<9; i++)
			{
				double temp_val[3] = {0};
				for(int j=0; j<3; j++)
				{
					for(int k=0; k<3; k++)	temp_val[j] += Gmatrix[j][k]*u_diff[i][k];
					ele_energy[i] += u_diff[i][j]*temp_val[j];
				}
			}
		}
		for(int i=0; i<9; i++) elements_energy[i] +=  0.5*ele_energy[i]*weight[count1];
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------
//计算RVE背景网格单元(接触力作用)等效应变能
void Global_Stiff_Matrix::Calculate_Contactforce_Energy(double elements_energy[], const double (*disp)[24], const vector<MatPro> &mats, const double &Jacobi, const double (*gauss_dfx)[3][8], const double weight[], const int &GS)const
{
	//--------------------------------------------------
	//初始化
	double element_stiff_matrix[24][24];
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++) 
			element_stiff_matrix[i][j] = 0;
	//--------------------------------------------------
	//循环高斯点计算积分
	for(int count=0; count<GS; count++)
	{
		//-----------------------------------------------------------------
		//计算此单元所对应的材料弹性矩阵
		double ele_elas[6][6];		
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				ele_elas[i][j] = mats[0].elas_matrix[i][j];
		
		//--------------------------------------------------------
		//B矩阵
		double B[6][24] = {{0}, {0}, {0}, {0}, {0}, {0}};

		for(int i=0; i<8; i++)
		{
			B[0][i*3+0]=gauss_dfx[count][0][i];
			B[1][i*3+1]=gauss_dfx[count][1][i];
			B[2][i*3+2]=gauss_dfx[count][2][i];
			B[3][i*3+0]=gauss_dfx[count][1][i];
			B[3][i*3+1]=gauss_dfx[count][0][i];
			B[4][i*3+1]=gauss_dfx[count][2][i];
			B[4][i*3+2]=gauss_dfx[count][1][i];
			B[5][i*3+0]=gauss_dfx[count][2][i];
			B[5][i*3+2]=gauss_dfx[count][0][i];
		}

		//--------------------------------------------------------
		//计算B_trans
		double B_trans[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
				B_trans[i][j] = B[j][i];
		//--------------------------------------------------------------------------------------------
		//求出B_trans矩阵与ele_elas矩阵的乘积array1
		double array1[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
			{
				array1[i][j]=0; 
				for(int k=0; k<6; k++)
					array1[i][j] += B_trans[i][k]*ele_elas[k][j];
			}
		//求出array1矩阵与B矩阵的乘积array2
		double array2[24][24];
		for(int i=0; i<24; i++)
			for(int j=0; j<24; j++)
			{
				array2[i][j]=0;
				for(int k=0; k<6; k++)
					array2[i][j] += array1[i][k]*B[k][j];
				element_stiff_matrix[i][j] += array2[i][j]*weight[count];
			}
	}

	//在高斯点循环外乘雅可比值减少乘法次数
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++) 
			element_stiff_matrix[i][j] = element_stiff_matrix[i][j]*Jacobi;

	//--------------------------------------------
	//组装应变能
	for(int i=0; i<9; i++)
		for(int j=0; j<24; j++)
			for(int k=0; k<24; k++)
				elements_energy[i] += disp[i][j]*element_stiff_matrix[j][k]*disp[i][k];
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//读入一行信息，并跳过注释行（以"%"开头）；
string Global_Stiff_Matrix::Get_Line(ifstream &infile)const
{
	string s;
	//读入信息一行
	getline(infile,s);
	//跳过注释行     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===============================================================
