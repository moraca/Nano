//===========================================================================
// Algorithm_FEM.cpp
// �㷨������Ԫʵ��ģ�����Ա����
// Member functions in classes to implement algorithm and finite element method
//===========================================================================
#include "Algorithm_FEM.h"

//--------------------------------------------------------------------------
//�㷨������Ԫ������
int Algorithm_FEM::Solve(ifstream &infile, vector<Node> &nodes, const vector<int>* peri_bnods, vector<Element> &elements,  const vector<MatPro> &mats, const struct Decay_Para &decay, const struct RVE_Geo &cell, 
											   const struct CNT_Geo &cnts, const vector<Point_3D> &cnps)
{
	clock_t ct0,ct1;  //���������¼ִ�еĿ�ʼ�ͽ���ʱ��

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���������������
	if(Import_computational_data(infile)==0) return 0;
	backup_file_name = "Displacement_Solution_Backup";  //����λ�ƽ⼰��ЧӦ���������ļ���
	if(wr_mod=="read_all")  return 1;   //����ȫ�����ݰ���λ�ƽ�͵�ЧӦ����
	
	Glosmat = new Global_Stiff_Matrix;		//��������նȾ����࣬ÿ�����������ɶȣ�����wr_mod��дģʽ��Ե��Ҫ��ǰ�ڴ˶���

	ct0 = clock();
	hout << "-_- ��ʼͳ�Ƶ�Ԫ����ص�Ԫ�ͽڵ����ؽڵ�" << endl;
	if(Deter_relative_nodes_elements(nodes, elements, decay.R, cell, com_mod)==0) return 0; //������Ԫ����ص�Ԫ�ͽڵ����ؽڵ�, ��Ԫ����ص�Ԫ��ϢҪ�ڹ��Ƶ�ЧӦ����ʱ�õ�

	ct1 = clock();
	hout << "    ͳ�Ƶ�Ԫ����ص�Ԫ�ͽڵ����ؽڵ��ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
	hout << "^_^ ͳ�Ƶ�Ԫ����ص�Ԫ�ͽڵ����ؽڵ����!" << endl << endl;

	//---------------------------------------------------------------------------
	//����洢����, ���ڵ�ЧӦ���ܵļ���
	vector<long> backup_iz;
	vector<int> backup_ig;
	vector<double> backup_matrix;
	vector<double> backup_equright[9];
	double backup_ege[9] = { 0 };

	if(wr_mod=="write")
	{
		//-----------------------------------------------------------------------------------------------------------------------------------------
		//���ɽ����洢�ڵ��ϵ����
		ct0 = clock();
		hout << "-_- ��ʼ�����洢�նȾ���" << endl;

		Solv = new SolveEqu;	//����ⷽ������
		vector<long> Iz;		//Iz,Ig�Ǳ����洢�������Ϣ
		vector<int> Ig;			//Iz,Ig�Ǳ����洢�������Ϣ
		if(Solv->izig(nodes, Iz, Ig)==0) return 0;
		//ע�ͣ��ڵ����ؽڵ���Ϣ�͵�Ԫ��������ص�Ԫ��ϢҪ�ں����õ�, �ڴ˲�����

		ct1 = clock();
		hout << "    �����洢�նȾ����ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
		hout << "^_^ �����洢�նȾ������!" << endl << endl;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//��������նȾ���
		ct0 = clock();
		hout << "-_- ��ʼ�����ܸ���" << endl;

		//�����ܸնȾ���洢�ռ�(���յ����ܸյķ�ʽ��������߶�)
		const long Total_Matrix_Size = 6*(long)nodes.size() + 9*Iz.back();
		vector<double> total_matrix(Total_Matrix_Size, 0);
		if(Total_Matrix_Size>(long)total_matrix.max_size()) { hout << "ע�⣡����ĸնȾ���ߴ��������vector������������ֵ��" << endl; return 0; }
		hout << "total_matrix.size: " << 	total_matrix.size() << "  " << "Iz.back: " << Iz.back() << endl;

		vector<double> nl_equright[9];		//���ڼ�¼���̷Ǿֲ����ڱ߽�㴦�������ڱ߽��������õ��Ҷ���
		for(int i=0; i<9; i++) nl_equright[i].assign(3*nodes.size(), 0.0);	//��ʼ��

		if(Glosmat->Gen_global_stiff_matrix(infile, cnps, nodes, elements, mats, decay, cell, cnts, com_mod, Iz, Ig, total_matrix, nl_equright, backup_ege)==0) return 0; //��������նȾ���
		//�����ڼ����ЧӦ���ܵ�ʱ����Ҫ�õ�����࣬��������û��ɾ��

		//------------------------------------------------------------------------
		//���ݶ������д��������ļ����Ա��Ժ󵥶�������Է�����ʱ�õ�
		string maeq_mod = "write";
//		string maeq_mod = "read";
		if(maeq_mod=="read") { Get_Line(infile); Get_Line(infile); }
		if(Write_read_matrix_equright_data(Total_Matrix_Size, total_matrix, nl_equright, backup_ege, maeq_mod)==0) return 0;

		//---------------------------------------------------------------------------
		backup_iz = Iz;									//��¼ԭʼIz
		backup_ig = Ig;									//��¼ԭʼIg
		backup_matrix = total_matrix;			//��¼ԭʼ�նȾ���
		for(int i=0; i<9; i++) backup_equright[i] = nl_equright[i];  //��¼ԭʼ�Ҷ�����
		//---------------------------------------------------------------------------

		ct1 = clock();
		hout << "    �����ܸ����ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
		hout << "^_^ �����ܸ������!" << endl << endl;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//����λ��Լ��
		ct0 = clock();
		hout << "-_- ��ʼ����λ��Լ��" << endl;

		int bnod_num = 0;		//���ڼ�¼�߽�ڵ�ĸ�����ע��Ҫ��ʼ����ֵ
		vector<int> ip;			//���ڼ�¼�߽�ڵ�������Ϣ
		vector<double> vp;	//���ڼ�¼�߽�ڵ��ϵ�ֵ
		if(Solv->Fixed_displacement_constraints(infile, nodes, bnod_num, ip, vp)==0) return 0;		//���ù̶�λ��Լ��
		ct1 = clock();
		hout << "    λ��Լ�������ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
		hout << "^_^ λ��Լ���������!" << endl << endl;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//�������ڱ߽�����
		ct0 = clock();
		hout << "-_- ��ʼ�������ڱ߽�����" << endl;

	//	Solv->Complete_matrix_equright_testing_periodical_constraints(peri_bnods, total_matrix, nl_equright, nodes, Iz, Ig);   //����ȫ�ĸնȾ����Ҷ��������ڱ߽������Ĵ���
		Solv->Periodical_boundary_constraints(nodes, peri_bnods, Iz, Ig, total_matrix, nl_equright);
		for(int i=0; i<(int)nodes.size(); i++) nodes[i].relative_nods.clear();  //�ڵ����ؽڵ���Ϣ�ڴ�����
	//	Solv->Export_complete_matrix_equright(total_matrix, nl_equright, nodes, Iz, Ig);		//�����ȫ�ĸնȾ����Ҷ������ڼ��
		ct1 = clock();
		hout << "    ���ڱ߽����������ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
		hout << "^_^ ���ڱ߽������������!" << endl << endl;

		//---------------------------------------------------------------------------
		//������ά�̶�λ��Լ�������λ��ֵ����
		ct0 = clock();
		hout << "-_- ��ʼ��ά�̶�λ��Լ�������λ��ֵ����" << endl;

		Solv->Deal_with_displacement_zero_value(bnod_num, (int)nodes.size(), Iz, Ig, ip, vp, nl_equright, total_matrix);
		ct1 = clock();
		hout << "    ��ά�̶�λ��Լ�������λ��ֵ���������ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
		hout << "^_^ ��ά�̶�λ��Լ�������λ��ֵ�����������!" << endl << endl;

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//ѭ������������Է�����
		ct0 = clock();
		hout << "-_- ��ʼѭ������������Է�����" << endl;

		vector<double> temp_solu;			//��ʱ����
		U_Solution.assign(9, temp_solu);	//��ʼ��
		//ִ��openmp
		#pragma omp parallel
		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<9; i++)
		{
			clock_t ct2, ct3;	//��ʱ��
			ct2 = clock();
			vector<double> solution(3*nodes.size(), 0.0); //����
			Solv->Solve_linear_equations(bnod_num, (int)nodes.size(), Iz, Ig, ip, vp, total_matrix, nl_equright[i], solution);	

			//-----------------------------------------------------------------------------------------------------------------------------------------
			//����������ڱ߽��Ͻڵ�Ľ�����ڱȶ�
//			Solv->Compared_periodical_bounday_solutions(i, peri_bnods, solution);
			//---------------------------------------------------------------------------
			//��¼�����
			#pragma omp critical
			{
				U_Solution[i] = solution;
			}
			ct3 = clock();
			hout << "    ^_^ ��" << i+1 <<"��������Է��������! ����ʱ��" << (double)(ct3 - ct2)/CLOCKS_PER_SEC << "��"  << endl;
		}

		//-----------------------------------------------------------------------------------------------------------------------------------------		
		//�������ڱ߽��������������ڵ�Ľ�
		Solv->Add_periodical_solution(peri_bnods, nodes, U_Solution);

		//-----------------------------------------------------------------------------------------------------------------------------------------		
		//���λ�ƽ����ڼ��
		//for(int i=0; i<9; i++)
		//{
		//	hout << endl << "//==============================================================" << endl;	
		//	hout << "    ^_^ ��" << i+1 <<"��λ�ƽ⣺"<< endl;		
		//	for(int j=0; j<(int)nodes.size(); j++)
		//	{
		//		hout << setwp(6) << j << "  ";
		//		hout << setwp(15,6) << U_Solution[i][3*j] << "  " << setwp(15,6) << U_Solution[i][3*j+1] << "  " << setwp(15,6) << U_Solution[i][3*j+2]<< endl;
		//	}
		//}
		ct1 = clock();
		hout << "    ѭ������������Է����鴦���ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
		hout << "^_^ ѭ������������Է����鴦�����!" << endl << endl;
		delete Solv;		//ɾ���ⷽ������
	}
	else if(wr_mod!="read_solution")	{ hout << "ע�⣡��дָ����Ƕ�read_all��read_solutionҲ����дָ����飡" << endl;	return 0; }

	//-----------------------------------------------------------------------------------------------------------------------------------------		
	//��ȡλ�ƽ�
	if(wr_mod=="read_solution") read_u_solution(backup_file_name);

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�����ЧӦ����, ע�Ȿ�����л���ҪGaussDats�������
	if(Glosmat->Calculate_equivalent_energy_simple(backup_iz, backup_ig, backup_matrix, backup_equright, backup_ege, U_Solution, equivalent_energy)==0) return 0; //��������նȾ���
//	if(Glosmat->Calculate_equivalent_energy(elements, mats, decay, U_Solution, cell, com_mod, equivalent_energy)==0) return 0; //��������նȾ���
	delete Glosmat;		//ɾ�������ܸն�����

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���λ�ƽ⼰��ЧӦ�����������ڼ�¼
	write_u_solution_equivalent_energy(backup_file_name);
	
	return 1;
}
//---------------------------------------------------------------------------
//���������������
int Algorithm_FEM::Import_computational_data(ifstream &infile)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�������ģ�����Ͳ���
	istringstream istr_com_mod(Get_Line(infile));
	istr_com_mod >> com_mod; 
	if(com_mod!="fem"&&com_mod!="nonlocal"&&com_mod!="hybrid")
	{ 
		hout <<"ע�⣡����ģ��ѡ�������������������������ģ�Ͳ�����" << endl; 
		return 0; 
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//λ�ƽ�Ķ�дģʽ����
	istringstream istr_wr_mod(Get_Line(infile));
	istr_wr_mod >> wr_mod;
	if(wr_mod=="read_all")	hout <<"ע�⣡����λ�ƽ⼰��ЧӦ���ܱ������ݣ�" << endl; 
	else if(wr_mod=="read_solution")	hout <<"ע�⣡����λ�ƽⱸ�����ݣ�" << endl; 
	else if(wr_mod!="write")
	{
		hout <<"ע�⣡λ�ƽ�Ķ�дģʽ����������������������дģʽ������" << endl; 
		return 0;
	}

	return 1;
}
//---------------------------------------------------------------------------
//�жϽڵ����ؽڵ�͵�Ԫ����ص�Ԫ
int Algorithm_FEM::Deter_relative_nodes_elements(vector<Node> &nodes, vector<Element> &elements, const double &dist, const struct RVE_Geo &cell, string mod)
{
	const int ES = (int)elements.size();
	const int NS = (int)nodes.size();
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//����(�������)
	for(int i=0; i<ES; i++) elements[i].relative_eles.clear();
	for(int i=0; i<NS; i++) nodes[i].relative_nods.clear();

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//��Ԫ����ص�Ԫ
	if(mod=="fem")
	{
		for(int i=0; i<ES; i++) 
			elements[i].relative_eles.push_back(i);
	}
	else if(mod=="hybrid"||mod=="nonlocal") 
	{
		//���ɵ�Ԫ�����ĵ�����
		vector<Point_3D> centrele(ES);
		for(int i=0; i<ES; i++)
		{
			Point_3D poi_tem(0,0,0);
			for(int j=0; j<8; j++)  //�����嵥Ԫ�ڵ����8
			{
				poi_tem.x += nodes[elements[i].nodes_id[j]].x;
				poi_tem.y += nodes[elements[i].nodes_id[j]].y;
				poi_tem.z += nodes[elements[i].nodes_id[j]].z;
			}
			centrele[i] = poi_tem/8.0;
		}

		//��֤��dist <= n*cell.delt_x(y,z)������ܶԵ�Ԫ��ÿ����˹�㶼ȡ����ȫ��dist����
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
	//�ڵ����ؽڵ�(����С��������)
	//����������֤�������·�ʽдOpenMP���룬ִ�����
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
							//���ַ�����
							int left = 0;
							int right = (int)rela_nods.size()-1;
							while(right>=left)
							{
								int middle = (left + right)/2;
								if(rela_nods[middle] == renod) goto Node_Num_Same_0; //�ڵ�����ͬ�����
								else if(rela_nods[middle] > renod) right = middle - 1;
								else left = middle + 1;
							}

							//�����ַ�ʽ����vector.insert�����Կ�һЩ
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
						//���ַ�����
						int right = (int)nodes[temp_num[i]].relative_nods.size()-1;
						while(right>=left)
						{
							int middle = (left + right)/2;
							if(nodes[temp_num[i]].relative_nods[middle] == temp_rnod[i][j]) goto Node_Num_Same_1; //�ڵ�����ͬ�����
							else if(nodes[temp_num[i]].relative_nods[middle] > temp_rnod[i][j]) right = middle - 1;
							else left = middle + 1;
						}

						//�����ַ�ʽ����vector.insert�����Կ�һЩ
						nodes[temp_num[i]].relative_nods.push_back(0);
						for(int n=(int)nodes[temp_num[i]].relative_nods.size()-1; n>left; n--) nodes[temp_num[i]].relative_nods[n] = nodes[temp_num[i]].relative_nods[n-1];
						nodes[temp_num[i]].relative_nods[left] = temp_rnod[i][j];

Node_Num_Same_1: ;					
					}
				}
			}
		}
	}

	//���³�����û��ִ��OpenMP�ĳ���
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
//						//���ַ�����
//						int left = 0;
//						int right = (int)nodes[nodn].relative_nods.size()-1;
//						while(right>=left)
//						{
//							int middle = (left + right)/2;
//							if(nodes[nodn].relative_nods[middle] == renod) goto Node_Num_Same; //�ڵ�����ͬ�����
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

	//��������ȽϺͼ��
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
//���ݶ������д��������ļ����Ա��Ժ󵥶�������Է�����ʱ�õ�
int Algorithm_FEM::Write_read_matrix_equright_data(const long &Total_Matrix_Size, vector<double> &total_matrix, vector<double>* nl_equright, double backup_ege[], string &maeq_mod)
{
	if(maeq_mod=="write")
	{
		//------------------------------------------------------------------------
		//����д��������ļ�
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
		//���ݶ���������ļ�, ���ڼ����ЧӦ����ʱ��
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
//������ȡλ�ƽ�
void Algorithm_FEM::read_u_solution(const string &output_file_name)
{
	//-------------------------------------------------------------------------------------------------
	//�����ƶ�ȡ����
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
//���λ�ƽ⼰��ЧӦ����
void Algorithm_FEM::write_u_solution_equivalent_energy(const string &output_file_name)const
{
	//-------------------------------------------------------------------------------------------------
	//�����ƴ洢����
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
//����һ����Ϣ��������ע���У���"%"��ͷ����
string Algorithm_FEM::Get_Line(ifstream &infile)const
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
