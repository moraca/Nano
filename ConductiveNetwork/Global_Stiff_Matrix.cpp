//===========================================================================
// Global_Stiff_Matrix.cpp
// ����ܸ������Ա����
// Member Functions in a Class of the Global Stiff Matrix
//===========================================================================

#include "Global_Stiff_Matrix.h"

//---------------------------------------------------------------------------
//��������նȾ���
int Global_Stiff_Matrix::Gen_global_stiff_matrix(ifstream &infile, const vector<Point_3D> &cnps, const vector<Node> &nodes, const vector<Element> &elements, const vector<MatPro> &mats, const struct Decay_Para &decay, 
																						const struct RVE_Geo &cell, const struct CNT_Geo &cnts, const string &com_mod, const vector<long> &Iz, const vector<int> &Ig, vector<double> &total_matrix, vector<double>* nl_equright, double backup_ege[])
{
	//---------------------------------------------------------------------------
	//ȡ��˹��
	Gauss gau;		//��˹������
	if(gau.Generate_gauss(infile)==0) return 0;

	clock_t ct0,ct1;
	//------------------------------------------------------------------------
	//���㵥Ԫÿ����˹������ݺ�������׹�Ȩ��
	ct0 = clock();
	hout << "-_- ��ʼ�����˹���������ݺ�Ȩ��" << endl;

	const int GS = (int)gau.gauss.size();
	const int ES = (int)elements.size();
	double (*gauss_ns)[8] = new double [GS][8];				//��¼��Ԫ��˹����κ���
	double (*gauss_nw)[8] = new double [GS][8];				//��¼��Ԫ��˹����κ���(��Ȩ��ϵ��)
	double (*gauss_po)[3] = new double [GS][3];				//��׼�����嵥Ԫ�ĸ�˹������
	double (*ele_cent)[3] = new double [ES][3];					//��¼��Ԫ���ĵ�λ������(x,y,z)�ֱ����[0],[1]��[2]�У��Ա�׼�����嵥Ԫ���ĵ�Ϊԭ�㣩
	double Jacobi = 0;																//��׼�����嵥Ԫ���ſɱ�ֵ
	double (*gauss_dfx)[3][8] = new double [GS][3][8];		//��¼��Ԫ��˹���κ����ĵ���
	Generate_element_gauss_data(elements, nodes, com_mod, gau.gauss, gau.weight, gauss_ns, gauss_nw, gauss_dfx, Jacobi, gauss_po, ele_cent);
	double (*gauss_wcnt)[10] = new double[ES*GS][10];	//���׹ܵ�Ȩ�أ�[0]��¼����Ȩ��; [1]-[9]�ֱ��¼��������ת�н�����; ����˳��Ϊx'����x, y, z��ļн� [cos(a1), cos(b1), cos(c1)];
																								//y'����x, y, z��ļн�[cos(a2), cos(b2), cos(c2)]; z'����x, y, z��ļн�[cos(a3), cos(b3), cos(c3)]; z'��Ϊ���׹ܷ���, x'��Ϊ��˹�㵽���׹ܴ�ֱ����, y'=z'���x'��
	if(Generate_element_gauss_weight_cnt(infile, decay.R, cell, cnts, cnps, elements, nodes, gau.gauss, gau.weight, gauss_po, ele_cent, gauss_wcnt, Jacobi)==0) return 0;

	//------------------------------------------------------------------------
	//����д��������ļ����Ա���������ЧӦ����ʱ�õ�
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
	hout << "    ���ɸ�˹���������ݺ�Ȩ�غ�ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
	hout << "^_^ �����˹���������ݺ�Ȩ�ز�����ϣ�" << endl << endl;

	//------------------------------------------------------------------------
	//ѭ�����е�Ԫ
	ct0 = clock();
	hout << "-_- ��ʼ���㵥�ղ��ӵ��ܸ�" << endl;
	//ִ��openmp
	#pragma omp parallel
	{
		//���嵥��
		double element_stiff_matrix1[24][24];
		double element_stiff_matrix2[24][24];
		double element_equright[9][24];
		double (*gwc_left)[10] = new double [GS][10];
		double (*gwc_right)[10] = new double [GS][10];

		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++)
		{
			//---------------------------------------------------------------------------
			//����RVE��������Ԫ(����������)�ĵ�������ӵ��ܸ�
			if(com_mod=="hybrid"||com_mod=="nonlocal")
			{
				//��ʼ��ʱ
				clock_t ctn1,ctn2;
				ctn1 = clock();

				//--------------------------------------------------
				//��ʼ������
				double elec_left[3] = { ele_cent[i][0], ele_cent[i][1], ele_cent[i][2] };

				//--------------------------------------------------
				//��Ԫ��ÿ����˹���������������׹ܵ�Ȩ��ֵ
				const int IGS = i*GS;
				for(int j=0; j<GS; j++)
					for(int k=0; k<10; k++)
						gwc_left[j][k] = gauss_wcnt[IGS+j][k];

				//------------------------------------------------------------------------
				//����������
				for(int j=0; j<(int)elements[i].relative_eles.size(); j++)
				{
					//--------------------------------------------------
					//��ʼ������
					const int ere = elements[i].relative_eles[j];
					double elec_right[3] = { ele_cent[ere][0], ele_cent[ere][1], ele_cent[ere][2] };

					//�ж���Щ�������ڱ߽�������������ص�Ԫ�����
					bool mark = false;
					double peri_dist[3] = {0};
					if(elec_right[0]-elec_left[0] > cell.len_x-decay.R-Zero) { elec_right[0] -= cell.len_x; peri_dist[0] = -cell.len_x; mark = true; }
					else if(elec_right[0]-elec_left[0] < -cell.len_x+decay.R+Zero) { elec_right[0] += cell.len_x; peri_dist[0] = cell.len_x; mark = true; }

					if(elec_right[1]-elec_left[1] > cell.wid_y-decay.R-Zero) { elec_right[1] -= cell.wid_y; peri_dist[1] = -cell.wid_y; mark = true; }
					else if(elec_right[1]-elec_left[1] < -cell.wid_y+decay.R+Zero) { elec_right[1] += cell.wid_y; peri_dist[1] = cell.wid_y; mark = true; }

					if(elec_right[2]-elec_left[2] > cell.hei_z-decay.R-Zero) { elec_right[2] -= cell.hei_z; peri_dist[2] = -cell.hei_z; mark = true; }
					else if(elec_right[2]-elec_left[2] < -cell.hei_z+decay.R+Zero) { elec_right[2] += cell.hei_z; peri_dist[2] = cell.hei_z; mark = true; }

					//--------------------------------------------------
					//��Ԫ��ÿ����˹���������������׹ܵ�Ȩ��ֵ
					const int EREGS = ere*GS;
					for(int k=0; k<GS; k++)
						for(int m=0; m<10; m++)
							gwc_right[k][m] = gauss_wcnt[EREGS+k][m];

					//------------------------------------------------------------------------
					//�߽絥Ԫ���Գ�������ص�Ԫ������ƽ�Ƶ����ⵥԪ
					if(mark)
					{
						double ege[9] = { 0 };
						Generate_Longforce_Equivalent_Equright(element_equright, Jacobi, gauss_nw, gau.weight, decay, gauss_po, elec_left, elec_right, gwc_left, gwc_right, peri_dist, ege);
						//������ӵ��Ҷ���
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
					//����RVE��������Ԫ(����������)������
					Generate_Longforce_Elestiff(element_stiff_matrix1, element_stiff_matrix2, gauss_ns, gauss_nw, gau.weight, decay, gauss_po, elec_left, elec_right, gwc_left, gwc_right);

					//������ӵ��ܸ�
					#pragma omp critical
					{
						Add_to_gsmatrix(element_stiff_matrix1, Iz, Ig, total_matrix, elements[i]);	//������openmp��, ���Բ��ܷ���0ֵ, ��ֹ����
						Add_to_gsmatrix(element_stiff_matrix2, Iz, Ig, total_matrix, elements[i], elements[elements[i].relative_eles[j]]);	//������openmp��, ���Բ��ܷ���0ֵ, ��ֹ����
					}
				}

				//------------------------------------------------------------------------
				ctn2 = clock();
				hout << "Total num of elements: " << (int)elements.size() << "; Element " << i << " took time: " << (double)(ctn2-ctn1)/CLOCKS_PER_SEC << "sec; " << endl;
			}

			//---------------------------------------------------------------------------
			//����RVE��������Ԫ(�Ӵ�������)�ĵ�������ӵ��ܸ�
			//����RVE����Ԫ(�Ӵ�������)�����󣨾�������Ԫ������
			if(com_mod=="hybrid"||com_mod=="fem")
			{
				if(com_mod=="fem")  //�ֲ�����, ���׹���ǿ����
				{
					//--------------------------------------------------
					//��Ԫ��ÿ����˹���������������׹ܵ�Ȩ��ֵ
					const int IGS = i*GS;
					for(int j=0; j<GS; j++)
						for(int k=0; k<10; k++)
							gwc_left[j][k] = gauss_wcnt[IGS+j][k];

					CNT_Reinforcement_Local_Continuum(element_stiff_matrix1, mats, Jacobi, gauss_dfx, gau.weight, gwc_left);
					
					//������ӵ��ܸ�
					#pragma omp critical
					{					
						Add_to_gsmatrix(element_stiff_matrix1, Iz, Ig, total_matrix, elements[i]);  //������openmp��, ���Բ��ܷ���0ֵ, ��ֹ����
					}
				}

				Generate_Contactforce_Elestiff(element_stiff_matrix1, mats, Jacobi, gauss_dfx, gau.weight);
				//������ӵ��ܸ�
				#pragma omp critical
				{
					Add_to_gsmatrix(element_stiff_matrix1, Iz, Ig, total_matrix, elements[i]);  //������openmp��, ���Բ��ܷ���0ֵ, ��ֹ����
				}
			}
		}
		//---------------------------------------------------------------------------
		//ɾ��ָ��
		delete[] gwc_left;
		delete[] gwc_right;
	}

	//---------------------------------------------------------------------------
	//ɾ��ָ��
	delete[] gauss_ns;
	delete[] gauss_nw;
	delete[] gauss_dfx;
	delete[] gauss_po;
	delete[] ele_cent;
	delete[] gauss_wcnt;

	//---------------------------------------------------------------------------
	ct1 = clock();
	hout << "    ���㵥�ղ��ӵ��ܸպ�ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
	hout << "^_^ ���㵥�ղ��ӵ��ܸղ�����ϣ�" << endl << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//����RVE��������Ԫ(�Ӵ�������)������(������)��CNT��ǿ�����ֲ�����ģ�Ͳ��֣�
void Global_Stiff_Matrix::CNT_Reinforcement_Local_Continuum(double (*element_stiff_matrix)[24], const vector<MatPro> &mats, const double &Jacobi, const double (*gauss_dfx)[3][8], 
																													 const vector<double> &weight, const double (*gwc_left)[10])const
{
	//--------------------------------------------------
	//��ʼ��
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++) 
			element_stiff_matrix[i][j] = 0;
	//--------------------------------------------------
	//ѭ����˹��������
	const int ws = (int)weight.size();
	for(int count=0; count<ws; count++)
	{
		if(gwc_left[count][0]==0) continue;
		//-----------------------------------------------------------------
		//�նȾ������ת�任����
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
		//����˵�Ԫ����Ӧ�Ĳ��ϵ��Ծ���
		double temp_ele[6][6];
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
			{
				temp_ele[i][j] = 0.0;
				for(int k=0; k<6; k++)
				temp_ele[i][j] += Ro[i][k]*mats[1].elas_matrix[k][j];
			}

		//-----------------------------------------------------------------
		//����Ro��ת��
		double ele_elas[6][6];
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
			{
				ele_elas[i][j] = 0.0;
				for(int k=0; k<6; k++)
				ele_elas[i][j] += temp_ele[i][k]*Ro[j][k];
			}

		//--------------------------------------------------------
		//B����
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
		//����B_trans
		double B_trans[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
				B_trans[i][j] = B[j][i];
		//--------------------------------------------------------------------------------------------
		//���B_trans������ele_elas����ĳ˻�array1
		double array1[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
			{
				array1[i][j]=0; 
				for(int k=0; k<6; k++)
					array1[i][j] += B_trans[i][k]*ele_elas[k][j];
			}
		//���array1������B����ĳ˻�array2
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

	//�ڸ�˹��ѭ������ſɱ�ֵ���ٳ˷�����
	for(int j=0; j<24; j++) 
		for (int k=0; k<24; k++) 
			element_stiff_matrix[j][k] = element_stiff_matrix[j][k]*Jacobi;

	//������������ڼ��
	//hout << "element stiffness matrix: " << endl;
	//for(int i=0; i<24; i++){
	//	for(int j=0; j<24; j++)
	//		hout << element_stiff_matrix[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl << endl;
}

//-----------------------------------------------------------------------------------------------
//����RVE��������Ԫ(�Ӵ�������)������(������)
void Global_Stiff_Matrix::Generate_Contactforce_Elestiff(double (*element_stiff_matrix)[24], const vector<MatPro> &mats, const double &Jacobi, const double (*gauss_dfx)[3][8], const vector<double> &weight)const
{
	//--------------------------------------------------
	//��ʼ��
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++) 
			element_stiff_matrix[i][j] = 0;
	//--------------------------------------------------
	//ѭ����˹��������
	const int ws = (int)weight.size();
	for(int count=0; count<ws; count++)
	{
		//-----------------------------------------------------------------
		//����˵�Ԫ����Ӧ�Ĳ��ϵ��Ծ���
		double ele_elas[6][6];		
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				ele_elas[i][j] = mats[0].elas_matrix[i][j];
		
		//--------------------------------------------------------
		//B����
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
		//����B_trans
		double B_trans[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
				B_trans[i][j] = B[j][i];
		//--------------------------------------------------------------------------------------------
		//���B_trans������ele_elas����ĳ˻�array1
		double array1[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
			{
				array1[i][j]=0; 
				for(int k=0; k<6; k++)
					array1[i][j] += B_trans[i][k]*ele_elas[k][j];
			}
		//���array1������B����ĳ˻�array2
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

	//�ڸ�˹��ѭ������ſɱ�ֵ���ٳ˷�����
	for(int j=0; j<24; j++) 
		for (int k=0; k<24; k++) 
			element_stiff_matrix[j][k] = element_stiff_matrix[j][k]*Jacobi;

	//������������ڼ��
	//hout << "element stiffness matrix: " << endl;
	//for(int i=0; i<24; i++){
	//	for(int j=0; j<24; j++)
	//		hout << element_stiff_matrix[i][j] << " ";
	//		hout << endl;
	//}
	//hout << endl << endl << endl;
}
//-----------------------------------------------------------------------------------------------
//����RVE��������Ԫ(����������)������(������)
void Global_Stiff_Matrix::Generate_Longforce_Elestiff(double (*element_stiff_matrix1)[24], double (*element_stiff_matrix2)[24], const double (*gauss_ns)[8], const double (*gauss_nw)[8], const vector<double> &weight, 
																								   const struct Decay_Para &decay, const double (*gauss_po)[3], const double elec_left[], const double elec_right[], const double (*gwc_left)[10], const double (*gwc_right)[10])const
{
	//--------------------------------------------------
	//��ʼ��
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++)
		{
			element_stiff_matrix1[i][j] = 0;
			element_stiff_matrix2[i][j] = 0;
		}

	//--------------------------------------------------	
	//ѭ������Ԫ��˹��������
	const int gau_num = (int)weight.size();
	for(int count1=0; count1<gau_num; count1++)
	{
		//--------------------------------------------------	
		double Gmat[6] = {0};
		double GNmat[6][8] = { {0}, {0}, {0}, {0}, {0}, {0} };
		//--------------------------------------------------	
		//��˸�˹������
		Point_3D gaupoi_left(0, 0, 0);
		gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
		gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
		gaupoi_left.z = gauss_po[count1][2] + elec_left[2];

		//------------------------------------------------------------------------------------------------------------------------
		//ѭ���ⵥԪ��˹��������
		for(int count2=0; count2<gau_num; count2++)
		{
			//--------------------------------------------
			//�Ҷ˸�˹������
			Point_3D gaupoi_right(0, 0, 0);
			gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
			gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
			gaupoi_right.z = gauss_po[count2][2] + elec_right[2];
	
			//--------------------------------------------
			//Ȩ��ֵ�ж�
			if(gwc_left[count1][0]<Zero&&gwc_right[count2][0]<Zero) continue;  //������ֵ��û�г���Ч��

			//--------------------------------------------
			//���㳤������˥������ֵ
			const double xleft = gaupoi_right.x-gaupoi_left.x;
			const double yleft = gaupoi_right.y-gaupoi_left.y;
			const double zleft = gaupoi_right.z-gaupoi_left.z;

			//--------------------------------------------
			//����bond�ĳ���
			const double dis_squr = xleft*xleft + yleft*yleft + zleft*zleft;
			const double poi_dis = sqrt(dis_squr);		//����֮��ľ���
			
			if(poi_dis>decay.R+Zero||poi_dis<Zero) continue;		//Բ�λ�����

			//--------------------------------------------
			//����left��Ȩ��
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
			//����right��Ȩ��
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
			//����Ȩ�غ���ֵ
			const double sum = 0.5*(gwc_left[count1][0]*sum_left+gwc_right[count2][0]*sum_right);
			if(fabs(sum)<Zero) continue;  //������ֵ��û�г���Ч��

			//--------------------------------------------
			//���㳤��������
			const double gv = exp(-poi_dis/decay.radius)*sum*weight[count2]; //˥���Լ�Ȩ�غ���ֵ

			const double temp_gmat[6] = {gv*xleft*xleft, gv*yleft*yleft, gv*zleft*zleft, gv*xleft*yleft, gv*yleft*zleft, gv*zleft*xleft}; //˥����������ĶԳ���(��xleft����xright�ĳ˻�����һ����)

			for(int i=0; i<6; i++) Gmat[i] += temp_gmat[i];  //����count2����

			for(int i=0; i<6; i++)
				for(int j=0; j<8; j++) GNmat[i][j] +=  temp_gmat[i]*gauss_ns[count2][j]; //����count2�����κ����˻�����
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

	////������������ڼ��
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
//���ɱ߽絥Ԫ(����������)������(������)��Ч�����Ҷ�����
void Global_Stiff_Matrix::Generate_Longforce_Equivalent_Equright(double (*element_equright)[24], const double &Jacobi, const double (*gauss_nw)[8], const vector<double> &weight, const struct Decay_Para &decay,
																														const double (*gauss_po)[3], const double elec_left[], const double elec_right[], const double (*gwc_left)[10], const double (*gwc_right)[10], const double peri_dist[], double ege[])const
{
	//--------------------------------------------------
	//��ʼ��
	for(int i=0; i<9; i++)
		for(int j=0; j<24; j++)
			element_equright[i][j] = 0;

	//--------------------------------------------------	
	//ѭ������Ԫ��˹��������
	double TGmat[6] = { 0 };
	double NGmat[8][6] = { {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0} };
	const int gau_num = (int)weight.size();
	for(int count1=0; count1<gau_num; count1++)
	{
		//--------------------------------------------------	
		double Gmat[6] = {0};

		//--------------------------------------------------	
		//��˸�˹������
		Point_3D gaupoi_left(0, 0, 0);
		gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
		gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
		gaupoi_left.z = gauss_po[count1][2] + elec_left[2];

		//------------------------------------------------------------------------------------------------------------------------
		//ѭ���ⵥԪ��˹��������
		for(int count2=0; count2<gau_num; count2++)
		{
			//--------------------------------------------
			//�Ҷ˸�˹������
			Point_3D gaupoi_right(0, 0, 0);
			gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
			gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
			gaupoi_right.z = gauss_po[count2][2] + elec_right[2];
	
			//--------------------------------------------
			//Ȩ��ֵ�ж�
			if(gwc_left[count1][0]<Zero&&gwc_right[count2][0]<Zero) continue;  //������ֵ��û�г���Ч��

			//--------------------------------------------
			//���㳤������˥������ֵ
			const double xleft = gaupoi_right.x-gaupoi_left.x;
			const double yleft = gaupoi_right.y-gaupoi_left.y;
			const double zleft = gaupoi_right.z-gaupoi_left.z;

			//--------------------------------------------
			//����bond�ĳ���
			const double dis_squr = xleft*xleft + yleft*yleft + zleft*zleft;
			const double poi_dis = sqrt(dis_squr);		//����֮��ľ���
			
			if(poi_dis>decay.R+Zero||poi_dis<Zero) continue;		//Բ�λ�����

			//--------------------------------------------
			//����left��Ȩ��
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
			//����right��Ȩ��
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
			//����Ȩ�غ���ֵ
			const double sum = 0.5*(gwc_left[count1][0]*sum_left+gwc_right[count2][0]*sum_right);
			if(fabs(sum)<Zero) continue;  //������ֵ��û�г���Ч��

			//--------------------------------------------
			//���㳤��������			
			const double gv = exp(-poi_dis/decay.radius)*sum*weight[count2]; //˥���Լ�Ȩ�غ���ֵ

			//------------------------------------------------------------------------------------------------
			//ע������ temp_gmat[6]�ĸ�ֵ˳��ͬ����������������ʱtemp_gmat[6]��˳��
			const double temp_gmat[6] = { gv*xleft*xleft, gv*xleft*yleft, gv*yleft*yleft, gv*xleft*zleft, gv*yleft*zleft, gv*zleft*zleft }; //˥����������ĶԳ���

			for(int i=0; i<6; i++) Gmat[i] += temp_gmat[i];  //����count2����
		}

		for(int i=0; i<8; i++)
			for(int j=0; j<6; j++)
				NGmat[i][j] += gauss_nw[count1][i]*Gmat[j];

		for(int i=0; i<6; i++) TGmat[i] += Gmat[i]*weight[count1];
	}

	for(int i=0; i<9; i++)
	{
		double E[3][3] = {{0}, {0}, {0}}; //�������ڱ߽���������ʱ�ı��ξ�������
		switch(i)	//���þ��Ȼ�Ӧ����Ϊ���ڱ߽�����
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
		default: hout << "���� �������ڱ߽�����Լ��ʱ��ѭ������ֵ����" << i << "С��0���ߴ���8�����飡" << endl;
		}

		//���Ȼ�λ�Ʋ�
		double uni_disp[3] = {0};
		for(int j=0; j<3; j++)
			for(int k=0; k<3; k++) 
				uni_disp[j] += E[j][k]*peri_dist[k];

		for(int j=0; j<3; j++) uni_disp[j] = uni_disp[j]*Jacobi; //Ϊ�˼��ٳ˷��Ĵ����������ȳ����ſɱȳ���

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
//��������ӵ��ܸ�
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
				//���ַ�����λ��
				bool mark = false;
				long left = Iz[row_node-1];
				long middle = 0;
				long right = Iz[row_node];
				while(right>=left)
				{
					middle = (left + right)/2;
					if(Ig[middle] == col_node) { mark = true; break; } //�ҵ���Ӧ�ڵ�
					else if(Ig[middle] > col_node) right = middle - 1;
					else left = middle + 1;
				}					
				if(!mark) { hout << "�����ڽڵ�" << row_node << "��Iz, Ig��û���ҵ��ڵ�" << col_node << "���飡" << endl; }
				
				const long Mt = 6*(long)row_node + 9*middle;	//�ҵ���Ӧһά�洢����ʼλ��
if( Mt>=(long)total_matrix.size() ) hout << "$$$1" <<  Mt << " row_now: " << row_node << " middle: " << middle << "  col_node: " << col_node << " Iz[row_node-1]: "  <<  Iz[row_node-1] << " Iz[row_node]: "  <<  Iz[row_node]<< endl;
				for(int k=0; k<=2; k++)
					for(int m=0; m<=2; m++)
						total_matrix[Mt+3*k+m] += element_stiff_matrix[3*i+k][3*j+m];
			}
			else
			{
				const long Mt = 6*(long)row_node + 9*Iz[row_node]; //�ҵ���Ӧһά�洢����ʼλ��
if( Mt>=(long)total_matrix.size() ) hout << "$$$2" <<  Mt << " row_now: " << row_node << "  col_node: " << col_node << " Iz[row_node-1]: "  <<  Iz[row_node-1] << " Iz[row_node]: "  <<  Iz[row_node] << endl;							
				for(int k=0; k<=2; k++)
					for(int m=0; m<=k; m++)	//ע��˴���С�ڵ���k, ��Ϊ����������
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
//��������ӵ��ܸ�(���е�Ԫ�ڵ㲻ͬ)
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
				//���ַ�����λ��
				bool mark = false;
				long left = Iz[row_node-1];
				long middle = 0;
				long right = Iz[row_node];
				while(right>=left)
				{
					middle = (left + right)/2;
					if(Ig[middle] == col_node) { mark = true; break; } //�ҵ���Ӧ�ڵ�
					else if(Ig[middle] > col_node) right = middle - 1;
					else left = middle + 1;
				}					
				if(!mark) { hout << "�����ڽڵ�" << row_node << "��Iz, Ig��û���ҵ��ڵ�" << col_node << "���飡" << endl; }
				
				const long Mt = 6*(long)row_node + 9*middle;	//�ҵ���Ӧһά�洢����ʼλ��

				for(int k=0; k<=2; k++)
					for(int m=0; m<=2; m++)
						total_matrix[Mt+3*k+m] += element_stiff_matrix[3*i+k][3*j+m];
			}
			else
			{
				const long Mt = 6*(long)row_node + 9*Iz[row_node]; //�ҵ���Ӧһά�洢����ʼλ��
				
				for(int k=0; k<=2; k++)
					for(int m=0; m<=k; m++)	//ע��˴���С�ڵ���k, ��Ϊ����������
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
//���㵥Ԫÿ����˹���������������׹ܵ�Ȩ��ֵ
int Global_Stiff_Matrix::Generate_element_gauss_weight_cnt(ifstream &infile, const double &dist, const struct RVE_Geo &cell, const struct CNT_Geo &cnts, const vector<Point_3D> &cnps, const vector<Element> &elements,
																											 const vector<Node> &nodes, const vector<Node> &gauss, const vector<double> &weight, const double (*gauss_po)[3], const double (*ele_cent)[3], double (*gauss_wcnt)[10], double &Jacobi)const
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//����Ȩ�غ��������þ����Լ�Ȩ�غ����Ľ���
	double w_dist, w_constant;
	int w_order;
	istringstream istrw(Get_Line(infile));
	istrw >> w_dist;
	if(w_dist<0&&w_dist!=-1) { hout << "ע�⣡Ȩ�غ�����������þ������" << w_dist << "��С��0�Ҳ�����-1����������������������룡" << endl; return 0; }
	if(w_dist>=0)
	{
		istrw >> w_order;
		if(w_order>3||w_order<-2) { hout << "ע�⣡Ȩ�غ������ݴε���" << w_order << "������3���ݻ��ߵ���-2���ݣ����������룡" << endl; return 0; }
		if(w_order==0)
		{
			istrw >> w_constant;		// �����0���ݱ�ʾ���������ﻹ��Ҫ���볣��ֵ
			if(w_constant>1.0||w_constant<0.0) { hout << "ע�⣡Ȩ�غ������ݴε���0�����ǳ���ϵ���������󣬴���1����С��0���߸���û�����룡" << endl; return 0; }
		}
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//������Ԫ����ص�Ԫ, ���ڼ����˹������׹ܵ���С����
	vector<vector<int> > relemd;
	Deter_relative_elements_min_dist(nodes, elements, cell, w_dist, relemd);

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�����˹�㵽������׹ܵľ���ͷ���
	double total_cnt_volume = 0;
	double total_interface_volume = 0;
	double total_value_gps = 0;
	#pragma omp parallel
	{
		//������ʱȨ�ر���
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
				temp_num.push_back(i*GS+count);  //��¼��˹�ڵ���һά�洢�е������ţ���omp��洢���õ�
				//�������, ÿһ����˹�㶼��Ȩ�ض���1������ת��Ϊ��λ��
				if(w_dist==-1) 
				{ 
					temp_gwcnt[0].push_back(1.0);
					//��֤��ת����ʱ��λ����
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
				Point_3D gaupoi(0, 0, 0);		//��˹���ֵ�
				gaupoi.x = gauss_po[count][0] + ele_cent[i][0];
				gaupoi.y = gauss_po[count][1] + ele_cent[i][1];
				gaupoi.z = gauss_po[count][2] + ele_cent[i][2];

				//--------------------------------------------
				//ѭ����Ԫ����ص�Ԫ
				double min_dist = w_dist;
				Point_3D min_ep[2];		//��С�����Ӧ�߶ε������˵�
				Point_3D min_fp;				//��С�����Ӧ����
				for(int j=0; j<(int)relemd[i].size(); j++) //����i��Ԫ����
				{
					int ere = relemd[i][j];
					//��֤��dist <= n*cell.delt_x(y,z)������ܶԵ�Ԫ��ÿ����˹�㶼ȡ����ȫ��w_dist����
					double distx, disty, distz;
					distx = (int(w_dist/cell.delt_x-Zero)+1)*cell.delt_x;
					disty = (int(w_dist/cell.delt_y-Zero)+1)*cell.delt_y;
					distz = (int(w_dist/cell.delt_z-Zero)+1)*cell.delt_z;
					//�ж���Щ�������ڱ߽�������������ص�Ԫ�����
					Point_3D delt_dist(0,0,0);
					if(ele_cent[ere][0]-ele_cent[i][0] > cell.len_x-distx-Zero) delt_dist.x = ele_cent[ere][0] - cell.len_x;
					else if(ele_cent[ere][0]-ele_cent[i][0] < distx+Zero-cell.len_x) delt_dist.x = ele_cent[ere][0] + cell.len_x;

					if(ele_cent[ere][1]-ele_cent[i][1] > cell.wid_y-disty-Zero) delt_dist.y = ele_cent[ere][1] - cell.wid_y;
					else if(ele_cent[ere][1]-ele_cent[i][1] < disty+Zero-cell.wid_y) delt_dist.y = ele_cent[ere][1] + cell.wid_y;

					if(ele_cent[ere][2]-ele_cent[i][2] > cell.hei_z-distz-Zero) delt_dist.z = ele_cent[ere][2] - cell.hei_z;
					else if(ele_cent[ere][2]-ele_cent[i][2] < distz+Zero-cell.hei_z) delt_dist.z = ele_cent[ere][2] + cell.hei_z;

					//ѭ����Ԫ��������׹��߶�
					for(int k=0; k<(int)elements[ere].relative_cnts.size(); k++)
					{
						int num = elements[ere].relative_cnts[k];
						//�߶ε������˵�
						Point_3D point[2] = { delt_dist+cnps[num], delt_dist+cnps[num+1] };

						//����ռ�ֱ����һ��Ĳ���
						double ratio_t = ((gaupoi.x-point[0].x)*(point[1].x-point[0].x)+(gaupoi.y-point[0].y)*(point[1].y-point[0].y)+(gaupoi.z-point[0].z)*(point[1].z-point[0].z))				
													/((point[1].x-point[0].x)*(point[1].x-point[0].x)+(point[1].y-point[0].y)*(point[1].y-point[0].y)+(point[1].z-point[0].z)*(point[1].z-point[0].z));

						//����ֱ�����������㣨���㣩
						Point_3D footpoi(point[0].x+ratio_t*(point[1].x-point[0].x), point[0].y+ratio_t*(point[1].y-point[0].y), point[0].z+ratio_t*(point[1].z-point[0].z));

						if(ratio_t>0.0&&ratio_t<1.0)
						{
							//�������
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
								//�������
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

						//���������˹��ʹ�����ͬһ��, ˵����˹�������׹�����, ������Ҫ�޸�min_fp
						if(gaupoi.distance_to(footpoi)<Zero)
						{
							//����point[1]-point[0]Ϊ����������gaupoiΪͨ�����ƽ�������ѡȡһ�㣬���㴹��㣨��x'�ᣩ
							//�������һ������
							int seed = rand()%MAX_INT;
							seed = (2053*seed + 13849)%MAX_INT;

							Point_3D gpoi(0,0,0);
							gpoi.x = 2.0*seed/MAX_INT-1.0;  //����[-1, 1]֮���һ����;
							seed = (2053*seed + 13849)%MAX_INT;
							gpoi.y = 2.0*seed/MAX_INT-1.0;  //����[-1, 1]֮���һ����;
							if(fabs(point[1].z-point[0].z)>Zero) gpoi.z =  ((point[1].x-point[0].x)*(gpoi.x-gaupoi.x) + (point[1].y-point[0].y)*(gpoi.y-gaupoi.y))/(point[0].z-point[1].z) + gaupoi.z;
							else gpoi.z = 0.0;  //����ֵ�����ԣ���Ϊϵ��(point[1].z-point[0].z)����0

							min_fp = gaupoi + gaupoi - gpoi;
						}
					}
				}
				//����Ȩ��ֵ(ǰ��ֻ��¼��w_distС��, ���Բ����ܴ���w_dist)
				if(min_dist<cnts.rad_max)	temp_cnt_volume += weight[count];  //��¼��˹���Ȩ��֮��, ʵ�����׹���ռ�������
				if(min_dist==w_dist)
				{
					temp_gwcnt[0].push_back(0.0);
					//��֤��ת����ʱ��λ����
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
					temp_interface_volume += weight[count];  //��¼��˹���Ȩ��֮��, ʵ�ʽ������ռ�������

					if(w_order==0) temp_gwcnt[0].push_back(w_constant);
					else if(w_order==1) temp_gwcnt[0].push_back(1.0-min_dist/w_dist);
					else if(w_order==3) temp_gwcnt[0].push_back(1.0+min_dist*min_dist*(2*min_dist-3*w_dist)/(w_dist*w_dist*w_dist));
					else if(w_order==-1) temp_gwcnt[0].push_back(2.0/(1.0+min_dist/w_dist)-1.0);
					else if(w_order==-2) temp_gwcnt[0].push_back(4.0/((1.0+min_dist/w_dist)*(1.0+min_dist/w_dist))-4.0/(1.0+min_dist/w_dist)+1.0);
					else	{	hout << "w_order=" << w_order << " û�ж�Ӧ���ݵ�Ȩ�غ������㹫ʽ�����������룡" << endl; }

					//������ת����ֵ
					double tgve[3][4] = {{0}, {0}, {0}};
					
					//x'��
					tgve[0][0] = gaupoi.x - min_fp.x;
					tgve[0][1] = gaupoi.y - min_fp.y;
					tgve[0][2] = gaupoi.z - min_fp.z;
					tgve[0][3] = sqrt(tgve[0][0]*tgve[0][0]+tgve[0][1]*tgve[0][1]+tgve[0][2]*tgve[0][2]);
					
					//z'��
					tgve[2][0] = min_ep[1].x - min_ep[0].x;
					tgve[2][1] = min_ep[1].y - min_ep[0].y;
					tgve[2][2] = min_ep[1].z - min_ep[0].z;
					tgve[2][3] = sqrt(tgve[2][0]*tgve[2][0]+tgve[2][1]*tgve[2][1]+tgve[2][2]*tgve[2][2]);
					
					//y'�� (z'���x')
					tgve[1][0] = tgve[2][1]*tgve[0][2] - tgve[2][2]*tgve[0][1];
					tgve[1][1] = tgve[2][2]*tgve[0][0] - tgve[2][0]*tgve[0][2];
					tgve[1][2] = tgve[2][0]*tgve[0][1] - tgve[2][1]*tgve[0][0];
					tgve[1][3] = sqrt(tgve[1][0]*tgve[1][0]+tgve[1][1]*tgve[1][1]+tgve[1][2]*tgve[1][2]);
					
					//������ת����
					for(int j=0; j<3; j++)
						for(int k=0; k<3; k++)
							temp_gwcnt[j*3+k+1].push_back(tgve[j][k]/tgve[j][3]);
				}
				else { hout << "ע�⣡������С��������趨���룬���飡" << endl; }
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
				//Ϊ���10��Ȩ��ֵ��ȷ���
				//if(gauss_wcnt[temp_num[i]][0]<-Zero||gauss_wcnt[temp_num[i]][0]>1.0+Zero) hout << "ע�����gauss_wcnt[" << temp_num[i] << "][0]=" << gauss_wcnt[temp_num[i]][0] << "�����飡" << endl;
				//for(int j=0; j<3; j++)
				//{
				//	double sum[2] = { 0 };
				//	for(int k=0; k<3; k++)	sum[0] +=  gauss_wcnt[temp_num[i]][3*j+k+1]*gauss_wcnt[temp_num[i]][3*j+k+1];
				//	for(int k=0; k<3; k++)	sum[1] +=  gauss_wcnt[temp_num[i]][j+3*k+1]*gauss_wcnt[temp_num[i]][j+3*k+1];
				//	if(fabs(sum[0]-1.0)>Zero||fabs(sum[1]-1.0)>Zero) hout << "ע�����sum[0]=" << sum[0] << "sum[1]=" << sum[1] << "�����飡" << endl;
				//}
			}
		}
	}

	hout << "    ʵ�����׹���ռ���������Ϊ��" << total_cnt_volume*Jacobi << endl;
	hout << "    ʵ�ʽ������ռ���������Ϊ��" << total_interface_volume*Jacobi << endl;
	hout << "    ���׹ܽ�����ƽ��Ȩ��ֵ��" << total_value_gps/((int)elements.size()*(int)gauss.size()) << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------
//������Ԫ����ص�Ԫ, ���ڼ����˹������׹ܵ���С����
void Global_Stiff_Matrix::Deter_relative_elements_min_dist(const vector<Node> &nodes, const vector<Element> &elements, const struct RVE_Geo &cell, const double &dist, vector<vector<int> > &relemd)const
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//��Ԫ����ص�Ԫ
	const int ES = (int)elements.size();
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
//���㵥Ԫÿ����˹�������
void Global_Stiff_Matrix::Generate_element_gauss_data(const vector<Element> &elements, const vector<Node> &nodes, const string &com_mod, const vector<Node> &gauss, const vector<double> &weight,
																									double (*gauss_ns)[8], double (*gauss_nw)[8], double (*gauss_dfx)[3][8], double &Jacobi,  double (*gauss_po)[3], double (*ele_cent)[3])const
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//ѭ����Ԫ
	Point_3D scenter(0,0,0);	//��׼�����嵥Ԫ�����ĵ�����
	int ie = 0; //��0�ŵ�ԪΪ��׼�����嵥Ԫ
	for(int count=0; count<(int)gauss.size(); count++)
	{	
		double Nshape[8] = {0};
		//--------------------------------------------
		//�����˹���ֵ���������
		Nshape[0]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
		Nshape[1]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
		Nshape[2]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
		Nshape[3]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
		Nshape[4]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
		Nshape[5]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
		Nshape[6]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);
		Nshape[7]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);

		//--------------------------------------------
		//�����˹���ֵ�����
		Point_3D gaupoi(0, 0, 0);		
		for(int j=0; j<8; j++) 
		{
			gaupoi.x += Nshape[j]*nodes[elements[ie].nodes_id[j]].x;
			gaupoi.y += Nshape[j]*nodes[elements[ie].nodes_id[j]].y;
			gaupoi.z += Nshape[j]*nodes[elements[ie].nodes_id[j]].z;
		}

		//��¼��˹���ֵ�����
		gauss_po[count][0] = gaupoi.x;
		gauss_po[count][1] = gaupoi.y;
		gauss_po[count][2] = gaupoi.z;

		//--------------------------------------------
		//����ʾ���
		//--------------------------------------------
		//�κ���N��gauss[count].x, gauss[count].y, gauss[count].z��ƫ������
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
		//��Ԫ�ڵ��������
		double elenode[8][3];
		for(int j=0; j<8; j++)
		{
			elenode[j][0]=nodes[elements[ie].nodes_id[j]].x;
			elenode[j][1]=nodes[elements[ie].nodes_id[j]].y;
			elenode[j][2]=nodes[elements[ie].nodes_id[j]].z;
		}
		//--------------------------------------------------
		//J����
		double Jmatrix[3][3];
		//������������Ļ�
		for(int j=0; j<3; j++)
			for(int k=0; k<3; k++)
			{
				Jmatrix[j][k]=0;
				for(int m=0; m<8; m++)
				Jmatrix[j][k] += diff[j][m]*elenode[m][k];
			}
		//--------------------------------------------------
		//���J���������ʽ
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
			//���J����������
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
			//���N��x,y,z��ƫ��
			double diffxy[3][8];
			for(int j=0; j<3; j++)
				for(int k=0; k<8; k++)
				{
					diffxy[j][k]=0;
					for(int m=0; m<3; m++)
						diffxy[j][k] += Jinverse[j][m]*diff[m][k];
				
					//��¼J���������ʽֵ
					gauss_dfx[count][j][k] = diffxy[j][k];
				}
		}
	}

	//����õ�Ԫ�����ĵ�����
	for(int j=0; j<8; j++) 
	{
		scenter.x += nodes[elements[ie].nodes_id[j]].x;
		scenter.y += nodes[elements[ie].nodes_id[j]].y;
		scenter.z += nodes[elements[ie].nodes_id[j]].z;
	}
	scenter = scenter/8;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//ѭ����Ԫ, ���������Ԫ���ĵ��������꣨�����0�ŵ�Ԫ��
	for(int i=0; i<(int)elements.size(); i++)
	{
		for(int j=0; j<3; j++) ele_cent[i][j] = 0.0;
		//����õ�Ԫ�����ĵ�����
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
//�����ЧӦ����
int Global_Stiff_Matrix::Calculate_equivalent_energy_simple(const vector<long> &Iz, const vector<int> &Ig, const vector<double> &total_matrix, const vector<double>* nl_equright, const double backup_ege[],
																											const vector<vector<double> > &disp_solu,  double equivalent_energy[])const
{
	//------------------------------------------------------------------------
	//������ʼ��
	for(int i=0; i<9; i++) equivalent_energy[i] = 0;

	clock_t ct0,ct1;
	//------------------------------------------------------------------------
	//ѭ�����нڵ����ؽڵ�
	ct0 = clock();
	hout << "-_- ��ʼ�����ЧӦ����" << endl;

	for(int i=0; i<(int)Iz.size(); i++)
	{
		//-----------------------------------------------------------
		//��ص�
		if(i!=0)
		{
			for(long j=Iz[i-1]; j<Iz[i]; j++)
			{
				const long Ij = 6*(long)i + 9*j;
				for(int k=0; k<(int)disp_solu.size(); k++)
					for(int m=0; m<=2; m++)
						for(int n=0; n<=2; n++)
						{
							equivalent_energy[k] += disp_solu[k][3*i+m]*total_matrix[Ij+3*m+n]*disp_solu[k][3*Ig[j]+n];   //����
							equivalent_energy[k] += disp_solu[k][3*Ig[j]+m]*total_matrix[Ij+m+3*n]*disp_solu[k][3*i+n];	 //ת��
						}
			}
		}
		//-----------------------------------------------------------
		//�Խ���
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
	//�Ҷ������ڱ߽���������
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
	hout << "    �����ЧӦ���ܺ�ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
	hout << "^_^ �����ЧӦ���ܲ�����ϣ�" << endl << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------
//�����ЧӦ����
int Global_Stiff_Matrix::Calculate_equivalent_energy(const vector<Element> &elements, const vector<MatPro> &mats, const struct Decay_Para &decay, const vector<vector<double> > &disp_solu, const struct RVE_Geo &cell,
																							 const string &com_mod,  double equivalent_energy[])const
{
	//------------------------------------------------------------------------
	//������ʼ��
	for(int i=0; i<9; i++) equivalent_energy[i] = 0;

	clock_t ct0,ct1;
	//------------------------------------------------------------------------
	//���ݶ���������ļ�, ���ڼ����ЧӦ����ʱ��
	fstream fin( "GaussDats", ios::in|ios::binary);
	int GS, ES; 	//��������	
	fin.read((char *)&GS, sizeof(int));	
	fin.read((char *)&ES, sizeof(int));
	double *weight = new double[GS];									//��¼��Ԫ��˹���Ȩ��ֵ
	fin.read((char *)weight, sizeof(double)*GS);

	double (*gauss_po)[3] = new double [GS][3];				//��׼�����嵥Ԫ�ĸ�˹������
	double (*gauss_ns)[8] = new double [GS][8];				//��¼��Ԫ��˹����κ���
	double (*ele_cent)[3] = new double [ES][3];					//��¼��Ԫ���ĵ�λ������(x,y,z)�ֱ����[0],[1]��[2]�У��Ա�׼�����嵥Ԫ���ĵ�Ϊԭ�㣩
	double (*gauss_wcnt)[10] = new double[ES*GS][10];	//���׹ܵ�Ȩ��
	double Jacobi;																	//��׼�����嵥Ԫ���ſɱ�ֵ
	double (*gauss_dfx)[3][8] = new double [GS][3][8];		//��¼��Ԫ��˹���κ����ĵ���

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
	//ѭ�����е�Ԫ
	ct0 = clock();
	hout << "-_- ��ʼ�����ЧӦ����" << endl;
	//ִ��openmp
	#pragma omp parallel
	{
		//���嵥��
		double elements_energy[9] = {0};
		double (*gwc_left)[10] = new double [GS][10];
		double (*gwc_right)[10] = new double [GS][10];

		#pragma omp for schedule(dynamic, CHUNKSIZE)
		for(int i=0; i<ES; i++)
		{
			//--------------------------------------------------
			//��Ԫ�����λ��ֵ
			double disp_left[9][24];
			for(int j=0; j<9; j++)
				for(int k=0; k<8; k++)
				{
					disp_left[j][3*k] = disp_solu[j][3*elements[i].nodes_id[k]];
					disp_left[j][3*k+1] = disp_solu[j][3*elements[i].nodes_id[k]+1];
					disp_left[j][3*k+2] = disp_solu[j][3*elements[i].nodes_id[k]+2];
				}

			//---------------------------------------------------------------------------
			//����RVE��������Ԫ(����������)�ĵ�ЧӦ����
			if(com_mod=="hybrid"||com_mod=="nonlocal")
			{
				//��ʼ��ʱ
				clock_t ctn1,ctn2;
				ctn1 = clock();

				//--------------------------------------------------
				//��ʼ������
				double elec_left[3] = { ele_cent[i][0], ele_cent[i][1], ele_cent[i][2] };

				//--------------------------------------------------
				//��Ԫ��ÿ����˹���������������׹ܵ�Ȩ��ֵ
				const int IGS = i*GS;
				for(int j=0; j<GS; j++)
					for(int k=0; k<10; k++)
						gwc_left[j][k] = gauss_wcnt[IGS+j][k];

				//------------------------------------------------------------------------
				//����������
				for(int j=0; j<(int)elements[i].relative_eles.size(); j++)
				{
					//--------------------------------------------------
					//��ʼ������
					const int ere = elements[i].relative_eles[j];
					double elec_right[3] = { ele_cent[ere][0], ele_cent[ere][1], ele_cent[ere][2] };

					//--------------------------------------------------
					//��Ԫ�����λ��ֵ
					double disp_right[9][24];
					for(int k=0; k<9; k++)
						for(int m=0; m<8; m++)
						{
							disp_right[k][3*m] = disp_solu[k][3*elements[ere].nodes_id[m]];
							disp_right[k][3*m+1] = disp_solu[k][3*elements[ere].nodes_id[m]+1];
							disp_right[k][3*m+2] = disp_solu[k][3*elements[ere].nodes_id[m]+2];
						}
					//--------------------------------------------------
					//�ж���Щ�������ڱ߽�������������ص�Ԫ����������޸����ֵ
					Verify_Periodical_Condition_Modify_Values(elec_left, elec_right, disp_right, cell, decay);

					//--------------------------------------------------
					//��Ԫ��ÿ����˹���������������׹ܵ�Ȩ��ֵ
					const int EREGS = ere*GS;
					for(int k=0; k<GS; k++)
						for(int m=0; m<10; m++)
							gwc_right[k][m] = gauss_wcnt[EREGS+k][m];

					//------------------------------------------------------------------------
					//����RVE��������Ԫ(����������)��ЧӦ����
					Calculate_Longforce_Energy(elements_energy, weight, gauss_ns, decay, gauss_po, elec_left, elec_right, gwc_left, gwc_right, disp_left, disp_right, GS);
				}
				//------------------------------------------------------------------------
				ctn2 = clock();
				hout << "Equivalent Energy: total num of elements: " << (int)elements.size() << "; Element " << i << " took time: " << (double)(ctn2-ctn1)/CLOCKS_PER_SEC << "sec;" << endl;
			}

			//---------------------------------------------------------------------------
			//����RVE��������Ԫ�Ӵ�����ЧӦ����
			if(com_mod=="hybrid"||com_mod=="fem")
			{
				Calculate_Contactforce_Energy(elements_energy, disp_left, mats, Jacobi, gauss_dfx, weight, GS);
			}
		}

		//��ӵ���Ӧ����
		#pragma omp critical
		{
			for(int i=0; i<9; i++)	equivalent_energy[i] += elements_energy[i];
		}

		//---------------------------------------------------------------------------
		//ɾ��ָ��
		delete[] gwc_left;
		delete[] gwc_right;
	}

	//---------------------------------------------------------------------------
	//ɾ��ָ��
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
	hout << "    �����ЧӦ���ܺ�ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
	hout << "^_^ �����ЧӦ���ܲ�����ϣ�" << endl << endl;

	return 1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------
//�ж���Щ�������ڱ߽�������������ص�Ԫ����������޸����ֵ
void Global_Stiff_Matrix::Verify_Periodical_Condition_Modify_Values(const double elec_left[], double elec_right[], double (*disp_right)[24], const struct RVE_Geo &cell, const struct Decay_Para &decay)const
{
	//---------------------------------------------------------------------------
	//�ж���Щ�������ڱ߽�������������ص�Ԫ
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
			double E[3][3] = {{0}, {0}, {0}}; //�������ڱ߽���������ʱ�ı��ξ�������
			switch(i)	//���þ��Ȼ�Ӧ����Ϊ���ڱ߽�����
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
			default: hout << "���� �������ڱ߽�����Լ��ʱ��ѭ������ֵ����" << i << "С��0���ߴ���8�����飡" << endl;
			}

			//---------------------------------
			//���Ȼ�λ�Ʋ�
			double uni_disp[3] = {0};
			for(int j=0; j<3; j++)
				for(int k=0; k<3; k++) 
					uni_disp[j] += E[j][k]*peri_dist[k];

			//---------------------------------
			//�������λ��ֵ
			for(int j=0; j<8; j++)
				for(int k=0; k<3; k++)
					disp_right[i][3*j+k] += uni_disp[k];
		}
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------
//����RVE��������Ԫ(����������)��ЧӦ����
void Global_Stiff_Matrix::Calculate_Longforce_Energy(double elements_energy[], const double weight[], const double (*gauss_ns)[8], const struct Decay_Para &decay, const double (*gauss_po)[3], const double elec_left[],
																								  const double elec_right[], const double (*gwc_left)[10], const double (*gwc_right)[10], const double (*disp_left)[24], const double (*disp_right)[24], const int &GS)const
{
	//--------------------------------------------------	
	//ѭ������Ԫ��˹��������
	for(int count1=0; count1<GS; count1++)
	{
		//��¼ÿ������Ԫ��˹���ϵĵ�ЧӦ�������ֵ
		double ele_energy[9] = { 0 };

		//--------------------------------------------------	
		double u_left[9][3] = { {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0} };  //��ʼ��
		for(int i=0; i<9; i++) 
			for(int j=0; j<3; j++)
				for(int k=0; k<8; k++)
					u_left[i][j] += gauss_ns[count1][k]*disp_left[i][3*k+j];

		//--------------------------------------------------	
		//��˸�˹������
		Point_3D gaupoi_left(0, 0, 0);
		gaupoi_left.x = gauss_po[count1][0] + elec_left[0];
		gaupoi_left.y = gauss_po[count1][1] + elec_left[1];
		gaupoi_left.z = gauss_po[count1][2] + elec_left[2];

		//------------------------------------------------------------------------------------------------------------------------
		//ѭ���ⵥԪ��˹��������
		for(int count2=0; count2<GS; count2++)
		{
			//--------------------------------------------------	
			double u_right[9][3] = { {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0} };  //��ʼ��
			for(int i=0; i<9; i++) 
				for(int j=0; j<3; j++)
					for(int k=0; k<8; k++)
						u_right[i][j] += gauss_ns[count2][k]*disp_right[i][3*k+j];

			double u_diff[9][3];
			for(int i=0; i<9; i++)
				for(int j=0; j<3; j++)
					u_diff[i][j] = u_right[i][j] - u_left[i][j];

			//--------------------------------------------
			//�Ҷ˸�˹������
			Point_3D gaupoi_right(0, 0, 0);
			gaupoi_right.x = gauss_po[count2][0] + elec_right[0];
			gaupoi_right.y = gauss_po[count2][1] + elec_right[1];
			gaupoi_right.z = gauss_po[count2][2] + elec_right[2];

			//--------------------------------------------
			//Ȩ��ֵ�ж�
			if(gwc_left[count1][0]<Zero&&gwc_right[count2][0]<Zero) continue;  //������ֵ��û�г���Ч��

			//--------------------------------------------
			//���㳤������˥������ֵ
			const double xleft = gaupoi_right.x-gaupoi_left.x;
			const double yleft = gaupoi_right.y-gaupoi_left.y;
			const double zleft = gaupoi_right.z-gaupoi_left.z;

			//--------------------------------------------
			//����bond�ĳ���
			const double dis_squr = xleft*xleft + yleft*yleft + zleft*zleft;
			const double poi_dis = sqrt(dis_squr);		//����֮��ľ���
			
			if(poi_dis>decay.R+Zero||poi_dis<Zero) continue;		//Բ�λ�����

			//--------------------------------------------
			//����left��Ȩ��
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
			//����right��Ȩ��
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
			//����Ȩ�غ���ֵ
			const double sum = 0.5*(gwc_left[count1][0]*sum_left+gwc_right[count2][0]*sum_right);
			if(fabs(sum)<Zero) continue;  //������ֵ��û�г���Ч��

			//--------------------------------------------
			//���㳤��������					
			const double gv = exp(-poi_dis/decay.radius)*sum*weight[count2]; //˥���Լ�Ȩ�غ���ֵ

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
//����RVE��������Ԫ(�Ӵ�������)��ЧӦ����
void Global_Stiff_Matrix::Calculate_Contactforce_Energy(double elements_energy[], const double (*disp)[24], const vector<MatPro> &mats, const double &Jacobi, const double (*gauss_dfx)[3][8], const double weight[], const int &GS)const
{
	//--------------------------------------------------
	//��ʼ��
	double element_stiff_matrix[24][24];
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++) 
			element_stiff_matrix[i][j] = 0;
	//--------------------------------------------------
	//ѭ����˹��������
	for(int count=0; count<GS; count++)
	{
		//-----------------------------------------------------------------
		//����˵�Ԫ����Ӧ�Ĳ��ϵ��Ծ���
		double ele_elas[6][6];		
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				ele_elas[i][j] = mats[0].elas_matrix[i][j];
		
		//--------------------------------------------------------
		//B����
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
		//����B_trans
		double B_trans[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
				B_trans[i][j] = B[j][i];
		//--------------------------------------------------------------------------------------------
		//���B_trans������ele_elas����ĳ˻�array1
		double array1[24][6];
		for(int i=0; i<24; i++)
			for(int j=0; j<6; j++)
			{
				array1[i][j]=0; 
				for(int k=0; k<6; k++)
					array1[i][j] += B_trans[i][k]*ele_elas[k][j];
			}
		//���array1������B����ĳ˻�array2
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

	//�ڸ�˹��ѭ������ſɱ�ֵ���ٳ˷�����
	for(int i=0; i<24; i++) 
		for (int j=0; j<24; j++) 
			element_stiff_matrix[i][j] = element_stiff_matrix[i][j]*Jacobi;

	//--------------------------------------------
	//��װӦ����
	for(int i=0; i<9; i++)
		for(int j=0; j<24; j++)
			for(int k=0; k<24; k++)
				elements_energy[i] += disp[i][j]*element_stiff_matrix[j][k]*disp[i][k];
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//����һ����Ϣ��������ע���У���"%"��ͷ����
string Global_Stiff_Matrix::Get_Line(ifstream &infile)const
{
	string s;
	//������Ϣһ��
	getline(infile,s);
	//����ע����     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===============================================================
