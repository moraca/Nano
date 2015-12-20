//===========================================================================
// MatBase.cpp
// ���Ͽ���ĳ�Ա����
// Member functions in a class of material database
//===========================================================================
#include "MatBase.h"

//---------------------------------------------------------------------------
//���ɲ������ݿ�
int MatBase::Generate_matbase(ifstream &infile)
{
	//---------------------------------------------------------------------------
	//��ȡ���ϵĸնȲ���
	MatPro Smat;
	//��ȡ���ϸն���Ϣ
	int mat_type;
	istringstream imat0(Get_Line(infile));
	imat0 >> mat_type;

	if(mat_type==0)//����ͬ�Բ���
	{
		double E, Nu;
		imat0 >> E >> Nu;
		Smat.set_ela_para(E, Nu);
	}
	else if(mat_type==1)//��۸���ͬ��
	{
		double E1, E2, Nu1, Nu2, G2;
		imat0 >> E1 >> E2 >> Nu1 >> Nu2 >> G2;
		Smat.set_ela_para(E1, E2, Nu1, Nu2, G2);
	}
	else if(mat_type==2)//������������
	{
		double E1, E2, E3, Nu12, Nu23, Nu13, G12, G23, G13;
		imat0 >> E1 >> E2 >> E3 >> Nu12 >> Nu23 >> Nu13 >> G12 >> G23 >> G13;
		Smat.set_ela_para(E1, E2, E3, Nu12, Nu23, Nu13, G12, G23, G13);
	}
	else
	{
		hout << "��ȡ������������" << mat_type << "����0��2֮�䡣" << endl;
		return 0;
	}
	if(Smat.Generate_local_elas_matrix()==0) return 0; //���ɲ��ϵĸնȾ���
	mats_vec.push_back(Smat);

	//---------------------------------------------------------------------------
	//��ȡ�����㳤�����ĵ�Ч�նȲ���
	MatPro known_mat;
	//��ȡ���ϸն���Ϣ
	istringstream imat(Get_Line(infile));
	imat >> mat_type;

	if(mat_type==0)//����ͬ�Բ���
	{
		double E;
		imat >> E;
		known_mat.set_ela_para(E);
	}
	else if(mat_type==1)//��۸���ͬ��(���벴�ɱ�)
	{
		double E1,E2,Nu2;
		imat >> E1 >> E2 >> Nu2;
		known_mat.set_ela_para(E1, E2, Nu2);
	}
	else if(mat_type==3)//��۸���ͬ��(�������ģ��)
	{
		double E1,E2,G2;
		imat >> E1 >> E2 >> G2;
		known_mat.set_ela_para_transverse(E1, E2, G2);
	}
	else if(mat_type==2)//������������
	{
		double E1,E2,E3,Nu12,Nu23,Nu13;
		imat >> E1 >> E2 >> E3 >> Nu12 >> Nu23 >> Nu13;
		known_mat.set_ela_para(E1, E2, E3, Nu12, Nu23, Nu13);
	}
	else
	{
		hout << "��ȡ��������Ч�նȲ�������" << mat_type << "����0��2֮�䡣" << endl;
		return 0;
	}

	if(known_mat.Generate_nonlocal_elas_matrix()==0) return 0; //���ɲ��ϵĸնȾ���

	//---------------------------------------------------------------------------
	//��ȡ˥����������
	istringstream in2(Get_Line(infile));
	in2 >> decay.R >> decay.radius;
	if(decay.R<0||decay.radius<=0) { hout << "ע�⣡˥����������С��0����С�ڵ���0�����������룡" << endl; return 0; }

	//�޸Ĳ���ϵ�������ݷǾֲ�������
	if(decay.R>0.0) 
	{
		MatPro mat;
		//---------------------------------------------------------------------------
		//���ɸ�������Ԥ���ֲ�ģ�͵�Ч�ն���
		Mesher grid;
		double grid_vol;
		if(grid.Generate_grids_for_effective_stiffness(infile, decay.R, grid_vol)==0) return 0;

		//ȡ��˹�㼰��Ȩϵ��
		Gauss gau;
		if(gau.Generate_gauss(infile)==0) return 0;
		
		//Ԥ���ֲ�ģ�͵�Ч�ն���
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

		//���¼���
		double Cm[6] = {	known_mat.elas_matrix[0][0],	known_mat.elas_matrix[1][1],	known_mat.elas_matrix[2][2], 
										known_mat.elas_matrix[0][1], known_mat.elas_matrix[1][2], known_mat.elas_matrix[0][2] };

		MathMatrix Dm(6,6);
		Dm = (Km.Transpose()*Km).Inverse()*Km.Transpose();
		for(int i=0; i<6; i++)
		{
			decay.acoe[i] = 0.0;
			for(int j=0; j<6; j++) decay.acoe[i] += Dm.element[i][j]*Cm[j];
		}

		if(Estimate_stiffness_long_range_interaction(grid.elements, grid.nodes, grid_vol, gau.gauss, gau.weight, decay.acoe, mat, 0)==0) return 0; //���һλ�β��������1����ʾ���error����

		//---------------------------------------------------------------------------
		//�������۹�ʽ����ϵ��������֤һ���ԣ�����Ҫ����Horizon�뾶������Horizon�뾶Ӧ�ô���1������������
//		mat.Compare_coef_by_analysis_formula(mat_type, decay);
		//---------------------------------------------------------------------------
		//�������۹�ʽ����նȾ��󣬲���֤һ����
//		mat.Compare_matrix_by_analysis_formula(decay);

		//---------------------------------------------------------------------------
		//���ݵ��Ծ��󣬷�����ϲ���
		mat.Get_ele_para_by_ela_matrix();

		mats_vec.push_back(mat); //�����������
	}
	else  istringstream istr_precision(Get_Line(infile));   //��������

	//---------------------------------------------------------------------------
	//������в��ϵĸնȾ��󼰶�Ӧ������ģ�������ɱ�
	Print_stiffness_matrix("Material_stiffness_matrix.dat");
	//���ݵ��Ծ��󣬷�����ϲ���
	known_mat.Get_ele_para_by_ela_matrix();
	Print_stiffness_matrix("Material_stiffness_matrix.dat", known_mat, 1);

	return 1;
}
//------------------------------------------------------------------------------
//������в��ϵĸնȾ����ļ���
void MatBase::Print_stiffness_matrix(string print_name)const
{
	ofstream opri(print_name.c_str());
	//���ϸնȾ������
	for(int num=0; num<(int)mats_vec.size(); num++)
	{
		opri << "Material " << num+1 << "  stiffness_matrix:" << endl;
		for(int i=0; i<6; i++)
		{
			for(int j=0; j<6; j++) opri << setprecision(18) << setw(24) << setiosflags(ios::right)  << mats_vec[num].elas_matrix[i][j] << "  ";
			opri << endl;
		}
		opri << endl;

		//���ϵ�����ģ�������ɱȺͼ���ģ�����
		opri << "%Young's Modulus��E11��E22��E33 ; Poisson's Ratio��Nu12��Nu23��Nu13 ; Shear Modulus��G12��G23��G13 " << endl;
		opri << mats_vec[num].E11 << " " << mats_vec[num].E22 << " " << mats_vec[num].E33 << "   ";
		opri << mats_vec[num].Nu12 << " " << mats_vec[num].Nu23 << " " << mats_vec[num].Nu13 << "   ";
		opri << mats_vec[num].G12 << " " << mats_vec[num].G23 << " " << mats_vec[num].G13 << "   ";

		opri << endl << endl;
	}

	opri.close();
}
//------------------------------------------------------------------------------
//������ϵĸնȾ����ļ���
void MatBase::Print_stiffness_matrix(string print_name, MatPro &mat, int key)const
{
	ofstream opri;
	if(key==0) opri.open(print_name.c_str());
	else opri.open(print_name.c_str(), ios_base::app);
	//���ϸնȾ������
	opri << "Compared Material  stiffness_matrix:" << endl;
	for(int i=0; i<6; i++)
	{
		for(int j=0; j<6; j++) opri << setprecision(18) << setw(24) << setiosflags(ios::right)  << mat.elas_matrix[i][j] << "  ";
		opri << endl;
	}
	opri << endl;

	//���ϵ�����ģ�������ɱȺͼ���ģ�����
	opri << "%Young's Modulus��E11��E22��E33 ; Poisson's Ratio��Nu12��Nu23��Nu13 ; Shear Modulus��G12��G23��G13 " << endl;
	opri << mat.E11 << " " << mat.E22 << " " << mat.E33 << "   ";
	opri << mat.Nu12 << " " << mat.Nu23 << " " << mat.Nu13 << "   ";
	opri << mat.G12 << " " << mat.G23 << " " << mat.G13 << "   ";

		opri << endl << endl;

	opri.close();
}
//------------------------------------------------------------------------------
//Ԥ���ֲ�ģ�͵�Ч�ն���
int MatBase::Estimate_stiffness_long_range_interaction(const vector<Element> &elements, const vector<Node> &nodes, const double &grid_vol,
						const vector<Node> &gauss, const vector<double> &weight, const double a[], MatPro &mat, const int &output_key)const
{
	//--------------------------------------------------	
	//��������Ԫ��ż���˹������
	if((int)elements.size()==0) { hout << "Ԥ���ն�����������Ԫ��Ϊ0, ���飡" << endl; return 0; }
	int nume = ((int)elements.size()-1)/2;

	//--------------------------------------------
	//��˹�����
	const int gau_num = (int)weight.size();

	//--------------------------------------------
	//��������Ԫ��˹������
	Point_3D gd_temp(0,0,0);
	vector<Point_3D> gd_main(gau_num, gd_temp);

	for(int count=0; count<gau_num; count++)
	{
		double Nshape[8] = {0};
		
		//�����˹���ֵ���������
		Nshape[0]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
		Nshape[1]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0-gauss[count].z);
		Nshape[2]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
		Nshape[3]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0-gauss[count].z);
		Nshape[4]=0.125*(1.0-gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
		Nshape[5]=0.125*(1.0+gauss[count].x)*(1.0-gauss[count].y)*(1.0+gauss[count].z);
		Nshape[6]=0.125*(1.0+gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);
		Nshape[7]=0.125*(1.0-gauss[count].x)*(1.0+gauss[count].y)*(1.0+gauss[count].z);

		//�����˹���ֵ�����
		for(int j=0; j<8; j++) 
		{
			gd_main[count].x += Nshape[j]*nodes[elements[nume].nodes_id[j]].x;
			gd_main[count].y += Nshape[j]*nodes[elements[nume].nodes_id[j]].y;
			gd_main[count].z += Nshape[j]*nodes[elements[nume].nodes_id[j]].z;
		}
	}

	//--------------------------------------------------	
	//��������Ч�ն���
	double (*Long_range_stiffness)[6][6] = new double [gau_num][6][6];
	//��������Ч�ն����ʼ��Ϊ��
	for(int count =0; count<gau_num; count++)		
		for(int i=0; i<6; i++)
			for(int j=0; j<6; j++)
				Long_range_stiffness[count][i][j] = 0.0;

	//--------------------------------------------------	
	//ѭ�����е�Ԫ
	for(int i=0; i<(int)elements.size(); i++)
	{
		//ѭ������Ԫ��˹��������
		for(int count1=0; count1<gau_num; count1++)
		{
			//ѭ���ⵥԪ��˹��������
			for(int count2=0; count2<gau_num; count2++)
			{
				double Nshape[8] = {0};
				//--------------------------------------------
				//�����˹���ֵ���������
				Nshape[0]=0.125*(1.0-gauss[count2].x)*(1.0-gauss[count2].y)*(1.0-gauss[count2].z);
				Nshape[1]=0.125*(1.0+gauss[count2].x)*(1.0-gauss[count2].y)*(1.0-gauss[count2].z);
				Nshape[2]=0.125*(1.0+gauss[count2].x)*(1.0+gauss[count2].y)*(1.0-gauss[count2].z);
				Nshape[3]=0.125*(1.0-gauss[count2].x)*(1.0+gauss[count2].y)*(1.0-gauss[count2].z);
				Nshape[4]=0.125*(1.0-gauss[count2].x)*(1.0-gauss[count2].y)*(1.0+gauss[count2].z);
				Nshape[5]=0.125*(1.0+gauss[count2].x)*(1.0-gauss[count2].y)*(1.0+gauss[count2].z);
				Nshape[6]=0.125*(1.0+gauss[count2].x)*(1.0+gauss[count2].y)*(1.0+gauss[count2].z);
				Nshape[7]=0.125*(1.0-gauss[count2].x)*(1.0+gauss[count2].y)*(1.0+gauss[count2].z);

				//--------------------------------------------
				//�����˹���ֵ�����
				Point_3D gd_assist(0,0,0);
				for(int j=0; j<8; j++) 
				{
					gd_assist.x += Nshape[j]*nodes[elements[i].nodes_id[j]].x;
					gd_assist.y += Nshape[j]*nodes[elements[i].nodes_id[j]].y;
					gd_assist.z += Nshape[j]*nodes[elements[i].nodes_id[j]].z;
				}				
			
				//--------------------------------------------
				//���㳤������˥������ֵ
				const double x = gd_assist.x-gd_main[count1].x;
				const double y = gd_assist.y-gd_main[count1].y;
				const double z = gd_assist.z-gd_main[count1].z;	

				const double dis_squr = x*x + y*y + z*z;
				const double poi_dis = sqrt(dis_squr);		//����֮��ľ���
			
				if(poi_dis>decay.R+Zero||poi_dis<Zero) continue;		//Բ�λ�����
			
				const double cos2sita = z*z/dis_squr;		//cos(sita)^2
				double cos2pha;
				if(fabs(x)<Zero&&fabs(y)<Zero) cos2pha = 1.0;
				else cos2pha = x*x/(x*x + y*y);					//cos(pha)^2

				double sum = a[0] + a[1]*0.5*(3*cos2sita-1) + a[2]*(2*cos2pha-1)*3*(1-cos2sita) + a[3]*0.125*(35*cos2sita*cos2sita - 30*cos2sita + 3)
										+ a[4]*(2*cos2pha-1)*7.5*(7*cos2sita-1)*(1-cos2sita) + a[5]*(8*cos2pha*cos2pha-8*cos2pha+1)*105*(1-cos2sita)*(1-cos2sita);

				const double gv = exp(-poi_dis/decay.radius)*sum*weight[count2]; //˥���Լ�Ȩ�غ���ֵ
//				const double gv = sum*weight[count2]/dis_squr; //˥���Լ�Ȩ�غ���ֵ  (���ڶԱȷ�����ʽ����ֵ��ʽ)

				const double r_comp[6] = {x*x, y*y, z*z, x*y, y*z, z*x}; //��Ч������˥������
				const double temp_gmat[6] = {gv*r_comp[0], gv*r_comp[1], gv*r_comp[2], gv*r_comp[3], gv*r_comp[4], gv*r_comp[5]}; //˥����������ĶԳ���

				//--------------------------------------------
				//���㳤������Ч�նȾ���
				for(int j=0; j<6; j++)
					for(int k=0; k<6; k++)
						Long_range_stiffness[count1][j][k] += temp_gmat[j]*r_comp[k]; //�൱��gv*r_comp[j]*r_comp[k]
			}
		}
	}

	//�����ſɱ�����ʽ��ֵ
	//--------------------------------------------
	//�κ���N�����ĵ��ƫ������
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
	//��Ԫ�ڵ��������
	double elenode[8][3];
	for(int j=0; j<8; j++)
	{
		elenode[j][0]=nodes[elements[nume].nodes_id[j]].x;
		elenode[j][1]=nodes[elements[nume].nodes_id[j]].y;
		elenode[j][2]=nodes[elements[nume].nodes_id[j]].z;
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
	double J_val = Jmatrix[0][0]*(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])
								-Jmatrix[0][1]*(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])
								+Jmatrix[0][2]*(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0]);

	//--------------------------------------------------
	//�������������С����
	double error_max[6][6], error_min[6][6];	

	//--------------------------------------------------
	//ѭ������Ԫ��˹����������˹���ϵĵ�Ч����������
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
		//����������
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

		//���������ڼ��
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
//����һ����Ϣ��������ע���У���"%"��ͷ����
string MatBase::Get_Line(ifstream &infile)const
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
