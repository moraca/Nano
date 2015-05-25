//===================================================
//Postprocessor.cpp
//�������
//===================================================
#include "Postprocessor.h"

//--------------------------------------------------------------------------
//ִ��ǰ�������
int Postprocessor::Treatment(ifstream &infile, vector<Node> &nodes, const vector<int>* peri_bnods, const vector<Element> &elements, const vector<MatPro> &mats, 
												   const vector<vector<double> > &u_solution, const string &wr_mod, const string &backup_file_name, double equivalent_energy[], const int &samples_count) //ִ�к������
{
	//��ֵ���ȡλ�ƽ⼰��ЧӦ��������
	if(wr_mod=="write"||wr_mod=="read_solution")	u = u_solution;
	else if(wr_mod=="read_all")  read_u_solution_equivalent_energy(backup_file_name, equivalent_energy);
	else {	hout << "ע�⣡��дָ����Ƕ�read_all��read_solutionҲ����дָ����飡" << endl;	return 0; }
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���λ�Ƴ���ͼ
	Export_Displacement_Contour("Displacement_Field_Contour.dat", nodes, elements);
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//����������ͼ
//	Export_Deformed_Mesh_Multimats("Deformed_Multimats_Mesh.dat", nodes, elements, mats);
	Export_Deformed_Mesh("Deformed_Mesh.dat", nodes, elements);
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���㵥Ԫ���Ĵ�Ӧ������
	Calculate_Elements_Strain(nodes, elements);
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//����ڵ㴦Ӧ������
	Calculate_Nodes_Strain(nodes, peri_bnods, elements);
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���Ӧ�䳡��ͼ
//	Export_Strain_Mesh("Strain_Field_Contour.dat", nodes, elements);
	Export_Strain_Contour("Strain_Field_Contour.dat", nodes, elements);
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���㲢������Ȼ�ϵ��
	Calculate_Export_Homogenized_Coefficients("Homogenized_Coefficitents.dat", elements, equivalent_energy, samples_count);

	return 1;
}
//---------------------------------------------------------------------------
//��ȡλ�ƽ⼰��ЧӦ��������
void Postprocessor::read_u_solution_equivalent_energy(const string &output_file_name, double equivalent_energy[])
{
	//-------------------------------------------------------------------------------------------------
	//�����ƶ�ȡ����
	int us, uis;
	fstream idata(output_file_name.c_str(), ios::in|ios::binary);
	//��ȡλ�ƽ�����
	idata.read((char *)&us, sizeof(int));
	idata.read((char *)&uis, sizeof(int));
	vector<double> temp_u(uis);
	u.assign(us, temp_u);
	for(int i=0; i<us; i++)
		for(int j=0; j<uis; j++)
			idata.read((char *)&u[i][j], sizeof(double));

	//��ȡ��ЧӦ��������
	idata.read((char *)equivalent_energy, sizeof(double)*9);

	idata.close();
}
//---------------------------------------------------------------------------
//���㵥Ԫ���Ĵ�Ӧ������
int Postprocessor::Calculate_Elements_Strain(const vector<Node> &nodes, const vector<Element> &elements)
{
	//---------------------------------------------------------------------------
	//���㵥Ԫ���ĵ㴦��Ӧ������
	vector<double> temp_strain(6*(int)u.size(),0.0);
	ele_strain.assign((int)elements.size(), temp_strain);
	for(int i=0; i<(int)elements.size(); i++)
	{
		//����B����
		//--------------------------------------------
		//�κ���N�������������ĵ�(0,0,0)��ƫ������	
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
			elenode[j][0]=nodes[elements[i].nodes_id[j]].x;
			elenode[j][1]=nodes[elements[i].nodes_id[j]].y;
			elenode[j][2]=nodes[elements[i].nodes_id[j]].z;
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
		//----------------------------------------------------
		//���J����������
		double Jinverse[3][3];
			
		Jinverse[0][0]=(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])/J_val;
		Jinverse[1][1]=(Jmatrix[0][0]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][0])/J_val;
		Jinverse[2][2]=(Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0])/J_val;

		Jinverse[0][1]=-(Jmatrix[0][1]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][1])/J_val;
		Jinverse[0][2]=(Jmatrix[0][1]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][1])/J_val;

		Jinverse[1][0]=-(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])/J_val;
		Jinverse[1][2]=-(Jmatrix[0][0]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][0])/J_val;

		Jinverse[2][0]=(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0])/J_val;
		Jinverse[2][1]=-(Jmatrix[0][0]*Jmatrix[2][1]-Jmatrix[0][1]*Jmatrix[2][0])/J_val;
		//-------------------------------------------------------
		//���N��x,y,z��ƫ��
		double diffxy[3][8];
		for(int j=0; j<3; j++)
			for(int k=0; k<8; k++)
			{
				diffxy[j][k]=0;
				for(int m=0; m<3; m++)
					diffxy[j][k] += Jinverse[j][m]*diff[m][k];
			}
		//--------------------------------------------------------
		//���B����
		double B[6][24];
		for(int j=0; j<6; j++)
			for(int k=0; k<24; k++)
				B[j][k]=0;
		for(int j=0; j<8; j++)
		{
			B[0][j*3+0]=diffxy[0][j];
			B[1][j*3+1]=diffxy[1][j];
			B[2][j*3+2]=diffxy[2][j];
			B[3][j*3+0]=diffxy[1][j];
			B[3][j*3+1]=diffxy[0][j];
			B[4][j*3+1]=diffxy[2][j];
			B[4][j*3+2]=diffxy[1][j];
			B[5][j*3+0]=diffxy[2][j];
			B[5][j*3+2]=diffxy[0][j];
		}		
		//--------------------------------------------
		//���㵥Ԫ���ĵ�Ӧ��ֵ
		vector<double> ele_u(24,0.0);
		for(int j=0; j<(int)u.size(); j++)
		{
			for(int k=0; k<8; k++)
			{
				ele_u[3*k] = u[j][3*elements[i].nodes_id[k]];
				ele_u[3*k+1] = u[j][3*elements[i].nodes_id[k]+1];
				ele_u[3*k+2] = u[j][3*elements[i].nodes_id[k]+2];
			}
			for(int k=0; k<6; k++)
			{
				for(int m=0; m<24; m++)
				{
					ele_strain[i][6*j+k] = ele_strain[i][6*j+k] + B[k][m]*ele_u[m];
				}
			}
		}
	}
		
	return 1;
}
//---------------------------------------------------------------------------
//����ڵ㴦Ӧ������
int Postprocessor::Calculate_Nodes_Strain(vector<Node> &nodes, const vector<int>* &peri_bnods, const vector<Element> &elements)
{
	//����ڵ㴦��Ӧ������
	vector<double> temp_strain(6*(int)u.size(),0.0);
	nod_strain.assign((int)nodes.size(), temp_strain);
	//---------------------------------------------------------------------------
	//ͳ�ƽڵ����ص�Ԫ��Ϣ
	for(int i=0; i<(int)nodes.size(); i++)  //��ص�Ԫ����
	{
		nodes[i].relative_eles.clear();
	}
	for(int i=0; i<(int)elements.size(); i++) //����ڵ����ص�Ԫ
	{
		int node_size = (int)elements[i].nodes_id.size();
		for(int j=0; j<node_size; j++)
		{
			nodes[elements[i].nodes_id[j]].relative_eles.push_back(i);
		}
	}
	//���ǵ������ԣ����ӱ߽�ڵ����ص�Ԫ
	int temp_num = 0;		//��֤����ʼ����
	vector<int> temp_relative_eles;
	for(int i=0; i<(int)peri_bnods[0].size(); i++)
	{
		for(int j=0; j<(int)nodes[peri_bnods[0][i]].relative_eles.size(); j++)
		{
			int nre = nodes[peri_bnods[0][i]].relative_eles[j];
			//���ַ����Ҳ�����
			bool mark = false;
			int left = 0;
			int middle = 0;
			int right = (int)temp_relative_eles.size()-1;
			while(right>=left)
			{
				middle = (left + right)/2;
				if(temp_relative_eles[middle] == nre) { mark = true; break; } //�ڵ�����ͬ�����
				else if(temp_relative_eles[middle] > nre) right = middle - 1;
				else left = middle + 1;
			}
			if(!mark) //û����ͬ��ŵĵ�Ԫ�Ͳ���
			{
				temp_relative_eles.push_back(0);
				for(int k=(int)temp_relative_eles.size()-1; k>left; k--) temp_relative_eles[k] = temp_relative_eles[k-1];
				temp_relative_eles[left] = nre;
			}
		}

		if(i==(int)peri_bnods[0].size()-1||peri_bnods[1][i+1]==0)
		{
			for(int j=temp_num; j<=i; j++)	nodes[peri_bnods[0][j]].relative_eles = temp_relative_eles;

			temp_num = i+1;
			temp_relative_eles.clear(); //��յ�Ԫ����
		}
	}

	//---------------------------------------------------------------------------
	//����ڵ㴦��Ӧ������
	for(int i=0; i<(int)nodes.size(); i++)
	{
		int nrev_size = (int)nodes[i].relative_eles.size();
		for(int j=0; j<6*(int)u.size(); j++)
		{
			for(int k=0; k<nrev_size; k++)
			{
				nod_strain[i][j] = nod_strain[i][j] + ele_strain[nodes[i].relative_eles[k]][j];
			}
			nod_strain[i][j] = nod_strain[i][j]/nrev_size;
		}
	}

	for(int i=0; i<(int)nodes.size(); i++) //��ص�Ԫ����
	{
		nodes[i].relative_eles.clear();
	}

	return 1;
}
//---------------------------------------------------------------------------
//���Tecplot���ӻ�λ�Ƴ���ͼ
int Postprocessor::Export_Displacement_Contour(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE =" << output_file_name << endl;
	otec << "VARIABLES = X, Y, Z";
	if((int)u.size()==1) otec << ", Ux, Uy, Uz" << endl;
	else if((int)u.size()>1)
	{
		for(int i=1; i<=(int)u.size(); i++) otec << ", U" << i << "x, U" << i << "y, U" << i << "z";
		otec << endl;
	}
	
	otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size() << ", F=FEPOINT, ET=BRICK" << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		otec << nodes[i].x << "  " << nodes[i].y << "  " << nodes[i].z;
		for(int j=0; j<(int)u.size(); j++)
		otec << "  " << u[j][3*i] << "  " << u[j][3*i+1] << "  " << u[j][3*i+2];
		otec << endl;
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
//����������ͼ
int Postprocessor::Export_Deformed_Mesh(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE =" << output_file_name << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	for(int i=0; i<(int)u.size(); i++)
	{
		otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size()  << ", F=FEPOINT, ET=BRICK" << endl;
		for(int j=0; j<(int)nodes.size(); j++)
		{
			otec << setprecision(15) << nodes[j].x+u[i][3*j] << "  " << nodes[j].y+u[i][3*j+1] << "  " << nodes[j].z+u[i][3*j+2] << endl;
		}
		for (int j=0; j<(int)elements.size(); j++)
		{
			otec	 << elements[j].nodes_id[0]+1 << "  " << elements[j].nodes_id[1]+1 << "  " 
					 << elements[j].nodes_id[2]+1 << "  " << elements[j].nodes_id[3]+1 << "  " 
					 << elements[j].nodes_id[4]+1 << "  " << elements[j].nodes_id[5]+1 << "  " 
					 << elements[j].nodes_id[6]+1 << "  " << elements[j].nodes_id[7]+1 << endl;
		}
		otec << endl;
	}
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//����������ͼ(���ֲ���)
int Postprocessor::Export_Deformed_Mesh_Multimats(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements, const vector<MatPro> &mats)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE =" << output_file_name << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	//���㵥Ԫ��
	vector<int> elenum(mats.size(),0);
	for (int j=0; j<(int)elements.size(); j++)	elenum[elements[j].mat]++;

	for(int i=0; i<(int)u.size(); i++)
	{
		for(int m=0; m<(int)mats.size(); m++)
		{
			otec << "ZONE N=" << (int)nodes.size() << ", E=" << elenum[m]  << ", F=FEPOINT, ET=BRICK" << endl;
			for(int j=0; j<(int)nodes.size(); j++)
			{
				otec << setprecision(15) << nodes[j].x+u[i][3*j] << "  " << nodes[j].y+u[i][3*j+1] << "  " << nodes[j].z+u[i][3*j+2] << endl;
			}
			for (int j=0; j<(int)elements.size(); j++)
			{
				if(elements[j].mat==m)
				{
					otec	 << elements[j].nodes_id[0]+1 << "  " << elements[j].nodes_id[1]+1 << "  " 
							 << elements[j].nodes_id[2]+1 << "  " << elements[j].nodes_id[3]+1 << "  " 
							 << elements[j].nodes_id[4]+1 << "  " << elements[j].nodes_id[5]+1 << "  " 
							 << elements[j].nodes_id[6]+1 << "  " << elements[j].nodes_id[7]+1 << endl;
				}
			}
			otec << endl;
		}
	}
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//���Ӧ�䳡��ͼ
int Postprocessor::Export_Strain_Mesh(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE =" << output_file_name << endl;
	otec << "VARIABLES = X, Y, Z";
	otec << ", <greek>e</greek><sub>xx</sub>, <greek>e</greek><sub>yy</sub>, <greek>e</greek><sub>zz</sub>";
	otec << ", <greek>e</greek><sub>xy</sub>, <greek>e</greek><sub>yz</sub>, <greek>e</greek><sub>xz</sub>" << endl;
	
	for(int i=0; i<(int)u.size(); i++)
	{
		otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size()  << ", F=FEPOINT, ET=BRICK" << endl;
		for(int j=0; j<(int)nodes.size(); j++)
		{
			otec << setprecision(15) << nodes[j].x+u[i][3*j] << "  " << nodes[j].y+u[i][3*j+1] << "  " << nodes[j].z+u[i][3*j+2];
			otec << setprecision(15) << "  " << nod_strain[j][6*i+0] << "  " << nod_strain[j][6*i+1] << "  " << nod_strain[j][6*i+2];
			otec << setprecision(15) << "  " << nod_strain[j][6*i+3] << "  " << nod_strain[j][6*i+4] << "  " << nod_strain[j][6*i+5] << endl;
		}
		for (int j=0; j<(int)elements.size(); j++)
		{
			otec	 << elements[j].nodes_id[0]+1 << "  " << elements[j].nodes_id[1]+1 << "  " 
					 << elements[j].nodes_id[2]+1 << "  " << elements[j].nodes_id[3]+1 << "  " 
					 << elements[j].nodes_id[4]+1 << "  " << elements[j].nodes_id[5]+1 << "  " 
					 << elements[j].nodes_id[6]+1 << "  " << elements[j].nodes_id[7]+1 << endl;
		}
		otec << endl;
	}
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//���Ӧ�䳡��ͼ(��׼����, ����û�б���)
int Postprocessor::Export_Strain_Contour(const string &output_file_name, const vector<Node> &nodes, const vector<Element> &elements)const
{
	ofstream otec(output_file_name.c_str());
	otec << "TITLE =" << output_file_name << endl;
	otec << "VARIABLES = X, Y, Z";
	otec << ", <greek>e</greek><sub>xx</sub>, <greek>e</greek><sub>yy</sub>, <greek>e</greek><sub>zz</sub>";
	otec << ", <greek>e</greek><sub>xy</sub>, <greek>e</greek><sub>yz</sub>, <greek>e</greek><sub>xz</sub>" << endl;
	
	for(int i=0; i<(int)u.size(); i++)
	{
		otec << "ZONE N=" << (int)nodes.size() << ", E=" << (int)elements.size()  << ", F=FEPOINT, ET=BRICK" << endl;
		for(int j=0; j<(int)nodes.size(); j++)
		{
			otec << setprecision(15) << nodes[j].x << "  " << nodes[j].y << "  " << nodes[j].z;
			otec << setprecision(15) << "  " << nod_strain[j][6*i+0] << "  " << nod_strain[j][6*i+1] << "  " << nod_strain[j][6*i+2];
			otec << setprecision(15) << "  " << nod_strain[j][6*i+3] << "  " << nod_strain[j][6*i+4] << "  " << nod_strain[j][6*i+5] << endl;
		}
		for (int j=0; j<(int)elements.size(); j++)
		{
			otec	 << elements[j].nodes_id[0]+1 << "  " << elements[j].nodes_id[1]+1 << "  " 
					 << elements[j].nodes_id[2]+1 << "  " << elements[j].nodes_id[3]+1 << "  " 
					 << elements[j].nodes_id[4]+1 << "  " << elements[j].nodes_id[5]+1 << "  " 
					 << elements[j].nodes_id[6]+1 << "  " << elements[j].nodes_id[7]+1 << endl;
		}
		otec << endl;
	}
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//���㲢������Ȼ�ϵ��
int Postprocessor::Calculate_Export_Homogenized_Coefficients(const string &output_file_name, const vector<Element> &elements, const double equivalent_energy[], const int &samples_count)const
{
	//---------------------------------------------------------------------------
	//������ļ�
	ofstream ohc;
	if(samples_count==1) ohc.open(output_file_name.c_str());
	else ohc.open(output_file_name.c_str(), ios_base::app);
	ohc << "==================================================================================" <<endl;
	ohc << "No. " << samples_count << " sample:" << endl;

	//---------------------------------------------------------------------------
	//���㵥ԪӦ���
	vector<double> aver_eles_strain(6*(int)u.size(),0.0);
	for(int i=0; i<6*(int)u.size(); i++)
	{
		for(int j=0; j<(int)elements.size(); j++)
		{
			aver_eles_strain[i] += ele_strain[j][i];
		}
		aver_eles_strain[i] = aver_eles_strain[i]/(int)elements.size();
	}

	//---------------------------------------------------------------------------
	//���ƽ����ԪӦ��
	//ohc << "//----------------------------------------------------------------------------------------------------------------------------" << endl;
	//ohc << "Average stains of elements��" << endl;
	//for(int i=0; i<(int)u.size(); i++)
	//{
	//	ohc << "The " << i+1 << " group solutions�� ";
	//	for(int j=0; j<6; j++)
	//	{
	//		if(fabs(aver_eles_strain[6*i+j])<Zero) aver_eles_strain[6*i+j] = 0;
	//		ohc << setw(8) << setprecision(4) << aver_eles_strain[6*i+j];
	//	}
	//	ohc << endl;
	//}

	//----------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------
	//������Ȼ���Ч�նȾ���Ԫ��ֵ

	//---------------------------------------------------------------------------
	//������Ȼ��նȾ���
	if((int)u.size()!=9) { hout << "��������λ�ƽ����������" << u.size() << "������Ԥ��ľ���⣬���飡" << endl; return 0; }
	if(fabs(aver_eles_strain[0])<Zero||fabs(aver_eles_strain[7])<Zero||fabs(aver_eles_strain[14])<Zero||
		fabs(aver_eles_strain[21])<Zero||fabs(aver_eles_strain[28])<Zero||fabs(aver_eles_strain[35])<Zero||
		fabs(aver_eles_strain[36])<Zero||fabs(aver_eles_strain[37])<Zero||fabs(aver_eles_strain[43])<Zero||
		fabs(aver_eles_strain[44])<Zero||fabs(aver_eles_strain[48])<Zero||fabs(aver_eles_strain[50])<Zero)
	{
		hout << "����ƽ��Ӧ��ֵ�б������������㣬���飡" << endl; 
		return 0;
	}

    MathMatrix Homo_Stiff(6,6);
	Homo_Stiff.element[0][0]=equivalent_energy[0]/(aver_eles_strain[0]*aver_eles_strain[0]);
	Homo_Stiff.element[1][1]=equivalent_energy[1]/(aver_eles_strain[7]*aver_eles_strain[7]);
	Homo_Stiff.element[2][2]=equivalent_energy[2]/(aver_eles_strain[14]*aver_eles_strain[14]);
	Homo_Stiff.element[3][3]=equivalent_energy[3]/(aver_eles_strain[21]*aver_eles_strain[21]);
	Homo_Stiff.element[4][4]=equivalent_energy[4]/(aver_eles_strain[28]*aver_eles_strain[28]);
	Homo_Stiff.element[5][5]=equivalent_energy[5]/(aver_eles_strain[35]*aver_eles_strain[35]);

	Homo_Stiff.element[0][1]=0.5*(equivalent_energy[6]-equivalent_energy[0]-equivalent_energy[1])/(aver_eles_strain[36]*aver_eles_strain[37]);
	Homo_Stiff.element[1][0]=Homo_Stiff.element[0][1];

	Homo_Stiff.element[1][2]=0.5*(equivalent_energy[7]-equivalent_energy[1]-equivalent_energy[2])/(aver_eles_strain[43]*aver_eles_strain[44]);
	Homo_Stiff.element[2][1]=Homo_Stiff.element[1][2];

	Homo_Stiff.element[0][2]=0.5*(equivalent_energy[8]-equivalent_energy[0]-equivalent_energy[2])/(aver_eles_strain[48]*aver_eles_strain[50]);
	Homo_Stiff.element[2][0]=Homo_Stiff.element[0][2];

	//������
    MathMatrix S(6,6);
	S=Homo_Stiff.Inverse();

	//����ϲ���
	double E11=1.0/S.element[0][0];
	double E22=1.0/S.element[1][1];
	double E33=1.0/S.element[2][2];
	double Nu12=-S.element[0][1]*E11;
	double Nu13=-S.element[0][2]*E11;
	double Nu23=-S.element[1][2]*E22;
	double G12=1.0/S.element[3][3];
	double G23=1.0/S.element[4][4];
	double G13=1.0/S.element[5][5];

	//---------------------------------------------------------------------------
	//������Ȼ��նȾ���
	ohc << endl << "Homogenized Stiffness Matrix:" << endl;
	ohc << Homo_Stiff.element << endl;

	//���ϵ�����ģ�������ɱȺͼ���ģ�����
	ohc << "%Young's Modulus��E11��E22��E33 ; Poisson's Ratio��Nu12��Nu13��Nu23 ; Shear Modulus��G12��G13��G23 "<<endl;
	ohc << E11<< " " << E22 << " " << E33 << "   ";
	ohc << Nu12<< " " << Nu13 << " " << Nu23 << "   ";
	ohc << G12 << " " <<  G13 << " " <<  G23 << "   ";

	ohc << endl << endl;

	//---------------------------------------------------------------------------
	//�ر�����ļ�
	ohc.close();

	return 1;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//��ȡ�����ļ�һ�У�����ע���У���%��ʼ��
char* Postprocessor::Get_Line(ifstream& input_stream, char* read_line)
{
	//�����һ��char read_line[200]     
	input_stream.getline(read_line,200);
	//����ע����
	while(!input_stream.eof() && read_line[0] == '%')
	{
		input_stream.getline(read_line,200);
	}

	return read_line;
}
//����һ����Ϣ��������ע���У���"%"��ͷ����
string Postprocessor::Get_Line(ifstream &infile)const
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
