//===========================================================================
// Mesher.cpp
// ��ά�����ʷֳ�Ա����
// Member functions in a class of 3D mesh generation
//===========================================================================
#include "Mesher.h"

//---------------------------------------------------------------------------
//�������ɺ�����
int Mesher::Generate_mesh(ifstream &infile, struct RVE_Geo &cell_geo, vector<Point_3D> &cnps)
{
	//���ɱ�������������
	if(Brick_background_mesh(infile, cell_geo)== 0) return 0;
	//һά�洢���ڱ߽�ڵ��ż���ʼ����Ϣ
	if(Record_periodic_boundary_nodes_number(cell_geo)==0) return 0;
	//��������Ԫ��������
	if(Element_material_property(infile)==0) return 0;
	//��¼��Ԫ��������׹��߶�
	if(Element_relative_cnts(cnps, cell_geo)==0) return 0;
	//�����ⵥԪ�Ƿ�ƥ������ȷ��������׹��߶���Ϣ
//	if(Test_element_relative_cnts(cnps)== 0) return 0;
	//���Tecplot���ӻ���������
//	if(Export_mesh_data_tecplot("Mesh_in_tecplot.dat")== 0) return 0;

	return 1;
}
//---------------------------------------------------------------------------
//������������������
int Mesher::Brick_background_mesh(ifstream &infile, struct RVE_Geo &cell_geo)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���뵥���ĳ�����͸ߵķָ�ϸ��
	istringstream istr2(Get_Line(infile));	
	istr2 >> dx >> dy >> dz;
	if(dx<0||dy<0||dz<0)
	{
		hout << "RVE����ߵķָ�ϸ��������� ���飡 ......... Brick_background_mesh" << endl;
		return 0;
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//��������(�ڵ�͵�Ԫ)

	//���ɸ�������ĵ���
	int x_sec_num = int(cell_geo.len_x/dx) +1;
	int y_sec_num = int(cell_geo.wid_y/dy)+1;
	int z_sec_num = int(cell_geo.hei_z/dz) +1;

	//΢�����ʷֲ�������֤�����ĳߴ粻��
	dx = cell_geo.len_x/(x_sec_num-1);
	dy = cell_geo.wid_y/(y_sec_num-1);
	dz = cell_geo.hei_z/(z_sec_num-1);

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//��¼�ָ�ϸ�ȵ�������Ϣ��
	cell_geo.delt_x = dx;
	cell_geo.delt_y = dy;
	cell_geo.delt_z = dz;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���ɽڵ�
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
				nd.type = Deter_node_type(i, j, k, z_sec_num-1, y_sec_num-1, x_sec_num-1);	//��ע�ڵ��λ��(�ǵ㡢�߽��ߵ㡢�߽���㡢�ڵ�)
				nodes.push_back(nd);
			}
		}
	}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���������嵥Ԫ

	//��������������
	for( int i=0; i<z_sec_num-1; i++ )
	{
		for( int j=0; j<y_sec_num-1; j++ )
		{
			for( int k=0; k<x_sec_num-1; k++ )
			{
				//������İ˸�����
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
//��������Ԫ��������
int Mesher::Element_material_property(ifstream &infile)
{
	istringstream istr(Get_Line(infile));		//����ע����
	int define_case;
	double fiber_radius;
	istr >> define_case;

	if(define_case!=0&&define_case!=1) { hout <<"ע�⣺����Ĳ��ϵ�Ԫ���������������....Element_material_property 1" << endl;	return 0; }
	if(define_case==1) istr >> fiber_radius;
	if(fiber_radius<0.0) { hout <<"ע�⣺��ά���ϵİ뾶�����������....Element_material_property 2" << endl;	return 0; }

	if(define_case==0)	//������
	{
		for(int i=0; i<(int)elements.size(); i++)
		{
			elements[i].type = 381;		//��Ӧ��Ԫ����״����(��������xyz�� x ά�ȣ�y��Ԫ�Ľڵ������z��Ԫ���κ����ݴ�)��
														//���磬121: һά���ڵ�(�߶�)�����κ�����231: ��ά���ڵ�(������)�����κ�����241: ��ά�Ľڵ�(�ı���)�����κ�����
														//341: ��ά�Ľڵ�(������)�����κ�����361: ��ά���ڵ�(������)�����κ�����381����ά�˽ڵ�(������)�����κ���
			elements[i].mat = 0;			//��Ӧ�������Ա��
		}
	}
	else if(define_case==1) //Z�ᵥ����ά
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
//��¼��Ԫ��������׹��߶�
int Mesher::Element_relative_cnts(vector<Point_3D> &cnps, const struct RVE_Geo &cell_geo)
{
	//���ɸ�������Ķ���(�뱳����������ɲ���һ��)
	int x_sec_num = int(cell_geo.len_x/dx)+1;
	int y_sec_num = int(cell_geo.wid_y/dy)+1;
	int z_sec_num = int(cell_geo.hei_z/dz)+1;

	for(int i=0; i<(int)cnps.size()-1; i++)
	{
		if(cnps[i+1].flag==0) continue;
		
		int n_temp;
		Point_3D point[2] = {cnps[i], cnps[i+1]};
		
		int nx[2] = { int((point[0].x-cell_geo.poi_min.x)/dx), int((point[1].x-cell_geo.poi_min.x)/dx) };		
		if(nx[0]>nx[1]) { n_temp = nx[0]; nx[0] = nx[1]; nx[1] = n_temp; } //����С����
		if(nx[1]==x_sec_num-1) nx[1] -= 1; //���߽��ϵĵ�����ǰһ����

		int ny[2] = { int((point[0].y-cell_geo.poi_min.y)/dy), int((point[1].y-cell_geo.poi_min.y)/dy) };
		if(ny[0]>ny[1]) { n_temp = ny[0]; ny[0] = ny[1]; ny[1] = n_temp; } //����С����
		if(ny[1]==y_sec_num-1) ny[1] -= 1; //���߽��ϵĵ�����ǰһ����

		int nz[2] = { int((point[0].z-cell_geo.poi_min.z)/dz), int((point[1].z-cell_geo.poi_min.z)/dz) };
		if(nz[0]>nz[1]) { n_temp = nz[0]; nz[0] = nz[1]; nz[1] = n_temp; } //����С����	  
		if(nz[1]==z_sec_num-1) nz[1] -= 1; //���߽��ϵĵ�����ǰһ����

		for(int jz=nz[0]; jz<=nz[1]; jz++)
			for(int jy=ny[0]; jy<=ny[1]; jy++)
				for(int jx=nx[0]; jx<=nx[1]; jx++)
				{
					int elen = jz*(x_sec_num-1)*(y_sec_num-1) + jy*(x_sec_num-1) + jx;

					Point_3D Hexc[2];  //��Ԫ����Сֵ������ֵ��
					Hexc[0].x = nodes[elements[elen].nodes_id[0]].x;
					Hexc[0].y = nodes[elements[elen].nodes_id[0]].y;
					Hexc[0].z = nodes[elements[elen].nodes_id[0]].z;
					Hexc[1].x = nodes[elements[elen].nodes_id[6]].x;
					Hexc[1].y = nodes[elements[elen].nodes_id[6]].y;
					Hexc[1].z = nodes[elements[elen].nodes_id[6]].z;

					//---------------------------------------------------------------------------
					//�߶ε�ĳ���˵��ڵ�Ԫ�ڲ�
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
					//�߶εĶ˵㵥Ԫ�ⲿ�����߶��뵥Ԫ�ཻ
					double t_ratio[6];
					//�߶�����ֱ����X�ᴹֱƽ���ཻ
					t_ratio[0] = (Hexc[0].x - point[0].x)/(point[1].x - point[0].x);
					t_ratio[1] = (Hexc[1].x - point[0].x)/(point[1].x - point[0].x);
					//�߶�����ֱ����Y�ᴹֱƽ���ཻ
					t_ratio[2] = (Hexc[0].y - point[0].y)/(point[1].y - point[0].y);
					t_ratio[3] = (Hexc[1].y - point[0].y)/(point[1].y - point[0].y);
					//�߶�����ֱ����Z�ᴹֱƽ���ཻ
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
					else if(count>2) { hout <<"����ֱ���뵥Ԫ����Ľ����������2�������飡 ...Element_relative_cnts" << endl;  return 0; }
				}
	}

	return 1;
}
//---------------------------------------------------------------------------
//�����ⵥԪ�Ƿ�ƥ������ȷ��������׹��߶���Ϣ
int Mesher::Test_element_relative_cnts(vector<Point_3D> &cnps)const
{
	ofstream otec("Test_element_cnts.dat");
	otec << "TITLE = Test_element_cnts" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	const int count = 10;
	//������N�������嵥Ԫ
	srand((unsigned int)time(NULL));	//����������ɵ���ʼʱ��
	for(int i=0; i<count; i++)
	{
		int ele_num = rand()%(int)elements.size();
		//---------------------------------------------------------------------------
		//����������
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
		//������׹�����ά����
		vector<int> veci_temp;
		vector<vector<int> > veci;
		for(int j=0; j<(int)elements[ele_num].relative_cnts.size(); j++)
		{
			veci_temp.push_back(elements[ele_num].relative_cnts[j]);  //����ڵ�

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
//���Tecplot���ӻ���������
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
//���ݸ�����i j k�Լ����i_max j_max k_max�����������ڵ��λ�ã��ǵ㡢�߽��ߵ㡢�߽���㡢�ڵ㣩
//����ֵ��0�ڵ㣬1�߽���㣬2�߽��ߵ㣬3�ǵ�
//i--z���ꣻ j--y���ꣻk--x����
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
//һά�洢���ڱ߽�ڵ��ż���ʼ����Ϣ
int Mesher::Record_periodic_boundary_nodes_number(const struct RVE_Geo &cell_geo)
{
	//���ڼ�¼��Щ�ڵ��Ѿ�����¼��(0:û�м�¼����1:��¼��)
	vector<int> nflag(nodes.size(),0);

	//���ɸ��������������
	int k_max = int(cell_geo.len_x/dx);
	int j_max = int(cell_geo.wid_y/dy);
	int i_max = int(cell_geo.hei_z/dz);

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���ڱ߽�ڵ��ż������Ϣ
	int count = 0;  //��¼����ڵ���
	for( int i=0; i<=i_max; i++ )
		for( int j=0; j<=j_max; j++ )
			for( int k=0; k<=k_max; k++ )
			{
				if(nodes[count].type!=0&&nflag[count]==0)  //�����ڵ㲢��û�м�¼��
				{
					int key = 0;
					for(int m=i; m<=i_max; m+=i_max)
						for(int n=j; n<=j_max; n+=j_max)
							for(int p=k; p<=k_max; p+=k_max)
							{
								int num = m*(j_max+1)*(k_max+1) + n*(k_max+1) + p;
								peri_bnods[0].push_back(num);
								nflag[num]=1;	//�˵��Ѿ���¼��
								peri_bnods[1].push_back(key);  //���λ��
								key++;
							}
				}
				count++;  //����ڵ��ż�1
			}

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//��һά�洢����ص�Ӵ�С��������
	int left = 0;
	int right = 0;
	for(int i=0; i<(int)peri_bnods[1].size(); i++)
	{
		if(i+1==(int)peri_bnods[1].size()||peri_bnods[1][i+1]==0) 
		{
			right = i;
			//---------------------------------------------------------------
			//����С���������Ϊ�Ӵ�С����
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

	//������
	//hout << "peri_bnods:" << endl << endl;
	//for(int i=0; i<(int)peri_bnods[0].size(); i++)
	//{
	//	hout << peri_bnods[0][i] << "  " << peri_bnods[1][i] << endl;
	//}

	return 1;
}
//---------------------------------------------------------------------------
//���ɸ�������Ԥ���ֲ�ģ�͵�Чģ��
int Mesher::Generate_grids_for_effective_stiffness(ifstream &infile, const double &decayR, double &grid_vol)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//������ӳߴ�(������͸�)
	istringstream istr(Get_Line(infile));	
	istr >> dx >> dy >> dz;
	if(dx<0||dy<0||dz<0)
	{
		hout << "Ԥ���նȾ���ʱ�����ɸ��ӳ���ߵķָ�ϸ��������� ���飡 ...Generate_grids_for_effective_stiffness" << endl;
		return 0;
	}

	//���ӵ����
	grid_vol = dx*dy*dz;

	//���ɸ�������ĵ���
	int x_sec_num = 2*(int(decayR/dx-Zero)+2);  //��һ��СֵΪ�˱���ȡ���Ĵ���
	int y_sec_num = 2*(int(decayR/dy-Zero)+2);  //+2��һ��+1�ǵ�ĸ�����һ��ȶ�����1��һ��+1�Ǳ�֤decayR<=n*dx(y,z)ʱ��������Ҫ��1
	int z_sec_num = 2*(int(decayR/dz-Zero)+2);

	Point_3D origin(0,0,0);
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���ɽڵ�
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
	//���������嵥Ԫ

	//��������������
	for( int i=0; i<z_sec_num-1; i++ )
		for( int j=0; j<y_sec_num-1; j++ )
			for( int k=0; k<x_sec_num-1; k++ )
			{
				//������İ˸�����
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
//����һ����Ϣ��������ע���У���"%"��ͷ����
string Mesher::Get_Line(ifstream &infile)const
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
