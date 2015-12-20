//===========================================================================
// SolveEqu.cpp
// �����Է����麯��
// Member functions in a class to solving linear equations
//===========================================================================
#include "SolveEqu.h"

//---------------------------------------------------------------------------
//�����洢�նȾ���
int SolveEqu::izig(const vector<Node> &nodes, vector<long> &Iz, vector<int> &Ig)const
{
	//iz[i-1], iz[i]�м�¼���i������ص�������ig�д洢�Ŀ�ʼλ�úͽ���λ��
	//��һ����ı����С(û�б�����Ż�С�Ľڵ�)������iz���ݴӵڶ�����(i==1)��ʼ
	Iz.push_back(0);
    
	//�ӵڶ����ڵ㿪ʼѭ�����еĽڵ�
	for(int i=1; i<(int)nodes.size(); i++)
	{
		for(int j=0; j<(int)nodes[i].relative_nods.size(); j++)	//ע������ڵ����ؽڵ��Ǵ�С�������еģ����Ҳ������ڵ�����ı��
		{
			if(nodes[i].relative_nods[j]>=i) break;						//һ��һ����Ŵ��ˣ����ڵ������ʵ�����ڣ���˵������ı�Ŷ�����
			else Ig.push_back(nodes[i].relative_nods[j]);
		}
		Iz.push_back((long)Ig.size());
	}
    
	//������ڼ��
	//hout << "Iz:" << endl;
	//for(int i=0; i<(int)Iz.size(); i++)
	//	hout << i << " " << Iz[i] << endl;
	//hout << endl << "Ig:" << endl;
	//for(int i=0; i<(long)Ig.size(); i++)
	//	hout << i << " " << Ig[i] << endl;
    
	return 1;
}
//---------------------------------------------------------------------------
//�������ڱ߽�����Լ��(�Ҷ���������)(�㷨�ο��ư�����˵����ģ�����alph=1, Q=(q1,q2,q3);)
void SolveEqu::Periodical_boundary_constraints(const vector<Node> &nodes, const vector<int>* &peri_bnods, vector<long> &Iz, vector<int> &Ig,	vector<double> &AK,  vector<double>* F)const
{
	//-----------------------------------------------------------
	//��ʼ������
	vector<vector<double> > temp_AKii;
	vector<vector<vector<double> > > temp_AK(nodes.size(), temp_AKii); //���������temp_AKII��ʼ��temp_AK
	vector<double> tempakii(6,0);
	vector<double> tempak(9,0);
    
	for(int i=0; i<6; i++) tempakii[i] = AK[i]; //����Žڵ�Խ�����������
	temp_AKii.push_back(tempakii);
	for(int i=1; i<(int)nodes.size(); i++) //��1�Žڵ㿪ʼѭ��
	{
		long count=0;
		for(long j=Iz[i]-1; j>=Iz[i-1]; j--) //��ؽڵ�Ӵ�С����(�����Ժ�����ɾ��)
		{
			count = 6*(long)i+9*j;
			for(int k=0; k<9; k++) tempak[k] = AK[count+k];
			temp_AK[i].push_back(tempak);
		}
		count = 6*(long)i+9*Iz[i];
		for(int j=0; j<6; j++) tempakii[j] = AK[count+j];
		temp_AKii.push_back(tempakii);
	}
    
	//-----------------------------------------------------------
	//���Iz, Ig��AK
	Iz.clear();
	Ig.clear();
	AK.clear();
    
	//-----------------------------------------------------------
	//��ʼ����ؽڵ�����
	vector<vector<int> > relative_small(nodes.size());
	vector<vector<int> > relative_large(nodes.size());
	//ѭ�����еĽڵ�
	for(int i=0; i<(int)nodes.size(); i++)
	{
		for(int j=(int)nodes[i].relative_nods.size()-1; j>=0; j--)	//ע������ڵ����ؽڵ��Ǵ�С�������еģ����Ҳ������ڵ�����ı��, ��������Ҫȡ�Ӵ�С����
		{
			if(nodes[i].relative_nods[j]>i) relative_large[i].push_back(nodes[i].relative_nods[j]);				//����ǴӴ�С����(�����Ժ�����ɾ��)
			else if(nodes[i].relative_nods[j]<i) relative_small[i].push_back(nodes[i].relative_nods[j]);    //����ǴӴ�С����(�����Ժ�����ɾ��)
		}
	}
    
	//-----------------------------------------------------------
	//ѭ�����ڱ߽�ڵ�����
	int num[2] = { 0, peri_bnods[0][0] }; //num[0]��num[1]�ֱ��Ӧ������Ԫ��������Ӧ�á��е�xj��xi
    
	for(int i=0; i<(int)peri_bnods[0].size()-1; i++)
	{
		if(peri_bnods[1][i+1]==0) { num[1] = peri_bnods[0][i+1]; continue; } //һ��������������ʼ��
		num[0] = peri_bnods[0][i+1];
		//-----------------------------------------------------------
		//�����Ҷ�����
		double dist[3] = { nodes[num[1]].x-nodes[num[0]].x, nodes[num[1]].y-nodes[num[0]].y, nodes[num[1]].z-nodes[num[0]].z };		//����ڵ�xi��xj������
		deal_with_F(num, dist, relative_small, relative_large, temp_AKii, temp_AK, F);
        
		//-----------------------------------------------------------
		//�����նȾ���
		//�޸���i,j�Žڵ�����������ز���
		//-----------------------------------------------------------
		//����AK[ii], AK[ij]��AK[jj]�Ĳ���, ����ڵ�i��ڵ�jԭ�����, ����Ҫɾ��i��j�������, ���޸�relative_small��relative_large
		deal_with_AK_ii_ij_jj(num, relative_small, relative_large, temp_AKii, temp_AK);
        
		//-----------------------------------------------------------
		//����AK[ki]�Ĳ���
		//����ڵ�k (k!=i,j), ֻҪk��j���, ��AK[ki]+=AK[kj](����ԭ��AK[ki]�Ƿ�Ϊ��) //ע�����AK[kj]�Ƿ������ת��
		//���AK[ki]ԭ��Ϊ�����, ����Ҫ���Ӿ���, ���޸�relative_small��relative_large
		deal_with_AK_ki_kj(num, relative_small, relative_large, temp_AK);
	}
    
 	//-----------------------------------------------------------
	//��������Iz, Ig��AK
    
	//iz[i-1], iz[i]�м�¼���i������ص�������ig�д洢�Ŀ�ʼλ�úͽ���λ��
	//��һ����ı����С(û�б�����Ż�С�Ľڵ�)������iz���ݴӵڶ�����(i==1)��ʼ
	Iz.push_back(0);
	AK = temp_AKii[0];			//������0�Žڵ����������
    
	//�ӵڶ����ڵ㿪ʼѭ�����еĽڵ�
	for(int i=1; i<(int)nodes.size(); i++)
	{
		for(int j=(int)relative_small[i].size()-1; j>=0; j--)		//ע��������ؽڵ��ǴӴ�С���еģ� ����Ҫ�������
		{
			Ig.push_back(relative_small[i][j]);
		}
		Iz.push_back((long)Ig.size());
        
		for(int j=(int)temp_AK[i].size()-1; j>=0; j--)		//ע����������ֵ����relative_small�Ľڵ��������Ӧ�� ����ҲҪ�������
		{
			for(int k=0; k<9; k++) AK.push_back(temp_AK[i][j][k]);
		}
		for(int j=0; j<6; j++) AK.push_back(temp_AKii[i][j]);
	}
    
	//������ڼ��
	//hout << "Iz:" << endl;
	//for(int i=0; i<(int)Iz.size(); i++)
	//hout << i << " " << Iz[i] << endl;
    
	//hout << endl << "Ig:" << endl;
	//for(int i=0; i<(long)Ig.size(); i++)
	//hout << i << " " << Ig[i] << endl;
    
	//hout << endl << "AK:" << endl;
	//for(int i=0; i<(long)AK.size(); i++)
	//hout << i << " " << AK[i] << endl;
    
	//hout << endl << "F:" << endl;
	//for(int i=0; i<(int)F.size(); i++)
	//hout << i << " " << F[i] << endl;
}
//-----------------------------------------------------------
//�����Ҷ�����F
void SolveEqu::deal_with_F(const int num[], const double dist[3], const vector<vector<int> > &relative_small, const vector<vector<int> > &relative_large,
                           const vector<vector<double> > &temp_AKii, const vector<vector<vector<double> > > &temp_AK, vector<double>* F)const
{
	double udiff[9][3] = { {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0} };
	for(int i=0; i<9; i++)
	{
		//-----------------------------------------------------------
		//���þ��Ȼ�Ӧ����Ϊ���ڱ߽�����
		double E[3][3] = {{0}, {0}, {0}}; //�������ڱ߽���������ʱ�ı��ξ�������
		switch(i)
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
        
		//-----------------------------------------------------------
		//��������Q (xj=alpha*xi - Q; alpha =1������Ԫ��������Ӧ�á�), dist��xi��xj������, ����udiffi�൱��Q
		for(int j=0; j<3; j++)
			for(int k=0; k<3; k++)
				udiff[i][j] += E[j][k]*dist[k];
	}
    
	//-----------------------------------------------------------
	//�ڵ���k<�ڵ���num[0]
	for(int j=0; j<(int)relative_small[num[0]].size(); j++)
	{
		for(int k=0; k<=2; k++)
		{
			const int rsn = 3*relative_small[num[0]][j]+k;
			for(int m=0; m<=2; m++)
			{
				const int km = k+3*m;
				for(int i=0; i<9; i++)
					F[i][rsn] += temp_AK[num[0]][j][km]*udiff[i][m];  //ע���ǳ˾����ת��
			}
		}
	}
    
	//-----------------------------------------------------------
	//�ڵ���k>�ڵ���num[0]
	for(int j=0; j<(int)relative_large[num[0]].size(); j++)
	{
		//���ַ�����j��k��relative_small�е�λ��
		int nk = relative_large[num[0]][j];
		if(nk==num[1]) continue; //k���ܵ���num[1], ����������������濼��
		bool mark = false;
		int left = 0;
		int middle = 0;
		int right = (int)relative_small[nk].size()-1;
		while(right>=left)
		{
			middle = (left + right)/2;
			if(relative_small[nk][middle] == num[0]) { mark = true; break; }				//�ҵ����λ��
			else if(relative_small[nk][middle] < num[0]) right = middle - 1;
			else left = middle + 1;
		}
		if(!mark) hout << "�����ڽڵ�" << nk << "�����С��Žڵ��в�û���ҵ��ڵ�" << num[0] << "�����飡" << endl;
		for(int k=0; k<=2; k++)
		{
			const int tnk = 3*nk+k;
			for(int m=0; m<=2; m++)
			{
				const int km = 3*k+m;
				for(int i=0; i<9; i++)
					F[i][tnk] += temp_AK[nk][middle][km]*udiff[i][m];  //ע���ǳ˾�����
			}
		}
	}
    
	//-----------------------------------------------------------
	//�ڵ���num[1]
	
	//���ַ���i����ؽڵ��в���j��λ�ã��������ؾ��������
	bool mark = false;
	int nm = 0;
	if(num[0]<=relative_small[num[1]][0]&&num[0]>=relative_small[num[1]].back())
	{
		int left = 0;
		int middle = 0;
		int right = (int)relative_small[num[1]].size()-1;
		while(right>=left)
		{
			middle = (left + right)/2;
			if(relative_small[num[1]][middle] == num[0]) { mark = true; nm = middle; break; }				//�ҵ����λ��
			else if(relative_small[num[1]][middle] < num[0]) right = middle - 1;
			else left = middle + 1;
		}
	}
    
	const int tnum0 = 3*num[0];
	const int tnum1 = 3*num[1];
	for(int i=0; i<9; i++)
	{
		if(mark) //num[1]��num[0]���
		{
			for(int j=0; j<=2; j++)
			{
				F[i][tnum1+j] += F[i][tnum0+j];
				for(int k=0; k<=2; k++)
				{
					if(j==2||k==2)
						F[i][tnum1+j] += (temp_AKii[num[0]][j+k+1]+temp_AK[num[1]][nm][3*j+k])*udiff[i][k];
					else
						F[i][tnum1+j] += (temp_AKii[num[0]][j+k]+temp_AK[num[1]][nm][3*j+k])*udiff[i][k];
				}
			}
		}
		else //num[1]��num[0]�����, AK[ij]=0
		{
			for(int j=0; j<=2; j++)
			{
				F[i][tnum1+j] += F[i][tnum0+j];
				for(int k=0; k<=2; k++)
				{
					if(j==2||k==2)
						F[i][tnum1+j] += temp_AKii[num[0]][j+k+1]*udiff[i][k];
					else
						F[i][tnum1+j] += temp_AKii[num[0]][j+k]*udiff[i][k];
				}
			}
		}
		//-----------------------------------------------------------
		//�ڵ���num[0]
		for(int j=0; j<=2; j++) F[i][tnum0+j] = 0;
	}
}
//-----------------------------------------------------------
//����AK[ii], AK[ij]��AK[jj]�Ĳ���
void SolveEqu::deal_with_AK_ii_ij_jj(const int num[], vector<vector<int> > &relative_small, vector<vector<int> > &relative_large, vector<vector<double> > &temp_AKii, vector<vector<vector<double> > > &temp_AK)const
{
	//��������㱾�����(��num[0]��num[1]����ؽڵ�����relative_small��, num[1]Ҳ��num[0]����ؽڵ�����relative_large��)��
	//��ô��num[0]��num[1]����ؽڵ�������ɾ��(�޸�relative_small), �ֽ�num[1]��num[0]����ؽڵ�������ɾ��(�޸�relative_large)
	
	bool mark[2] = { false, false };		//��¼��ػ��߲����
	int nodij = 0;										//��¼j����i����ؽڵ��е�λ��
    
	//-----------------------------------------------------------
	//��num[1]��relative_small�в��Ҳ�ɾ��num[0]���˴�num[1]�����ܵ���0����Ϊ������num[0]��num[1]С��
	if(num[0]<=relative_small[num[1]][0]&&num[0]>=relative_small[num[1]].back())
	{
		//���ַ����Ҳ�ɾ��
		int left = 0;
		int middle = 0;
		int right = (int)relative_small[num[1]].size()-1;
		while(right>=left)
		{
			middle = (left + right)/2;
			if(relative_small[num[1]][middle] == num[0]) { mark[0]=true; nodij =  middle; break; } //�ڵ�����ͬ�����
			else if(relative_small[num[1]][middle] < num[0]) right = middle - 1;
			else left = middle + 1;
		}
		if(mark[0])	//ɾ��
		{
			for(int i=middle; i<(int)relative_small[num[1]].size()-1; i++) relative_small[num[1]][i] = relative_small[num[1]][i+1];
			relative_small[num[1]].pop_back();	//ɾ�����һ��Ԫ��
		}
	}
    
	//-----------------------------------------------------------
	//��num[0]��relative_large�в��Ҳ�ɾ��num[1]���˴�num[0]�����ܵ������һ���ڵ�ı��, ��Ϊ������num[1]��num[0]��
	if(num[1]<=relative_large[num[0]][0]&&num[1]>=relative_large[num[0]].back())
	{
		//���ַ����Ҳ�ɾ��
		int left = 0;
		int middle = 0;
		int right = (int)relative_large[num[0]].size()-1;
		while(right>=left)
		{
			middle = (left + right)/2;
			if(relative_large[num[0]][middle] == num[1]) { mark[1]=true; break; } //�ڵ�����ͬ�����
			else if(relative_large[num[0]][middle] < num[1]) right = middle - 1;
			else left = middle + 1;
		}
		if(mark[0]!=mark[1]) hout << "���󣡽ڵ�" << num[0] << "�ͽڵ�" << num[1] << "֮���relative_small��relative_large��Ӧ��ϵ�������飡" << endl;
		if(mark[1])	//ɾ��
		{
			for(int i=middle; i<(int)relative_large[num[0]].size()-1; i++) relative_large[num[0]][i] = relative_large[num[0]][i+1];
			relative_large[num[0]].pop_back();	//ɾ�����һ��Ԫ��
		}
	}
    
	//-----------------------------------------------------------
	//�޸�AK[ii]��AK[ij]
	if(mark[0]) //�ڵ�i��ڵ�j���
	{
		//-----------------------------------------------------------
		//�޸��ܸ���AK[ii], ����AK[ii]��AK[ij]�Ĺ�����, ����Ҫ�����޸�AK[ii]
		temp_AKii[num[1]][0] += temp_AKii[num[0]][0] + 2*temp_AK[num[1]][nodij][0];
		temp_AKii[num[1]][1] += temp_AKii[num[0]][1] + temp_AK[num[1]][nodij][1]+temp_AK[num[1]][nodij][3];
		temp_AKii[num[1]][2] += temp_AKii[num[0]][2] + 2*temp_AK[num[1]][nodij][4];
		temp_AKii[num[1]][3] += temp_AKii[num[0]][3] + temp_AK[num[1]][nodij][2]+temp_AK[num[1]][nodij][6];
		temp_AKii[num[1]][4] += temp_AKii[num[0]][4] + temp_AK[num[1]][nodij][5]+temp_AK[num[1]][nodij][7];
		temp_AKii[num[1]][5] += temp_AKii[num[0]][5] + 2*temp_AK[num[1]][nodij][8];
        
		//-----------------------------------------------------------
		//ɾ��AK[ij]
		for(int i=nodij; i<(int)temp_AK[num[1]].size()-1; i++) temp_AK[num[1]][i] = temp_AK[num[1]][i+1];
		temp_AK[num[1]].pop_back();	//ɾ�����һ��Ԫ��
	}
	else	//�ڵ�i��ڵ�j�����
	{
		//-----------------------------------------------------------
		//�޸��ܸ���AK[ii]
		temp_AKii[num[1]][0] += temp_AKii[num[0]][0];
		temp_AKii[num[1]][1] += temp_AKii[num[0]][1];
		temp_AKii[num[1]][2] += temp_AKii[num[0]][2];
		temp_AKii[num[1]][3] += temp_AKii[num[0]][3];
		temp_AKii[num[1]][4] += temp_AKii[num[0]][4];
		temp_AKii[num[1]][5] += temp_AKii[num[0]][5];
	}
    
	//-----------------------------------------------------------
	//�޸��ܸ���AK[jj]
	temp_AKii[num[0]][0] = temp_AKii[num[0]][2] = temp_AKii[num[0]][5] = 1;
	temp_AKii[num[0]][1] = temp_AKii[num[0]][3] = temp_AKii[num[0]][4] = 0;
    
}
//-----------------------------------------------------------
//����AK[ki]��AK[kj]�Ĳ���
void SolveEqu::deal_with_AK_ki_kj(const int num[], vector<vector<int> > &relative_small, vector<vector<int> > &relative_large, vector<vector<vector<double> > > &temp_AK)const
{
	//-----------------------------------------------------------
	//ѭ��relative_large[num[0]]
	for(int i=0; i<(int)relative_large[num[0]].size(); i++)
	{
		int nodkj = 0; //��¼j��k����ؽڵ��е�λ�ñ��
		//--------------------------------------------------------------------------------------------------------------------------------------------
		int nk = relative_large[num[0]][i];
		//��relative_small[nk]��Ѱ��num[0]��λ��
		if(num[0]<=relative_small[nk][0]&&num[0]>=relative_small[nk].back())
		{
			//���ַ����Ҳ�ɾ��
			bool mark = false;
			int left = 0;
			int middle = 0;
			int right = (int)relative_small[nk].size()-1;
			while(right>=left)
			{
				middle = (left + right)/2;
				if(relative_small[nk][middle] == num[0]) { mark = true; nodkj =  middle; break; } //�ڵ�����ͬ�����
				else if(relative_small[nk][middle] < num[0]) right = middle - 1;
				else left = middle + 1;
			}
			if(!mark) hout << "���󣡽ڵ�" << nk << "��relative_small��ؽڵ���û�ҵ�" << num[0] << "�����飡" << endl;
		}
		else hout << "���󣡽ڵ�" << nk << "��relative_small��ؽڵ��в�����" << num[0] << "�����飡" << endl;
        
		//------------------------------------------------------------------------------------------------------------------------------------------------
		//��Ѱnk��num[1]�Ĺ�ϵ���޸�
		if(nk == num[1]) hout << "���󣡽ڵ�" << num[0] << "����ؽڵ������нڵ�" << num[1] << "���ڣ����飡" << endl;
		else if(nk > num[1])
		{
			bool mark[2] = {false, false};
			//��nk��relative_small�в���num[1], ���޸�temp_AK, �ڲ���num[1]��ص�ͬʱɾ��num[0]
			if(num[1] > relative_small[nk][0])
			{
				for(int j=nodkj; j>0; j--) relative_small[nk][j] = relative_small[nk][j-1];
				relative_small[nk][0] = num[1];
                
				vector<double> tempak(temp_AK[nk][nodkj]);
				for(int j=nodkj; j>0; j--) temp_AK[nk][j] = temp_AK[nk][j-1];
				temp_AK[nk][0] = tempak;
			}
			else if(num[1]<relative_small[nk].back())
			{
				for(int i=nodkj; i<(int)relative_small[nk].size()-1; i++) relative_small[nk][i] = relative_small[nk][i+1];
				relative_small[nk].back() = num[1];
                
				vector<double> tempak(temp_AK[nk][nodkj]);
				for(int i=nodkj; i<(int)temp_AK[nk].size()-1; i++) temp_AK[nk][i] = temp_AK[nk][i+1];
				temp_AK[nk].back() = tempak;
			}
			else //num[1]��relative_small[nk]֮��
			{
				//���ַ�����
				int left = 0;
				int middle = 0;
				int right = (int)relative_small[nk].size()-1;
				while(right>=left)
				{
					middle = (left + right)/2;
					if(relative_small[nk][middle] == num[1]) { mark[0] = true; break; } //�ڵ�����ͬ�����
					else if(relative_small[nk][middle] < num[1]) right = middle - 1;
					else left = middle + 1;
				}
				if(mark[0]) //num[1]�Ѿ���nk����ؽڵ�
				{
					for(int i=nodkj; i<(int)relative_small[nk].size()-1; i++)	relative_small[nk][i] = relative_small[nk][i+1];
					relative_small[nk].pop_back();	//ɾ�����һ��Ԫ��
                    
					for(int j=0; j<9; j++)	temp_AK[nk][middle][j] += temp_AK[nk][nodkj][j];  //����
                    
					for(int i=nodkj; i<(int)temp_AK[nk].size()-1; i++) temp_AK[nk][i] = temp_AK[nk][i+1];
					temp_AK[nk].pop_back();	//ɾ�����һ��Ԫ��
				}
				else
				{
					if(nodkj<left) hout << "���󣡽ڵ�" << num[0] << "��" << num[1] << "���ڽڵ�" << nk << "��relative_small�ڣ�������λ�ó������飡" << endl;
					for(int j=nodkj; j>left; j--) relative_small[nk][j] = relative_small[nk][j-1];
					relative_small[nk][left] = num[1];
                    
					vector<double> tempak(temp_AK[nk][nodkj]);
					for(int j=nodkj; j>left; j--) temp_AK[nk][j] = temp_AK[nk][j-1];
					temp_AK[nk][left] = tempak;
				}
			}
            
			//��num[1]��relative_large�в���nk���޸�
			if(nk>relative_large[num[1]][0])
			{
				relative_large[num[1]].push_back(0);
				for(int j=(int)relative_large[num[1]].size()-1; j>0; j--) relative_large[num[1]][j] = relative_large[num[1]][j-1];
				relative_large[num[1]][0] = nk;
			}
			else if(nk<relative_large[num[1]].back())	relative_large[num[1]].push_back(nk);
			else //nk��relative_large[num[1]]֮��
			{
				//���ַ�����
				int left = 0;
				int middle = 0;
				int right = (int)relative_large[num[1]].size()-1;
				while(right>=left)
				{
					middle = (left + right)/2;
					if(relative_large[num[1]][middle] == nk) { mark[1] = true; break; } //�ڵ�����ͬ�����
					else if(relative_large[num[1]][middle] < nk) right = middle - 1;
					else left = middle + 1;
				}
				if(mark[0]!=mark[1]) hout << "���󣡽ڵ�" << num[0] << "�ͽڵ�" << nk << "������Բ�һ�£����飡" << endl;
				if(!mark[1]) //nk����num[1]����ؽڵ�
				{
					relative_large[num[1]].push_back(0);
					for(int j=(int)relative_large[num[1]].size()-1; j>left; j--) relative_large[num[1]][j] = relative_large[num[1]][j-1];
					relative_large[num[1]][left] = nk;
				}
			}
		}
		else  //nk < num[1]
		{
			bool mark[2] = {false, false};
			//��num[1]��relative_small�в���nk, ���޸�temp_AK
			if(nk > relative_small[num[1]][0])
			{
				relative_small[num[1]].push_back(0);
				for(int j=(int)relative_small[num[1]].size()-1; j>0; j--) relative_small[num[1]][j] = relative_small[num[1]][j-1];
				relative_small[num[1]][0] = nk;
				
				vector<double> tempak(9,0);
				temp_AK[num[1]].push_back(tempak);
				for(int j=(int)temp_AK[num[1]].size()-1; j>0; j--) temp_AK[num[1]][j] = temp_AK[num[1]][j-1];
				for(int j=0; j<=2; j++)											//��j<k<iʱҪ����AK[kj]��ת�þ���
					for(int m=0; m<=2; m++)
						temp_AK[num[1]][0][3*j+m] = temp_AK[nk][nodkj][j+3*m];
			}
			else if(nk < relative_small[num[1]].back())
			{
				relative_small[num[1]].push_back(nk);
                
				vector<double> tempak(9,0);
				temp_AK[num[1]].push_back(tempak);
				for(int j=0; j<=2; j++)											//��j<k<iʱҪ����AK[kj]��ת�þ���
					for(int m=0; m<=2; m++)
						temp_AK[num[1]].back()[3*j+m] = temp_AK[nk][nodkj][j+3*m];
			}
			else //nk��relative_small[num[1]]֮��
			{
				//���ַ�����
				int left = 0;
				int middle = 0;
				int right = (int)relative_small[num[1]].size()-1;
				while(right>=left)
				{
					middle = (left + right)/2;
					if(relative_small[num[1]][middle] == nk) { mark[0] = true; break; } //�ڵ�����ͬ�����
					else if(relative_small[num[1]][middle] < nk) right = middle - 1;
					else left = middle + 1;
				}
				if(mark[0]) //nk�Ѿ���num[1]����ؽڵ�
				{
					for(int j=0; j<=2; j++)											//��j<k<iʱҪ����AK[kj]��ת�þ���
						for(int m=0; m<=2; m++)
							temp_AK[num[1]][middle][3*j+m] += temp_AK[nk][nodkj][j+3*m];
				}
				else
				{
					relative_small[num[1]].push_back(0);
					for(int j=(int)relative_small[num[1]].size()-1; j>left; j--) relative_small[num[1]][j] = relative_small[num[1]][j-1];
					relative_small[num[1]][left] = nk;
                    
					vector<double> tempak(9,0);
					temp_AK[num[1]].push_back(tempak);
					for(int j=(int)temp_AK[num[1]].size()-1; j>left; j--) temp_AK[num[1]][j] = temp_AK[num[1]][j-1];
					for(int j=0; j<=2; j++)											//��j<k<iʱҪ����AK[kj]��ת�þ���
						for(int m=0; m<=2; m++)
							temp_AK[num[1]][left][3*j+m] = temp_AK[nk][nodkj][j+3*m];
				}
			}
            
			//��nk��relative_large�в���num[1]���޸�
			if(num[1] > relative_large[nk][0])
			{
				relative_large[nk].push_back(0);
				for(int j=(int)relative_large[nk].size()-1; j>0; j--) relative_large[nk][j] = relative_large[nk][j-1];
				relative_large[nk][0] = num[1];
			}
			else if(num[1] < relative_large[nk].back())	relative_large[nk].push_back(num[1]);
			else //num[1]��relative_large[nk]֮��
			{
				//���ַ�����
				int left = 0;
				int middle = 0;
				int right = (int)relative_large[nk].size()-1;
				while(right>=left)
				{
					middle = (left + right)/2;
					if(relative_large[nk][middle] == num[1]) { mark[1] = true; break; } //�ڵ�����ͬ�����
					else if(relative_large[nk][middle] < num[1]) right = middle - 1;
					else left = middle + 1;
				}
				if(mark[0]!=mark[1]) hout << "���󣡽ڵ�" << num[0] << "�ͽڵ�" << nk << "������Բ�һ�£����飡" << endl;
				if(!mark[1]) //num[1]����nk����ؽڵ�,
				{
					relative_large[nk].push_back(0);
					for(int j=(int)relative_large[nk].size()-1; j>left; j--) relative_large[nk][j] = relative_large[nk][j-1];
					relative_large[nk][left] = num[1];
				}
			}
            
			//ɾ��relative_small[nk]�е�num[0]
			for(int j=nodkj; j<(int)relative_small[nk].size()-1; j++) relative_small[nk][j] = relative_small[nk][j+1];
			relative_small[nk].pop_back();
            
			//ɾ��temp_AK[nk]�е�num[0]����ؾ���
			for(int j=nodkj; j<(int)temp_AK[nk].size()-1; j++) temp_AK[nk][j] = temp_AK[nk][j+1];
			temp_AK[nk].pop_back();
		}
	}
	//���relative_large[num[0]]
	relative_large[num[0]].clear();
    
	//-----------------------------------------------------------
	//ѭ��relative_small[num[0]], �����һ��Ҳ��num[1]С
	for(int i=0; i<(int)relative_small[num[0]].size(); i++)
	{
		int nk = relative_small[num[0]][i];
        
		bool mark[2] = {false, false};
		//��num[1]��relative_small�в���nk, ���޸�temp_AK
		if(nk > relative_small[num[1]][0])
		{
			relative_small[num[1]].push_back(0);
			for(int j=(int)relative_small[num[1]].size()-1; j>0; j--) relative_small[num[1]][j] = relative_small[num[1]][j-1];
			relative_small[num[1]][0] = nk;
            
			vector<double> tempak(9,0);
			temp_AK[num[1]].push_back(tempak);
			for(int j=(int)temp_AK[num[1]].size()-1; j>0; j--) temp_AK[num[1]][j] = temp_AK[num[1]][j-1];
			temp_AK[num[1]][0] = temp_AK[num[0]][i];
		}
		else if(nk < relative_small[num[1]].back())
		{
			relative_small[num[1]].push_back(nk);
            
			temp_AK[num[1]].push_back(temp_AK[num[0]][i]);
		}
		else //nk��relative_small[num[1]]֮��
		{
			//���ַ�����
			int left = 0;
			int middle = 0;
			int right = (int)relative_small[num[1]].size()-1;
			while(right>=left)
			{
				middle = (left + right)/2;
				if(relative_small[num[1]][middle] == nk) { mark[0] = true; break; } //�ڵ�����ͬ�����
				else if(relative_small[num[1]][middle] < nk) right = middle - 1;
				else left = middle + 1;
			}
			if(mark[0]) //nk�Ѿ���num[1]����ؽڵ�
			{
				for(int j=0; j<9; j++)	temp_AK[num[1]][middle][j] += temp_AK[num[0]][i][j];
			}
			else
			{
				relative_small[num[1]].push_back(0);
				for(int j=(int)relative_small[num[1]].size()-1; j>left; j--) relative_small[num[1]][j] = relative_small[num[1]][j-1];
				relative_small[num[1]][left] = nk;
                
				vector<double> tempak(9,0);
				temp_AK[num[1]].push_back(tempak);
				for(int j=(int)temp_AK[num[1]].size()-1; j>left; j--) temp_AK[num[1]][j] = temp_AK[num[1]][j-1];
				temp_AK[num[1]][left] = temp_AK[num[0]][i];
			}
		}
        
		int nodkj = 0; //��¼j��k����ؽڵ��е�λ�ñ��
		//��relative_large[nk]��Ѱ��num[0]��λ��
		if(num[0]<=relative_large[nk][0]&&num[0]>=relative_large[nk].back())
		{
			//���ַ����Ҳ�ɾ��
			bool mark = false;
			int left = 0;
			int middle = 0;
			int right = (int)relative_large[nk].size()-1;
			while(right>=left)
			{
				middle = (left + right)/2;
				if(relative_large[nk][middle] == num[0]) { mark = true; nodkj =  middle; break; } //�ڵ�����ͬ�����
				else if(relative_large[nk][middle] < num[0]) right = middle - 1;
				else left = middle + 1;
			}
			if(!mark) hout << "���󣡽ڵ�" << nk << "��relative_large��ؽڵ���û�ҵ�" << num[0] << "�����飡" << endl;
		}
		else hout << "���󣡽ڵ�" << nk << "��relative_large��ؽڵ��в�����" << num[0] << "�����飡" << endl;
        
        
		//��nk��relative_large��Ѱ��num[1]������num[1], ɾ��num[0]
		if(num[1]>relative_large[nk][0])
		{
			for(int j=nodkj; j>0; j--) relative_large[nk][j] = relative_large[nk][j-1];
			relative_large[nk][0] = num[1];
		}
		else //num[1]��relative_large[nk]��0��nodkj֮���ң�num[1]������С��relative_large[nk].back(), ��Ϊ����num[0]<num[1]
		{
			//���ַ�����
			int left = 0;
			int middle = 0;
			int right = nodkj-1;
			while(right>=left)
			{
				middle = (left + right)/2;
				if(relative_large[nk][middle] == num[1]) { mark[1] = true; break; } //�ڵ�����ͬ�����
				else if(relative_large[nk][middle] < num[1]) right = middle - 1;
				else left = middle + 1;
			}
			if(mark[0]!=mark[1]) hout << "���󣡽ڵ�" << num[1] << "�ͽڵ�" << nk << "������Բ�һ�£����飡" << endl;
			if(!mark[1])  //num[1]����nk����ؽڵ�
			{
				for(int j=nodkj; j>left; j--) relative_large[nk][j] = relative_large[nk][j-1];
				relative_large[nk][left] = num[1];
			}
			else //num[1]��nk����ؽڵ㣬ɾ��num[0]
			{
				for(int j=nodkj; j<(int)relative_large[nk].size()-1; j++) relative_large[nk][j] = relative_large[nk][j+1];
				relative_large[nk].pop_back();
			}
		}
	}
	//���relative_small[num[0]]
	relative_small[num[0]].clear();
	//���temp_AK[num[0]]
	temp_AK[num[0]].clear();
}
//---------------------------------------------------------------------------
//���ù̶�λ��Լ��
int SolveEqu::Fixed_displacement_constraints(ifstream &infile, const vector<Node> &nodes, int &bnod_num, vector<int> &ip, vector<double> &vp)const
{
	//-----------------------------------------------------------
	//����λ������
	int fixed_displacement_num;
	istringstream istn(Get_Line(infile));
	istn >> fixed_displacement_num;		//����λ������
	while(fixed_displacement_num)
	{
		//Լ����������һ
		fixed_displacement_num--;
		//����Լ������
		istringstream istr(Get_Line(infile));
		//Լ����������
		Point_3D fpoi;
		Plane_3D fpla;
		string Boundary_case;
		istr >> Boundary_case;
		if(Boundary_case=="Point") istr >> fpoi.x >> fpoi.y >> fpoi.z;
		else if(Boundary_case=="Surface") istr >> fpla.coef[0] >> fpla.coef[1] >> fpla.coef[2] >> fpla.coef[3];
		else { hout << "ע�⣡Լ���������ǵ�����棬����������λ��Լ�����ݣ�" << endl; return 0; }
		//Լ��ֵ
		vector<double> value;
		int sign = 0;
		double val_temp;
		while(!istr.eof())
		{
			string str_temp;
			istr >> str_temp;
			if(str_temp.empty()) continue;
			if(str_temp=="Fixed_displacement_x")
			{
				if(sign>=1) { hout << "ע�⣡λ��Լ�������б�ʶ�ظ������������˳��ߵ�������������Լ�����ݣ�" << endl; return 0; }
				sign += 1;
				istr >> val_temp;		//��������
				value.push_back(val_temp);
			}
			else if(str_temp=="Fixed_displacement_y")
			{
				if(sign>=2) { hout << "ע�⣡λ��Լ�������б�ʶ�ظ������������˳��ߵ�������������Լ�����ݣ�" << endl; return 0; }
				sign += 2;
				istr >> val_temp;		//��������
				value.push_back(val_temp);
			}
			else if(str_temp=="Fixed_displacement_z")
			{
				if(sign>=4) { hout << "ע�⣡λ��Լ�������б�ʶ�ظ������������˳��ߵ�������������Լ�����ݣ�" << endl; return 0; }
				sign += 4;
				istr >> val_temp;		//��������
				value.push_back(val_temp);
			}
			else { hout << "ע�⣡Լ�������й̶�λ��������ʶ�д������Ǵ���ĩβ�пո�����������λ��Լ�����ݣ�" << endl; return 0; }
		}
		if((int)value.size()==0||(int)value.size()>3)
		{
			hout << "ע�⣡Լ�������й̶�λ������û���������������࣬����������Լ�����ݣ�" << endl;
			return 0;
		}
        
		//����λ��Լ������
		int count = 0;
		for(int i=0; i<(int)nodes.size(); i++)
		{
			if((Boundary_case=="Point"&&fpoi.distance_to(nodes[i].x, nodes[i].y, nodes[i].z)==0.0)||			//���̶���ľ������0
               (Boundary_case=="Surface"&&fpla.contain(nodes[i].x, nodes[i].y, nodes[i].z)==1))					//����һ���ռ�ƽ����
			{
				count ++;	//��¼�ҵ�һ����
				ip.push_back(i);
				ip.push_back(sign);
				//��ʼ���߽�λ��Լ������Ϣ
				//--------------------------------------------------------------------------------------------------------
				//Type	u		v		w		�ձ�ʾ��һά�����ɱ߽�(�ڴ˴�ͳһ����ֵ�����ǲ��ᱻ����)
				//   7		0		0		0		����0����ʾ����ֵ�����ֵ ��������ֵ
				//	 6		��	0		0
				//   5		0		��	0
				//	 4		��	��	0
				//   3		0		0		��
				//   2		��	0		��
				//   1		0		��	��
				//--------------------------------------------------------------------------------------------------------
				switch(sign)
				{
					case 1:	vp.push_back(value[0]);	vp.push_back(0);				vp.push_back(0);	break;
					case 2:	vp.push_back(0);				vp.push_back(value[0]);	vp.push_back(0);	break;
					case 3:	vp.push_back(value[0]);	vp.push_back(value[1]);	vp.push_back(0);	break;
					case 4:	vp.push_back(0);				vp.push_back(0);				vp.push_back(value[0]);	break;
					case 5:	vp.push_back(value[0]);	vp.push_back(0);				vp.push_back(value[1]);	break;
					case 6:	vp.push_back(0);				vp.push_back(value[0]);	vp.push_back(value[1]);	break;
					case 7:	vp.push_back(value[0]);	vp.push_back(value[1]);	vp.push_back(value[2]);	break;
					default:	hout << "ע�⣡λ��Լ�������ͱ����������������λ��Լ��������" << endl;	return 0;
				}
				bnod_num++;  //�߽�ڵ������һ
			}
		}
		if(count==0) { hout << "ע�⣡λ��Լ����λ���������󣨿��ܲ������κ�����ڵ㣩������������λ��Լ��������" << endl; return 0; }
	}
    
	return 1;
}
//---------------------------------------------------------------------------
//������ά�̶�λ��Լ�������λ��ֵ
void SolveEqu::Deal_with_displacement_zero_value(const int &bnod_num, const int &nod_size, const vector<long> &Iz, const vector<int> &Ig, const vector<int> &ip,
                                                 const vector<double> &vp, vector <double>* F, vector <double> &AK)const
//--------------------------------------------------------------------------------------------------------
//Type	u		v		w		�ձ�ʾ��һά�����ɱ߽�(�ڴ˴�ͳһ����ֵ�����ǲ��ᱻ����)
//   7		0		0		0		����0����ʾ����ֵ�����ֵ ��������ֵ
//	 6		��	0		0
//   5		0		��	0
//	 4		��	��	0
//   3		0		0		��
//   2		��	0		��
//   1		0		��	��
//--------------------------------------------------------------------------------------------------------
{
	for(int i=0; i<bnod_num; i++)
	{
		int key[3] = {1, 2, 4};
		int Ia = 2*i;
		int Ih = 3*i;
		int Mt = ip[Ia];
		int Lkt = ip[Ia+1];
		for(int j=2; j>=0; j--)	 //����һС��ѭ������Ҫ��key[0]~key[2]�зֱ��¼u,v,w������ֵ(key==1)������Ϊ��(key==0)
		{
			int item = Lkt%key[j];	//ȡ����
			key[j] = Lkt/key[j];		//ȡ��
			Lkt = item;
		}
		for(int j=0; j<=2; j++)
		{
			if(key[j]==1&&fabs(vp[Ih+j])<Zero)
			{
				long Kmt = 0;
				//---------------------------------------------------------------------
				for(int k=0; k<9; k++)	F[k][3*Mt+j]=0.0;		//�Ҷ���
				//---------------------------------------------------------------------
				if(Mt!=0)				//Mt����ؽڵ��б�Mt���С��
				{
					for(long k=Iz[Mt-1]; k<Iz[Mt]; k++)
						for(int m=0; m<=2; m++)
						{
							Kmt = 6*(long)Mt+9*k+3*j+m;
							AK[Kmt]=0.0;
						}
				}
				//---------------------------------------------------------------------
				for(int k=Mt+1; k<nod_size; k++)		//Mt����ؽڵ��б�Mt��Ŵ��
				{
					if(Mt>=Ig[Iz[k-1]]&&Mt<=Ig[Iz[k]-1])
					{
						//���ַ�����
						bool mark = false;
						long left = Iz[k-1];
						long middle = 0;
						long right = Iz[k]-1;
						while(right>=left)
						{
							middle = (left + right)/2;
							if(Ig[middle] == Mt) { mark = true; break; } //�ڵ�����ͬ�����
							else if(Ig[middle] > Mt) right = middle - 1;
							else left = middle + 1;
						}
						if(mark)  //�ҵ���Ӧ��
						{
							for(int n=0; n<=6; n+=3)
							{
								Kmt = 6*(long)k+9*middle+j+n;
								AK[Kmt]=0.0;
							}
						}
					}
				}
				//---------------------------------------------------------------------
				for(int k=0; k<=2; k++)  //����ԽǾ���
				{
					if(j==2||k==2) Kmt=6*(long)Mt+9*Iz[Mt]+j+k+1;				//�൱��	[0,1,	3]		��j==0ʱ, �� 0,	1,	3;
					else Kmt=6*(long)Mt+9*Iz[Mt]+j+k;									//				[1,2,	4]		��j==1ʱ, ��	1,	2,	4;
 					if(k==j) AK[Kmt]=1.0;																//				[3,4,	5]		�� j==2ʱ �� 3,	4,	5;
					else AK[Kmt]=0.0;
				}
			}
		}
	}
}
//---------------------------------------------------------------------------
//������Է����麯��(�����ݶ�CONJUGATE GRADIENT METHOD)
void SolveEqu::Solve_linear_equations(const int &bnod_num, const int &N, const vector<long> &Iz, const vector<int> &Ig, const vector<int> &ip,
                                      const vector<double> &vp, const vector<double> &A, const vector<double> &B, vector<double> &X)const
{
	//---------------------------------------------------------------------
	vector<double> P(3*N, 0);
	vector<double> R(3*N, 0);
	vector<double> S(3*N, 0);
	vector<double> V(3*N, 0);
	double Rk,R0,RR0,APP,AK,RRk,BK;
	int N1=3*N;  //�ڵ������ɶ�
    //	int Nup=240*N1; //����������
	int Nup=N1; //����������
	int K=0;
	//---------------------------------------------------------------------
	if(bnod_num!=0) displacement_nonzero_value(1,1,bnod_num,ip,vp,X,V);  //����λ�Ʒ���ֵԼ������
	//---------------------------------------------------------------------
	mabvm(N,N1,Iz,Ig,A,X,V);		//CALCULATE PRODUCT A*X0=> AP
	//---------------------------------------------------------------------
	for(int i=0; i<N1; i++)		//CALCULATE R0=P0=B-A*X0
	{
		R[i]=B[i]-V[i];
		P[i]=R[i];
	}
	//---------------------------------------------------------------------
	if(bnod_num!= 0) displacement_nonzero_value(0,1,bnod_num,ip,vp,P,R);		//����λ�Ʒ���ֵԼ������
	//---------------------------------------------------------------------
	RR0=0.0;
	for(int i=0; i<N1; i++)			//COMPUTE R(K)*R(K)
	{
		RR0=RR0+R[i]*R[i];
	}
	R0=sqrt(RR0);
	//---------------------------------------------------------------------
	//ENTER TO ITERATE
	//---------------------------------------------------------------------
    int test_sum = 0;
	while(K<Nup)
	{
		K=K+1;
		//---------------------------------------------------------------------
		mabvm(N,N1,Iz,Ig,A,P,S);		//CALCULATE AP=A*P=>S
   		//---------------------------------------------------------------------
		if(bnod_num!=0) displacement_nonzero_value(0,0,bnod_num,ip,vp,S,P);	//����λ�Ʒ���ֵԼ������
        
		APP=0.0;
		for(int i=0; i<N1; i++)	//CALCULATE INNER PRODUCT (AP,P)
		{
			APP=APP+S[i]*P[i];
		}
		AK=-RR0/APP;
		//---------------------------------------------------------------------
		//CALCULATE X(K)=X(K-1)-AK*P(K)
		//CALCULATE R(K)=R(K-1)+AK*AP(K)
		RRk=0.0;
		for(int i=0; i<N1; i++)
		{
			X[i]=X[i]-AK*P[i];
			R[i]=R[i]+AK*S[i];
			RRk=RRk+R[i]*R[i];
		}
		Rk=sqrt(RRk);
        
        if(K==test_sum+1)
        {
            hout << "test_sum=" << test_sum;
            if(fabs(R0)!=0.0) hout << "  error=" <<  fabs(Rk)/fabs(R0) << endl;
            test_sum +=100;
        }
        
        //		if(fabs(Rk)>=Zero*fabs(R0))
		if(fabs(Rk)>=1.0E-3*fabs(R0))
		{
			BK=RRk/RR0;
			RR0=RRk;
			for(int i=0; i<N1; i++)
			{
				P[i]=R[i]+BK*P[i];
			}
		}
		else	break;
	}
	if(K==Nup) hout << "ע�⣡�����ݶȷ������Է��������" <<K<< "����δ���������飡" << endl;
    
}
//---------------------------------------------------------------------------
//������ά�̶�λ��Լ����ķ���λ��ֵ
void SolveEqu::displacement_nonzero_value(const int &Kw, const int &Kg, const int &bnod_num, const vector<int> &ip, const vector<double> &vp, vector<double> &W, vector<double> &G)const
//Kw	��Kw!=0, ���̶�ֵ����W, ����Kg��ֵ, ������ֵ����W;
//Kg		��Kg!=0, ����ֵ����G.
//--------------------------------------------------------------------------------------------------------
//Type	u		v		w		�ձ�ʾ��һά�����ɱ߽�(�ڴ˴�ͳһ����ֵ�����ǲ��ᱻ����)
//   7		0		0		0		����0����ʾ����ֵ�����ֵ ��������ֵ
//	 6		��	0		0
//   5		0		��	0
//	 4		��	��	0
//   3		0		0		��
//   2		��	0		��
//   1		0		��	��
//--------------------------------------------------------------------------------------------------------
{
	for(int i=0; i<bnod_num; i++)
	{
		int key[3] = {1, 2, 4};
		int Ia = 2*i;
		int Ib = 3*i;
		int Ih = 3*ip[Ia];
		int Lkt = ip[Ia+1];
		for(int j=2; j>=0; j--)	//����һС��ѭ������Ҫ��key[0]~key[2]�зֱ��¼u,v,w������ֵ(key==1)������Ϊ��(key==0)
		{
			int item = Lkt%key[j];	//ȡ����
			key[j] = Lkt/key[j];		//ȡ��
			Lkt = item;
		}
		for(int j=0; j<=2; j++)
		{
			if(key[j]==1&&fabs(vp[Ib+j])>=Zero)
			{
				if(Kw!=0) W[Ih+j]=vp[Ib+j];
				else
				{
					if(Kg!=0) G[Ih+j]=0.0;
					W[Ih+j]=0.0;
				}
			}
		}
	}
}
//---------------------------------------------------------------------------
//����AK������U
void SolveEqu::mabvm(const int &N, const int &N1, const vector<long> &Iz, const vector<int> &Ig, const vector<double> &AK, const vector<double> &U, vector<double> &V)const
//AK	�նȾ���
//U		��������
//V		�������
{
	for(int i=0; i< N1; i++) V[i]=0.0;		//����ֵ
    
	for(int i=0; i<N; i++)
	{
		long Ik=6*(long)i+9*Iz[i];
		int Kh=3*i;
		for(int j=0; j<=2; j++)
		{
			for(int k=0; k<=2; k++)
			{
				if(j==2||k==2) V[Kh+j]=V[Kh+j]+AK[Ik+j+k+1]*U[Kh+k];
				else	V[Kh+j]=V[Kh+j]+AK[Ik+j+k]*U[Kh+k];
			}
		}
		if(i==0) continue;
		for(long m=Iz[i-1]; m<Iz[i]; m++)
		{
			Ik=6*(long)i+9*m;
			int Ia=3*Ig[m];
			for(int j=0; j<=2; j++)
			{
				for(int k=0; k<=2; k++)
				{
					V[Ia+j]=V[Ia+j]+AK[Ik+3*k+j]*U[Kh+k];
					V[Kh+j]=V[Kh+j]+AK[Ik+3*j+k]*U[Ia+k];
				}
			}
		}
	}
}
//---------------------------------------------------------------------------
//�������ڱ߽��������������ڵ�Ľ�
void  SolveEqu::Add_periodical_solution(const vector<int>* &peri_bnods, const vector<Node> nodes, vector<vector<double> > &solution)const
{
	//-----------------------------------------------------------
	//ѭ�����ڱ߽�ڵ�����
	int num[2] = { 0, peri_bnods[0][0] }; //num[0]��num[1]�ֱ��Ӧ������Ԫ��������Ӧ�á��е�xj��xi
	for(int i=0; i<(int)peri_bnods[0].size()-1; i++)
	{
		if(peri_bnods[1][i+1]==0) { num[1] = peri_bnods[0][i+1]; continue;} //һ��������������ʼ��
		num[0] = peri_bnods[0][i+1];
        
		//-----------------------------------------------------------
		//���������ڵ�ľ���
		double dist[3] = { nodes[num[1]].x-nodes[num[0]].x, nodes[num[1]].y-nodes[num[0]].y, nodes[num[1]].z-nodes[num[0]].z };		//����ڵ�xi��xj������
        
		for(int g=0; g<9; g++)
		{
			//-----------------------------------------------------------
			//���þ��Ȼ�Ӧ����Ϊ���ڱ߽�����
			double E[3][3] = {{0}, {0}, {0}}; //�������ڱ߽���������ʱ�ı��ξ�������
			switch(g)
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
                default: hout << "���� �������ڱ߽�����Լ��ʱ��ѭ������ֵ����" << g << "С��0���ߴ���8�����飡" << endl;
			}
			//-----------------------------------------------------------
			//����������������Q (xj=alpha*xi - Q; alpha =1������Ԫ��������Ӧ�á�), dist��xi��xj������, ����udiffi�൱��Q
			double udiff[3] = {0, 0, 0};
			for(int j=0; j<3; j++)
				for(int k=0; k<3; k++)
					udiff[j] += E[j][k]*dist[k];
            
			//-----------------------------------------------------------
			//�������ڱ߽�����������
			for(int j=0; j<3; j++)	solution[g][3*num[0]+j] = solution[g][3*num[1]+j] - udiff[j];
		}
	}
    
}
//---------------------------------------------------------------------------
//�����ȫ�ĸնȾ����Ҷ������ڼ��
void SolveEqu::Export_complete_matrix_equright(const vector<double> &A, const vector <double>* F,  const vector<Node> &nodes, const vector<long> &Iz, const vector<int> &Ig)const
{
	//---------------------------------------------------------------------------
	//������ȫ�ĸնȾ���
	vector<double> temp_vec(3*(int)nodes.size(),0.0);
	vector<vector<double> > full_matrix(3*(int)nodes.size(), temp_vec);
	for(int i=0; i<(int)Iz.size(); i++)
	{
		//�ԽǾ���
		long Mt =  6*(long)i + 9*Iz[i];
		for(int j=0; j<=2; j++)
			for(int k=0; k<=2; k++)
			{
				if(j==2||k==2) full_matrix[3*i+j][3*i+k] = A[Mt+j+k+1];
				else full_matrix[3*i+j][3*i+k] = A[Mt+j+k];
			}
		if(i!=0)
		{
			for(long j=Iz[i-1]; j<Iz[i]; j++)
			{
				Mt =  6*(long)i + 9*j;
				for(int k=0; k<=2; k++)
					for(int m=0; m<=2; m++)
					{
						full_matrix[3*i+k][3*Ig[j]+m] = A[Mt+3*k+m];
						full_matrix[3*Ig[j]+m][3*i+k] = A[Mt+3*k+m];
					}
			}
		}
	}
	//---------------------------------------------------------------------------
	//�����ȫ�ĸնȾ���鿴
	ofstream odata;
	odata.open("Complete_matrix_equright.dat");
    
	odata << "Complete_matrix:" << endl;
	odata << "//-------------------------------------------------------------------------------------------------" << endl;
	odata << setw(8) << setprecision(4)  << " ";
	for(int i=0; i<(int)nodes.size(); i++)
	{
		for(int j=0; j<3; j++)	odata << setw(8) << setprecision(4) << i;
		odata << " | ";
	}
	odata << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		for(int k=0; k<3; k++)
		{
			odata << setw(8) << setprecision(4) << i;
			for(int j=0; j<(int)nodes.size(); j++)
			{
				for(int m=0; m<3; m++)
				{
					if(fabs(full_matrix[3*i+k][3*j+m])<=Zero) odata << setw(8) << setprecision(4) << 0;
					else odata << setw(8) << setprecision(4) << full_matrix[3*i+k][3*j+m];
				}
				odata << " | ";
			}
			odata << endl;
		}
		odata << endl;
	}
	odata << endl << endl;
    
	//---------------------------------------------------------------------------
	//�����ȫ���Ҷ���鿴
	odata << "Complete_equright:" << endl;
	odata << "//-------------------------------------------------------------------------------------------------" << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		for(int k=0; k<3; k++)
		{
			odata << setw(10) << setprecision(4) << i;
			for(int j=0; j<9; j++)
			{
				if(fabs(F[j][3*i+k])<=Zero) odata << setw(15) << setprecision(4) << 0;
				else odata << setw(15) << setprecision(4) << F[j][3*i+k];
			}
			odata << endl;
		}
		odata << endl;
	}
	odata << endl << endl;
}
//---------------------------------------------------------------------------
//����ȫ�ĸնȾ����Ҷ��������ڱ߽������Ĵ���
void SolveEqu::Complete_matrix_equright_testing_periodical_constraints(const vector<int>* &peri_bnods, const vector<double> &A, const vector <double>* F, const vector<Node> &nodes, const vector<long> &Iz, const vector<int> &Ig)const
{
	//---------------------------------------------------------------------------
	//������ȫ�ĸնȾ�����Ҷ���
	vector<double> temp_vec(3*(int)nodes.size(),0.0);
	vector<vector<double> > full_matrix(3*(int)nodes.size(), temp_vec);
	for(int i=0; i<(int)Iz.size(); i++)
	{
		//�ԽǾ���
		long Mt =  6*(long)i + 9*Iz[i];
		for(int j=0; j<=2; j++)
			for(int k=0; k<=2; k++)
			{
				if(j==2||k==2) full_matrix[3*i+j][3*i+k] = A[Mt+j+k+1];
				else full_matrix[3*i+j][3*i+k] = A[Mt+j+k];
			}
		if(i!=0)
		{
			for(long j=Iz[i-1]; j<Iz[i]; j++)
			{
				Mt =  6*(long)i + 9*j;
				for(int k=0; k<=2; k++)
					for(int m=0; m<=2; m++)
					{
						full_matrix[3*i+k][3*Ig[j]+m] = A[Mt+3*k+m];
						full_matrix[3*Ig[j]+m][3*i+k] = A[Mt+3*k+m];
					}
			}
		}
	}
	vector<double> fullright[9];
	for(int i=0; i<9; i++) fullright[i] = F[i];
    
	//-----------------------------------------------------------
	//ѭ�����ڱ߽�ڵ�����
	int num[2] = { 0, peri_bnods[0][0] }; //num[0]��num[1]�ֱ��Ӧ������Ԫ��������Ӧ�á��е�xj��xi
	for(int i=0; i<(int)peri_bnods[0].size()-1; i++)
	{
		if(peri_bnods[1][i+1]==0) { num[1] = peri_bnods[0][i+1]; continue;} //һ��������������ʼ��
		num[0] = peri_bnods[0][i+1];
        
		//---------------------------------------------------------------------------
		//�����б任
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<9; k++) fullright[k][3*num[1]+j] += fullright[k][3*num[0]+j];
			for(int k=0; k<3*(int)nodes.size(); k++) full_matrix[3*num[1]+j][k] += full_matrix[3*num[0]+j][k];
		}
        
		//---------------------------------------------------------------------------
		//�����б任
		double dist[3] = { nodes[num[1]].x-nodes[num[0]].x, nodes[num[1]].y-nodes[num[0]].y, nodes[num[1]].z-nodes[num[0]].z };		//����ڵ�xi��xj������
		double udiff[9][3] = { {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0} };
        
		for(int g=0; g<9; g++)
		{
			//---------------------------------------------------------------------------
			//����Լ���Ĵ���
			//-----------------------------------------------------------
			//���þ��Ȼ�Ӧ����Ϊ���ڱ߽�����
			double E[3][3] = {{0}, {0}, {0}}; //�������ڱ߽���������ʱ�ı��ξ�������
			switch(g)
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
                default: hout << "���� �������ڱ߽�����Լ��ʱ��ѭ������ֵ����" << g << "С��0���ߴ���8�����飡" << endl;
			}
            
			for(int j=0; j<3; j++)
				for(int k=0; k<3; k++)
					udiff[g][j] += E[j][k]*dist[k];
		}
        
		for(int j=0; j<3*(int)nodes.size(); j++)
			for(int k=0; k<3; k++)
			{
				full_matrix[j][3*num[1]+k] += full_matrix[j][3*num[0]+k];
				for(int g=0; g<9; g++) fullright[g][j] += full_matrix[j][3*num[0]+k]*udiff[g][k];
			}
        
		//---------------------------------------------------------------------------
		//���㼰��λ��
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3*(int)nodes.size(); k++)
			{
				if(3*num[0]+j==k) full_matrix[3*num[0]+j][k] = 1.0;
				else
				{
					//������
					full_matrix[3*num[0]+j][k] = 0;
					//������
					full_matrix[k][3*num[0]+j] = 0;
				}
			}
			//�Ҷ�������
			for(int g=0; g<9; g++) fullright[g][3*num[0]+j] = 0;
		}
	}
    
	//---------------------------------------------------------------------------
	//�����ȫ�ĸնȾ���鿴
	ofstream odata;
	odata.open("Testing_complete_matrix_equright.dat");
    
	odata << "Complete_matrix:" << endl;
	odata << "//-------------------------------------------------------------------------------------------------" << endl;
	odata << setw(8) << setprecision(4)  << " ";
	for(int i=0; i<(int)nodes.size(); i++)
	{
		for(int j=0; j<3; j++)	odata << setw(8) << setprecision(4) << i;
		odata << " | ";
	}
	odata << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		for(int k=0; k<3; k++)
		{
			odata << setw(8) << setprecision(4) << i;
			for(int j=0; j<(int)nodes.size(); j++)
			{
				for(int m=0; m<3; m++)
				{
					if(fabs(full_matrix[3*i+k][3*j+m])<=Zero) odata << setw(8) << setprecision(4) << 0;
					else odata << setw(8) << setprecision(4) << full_matrix[3*i+k][3*j+m];
				}
				odata << " | ";
			}
			odata << endl;
		}
		odata << endl;
	}
	odata << endl << endl;
    
	//---------------------------------------------------------------------------
	//�����ȫ���Ҷ���鿴
	odata << "Complete_equright:" << endl;
	odata << "//-------------------------------------------------------------------------------------------------" << endl;
	for(int i=0; i<(int)nodes.size(); i++)
	{
		for(int k=0; k<3; k++)
		{
			odata << setw(10) << setprecision(4) << i;
			for(int j=0; j<9; j++)
			{
				if(fabs(fullright[j][3*i+k])<=Zero) odata << setw(15) << setprecision(4) << 0;
				else odata << setw(15) << setprecision(4) << fullright[j][3*i+k];
			}
			odata << endl;
		}
		odata << endl;
	}
	odata << endl << endl;
    
}
//---------------------------------------------------------------------------
//����������ڱ߽��Ͻڵ�Ľ�����ڱȶ�
void SolveEqu::Compared_periodical_bounday_solutions(const int &loop_num, const vector<int>* &peri_bnods, vector<double> &solution)const
{
	ofstream osol;
	if(loop_num==0) osol.open("Solution_in_periodical_boundary.dat");
	else osol.open("Solution_in_periodical_boundary.dat", ios::app);
	osol << endl << endl<< "//------------------------------------------------------------------------------------------------- " << endl;
	osol << "Solution "<< loop_num+1 << "  " << endl;
	for(int i=0; i<(int)peri_bnods[0].size(); i++)
	{
		if(i!=0&&peri_bnods[1][i]==0) osol << endl << endl;
		osol << setw(8) << peri_bnods[0][i]  << "  ";
		osol << setw(15) << setprecision(6)  << solution[3*peri_bnods[0][i]] << "  ";
		osol << setw(15) << setprecision(6) << solution[3*peri_bnods[0][i]+1] << "  ";
		osol << setw(15) << setprecision(6) << solution[3*peri_bnods[0][i]+2]<< endl;
	}
	osol.close();
    
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//����һ����Ϣ��������ע���У���"%"��ͷ����
string SolveEqu::Get_Line(ifstream &infile)const
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
