//===========================================================================
// SolveEqu.cpp
// 解线性方程组函数
// Member functions in a class to solving linear equations
//===========================================================================
#include "SolveEqu.h"

//---------------------------------------------------------------------------
//紧缩存储刚度矩阵
int SolveEqu::izig(const vector<Node> &nodes, vector<long> &Iz, vector<int> &Ig)const
{
	//iz[i-1], iz[i]中记录与第i个点相关的数据在ig中存储的开始位置和结束位置
	//第一个点的编号最小(没有比它编号还小的节点)，所以iz数据从第二个点(i==1)开始
	Iz.push_back(0);
    
	//从第二个节点开始循环所有的节点
	for(int i=1; i<(int)nodes.size(); i++)
	{
		for(int j=0; j<(int)nodes[i].relative_nods.size(); j++)	//注意这里节点的相关节点是从小到大排列的，并且不包括节点自身的编号
		{
			if(nodes[i].relative_nods[j]>=i) break;						//一旦一个编号大了（等于的情况其实不存在），说明后面的编号都大了
			else Ig.push_back(nodes[i].relative_nods[j]);
		}
		Iz.push_back((long)Ig.size());
	}
    
	//输出用于检查
	//hout << "Iz:" << endl;
	//for(int i=0; i<(int)Iz.size(); i++)
	//	hout << i << " " << Iz[i] << endl;
	//hout << endl << "Ig:" << endl;
	//for(int i=0; i<(long)Ig.size(); i++)
	//	hout << i << " " << Ig[i] << endl;
    
	return 1;
}
//---------------------------------------------------------------------------
//处理周期边界条件约束(右端向量部分)(算法参考黄艾香等人的论文，其中alph=1, Q=(q1,q2,q3);)
void SolveEqu::Periodical_boundary_constraints(const vector<Node> &nodes, const vector<int>* &peri_bnods, vector<long> &Iz, vector<int> &Ig,	vector<double> &AK,  vector<double>* F)const
{
	//-----------------------------------------------------------
	//初始化矩阵
	vector<vector<double> > temp_AKii;
	vector<vector<vector<double> > > temp_AK(nodes.size(), temp_AKii); //这里借用了temp_AKII初始化temp_AK
	vector<double> tempakii(6,0);
	vector<double> tempak(9,0);
    
	for(int i=0; i<6; i++) tempakii[i] = AK[i]; //把零号节点对角阵首先输入
	temp_AKii.push_back(tempakii);
	for(int i=1; i<(int)nodes.size(); i++) //从1号节点开始循环
	{
		long count=0;
		for(long j=Iz[i]-1; j>=Iz[i-1]; j--) //相关节点从大到小排列(便于以后插入和删除)
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
	//清空Iz, Ig和AK
	Iz.clear();
	Ig.clear();
	AK.clear();
    
	//-----------------------------------------------------------
	//初始化相关节点向量
	vector<vector<int> > relative_small(nodes.size());
	vector<vector<int> > relative_large(nodes.size());
	//循环所有的节点
	for(int i=0; i<(int)nodes.size(); i++)
	{
		for(int j=(int)nodes[i].relative_nods.size()-1; j>=0; j--)	//注意这里节点的相关节点是从小到大排列的，并且不包括节点自身的编号, 但是我们要取从大到小排列
		{
			if(nodes[i].relative_nods[j]>i) relative_large[i].push_back(nodes[i].relative_nods[j]);				//结果是从大到小排列(便于以后插入和删除)
			else if(nodes[i].relative_nods[j]<i) relative_small[i].push_back(nodes[i].relative_nods[j]);    //结果是从大到小排列(便于以后插入和删除)
		}
	}
    
	//-----------------------------------------------------------
	//循环周期边界节点向量
	int num[2] = { 0, peri_bnods[0][0] }; //num[0]和num[1]分别对应《有限元方法及其应用》中的xj和xi
    
	for(int i=0; i<(int)peri_bnods[0].size()-1; i++)
	{
		if(peri_bnods[1][i+1]==0) { num[1] = peri_bnods[0][i+1]; continue; } //一组周期条件的起始点
		num[0] = peri_bnods[0][i+1];
		//-----------------------------------------------------------
		//调整右端向量
		double dist[3] = { nodes[num[1]].x-nodes[num[0]].x, nodes[num[1]].y-nodes[num[0]].y, nodes[num[1]].z-nodes[num[0]].z };		//计算节点xi到xj的向量
		deal_with_F(num, dist, relative_small, relative_large, temp_AKii, temp_AK, F);
        
		//-----------------------------------------------------------
		//调整刚度矩阵
		//修改与i,j号节点所在行列相关部分
		//-----------------------------------------------------------
		//处理AK[ii], AK[ij]和AK[jj]的部分, 如果节点i与节点j原本相关, 则需要删除i与j的相关性, 并修改relative_small和relative_large
		deal_with_AK_ii_ij_jj(num, relative_small, relative_large, temp_AKii, temp_AK);
        
		//-----------------------------------------------------------
		//处理AK[ki]的部分
		//处理节点k (k!=i,j), 只要k与j相关, 则AK[ki]+=AK[kj](无论原本AK[ki]是否为零) //注意矩阵AK[kj]是否采用其转置
		//如果AK[ki]原本为零矩阵, 则还需要增加矩阵, 并修改relative_small和relative_large
		deal_with_AK_ki_kj(num, relative_small, relative_large, temp_AK);
	}
    
 	//-----------------------------------------------------------
	//重新生成Iz, Ig和AK
    
	//iz[i-1], iz[i]中记录与第i个点相关的数据在ig中存储的开始位置和结束位置
	//第一个点的编号最小(没有比它编号还小的节点)，所以iz数据从第二个点(i==1)开始
	Iz.push_back(0);
	AK = temp_AKii[0];			//先输入0号节点的下三角阵
    
	//从第二个节点开始循环所有的节点
	for(int i=1; i<(int)nodes.size(); i++)
	{
		for(int j=(int)relative_small[i].size()-1; j>=0; j--)		//注意这里相关节点是从大到小排列的， 所以要倒着输出
		{
			Ig.push_back(relative_small[i][j]);
		}
		Iz.push_back((long)Ig.size());
        
		for(int j=(int)temp_AK[i].size()-1; j>=0; j--)		//注意这里向量值是与relative_small的节点编号排序对应， 所以也要倒着输出
		{
			for(int k=0; k<9; k++) AK.push_back(temp_AK[i][j][k]);
		}
		for(int j=0; j<6; j++) AK.push_back(temp_AKii[i][j]);
	}
    
	//输出用于检查
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
//处理右端向量F
void SolveEqu::deal_with_F(const int num[], const double dist[3], const vector<vector<int> > &relative_small, const vector<vector<int> > &relative_large,
                           const vector<vector<double> > &temp_AKii, const vector<vector<vector<double> > > &temp_AK, vector<double>* F)const
{
	double udiff[9][3] = { {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0} };
	for(int i=0; i<9; i++)
	{
		//-----------------------------------------------------------
		//设置均匀化应变做为周期边界条件
		double E[3][3] = {{0}, {0}, {0}}; //用于周期边界条件处理时的变形矩阵条件
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
            default: hout << "错误！ 处理周期边界条件约束时，循环次数值等于" << i << "小于0或者大于8！请检查！" << endl;
		}
        
		//-----------------------------------------------------------
		//计算向量Q (xj=alpha*xi - Q; alpha =1《有限元方法及其应用》), dist是xi到xj的向量, 所以udiffi相当于Q
		for(int j=0; j<3; j++)
			for(int k=0; k<3; k++)
				udiff[i][j] += E[j][k]*dist[k];
	}
    
	//-----------------------------------------------------------
	//节点编号k<节点编号num[0]
	for(int j=0; j<(int)relative_small[num[0]].size(); j++)
	{
		for(int k=0; k<=2; k++)
		{
			const int rsn = 3*relative_small[num[0]][j]+k;
			for(int m=0; m<=2; m++)
			{
				const int km = k+3*m;
				for(int i=0; i<9; i++)
					F[i][rsn] += temp_AK[num[0]][j][km]*udiff[i][m];  //注意是乘矩阵的转置
			}
		}
	}
    
	//-----------------------------------------------------------
	//节点编号k>节点编号num[0]
	for(int j=0; j<(int)relative_large[num[0]].size(); j++)
	{
		//二分法查找j在k的relative_small中的位置
		int nk = relative_large[num[0]][j];
		if(nk==num[1]) continue; //k不能等于num[1], 这种情况单独在下面考虑
		bool mark = false;
		int left = 0;
		int middle = 0;
		int right = (int)relative_small[nk].size()-1;
		while(right>=left)
		{
			middle = (left + right)/2;
			if(relative_small[nk][middle] == num[0]) { mark = true; break; }				//找到相关位置
			else if(relative_small[nk][middle] < num[0]) right = middle - 1;
			else left = middle + 1;
		}
		if(!mark) hout << "错误！在节点" << nk << "的相关小编号节点中并没有找到节点" << num[0] << "，请检查！" << endl;
		for(int k=0; k<=2; k++)
		{
			const int tnk = 3*nk+k;
			for(int m=0; m<=2; m++)
			{
				const int km = 3*k+m;
				for(int i=0; i<9; i++)
					F[i][tnk] += temp_AK[nk][middle][km]*udiff[i][m];  //注意是乘矩阵本身
			}
		}
	}
    
	//-----------------------------------------------------------
	//节点编号num[1]
	
	//二分法在i的相关节点中查找j的位置，如果不相关就是零矩阵
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
			if(relative_small[num[1]][middle] == num[0]) { mark = true; nm = middle; break; }				//找到相关位置
			else if(relative_small[num[1]][middle] < num[0]) right = middle - 1;
			else left = middle + 1;
		}
	}
    
	const int tnum0 = 3*num[0];
	const int tnum1 = 3*num[1];
	for(int i=0; i<9; i++)
	{
		if(mark) //num[1]与num[0]相关
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
		else //num[1]与num[0]不相关, AK[ij]=0
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
		//节点编号num[0]
		for(int j=0; j<=2; j++) F[i][tnum0+j] = 0;
	}
}
//-----------------------------------------------------------
//处理AK[ii], AK[ij]和AK[jj]的部分
void SolveEqu::deal_with_AK_ii_ij_jj(const int num[], vector<vector<int> > &relative_small, vector<vector<int> > &relative_large, vector<vector<double> > &temp_AKii, vector<vector<vector<double> > > &temp_AK)const
{
	//如果两个点本身相关(即num[0]在num[1]的相关节点向量relative_small中, num[1]也在num[0]的相关节点向量relative_large中)，
	//那么将num[0]从num[1]的相关节点向量中删除(修改relative_small), 又将num[1]从num[0]的相关节点向量中删除(修改relative_large)
	
	bool mark[2] = { false, false };		//记录相关或者不相关
	int nodij = 0;										//记录j点在i的相关节点中的位置
    
	//-----------------------------------------------------------
	//在num[1]的relative_small中查找并删除num[0]（此处num[1]不可能等于0，因为至少有num[0]比num[1]小）
	if(num[0]<=relative_small[num[1]][0]&&num[0]>=relative_small[num[1]].back())
	{
		//二分法查找并删除
		int left = 0;
		int middle = 0;
		int right = (int)relative_small[num[1]].size()-1;
		while(right>=left)
		{
			middle = (left + right)/2;
			if(relative_small[num[1]][middle] == num[0]) { mark[0]=true; nodij =  middle; break; } //节点编号相同的情况
			else if(relative_small[num[1]][middle] < num[0]) right = middle - 1;
			else left = middle + 1;
		}
		if(mark[0])	//删除
		{
			for(int i=middle; i<(int)relative_small[num[1]].size()-1; i++) relative_small[num[1]][i] = relative_small[num[1]][i+1];
			relative_small[num[1]].pop_back();	//删除最后一个元素
		}
	}
    
	//-----------------------------------------------------------
	//在num[0]的relative_large中查找并删除num[1]（此处num[0]不可能等于最后一个节点的编号, 因为至少有num[1]比num[0]大）
	if(num[1]<=relative_large[num[0]][0]&&num[1]>=relative_large[num[0]].back())
	{
		//二分法查找并删除
		int left = 0;
		int middle = 0;
		int right = (int)relative_large[num[0]].size()-1;
		while(right>=left)
		{
			middle = (left + right)/2;
			if(relative_large[num[0]][middle] == num[1]) { mark[1]=true; break; } //节点编号相同的情况
			else if(relative_large[num[0]][middle] < num[1]) right = middle - 1;
			else left = middle + 1;
		}
		if(mark[0]!=mark[1]) hout << "错误！节点" << num[0] << "和节点" << num[1] << "之间的relative_small与relative_large对应关系有误，请检查！" << endl;
		if(mark[1])	//删除
		{
			for(int i=middle; i<(int)relative_large[num[0]].size()-1; i++) relative_large[num[0]][i] = relative_large[num[0]][i+1];
			relative_large[num[0]].pop_back();	//删除最后一个元素
		}
	}
    
	//-----------------------------------------------------------
	//修改AK[ii]和AK[ij]
	if(mark[0]) //节点i与节点j相关
	{
		//-----------------------------------------------------------
		//修改总刚阵AK[ii], 由于AK[ii]与AK[ij]的关联性, 所以要首先修改AK[ii]
		temp_AKii[num[1]][0] += temp_AKii[num[0]][0] + 2*temp_AK[num[1]][nodij][0];
		temp_AKii[num[1]][1] += temp_AKii[num[0]][1] + temp_AK[num[1]][nodij][1]+temp_AK[num[1]][nodij][3];
		temp_AKii[num[1]][2] += temp_AKii[num[0]][2] + 2*temp_AK[num[1]][nodij][4];
		temp_AKii[num[1]][3] += temp_AKii[num[0]][3] + temp_AK[num[1]][nodij][2]+temp_AK[num[1]][nodij][6];
		temp_AKii[num[1]][4] += temp_AKii[num[0]][4] + temp_AK[num[1]][nodij][5]+temp_AK[num[1]][nodij][7];
		temp_AKii[num[1]][5] += temp_AKii[num[0]][5] + 2*temp_AK[num[1]][nodij][8];
        
		//-----------------------------------------------------------
		//删除AK[ij]
		for(int i=nodij; i<(int)temp_AK[num[1]].size()-1; i++) temp_AK[num[1]][i] = temp_AK[num[1]][i+1];
		temp_AK[num[1]].pop_back();	//删除最后一个元素
	}
	else	//节点i与节点j不相关
	{
		//-----------------------------------------------------------
		//修改总刚阵AK[ii]
		temp_AKii[num[1]][0] += temp_AKii[num[0]][0];
		temp_AKii[num[1]][1] += temp_AKii[num[0]][1];
		temp_AKii[num[1]][2] += temp_AKii[num[0]][2];
		temp_AKii[num[1]][3] += temp_AKii[num[0]][3];
		temp_AKii[num[1]][4] += temp_AKii[num[0]][4];
		temp_AKii[num[1]][5] += temp_AKii[num[0]][5];
	}
    
	//-----------------------------------------------------------
	//修改总刚阵AK[jj]
	temp_AKii[num[0]][0] = temp_AKii[num[0]][2] = temp_AKii[num[0]][5] = 1;
	temp_AKii[num[0]][1] = temp_AKii[num[0]][3] = temp_AKii[num[0]][4] = 0;
    
}
//-----------------------------------------------------------
//处理AK[ki]和AK[kj]的部分
void SolveEqu::deal_with_AK_ki_kj(const int num[], vector<vector<int> > &relative_small, vector<vector<int> > &relative_large, vector<vector<vector<double> > > &temp_AK)const
{
	//-----------------------------------------------------------
	//循环relative_large[num[0]]
	for(int i=0; i<(int)relative_large[num[0]].size(); i++)
	{
		int nodkj = 0; //记录j在k的相关节点中的位置编号
		//--------------------------------------------------------------------------------------------------------------------------------------------
		int nk = relative_large[num[0]][i];
		//在relative_small[nk]中寻找num[0]的位置
		if(num[0]<=relative_small[nk][0]&&num[0]>=relative_small[nk].back())
		{
			//二分法查找并删除
			bool mark = false;
			int left = 0;
			int middle = 0;
			int right = (int)relative_small[nk].size()-1;
			while(right>=left)
			{
				middle = (left + right)/2;
				if(relative_small[nk][middle] == num[0]) { mark = true; nodkj =  middle; break; } //节点编号相同的情况
				else if(relative_small[nk][middle] < num[0]) right = middle - 1;
				else left = middle + 1;
			}
			if(!mark) hout << "错误！节点" << nk << "的relative_small相关节点中没找到" << num[0] << "，请检查！" << endl;
		}
		else hout << "错误！节点" << nk << "的relative_small相关节点中不包含" << num[0] << "，请检查！" << endl;
        
		//------------------------------------------------------------------------------------------------------------------------------------------------
		//找寻nk和num[1]的关系并修改
		if(nk == num[1]) hout << "错误！节点" << num[0] << "的相关节点中仍有节点" << num[1] << "存在，请检查！" << endl;
		else if(nk > num[1])
		{
			bool mark[2] = {false, false};
			//在nk的relative_small中查找num[1], 并修改temp_AK, 在插入num[1]相关的同时删除num[0]
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
			else //num[1]在relative_small[nk]之间
			{
				//二分法查找
				int left = 0;
				int middle = 0;
				int right = (int)relative_small[nk].size()-1;
				while(right>=left)
				{
					middle = (left + right)/2;
					if(relative_small[nk][middle] == num[1]) { mark[0] = true; break; } //节点编号相同的情况
					else if(relative_small[nk][middle] < num[1]) right = middle - 1;
					else left = middle + 1;
				}
				if(mark[0]) //num[1]已经是nk的相关节点
				{
					for(int i=nodkj; i<(int)relative_small[nk].size()-1; i++)	relative_small[nk][i] = relative_small[nk][i+1];
					relative_small[nk].pop_back();	//删除最后一个元素
                    
					for(int j=0; j<9; j++)	temp_AK[nk][middle][j] += temp_AK[nk][nodkj][j];  //叠加
                    
					for(int i=nodkj; i<(int)temp_AK[nk].size()-1; i++) temp_AK[nk][i] = temp_AK[nk][i+1];
					temp_AK[nk].pop_back();	//删除最后一个元素
				}
				else
				{
					if(nodkj<left) hout << "错误！节点" << num[0] << "和" << num[1] << "都在节点" << nk << "的relative_small内，但排列位置出错，请检查！" << endl;
					for(int j=nodkj; j>left; j--) relative_small[nk][j] = relative_small[nk][j-1];
					relative_small[nk][left] = num[1];
                    
					vector<double> tempak(temp_AK[nk][nodkj]);
					for(int j=nodkj; j>left; j--) temp_AK[nk][j] = temp_AK[nk][j-1];
					temp_AK[nk][left] = tempak;
				}
			}
            
			//在num[1]的relative_large中查找nk并修改
			if(nk>relative_large[num[1]][0])
			{
				relative_large[num[1]].push_back(0);
				for(int j=(int)relative_large[num[1]].size()-1; j>0; j--) relative_large[num[1]][j] = relative_large[num[1]][j-1];
				relative_large[num[1]][0] = nk;
			}
			else if(nk<relative_large[num[1]].back())	relative_large[num[1]].push_back(nk);
			else //nk在relative_large[num[1]]之间
			{
				//二分法查找
				int left = 0;
				int middle = 0;
				int right = (int)relative_large[num[1]].size()-1;
				while(right>=left)
				{
					middle = (left + right)/2;
					if(relative_large[num[1]][middle] == nk) { mark[1] = true; break; } //节点编号相同的情况
					else if(relative_large[num[1]][middle] < nk) right = middle - 1;
					else left = middle + 1;
				}
				if(mark[0]!=mark[1]) hout << "错误！节点" << num[0] << "和节点" << nk << "的相关性不一致，请检查！" << endl;
				if(!mark[1]) //nk不是num[1]的相关节点
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
			//在num[1]的relative_small中查找nk, 并修改temp_AK
			if(nk > relative_small[num[1]][0])
			{
				relative_small[num[1]].push_back(0);
				for(int j=(int)relative_small[num[1]].size()-1; j>0; j--) relative_small[num[1]][j] = relative_small[num[1]][j-1];
				relative_small[num[1]][0] = nk;
				
				vector<double> tempak(9,0);
				temp_AK[num[1]].push_back(tempak);
				for(int j=(int)temp_AK[num[1]].size()-1; j>0; j--) temp_AK[num[1]][j] = temp_AK[num[1]][j-1];
				for(int j=0; j<=2; j++)											//当j<k<i时要加上AK[kj]的转置矩阵
					for(int m=0; m<=2; m++)
						temp_AK[num[1]][0][3*j+m] = temp_AK[nk][nodkj][j+3*m];
			}
			else if(nk < relative_small[num[1]].back())
			{
				relative_small[num[1]].push_back(nk);
                
				vector<double> tempak(9,0);
				temp_AK[num[1]].push_back(tempak);
				for(int j=0; j<=2; j++)											//当j<k<i时要加上AK[kj]的转置矩阵
					for(int m=0; m<=2; m++)
						temp_AK[num[1]].back()[3*j+m] = temp_AK[nk][nodkj][j+3*m];
			}
			else //nk在relative_small[num[1]]之间
			{
				//二分法查找
				int left = 0;
				int middle = 0;
				int right = (int)relative_small[num[1]].size()-1;
				while(right>=left)
				{
					middle = (left + right)/2;
					if(relative_small[num[1]][middle] == nk) { mark[0] = true; break; } //节点编号相同的情况
					else if(relative_small[num[1]][middle] < nk) right = middle - 1;
					else left = middle + 1;
				}
				if(mark[0]) //nk已经是num[1]的相关节点
				{
					for(int j=0; j<=2; j++)											//当j<k<i时要加上AK[kj]的转置矩阵
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
					for(int j=0; j<=2; j++)											//当j<k<i时要加上AK[kj]的转置矩阵
						for(int m=0; m<=2; m++)
							temp_AK[num[1]][left][3*j+m] = temp_AK[nk][nodkj][j+3*m];
				}
			}
            
			//在nk的relative_large中查找num[1]并修改
			if(num[1] > relative_large[nk][0])
			{
				relative_large[nk].push_back(0);
				for(int j=(int)relative_large[nk].size()-1; j>0; j--) relative_large[nk][j] = relative_large[nk][j-1];
				relative_large[nk][0] = num[1];
			}
			else if(num[1] < relative_large[nk].back())	relative_large[nk].push_back(num[1]);
			else //num[1]在relative_large[nk]之间
			{
				//二分法查找
				int left = 0;
				int middle = 0;
				int right = (int)relative_large[nk].size()-1;
				while(right>=left)
				{
					middle = (left + right)/2;
					if(relative_large[nk][middle] == num[1]) { mark[1] = true; break; } //节点编号相同的情况
					else if(relative_large[nk][middle] < num[1]) right = middle - 1;
					else left = middle + 1;
				}
				if(mark[0]!=mark[1]) hout << "错误！节点" << num[0] << "和节点" << nk << "的相关性不一致，请检查！" << endl;
				if(!mark[1]) //num[1]不是nk的相关节点,
				{
					relative_large[nk].push_back(0);
					for(int j=(int)relative_large[nk].size()-1; j>left; j--) relative_large[nk][j] = relative_large[nk][j-1];
					relative_large[nk][left] = num[1];
				}
			}
            
			//删除relative_small[nk]中的num[0]
			for(int j=nodkj; j<(int)relative_small[nk].size()-1; j++) relative_small[nk][j] = relative_small[nk][j+1];
			relative_small[nk].pop_back();
            
			//删除temp_AK[nk]中的num[0]的相关矩阵
			for(int j=nodkj; j<(int)temp_AK[nk].size()-1; j++) temp_AK[nk][j] = temp_AK[nk][j+1];
			temp_AK[nk].pop_back();
		}
	}
	//清空relative_large[num[0]]
	relative_large[num[0]].clear();
    
	//-----------------------------------------------------------
	//循环relative_small[num[0]], 这个量一定也比num[1]小
	for(int i=0; i<(int)relative_small[num[0]].size(); i++)
	{
		int nk = relative_small[num[0]][i];
        
		bool mark[2] = {false, false};
		//在num[1]的relative_small中查找nk, 并修改temp_AK
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
		else //nk在relative_small[num[1]]之间
		{
			//二分法查找
			int left = 0;
			int middle = 0;
			int right = (int)relative_small[num[1]].size()-1;
			while(right>=left)
			{
				middle = (left + right)/2;
				if(relative_small[num[1]][middle] == nk) { mark[0] = true; break; } //节点编号相同的情况
				else if(relative_small[num[1]][middle] < nk) right = middle - 1;
				else left = middle + 1;
			}
			if(mark[0]) //nk已经是num[1]的相关节点
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
        
		int nodkj = 0; //记录j在k的相关节点中的位置编号
		//在relative_large[nk]中寻找num[0]的位置
		if(num[0]<=relative_large[nk][0]&&num[0]>=relative_large[nk].back())
		{
			//二分法查找并删除
			bool mark = false;
			int left = 0;
			int middle = 0;
			int right = (int)relative_large[nk].size()-1;
			while(right>=left)
			{
				middle = (left + right)/2;
				if(relative_large[nk][middle] == num[0]) { mark = true; nodkj =  middle; break; } //节点编号相同的情况
				else if(relative_large[nk][middle] < num[0]) right = middle - 1;
				else left = middle + 1;
			}
			if(!mark) hout << "错误！节点" << nk << "的relative_large相关节点中没找到" << num[0] << "，请检查！" << endl;
		}
		else hout << "错误！节点" << nk << "的relative_large相关节点中不包含" << num[0] << "，请检查！" << endl;
        
        
		//在nk的relative_large中寻找num[1]并插入num[1], 删除num[0]
		if(num[1]>relative_large[nk][0])
		{
			for(int j=nodkj; j>0; j--) relative_large[nk][j] = relative_large[nk][j-1];
			relative_large[nk][0] = num[1];
		}
		else //num[1]在relative_large[nk]的0和nodkj之间找，num[1]不可能小于relative_large[nk].back(), 因为总有num[0]<num[1]
		{
			//二分法查找
			int left = 0;
			int middle = 0;
			int right = nodkj-1;
			while(right>=left)
			{
				middle = (left + right)/2;
				if(relative_large[nk][middle] == num[1]) { mark[1] = true; break; } //节点编号相同的情况
				else if(relative_large[nk][middle] < num[1]) right = middle - 1;
				else left = middle + 1;
			}
			if(mark[0]!=mark[1]) hout << "错误！节点" << num[1] << "和节点" << nk << "的相关性不一致，请检查！" << endl;
			if(!mark[1])  //num[1]不是nk的相关节点
			{
				for(int j=nodkj; j>left; j--) relative_large[nk][j] = relative_large[nk][j-1];
				relative_large[nk][left] = num[1];
			}
			else //num[1]是nk的相关节点，删除num[0]
			{
				for(int j=nodkj; j<(int)relative_large[nk].size()-1; j++) relative_large[nk][j] = relative_large[nk][j+1];
				relative_large[nk].pop_back();
			}
		}
	}
	//清空relative_small[num[0]]
	relative_small[num[0]].clear();
	//清空temp_AK[num[0]]
	temp_AK[num[0]].clear();
}
//---------------------------------------------------------------------------
//设置固定位移约束
int SolveEqu::Fixed_displacement_constraints(ifstream &infile, const vector<Node> &nodes, int &bnod_num, vector<int> &ip, vector<double> &vp)const
{
	//-----------------------------------------------------------
	//读入位移数据
	int fixed_displacement_num;
	istringstream istn(Get_Line(infile));
	istn >> fixed_displacement_num;		//读入位移数据
	while(fixed_displacement_num)
	{
		//约束条件数减一
		fixed_displacement_num--;
		//读入约束条件
		istringstream istr(Get_Line(infile));
		//约束条件类型
		Point_3D fpoi;
		Plane_3D fpla;
		string Boundary_case;
		istr >> Boundary_case;
		if(Boundary_case=="Point") istr >> fpoi.x >> fpoi.y >> fpoi.z;
		else if(Boundary_case=="Surface") istr >> fpla.coef[0] >> fpla.coef[1] >> fpla.coef[2] >> fpla.coef[3];
		else { hout << "注意！约束条件不是点或者面，请重新输入位移约束数据！" << endl; return 0; }
		//约束值
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
				if(sign>=1) { hout << "注意！位移约束条件中标识重复输入或者输入顺序颠倒，请重新输入约束数据！" << endl; return 0; }
				sign += 1;
				istr >> val_temp;		//读入数据
				value.push_back(val_temp);
			}
			else if(str_temp=="Fixed_displacement_y")
			{
				if(sign>=2) { hout << "注意！位移约束条件中标识重复输入或者输入顺序颠倒，请重新输入约束数据！" << endl; return 0; }
				sign += 2;
				istr >> val_temp;		//读入数据
				value.push_back(val_temp);
			}
			else if(str_temp=="Fixed_displacement_z")
			{
				if(sign>=4) { hout << "注意！位移约束条件中标识重复输入或者输入顺序颠倒，请重新输入约束数据！" << endl; return 0; }
				sign += 4;
				istr >> val_temp;		//读入数据
				value.push_back(val_temp);
			}
			else { hout << "注意！约束条件中固定位移条件标识有错，或者是此行末尾有空格，请重新输入位移约束数据！" << endl; return 0; }
		}
		if((int)value.size()==0||(int)value.size()>3)
		{
			hout << "注意！约束条件中固定位移条件没有输入或者输入过多，请重新输入约束数据！" << endl;
			return 0;
		}
        
		//建立位移约束条件
		int count = 0;
		for(int i=0; i<(int)nodes.size(); i++)
		{
			if((Boundary_case=="Point"&&fpoi.distance_to(nodes[i].x, nodes[i].y, nodes[i].z)==0.0)||			//到固定点的距离等于0
               (Boundary_case=="Surface"&&fpla.contain(nodes[i].x, nodes[i].y, nodes[i].z)==1))					//点在一个空间平面上
			{
				count ++;	//记录找到一个点
				ip.push_back(i);
				ip.push_back(sign);
				//初始化边界位移约束点信息
				//--------------------------------------------------------------------------------------------------------
				//Type	u		v		w		空表示这一维是自由边界(在此处统一赋零值，但是不会被调用)
				//   7		0		0		0		其它0处表示是零值或非零值 但必须有值
				//	 6		空	0		0
				//   5		0		空	0
				//	 4		空	空	0
				//   3		0		0		空
				//   2		空	0		空
				//   1		0		空	空
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
					default:	hout << "注意！位移约束点类型标记有误，请重新输入位移约束条件！" << endl;	return 0;
				}
				bnod_num++;  //边界节点个数加一
			}
		}
		if(count==0) { hout << "注意！位移约束点位置输入有误（可能不经过任何网格节点），请重新输入位移约束条件！" << endl; return 0; }
	}
    
	return 1;
}
//---------------------------------------------------------------------------
//处理三维固定位移约束点的零位移值
void SolveEqu::Deal_with_displacement_zero_value(const int &bnod_num, const int &nod_size, const vector<long> &Iz, const vector<int> &Ig, const vector<int> &ip,
                                                 const vector<double> &vp, vector <double>* F, vector <double> &AK)const
//--------------------------------------------------------------------------------------------------------
//Type	u		v		w		空表示这一维是自由边界(在此处统一赋零值，但是不会被调用)
//   7		0		0		0		其它0处表示是零值或非零值 但必须有值
//	 6		空	0		0
//   5		0		空	0
//	 4		空	空	0
//   3		0		0		空
//   2		空	0		空
//   1		0		空	空
//--------------------------------------------------------------------------------------------------------
{
	for(int i=0; i<bnod_num; i++)
	{
		int key[3] = {1, 2, 4};
		int Ia = 2*i;
		int Ih = 3*i;
		int Mt = ip[Ia];
		int Lkt = ip[Ia+1];
		for(int j=2; j>=0; j--)	 //以下一小段循环体是要在key[0]~key[2]中分别记录u,v,w哪项有值(key==1)而哪项为空(key==0)
		{
			int item = Lkt%key[j];	//取余数
			key[j] = Lkt/key[j];		//取商
			Lkt = item;
		}
		for(int j=0; j<=2; j++)
		{
			if(key[j]==1&&fabs(vp[Ih+j])<Zero)
			{
				long Kmt = 0;
				//---------------------------------------------------------------------
				for(int k=0; k<9; k++)	F[k][3*Mt+j]=0.0;		//右端项
				//---------------------------------------------------------------------
				if(Mt!=0)				//Mt的相关节点中比Mt编号小的
				{
					for(long k=Iz[Mt-1]; k<Iz[Mt]; k++)
						for(int m=0; m<=2; m++)
						{
							Kmt = 6*(long)Mt+9*k+3*j+m;
							AK[Kmt]=0.0;
						}
				}
				//---------------------------------------------------------------------
				for(int k=Mt+1; k<nod_size; k++)		//Mt的相关节点中比Mt编号大的
				{
					if(Mt>=Ig[Iz[k-1]]&&Mt<=Ig[Iz[k]-1])
					{
						//二分法查找
						bool mark = false;
						long left = Iz[k-1];
						long middle = 0;
						long right = Iz[k]-1;
						while(right>=left)
						{
							middle = (left + right)/2;
							if(Ig[middle] == Mt) { mark = true; break; } //节点编号相同的情况
							else if(Ig[middle] > Mt) right = middle - 1;
							else left = middle + 1;
						}
						if(mark)  //找到对应点
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
				for(int k=0; k<=2; k++)  //计算对角矩阵
				{
					if(j==2||k==2) Kmt=6*(long)Mt+9*Iz[Mt]+j+k+1;				//相当于	[0,1,	3]		当j==0时, 消 0,	1,	3;
					else Kmt=6*(long)Mt+9*Iz[Mt]+j+k;									//				[1,2,	4]		当j==1时, 消	1,	2,	4;
 					if(k==j) AK[Kmt]=1.0;																//				[3,4,	5]		当 j==2时 消 3,	4,	5;
					else AK[Kmt]=0.0;
				}
			}
		}
	}
}
//---------------------------------------------------------------------------
//求解线性方程组函数(共轭梯度CONJUGATE GRADIENT METHOD)
void SolveEqu::Solve_linear_equations(const int &bnod_num, const int &N, const vector<long> &Iz, const vector<int> &Ig, const vector<int> &ip,
                                      const vector<double> &vp, const vector<double> &A, const vector<double> &B, vector<double> &X)const
{
	//---------------------------------------------------------------------
	vector<double> P(3*N, 0);
	vector<double> R(3*N, 0);
	vector<double> S(3*N, 0);
	vector<double> V(3*N, 0);
	double Rk,R0,RR0,APP,AK,RRk,BK;
	int N1=3*N;  //节点总自由度
    //	int Nup=240*N1; //最大迭代步数
	int Nup=N1; //最大迭代步数
	int K=0;
	//---------------------------------------------------------------------
	if(bnod_num!=0) displacement_nonzero_value(1,1,bnod_num,ip,vp,X,V);  //处理位移非零值约束条件
	//---------------------------------------------------------------------
	mabvm(N,N1,Iz,Ig,A,X,V);		//CALCULATE PRODUCT A*X0=> AP
	//---------------------------------------------------------------------
	for(int i=0; i<N1; i++)		//CALCULATE R0=P0=B-A*X0
	{
		R[i]=B[i]-V[i];
		P[i]=R[i];
	}
	//---------------------------------------------------------------------
	if(bnod_num!= 0) displacement_nonzero_value(0,1,bnod_num,ip,vp,P,R);		//处理位移非零值约束条件
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
		if(bnod_num!=0) displacement_nonzero_value(0,0,bnod_num,ip,vp,S,P);	//处理位移非零值约束条件
        
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
	if(K==Nup) hout << "注意！共轭梯度法解线性方程组迭代" <<K<< "次仍未收敛，请检查！" << endl;
    
}
//---------------------------------------------------------------------------
//处理三维固定位移约束点的非零位移值
void SolveEqu::displacement_nonzero_value(const int &Kw, const int &Kg, const int &bnod_num, const vector<int> &ip, const vector<double> &vp, vector<double> &W, vector<double> &G)const
//Kw	当Kw!=0, 将固定值赋予W, 无论Kg的值, 否则将零值赋予W;
//Kg		当Kg!=0, 将零值赋予G.
//--------------------------------------------------------------------------------------------------------
//Type	u		v		w		空表示这一维是自由边界(在此处统一赋零值，但是不会被调用)
//   7		0		0		0		其它0处表示是零值或非零值 但必须有值
//	 6		空	0		0
//   5		0		空	0
//	 4		空	空	0
//   3		0		0		空
//   2		空	0		空
//   1		0		空	空
//--------------------------------------------------------------------------------------------------------
{
	for(int i=0; i<bnod_num; i++)
	{
		int key[3] = {1, 2, 4};
		int Ia = 2*i;
		int Ib = 3*i;
		int Ih = 3*ip[Ia];
		int Lkt = ip[Ia+1];
		for(int j=2; j>=0; j--)	//以下一小段循环体是要在key[0]~key[2]中分别记录u,v,w哪项有值(key==1)而哪项为空(key==0)
		{
			int item = Lkt%key[j];	//取余数
			key[j] = Lkt/key[j];		//取商
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
//矩阵AK乘向量U
void SolveEqu::mabvm(const int &N, const int &N1, const vector<long> &Iz, const vector<int> &Ig, const vector<double> &AK, const vector<double> &U, vector<double> &V)const
//AK	刚度矩阵
//U		输入向量
//V		结果向量
{
	for(int i=0; i< N1; i++) V[i]=0.0;		//赋初值
    
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
//根据周期边界条件附加其他节点的解
void  SolveEqu::Add_periodical_solution(const vector<int>* &peri_bnods, const vector<Node> nodes, vector<vector<double> > &solution)const
{
	//-----------------------------------------------------------
	//循环周期边界节点向量
	int num[2] = { 0, peri_bnods[0][0] }; //num[0]和num[1]分别对应《有限元方法及其应用》中的xj和xi
	for(int i=0; i<(int)peri_bnods[0].size()-1; i++)
	{
		if(peri_bnods[1][i+1]==0) { num[1] = peri_bnods[0][i+1]; continue;} //一组周期条件的起始点
		num[0] = peri_bnods[0][i+1];
        
		//-----------------------------------------------------------
		//计算两个节点的距离
		double dist[3] = { nodes[num[1]].x-nodes[num[0]].x, nodes[num[1]].y-nodes[num[0]].y, nodes[num[1]].z-nodes[num[0]].z };		//计算节点xi到xj的向量
        
		for(int g=0; g<9; g++)
		{
			//-----------------------------------------------------------
			//设置均匀化应变做为周期边界条件
			double E[3][3] = {{0}, {0}, {0}}; //用于周期边界条件处理时的变形矩阵条件
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
                default: hout << "错误！ 处理周期边界条件约束时，循环次数值等于" << g << "小于0或者大于8！请检查！" << endl;
			}
			//-----------------------------------------------------------
			//计算周期条件向量Q (xj=alpha*xi - Q; alpha =1《有限元方法及其应用》), dist是xi到xj的向量, 所以udiffi相当于Q
			double udiff[3] = {0, 0, 0};
			for(int j=0; j<3; j++)
				for(int k=0; k<3; k++)
					udiff[j] += E[j][k]*dist[k];
            
			//-----------------------------------------------------------
			//根据周期边界条件修正解
			for(int j=0; j<3; j++)	solution[g][3*num[0]+j] = solution[g][3*num[1]+j] - udiff[j];
		}
	}
    
}
//---------------------------------------------------------------------------
//输出完全的刚度矩阵及右端项用于检测
void SolveEqu::Export_complete_matrix_equright(const vector<double> &A, const vector <double>* F,  const vector<Node> &nodes, const vector<long> &Iz, const vector<int> &Ig)const
{
	//---------------------------------------------------------------------------
	//计算完全的刚度矩阵
	vector<double> temp_vec(3*(int)nodes.size(),0.0);
	vector<vector<double> > full_matrix(3*(int)nodes.size(), temp_vec);
	for(int i=0; i<(int)Iz.size(); i++)
	{
		//对角矩阵
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
	//输出完全的刚度矩阵查看
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
	//输出完全的右端项查看
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
//用完全的刚度矩阵及右端项检测周期边界条件的处理
void SolveEqu::Complete_matrix_equright_testing_periodical_constraints(const vector<int>* &peri_bnods, const vector<double> &A, const vector <double>* F, const vector<Node> &nodes, const vector<long> &Iz, const vector<int> &Ig)const
{
	//---------------------------------------------------------------------------
	//计算完全的刚度矩阵和右端项
	vector<double> temp_vec(3*(int)nodes.size(),0.0);
	vector<vector<double> > full_matrix(3*(int)nodes.size(), temp_vec);
	for(int i=0; i<(int)Iz.size(); i++)
	{
		//对角矩阵
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
	//循环周期边界节点向量
	int num[2] = { 0, peri_bnods[0][0] }; //num[0]和num[1]分别对应《有限元方法及其应用》中的xj和xi
	for(int i=0; i<(int)peri_bnods[0].size()-1; i++)
	{
		if(peri_bnods[1][i+1]==0) { num[1] = peri_bnods[0][i+1]; continue;} //一组周期条件的起始点
		num[0] = peri_bnods[0][i+1];
        
		//---------------------------------------------------------------------------
		//初等行变换
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<9; k++) fullright[k][3*num[1]+j] += fullright[k][3*num[0]+j];
			for(int k=0; k<3*(int)nodes.size(); k++) full_matrix[3*num[1]+j][k] += full_matrix[3*num[0]+j][k];
		}
        
		//---------------------------------------------------------------------------
		//初等列变换
		double dist[3] = { nodes[num[1]].x-nodes[num[0]].x, nodes[num[1]].y-nodes[num[0]].y, nodes[num[1]].z-nodes[num[0]].z };		//计算节点xi到xj的向量
		double udiff[9][3] = { {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0} };
        
		for(int g=0; g<9; g++)
		{
			//---------------------------------------------------------------------------
			//周期约束的处理
			//-----------------------------------------------------------
			//设置均匀化应变做为周期边界条件
			double E[3][3] = {{0}, {0}, {0}}; //用于周期边界条件处理时的变形矩阵条件
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
                default: hout << "错误！ 处理周期边界条件约束时，循环次数值等于" << g << "小于0或者大于8！请检查！" << endl;
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
		//清零及单位化
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3*(int)nodes.size(); k++)
			{
				if(3*num[0]+j==k) full_matrix[3*num[0]+j][k] = 1.0;
				else
				{
					//行清零
					full_matrix[3*num[0]+j][k] = 0;
					//列清零
					full_matrix[k][3*num[0]+j] = 0;
				}
			}
			//右端项清零
			for(int g=0; g<9; g++) fullright[g][3*num[0]+j] = 0;
		}
	}
    
	//---------------------------------------------------------------------------
	//输出完全的刚度矩阵查看
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
	//输出完全的右端项查看
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
//用于输出周期边界上节点的结果用于比对
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
//读入一行信息，并跳过注释行（以"%"开头）；
string SolveEqu::Get_Line(ifstream &infile)const
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
