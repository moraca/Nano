//===========================================================================
// SolveEqu.h
// �����Է�����ͷ�ļ�
// A class to solving linear equations
//===========================================================================
#ifndef SOVLEEQU_H
#define SOVLEEQU_H

#include<assert.h>
#include<iostream>
#include<cmath>
#include<vector>
#include "Geometry_3D.h"
#include "Fem_3D.h"
#include "Hns.h"
using namespace hns;

//---------------------------------------------------------------------------
//����ⷽ����
class SolveEqu
{
public:
	//�����洢�նȾ���(�����洢�ڵ��Ӧ��ϵ)
	int izig(const vector<Node> &nodes, vector<long> &Iz, vector<int> &Ig)const;

	//�������ڱ߽�����Լ��
	void Periodical_boundary_constraints(const vector<Node> &nodes, const vector<int>* &peri_bnods, vector<long> &Iz, vector<int> &Ig, vector<double> &AK,  vector<double>* F)const;

	//���ù̶�λ��Լ��
	int Fixed_displacement_constraints(ifstream &infile, const vector<Node> &nodes, int &bnod_num, vector<int> &ip, vector<double> &vp)const;
	//������ά�̶�λ��Լ�������λ��ֵ
	void Deal_with_displacement_zero_value(const int &bnod_num, const int &nod_size, const vector<long> &Iz, const vector<int> &Ig, const vector<int> &Ip,
																			 const vector<double> &vp, vector <double>* F, vector <double> &AK)const;
	
	//������Է�����
	void Solve_linear_equations(const int &bnod_num, const int &N, const vector<long> &Iz, const vector<int> &Ig, const vector<int> &ip, const vector<double> &vp, 
													  const vector<double> &A, const vector<double> &B, vector<double> &X)const;

	//�������ڱ߽��������������ڵ�Ľ�
	void Add_periodical_solution(const vector<int>* &peri_bnods, const vector<Node> nodes, vector<vector<double> > &solution)const;

	//�����ȫ�ĸնȾ����Ҷ������ڼ��
	void Export_complete_matrix_equright(const vector<double> &A, const vector <double>* F,  const vector<Node> &nodes, const vector<long> &Iz, const vector<int> &Ig)const;

	//����ȫ�ĸնȾ����Ҷ��������ڱ߽������Ĵ���
	void Complete_matrix_equright_testing_periodical_constraints(const vector<int>* &peri_bnods, const vector<double> &A, const vector <double>* F, const vector<Node> &nodes, const vector<long> &Iz, const vector<int> &Ig)const;

	//����������ڱ߽��Ͻڵ�Ľ�����ڱȶ�
	void Compared_periodical_bounday_solutions(const int &loop_num, const vector<int>* &peri_bnods, vector<double> &solution)const;
private:
	//�����Ҷ�����
	void deal_with_F(const int num[], const double dist[3], const vector<vector<int> > &relative_small, const vector<vector<int> > &relative_large,
								   const vector<vector<double> > &temp_AKii, const vector<vector<vector<double> > > &temp_AK, vector<double>* F)const;
	
	//����AK[ii], AK[ij]��AK[jj]�Ĳ���
	void deal_with_AK_ii_ij_jj(const int num[], vector<vector<int> > &relative_small, vector<vector<int> > &relative_large, vector<vector<double> > &temp_AKii, vector<vector<vector<double> > > &temp_AK)const;

	//����AK[ki]��AK[kj]�Ĳ���
	void deal_with_AK_ki_kj(const int num[], vector<vector<int> > &relative_small, vector<vector<int> > &relative_large, vector<vector<vector<double> > > &temp_AK)const;	
	
	//������ά�̶�λ��Լ����ķ���λ��ֵ
	void displacement_nonzero_value(const int &Kw, const int &Kg, const int &bnod_num, const vector<int> &Ip, const vector<double> &vp, vector<double> &W, vector<double> &G)const;

	//����AK������U
	void mabvm(const int &N, const int &N1, const vector<long> &Iz, const vector<int> &Ig, const vector<double> &AK, const vector<double> &U, vector<double> &V)const;

	//������Ϣһ�У�����ע���У���%��ͷ��
	string Get_Line(ifstream &infile)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
