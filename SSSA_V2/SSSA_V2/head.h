#pragma once
#ifndef HEAD_H
#define HEAD_H

#include<iostream>
#include<fstream>
#include<sstream>		
#include<iomanip>
#include<cmath>
#include<cstring>
#include<complex>
#include<ctime>
#include<vector>

using namespace std;

#define PI 3.14159265358979
#define SITA 57.295779513082323		//180 / pi
#define base 100					//����ֵ��׼

struct NODE				//�ڵ����ݣ��������нڵ㣩
{
	int no;				//�ڵ��
	string name;		//�ڵ�����
	double vol;			//�ڵ��׼��ѹ��kV��
};

struct NODE_DATA		//�ڵ㳱�����ݣ��������нڵ㣩
{
	int no;				//�����no�ǳ�������� �Ż���ź��
	int type;			//�ڵ����� -1=������ڵ� 0=����ڵ� 1=���ɽڵ�
	string name;
	double vol_base, vol_pu, vol, vol_ang;
	double P, Q;		//����Ϊ��
};


struct GENE_INIT		//�������ֵ
{
	double Delta0, Vd0, Vq0, Id0, Iq0;
	double Efq0, Eq10, Eq20, Ed10, Ed20;
};


struct GENE_PARA		//���������
{
	//ͬ��������
	string card_type;		//�����������
	string gene_model;		//�����ģ�� ��four_winding, three_winding, two_axle, Eq_change, Eq_constant, classic
	string name;			//���������
	double vol;				//ĸ�߻�׼��ѹ(kV)
	double capacity;		//��������(MVA)
	double power_factor;	//�����������
	string type;			//�������
	double Ek;				//���������(MW s��
	double Sh;				//��׼������MVA��

	double Xd, Xd1, Xd2;
	double Xq, Xq1, Xq2;
	double Ra;
	double Td01, Td02;
	double Tq01, Tq02;
	double Tj;
	double D;				//����ϵ��
	double XL;				//���������©��
	double A, B, N;			//���������ϵ��
};//XL, A, B, N��������δ�õ�


struct LOAD_PARA		//���ɲ���
{
	string card_type;		//���ɿ�����
	string load_model;		//���ɿ�ģ��
	string name;
	double vol;
	double ap, bp, cp;
	double aq, bq, cq;
};




//////////////////////////����
void nrerror(char error_text[]);//��ֵ�������ʱ���������
int **imatrix(long nrow, long ncol);//����һ��int�͵ľ���nrow��ncol�ף��õ��˶�ά���齵ά����
double **dmatrix(long nrow, long ncol);//����һ��double�͵ľ���nrow��ncol�ף��õ��˶�ά���齵ά����
complex<double> **cmatrix(long nrow, long ncol);//����һ��double�͵ĸ�����nrow��ncol�ף��õ��˶�ά���齵ά����

#endif
