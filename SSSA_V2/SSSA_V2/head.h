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
#define base 100					//标幺值基准

struct NODE				//节点数据（包括所有节点）
{
	int no;				//节点号
	string name;		//节点名字
	double vol;			//节点基准电压（kV）
};

struct NODE_DATA		//节点潮流数据（包括所有节点）
{
	int no;				//这里的no是潮流结果中 优化编号后的
	int type;			//节点类型 -1=发电机节点 0=联络节点 1=负荷节点
	string name;
	double vol_base, vol_pu, vol, vol_ang;
	double P, Q;		//流出为正
};


struct GENE_INIT		//发电机初值
{
	double Delta0, Vd0, Vq0, Id0, Iq0;
	double Efq0, Eq10, Eq20, Ed10, Ed20;
};


struct GENE_PARA		//发电机参数
{
	//同步机参数
	string card_type;		//发电机卡类型
	string gene_model;		//发电机模型 有four_winding, three_winding, two_axle, Eq_change, Eq_constant, classic
	string name;			//发电机名称
	double vol;				//母线基准电压(kV)
	double capacity;		//电机额定容量(MVA)
	double power_factor;	//电机功率因数
	string type;			//电机类型
	double Ek;				//发电机动能(MW s）
	double Sh;				//基准容量（MVA）

	double Xd, Xd1, Xd2;
	double Xq, Xq1, Xq2;
	double Ra;
	double Td01, Td02;
	double Tq01, Tq02;
	double Tj;
	double D;				//阻尼系数
	double XL;				//发电机定子漏抗
	double A, B, N;			//发电机饱和系数
};//XL, A, B, N本程序中未用到


struct LOAD_PARA		//负荷参数
{
	string card_type;		//负荷卡类型
	string load_model;		//负荷卡模型
	string name;
	double vol;
	double ap, bp, cp;
	double aq, bq, cq;
};




//////////////////////////函数
void nrerror(char error_text[]);//数值计算错误时，输出报告
int **imatrix(long nrow, long ncol);//分配一个int型的矩阵nrow×ncol阶，用到了二维数组降维处理
double **dmatrix(long nrow, long ncol);//分配一个double型的矩阵nrow×ncol阶，用到了二维数组降维处理
complex<double> **cmatrix(long nrow, long ncol);//分配一个double型的复矩阵nrow×ncol阶，用到了二维数组降维处理

#endif
