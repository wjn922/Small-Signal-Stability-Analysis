//励磁系统模型
#include"head.h"

//EA模型
void exci_EA(vector<EXCI_INIT> init_exci, vector<EXCI_PARA> para_exci, double **Tgb, double **Jgb, int k)
{
	double Tr, Ta, Ta1, Te, Tf;
	double Ka, Ke, Kf, Se;

	//参数
	Tr = para_exci[k].Tr;	
	Ta = para_exci[k].Ta;	Te = para_exci[k].Te;	Tf = para_exci[k].Tf;	Ta1 = para_exci[k].Ta1;
	Ka = para_exci[k].Ka;	Ke = para_exci[k].Ke;	Kf = para_exci[k].Kf;
	Se = para_exci[k].Se;

	//形成模型
	/*
	Tgb[1][1] = Te;
	Tgb[2][2] = Ta;
	Tgb[3][3] = Ta1;
	Tgb[4][1] = Kf;
	Tgb[4][4] = -Tf;
	Tgb[5][5] = Tr;

	Jgb[1][1] = -(Ke + Se);
	Jgb[1][3] = 1.0;
	Jgb[2][2] = -1.0;
	Jgb[2][4] = -Ka;
	Jgb[2][5] = -Ka;
	Jgb[3][2] = 1.0;
	Jgb[3][3] = -1.0;
	Jgb[4][4] = 1.0;
	Jgb[5][5] = -1.0;
	*/

	Tgb[1][1] = Te;
	Tgb[2][2] = Ta + Ta1;
	Tgb[2][4] = -0.1 * Ka;
	Tgb[3][3] = Tf;
	Tgb[4][3] = 1.0;
	Tgb[4][4] = 0.1;
	Tgb[5][5] = Tr;

	Jgb[1][1] = -(Ke + Se);
	Jgb[1][2] = 1.0;
	Jgb[2][2] = -1.0;
	Jgb[2][4] = Ka;
	Jgb[3][1] = Kf;
	Jgb[3][3] = -1.0;
	Jgb[4][4] = -1.0;
	Jgb[4][5] = -1.0;
	Jgb[5][5] = -1.0;

}