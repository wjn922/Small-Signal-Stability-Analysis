//PSS系统模型
#include"head.h"

//SS模型
void pss_SS(vector<PSS_PARA> para_pss, double **Tgb, double **Jgb, int k)
{
	double Tqs, Tq, Tq1, Tq11, Tq2, Tq21, Tq3, Tq31;

	//参数
	Tqs = para_pss[k].Tqs;	Tq = para_pss[k].Tq;
	Tq1 = para_pss[k].Tq1;	Tq11 = para_pss[k].Tq11;
	Tq2 = para_pss[k].Tq2;	Tq21 = para_pss[k].Tq21;
	Tq3 = para_pss[k].Tq3;	Tq31 = para_pss[k].Tq31;

	//形成模型
	Tgb[1][1] = Tqs;
	Tgb[2][1] = Tq;
	Tgb[2][2] = -Tq;
	Tgb[3][2] = -Tq11;
	Tgb[3][3] = Tq1;
	Tgb[4][3] = -Tq21;
	Tgb[4][4] = Tq2;
	Tgb[5][4] = -Tq31;
	Tgb[5][5] = Tq3;

	Jgb[1][1] = -1.0;
	Jgb[2][2] = 1.0;
	Jgb[3][2] = -1.0;
	Jgb[3][3] = 1.0;
	Jgb[4][3] = -1.0;
	Jgb[4][4] = 1.0;
	Jgb[5][4] = -1.0;
	Jgb[5][5] = 1.0;
}