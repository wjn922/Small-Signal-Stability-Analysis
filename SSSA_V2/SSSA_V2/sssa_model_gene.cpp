//同步机模型
#include"head.h"

//////////////////////////////////////MF/MC卡//////////////////////////////////////////
// MF 四绕组模型
void gene_four_winding(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,		
	double **Tgb, double **Jgb, int k)
{
	double Delta0, Vd0, Vq0, Id0, Iq0, Efq0, Eq10, Eq20, Ed10, Ed20;		//发电机初值
	double Xd, Xd1, Xd2, Xq, Xq1, Xq2, Ra, Td01, Td02, Tq01, Tq02, Tj, D;	//发电机参数

	//开始
	double rate = base / para_gene[k].Sh;

	//参数
	Xd = rate * para_gene[k].Xd;  Xd1 = rate * para_gene[k].Xd1;  Xd2 = rate * para_gene[k].Xd2;
	Xq = rate * para_gene[k].Xq;  Xq1 = rate * para_gene[k].Xq1;  Xq2 = rate * para_gene[k].Xq2;
	Ra = rate * para_gene[k].Ra;

	Td01 = para_gene[k].Td01 / rate;  Td02 = para_gene[k].Td02 / rate;
	Tq01 = para_gene[k].Tq01 / rate;  Tq02 = para_gene[k].Tq02 / rate;

	Tj = para_gene[k].Tj / rate;  D = para_gene[k].D / rate;

	//初值
	Delta0 = init_gene[k].Delta0;
	Vd0 = init_gene[k].Vd0;      Vq0 = init_gene[k].Vq0;
	Id0 = init_gene[k].Id0;	     Iq0 = init_gene[k].Iq0;	  Efq0 = init_gene[k].Efq0;
	Eq10 = init_gene[k].Eq10;    Eq20 = init_gene[k].Eq20;
	Ed10 = init_gene[k].Ed10;    Ed20 = init_gene[k].Ed20;

	//形成模型
	Tgb[1][1] = 1.0;
	Tgb[2][2] = Tj;
	Tgb[3][3] = Td01;
	Tgb[4][4] = Td02;
	Tgb[5][5] = Tq01;
	Tgb[6][6] = Tq02;

	Jgb[1][2] = 100.0 * PI;
	Jgb[2][2] = -D;
	Jgb[2][7] = -1.0;
	Jgb[3][3] = -(Xd - Xd2) / (Xd1 - Xd2);
	Jgb[3][4] = (Xd - Xd1) / (Xd1 - Xd2);
	Jgb[4][3] = 1.0;
	Jgb[4][4] = -1.0;
	Jgb[4][8] = -(Xd1 - Xd2);
	Jgb[5][5] = -(Xq - Xq2) / (Xq1 - Xq2);
	Jgb[5][6] = (Xq - Xq1) / (Xq1 - Xq2);
	Jgb[6][5] = 1.0;
	Jgb[6][6] = -1.0;
	Jgb[6][9] = Xq1 - Xq2;

	//Jgb[7][4] = Iq0;
	//Jgb[7][6] = Id0;
	//Jgb[7][7] = -1.0;
	//Jgb[7][8] = Ed20 + (Xq2 - Xd2) * Iq0;
	//Jgb[7][9] = Eq20 + (Xq2 - Xd2) * Id0;
	Jgb[7][7] = -1.0;
	Jgb[7][8] = Vd0;
	Jgb[7][9] = Vq0;
	Jgb[7][10] = Id0;
	Jgb[7][11] = Iq0;
	Jgb[8][6] = 1.0;	//△Ed2 - Ra△Id + Xq2△Iq - △Vd = 0;
	Jgb[8][8] = -Ra;
	Jgb[8][9] = Xq2;
	Jgb[8][10] = -1.0;
	Jgb[9][4] = 1.0;	//△Eq2 - Ra△Iq - Xd2△Id - △Vq = 0;
	Jgb[9][8] = -Xd2;
	Jgb[9][9] = -Ra;
	Jgb[9][11] = -1.0;
	Jgb[10][1] = Vq0;	//sin(delta0)△Vx - cos(delta0)△Vy + Vq0△omega - △Vd = 0
	Jgb[10][10] = -1.0;
	Jgb[11][1] = -Vd0;	//cos(delta0)△Vx + sin(delta0)△Vy - Vd0△omega - △Vq = 0
	Jgb[11][11] = -1.0;
	Jgb[12][1] = Iq0;	//sin(delta0)△Ix - cos(delta0)△Iy + Iq0△omega - △Id = 0
	Jgb[12][8] = -1.0;
	Jgb[12][12] = sin(Delta0);
	Jgb[12][13] = -cos(Delta0);
	Jgb[13][1] = -Id0;	//cos(delta0)△Ix + sin(delta0)△Iy - Id0△omega - △Iq = 0
	Jgb[13][9] = -1.0;
	Jgb[13][12] = cos(Delta0);
	Jgb[13][13] = sin(Delta0);
}


//MF 三绕组模型
void gene_three_winding(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,		
	double **Tgb, double **Jgb, int k)
{
	double Delta0, Vd0, Vq0, Id0, Iq0, Efq0, Eq10, Eq20, Ed10, Ed20;		//发电机初值
	double Xd, Xd1, Xd2, Xq, Xq1, Xq2, Ra, Td01, Td02, Tq01, Tq02, Tj, D;	//发电机参数

	//开始
	double rate = base / para_gene[k].Sh;

	//参数
	Xd = rate * para_gene[k].Xd;  Xd1 = rate * para_gene[k].Xd1;  Xd2 = rate * para_gene[k].Xd2;
	Xq = rate * para_gene[k].Xq;  Xq1 = rate * para_gene[k].Xq1;  Xq2 = rate * para_gene[k].Xq2;
	Ra = rate * para_gene[k].Ra;

	Td01 = para_gene[k].Td01 / rate;  Td02 = para_gene[k].Td02 / rate;
	Tq01 = para_gene[k].Tq01 / rate;  Tq02 = para_gene[k].Tq02 / rate;

	Tj = para_gene[k].Tj / rate;  D = para_gene[k].D / rate;


	//初值
	Delta0 = init_gene[k].Delta0;
	Vd0 = init_gene[k].Vd0;      Vq0 = init_gene[k].Vq0;
	Id0 = init_gene[k].Id0;		 Iq0 = init_gene[k].Iq0;	  Efq0 = init_gene[k].Efq0;
	Eq10 = init_gene[k].Eq10;    Eq20 = init_gene[k].Eq20;
	Ed10 = init_gene[k].Ed10;    Ed20 = init_gene[k].Ed20;

	//形成模型
	Tgb[1][1] = 1.0;
	Tgb[2][2] = Tj;
	Tgb[3][3] = Td01;
	Tgb[4][4] = Td02;
	Tgb[5][5] = Tq02;

	Jgb[1][2] = 100.0 * PI;
	Jgb[2][2] = -D;
	Jgb[2][6] = -1.0;
	Jgb[3][3] = -(Xd - Xd2) / (Xd1 - Xd2);
	Jgb[3][4] = (Xd - Xd1) / (Xd1 - Xd2);
	Jgb[4][3] = 1.0;
	Jgb[4][4] = -1.0;
	Jgb[4][7] = -(Xd1 - Xd2);
	Jgb[5][5] = -1.0;
	Jgb[5][8] = Xq1 - Xq2;

	//Jgb[6][4] = Iq0;
	//Jgb[6][5] = Id0;
	//Jgb[6][6] = -1.0;
	//Jgb[6][7] = Ed20 + (Xq2 - Xd2) * Iq0;
	//Jgb[6][8] = Eq20 + (Xq2 - Xd2) * Id0;
	Jgb[6][6] = -1.0;
	Jgb[6][7] = Vd0;
	Jgb[6][8] = Vq0;
	Jgb[6][9] = Id0;
	Jgb[6][10] = Iq0;
	Jgb[7][5] = 1.0;	//△Ed2 - Ra△Id + Xq2△Iq - △Vd = 0;
	Jgb[7][7] = -Ra;
	Jgb[7][8] = Xq2;
	Jgb[7][9] = -1.0;
	Jgb[8][4] = 1.0;	//△Eq2 - Ra△Iq - Xd2△Id - △Vq = 0;
	Jgb[8][7] = -Xd2;
	Jgb[8][8] = -Ra;
	Jgb[8][10] = -1.0;
	Jgb[9][1] = Vq0;	//sin(delta0)△Vx - cos(delta0)△Vy + Vq0△omega - △Vd = 0
	Jgb[9][9] = -1.0;
	Jgb[10][1] = -Vd0;	//cos(delta0)△Vx + sin(delta0)△Vy - Vd0△omega - △Vq = 0
	Jgb[10][10] = -1.0;
	Jgb[11][1] = Iq0;	//sin(delta0)△Ix - cos(delta0)△Iy + Iq0△omega - △Id = 0
	Jgb[11][7] = -1.0;
	Jgb[11][11] = sin(Delta0);
	Jgb[11][12] = -cos(Delta0);
	Jgb[12][1] = -Id0;	//cos(delta0)△Ix + sin(delta0)△Iy - Id0△omega - △Iq = 0
	Jgb[12][8] = -1.0;
	Jgb[12][11] = cos(Delta0);
	Jgb[12][12] = sin(Delta0);
}

//MF 双轴模型
void gene_two_axle(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,
	double **Tgb, double **Jgb, int k)
{
	double Delta0, Vd0, Vq0, Id0, Iq0, Efq0, Eq10, Eq20, Ed10, Ed20;		//发电机初值
	double Xd, Xd1, Xd2, Xq, Xq1, Xq2, Ra, Td01, Td02, Tq01, Tq02, Tj, D;	//发电机参数

	//开始
	double rate = base / para_gene[k].Sh;

	//参数
	Xd = rate * para_gene[k].Xd;  Xd1 = rate * para_gene[k].Xd1;  Xd2 = rate * para_gene[k].Xd2;
	Xq = rate * para_gene[k].Xq;  Xq1 = rate * para_gene[k].Xq1;  Xq2 = rate * para_gene[k].Xq2;
	Ra = rate * para_gene[k].Ra;

	Td01 = para_gene[k].Td01 / rate;  Td02 = para_gene[k].Td02 / rate;
	Tq01 = para_gene[k].Tq01 / rate;  Tq02 = para_gene[k].Tq02 / rate;

	Tj = para_gene[k].Tj / rate;  D = para_gene[k].D / rate;


	//初值
	Delta0 = init_gene[k].Delta0;
	Vd0 = init_gene[k].Vd0;      Vq0 = init_gene[k].Vq0;
	Id0 = init_gene[k].Id0;	     Iq0 = init_gene[k].Iq0;	  Efq0 = init_gene[k].Efq0;
	Eq10 = init_gene[k].Eq10;    Eq20 = init_gene[k].Eq20;
	Ed10 = init_gene[k].Ed10;    Ed20 = init_gene[k].Ed20;

	//形成模型
	Tgb[1][1] = 1.0;
	Tgb[2][2] = Tj;
	Tgb[3][3] = Td01;
	Tgb[4][4] = Tq01;

	Jgb[1][2] = 100.0 * PI;
	Jgb[2][2] = -D;
	Jgb[2][5] = -1.0;
	Jgb[3][3] = -1.0;
	Jgb[3][6] = -(Xd - Xd1);
	Jgb[4][4] = -1.0;
	Jgb[4][7] = Xq - Xq1;

	//Jgb[5][3] = Iq0;
	//Jgb[5][4] = Id0;
	//Jgb[5][5] = -1.0;
	//Jgb[5][6] = Ed10 + (Xq1 - Xd1) * Iq0;
	//Jgb[5][7] = Eq10 + (Xq1 - Xd1) * Id0;
	Jgb[5][5] = -1.0;
	Jgb[5][6] = Vd0;
	Jgb[5][7] = Vq0;
	Jgb[5][8] = Id0;
	Jgb[5][9] = Iq0;
	Jgb[6][4] = 1.0;	//△Ed1 - Ra△Id + Xq1△Iq - △Vd = 0;
	Jgb[6][6] = -Ra;
	Jgb[6][7] = Xq1;
	Jgb[6][8] = -1.0;
	Jgb[7][3] = 1.0;	//△Eq1 - Ra△Iq - Xd1△Id - △Vq = 0;
	Jgb[7][6] = -Xd1;
	Jgb[7][7] = -Ra;
	Jgb[7][9] = -1.0;
	Jgb[8][1] = Vq0;	//sin(delta0)△Vx - cos(delta0)△Vy + Vq0△omega - △Vd = 0
	Jgb[8][8] = -1.0;
	Jgb[9][1] = -Vd0;	//cos(delta0)△Vx + sin(delta0)△Vy - Vd0△omega - △Vq = 0
	Jgb[9][9] = -1.0;
	Jgb[10][1] = Iq0;	//sin(delta0)△Ix - cos(delta0)△Iy + Iq0△omega - △Id = 0
	Jgb[10][6] = -1.0;
	Jgb[10][10] = sin(Delta0);
	Jgb[10][11] = -cos(Delta0);
	Jgb[11][1] = -Id0;	//cos(delta0)△Ix + sin(delta0)△Iy - Id0△omega - △Iq = 0
	Jgb[11][7] = -1.0;
	Jgb[11][10] = cos(Delta0);
	Jgb[11][11] = sin(Delta0);
}

//MF Eq'变化
void gene_Eq_change(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,
	double **Tgb, double **Jgb, int k)
{
	double Delta0, Vd0, Vq0, Id0, Iq0, Efq0, Eq10, Eq20, Ed10, Ed20;		//发电机初值
	double Xd, Xd1, Xd2, Xq, Xq1, Xq2, Ra, Td01, Td02, Tq01, Tq02, Tj, D;	//发电机参数

	//开始
	double rate = base / para_gene[k].Sh;

	//参数
	Xd = rate * para_gene[k].Xd;  Xd1 = rate * para_gene[k].Xd1;  Xd2 = rate * para_gene[k].Xd2;
	Xq = rate * para_gene[k].Xq;  Xq1 = rate * para_gene[k].Xq1;  Xq2 = rate * para_gene[k].Xq2;
	Ra = rate * para_gene[k].Ra;

	Td01 = para_gene[k].Td01 / rate;  Td02 = para_gene[k].Td02 / rate;
	Tq01 = para_gene[k].Tq01 / rate;  Tq02 = para_gene[k].Tq02 / rate;

	Tj = para_gene[k].Tj / rate;  D = para_gene[k].D / rate;


	//初值
	Delta0 = init_gene[k].Delta0;
	Vd0 = init_gene[k].Vd0;      Vq0 = init_gene[k].Vq0;
	Id0 = init_gene[k].Id0;	     Iq0 = init_gene[k].Iq0;	  Efq0 = init_gene[k].Efq0;
	Eq10 = init_gene[k].Eq10;    Eq20 = init_gene[k].Eq20;
	Ed10 = init_gene[k].Ed10;    Ed20 = init_gene[k].Ed20;

	//形成模型
	Tgb[1][1] = 1.0;
	Tgb[2][2] = Tj;
	Tgb[3][3] = Td01;

	Jgb[1][2] = 100.0 * PI;
	Jgb[2][2] = -D;
	Jgb[2][4] = -1.0;
	Jgb[3][3] = -1.0;
	Jgb[3][5] = -(Xd - Xd1);

	//Jgb[4][3] = Iq0;
	//Jgb[4][4] = -1.0;
	//Jgb[4][5] = (Xq - Xd1) * Iq0;
	//Jgb[4][6] = Eq10 + (Xq - Xd1) * Id0;
	Jgb[4][4] = -1.0;
	Jgb[4][5] = Vd0;
	Jgb[4][6] = Vq0;
	Jgb[4][7] = Id0;
	Jgb[4][8] = Iq0;
	Jgb[5][5] = -Ra;	//- Ra△Id + Xq△Iq - △Vd = 0;
	Jgb[5][6] = Xq;
	Jgb[5][7] = -1.0;
	Jgb[6][3] = 1.0;	//△Eq1 - Ra△Iq - Xd1△Id - △Vq = 0;
	Jgb[6][5] = -Xd1;
	Jgb[6][6] = -Ra;
	Jgb[6][8] = -1.0;
	Jgb[7][1] = Vq0;	//sin(delta0)△Vx - cos(delta0)△Vy + Vq0△omega - △Vd = 0
	Jgb[7][7] = -1.0;
	Jgb[8][1] = -Vd0;	//cos(delta0)△Vx + sin(delta0)△Vy - Vd0△omega - △Vq = 0
	Jgb[8][8] = -1.0;
	Jgb[9][1] = Iq0;	//sin(delta0)△Ix - cos(delta0)△Iy + Iq0△omega - △Id = 0
	Jgb[9][5] = -1.0;
	Jgb[9][9] = sin(Delta0);
	Jgb[9][10] = -cos(Delta0);
	Jgb[10][1] = -Id0;	//cos(delta0)△Ix + sin(delta0)△Iy - Id0△omega - △Iq = 0
	Jgb[10][6] = -1.0;
	Jgb[10][9] = cos(Delta0);
	Jgb[10][10] = sin(Delta0);
}


//MF Eq'恒定
void gene_Eq_constant(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,
	double **Tgb, double **Jgb, int k)
{
	double Delta0, Vd0, Vq0, Id0, Iq0, Efq0, Eq10, Eq20, Ed10, Ed20;		//发电机初值
	double Xd, Xd1, Xd2, Xq, Xq1, Xq2, Ra, Td01, Td02, Tq01, Tq02, Tj, D;	//发电机参数

	//开始
	double rate = base / para_gene[k].Sh;

	//参数
	Xd = rate * para_gene[k].Xd;  Xd1 = rate * para_gene[k].Xd1;  Xd2 = rate * para_gene[k].Xd2;
	Xq = rate * para_gene[k].Xq;  Xq1 = rate * para_gene[k].Xq1;  Xq2 = rate * para_gene[k].Xq2;
	Ra = rate * para_gene[k].Ra;

	Td01 = para_gene[k].Td01 / rate;  Td02 = para_gene[k].Td02 / rate;
	Tq01 = para_gene[k].Tq01 / rate;  Tq02 = para_gene[k].Tq02 / rate;

	Tj = para_gene[k].Tj / rate;  D = para_gene[k].D / rate;


	//初值
	Delta0 = init_gene[k].Delta0;
	Vd0 = init_gene[k].Vd0;      Vq0 = init_gene[k].Vq0;
	Id0 = init_gene[k].Id0;	     Iq0 = init_gene[k].Iq0;	  Efq0 = init_gene[k].Efq0;
	Eq10 = init_gene[k].Eq10;    Eq20 = init_gene[k].Eq20;
	Ed10 = init_gene[k].Ed10;    Ed20 = init_gene[k].Ed20;

	//形成模型
	Tgb[1][1] = 1.0;
	Tgb[2][2] = Tj;

	Jgb[1][2] = 100.0 * PI;
	Jgb[2][2] = -D;
	Jgb[2][3] = -1.0;

	//Jgb[3][3] = -1.0;
	//Jgb[3][4] = (Xq - Xd1) * Iq0;
	//Jgb[3][5] = Eq10 + (Xq - Xd1) * Id0;
	Jgb[3][3] = -1.0;
	Jgb[3][4] = Vd0;
	Jgb[3][5] = Vq0;
	Jgb[3][6] = Id0;
	Jgb[3][7] = Iq0;
	Jgb[4][4] = -Ra;	//- Ra△Id + Xq△Iq - △Vd = 0;
	Jgb[4][5] = Xq;
	Jgb[4][6] = -1.0;
	Jgb[5][4] = -Xd1;	//- Ra△Iq - Xd1△Id - △Vq = 0;
	Jgb[5][5] = -Ra;
	Jgb[5][7] = -1.0;
	Jgb[6][1] = Vq0;	//sin(delta0)△Vx - cos(delta0)△Vy + Vq0△omega - △Vd = 0
	Jgb[6][6] = -1.0;
	Jgb[7][1] = -Vd0;	//cos(delta0)△Vx + sin(delta0)△Vy - Vd0△omega - △Vq = 0
	Jgb[7][7] = -1.0;
	Jgb[8][1] = Iq0;	//sin(delta0)△Ix - cos(delta0)△Iy + Iq0△omega - △Id = 0
	Jgb[8][4] = -1.0;
	Jgb[8][8] = sin(Delta0);
	Jgb[8][9] = -cos(Delta0);
	Jgb[9][1] = -Id0;	//cos(delta0)△Ix + sin(delta0)△Iy - Id0△omega - △Iq = 0
	Jgb[9][5] = -1.0;
	Jgb[9][8] = cos(Delta0);
	Jgb[9][9] = sin(Delta0);
}

//MC 经典模型
void gene_classic(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,
	double **Tgb, double **Jgb, int k)
{
	double Delta0, Vd0, Vq0, Id0, Iq0, Efq0, Eq10, Eq20, Ed10, Ed20;		//发电机初值
	double Xd, Xd1, Xd2, Xq, Xq1, Xq2, Ra, Td01, Td02, Tq01, Tq02, Tj, D;	//发电机参数

	//开始
	double rate = base / para_gene[k].Sh;

	//参数 注意Xq = Xd1
	Xd = rate * para_gene[k].Xd;  Xd1 = rate * para_gene[k].Xd1;  Xd2 = rate * para_gene[k].Xd2;
	Xq = rate * para_gene[k].Xd1;  Xq1 = rate * para_gene[k].Xq1;  Xq2 = rate * para_gene[k].Xq2;
	Ra = rate * para_gene[k].Ra;

	Td01 = para_gene[k].Td01 / rate;  Td02 = para_gene[k].Td02 / rate;
	Tq01 = para_gene[k].Tq01 / rate;  Tq02 = para_gene[k].Tq02 / rate;

	Tj = para_gene[k].Tj / rate;  D = para_gene[k].D / rate;


	//初值
	Delta0 = init_gene[k].Delta0;
	Vd0 = init_gene[k].Vd0;      Vq0 = init_gene[k].Vq0;
	Id0 = init_gene[k].Id0;	     Iq0 = init_gene[k].Iq0;	  Efq0 = init_gene[k].Efq0;
	Eq10 = init_gene[k].Eq10;    Eq20 = init_gene[k].Eq20;
	Ed10 = init_gene[k].Ed10;    Ed20 = init_gene[k].Ed20;

	//形成模型
	Tgb[1][1] = 1.0;
	Tgb[2][2] = Tj;

	Jgb[1][2] = 100.0 * PI;
	Jgb[2][2] = -D;
	Jgb[2][3] = -1.0;

	Jgb[3][3] = -1.0;
	//Jgb[3][5] = Eq10;
	Jgb[3][4] = Vd0;
	Jgb[3][5] = Vq0;
	Jgb[3][6] = Id0;
	Jgb[3][7] = Iq0;
	Jgb[4][4] = -Ra;	//- Ra△Id + Xd1△Iq - △Vd = 0;
	Jgb[4][5] = Xd1;
	Jgb[4][6] = -1.0;
	Jgb[5][4] = -Xd1;	//- Ra△Iq - Xd1△Id - △Vq = 0;
	Jgb[5][5] = -Ra;
	Jgb[5][7] = -1.0;
	Jgb[6][1] = Vq0;	//sin(delta0)△Vx - cos(delta0)△Vy + Vq0△omega - △Vd = 0
	Jgb[6][6] = -1.0;
	Jgb[7][1] = -Vd0;	//cos(delta0)△Vx + sin(delta0)△Vy - Vd0△omega - △Vq = 0
	Jgb[7][7] = -1.0;
	Jgb[8][1] = Iq0;	//sin(delta0)△Ix - cos(delta0)△Iy + Iq0△omega - △Id = 0
	Jgb[8][4] = -1.0;
	Jgb[8][8] = sin(Delta0);
	Jgb[8][9] = -cos(Delta0);
	Jgb[9][1] = -Id0;	//cos(delta0)△Ix + sin(delta0)△Iy - Id0△omega - △Iq = 0
	Jgb[9][5] = -1.0;
	Jgb[9][8] = cos(Delta0);
	Jgb[9][9] = sin(Delta0);
}


////////////////////////////////////////////MG卡///////////////////////////////////////////////
//6型，Eq',Eq'',Ed',Ed''变化的6阶模型
void gene_type6(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,
	double **Tgb, double **Jgb, int k)
{
	double Delta0, Vd0, Vq0, Id0, Iq0, Efq0, Eq10, Eq20, Ed10, Ed20;		//发电机初值
	double Xd, Xd1, Xd2, Xq, Xq1, Xq2, Ra, Td01, Td02, Tq01, Tq02, Tj, D;	//发电机参数
	double A, B, N, KG;

	//开始
	double rate = base / para_gene[k].Sh;

	//参数
	Xd = rate * para_gene[k].Xd;  Xd1 = rate * para_gene[k].Xd1;  Xd2 = rate * para_gene[k].Xd2;
	Xq = rate * para_gene[k].Xq;  Xq1 = rate * para_gene[k].Xq1;  Xq2 = rate * para_gene[k].Xq2;
	Ra = rate * para_gene[k].Ra;

	Td01 = para_gene[k].Td01 / rate;  Td02 = para_gene[k].Td02 / rate;
	Tq01 = para_gene[k].Tq01 / rate;  Tq02 = para_gene[k].Tq02 / rate;

	Tj = para_gene[k].Tj / rate;  D = para_gene[k].D / rate;

	A = para_gene[k].A;		B = para_gene[k].B;		N = para_gene[k].N;


	//初值
	Delta0 = init_gene[k].Delta0;
	Vd0 = init_gene[k].Vd0;      Vq0 = init_gene[k].Vq0;
	Id0 = init_gene[k].Id0;	     Iq0 = init_gene[k].Iq0;	  Efq0 = init_gene[k].Efq0;
	Eq10 = init_gene[k].Eq10;    Eq20 = init_gene[k].Eq20;
	Ed10 = init_gene[k].Ed10;    Ed20 = init_gene[k].Ed20;

	if (B == 0)	KG = 1.0;
	else
		KG = 1 + B * pow(Eq10, N - 1) / A;

	//形成模型
	Tgb[1][1] = 1.0;
	Tgb[2][2] = Tj;
	Tgb[3][3] = Td01;
	Tgb[4][3] = -Td02;
	Tgb[4][4] = Td02;
	Tgb[5][5] = Tq01;
	Tgb[6][5] = -Tq02;
	Tgb[6][6] = Tq02;

	Jgb[1][2] = 100.0 * PI;
	Jgb[2][2] = -D;
	Jgb[2][7] = -1.0;
	Jgb[3][3] = -KG;
	Jgb[3][8] = -(Xd - Xd1);
	Jgb[4][3] = 1.0;
	Jgb[4][4] = -1.0;
	Jgb[4][8] = -(Xd1 - Xd2);
	Jgb[5][5] = -1.0;
	Jgb[5][9] = Xq - Xq1;
	Jgb[6][5] = 1.0;
	Jgb[6][6] = -1.0;
	Jgb[6][9] = Xq1 - Xq2;


	Jgb[7][7] = -1.0;
	Jgb[7][8] = Vd0;
	Jgb[7][9] = Vq0;
	Jgb[7][10] = Id0;
	Jgb[7][11] = Iq0;
	Jgb[8][6] = 1.0;	//△Ed2 - Ra△Id + Xq2△Iq - △Vd = 0;
	Jgb[8][8] = -Ra;
	Jgb[8][9] = Xq2;
	Jgb[8][10] = -1.0;
	Jgb[9][4] = 1.0;	//△Eq2 - Ra△Iq - Xd2△Id - △Vq = 0;
	Jgb[9][8] = -Xd2;
	Jgb[9][9] = -Ra;
	Jgb[9][11] = -1.0;
	Jgb[10][1] = Vq0;	//sin(delta0)△Vx - cos(delta0)△Vy + Vq0△omega - △Vd = 0
	Jgb[10][10] = -1.0;
	Jgb[11][1] = -Vd0;	//cos(delta0)△Vx + sin(delta0)△Vy - Vd0△omega - △Vq = 0
	Jgb[11][11] = -1.0;
	Jgb[12][1] = Iq0;	//sin(delta0)△Ix - cos(delta0)△Iy + Iq0△omega - △Id = 0
	Jgb[12][8] = -1.0;
	Jgb[12][12] = sin(Delta0);
	Jgb[12][13] = -cos(Delta0);
	Jgb[13][1] = -Id0;	//cos(delta0)△Ix + sin(delta0)△Iy - Id0△omega - △Iq = 0
	Jgb[13][9] = -1.0;
	Jgb[13][12] = cos(Delta0);
	Jgb[13][13] = sin(Delta0);
}


//5型，Eq',Ed'变化的4阶模型
void gene_type5(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,
	double **Tgb, double **Jgb, int k)
{
	double Delta0, Vd0, Vq0, Id0, Iq0, Efq0, Eq10, Eq20, Ed10, Ed20;		//发电机初值
	double Xd, Xd1, Xd2, Xq, Xq1, Xq2, Ra, Td01, Td02, Tq01, Tq02, Tj, D;	//发电机参数
	double A, B, N, KG;														//饱和系数

	//开始
	double rate = base / para_gene[k].Sh;

	//参数
	Xd = rate * para_gene[k].Xd;  Xd1 = rate * para_gene[k].Xd1;  Xd2 = rate * para_gene[k].Xd2;
	Xq = rate * para_gene[k].Xq;  Xq1 = rate * para_gene[k].Xq1;  Xq2 = rate * para_gene[k].Xq2;
	Ra = rate * para_gene[k].Ra;

	Td01 = para_gene[k].Td01 / rate;  Td02 = para_gene[k].Td02 / rate;
	Tq01 = para_gene[k].Tq01 / rate;  Tq02 = para_gene[k].Tq02 / rate;

	Tj = para_gene[k].Tj / rate;  D = para_gene[k].D / rate;

	A = para_gene[k].A;		B = para_gene[k].B;		N = para_gene[k].N;

	//初值
	Delta0 = init_gene[k].Delta0;
	Vd0 = init_gene[k].Vd0;      Vq0 = init_gene[k].Vq0;
	Id0 = init_gene[k].Id0;	     Iq0 = init_gene[k].Iq0;	  Efq0 = init_gene[k].Efq0;
	Eq10 = init_gene[k].Eq10;    Eq20 = init_gene[k].Eq20;
	Ed10 = init_gene[k].Ed10;    Ed20 = init_gene[k].Ed20;

	if (B == 0)	KG = 1.0;
	else
		KG = 1 + B * pow(Eq10, N - 1) / A;

	//形成模型
	Tgb[1][1] = 1.0;
	Tgb[2][2] = Tj;
	Tgb[3][3] = Td01;
	Tgb[4][4] = Tq01;

	Jgb[1][2] = 100.0 * PI;
	Jgb[2][2] = -D;
	Jgb[2][5] = -1.0;
	Jgb[3][3] = -KG;
	Jgb[3][6] = -(Xd - Xd1);
	Jgb[4][4] = -1.0;
	Jgb[4][7] = Xq - Xq1;

	Jgb[5][5] = -1.0;
	Jgb[5][6] = Vd0;
	Jgb[5][7] = Vq0;
	Jgb[5][8] = Id0;
	Jgb[5][9] = Iq0;
	Jgb[6][4] = 1.0;	//△Ed1 - Ra△Id + Xq1△Iq - △Vd = 0;
	Jgb[6][6] = -Ra;
	Jgb[6][7] = Xq1;
	Jgb[6][8] = -1.0;
	Jgb[7][3] = 1.0;	//△Eq1 - Ra△Iq - Xd1△Id - △Vq = 0;
	Jgb[7][6] = -Xd1;
	Jgb[7][7] = -Ra;
	Jgb[7][9] = -1.0;
	Jgb[8][1] = Vq0;	//sin(delta0)△Vx - cos(delta0)△Vy + Vq0△omega - △Vd = 0
	Jgb[8][8] = -1.0;
	Jgb[9][1] = -Vd0;	//cos(delta0)△Vx + sin(delta0)△Vy - Vd0△omega - △Vq = 0
	Jgb[9][9] = -1.0;
	Jgb[10][1] = Iq0;	//sin(delta0)△Ix - cos(delta0)△Iy + Iq0△omega - △Id = 0
	Jgb[10][6] = -1.0;
	Jgb[10][10] = sin(Delta0);
	Jgb[10][11] = -cos(Delta0);
	Jgb[11][1] = -Id0;	//cos(delta0)△Ix + sin(delta0)△Iy - Id0△omega - △Iq = 0
	Jgb[11][7] = -1.0;
	Jgb[11][10] = cos(Delta0);
	Jgb[11][11] = sin(Delta0);
}


//3型，Eq',Eq'',Ed''变化的5阶模型
void gene_type3(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,
	double **Tgb, double **Jgb, int k)
{
	double Delta0, Vd0, Vq0, Id0, Iq0, Efq0, Eq10, Eq20, Ed10, Ed20;		//发电机初值
	double Xd, Xd1, Xd2, Xq, Xq1, Xq2, Ra, Td01, Td02, Tq01, Tq02, Tj, D;	//发电机参数
	double A, B, N, KG;

	//开始
	double rate = base / para_gene[k].Sh;

	//参数
	Xd = rate * para_gene[k].Xd;  Xd1 = rate * para_gene[k].Xd1;  Xd2 = rate * para_gene[k].Xd2;
	Xq = rate * para_gene[k].Xq;  Xq1 = rate * para_gene[k].Xq1;  Xq2 = rate * para_gene[k].Xq2;
	Ra = rate * para_gene[k].Ra;

	Td01 = para_gene[k].Td01 / rate;  Td02 = para_gene[k].Td02 / rate;
	Tq01 = para_gene[k].Tq01 / rate;  Tq02 = para_gene[k].Tq02 / rate;

	Tj = para_gene[k].Tj / rate;  D = para_gene[k].D / rate;

	A = para_gene[k].A;		B = para_gene[k].B;		N = para_gene[k].N;


	//初值
	Delta0 = init_gene[k].Delta0;
	Vd0 = init_gene[k].Vd0;      Vq0 = init_gene[k].Vq0;
	Id0 = init_gene[k].Id0;		 Iq0 = init_gene[k].Iq0;	  Efq0 = init_gene[k].Efq0;
	Eq10 = init_gene[k].Eq10;    Eq20 = init_gene[k].Eq20;
	Ed10 = init_gene[k].Ed10;    Ed20 = init_gene[k].Ed20;

	if (B == 0)	KG = 1.0;
	else
		KG = 1 + B * pow(Eq10, N - 1) / A;

	//形成模型
	Tgb[1][1] = 1.0;
	Tgb[2][2] = Tj;
	Tgb[3][3] = Td01;
	Tgb[4][3] = -Td02;
	Tgb[4][4] = Td02;
	Tgb[5][5] = Tq02;

	Jgb[1][2] = 100.0 * PI;
	Jgb[2][2] = -D;
	Jgb[2][6] = -1.0;
	Jgb[3][3] = -KG;
	Jgb[3][7] = -(Xd - Xd1);
	Jgb[4][3] = 1.0;
	Jgb[4][4] = -1.0;
	Jgb[4][7] = -(Xd1 - Xd2);
	Jgb[5][5] = -1.0;
	Jgb[5][8] = Xq1 - Xq2;

	Jgb[6][6] = -1.0;
	Jgb[6][7] = Vd0;
	Jgb[6][8] = Vq0;
	Jgb[6][9] = Id0;
	Jgb[6][10] = Iq0;
	Jgb[7][5] = 1.0;	//△Ed2 - Ra△Id + Xq2△Iq - △Vd = 0;
	Jgb[7][7] = -Ra;
	Jgb[7][8] = Xq2;
	Jgb[7][9] = -1.0;
	Jgb[8][4] = 1.0;	//△Eq2 - Ra△Iq - Xd2△Id - △Vq = 0;
	Jgb[8][7] = -Xd2;
	Jgb[8][8] = -Ra;
	Jgb[8][10] = -1.0;
	Jgb[9][1] = Vq0;	//sin(delta0)△Vx - cos(delta0)△Vy + Vq0△omega - △Vd = 0
	Jgb[9][9] = -1.0;
	Jgb[10][1] = -Vd0;	//cos(delta0)△Vx + sin(delta0)△Vy - Vd0△omega - △Vq = 0
	Jgb[10][10] = -1.0;
	Jgb[11][1] = Iq0;	//sin(delta0)△Ix - cos(delta0)△Iy + Iq0△omega - △Id = 0
	Jgb[11][7] = -1.0;
	Jgb[11][11] = sin(Delta0);
	Jgb[11][12] = -cos(Delta0);
	Jgb[12][1] = -Id0;	//cos(delta0)△Ix + sin(delta0)△Iy - Id0△omega - △Iq = 0
	Jgb[12][8] = -1.0;
	Jgb[12][11] = cos(Delta0);
	Jgb[12][12] = sin(Delta0);
}



//2型，Eq'变化的3阶模型
void gene_type2(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,
	double **Tgb, double **Jgb, int k)
{
	double Delta0, Vd0, Vq0, Id0, Iq0, Efq0, Eq10, Eq20, Ed10, Ed20;		//发电机初值
	double Xd, Xd1, Xd2, Xq, Xq1, Xq2, Ra, Td01, Td02, Tq01, Tq02, Tj, D;	//发电机参数
	double A, B, N, KG;														//饱和系数

	//开始
	double rate = base / para_gene[k].Sh;

	//参数
	Xd = rate * para_gene[k].Xd;  Xd1 = rate * para_gene[k].Xd1;  Xd2 = rate * para_gene[k].Xd2;
	Xq = rate * para_gene[k].Xq;  Xq1 = rate * para_gene[k].Xq1;  Xq2 = rate * para_gene[k].Xq2;
	Ra = rate * para_gene[k].Ra;

	Td01 = para_gene[k].Td01 / rate;  Td02 = para_gene[k].Td02 / rate;
	Tq01 = para_gene[k].Tq01 / rate;  Tq02 = para_gene[k].Tq02 / rate;

	Tj = para_gene[k].Tj / rate;  D = para_gene[k].D / rate;

	A = para_gene[k].A;		B = para_gene[k].B;		N = para_gene[k].N;


	//初值
	Delta0 = init_gene[k].Delta0;
	Vd0 = init_gene[k].Vd0;      Vq0 = init_gene[k].Vq0;
	Id0 = init_gene[k].Id0;	     Iq0 = init_gene[k].Iq0;	  Efq0 = init_gene[k].Efq0;
	Eq10 = init_gene[k].Eq10;    Eq20 = init_gene[k].Eq20;
	Ed10 = init_gene[k].Ed10;    Ed20 = init_gene[k].Ed20;

	if (B == 0)	KG = 1.0;
	else
		KG = 1 + B * pow(Eq10, N - 1) / A;
	
	//形成模型
	Tgb[1][1] = 1.0;
	Tgb[2][2] = Tj;
	Tgb[3][3] = Td01;

	Jgb[1][2] = 100.0 * PI;
	Jgb[2][2] = -D;
	Jgb[2][4] = -1.0;
	Jgb[3][3] = -KG;
	Jgb[3][5] = -(Xd - Xd1);

	Jgb[4][4] = -1.0;
	Jgb[4][5] = Vd0;
	Jgb[4][6] = Vq0;
	Jgb[4][7] = Id0;
	Jgb[4][8] = Iq0;
	Jgb[5][5] = -Ra;	//- Ra△Id + Xq△Iq - △Vd = 0;
	Jgb[5][6] = Xq;
	Jgb[5][7] = -1.0;
	Jgb[6][3] = 1.0;	//△Eq1 - Ra△Iq - Xd1△Id - △Vq = 0;
	Jgb[6][5] = -Xd1;
	Jgb[6][6] = -Ra;
	Jgb[6][8] = -1.0;
	Jgb[7][1] = Vq0;	//sin(delta0)△Vx - cos(delta0)△Vy + Vq0△omega - △Vd = 0
	Jgb[7][7] = -1.0;
	Jgb[8][1] = -Vd0;	//cos(delta0)△Vx + sin(delta0)△Vy - Vd0△omega - △Vq = 0
	Jgb[8][8] = -1.0;
	Jgb[9][1] = Iq0;	//sin(delta0)△Ix - cos(delta0)△Iy + Iq0△omega - △Id = 0
	Jgb[9][5] = -1.0;
	Jgb[9][9] = sin(Delta0);
	Jgb[9][10] = -cos(Delta0);
	Jgb[10][1] = -Id0;	//cos(delta0)△Ix + sin(delta0)△Iy - Id0△omega - △Iq = 0
	Jgb[10][6] = -1.0;
	Jgb[10][9] = cos(Delta0);
	Jgb[10][10] = sin(Delta0);
}

//1型，Eq'恒定的2阶模型
void gene_type1(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,
	double **Tgb, double **Jgb, int k)
{
	double Delta0, Vd0, Vq0, Id0, Iq0, Efq0, Eq10, Eq20, Ed10, Ed20;		//发电机初值
	double Xd, Xd1, Xd2, Xq, Xq1, Xq2, Ra, Td01, Td02, Tq01, Tq02, Tj, D;	//发电机参数

	//开始
	double rate = base / para_gene[k].Sh;

	//参数
	Xd = rate * para_gene[k].Xd;  Xd1 = rate * para_gene[k].Xd1;  Xd2 = rate * para_gene[k].Xd2;
	Xq = rate * para_gene[k].Xq;  Xq1 = rate * para_gene[k].Xq1;  Xq2 = rate * para_gene[k].Xq2;
	Ra = rate * para_gene[k].Ra;

	Td01 = para_gene[k].Td01 / rate;  Td02 = para_gene[k].Td02 / rate;
	Tq01 = para_gene[k].Tq01 / rate;  Tq02 = para_gene[k].Tq02 / rate;

	Tj = para_gene[k].Tj / rate;  D = para_gene[k].D / rate;


	//初值
	Delta0 = init_gene[k].Delta0;
	Vd0 = init_gene[k].Vd0;      Vq0 = init_gene[k].Vq0;
	Id0 = init_gene[k].Id0;	     Iq0 = init_gene[k].Iq0;	  Efq0 = init_gene[k].Efq0;
	Eq10 = init_gene[k].Eq10;    Eq20 = init_gene[k].Eq20;
	Ed10 = init_gene[k].Ed10;    Ed20 = init_gene[k].Ed20;

	//形成模型
	Tgb[1][1] = 1.0;
	Tgb[2][2] = Tj;

	Jgb[1][2] = 100.0 * PI;
	Jgb[2][2] = -D;
	Jgb[2][3] = -1.0;

	Jgb[3][3] = -1.0;
	Jgb[3][4] = Vd0;
	Jgb[3][5] = Vq0;
	Jgb[3][6] = Id0;
	Jgb[3][7] = Iq0;
	Jgb[4][4] = -Ra;	//- Ra△Id + Xq△Iq - △Vd = 0;
	Jgb[4][5] = Xq;
	Jgb[4][6] = -1.0;
	Jgb[5][4] = -Xd1;	//- Ra△Iq - Xd1△Id - △Vq = 0;
	Jgb[5][5] = -Ra;
	Jgb[5][7] = -1.0;
	Jgb[6][1] = Vq0;	//sin(delta0)△Vx - cos(delta0)△Vy + Vq0△omega - △Vd = 0
	Jgb[6][6] = -1.0;
	Jgb[7][1] = -Vd0;	//cos(delta0)△Vx + sin(delta0)△Vy - Vd0△omega - △Vq = 0
	Jgb[7][7] = -1.0;
	Jgb[8][1] = Iq0;	//sin(delta0)△Ix - cos(delta0)△Iy + Iq0△omega - △Id = 0
	Jgb[8][4] = -1.0;
	Jgb[8][8] = sin(Delta0);
	Jgb[8][9] = -cos(Delta0);
	Jgb[9][1] = -Id0;	//cos(delta0)△Ix + sin(delta0)△Iy - Id0△omega - △Iq = 0
	Jgb[9][5] = -1.0;
	Jgb[9][8] = cos(Delta0);
	Jgb[9][9] = sin(Delta0);
}