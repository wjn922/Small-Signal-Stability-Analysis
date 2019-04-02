//发电机部分线性化模型
#include"head.h"
#include"sssa_initvalue_gene.cpp"
#include"sssa_form_gene.cpp"
#include"sssa_initvalue_exci.cpp"
#include"sssa_form_exci.cpp"
#include"sssa_form_pss.cpp"

void sssa_genesystem(char name_st[], vector<NODE_DATA> data_node, vector<GENE_PARA> para_gene, vector<GENE_INIT> &init_gene,
	vector<EXCI_PARA> para_exci, vector<EXCI_INIT> &init_exci, vector<PSS_PARA> para_pss,
	vector<int> gdim, vector<int> gl_gene, vector<int> gl_exci, vector<int> gl_pss, double **Tgsys, double **Jgsys)
{
	//初值计算
	sssa_initvalue_gene(name_st, data_node, para_gene, init_gene);
	sssa_initvalue_exci(name_st, data_node, para_exci, init_gene, init_exci);
	
	//线性化模型
	int i, j;
	int num_gene = para_gene.size();					//发电机数目

	string card_type_gene;											//card_type_exci, card_type_pss...
	string gene_model, exci_model, pss_model;						//exci_model, pss_model,...
	int dimg, dimg_gene, dimg_exci, dimg_pss;						//每台发电机总维数,同步机维数  dim_ex, dim_pss, dim_go...
	int start;														//暂存每台发电机的起始行数

	for (int k = 0; k < num_gene; k++)					//每次形成一台发电机的线性化模型
	{
		card_type_gene = para_gene[k].card_type;
		gene_model = para_gene[k].gene_model;
		exci_model = para_exci[k].exci_model;
		pss_model = para_pss[k].pss_model;
		//exci,pss,....

		dimg = gl_gene[k] + gl_exci[k] + gl_pss[k];;		//...+gl_ex[k]+gl_pss[k]+gl_go[k]
		dimg_gene = gl_gene[k];
		dimg_exci = gl_exci[k];
		dimg_pss = gl_pss[k];

		//////////////////////////////////////////////////////////
		//初始化矩阵存储该发电机,包括同步机，励磁，PSS，调速器
		double **Jgb = dmatrix(dimg + 1, dimg + 1);
		double **Tgb = dmatrix(dimg + 1, dimg + 1);
		for (i = 1; i <= dimg; i++) {
			for (j = 1; j <= dimg; j++) {
				Tgb[i][j] = Jgb[i][j] = 0.0;
			}
		}

		//初始化发电机中的同步机部分矩阵
		double **Jgb_gene = dmatrix(dimg_gene + 1, dimg_gene + 1);
		double **Tgb_gene = dmatrix(dimg_gene + 1, dimg_gene + 1);
		for (i = 1; i <= dimg_gene; i++) {
			for (j = 1; j <= dimg_gene; j++) {
				Tgb_gene[i][j] = Jgb_gene[i][j] = 0.0;
			}
		}

		//初始化发电机中励磁部分矩阵
		double **Jgb_exci = dmatrix(dimg_exci + 1, dimg_exci + 1);
		double **Tgb_exci = dmatrix(dimg_exci + 1, dimg_exci + 1);
		for (i = 1; i <= dimg_exci; i++) {
			for (j = 1; j <= dimg_exci; j++) {
				Tgb_exci[i][j] = Jgb_exci[i][j] = 0.0;
			}
		}

		//初始化发电机中的PSS矩阵 
		double **Jgb_pss = dmatrix(dimg_pss + 1, dimg_pss + 1);
		double **Tgb_pss = dmatrix(dimg_pss + 1, dimg_pss + 1);
		for (i = 1; i <= dimg_pss; i++) {
			for (j = 1; j <= dimg_pss; j++) {
				Tgb_pss[i][j] = Jgb_pss[i][j] = 0.0;
			}
		}
		//初始化调速...部分矩阵


		////////////////////////////////////////////////////////////////
		//同步机线性化模型
		sssa_form_gene(init_gene, para_gene, card_type_gene, gene_model, Tgb_gene, Jgb_gene, k);
		//励磁线性化模型
		sssa_form_exci(init_exci, para_exci, exci_model, Tgb_exci, Jgb_exci, k);
		//PSS线性化模型
		sssa_form_pss(para_pss, pss_model, Tgb_pss, Jgb_pss, k);
		

		//PSS线性化模型
		//调速器线性化模型...

		

		////////////////////////////////////////////////////////////////////
		//将同步机的矩阵放入 发电机矩阵相应位置
		for (i = 1; i <= dimg_gene; i++) {
			for (j = 1; j <= dimg_gene; j++) {
				Tgb[i][j] = Tgb_gene[i][j];
				Jgb[i][j] = Jgb_gene[i][j];
			}
		}

		//将励磁的矩阵放入 发电机矩阵相应位置
		if (exci_model.size() != 0)
		{
			for (i = 1; i <= dimg_exci; i++){
				for (j = 1; j <= dimg_exci; j++){
					Tgb[dimg_gene + i][dimg_gene + j] = Tgb_exci[i][j];
					Jgb[dimg_gene + i][dimg_gene + j] = Jgb_exci[i][j];
				}
			}
		}

		//将PSS的矩阵放入 发电机矩阵相应位置
		if (pss_model.size() != 0)
		{
			for (i = 1; i <= dimg_pss; i++) {
				for (j = 1; j <= dimg_pss; j++) {
					Tgb[dimg_gene + dimg_exci + i][dimg_gene + dimg_exci + j] = Tgb_pss[i][j];
					Jgb[dimg_gene + dimg_exci + i][dimg_gene + dimg_exci + j] = Jgb_pss[i][j];
				}
			}
		}

		//将调速器的矩阵放入 发电机矩阵相应位置...



		//////////////////////////////////////////////////////////////////
		//发电机系统里 励磁，PSS，调速器的连接....
		if (exci_model.size() != 0)		//有励磁
		{
			//补同步机内元素
			if (gene_model == "four_winding" || gene_model == "three_winding" || gene_model == "two_axle" || gene_model == "Eq_change" 
				|| gene_model == "type6" || gene_model == "type5" || gene_model == "type3" || gene_model == "type2")
			{
				//Eq'方程里含有Efq
				Jgb[3][dimg_gene + 1] = 1.0;
			}
			//补励磁内元素
			if (exci_model == "EA")
			{
				double Vd0 = init_gene[k].Vd0;
				double Vq0 = init_gene[k].Vq0;
				double V0 = sqrt(Vd0 * Vd0 + Vq0 * Vq0);
				//Vm方程里有Vd，Vq
				Jgb[dimg_gene + dimg_exci][dimg_gene - 3] = Vd0 / V0;
				Jgb[dimg_gene + dimg_exci][dimg_gene - 2] = Vq0 / V0;
			}
		}
		if (pss_model.size() != 0)		//有PSS
		{
			//补励磁内元素
			if (exci_model == "EA")
			{
				double Ka = para_exci[k].Ka;
				//Vr方程里有Vs
				Jgb[dimg_gene + 2][dimg_gene + dimg_exci + dimg_pss] = Ka;
			}
			//补PSS内元素
			if (pss_model == "SF" || pss_model == "SP" || pss_model == "SS" || pss_model == "SG")
			{
				double Kqs = para_pss[k].Kqs;
				//V1方程里有Omega
				Jgb[dimg_gene + dimg_exci + 1][2] = Kqs;
			}
		}


		//////////////////////////////////////////////////////////////
		//将该发电机矩阵放入发电机系统矩阵的相应位置
		start = gdim[k];
		for (i = 1; i <= dimg; i++)
		{
			for (j = 1; j <= dimg; j++)
			{
				Tgsys[start + i][start + j] = Tgb[i][j];
				Jgsys[start + i][start + j] = Jgb[i][j];
			}
		}
	}


	/*
	/////////////////////////////////输出发电机模型
	int dim_gene_all = gdim[gdim.size() - 1];
	char temp[20];
	strcpy_s(temp, name_st);
	strcat_s(temp, "_geneT.model");
	ofstream out1;
	out1.open(temp);

	for (i = 1; i <= dim_gene_all; i++)
	{
		for (j = 1; j <= dim_gene_all; j++)
			out1 << Tgsys[i][j] << " ";
		out1 << endl;
	}
	out1.close();

	strcpy_s(temp, name_st);
	strcat_s(temp, "_geneJ.model");
	ofstream out2;
	out2.open(temp);

	for (i = 1; i <= dim_gene_all; i++)
	{
		for (j = 1; j <= dim_gene_all; j++)
			out2 << Jgsys[i][j] << " ";
		out2 << endl;
	}
	out2.close();
	*/
	
}