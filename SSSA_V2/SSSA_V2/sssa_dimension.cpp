#include"head.h"

void sssa_dimension(vector<int> &gdim, vector<string> &gname, vector<int> &gl_gene,
	vector<GENE_PARA> para_gene)
{
	// 寻找各发电机的变量数
	int dim_sum = 0;  //各发电机变量数之和

	gdim.push_back(0);		//保存每台发电机的变量数,可用来索引起始行数	

	int num_gene = para_gene.size();
	for (int k = 0; k < num_gene; k++)	//每次进行一台发电机的维数计算
	{
		string card_type_gene;
		string gene_model;
		int count = 0;			//维数计数

		card_type_gene = para_gene[k].card_type;
		gene_model = para_gene[k].gene_model;
		//exci_model, pss_model..

		if (card_type_gene == "MF")
		{
			if(gene_model == "four_winding")
			{gl_gene.push_back(13); gname.push_back(para_gene[k].name); count += 13; }		//四绕组模型
			if(gene_model == "three_winding")
			{gl_gene.push_back(12); gname.push_back(para_gene[k].name); count += 12;  }		//三绕组模型
			if (gene_model == "two_axle")
			{gl_gene.push_back(11); gname.push_back(para_gene[k].name); count += 11; }		//双轴模型
			if (gene_model == "Eq_change")
			{gl_gene.push_back(10); gname.push_back(para_gene[k].name); count += 10; }		//Eq'变化
			if (gene_model == "Eq_constant")
			{gl_gene.push_back(9); gname.push_back(para_gene[k].name); count += 9; }		//Eq'恒定
		}
		if (card_type_gene == "MC")
		{gl_gene.push_back(9); gname.push_back(para_gene[k].name); count += 9;}				//经典模型
		//if (card_type_gene == "MG")...
		

		//switch(exci_model)....switch(pss_model)...	
		//要注意这时 没有exci，pss..的情况

		dim_sum += count;	//每台发电机的总维数（包括同步机，励磁，PSS,调速器，且包括状态和代数变量）
		gdim.push_back(dim_sum);
	}
}