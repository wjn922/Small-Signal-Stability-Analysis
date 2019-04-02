#include"head.h"

void sssa_dimension(vector<int> &gdim, vector<string> &gname, vector<int> &gl_gene, vector<int> &gl_exci, vector<int> &gl_pss,
	vector<GENE_PARA> para_gene, vector<EXCI_PARA> para_exci, vector<PSS_PARA> para_pss)
{
	// 寻找各发电机的变量数
	int dim_sum = 0;		//各发电机变量数之和

	gdim.push_back(0);		//保存每台发电机的变量数,可用来索引起始行数	

	int num_gene = para_gene.size();
	for (int k = 0; k < num_gene; k++)	//每次进行一台发电机的维数计算
	{
		string card_type_gene;
		string gene_model, exci_model, pss_model;
		int count = 0;			//维数计数

		card_type_gene = para_gene[k].card_type;
		gene_model = para_gene[k].gene_model;
		exci_model = para_exci[k].exci_model;
		pss_model = para_pss[k].pss_model;

		//exci_model, pss_model..

		//同步机维数
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
		if (card_type_gene == "MG")
		{
			if (gene_model == "type6")
			{gl_gene.push_back(13); gname.push_back(para_gene[k].name); count += 13;}		//6型，Eq',Eq'',Ed',Ed''变化的6阶模型
			if (gene_model == "type3")
			{gl_gene.push_back(12); gname.push_back(para_gene[k].name); count += 12;}		//3型，Eq',Eq'',Ed''变化的5阶模型
			if (gene_model == "type5")
			{gl_gene.push_back(11); gname.push_back(para_gene[k].name); count += 11;}		//5型，Eq',Ed'变化的4阶模型
			if (gene_model == "type2")
			{gl_gene.push_back(10); gname.push_back(para_gene[k].name); count += 10;}		//2型，Eq'变化的3阶模型
			if (gene_model == "type1")
			{gl_gene.push_back(9); gname.push_back(para_gene[k].name); count += 9;}			//1型，Eq'恒定的2阶模型
		}

		//励磁系统维数
		if (exci_model == "EA")
		{gl_exci.push_back(5);	count += 5;}
		else       //没有励磁
		{gl_exci.push_back(0);	count += 0;}

		//PSS系统维数
		if(pss_model == "SF" || pss_model == "SP" || pss_model == "SS" || pss_model == "SG")
		{gl_pss.push_back(5);	count += 5;}
		else      //没有PSS
		{gl_pss.push_back(0);	count += 0;}
	
		

		//switch(exci_model)....switch(pss_model)...	
		//要注意这时 没有exci，pss..的情况

		dim_sum += count;	//每台发电机的总维数（包括同步机，励磁，PSS,调速器，且包括状态和代数变量）
		gdim.push_back(dim_sum);
	}
}