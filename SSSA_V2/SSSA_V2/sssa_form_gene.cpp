#include"head.h"
#include"sssa_model_gene.cpp"

void sssa_form_gene(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,
	string card_type, string gene_model, double **Tgb_gene, double **Jgb_gene, int k)
{
	if (card_type == "MF")
	{
		if(gene_model == "four_winding")
			gene_four_winding(init_gene, para_gene, Tgb_gene, Jgb_gene, k); 		//四绕组模型
		if (gene_model == "three_winding") 
			gene_three_winding(init_gene, para_gene, Tgb_gene, Jgb_gene, k); 		//三绕组模型
		if (gene_model == "two_axle")
			gene_two_axle(init_gene, para_gene, Tgb_gene, Jgb_gene, k);				//双轴模型
		if (gene_model == "Eq_change")
			gene_Eq_change(init_gene, para_gene, Tgb_gene, Jgb_gene, k);			//Eq'变化模型
		if (gene_model == "Eq_constant")
			gene_Eq_constant(init_gene, para_gene, Tgb_gene, Jgb_gene, k);			//Eq'恒定模型

	}
	if (card_type == "MC")
		gene_classic(init_gene, para_gene, Tgb_gene, Jgb_gene, k);					//经典模型
	if (card_type == "MG")
	{
		if (gene_model == "type6")
			gene_type6(init_gene, para_gene, Tgb_gene, Jgb_gene, k);				//6型，Eq',Eq'',Ed',Ed''变化的6阶模型
		if (gene_model == "type5")
			gene_type5(init_gene, para_gene, Tgb_gene, Jgb_gene, k);				//5型，Eq',Ed'变化的4阶模型
		if (gene_model == "type3")
			gene_type3(init_gene, para_gene, Tgb_gene, Jgb_gene, k);				//3型，Eq',Eq'',Ed''变化的5阶模型
		if (gene_model == "type2")
			gene_type2(init_gene, para_gene, Tgb_gene, Jgb_gene, k);				//2型，Eq'变化的3阶模型
		if (gene_model == "type1")
			gene_type1(init_gene, para_gene, Tgb_gene, Jgb_gene, k);				//1型，Eq'恒定的2阶模型
		//type0 经典模型依然使用MC卡
	}
}