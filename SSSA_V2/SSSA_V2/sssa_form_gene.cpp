#include"head.h"
#include"sssa_model_gene.cpp"

void sssa_form_gene(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,
	string card_type, string gene_model, double **Tgb_gene, double **Jgb_gene, int k)
{
	if (card_type == "MF")
	{
		if(gene_model == "four_winding")
			gene_four_winding(init_gene, para_gene, Tgb_gene, Jgb_gene, k); 		//������ģ��
		if (gene_model == "three_winding") 
			gene_three_winding(init_gene, para_gene, Tgb_gene, Jgb_gene, k); 		//������ģ��
		if (gene_model == "two_axle")
			gene_two_axle(init_gene, para_gene, Tgb_gene, Jgb_gene, k);				//˫��ģ��
		if (gene_model == "Eq_change")
			gene_Eq_change(init_gene, para_gene, Tgb_gene, Jgb_gene, k);			//Eq'�仯ģ��
		if (gene_model == "Eq_constant")
			gene_Eq_constant(init_gene, para_gene, Tgb_gene, Jgb_gene, k);			//Eq'�㶨ģ��

	}
	if (card_type == "MC")
		gene_classic(init_gene, para_gene, Tgb_gene, Jgb_gene, k);					//����ģ��
	if (card_type == "MG")
	{
		if (gene_model == "type6")
			gene_type6(init_gene, para_gene, Tgb_gene, Jgb_gene, k);				//6�ͣ�Eq',Eq'',Ed',Ed''�仯��6��ģ��
		if (gene_model == "type5")
			gene_type5(init_gene, para_gene, Tgb_gene, Jgb_gene, k);				//5�ͣ�Eq',Ed'�仯��4��ģ��
		if (gene_model == "type3")
			gene_type3(init_gene, para_gene, Tgb_gene, Jgb_gene, k);				//3�ͣ�Eq',Eq'',Ed''�仯��5��ģ��
		if (gene_model == "type2")
			gene_type2(init_gene, para_gene, Tgb_gene, Jgb_gene, k);				//2�ͣ�Eq'�仯��3��ģ��
		if (gene_model == "type1")
			gene_type1(init_gene, para_gene, Tgb_gene, Jgb_gene, k);				//1�ͣ�Eq'�㶨��2��ģ��
		//type0 ����ģ����Ȼʹ��MC��
	}
}