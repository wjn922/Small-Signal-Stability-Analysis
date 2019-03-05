#include"head.h"
#include"sssa_model_gene.cpp"

void sssa_form_gene(vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,
	string card_type, string gene_model, double **Tgb_gene, double **Jgb_gene, int k)
{
	if (card_type == "MF")
	{
		if(gene_model == "four_winding")
			gene_four_winding(init_gene, para_gene, Tgb_gene, Jgb_gene, k); 		//M������ģ��
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
}