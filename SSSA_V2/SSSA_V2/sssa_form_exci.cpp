#include"head.h"
#include"sssa_model_exci.cpp"

void sssa_form_exci(vector<EXCI_INIT> init_exci, vector<EXCI_PARA> para_exci,
	string exci_model, double **Tgb_exci, double **Jgb_exci, int k)
{
	if (exci_model == "EA")
		exci_EA(init_exci, para_exci, Tgb_exci, Jgb_exci, k);
	//if(exci_model == "F")			//其他模型
	
}