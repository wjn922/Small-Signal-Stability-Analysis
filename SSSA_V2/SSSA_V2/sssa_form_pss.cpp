#include"head.h"
#include"sssa_model_pss.cpp"

void sssa_form_pss(vector<PSS_PARA> para_pss, string pss_model, double **Tgb_pss, double **Jgb_pss, int k)
{
	if (pss_model == "SS" )
		pss_SS(para_pss, Tgb_pss, Jgb_pss, k);
	//if(pss_model == "SI")			//其他模型

}