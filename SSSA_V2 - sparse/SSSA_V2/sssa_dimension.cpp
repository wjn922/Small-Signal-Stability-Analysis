#include"head.h"

void sssa_dimension(vector<int> &gdim, vector<string> &gname, vector<int> &gl_gene, vector<int> &gl_exci, vector<int> &gl_pss,
	vector<GENE_PARA> para_gene, vector<EXCI_PARA> para_exci, vector<PSS_PARA> para_pss)
{
	// Ѱ�Ҹ�������ı�����
	int dim_sum = 0;		//�������������֮��

	gdim.push_back(0);		//����ÿ̨������ı�����,������������ʼ����	

	int num_gene = para_gene.size();
	for (int k = 0; k < num_gene; k++)	//ÿ�ν���һ̨�������ά������
	{
		string card_type_gene;
		string gene_model, exci_model, pss_model;
		int count = 0;			//ά������

		card_type_gene = para_gene[k].card_type;
		gene_model = para_gene[k].gene_model;
		exci_model = para_exci[k].exci_model;
		pss_model = para_pss[k].pss_model;

		//exci_model, pss_model..

		//ͬ����ά��
		if (card_type_gene == "MF")
		{
			if(gene_model == "four_winding")
			{gl_gene.push_back(13); gname.push_back(para_gene[k].name); count += 13; }		//������ģ��
			if(gene_model == "three_winding")
			{gl_gene.push_back(12); gname.push_back(para_gene[k].name); count += 12;  }		//������ģ��
			if (gene_model == "two_axle")
			{gl_gene.push_back(11); gname.push_back(para_gene[k].name); count += 11; }		//˫��ģ��
			if (gene_model == "Eq_change")
			{gl_gene.push_back(10); gname.push_back(para_gene[k].name); count += 10; }		//Eq'�仯
			if (gene_model == "Eq_constant")
			{gl_gene.push_back(9); gname.push_back(para_gene[k].name); count += 9; }		//Eq'�㶨
		}
		if (card_type_gene == "MC")
		{gl_gene.push_back(9); gname.push_back(para_gene[k].name); count += 9;}				//����ģ��
		if (card_type_gene == "MG")
		{
			if (gene_model == "type6")
			{gl_gene.push_back(13); gname.push_back(para_gene[k].name); count += 13;}		//6�ͣ�Eq',Eq'',Ed',Ed''�仯��6��ģ��
			if (gene_model == "type3")
			{gl_gene.push_back(12); gname.push_back(para_gene[k].name); count += 12;}		//3�ͣ�Eq',Eq'',Ed''�仯��5��ģ��
			if (gene_model == "type5")
			{gl_gene.push_back(11); gname.push_back(para_gene[k].name); count += 11;}		//5�ͣ�Eq',Ed'�仯��4��ģ��
			if (gene_model == "type2")
			{gl_gene.push_back(10); gname.push_back(para_gene[k].name); count += 10;}		//2�ͣ�Eq'�仯��3��ģ��
			if (gene_model == "type1")
			{gl_gene.push_back(9); gname.push_back(para_gene[k].name); count += 9;}			//1�ͣ�Eq'�㶨��2��ģ��
		}

		//����ϵͳά��
		if (exci_model == "EA")
		{gl_exci.push_back(5);	count += 5;}
		else       //û������
		{gl_exci.push_back(0);	count += 0;}

		//PSSϵͳά��
		if(pss_model == "SF" || pss_model == "SP" || pss_model == "SS" || pss_model == "SG")
		{gl_pss.push_back(5);	count += 5;}
		else      //û��PSS
		{gl_pss.push_back(0);	count += 0;}
	
		

		//switch(exci_model)....switch(pss_model)...	
		//Ҫע����ʱ û��exci��pss..�����

		dim_sum += count;	//ÿ̨���������ά��������ͬ���������ţ�PSS,���������Ұ���״̬�ʹ���������
		gdim.push_back(dim_sum);
	}
}