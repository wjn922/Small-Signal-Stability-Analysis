//������������Ի�ģ��
#include"head.h"
#include"sssa_initvalue_gene.cpp"
#include"sssa_form_gene.cpp"
#include"sssa_initvalue_exci.cpp"
#include"sssa_form_exci.cpp"
#include"sssa_form_pss.cpp"

//ÿ̨�����������ͬ���������ţ�PSS�������������ö�ά����洢
void sssa_genesystem(char name_st[], vector<NODE> node_table, vector<NODE_DATA> data_node, vector<GENE_PARA> para_gene, vector<GENE_INIT> &init_gene,
	vector<EXCI_PARA> para_exci, vector<EXCI_INIT> &init_exci, vector<PSS_PARA> para_pss,
	vector<int> gdim, vector<int> gl_gene, vector<int> gl_exci, vector<int> gl_pss,
	vector<int> &Trow, vector<int> &Tcol, vector<double> &Tvalue,
	vector<int> &Jrow, vector<int> &Jcol, vector<double> &Jvalue)
{
	//��ֵ����
	sssa_initvalue_gene(name_st, data_node, para_gene, init_gene);
	sssa_initvalue_exci(name_st, data_node, para_exci, init_gene, init_exci);
	
	//���Ի�ģ��
	int i, j;
	int num_gene = para_gene.size();					//�������Ŀ

	string card_type_gene;								//card_type_exci, card_type_pss...
	string gene_model, exci_model, pss_model;			//exci_model, pss_model,...
	int dimg, dimg_gene, dimg_exci, dimg_pss;			//ÿ̨�������ά��,ͬ����ά��  dim_ex, dim_pss, dim_go...
	int start;											//�ݴ�ÿ̨���������ʼ����

	int dim_gene_all = gdim[gdim.size() - 1];			//�����ϵͳȫ��ά��
	int dim_sys_all = dim_gene_all;						//����ֱ��ϵͳ�ȣ�ά�������ϼ�

	int count1 = 1;										//T������ѹ���洢��ʼֵ
	int count2 = 1;										//J������ѹ���洢��ʼֵ
	Trow.push_back(count1);								//��1�е����׵�ַԪ��Ϊ1
	Jrow.push_back(count2);

	for (int k = 0; k < num_gene; k++)					//ÿ���γ�һ̨����������Ի�ģ��
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
		//��ʼ������洢�÷����,����ͬ���������ţ�PSS��������
		double **Jgb = dmatrix(dimg + 1, dimg + 1);
		double **Tgb = dmatrix(dimg + 1, dimg + 1);
		for (i = 1; i <= dimg; i++) {
			for (j = 1; j <= dimg; j++) {
				Tgb[i][j] = Jgb[i][j] = 0.0;
			}
		}

		//��ʼ��������е�ͬ�������־���
		double **Jgb_gene = dmatrix(dimg_gene + 1, dimg_gene + 1);
		double **Tgb_gene = dmatrix(dimg_gene + 1, dimg_gene + 1);
		for (i = 1; i <= dimg_gene; i++) {
			for (j = 1; j <= dimg_gene; j++) {
				Tgb_gene[i][j] = Jgb_gene[i][j] = 0.0;
			}
		}

		//��ʼ������������Ų��־���
		double **Jgb_exci = dmatrix(dimg_exci + 1, dimg_exci + 1);
		double **Tgb_exci = dmatrix(dimg_exci + 1, dimg_exci + 1);
		for (i = 1; i <= dimg_exci; i++) {
			for (j = 1; j <= dimg_exci; j++) {
				Tgb_exci[i][j] = Jgb_exci[i][j] = 0.0;
			}
		}

		//��ʼ��������е�PSS���� 
		double **Jgb_pss = dmatrix(dimg_pss + 1, dimg_pss + 1);
		double **Tgb_pss = dmatrix(dimg_pss + 1, dimg_pss + 1);
		for (i = 1; i <= dimg_pss; i++) {
			for (j = 1; j <= dimg_pss; j++) {
				Tgb_pss[i][j] = Jgb_pss[i][j] = 0.0;
			}
		}
		//��ʼ������...���־���


		////////////////////////////////////////////////////////////////
		//ͬ�������Ի�ģ��
		sssa_form_gene(init_gene, para_gene, card_type_gene, gene_model, Tgb_gene, Jgb_gene, k);
		//�������Ի�ģ��
		sssa_form_exci(init_exci, para_exci, exci_model, Tgb_exci, Jgb_exci, k);
		//PSS���Ի�ģ��
		sssa_form_pss(para_pss, pss_model, Tgb_pss, Jgb_pss, k);
		

		//PSS���Ի�ģ��
		//���������Ի�ģ��...

		

		////////////////////////////////////////////////////////////////////
		//��ͬ�����ľ������ �����������Ӧλ��
		for (i = 1; i <= dimg_gene; i++) {
			for (j = 1; j <= dimg_gene; j++) {
				Tgb[i][j] = Tgb_gene[i][j];
				Jgb[i][j] = Jgb_gene[i][j];
			}
		}

		//�����ŵľ������ �����������Ӧλ��
		if (exci_model.size() != 0)
		{
			for (i = 1; i <= dimg_exci; i++){
				for (j = 1; j <= dimg_exci; j++){
					Tgb[dimg_gene + i][dimg_gene + j] = Tgb_exci[i][j];
					Jgb[dimg_gene + i][dimg_gene + j] = Jgb_exci[i][j];
				}
			}
		}

		//��PSS�ľ������ �����������Ӧλ��
		if (pss_model.size() != 0)
		{
			for (i = 1; i <= dimg_pss; i++) {
				for (j = 1; j <= dimg_pss; j++) {
					Tgb[dimg_gene + dimg_exci + i][dimg_gene + dimg_exci + j] = Tgb_pss[i][j];
					Jgb[dimg_gene + dimg_exci + i][dimg_gene + dimg_exci + j] = Jgb_pss[i][j];
				}
			}
		}

		//���������ľ������ �����������Ӧλ��...



		//////////////////////////////////////////////////////////////////
		//�����ϵͳ�� ���ţ�PSS��������������....
		if (exci_model.size() != 0)		//������
		{
			//��ͬ������Ԫ��
			if (gene_model == "four_winding" || gene_model == "three_winding" || gene_model == "two_axle" || gene_model == "Eq_change"
				|| gene_model == "type6" || gene_model == "type5" || gene_model == "type3" || gene_model == "type2")
			{
				//Eq'�����ﺬ��Efq
				Jgb[3][dimg_gene + 1] = 1.0;
			}
			//��������Ԫ��
			if (exci_model == "EA")
			{
				double Vd0 = init_gene[k].Vd0;
				double Vq0 = init_gene[k].Vq0;
				double V0 = sqrt(Vd0 * Vd0 + Vq0 * Vq0);
				//Vm��������Vd��Vq
				Jgb[dimg_gene + dimg_exci][dimg_gene - 3] = Vd0 / V0;
				Jgb[dimg_gene + dimg_exci][dimg_gene - 2] = Vq0 / V0;
			}
		}
		if (pss_model.size() != 0)		//��PSS
		{
			//��������Ԫ��
			if (exci_model == "EA")
			{
				double Ka = para_exci[k].Ka;
				//Vr��������Vs
				Jgb[dimg_gene + 2][dimg_gene + dimg_exci + dimg_pss] = Ka;
			}
			//��PSS��Ԫ��
			if (pss_model == "SF" || pss_model == "SP" || pss_model == "SS" || pss_model == "SG")
			{
				double Kqs = para_pss[k].Kqs;
				//V1��������Omega
				Jgb[dimg_gene + dimg_exci + 1][2] = Kqs;
			}
		}


		//////////////////////////////////////////////////////////////
		//��ѹ����ʽ�洢�÷��������
		start = gdim[k];
		//cout << dimg << endl;
		for (i = 1; i <= dimg; i++)
		{
			for (j = 1; j <= dimg; j++)
			{
				if (Tgb[i][j] != 0.0)	//��Ԫ��
				{
					count1++;
					Tcol.push_back(start+j);
					Tvalue.push_back(Tgb[i][j]);
				}
				if (Jgb[i][j] != 0.0)	//��Ԫ��
				{
					count2++;
					Jcol.push_back(start + j);
					Jvalue.push_back(Jgb[i][j]);
				}
			}

		
			//��Ϊ������Ķ�ӦVd��Vq�У���Ҫ�����ѹ������ת�任����
			//Vd
			int i1 = gl_gene[k] - 3;	//��Ӧ��k̨�������Vd����
			if (i == i1)
			{	
				int no;
				vector<NODE>::iterator kn;
				for (kn = node_table.begin(); kn != node_table.end(); kn++)
				{
					if (kn->name == para_gene[k].name)
					{
						no = kn->no;
						break;
					}
				}
				int i2 = no + no;
				int col = i2 - 1;

				//������ת����;
				double Delta0 = init_gene[k].Delta0;
				Jcol.push_back(dim_sys_all + col);
				Jvalue.push_back(sin(Delta0));
				count2++;
				Jcol.push_back(dim_sys_all + col + 1);
				Jvalue.push_back(-cos(Delta0));
				count2++;
			}
			
			//Vq
			i1 = gl_gene[k] - 2;	//��Ӧ��k̨�������Vq
			if (i == i1)
			{
				int no;
				vector<NODE>::iterator kn;
				for (kn = node_table.begin(); kn != node_table.end(); kn++)
				{
					if (kn->name == para_gene[k].name)
					{
						no = kn->no;
						break;
					}
				}
				int i2 = no + no;
				int col = i2 - 1;

				//������ת����;
				double Delta0 = init_gene[k].Delta0;
				Jcol.push_back(dim_sys_all + col);
				Jvalue.push_back(cos(Delta0));
				count2++;
				Jcol.push_back(dim_sys_all + col + 1);
				Jvalue.push_back(sin(Delta0));
				count2++;
			}

			Trow.push_back(count1);		//���׵�ַ ��ά��Ϊ��������+1
			Jrow.push_back(count2);
		}
		

		delete Tgb;		delete Tgb_gene;	delete Tgb_exci;	delete Tgb_pss;
		delete Jgb;		delete Jgb_gene;	delete Jgb_exci;	delete Jgb_pss;
	}


	/*
	/////////////////////////////////��������ģ��
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