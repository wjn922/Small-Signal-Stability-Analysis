//������������Ի�ģ��
#include"head.h"
#include"sssa_initvalue_gene.cpp"
#include"sssa_form_gene.cpp"

void sssa_genesystem(char name_st[], vector<NODE_DATA> data_node, vector<GENE_PARA> para_gene, vector<GENE_INIT> &init_gene,
	vector<int> gdim, vector<int> gl_gene, double **Tgsys, double **Jgsys)
{
	//��ֵ����
	sssa_initvalue_gene(name_st, data_node, para_gene, init_gene);
	//sssa_initvalue_exci...
	
	//���Ի�ģ��
	int i, j;
	int num_gene = para_gene.size();		//�������Ŀ

	string card_type_gene;					//card_type_exci, card_type_pss...
	string gene_model;						//exci_model, pss_model,...
	int dimg, dimg_gene;					//ÿ̨�������ά��,ͬ����ά��  dim_ex, dim_pss, dim_go...
	int start;								//�ݴ�ÿ̨���������ʼ����

	for (int k = 0; k < num_gene; k++)			//ÿ���γ�һ̨����������Ի�ģ��
	{
		card_type_gene = para_gene[k].card_type;
		gene_model = para_gene[k].gene_model;
		//exci,pss,....

		dimg = gl_gene[k];		//...+gl_ex[k]+gl_pss[k]+gl_go[k]
		dimg_gene = gl_gene[k];

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

		//��ʼ������������ţ�PSS�� ����...���־���


		//ͬ�������Ի�ģ��
		sssa_form_gene(init_gene, para_gene, card_type_gene, gene_model, Tgb_gene, Jgb_gene, k);


		//�������Ի�ģ��
		//PSS���Ի�ģ��
		//���������Ի�ģ��...

		//�����ϵͳ�� ���ţ�PSS��������������....


		//��ͬ�����ľ������ �����������Ӧλ��
		for (i = 1; i <= dimg_gene; i++) {
			for (j = 1; j <= dimg_gene; j++) {
				Tgb[i][j] = Tgb_gene[i][j];
				Jgb[i][j] = Jgb_gene[i][j];
			}
		}


		//�����ŵľ������ �����������Ӧλ��
		//��PSS�ľ������ �����������Ӧλ��...
		//���������ľ������ �����������Ӧλ��...


		//���÷����������뷢���ϵͳ�������Ӧλ��
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