#include"head.h"

void sssa_connect_system(char name_st[], double **Tgsys, double **Jgsys, double **Ysystem, 
	vector<NODE> node_table, vector<GENE_INIT> init_gene, vector<GENE_PARA> para_gene,  
	vector<int> gdim, vector<int> gl_gene, double **T, double **J)
{
	//�����з�����ľ��������T��J��Ӧλ��
	int dim_gene_all = gdim[gdim.size() - 1];
	for (int i = 1; i <= dim_gene_all; i++) 
	{
		for (int j = 1; j <= dim_gene_all; j++)
		{
			T[i][j] = Tgsys[i][j];
			J[i][j] = Jgsys[i][j];
		}
	}


	//�����緽�̷�����T��J��Ӧλ��
	int num_node = node_table.size();
	int nd = 2 * num_node;
	for (int i = 1; i <= nd; i++)
	{
		for (int j = 1; j <= nd; j++)
		{
			J[dim_gene_all + i][dim_gene_all + j] = Ysystem[i][j];
		}
	}


	//��ȫ�������ѹ������ת �� ��������緽��
	vector<NODE>::iterator kn;
	int num_gene = init_gene.size();
	double Delta0;
	for (int k = 0; k < num_gene; k++)
	{
		Delta0 = init_gene[k].Delta0;
		int start = gdim[k];		//��k̨���������ʼ��

		//�ҷ����Vd��Vq
		int i1 = gl_gene[k] - 3;	//��Ӧ��k̨�������Vd��������i1+1��ӦVq����
		int row = start + i1;

		int no;
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

		//������ת����
		//cout << row << " " << dim_gene_all + col << endl;
		J[row][dim_gene_all + col] = sin(Delta0);
		J[row][dim_gene_all + col + 1] = -cos(Delta0);
		J[row + 1][dim_gene_all + col] = cos(Delta0);
		J[row + 1][dim_gene_all + col + 1] = sin(Delta0);

		//�ҷ����Ix��Iy
		row = col;
		i1 = gl_gene[k] - 1;		//��Ӧ��k̨�������Ix��������i1+1��ӦIy����
		col = start + i1;

		//��������緽��
		//cout << dim_gene_all + row << " " << col << endl;
		J[dim_gene_all + row][col] = -1;
		J[dim_gene_all + row + 1][col + 1] = -1;
	}


	///////////////////////////���ȫ�����Ի�ģ��
	char temp[20];
	strcpy_s(temp, name_st);
	strcat_s(temp, "_T.txt");
	ofstream out1;
	out1.open(temp);

	out1.setf(ios::fixed);
	out1.precision(16);


	for (int i = 1; i <= dim_gene_all + nd; i++)
	{
		for (int j = 1; j <= dim_gene_all + nd; j++)
			out1 << T[i][j] <<'\t';
		out1 << endl;
	}
	out1.close();

	strcpy_s(temp, name_st);
	strcat_s(temp, "_J.txt");
	ofstream out2;
	out2.open(temp);

	out2.setf(ios::fixed);
	out2.precision(16);

	for (int i = 1; i <= dim_gene_all + nd; i++)
	{
		for (int j = 1; j <= dim_gene_all + nd; j++)
			out2 << J[i][j] << '\t';
		out2 << endl;
	}
	out2.close();
}