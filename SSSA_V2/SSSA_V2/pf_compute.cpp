//�������㣬�ɵ��ɾ���͵�ѹ������㹦��
#include"head.h"

void pf_compute(vector<NODE> node_table, vector<NODE_DATA> &data_node, double **Ymatrix)
{
	int num_node = node_table.size();
	int nd = 2 * num_node;
	vector<NODE>::iterator kn;

	double **Umatrix = dmatrix(nd + 1, 1 + 1);
	double **Imatrix = dmatrix(nd + 1, 1 + 1);

	for (int k = 0; k < num_node; k++)		//ÿ���γ�һ���ڵ��Vx��Vy
	{
		double Vx, Vy;
		int no, i1, i2;

		Vx = data_node[k].vol_pu * cos(data_node[k].vol_ang / SITA);
		Vy = data_node[k].vol_pu * sin(data_node[k].vol_ang / SITA);

		//�ҵ��ɾ����Ӧ��λ��
		for (kn = node_table.begin(); kn != node_table.end(); kn++)
		{
			if (kn->name == data_node[k].name)
			{
				no = kn->no;
				break;
			}
		}
		i2 = no + no;
		i1 = i2 - 1;

		Umatrix[i1][1] = Vx;
		Umatrix[i2][1] = Vy;
	}

	matrix_mul(nd, nd, 1, Ymatrix, Umatrix, Imatrix);	//����õ�Ix, Iy

	for (int k = 0; k < num_node; k++)		//ÿ�μ���һ���ڵ��P,Q
	{
		double Vx, Vy, Ix, Iy;
		int no, i1, i2;
		double P, Q;

		Vx = data_node[k].vol_pu * cos(data_node[k].vol_ang / SITA);
		Vy = data_node[k].vol_pu * sin(data_node[k].vol_ang / SITA);

		//�Ҷ�Ӧ��Ix, Iy
		for (kn = node_table.begin(); kn != node_table.end(); kn++)
		{
			if (kn->name == data_node[k].name)
			{
				no = kn->no;
				break;
			}
		}
		i2 = no + no;
		i1 = i2 - 1;

		Ix = Imatrix[i1][1];
		Iy = Imatrix[i2][1];

		//���㹦��
		P = Vx * Ix + Vy * Iy;
		Q = Vy * Ix - Vx * Iy;

		data_node[k].P = P * 100.0;
		data_node[k].Q = Q * 100.0;
	}

	/*
	//�����ʾ����
	for (int i = 1; i <= 18; i++)
	{
		cout << Imatrix[i][1] << "\t";
	}
	cout << endl;
	*/
	
	
	
}