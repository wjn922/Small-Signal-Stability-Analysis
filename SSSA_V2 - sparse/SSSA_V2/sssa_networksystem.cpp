//�����緽�̵ĸ��ɽ�������
#include"head.h"

void sssa_networksystem(char name_st[],  vector<NODE> node_table, vector<NODE_DATA> data_node, 
	vector<LOAD_PARA> para_load, vector<GENE_PARA> para_gene, vector<int> gdim, vector<int> gl_gene,
	vector<int> &Trow, vector<int> &Tcol, vector<double> &Tvalue,
	vector<int> &Jrow, vector<int> &Jcol, vector<double> &Jvalue)
{
	//////////////////////////////���緽�̸��ɽڵ�����///////////////////////////////////
	//��Ԫ��洢���ɾ���
	vector<int> Yrow;
	vector<int> Ycol;
	vector<double> Yvalue;

	//��ȡ���ɾ���
	char temp[50];
	int ry, cy, row, col;
	double G, B;
	char dot;

	strcpy_s(temp, name_st);
	strcat_s(temp, "���ɾ���.txt");
	ifstream input;
	input.open(temp);
	if (!input)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}
	while (input.peek() != EOF)
	{
		input >> ry >> dot >> cy >> dot >> G >> dot >> B;
		//cout << ry << " " << cy << " " << G << " " << B << endl;
		row = (ry - 1) * 2;
		col = (cy - 1) * 2;
		Yrow.push_back(row + 1);
		Ycol.push_back(col + 1);
		Yvalue.push_back(G);
		Yrow.push_back(row + 1);
		Ycol.push_back(col + 2);
		Yvalue.push_back(-B);
		Yrow.push_back(row + 2);
		Ycol.push_back(col + 1);
		Yvalue.push_back(B);
		Yrow.push_back(row + 2);
		Ycol.push_back(col + 2);
		Yvalue.push_back(G);

	}

	//�������ɽڵ�
	int num_load = para_load.size();
	vector<NODE>::iterator kn;
	vector<NODE_DATA>::iterator knb;

	double P, Q;
	double Vx, Vy, Ix, Iy;
	double V, V2;
	double ap, bp, cp;
	double aq, bq, cq;
	double pu, qu, gx, gy, bx, by;
	int no, i1, i2;

	for (int k = 0; k < num_load; k++)		//ÿ������һ�����ɽڵ�
	{

		for (knb = data_node.begin(); knb != data_node.end(); knb++)
		{
			if (knb->name == para_load[k].name)
			{
				P = -knb->P / base;
				Q = -knb->Q / base;

				/*
				if (para_load[k].name == "ĸ��A")
					Q = 0.5;
				if (para_load[k].name == "ĸ��B")
					Q = 0.3;
				if (para_load[k].name == "ĸ��C")
					Q = 0.35;
				if (para_load[k].name == "ĸ��2")
					Q = 0.1;
				*/
				
				
				V = knb->vol_pu;
				Vx = knb->vol_pu * cos(knb->vol_ang / SITA);
				Vy = knb->vol_pu * sin(knb->vol_ang / SITA);
				break;
			}
		}
		//cout << para_load[k].name << '\t' << setw(20)<< P << setw(20) << Q << setw(20) << Vx << setw(20) << Vy  << endl;

		V2 = V * V;
		ap = para_load[k].ap;
		bp = para_load[k].bp;
		cp = para_load[k].cp;
		aq = para_load[k].aq;
		bq = para_load[k].bq;
		cq = para_load[k].cq;

		/*
		Ix = (Vx * P + Vy * Q) / V2;
		Iy = (Vy * P - Vx * Q) / V2;
		pu = (2.0 * ap + bp)* P / V;
		qu = (2.0 * aq + bq)* Q / V;
		gx = ((pu * Vx + qu * Vy) * Vx / V - Vx * Ix + Vy * Iy) / V2;
		bx = ((pu * Vx + qu * Vy) * Vy / V - Vx * Iy - Vy * Ix) / V2;
		by = ((pu * Vy - qu * Vx) * Vx / V - Vx * Iy - Vy * Ix) / V2;
		gy = ((pu * Vy - qu * Vx) * Vy / V + Vx * Ix - Vy * Iy) / V2;
		*/
		
		/*
		gx = (P * Vx * Vx * (bp + 2 * cp) + Q * Vx * Vy * (bq + 2 * cq)) / (V2 * V2) - P / V2;
		bx = (Q * Vy * Vy * (bq + 2 * cq) + P * Vx * Vy * (bp + 2 * cp)) / (V2 * V2) - Q / V2;	bx = -bx;
		by = (Q * Vx * Vx * (bq + 2 * cq) - P * Vx * Vy * (bp + 2 * cp)) / (V2 * V2) - Q / V2;
		gy = (P * Vy * Vy * (bp + 2 * cp) + Q * Vx * Vy * (bq + 2 * cq)) / (V2 * V2) - P / V2;
		*/

		//���迹ģ��
		gx = P / V2;
		bx = Q / V2;
		by = -Q / V2;
		gy = P / V2;
		

		for (kn = node_table.begin(); kn != node_table.end(); kn++)
		{
			if (kn->name == para_load[k].name)
			{
				no = kn->no;
				break;
			}
		}
		
		i2 = no + no;
		i1 = i2 - 1;

		//����
		for (int k = 0; k < Yrow.size(); k++)
		{
			if (Yrow[k] == i1 && Ycol[k] == i1)
				Yvalue[k] += gx;
			if (Yrow[k] == i1 && Ycol[k] == i2)
				Yvalue[k] += bx;
			if (Yrow[k] == i2 && Ycol[k] == i1)
				Yvalue[k] += by;
			if (Yrow[k] == i2 && Ycol[k] == i2)
				Yvalue[k] += gy;
		}
		
	}

	///////////////////////////////////////��ѹ���洢���緽�̾���////////////////////////////////
	int dim_gene_all = gdim[gdim.size() - 1];			//�����ϵͳȫ��ά��
	int dim_sys_all = dim_gene_all;						//����ֱ��ϵͳ�ȣ�ά�������ϼ�+...

	int num_node = node_table.size();					//�ڵ���
	int num_gene = para_gene.size();					//������ڵ���

	int count = Jrow[Jrow.size() - 1];				//���緽����ʼ���׵�ַ
	//cout << count << endl;


	for (int k = 0; k < num_node; k++)		//һ��һ�д洢���緽�̣����緽��ά��Ϊ2*num_node����һ��ѭ���γ�����
	{
		int flag = 0;		//��Ǹýڵ��Ƿ�Ϊ������ڵ�
		int index;
		for (int kk = 0; kk < num_gene; kk++)
		{
			if (para_gene[kk].name == node_table[k].name)	//�ҵ���Ӧ�ڵ�ķ��������
			{
				flag = 1;
				index = kk;
				break;
			}
		}

		if (flag)	//�Ƿ�����ڵ㣬��ע�����
		{
			//////////Ix��///////////////
			int start = gdim[index];
			col = start + gl_gene[k] - 1;		//��ӦIx
			Jcol.push_back(col);
			Jvalue.push_back(-1.0);
			count++;

			//�������ɾ���ȫ��Ԫ�أ����ڸ��е�����
			for (int i = 0; i < Yrow.size(); i++)
			{
				if (Yrow[i] == 2 * k + 1)
				{
					Jcol.push_back(dim_sys_all + Ycol[i]);
					Jvalue.push_back(Yvalue[i]);
					count++;
				}
			}
			Jrow.push_back(count);

			//////////Iy��/////////////////
			col = start + gl_gene[k];			///��ӦIy
			Jcol.push_back(col);
			Jvalue.push_back(-1.0);
			count++;

			//�������ɾ���ȫ��Ԫ�أ����ڸ��е�����
			for (int i = 0; i < Yrow.size(); i++)
			{
				if (Yrow[i] == 2 * k + 2)
				{
					Jcol.push_back(dim_sys_all + Ycol[i]);
					Jvalue.push_back(Yvalue[i]);
					count++;
				}
			}
			Jrow.push_back(count);

		}

		else          //����ڵ㣬��ע�����
		{
			//�������ɾ���ȫ��Ԫ�أ����ڸ��е�����
			for (int i = 0; i < Yrow.size(); i++)
			{
				if (Yrow[i] == 2 * k + 1)
				{
					Jcol.push_back(dim_sys_all + Ycol[i]);
					Jvalue.push_back(Yvalue[i]);
					count++;
				}
			}
			Jrow.push_back(count);

			//�������ɾ���ȫ��Ԫ�أ����ڸ��е�����
			for (int i = 0; i < Yrow.size(); i++)
			{
				if (Yrow[i] == 2 * k + 2)
				{
					Jcol.push_back(dim_sys_all + Ycol[i]);
					Jvalue.push_back(Yvalue[i]);
					count++;
				}
			}
			Jrow.push_back(count);
		}
	}

}