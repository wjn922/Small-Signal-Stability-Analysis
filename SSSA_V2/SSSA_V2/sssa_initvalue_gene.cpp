//��ֵ���� ͬ����
#include"head.h"

void sssa_initvalue_gene(char name_st[], vector<NODE_DATA> data_node, vector<GENE_PARA> para_gene, vector<GENE_INIT> &init_gene)
{
	double P_pu, Q_pu;
	double Ra, Xq;	
	double V, Vx, Vy, Ix, Iy;
	double delta, sd, cd;			//ת�ӽǶȣ����Ǻ���ֵ
	double Vd0, Vq0, Id0, Iq0;		//d��q���ѹ������ֵ

	int num_gene = para_gene.size();
	vector<NODE_DATA>::iterator knd;
	vector<GENE_PARA>::iterator kgp;
	GENE_INIT gi;

	for (int k = 0; k < num_gene; k++)	//ÿ�ν���һ̨������ĳ�ֵ����
	{
		double rate;
		for (knd = data_node.begin(); knd != data_node.end(); knd++)
		{
			if (knd->name == para_gene[k].name)
			{
				rate = base / para_gene[k].Sh;
				P_pu = knd->P / base;
				Q_pu = knd->Q / base;
				V = knd->vol_pu;
				Vx = knd->vol_pu * cos(knd->vol_ang / SITA);
				Vy = knd->vol_pu * sin(knd->vol_ang / SITA);
				break;
			}
		}


		Ra = rate * para_gene[k].Ra;
		if (para_gene[k].gene_model == "classic")
			Xq = rate * para_gene[k].Xd1;
		else
			Xq = rate * para_gene[k].Xq;

		
		Ix = (P_pu * Vx + Q_pu * Vy) / (V * V);
		Iy = (P_pu * Vy - Q_pu * Vx) / (V * V);

		//cout << para_gene[k].name << '\t' << setw(20) << P_pu << setw(20) << Q_pu << setw(20) <<
			//Vx << setw(20) << Vy << setw(20) << Ix << setw(20) << Iy << endl;


		//��ʼ��ֵ����
		gi.Delta0 = delta = atan2(Vy + Ra * Iy + Xq * Ix, Vx + Ra * Ix - Xq * Iy);
		sd = sin(delta);
		cd = cos(delta);

		gi.Vd0 = Vd0 = sd * Vx - cd * Vy;		//������ת��ʽ
		gi.Vq0 = Vq0 = cd * Vx + sd * Vy;		//d��q���ѹ������ֵ
		gi.Id0 = Id0 = sd * Ix - cd * Iy;
		gi.Iq0 = Iq0 = cd * Ix + sd * Iy;

		gi.Efq0 = Vq0 + Ra * Iq0 + rate * para_gene[k].Xd * Id0;		//���ص��Ƴ�ֵ
		gi.Eq10 = Vq0 + Ra * Iq0 + rate * para_gene[k].Xd1 * Id0;		//��̬���Ƴ�ֵ
		gi.Eq20 = Vq0 + Ra * Iq0 + rate * para_gene[k].Xd2 * Id0;
		gi.Ed10 = Vd0 + Ra * Id0 - rate * para_gene[k].Xq1 * Iq0;
		gi.Ed20 = Vd0 + Ra * Id0 - rate * para_gene[k].Xq2 * Iq0;

		init_gene.push_back(gi);
	}

	/*
	//�洢��ֵ������
	char temp[50];
	strcpy_s(temp, name_st);
	strcat_s(temp, "ͬ������ֵ���.init");
	ofstream out;
	out.open(temp);

	vector<GENE_INIT>::iterator kgi;
	for (kgi = init_gene.begin(); kgi != init_gene.end(); kgi++) {
		out << "Delta0 = " << kgi->Delta0 << endl;
		out << "Vd0 = " << kgi->Vd0 << endl;
		out << "Vq0 = " << kgi->Vq0 << endl;
		out << "Id0 = " << kgi->Id0 << endl;
		out << "Iq0 = " << kgi->Iq0 << endl;
		out << "Efq0 = " << kgi->Efq0 << endl;
		out << "Eq10 = " << kgi->Eq10 << endl;
		out << "Eq20 = " << kgi->Eq20 << endl;
		out << "Ed10 = " << kgi->Ed10 << endl;
		out << "Ed20 = " << kgi->Ed20 << endl << endl;
	}
	*/
	
}