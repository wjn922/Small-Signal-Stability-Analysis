//初值计算
#include"head.h"

void sssa_initvalue_gene(char name_st[], vector<NODE_DATA> data_node, vector<GENE_PARA> para_gene, vector<GENE_INIT> &init_gene)
{
	double P_pu, Q_pu;
	double Ra, Xq;	
	double Vx, Vy, Ix, Iy;
	double delta, sd, cd;			//转子角度，三角函数值
	double Vd0, Vq0, Id0, Iq0;		//d、q轴电压电流初值

	int num_gene = para_gene.size();
	vector<NODE_DATA>::iterator knb;
	vector<GENE_PARA>::iterator kgp;
	GENE_INIT gi;

	for (int k = 0; k < num_gene; k++)	//每次进行一台发电机的初值计算
	{
		for (knb = data_node.begin(); knb != data_node.end(); knb++)
		{
			if (knb->name == para_gene[k].name)
			{
				
				P_pu = knb->P / para_gene[k].Sh;
				Q_pu = knb->Q / para_gene[k].Sh;
				Vx = knb->vol_pu * cos(knb->vol_ang / SITA);
				Vy = knb->vol_pu * sin(knb->vol_ang / SITA);
				break;
			}
		}

		Ix = (P_pu * Vx + Q_pu * Vy) / (Vx * Vx + Vy * Vy);
		Iy = (P_pu * Vy - Q_pu * Vx) / (Vx * Vx + Vy * Vy);
		Ra = para_gene[k].Ra;	
		Xq = para_gene[k].Xq;

		//开始初值计算
		gi.Delta0 = delta = atan2(Vy + Ra * Iy + Xq * Ix, Vx + Ra * Ix - Xq * Iy);
		sd = sin(delta);
		cd = cos(delta);

		gi.Vd0 = Vd0 = sd * Vx - cd * Vy;		//坐标旋转公式
		gi.Vq0 = Vq0 = cd * Vx + sd * Vy;		//d、q轴电压电流初值
		gi.Id0 = Id0 = sd * Ix - cd * Iy;
		gi.Iq0 = Iq0 = cd * Ix + sd * Iy;

		gi.Efq0 = Vq0 + Ra * Iq0 + para_gene[k].Xd * Id0;	//空载电势初值
		gi.Eq10 = Vq0 + Ra * Iq0 + para_gene[k].Xd1 * Id0;	//暂态电势初值
		gi.Eq20 = Vq0 + Ra * Iq0 + para_gene[k].Xd2 * Id0;
		gi.Ed10 = Vd0 + Ra * Id0 - para_gene[k].Xq1 * Iq0;
		gi.Ed20 = Vd0 + Ra * Id0 - para_gene[k].Xq2 * Iq0;

		init_gene.push_back(gi);
	}

	 
	//存储初值计算结果
	char temp[50];
	strcpy_s(temp, name_st);
	strcat_s(temp, "发电机初值结果.init");
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
}