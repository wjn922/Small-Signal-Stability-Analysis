//网络方程
#include"head.h"

void sssa_networksystem(double **Ysystem, vector<NODE> node_table, vector<NODE_DATA> data_node, vector<LOAD_PARA> para_load)
{
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

	for (int k = 0; k < num_load; k++)		//每次修正一个负荷节点
	{

		for (knb = data_node.begin(); knb != data_node.end(); knb++)
		{
			if (knb->name == para_load[k].name)
			{
				P = -knb->P / base;
				Q = -knb->Q / base;

				/*
				if (para_load[k].name == "母线A")
					Q = 0.5;
				if (para_load[k].name == "母线B")
					Q = 0.3;
				if (para_load[k].name == "母线C")
					Q = 0.35;
				if (para_load[k].name == "母线2")
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

		//恒阻抗模型
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
		Ysystem[i1][i1] += gx;
		Ysystem[i2][i2] += gy;
		Ysystem[i1][i2] += bx;
		Ysystem[i2][i1] += by;
	}



}