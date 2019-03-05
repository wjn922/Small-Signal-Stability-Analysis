//网络方程
#include"head.h"

void sssa_networksystem(double **Ysystem, vector<NODE> node_table, vector<NODE_DATA> data_node, vector<LOAD_PARA> para_load)
{
	double P, Q;
	double Vx, Vy, Ix, Iy;
	double V, V2;
	double ap, bp, cp;
	double aq, bq, cq;
	double pu, qu, gx, gy, bx, by;
	int no, i1, i2;

	int num_load = para_load.size();
	vector<NODE>::iterator kn;
	vector<NODE_DATA>::iterator knb;

	
	for (int k = 0; k < num_load; k++)		//每次修正一个负荷节点
	{
		for (knb = data_node.begin(); knb != data_node.end(); knb++)
		{
			if (knb->name == para_load[k].name)
			{
				P = -knb->P / base;
				Q = -knb->Q / base;
				V = knb->vol_pu;
				Vx = knb->vol_pu * cos(knb->vol_ang / SITA);
				Vy = knb->vol_pu * sin(knb->vol_ang / SITA);
				break;
			}
		}
		V2 = Vx * Vx + Vy * Vy;
		ap = para_load[k].ap;
		bp = para_load[k].bp;
		cp = para_load[k].cp;
		aq = para_load[k].aq;
		bq = para_load[k].bq;
		cq = para_load[k].cq;

		Ix = (Vx * P + Vy * Q) / V2;
		Iy = (Vy * P - Vx * Q) / V2;
		pu = (2.0 * ap + bp)* P / V;
		qu = (2.0 * aq + bq)* Q / V;
		gx = ((pu * Vx + qu * Vy)*Vx / V - Vx * Ix + Vy * Iy) / V2;
		bx = ((pu * Vx + qu * Vy)*Vy / V - Vx * Iy - Vy * Ix) / V2;
		by = ((pu * Vy - qu * Vx)*Vx / V - Vx * Iy - Vy * Ix) / V2;
		gy = ((pu * Vy - qu * Vx)*Vy / V + Vx * Ix - Vy * Iy) / V2;

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
		Ysystem[i1][i1] -= gx;
		Ysystem[i2][i2] -= gy;
		Ysystem[i1][i2] -= bx;
		Ysystem[i2][i1] -= by;
	}

}