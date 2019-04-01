//对网络方程的负荷进行修正
#include"head.h"

void sssa_networksystem(char name_st[],  vector<NODE> node_table, vector<NODE_DATA> data_node, 
	vector<LOAD_PARA> para_load, vector<GENE_PARA> para_gene, vector<int> gdim, vector<int> gl_gene,
	vector<int> &Trow, vector<int> &Tcol, vector<double> &Tvalue,
	vector<int> &Jrow, vector<int> &Jcol, vector<double> &Jvalue)
{
	//////////////////////////////网络方程负荷节点修正///////////////////////////////////
	//三元组存储导纳矩阵
	vector<int> Yrow;
	vector<int> Ycol;
	vector<double> Yvalue;

	//读取导纳矩阵
	char temp[50];
	int ry, cy, row, col;
	double G, B;
	char dot;

	strcpy_s(temp, name_st);
	strcat_s(temp, "导纳矩阵.txt");
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

	//修正负荷节点
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

		//修正
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

	///////////////////////////////////////行压缩存储网络方程矩阵////////////////////////////////
	int dim_gene_all = gdim[gdim.size() - 1];			//发电机系统全部维数
	int dim_sys_all = dim_gene_all;						//加上直流系统等，维数再往上加+...

	int num_node = node_table.size();					//节点数
	int num_gene = para_gene.size();					//发电机节点数

	int count = Jrow[Jrow.size() - 1];				//网络方程起始行首地址
	//cout << count << endl;


	for (int k = 0; k < num_node; k++)		//一行一行存储网络方程（网络方程维数为2*num_node），一次循环形成两行
	{
		int flag = 0;		//标记该节点是否为发电机节点
		int index;
		for (int kk = 0; kk < num_gene; kk++)
		{
			if (para_gene[kk].name == node_table[k].name)	//找到对应节点的发电机索引
			{
				flag = 1;
				index = kk;
				break;
			}
		}

		if (flag)	//是发电机节点，有注入电流
		{
			//////////Ix行///////////////
			int start = gdim[index];
			col = start + gl_gene[k] - 1;		//对应Ix
			Jcol.push_back(col);
			Jvalue.push_back(-1.0);
			count++;

			//遍历导纳矩阵全部元素，属于该行的填入
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

			//////////Iy行/////////////////
			col = start + gl_gene[k];			///对应Iy
			Jcol.push_back(col);
			Jvalue.push_back(-1.0);
			count++;

			//遍历导纳矩阵全部元素，属于该行的填入
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

		else          //联络节点，无注入电流
		{
			//遍历导纳矩阵全部元素，属于该行的填入
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

			//遍历导纳矩阵全部元素，属于该行的填入
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