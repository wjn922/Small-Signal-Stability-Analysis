#include"head.h"

//读取导纳矩阵
void sssa_inputY(char name_st[], double **Ymatrix)
{
	char temp[50];
	int ry, cy, row, col;
	double G, B;
	char dot;

	//读导纳矩阵
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
		Ymatrix[row + 1][col + 1] = G;
		Ymatrix[row + 1][col + 2] = -B;
		Ymatrix[row + 2][col + 1] = B;
		Ymatrix[row + 2][col + 2] = G;
	}
}


//读取节点对应表
void sssa_inputnode(char name_st[], vector<NODE> &node_table)
{
	char temp[50];

	///////////////////////////读节点对应表
	strcpy_s(temp, name_st);
	strcat_s(temp, "节点对应表.txt");
	ifstream input;

	input.open(temp);
	if (!input)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}
	while (!input.eof())			//文件未读取完	
	{
		NODE node;
		char dot;
		input >> node.no >> dot >> node.name >> node.vol;
		node_table.push_back(node);
	}
	node_table.pop_back();		//若直接!input.eof()最后一行会输出两遍，弹出最后一个元素
	input.close();
}


//读取各节点潮流结果
void sssa_inputPF(char name_st[], vector<NODE> node_table, vector<NODE_DATA> &data_node)
{
	
	char temp[50];

	int num_node = node_table.size();
	int num_gene = 0;
	vector<NODE>::iterator kn;
	for (kn = node_table.begin(); kn != node_table.end(); kn++)
	{
		string name = kn->name;
		if (name.find("发电机") != name.npos)		//是发电机
			num_gene++;
	}

	strcpy_s(temp, name_st);
	strcat_s(temp, "潮流结果.txt");
	ifstream input;

	input.open(temp);
	if (!input)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}

	string line;
	//节点电压信息
	getline(input, line);
	getline(input, line);		//前两行
	for(int i = 0; i<num_node; i++)
	{
		getline(input, line);
		stringstream word(line);
		NODE_DATA nb;
		word >> nb.no >> nb.name >> nb.vol_base >> nb.vol_pu >> nb.vol >> nb.vol_ang;
		nb.P = nb.Q = 0.0;
		data_node.push_back(nb);
	}

	//发电机功率信息
	getline(input, line);
	getline(input, line);		//再两行
	for(int i = 0; i < num_gene; i++)
		getline(input, line);

	//线路/变压器潮流信息
	getline(input, line);
	getline(input, line);		//再两行
	while (getline(input, line))
	{
		int Fno, Tno;						//送端编号，受端编号
		string Fname, Tname;				//送端母线名，受端母线名
		double Fvol, Tvol;					//送端母线基准电压， 受端母线基准电压
		double FP, FQ, TP, TQ;				//送端有功，送端无功，受端有功，受端无功

		stringstream word(line);
		word >> Fno >> Fname >> Fvol >> Tno >> Tname >> Tvol >> FP >> FQ >> TP >> TQ;

		vector<NODE_DATA>::iterator knb;
		for (knb = data_node.begin(); knb != data_node.end(); knb++)
		{
			if (knb->name == Fname)	 //发电机母线
			{
				knb->P += FP;
				knb->Q += FQ;
			}
			if (knb->name == Tname)	 //发电机母线
			{
				knb->P += TP;
				knb->Q += TQ;
			}
		}
	}
	input.close();
	
	vector<NODE_DATA>::iterator knb;
	for (knb = data_node.begin(); knb != data_node.end(); knb++)
	{
		if (knb->P > 0)		knb->type = -1;			//发电机节点
		else if (knb->P < 0)	knb->type = 1;		//负荷节点
		else	knb->type = 0;						//联络节点
	}
	
	
	/*
	/////////////////////////读潮流结果所有节点母线电压
	strcpy_s(temp, name_st);
	strcat_s(temp, "潮流结果_Bus.txt");
	ifstream input1;
	
	input1.open(temp);
	if (!input1)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}
	while (!input1.eof())			//文件未读取完	
	{
		NODE_DATA nb;
		input1 >> nb.no >> nb.name >> nb.vol_base >> nb.vol_pu >> nb.vol >> nb.vol_ang;
		nb.P = nb.Q = 0.0;
		data_node.push_back(nb);
	}
	data_node.pop_back();		//若直接!input.eof()最后一行会输出两遍，弹出最后一个元素
	input1.close();


	//////////////////////////读潮流结果各节点功率
	strcpy_s(temp, name_st);
	strcat_s(temp, "潮流结果_Line.txt");
	ifstream input2;

	input2.open(temp);
	if (!input2)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}
	string line;
	while (getline(input2, line))			//按行读取文件
	{
		int Fno, Tno;						//送端编号，受端编号
		string Fname, Tname;				//送端母线名，受端母线名
		double Fvol, Tvol;					//送端母线基准电压， 受端母线基准电压
		double FP, FQ, TP, TQ;				//送端有功，送端无功，受端有功，受端无功

		stringstream word(line);
		word >> Fno >> Fname >> Fvol >> Tno >> Tname >> Tvol >> FP >> FQ >> TP >> TQ;

		vector<NODE_DATA>::iterator knb;
		for (knb = data_node.begin(); knb != data_node.end(); knb++)
		{
			if (knb->name == Fname)	 //发电机母线
			{
				knb->P += FP;
				knb->Q += FQ;
			}
			if (knb->name == Tname)	 //发电机母线
			{
				knb->P += TP;
				knb->Q += TQ;
			}
		}
	}
	input2.close;

	vector<NODE_DATA>::iterator knb;
	for (knb = data_node.begin(); knb != data_node.end(); knb++)
	{
		if (knb->P > 0)		knb->type = -1;			//发电机节点
		else if (knb->P < 0)	knb->type = 1;		//负荷节点
		else	knb->type = 0;						//联络节点
	}
	*/

}


//读取发电机参数（包括同步机，励磁，PSS，调速器）
void sssa_inputgene(char name_st[], vector<GENE_PARA> &para_gene)
{
	char temp[50];

	/////////////////////////读各发电机参数
	strcpy_s(temp, name_st);
	strcat_s(temp, ".swi");
	ifstream input1;

	input1.open(temp);
	if (!input1)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}

	char ch[100];
	string sh;
	char *zh;
	while (input1.getline(ch, 100))			//按行读取文件
	{
		sh = ch;
		//cout << sh << endl;
									
		string card_type = sh.substr(0, 3);
		card_type.erase(0, card_type.find_first_not_of(" "));	//清除头尾空格
		card_type.erase(card_type.find_last_not_of(" ") + 1);

		//同步机参数
		GENE_PARA gp;

		if (card_type == "M")					//发电机次暂态参数模型卡，排在发电机的第一个
		{
			string name = sh.substr(3, 8);						//电机母线名
			name.erase(0, name.find_first_not_of(" "));
			name.erase(name.find_last_not_of(" ") + 1);
			gp.name = name;
			string vol_str = sh.substr(11, 4);					//电机母线基准电压
			gp.vol = atof(vol_str.c_str());
			string capacity_str = sh.substr(16, 5);				//电机额定容量
			gp.capacity = atof(capacity_str.c_str());
			string power_factor_str = sh.substr(22, 3);			//电机功率因数
			gp.power_factor = atof(power_factor_str.c_str());
			string type = sh.substr(30, 2);						//电机类型
			type.erase(0, type.find_first_not_of(" "));
			type.erase(type.find_last_not_of(" ") + 1);
			gp.type = type;

			zh = &ch[37];										//Xd2
			string Xd2_str = zh;
			Xd2_str = Xd2_str.substr(0, 5);
			gp.Xd2 = atof(Xd2_str.c_str());
			zh = &ch[42];										//Xq2
			string Xq2_str = zh;
			Xq2_str = Xq2_str.substr(0, 5);
			gp.Xq2 = atof(Xq2_str.c_str());
			zh = &ch[47];										//Td02
			string Td02_str = zh;
			Td02_str = Td02_str.substr(0, 4);
			gp.Td02 = atof(Td02_str.c_str());
			zh = &ch[51];										//Tq02
			string Tq02_str = zh;
			Tq02_str = Tq02_str.substr(0, 4);
			gp.Tq02 = atof(Tq02_str.c_str());
			
		}

		if (card_type == "MF")					//发电机模型卡
		{
			gp.card_type = card_type;							//发电机卡片名
			string name = sh.substr(3, 8);						//电机母线名
			name.erase(0, name.find_first_not_of(" "));
			name.erase(name.find_last_not_of(" ") + 1);
			gp.name = name;

			zh = &ch[11];										//电机母线基准电压
			string vol_str = zh;
			vol_str = vol_str.substr(0, 4);
			gp.vol = atof(vol_str.c_str());
			zh = &ch[16];										//发电机动能
			string Ek_str = zh;
			Ek_str = Ek_str.substr(0, 6);
			gp.Ek = atof(Ek_str.c_str());
			zh = &ch[28];										//标幺参数的基准容量
			string Sh_str = zh;
			Sh_str = Sh_str.substr(0, 4);
			gp.Sh = atof(Sh_str.c_str());
			zh = &ch[32];										//Ra
			string Ra_str = zh;
			Ra_str = Ra_str.substr(0, 4);
			gp.Ra = atof(Ra_str.c_str());
			zh = &ch[36];										//Xd1
			string Xd1_str = zh;
			Xd1_str = Xd1_str.substr(0, 5);
			gp.Xd1 = atof(Xd1_str.c_str());
			zh = &ch[41];										//Xq1
			string Xq1_str = zh;
			Xq1_str = Xq1_str.substr(0, 5);
			gp.Xq1 = atof(Xq1_str.c_str());
			zh = &ch[46];										//Xd
			string Xd_str = zh;
			Xd_str = Xd_str.substr(0, 5);
			gp.Xd = atof(Xd_str.c_str());
			zh = &ch[51];										//Xq
			string Xq_str = zh;
			Xq_str = Xq_str.substr(0, 5);
			gp.Xq = atof(Xq_str.c_str());
			zh = &ch[56];										//Td01
			string Td01_str = zh;	
			Td01_str = Td01_str.substr(0, 4);
			gp.Td01 = atof(Td01_str.c_str());
			zh = &ch[60];										//Tq01	
			string Tq01_str = zh;
			Tq01_str = Tq01_str.substr(0, 3);
			gp.Tq01 = atof(Tq01_str.c_str());
			zh = &ch[63];										//XL
			string XL_str = zh;
			XL_str = XL_str.substr(0, 5);
			gp.XL = atof(XL_str.c_str());
			zh = &ch[77];										//D
			string D_str = zh;
			D_str = D_str.substr(0, 3);
			gp.D = atof(D_str.c_str());

			gp.Tj = 2 * gp.Ek / gp.Sh;
			//判断同步机模型
			if (gp.Td02 != 0.0 && gp.Tq02 != 0.0)		//含有次暂态参数
			{
				if (gp.Tq01 != 0 && gp.Xq != gp.Xq1)			//四绕组模型
					gp.gene_model = "four_winding";
				if (gp.Tq01 == 0 || gp.Xq == gp.Xq1)			//三绕组模型
					gp.gene_model = "three_winding";
			}
			else
			{
				if (gp.Tq01 != 0 && gp.Xq != gp.Xq1)	//不含次暂态参数
					gp.gene_model = "two_axle";					//双轴模型
				if (gp.Tq01 == 0 || gp.Xq == gp.Xq1)
					gp.gene_model = "Eq_change";				//Eq'变化模型
			}
			if (gp.Td01 == 999.0)
				gp.gene_model = "Eq_constant";					//Eq'恒定模型

			para_gene.push_back(gp);		//将同步机参数压栈
		}

		if (card_type == "MC")					//发电机模型卡(经典模型)
		{
			gp.card_type = card_type;							//发电机卡片名
			string name = sh.substr(3, 8);						//电机母线名
			name.erase(0, name.find_first_not_of(" "));
			name.erase(name.find_last_not_of(" ") + 1);
			gp.name = name;
			
			zh = &ch[11];										//电机母线基准电压
			string vol_str = zh;
			vol_str = vol_str.substr(0, 4);
			gp.vol = atof(vol_str.c_str());
			zh = &ch[16];										//发电机动能
			string Ek_str = zh;
			Ek_str = Ek_str.substr(0, 6);
			gp.Ek = atof(Ek_str.c_str());
			zh = &ch[28];										//标幺参数的基准容量
			string Sh_str = zh;
			Sh_str = Sh_str.substr(0, 4);
			gp.Sh = atof(Sh_str.c_str());
			zh = &ch[32];										//Ra
			string Ra_str = zh;
			Ra_str = Ra_str.substr(0, 4);
			gp.Ra = atof(Ra_str.c_str());
			zh = &ch[36];										//Xd1
			string Xd1_str = zh;
			Xd1_str = Xd1_str.substr(0, 5);
			gp.Xd1 = atof(Xd1_str.c_str());
			zh = &ch[41];										//Xq1
			string Xq1_str = zh;
			Xq1_str = Xq1_str.substr(0, 5);
			gp.Xq1 = atof(Xq1_str.c_str());
			zh = &ch[46];										//Xd
			string Xd_str = zh;
			Xd_str = Xd_str.substr(0, 5);
			gp.Xd = atof(Xd_str.c_str());
			zh = &ch[51];										//Xq
			string Xq_str = zh;
			Xq_str = Xq_str.substr(0, 5);
			gp.Xq = atof(Xq_str.c_str());
			zh = &ch[56];										//Td01
			string Td01_str = zh;
			Td01_str = Td01_str.substr(0, 4);
			gp.Td01 = atof(Td01_str.c_str());
			zh = &ch[60];										//Tq01	
			string Tq01_str = zh;
			Tq01_str = Tq01_str.substr(0, 3);
			gp.Tq01 = atof(Tq01_str.c_str());
			zh = &ch[63];										//XL
			string XL_str = zh;
			XL_str = XL_str.substr(0, 5);
			gp.XL = atof(XL_str.c_str());
			zh = &ch[77];										//D
			string D_str = zh;
			D_str = D_str.substr(0, 3);
			gp.D = atof(D_str.c_str());

			gp.Tj = 2 * gp.Ek / gp.Sh;
			//判断同步机模型
			gp.gene_model == "classic";							//经典模型

			para_gene.push_back(gp);		//将同步机参数压栈
		}

		//if (card_type == "MG")					//发电机模型卡...
	}

	//int num_gene = para_gene.size();
	//while (input1.peek() != EOF)			//继续读励磁，PSS，调速器等参数....
	//注意要根据发电机名，vector同个索引要指向同一个发电机

	input1.close();
}



//读取负荷节点参数
void sssa_inputload(char name_st[], vector<LOAD_PARA> &para_load)
{
	char temp[50];

	/////////////////////////读各发电机参数
	strcpy_s(temp, name_st);
	strcat_s(temp, ".swi");
	ifstream input1;

	input1.open(temp);
	if (!input1)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}
	char ch[100];
	string sh;
	char *zh;
	while (input1.getline(ch, 100))			//按行读取文件
	{
		sh = ch;
		//cout << sh << endl;

		string card_type = sh.substr(0, 3);
		card_type.erase(0, card_type.find_first_not_of(" "));	//清除头尾空格
		card_type.erase(card_type.find_last_not_of(" ") + 1);

		//负荷参数
		LOAD_PARA lp;
		if (card_type == "LB")				//负荷模型卡（zip模型）
		{
			lp.card_type = card_type;
			string name = sh.substr(3, 8);						//负荷母线名
			name.erase(0, name.find_first_not_of(" "));
			name.erase(name.find_last_not_of(" ") + 1);
			lp.name = name;
			string vol_str = sh.substr(11, 4);					//负荷母线基准电压
			lp.vol = atof(vol_str.c_str());

			zh = &ch[27];										//ap
			string ap_str = zh; 
			ap_str = ap_str.substr(0, 5);
			lp.ap = atof(ap_str.c_str());
			zh = &ch[32];										//aq
			string aq_str = zh;
			aq_str = aq_str.substr(0, 5);
			lp.aq = atof(aq_str.c_str());
			zh = &ch[37];										//bp
			string bp_str = zh;
			bp_str = bp_str.substr(0, 5);
			lp.bp = atof(bp_str.c_str());
			zh = &ch[42];										//bq
			string bq_str = zh;
			bq_str = bq_str.substr(0, 5);
			lp.bq = atof(bq_str.c_str());
			zh = &ch[47];										//cp
			string cp_str = zh;
			cp_str = cp_str.substr(0, 5);
			lp.cp = atof(cp_str.c_str());
			zh = &ch[52];										//cq
			string cq_str = zh;
			cq_str = cq_str.substr(0, 5);
			lp.cq = atof(cq_str.c_str());

			
			para_load.push_back(lp);
		}	
	}
	LOAD_PARA lp;					//这一段强制给母线2赋值!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	lp.card_type = "LB";										
	lp.name = "母线2";
	lp.ap = lp.aq = 1.0;
	lp.bp = lp.bq = lp.cp = lp.cq = 0.0;
	para_load.push_back(lp);
	input1.close();
}



//读取其他文件
void sssa_inputfile(char name_st[], vector<NODE> &node_table, vector<NODE_DATA> &data_node, vector<GENE_PARA> &para_gene, vector<LOAD_PARA> &para_load)
{
	sssa_inputnode(name_st, node_table);			//读取节点对应表
	sssa_inputPF(name_st, node_table, data_node);	//读取节点潮流结果
	sssa_inputgene(name_st, para_gene);				//读取发电机参数（包括同步机，励磁，PSS，调速器）	
	sssa_inputload(name_st, para_load);				//读取负荷参数
	
}



