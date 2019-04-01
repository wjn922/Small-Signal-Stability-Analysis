#include"head.h"

//=====================读取导纳矩阵=======================
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


//==========================读取节点对应表========================
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


//============================读取各节点潮流结果=================================
void sssa_inputPF(char name_st[], int num_gene, int num_node, vector<NODE_DATA> &data_node)
{
	
	char temp[50];

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
			if (knb->name == Fname)	 //送端母线
			{
				knb->P += FP;
				knb->Q += FQ;
			}
			if (knb->name == Tname)	 //受端母线
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

}


//==========================读取发电机参数（包括同步机，励磁，PSS，调速器）===============================
void sssa_inputgene(char name_st[], vector<GENE_PARA> &para_gene, vector<EXCI_PARA> &para_exci, vector<PSS_PARA> &para_pss)
{
	char temp[50];
	char ch[100];
	string sh;
	char *zh;

	/////////////////////////////////////////读各同步机参数///////////////////////////////////////
	strcpy_s(temp, name_st);
	strcat_s(temp, ".swi");
	ifstream input1;

	input1.open(temp);
	if (!input1)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}

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

		if (card_type == "MF" || card_type == "MC")					//发电机模型卡
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

			if (gp.Sh == 0.0)
				gp.Sh = base;
			gp.Tj = 2 * gp.Ek / gp.Sh;

			if (card_type == "MC")
				gp.gene_model = "classic";
			else        //MF卡
			{
				//判断同步机模型
				if (gp.Td02 != 0.0 && gp.Tq02 != 0.0)		//含有次暂态参数
				{
					if (gp.Tq01 != 0 && gp.Xq != gp.Xq1)			//四绕组模型
						gp.gene_model = "four_winding";
					if (gp.Tq01 == 0 || gp.Xq == gp.Xq1)			//三绕组模型
						gp.gene_model = "three_winding";
				}
				else										//不含次暂态参数
				{
					if (gp.Tq01 != 0 && gp.Xq != gp.Xq1)	
						gp.gene_model = "two_axle";					//双轴模型
					if (gp.Tq01 == 0 || gp.Xq == gp.Xq1)
						gp.gene_model = "Eq_change";				//Eq'变化模型
				}
				if (gp.Td01 == 999.0)
					gp.gene_model = "Eq_constant";					//Eq'恒定模型
			}
			

			para_gene.push_back(gp);		//将同步机参数压栈
		}

		if (card_type == "MG")					//发电机模型卡...
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
			zh = &ch[63];										//N
			string N_str = zh;
			N_str = N_str.substr(0, 5);
			gp.N = atof(N_str.c_str());
			zh = &ch[68];										//A
			string A_str = zh;
			A_str = A_str.substr(0, 5);
			gp.A = atof(A_str.c_str());
			zh = &ch[73];										//B
			string B_str = zh;
			B_str = B_str.substr(0, 5);
			gp.B = atof(B_str.c_str());
			zh = &ch[77];										//D
			string D_str = zh;
			D_str = D_str.substr(0, 3);
			gp.D = atof(D_str.c_str());

			if (gp.Sh == 0.0)
				gp.Sh = base;
			gp.Tj = 2 * gp.Ek / gp.Sh;

			//判断同步机模型
			if (gp.Td02 != 0.0 && gp.Tq02 != 0.0)		//含有次暂态参数
			{
				if (gp.Tq01 != 0 && gp.Xq != gp.Xq1)			//6型，Eq',Eq'',Ed',Ed''变化的6阶模型
					gp.gene_model = "type6";
				if (gp.Tq01 == 0 || gp.Xq == gp.Xq1)			//3型，Eq',Eq'',Ed''变化的5阶模型
					gp.gene_model = "type3";
			}
			else										//不含次暂态参数
			{
				if (gp.Tq01 != 0 && gp.Xq != gp.Xq1)
					gp.gene_model = "type5";					//5型，Eq',Ed'变化的4阶模型
				if (gp.Tq01 == 0 || gp.Xq == gp.Xq1)
					gp.gene_model = "type2";					//2型，Eq'变化的3阶模型
			}
			if (gp.Td01 == 999.0)
				gp.gene_model = "type1";						//1型，Eq'恒定的2阶模型
			//经典模型 type0 还是用MC卡

			para_gene.push_back(gp);		//将同步机参数压栈

		}

		memset(ch, 0, sizeof(ch));	//清空ch

		
	}
	input1.close();


	/////////////////////////////////////读各励磁系统参数////////////////////////////////////////////
	int num_gene = para_gene.size();

	input1.open(temp);
	if (!input1)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}

	while (input1.getline(ch, 100))			//按行读取文件
	{
		sh = ch;
		//cout << sh << endl;

		string card_type = sh.substr(0, 3);
		card_type.erase(0, card_type.find_first_not_of(" "));	//清除头尾空格
		card_type.erase(card_type.find_last_not_of(" ") + 1);

		if (card_type[0] == 'E')
		{
			string name = sh.substr(3, 8);						//电机母线名
			name.erase(0, name.find_first_not_of(" "));
			name.erase(name.find_last_not_of(" ") + 1);

			int index = 0;	
			for (int k = 0; k < num_gene; k++)			//找同台发电机的索引
			{
				if (para_gene[k].name == name)
					break;
				index++;
			}
			for (int k = 0; k < index - para_exci.size(); k++)
			{
				EXCI_PARA ep;
				para_exci.push_back(ep);
			}

			EXCI_PARA ep;
			ep.card_type = card_type;							//励磁卡片名
			ep.exci_model = card_type;							//励磁模型名
			ep.name = name;										//母线名
			string vol_str = sh.substr(11, 4);					//电机母线基准电压
			ep.vol = atof(vol_str.c_str());
			zh = &ch[16];										//Tr
			string Tr_str = zh;
			Tr_str = Tr_str.substr(0, 4);
			ep.Tr = atof(Tr_str.c_str());
			zh = &ch[20];										//Ka
			string Ka_str = zh;
			Ka_str = Ka_str.substr(0, 5);
			ep.Ka = atof(Ka_str.c_str());
			zh = &ch[25];										//Ta
			string Ta_str = zh;
			Ta_str = Ta_str.substr(0, 4);
			ep.Ta = atof(Ta_str.c_str());
			zh = &ch[29];										//Ta1
			string Ta1_str = zh;
			Ta1_str = Ta1_str.substr(0, 4);
			ep.Ta1 = atof(Ta1_str.c_str());
			zh = &ch[33];										//Vmult
			string Vmult_str = zh;
			Vmult_str = Vmult_str.substr(0, 4);
			ep.Vmult = atof(Vmult_str.c_str());
			zh = &ch[37];										//Ke
			string Ke_str = zh;
			Ke_str = Ke_str.substr(0, 4);
			ep.Ke = atof(Ke_str.c_str());
			zh = &ch[41];										//Te
			string Te_str = zh;
			Te_str = Te_str.substr(0, 4);
			ep.Te = atof(Te_str.c_str());
			zh = &ch[45];										//Se0.75
			string Sei_str = zh;
			Sei_str = Sei_str.substr(0, 4);
			ep.Sei = atof(Sei_str.c_str());
			zh = &ch[49];										//Se
			string Se_str = zh;
			Se_str = Se_str.substr(0, 4);
			ep.Se = atof(Se_str.c_str());
			zh = &ch[53];										//Efqmin
			string Efqmin_str = zh;
			Efqmin_str = Efqmin_str.substr(0, 5);
			ep.Efqmin = atof(Efqmin_str.c_str());
			zh = &ch[58];										//Efqmax
			string Efqmax_str = zh;
			Efqmax_str = Efqmax_str.substr(0, 4);
			ep.Efqmax = atof(Efqmax_str.c_str());
			zh = &ch[62];										//Kf
			string Kf_str = zh;
			Kf_str = Kf_str.substr(0, 4);
			ep.Kf = atof(Kf_str.c_str());
			zh = &ch[66];										//Tf
			string Tf_str = zh;
			Tf_str = Tf_str.substr(0, 4);
			ep.Tf = atof(Tf_str.c_str());
			zh = &ch[70];										//XL
			string XL_str = zh;
			XL_str = XL_str.substr(0, 5);
			ep.XL = atof(XL_str.c_str());
			zh = &ch[75];										//Tf1
			string Tf1_str = zh;
			Tf1_str = Tf1_str.substr(0, 4);
			ep.Tf1 = atof(Tf1_str.c_str());
			
			//ep.Ta = ep.Ta - ep.Ta1;
			para_exci.push_back(ep);
		}
		
		//if (card_type[0] == 'F')...其他模型

		memset(ch, 0, sizeof(ch));	//清空ch
	}
	input1.close();

	int num_exci = para_exci.size();
	//保证发电机各部分元件vector的size相同
	for (int k = 0; k < num_gene - num_exci; k++)
	{
		EXCI_PARA ep;
		para_exci.push_back(ep);
	}


	/////////////////////////////////////读各PSS系统参数////////////////////////////////////////////
	input1.open(temp);
	if (!input1)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}

	while (input1.getline(ch, 100))			//按行读取文件
	{
		sh = ch;
		//cout << sh << endl;

		string card_type = sh.substr(0, 3);
		card_type.erase(0, card_type.find_first_not_of(" "));	//清除头尾空格
		card_type.erase(card_type.find_last_not_of(" ") + 1);

		if (card_type == "SF" || card_type == "SP" || card_type == "SS" || card_type == "SG")
		{
			string name = sh.substr(3, 8);						//电机母线名
			name.erase(0, name.find_first_not_of(" "));
			name.erase(name.find_last_not_of(" ") + 1);

			int index = 0;
			for (int k = 0; k < num_gene; k++)					//找同台发电机的索引
			{
				if (para_gene[k].name == name)
					break;
				index++;
			}
			for (int k = 0; k < index - para_pss.size(); k++)
			{
				PSS_PARA pp;
				para_pss.push_back(pp);
			}

			PSS_PARA pp;
			pp.card_type = card_type;							//PSS卡片名
			pp.pss_model = card_type;							//PSS模型名
			pp.name = name;										//母线名
			string vol_str = sh.substr(11, 4);					//电机母线基准电压
			pp.vol = atof(vol_str.c_str());
			zh = &ch[16];										//Kqv
			string Kqv_str = zh;
			Kqv_str = Kqv_str.substr(0, 4);
			pp.Kqv = atof(Kqv_str.c_str());
			zh = &ch[20];										//Tqv
			string Tqv_str = zh;
			Tqv_str = Tqv_str.substr(0, 3);
			pp.Tqv = atof(Tqv_str.c_str());
			zh = &ch[23];										//Kqs
			string Kqs_str = zh;
			Kqs_str = Kqs_str.substr(0, 4);
			pp.Kqs = atof(Kqs_str.c_str());
			zh = &ch[27];										//Tqs
			string Tqs_str = zh;
			Tqs_str = Tqs_str.substr(0, 3);
			pp.Tqs = atof(Tqs_str.c_str());
			zh = &ch[30];										//Tq
			string Tq_str = zh;
			Tq_str = Tq_str.substr(0, 4);
			pp.Tq = atof(Tq_str.c_str());
			zh = &ch[34];										//Tq1
			string Tq1_str = zh;
			Tq1_str = Tq1_str.substr(0, 4);
			pp.Tq1 = atof(Tq1_str.c_str());
			zh = &ch[38];										//Tq11
			string Tq11_str = zh;
			Tq11_str = Tq11_str.substr(0, 4);
			pp.Tq11 = atof(Tq11_str.c_str());
			zh = &ch[42];										//Tq2
			string Tq2_str = zh;
			Tq2_str = Tq2_str.substr(0, 4);
			pp.Tq2 = atof(Tq2_str.c_str());
			zh = &ch[46];										//Tq21
			string Tq21_str = zh;
			Tq21_str = Tq21_str.substr(0, 4);
			pp.Tq21 = atof(Tq21_str.c_str());
			zh = &ch[50];										//Tq3
			string Tq3_str = zh;
			Tq3_str = Tq3_str.substr(0, 4);
			pp.Tq3 = atof(Tq3_str.c_str());
			zh = &ch[54];										//Tq31
			string Tq31_str = zh;
			Tq31_str = Tq31_str.substr(0, 4);
			pp.Tq31 = atof(Tq31_str.c_str());
			zh = &ch[58];										//Vsmax
			string Vsmax_str = zh;
			Vsmax_str = Vsmax_str.substr(0, 4);
			pp.Vsmax = atof(Vsmax_str.c_str());
			zh = &ch[62];										//Vcutoff
			string Vcutoff_str = zh;
			Vcutoff_str = Vcutoff_str.substr(0, 4);
			pp.Vcutoff = atof(Vcutoff_str.c_str());
			zh = &ch[66];										//Vslow
			string Vslow_str = zh;
			Vslow_str = Vslow_str.substr(0, 2);
			pp.Vslow = atof(Vslow_str.c_str());
			zh = &ch[76];										//Sh(Kqs的基准容量) 仅适用于SP/SG
			string Sh_str = zh;
			Sh_str = Sh_str.substr(0, 4);
			pp.Sh = atof(Sh_str.c_str());

			para_pss.push_back(pp);
		}
		memset(ch, 0, sizeof(ch));	//清空ch
	}
	input1.close();

	//保证发电机各部分元件vector的size相同
	int num_pss = para_pss.size();
	for (int k = 0; k < num_gene - num_pss; k++)
	{
		PSS_PARA pp;
		para_pss.push_back(pp);
	}
	//while (input1.peek() != EOF)			//继续读PSS，调速器等参数....
	//注意要根据发电机名，vector同个索引要指向同一个发电机

	
}



//=========================读取负荷节点参数=============================
void sssa_inputload(char name_st[], vector<LOAD_PARA> &para_load)
{
	char temp[50];

	/////////////////////////读各负荷节点参数
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
			string vol_str = sh.substr(11,  4);					//负荷母线基准电压
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
		memset(ch, 0, sizeof(ch));	//清空ch
	}
	input1.close();

	string name_temp = name_st;
	if(name_temp == "IEEE90")
	{
		LOAD_PARA lp;					//这一段强制给母线2赋值!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		lp.card_type = "LB";
		lp.name = "母线2";
		lp.ap = lp.aq = 1.0;
		lp.bp = lp.bq = lp.cp = lp.cq = 0.0;
		para_load.push_back(lp);
	}

}



//读取其他文件
void sssa_inputfile(char name_st[], vector<NODE> &node_table, vector<NODE_DATA> &data_node, 
	vector<GENE_PARA> &para_gene, vector<EXCI_PARA> &para_exci, vector<PSS_PARA> &para_pss, vector<LOAD_PARA> &para_load)
{
	sssa_inputnode(name_st, node_table);								//读取节点对应表
	sssa_inputgene(name_st, para_gene, para_exci, para_pss);			//读取发电机参数（包括同步机，励磁，PSS，调速器）	
	sssa_inputload(name_st, para_load);									//读取负荷参数
	int num_gene = para_gene.size();
	int num_node = node_table.size();
	sssa_inputPF(name_st, num_gene, num_node, data_node);				//读取节点潮流结果
	
}



