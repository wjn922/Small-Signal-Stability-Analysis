#include"head.h"

//��ȡ���ɾ���
void sssa_inputY(char name_st[], double **Ymatrix)
{
	char temp[50];
	int ry, cy, row, col;
	double G, B;
	char dot;

	//�����ɾ���
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
		Ymatrix[row + 1][col + 1] = G;
		Ymatrix[row + 1][col + 2] = -B;
		Ymatrix[row + 2][col + 1] = B;
		Ymatrix[row + 2][col + 2] = G;
	}
}


//��ȡ�ڵ��Ӧ��
void sssa_inputnode(char name_st[], vector<NODE> &node_table)
{
	char temp[50];

	///////////////////////////���ڵ��Ӧ��
	strcpy_s(temp, name_st);
	strcat_s(temp, "�ڵ��Ӧ��.txt");
	ifstream input;

	input.open(temp);
	if (!input)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}
	while (!input.eof())			//�ļ�δ��ȡ��	
	{
		NODE node;
		char dot;
		input >> node.no >> dot >> node.name >> node.vol;
		node_table.push_back(node);
	}
	node_table.pop_back();		//��ֱ��!input.eof()���һ�л�������飬�������һ��Ԫ��
	input.close();
}


//��ȡ���ڵ㳱�����
void sssa_inputPF(char name_st[], vector<NODE> node_table, vector<NODE_DATA> &data_node)
{
	
	char temp[50];

	int num_node = node_table.size();
	int num_gene = 0;
	vector<NODE>::iterator kn;
	for (kn = node_table.begin(); kn != node_table.end(); kn++)
	{
		string name = kn->name;
		if (name.find("�����") != name.npos)		//�Ƿ����
			num_gene++;
	}

	strcpy_s(temp, name_st);
	strcat_s(temp, "�������.txt");
	ifstream input;

	input.open(temp);
	if (!input)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}

	string line;
	//�ڵ��ѹ��Ϣ
	getline(input, line);
	getline(input, line);		//ǰ����
	for(int i = 0; i<num_node; i++)
	{
		getline(input, line);
		stringstream word(line);
		NODE_DATA nb;
		word >> nb.no >> nb.name >> nb.vol_base >> nb.vol_pu >> nb.vol >> nb.vol_ang;
		nb.P = nb.Q = 0.0;
		data_node.push_back(nb);
	}

	//�����������Ϣ
	getline(input, line);
	getline(input, line);		//������
	for(int i = 0; i < num_gene; i++)
		getline(input, line);

	//��·/��ѹ��������Ϣ
	getline(input, line);
	getline(input, line);		//������
	while (getline(input, line))
	{
		int Fno, Tno;						//�Ͷ˱�ţ��ܶ˱��
		string Fname, Tname;				//�Ͷ�ĸ�������ܶ�ĸ����
		double Fvol, Tvol;					//�Ͷ�ĸ�߻�׼��ѹ�� �ܶ�ĸ�߻�׼��ѹ
		double FP, FQ, TP, TQ;				//�Ͷ��й����Ͷ��޹����ܶ��й����ܶ��޹�

		stringstream word(line);
		word >> Fno >> Fname >> Fvol >> Tno >> Tname >> Tvol >> FP >> FQ >> TP >> TQ;

		vector<NODE_DATA>::iterator knb;
		for (knb = data_node.begin(); knb != data_node.end(); knb++)
		{
			if (knb->name == Fname)	 //�����ĸ��
			{
				knb->P += FP;
				knb->Q += FQ;
			}
			if (knb->name == Tname)	 //�����ĸ��
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
		if (knb->P > 0)		knb->type = -1;			//������ڵ�
		else if (knb->P < 0)	knb->type = 1;		//���ɽڵ�
		else	knb->type = 0;						//����ڵ�
	}
	
	
	/*
	/////////////////////////������������нڵ�ĸ�ߵ�ѹ
	strcpy_s(temp, name_st);
	strcat_s(temp, "�������_Bus.txt");
	ifstream input1;
	
	input1.open(temp);
	if (!input1)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}
	while (!input1.eof())			//�ļ�δ��ȡ��	
	{
		NODE_DATA nb;
		input1 >> nb.no >> nb.name >> nb.vol_base >> nb.vol_pu >> nb.vol >> nb.vol_ang;
		nb.P = nb.Q = 0.0;
		data_node.push_back(nb);
	}
	data_node.pop_back();		//��ֱ��!input.eof()���һ�л�������飬�������һ��Ԫ��
	input1.close();


	//////////////////////////������������ڵ㹦��
	strcpy_s(temp, name_st);
	strcat_s(temp, "�������_Line.txt");
	ifstream input2;

	input2.open(temp);
	if (!input2)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}
	string line;
	while (getline(input2, line))			//���ж�ȡ�ļ�
	{
		int Fno, Tno;						//�Ͷ˱�ţ��ܶ˱��
		string Fname, Tname;				//�Ͷ�ĸ�������ܶ�ĸ����
		double Fvol, Tvol;					//�Ͷ�ĸ�߻�׼��ѹ�� �ܶ�ĸ�߻�׼��ѹ
		double FP, FQ, TP, TQ;				//�Ͷ��й����Ͷ��޹����ܶ��й����ܶ��޹�

		stringstream word(line);
		word >> Fno >> Fname >> Fvol >> Tno >> Tname >> Tvol >> FP >> FQ >> TP >> TQ;

		vector<NODE_DATA>::iterator knb;
		for (knb = data_node.begin(); knb != data_node.end(); knb++)
		{
			if (knb->name == Fname)	 //�����ĸ��
			{
				knb->P += FP;
				knb->Q += FQ;
			}
			if (knb->name == Tname)	 //�����ĸ��
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
		if (knb->P > 0)		knb->type = -1;			//������ڵ�
		else if (knb->P < 0)	knb->type = 1;		//���ɽڵ�
		else	knb->type = 0;						//����ڵ�
	}
	*/

}


//��ȡ���������������ͬ���������ţ�PSS����������
void sssa_inputgene(char name_st[], vector<GENE_PARA> &para_gene)
{
	char temp[50];

	/////////////////////////�������������
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
	while (input1.getline(ch, 100))			//���ж�ȡ�ļ�
	{
		sh = ch;
		//cout << sh << endl;
									
		string card_type = sh.substr(0, 3);
		card_type.erase(0, card_type.find_first_not_of(" "));	//���ͷβ�ո�
		card_type.erase(card_type.find_last_not_of(" ") + 1);

		//ͬ��������
		GENE_PARA gp;

		if (card_type == "M")					//���������̬����ģ�Ϳ������ڷ�����ĵ�һ��
		{
			string name = sh.substr(3, 8);						//���ĸ����
			name.erase(0, name.find_first_not_of(" "));
			name.erase(name.find_last_not_of(" ") + 1);
			gp.name = name;
			string vol_str = sh.substr(11, 4);					//���ĸ�߻�׼��ѹ
			gp.vol = atof(vol_str.c_str());
			string capacity_str = sh.substr(16, 5);				//��������
			gp.capacity = atof(capacity_str.c_str());
			string power_factor_str = sh.substr(22, 3);			//�����������
			gp.power_factor = atof(power_factor_str.c_str());
			string type = sh.substr(30, 2);						//�������
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

		if (card_type == "MF")					//�����ģ�Ϳ�
		{
			gp.card_type = card_type;							//�������Ƭ��
			string name = sh.substr(3, 8);						//���ĸ����
			name.erase(0, name.find_first_not_of(" "));
			name.erase(name.find_last_not_of(" ") + 1);
			gp.name = name;

			zh = &ch[11];										//���ĸ�߻�׼��ѹ
			string vol_str = zh;
			vol_str = vol_str.substr(0, 4);
			gp.vol = atof(vol_str.c_str());
			zh = &ch[16];										//���������
			string Ek_str = zh;
			Ek_str = Ek_str.substr(0, 6);
			gp.Ek = atof(Ek_str.c_str());
			zh = &ch[28];										//���۲����Ļ�׼����
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
			//�ж�ͬ����ģ��
			if (gp.Td02 != 0.0 && gp.Tq02 != 0.0)		//���д���̬����
			{
				if (gp.Tq01 != 0 && gp.Xq != gp.Xq1)			//������ģ��
					gp.gene_model = "four_winding";
				if (gp.Tq01 == 0 || gp.Xq == gp.Xq1)			//������ģ��
					gp.gene_model = "three_winding";
			}
			else
			{
				if (gp.Tq01 != 0 && gp.Xq != gp.Xq1)	//��������̬����
					gp.gene_model = "two_axle";					//˫��ģ��
				if (gp.Tq01 == 0 || gp.Xq == gp.Xq1)
					gp.gene_model = "Eq_change";				//Eq'�仯ģ��
			}
			if (gp.Td01 == 999.0)
				gp.gene_model = "Eq_constant";					//Eq'�㶨ģ��

			para_gene.push_back(gp);		//��ͬ��������ѹջ
		}

		if (card_type == "MC")					//�����ģ�Ϳ�(����ģ��)
		{
			gp.card_type = card_type;							//�������Ƭ��
			string name = sh.substr(3, 8);						//���ĸ����
			name.erase(0, name.find_first_not_of(" "));
			name.erase(name.find_last_not_of(" ") + 1);
			gp.name = name;
			
			zh = &ch[11];										//���ĸ�߻�׼��ѹ
			string vol_str = zh;
			vol_str = vol_str.substr(0, 4);
			gp.vol = atof(vol_str.c_str());
			zh = &ch[16];										//���������
			string Ek_str = zh;
			Ek_str = Ek_str.substr(0, 6);
			gp.Ek = atof(Ek_str.c_str());
			zh = &ch[28];										//���۲����Ļ�׼����
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
			//�ж�ͬ����ģ��
			gp.gene_model == "classic";							//����ģ��

			para_gene.push_back(gp);		//��ͬ��������ѹջ
		}

		//if (card_type == "MG")					//�����ģ�Ϳ�...
	}

	//int num_gene = para_gene.size();
	//while (input1.peek() != EOF)			//���������ţ�PSS���������Ȳ���....
	//ע��Ҫ���ݷ��������vectorͬ������Ҫָ��ͬһ�������

	input1.close();
}



//��ȡ���ɽڵ����
void sssa_inputload(char name_st[], vector<LOAD_PARA> &para_load)
{
	char temp[50];

	/////////////////////////�������������
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
	while (input1.getline(ch, 100))			//���ж�ȡ�ļ�
	{
		sh = ch;
		//cout << sh << endl;

		string card_type = sh.substr(0, 3);
		card_type.erase(0, card_type.find_first_not_of(" "));	//���ͷβ�ո�
		card_type.erase(card_type.find_last_not_of(" ") + 1);

		//���ɲ���
		LOAD_PARA lp;
		if (card_type == "LB")				//����ģ�Ϳ���zipģ�ͣ�
		{
			lp.card_type = card_type;
			string name = sh.substr(3, 8);						//����ĸ����
			name.erase(0, name.find_first_not_of(" "));
			name.erase(name.find_last_not_of(" ") + 1);
			lp.name = name;
			string vol_str = sh.substr(11, 4);					//����ĸ�߻�׼��ѹ
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
	LOAD_PARA lp;					//��һ��ǿ�Ƹ�ĸ��2��ֵ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	lp.card_type = "LB";										
	lp.name = "ĸ��2";
	lp.ap = lp.aq = 1.0;
	lp.bp = lp.bq = lp.cp = lp.cq = 0.0;
	para_load.push_back(lp);
	input1.close();
}



//��ȡ�����ļ�
void sssa_inputfile(char name_st[], vector<NODE> &node_table, vector<NODE_DATA> &data_node, vector<GENE_PARA> &para_gene, vector<LOAD_PARA> &para_load)
{
	sssa_inputnode(name_st, node_table);			//��ȡ�ڵ��Ӧ��
	sssa_inputPF(name_st, node_table, data_node);	//��ȡ�ڵ㳱�����
	sssa_inputgene(name_st, para_gene);				//��ȡ���������������ͬ���������ţ�PSS����������	
	sssa_inputload(name_st, para_load);				//��ȡ���ɲ���
	
}



