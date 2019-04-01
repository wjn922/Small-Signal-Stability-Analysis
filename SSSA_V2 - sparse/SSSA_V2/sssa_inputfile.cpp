#include"head.h"

//=====================��ȡ���ɾ���=======================
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


//==========================��ȡ�ڵ��Ӧ��========================
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


//============================��ȡ���ڵ㳱�����=================================
void sssa_inputPF(char name_st[], int num_gene, int num_node, vector<NODE_DATA> &data_node)
{
	
	char temp[50];

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
			if (knb->name == Fname)	 //�Ͷ�ĸ��
			{
				knb->P += FP;
				knb->Q += FQ;
			}
			if (knb->name == Tname)	 //�ܶ�ĸ��
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

}


//==========================��ȡ���������������ͬ���������ţ�PSS����������===============================
void sssa_inputgene(char name_st[], vector<GENE_PARA> &para_gene, vector<EXCI_PARA> &para_exci, vector<PSS_PARA> &para_pss)
{
	char temp[50];
	char ch[100];
	string sh;
	char *zh;

	/////////////////////////////////////////����ͬ��������///////////////////////////////////////
	strcpy_s(temp, name_st);
	strcat_s(temp, ".swi");
	ifstream input1;

	input1.open(temp);
	if (!input1)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}

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

		if (card_type == "MF" || card_type == "MC")					//�����ģ�Ϳ�
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

			if (gp.Sh == 0.0)
				gp.Sh = base;
			gp.Tj = 2 * gp.Ek / gp.Sh;

			if (card_type == "MC")
				gp.gene_model = "classic";
			else        //MF��
			{
				//�ж�ͬ����ģ��
				if (gp.Td02 != 0.0 && gp.Tq02 != 0.0)		//���д���̬����
				{
					if (gp.Tq01 != 0 && gp.Xq != gp.Xq1)			//������ģ��
						gp.gene_model = "four_winding";
					if (gp.Tq01 == 0 || gp.Xq == gp.Xq1)			//������ģ��
						gp.gene_model = "three_winding";
				}
				else										//��������̬����
				{
					if (gp.Tq01 != 0 && gp.Xq != gp.Xq1)	
						gp.gene_model = "two_axle";					//˫��ģ��
					if (gp.Tq01 == 0 || gp.Xq == gp.Xq1)
						gp.gene_model = "Eq_change";				//Eq'�仯ģ��
				}
				if (gp.Td01 == 999.0)
					gp.gene_model = "Eq_constant";					//Eq'�㶨ģ��
			}
			

			para_gene.push_back(gp);		//��ͬ��������ѹջ
		}

		if (card_type == "MG")					//�����ģ�Ϳ�...
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

			//�ж�ͬ����ģ��
			if (gp.Td02 != 0.0 && gp.Tq02 != 0.0)		//���д���̬����
			{
				if (gp.Tq01 != 0 && gp.Xq != gp.Xq1)			//6�ͣ�Eq',Eq'',Ed',Ed''�仯��6��ģ��
					gp.gene_model = "type6";
				if (gp.Tq01 == 0 || gp.Xq == gp.Xq1)			//3�ͣ�Eq',Eq'',Ed''�仯��5��ģ��
					gp.gene_model = "type3";
			}
			else										//��������̬����
			{
				if (gp.Tq01 != 0 && gp.Xq != gp.Xq1)
					gp.gene_model = "type5";					//5�ͣ�Eq',Ed'�仯��4��ģ��
				if (gp.Tq01 == 0 || gp.Xq == gp.Xq1)
					gp.gene_model = "type2";					//2�ͣ�Eq'�仯��3��ģ��
			}
			if (gp.Td01 == 999.0)
				gp.gene_model = "type1";						//1�ͣ�Eq'�㶨��2��ģ��
			//����ģ�� type0 ������MC��

			para_gene.push_back(gp);		//��ͬ��������ѹջ

		}

		memset(ch, 0, sizeof(ch));	//���ch

		
	}
	input1.close();


	/////////////////////////////////////��������ϵͳ����////////////////////////////////////////////
	int num_gene = para_gene.size();

	input1.open(temp);
	if (!input1)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}

	while (input1.getline(ch, 100))			//���ж�ȡ�ļ�
	{
		sh = ch;
		//cout << sh << endl;

		string card_type = sh.substr(0, 3);
		card_type.erase(0, card_type.find_first_not_of(" "));	//���ͷβ�ո�
		card_type.erase(card_type.find_last_not_of(" ") + 1);

		if (card_type[0] == 'E')
		{
			string name = sh.substr(3, 8);						//���ĸ����
			name.erase(0, name.find_first_not_of(" "));
			name.erase(name.find_last_not_of(" ") + 1);

			int index = 0;	
			for (int k = 0; k < num_gene; k++)			//��̨ͬ�����������
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
			ep.card_type = card_type;							//���ſ�Ƭ��
			ep.exci_model = card_type;							//����ģ����
			ep.name = name;										//ĸ����
			string vol_str = sh.substr(11, 4);					//���ĸ�߻�׼��ѹ
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
		
		//if (card_type[0] == 'F')...����ģ��

		memset(ch, 0, sizeof(ch));	//���ch
	}
	input1.close();

	int num_exci = para_exci.size();
	//��֤�����������Ԫ��vector��size��ͬ
	for (int k = 0; k < num_gene - num_exci; k++)
	{
		EXCI_PARA ep;
		para_exci.push_back(ep);
	}


	/////////////////////////////////////����PSSϵͳ����////////////////////////////////////////////
	input1.open(temp);
	if (!input1)
	{
		cout << "can't open the file:" << temp << endl;
		exit(1);
	}

	while (input1.getline(ch, 100))			//���ж�ȡ�ļ�
	{
		sh = ch;
		//cout << sh << endl;

		string card_type = sh.substr(0, 3);
		card_type.erase(0, card_type.find_first_not_of(" "));	//���ͷβ�ո�
		card_type.erase(card_type.find_last_not_of(" ") + 1);

		if (card_type == "SF" || card_type == "SP" || card_type == "SS" || card_type == "SG")
		{
			string name = sh.substr(3, 8);						//���ĸ����
			name.erase(0, name.find_first_not_of(" "));
			name.erase(name.find_last_not_of(" ") + 1);

			int index = 0;
			for (int k = 0; k < num_gene; k++)					//��̨ͬ�����������
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
			pp.card_type = card_type;							//PSS��Ƭ��
			pp.pss_model = card_type;							//PSSģ����
			pp.name = name;										//ĸ����
			string vol_str = sh.substr(11, 4);					//���ĸ�߻�׼��ѹ
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
			zh = &ch[76];										//Sh(Kqs�Ļ�׼����) ��������SP/SG
			string Sh_str = zh;
			Sh_str = Sh_str.substr(0, 4);
			pp.Sh = atof(Sh_str.c_str());

			para_pss.push_back(pp);
		}
		memset(ch, 0, sizeof(ch));	//���ch
	}
	input1.close();

	//��֤�����������Ԫ��vector��size��ͬ
	int num_pss = para_pss.size();
	for (int k = 0; k < num_gene - num_pss; k++)
	{
		PSS_PARA pp;
		para_pss.push_back(pp);
	}
	//while (input1.peek() != EOF)			//������PSS���������Ȳ���....
	//ע��Ҫ���ݷ��������vectorͬ������Ҫָ��ͬһ�������

	
}



//=========================��ȡ���ɽڵ����=============================
void sssa_inputload(char name_st[], vector<LOAD_PARA> &para_load)
{
	char temp[50];

	/////////////////////////�������ɽڵ����
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
			string vol_str = sh.substr(11,  4);					//����ĸ�߻�׼��ѹ
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
		memset(ch, 0, sizeof(ch));	//���ch
	}
	input1.close();

	string name_temp = name_st;
	if(name_temp == "IEEE90")
	{
		LOAD_PARA lp;					//��һ��ǿ�Ƹ�ĸ��2��ֵ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		lp.card_type = "LB";
		lp.name = "ĸ��2";
		lp.ap = lp.aq = 1.0;
		lp.bp = lp.bq = lp.cp = lp.cq = 0.0;
		para_load.push_back(lp);
	}

}



//��ȡ�����ļ�
void sssa_inputfile(char name_st[], vector<NODE> &node_table, vector<NODE_DATA> &data_node, 
	vector<GENE_PARA> &para_gene, vector<EXCI_PARA> &para_exci, vector<PSS_PARA> &para_pss, vector<LOAD_PARA> &para_load)
{
	sssa_inputnode(name_st, node_table);								//��ȡ�ڵ��Ӧ��
	sssa_inputgene(name_st, para_gene, para_exci, para_pss);			//��ȡ���������������ͬ���������ţ�PSS����������	
	sssa_inputload(name_st, para_load);									//��ȡ���ɲ���
	int num_gene = para_gene.size();
	int num_node = node_table.size();
	sssa_inputPF(name_st, num_gene, num_node, data_node);				//��ȡ�ڵ㳱�����
	
}



