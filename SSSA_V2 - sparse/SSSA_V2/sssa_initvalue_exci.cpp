//��ֵ���� ����ϵͳ
#include"head.h"

void sssa_initvalue_exci(char name_st[], vector<NODE_DATA> data_node, vector<EXCI_PARA> para_exci, vector<GENE_INIT> init_gene, vector<EXCI_INIT> &init_exci)
{
	double Se, Ke, Ka;		//����ϵͳ����
	double Efq0, V0;

	int num_gene = para_exci.size();
	vector<NODE_DATA>::iterator knd;

	for (int k = 0; k < num_gene; k++)	//ÿ�ν���һ̨������ĳ�ֵ����
	{
		string exci_model = para_exci[k].exci_model;

		if (exci_model == "EA")
		{
			EXCI_INIT ei;

			for (knd = data_node.begin(); knd != data_node.end(); knd++)
			{
				if (knd->name == para_exci[k].name)
				{
					V0 = knd->vol_pu;
					break;
				}
			}

			Se = para_exci[k].Se;
			Ke = para_exci[k].Ke;
			Ka = para_exci[k].Ka;

			Efq0 = init_gene[k].Efq0;

			//��ʼ��ֵ����
			ei.Vr0 = (Se + Ke) * Efq0;
			ei.Vr10 = ei.Vr0;
			ei.Vf0 = 0.0;
			ei.Vm0 = V0;
			ei.Vref0 = ei.Vm0 + ei.Vr0 / Ka;

			init_exci.push_back(ei);

		}

		//if(exci_model == "F")   //��������ϵͳģ�ͳ�ֵ����

		else     //û������
		{
			EXCI_INIT ei;
			init_exci.push_back(ei);
		}
	}

	/*
	//�洢��ֵ������
	char temp[50];
	strcpy_s(temp, name_st);
	strcat_s(temp, "���ų�ֵ���.init");
	ofstream out;
	out.open(temp);

	vector<EXCI_INIT>::iterator kei;
	for (kei = init_exci.begin(); kei != init_exci.end(); kei++) {
		out << "Vr0 = " << kei->Vr0 << endl;
		out << "Vr10 = " << kei->Vr10 << endl;
		out << "Vf0 = " << kei->Vf0 << endl;
		out << "Vm0 = " << kei->Vm0 << endl;
		out << "Vref0 = " << kei->Vref0 << endl << endl;
	}
	*/
	
}