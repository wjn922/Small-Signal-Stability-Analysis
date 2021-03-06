// SSSA_V2.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include <iostream>
#include"head.h"
#include"head.cpp"
#include"matrix_operator.cpp"
#include"sssa_inputfile.cpp"
#include"sssa_dimension.cpp"
#include"sssa_genesystem.cpp"
#include"sssa_networksystem.cpp"
//#include"sssa_connect_system.cpp"
#include"pf_compute.cpp"


int main()
{
	char name_st[50];
	strcpy_s(name_st, "IEEE90");			//系统名字

	cout << "start" << endl << endl;

	//======================读取系统参数数据=============================
	vector<NODE> node_table;				//存储节点对应表
	vector<NODE_DATA> data_node;			//存储所有节点的潮流结果
	vector<GENE_PARA> para_gene;			//存储各台发电机同步机参数
	vector<EXCI_PARA> para_exci;			//存储各台发电机励磁参数
	vector<PSS_PARA> para_pss;				//存储各台发电机PSS参数
	vector<LOAD_PARA> para_load;			//存储负荷参数
	sssa_inputfile(name_st, node_table, data_node, para_gene, para_exci, para_pss, para_load);

	int num_node, num_gene, nd, ng;	//系统节点数，发电机节点数，2*系统节点数，2*发电机节点数
	num_gene = para_gene.size();
	num_node = node_table.size();
	nd = 2 * num_node;
	ng = 2 * num_gene;

	//同步机参数
	cout << "节点数：" << num_node << endl;
	cout << "发电机节点数：" << num_gene << endl << endl;
	vector<GENE_PARA>::iterator kgp;
	for (kgp = para_gene.begin(); kgp != para_gene.end(); kgp++)
	{
		cout << kgp->card_type << " " << kgp->name << " " << kgp->vol << " " << kgp->capacity << " " << kgp->power_factor
			<< " " << kgp->type << " " << kgp->Xd2 << " " << kgp->Xq2 << " " << kgp->Td02 << " " << kgp->Tq02
			<< " " << kgp->Ek << " " << kgp->Sh << " " << kgp->Ra << " " << kgp->Xd1 << " " << kgp->Xq1 << " " << kgp->Xd << " " << kgp->Xq
			<< " " << kgp->Td01 << " " << kgp->Tq01 << " "  << kgp->Tj << " " << kgp->D << " " << kgp->gene_model;
		cout << endl;
	}
	cout << endl;


	//节点潮流数据
	vector<NODE_DATA>::iterator knb;
	for (knb = data_node.begin(); knb != data_node.end(); knb++)
	{
		cout << knb->name << "\t" << setw(20) << knb->vol_pu << setw(20) << knb->vol_ang << setw(20) << knb->P << setw(20) << knb->Q << endl;
	}
	cout << endl;


	//读取导纳矩阵
	double **Ymatrix = dmatrix(nd + 1, nd + 1);
	for (int i = 1; i <= nd; i++) {
		for (int j = 1; j <= nd; j++)	Ymatrix[i][j] = 0.0;		//初始化
	}
	sssa_inputY(name_st, Ymatrix);

	pf_compute(node_table, data_node, Ymatrix);			//这里重新计算了潮流功率！！！！！！！1

	//节点潮流数据
	for (knb = data_node.begin(); knb != data_node.end(); knb++)
	{
		cout << knb->name << "\t" << setw(20) << knb->vol_pu << setw(20) << knb->vol_ang << setw(20) << knb->P << setw(20) << knb->Q << endl;
	}
	cout << endl;
	

	//=======================确定发电机系统模型维数===============================
	vector<int> gdim, gl_gene, gl_exci, gl_pss;		//存储系统维数, gl_exci, gl_pss...
	vector<string> gname;
	sssa_dimension(gdim, gname, gl_gene, gl_exci, gl_pss, para_gene, para_exci, para_pss);
	

	//=====================形成发电机系统线性化模型=========================
	//int dim_gene_all = gdim[gdim.size() - 1];						//所有发电机的总维数
	//行压缩格式存储矩阵
	vector<int> Trow;
	vector<int> Tcol;
	vector<double> Tvalue;
	vector<int> Jrow;
	vector<int> Jcol;
	vector<double> Jvalue;

	vector<GENE_INIT> init_gene;
	vector<EXCI_INIT> init_exci;
	sssa_genesystem(name_st, node_table, data_node, para_gene, init_gene, para_exci, init_exci, para_pss,
		gdim, gl_gene, gl_exci, gl_pss, Trow, Tcol, Tvalue, Jrow, Jcol, Jvalue);
	
	//===========================形成网络方程=================================
	sssa_networksystem(name_st, node_table, data_node, para_load, para_gene, gdim, gl_gene, Trow, Tcol, Tvalue, Jrow, Jcol, Jvalue);

	//////=====================================输出行压缩存储格式矩阵=============================/////////
	char temp[20];
	strcpy_s(temp, name_st);
	strcat_s(temp, "_Trow.txt");
	ofstream out1;
	out1.open(temp);
	for (int i = 0; i < Trow.size(); i++)
	{
		out1 << Trow[i] << '\t';
	}
	out1 << endl;
	out1.close();

	strcpy_s(temp, name_st);
	strcat_s(temp, "_Tcol.txt");
	ofstream out2;
	out2.open(temp);
	for (int j = 0; j < Tcol.size(); j++)
	{
		out2 << Tcol[j] << '\t';
	}
	out2 << endl;
	out2.close();

	strcpy_s(temp, name_st);
	strcat_s(temp, "_Tvalue.txt");
	ofstream out3;
	out3.open(temp);
	out3.setf(ios::fixed);
	out3.precision(16);
	for (int k = 0; k < Tvalue.size(); k++)
	{
		out3 << Tvalue[k] << '\t';
	}
	out3 << endl;
	out3.close();


	/////////////////
	strcpy_s(temp, name_st);
	strcat_s(temp, "_Jrow.txt");
	ofstream out4;
	out4.open(temp);
	for (int i = 0; i < Jrow.size(); i++)
	{
		out4 << Jrow[i] << '\t';
	}
	out4 << endl;
	out4.close();

	strcpy_s(temp, name_st);
	strcat_s(temp, "_Jcol.txt");
	ofstream out5;
	out5.open(temp);
	for (int j = 0; j < Jcol.size(); j++)
	{
		out5 << Jcol[j] << '\t';
	}
	out5 << endl;

	strcpy_s(temp, name_st);
	strcat_s(temp, "_Jvalue.txt");
	ofstream out6;
	out6.open(temp);
	out6.setf(ios::fixed);
	out6.precision(16);
	for (int k = 0; k < Jvalue.size(); k++)
	{
		out6 << Jvalue[k] << '\t';
	}
	out6 << endl;
	out6.close();

	
	cout << "end" << endl;

	/*
	//==========================连接全系统线性化模型============================
	double **T = dmatrix(dim_gene_all + nd + 1, dim_gene_all + nd + 1);
	double **J = dmatrix(dim_gene_all + nd + 1, dim_gene_all + nd + 1);	//全系统的线性化模型T，J矩阵
	for (int i = 1; i <= dim_gene_all + nd; i++) {
		for (int j = 1; j <= dim_gene_all + nd; j++)	T[i][j] = J[i][j] = 0.0;		//初始化
	}
	sssa_connect_system(name_st, Tgsys, Jgsys, Ysystem, node_table, init_gene, para_gene, gdim, gl_gene, T, J);

	cout << "end" << endl;
	*/
	



	/*   测试输出
	//同步机参数
	cout << num_node << " " << num_gene << endl;
	vector<GENE_PARA>::iterator kgp;
	for (kgp = para_gene.begin(); kgp != para_gene.end(); kgp++)
	{
		cout << kgp->card_type << " " << kgp->name << " " << kgp->vol << " " << kgp->capacity << " " << kgp->power_factor
			<< " " << kgp->type << " " << kgp->Xd2 << " " << kgp->Xq2 << " " << kgp->Td02 << " " << kgp->Tq02
			<< " " << kgp->Ek << " " << kgp->Sh << " " << kgp->Ra << " " << kgp->Xd1 << " " << kgp->Xq1 << " " << kgp->Xd << " " << kgp->Xq
			<< " " << kgp->Td01 << " " << kgp->Tq01 << " " << kgp->XL << " " << kgp->Tj << " " << kgp->D << " " << kgp->gene_model;
		cout << endl;
	}
	cout << endl;

	//励磁参数
	vector<EXCI_PARA>::iterator kei;
	for (kei = para_exci.begin(); kei != para_exci.end(); kei++)
	{
		cout << kei->exci_model << " " << kei->name << " " << kei->Tr << " " << kei->Ta << " " << kei->Te << " " << kei->Tf << endl;
	}
	cout << endl;

	//PSS参数
	vector<PSS_PARA>::iterator kpi;
	for (kpi = para_pss.begin(); kpi != para_pss.end(); kpi++)
	{
		cout << kpi->pss_model << " " << kpi->name << " " << kpi->Kqs << " " << kpi->Tqs << " "
			<< kpi->Tq << " " << kpi->Tq1 << " " << kpi->Tq11 << " " << kpi->Tq2 << " " << kpi->Tq21 <<endl;
	}
	cout << endl;

	//负荷参数
	vector<LOAD_PARA>::iterator kli;
	for(kli = para_load.begin(); kli != para_load.end(); kli++)
	{
		cout << kli->ap << " " << kli->aq  << " " << kli->bp  << " " << kli->bq  << " " << kli->cp  << " " << kli->cq << endl;
	}
	cout << endl;

	//节点潮流数据
	vector<NODE_DATA>::iterator knb;
	for (knb = data_node.begin(); knb != data_node.end(); knb++)
	{
		cout << knb->name << "\t" << setw(20) << knb->vol_pu << setw(20) << knb->vol_ang << setw(20) << knb->P << setw(20) << knb->Q << endl;
	}
	cout << endl;

	*/

	

	return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门提示: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
