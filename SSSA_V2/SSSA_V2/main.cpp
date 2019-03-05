// SSSA_V2.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include <iostream>
#include"head.h"
#include"head.cpp"
#include"sssa_inputfile.cpp"
#include"sssa_dimension.cpp"
#include"sssa_genesystem.cpp"
#include"sssa_networksystem.cpp"
#include"sssa_connect_system.cpp"


int main()
{
	char name_st[50];
	strcpy_s(name_st, "IEEE90");			//系统名字

	cout << "start" << endl;

	//======================读取系统参数数据=============================
	vector<NODE> node_table;				//存储节点对应表
	vector<NODE_DATA> data_node;			//存储所有节点的潮流结果
	vector<GENE_PARA> para_gene;			//存储各台发电机参数
	vector<LOAD_PARA> para_load;			//存储负荷参数
	sssa_inputfile(name_st, node_table, data_node, para_gene, para_load);

	int num_node, num_gene, nd, ng;	//系统节点数，发电机节点数，2*系统节点数，2*发电机节点数
	num_gene = para_gene.size();
	num_node = node_table.size();
	nd = 2 * num_node;
	ng = 2 * num_gene;
	
	//读取导纳矩阵
	double **Ymatrix = dmatrix(nd + 1, nd + 1);
	for (int i = 1; i <= nd; i++) {
		for (int j = 1; j <= nd; j++)	Ymatrix[i][j] = 0.0;		//初始化
	}
	sssa_inputY(name_st, Ymatrix);
	

	//=======================确定发电机系统模型维数===============================
	vector<int> gdim, gl_gene;		//存储系统维数, gl_exci, gl_pss...
	vector<string> gname;
	sssa_dimension(gdim, gname, gl_gene, para_gene);
	

	//=====================形成发电机系统线性化模型=========================
	int dim_gene_all = gdim[gdim.size() - 1];						//所有发电机的总维数
	double **Tgsys = dmatrix(dim_gene_all + 1, dim_gene_all + 1);	//声明发电机系统T矩阵			
	double **Jgsys = dmatrix(dim_gene_all + 1, dim_gene_all + 1);	//声明发电机系统J矩阵
	for (int i = 1; i <= dim_gene_all; i++) {	//初始化
		for (int j = 1; j <= dim_gene_all; j++)		Tgsys[i][j] = Jgsys[i][j] = 0.0;
	}
	vector<GENE_INIT> init_gene;
	sssa_genesystem(name_st, data_node, para_gene, init_gene, gdim, gl_gene, Tgsys, Jgsys);

	//===========================形成网络方程=================================
	double **Ysystem = dmatrix(nd + 1, nd + 1);			//修正后的网络方程
	for (int i = 1; i <= nd; i++) {
		for (int j = 1; j <= nd; j++)	Ysystem[i][j] = Ymatrix[i][j];		//初始化
	}
	sssa_networksystem(Ysystem, node_table, data_node, para_load);

	
	//==========================连接全系统线性化模型============================
	double **T = dmatrix(dim_gene_all + nd + 1, dim_gene_all + nd + 1);
	double **J = dmatrix(dim_gene_all + nd + 1, dim_gene_all + nd + 1);	//全系统的线性化模型T，J矩阵
	for (int i = 1; i <= dim_gene_all + nd; i++) {
		for (int j = 1; j <= dim_gene_all + nd; j++)	T[i][j] = J[i][j] = 0.0;		//初始化
	}
	sssa_connect_system(name_st, Tgsys, Jgsys, Ysystem, node_table, init_gene, para_gene, gdim, gl_gene, T, J);
	
	cout << "end" << endl;



	/*   测试输出
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

	vector<NODE_DATA>::iterator knb;
	for (knb = data_node.begin(); knb != data_node.end(); knb++)
	{
		cout << knb->no << " " << knb->name << " " << knb->P << " " << knb->Q << endl;
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
