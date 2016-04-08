#if 0
#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<fstream>
#include<time.h>
#include<limits.h>
#include <iomanip>
#include<algorithm>
#include<Windows.h>
#include <sstream>
#include"TraverseDir.h"
#include<map>
#define N 5000
using namespace std;
void analyze_result(string fnr, string fnw, int whe_sol_seq, int whe_of_sol_each_run, int whe_curf_tm)
{	//for current f with respect to iteration, genteration, and time.
	FILE *fpr = fopen(fnr.c_str(), "r");
	if (fpr == NULL)
	{
		cout << fpr << endl;	perror("fpr");
		exit(0);
	}
	ofstream fpw;
	if (whe_curf_tm == 0)
		fpw.open(fnw, ios::app | ios::out);
	else
		fpw.open(fnw, ios::out);

	if (!fpw.is_open())
	{
		cout << " open failed" << endl;	perror("fpw");
		exit(0);
	}
	char strline1[2 * N];
	char strline2[2 * N];
	int a[N + 3][8];//run_cnt,A,B,C,  f,iter,gen,time(s)
	int b[N][5];
	int rc;
	memset(b, 0, sizeof(b));
	cout << fnr << endl;
	if (whe_curf_tm == 0)
	{
		if (whe_of_sol_each_run == 1)
			fpw << fnr << endl;
		else
			fpw << fnr << "\t";
	}
	int line_cnt = 0;
	while (fgets(strline1, 2 * N, fpr))
	{
		line_cnt += 1;
		if (whe_sol_seq == 1 && line_cnt % 2 == 0)
			continue;
		//fputs(strline1,stdout);
		sscanf(strline1, "%d", &rc);
		rc -= 1;
		sscanf(strline1, "%d %d %d %d %d %d %d %d", &a[rc][0], &a[rc][1]
			, &a[rc][2], &a[rc][3], &a[rc][4], &a[rc][5], &a[rc][6], &a[rc][7]);
		b[a[rc][4]][0] += 1;
		b[a[rc][4]][1] += a[rc][5];
		b[a[rc][4]][2] += a[rc][6];
		b[a[rc][4]][3] += a[rc][7];
	}
	if (whe_curf_tm)
	{
		for (int i = N - 1; i >= 1; i--)
		{
			if (b[i - 1][3] == b[i][3])//time
				continue;
			cout << i << "\t";//obj
			fpw << i << "\t";
			for (int j = 0; j<4; j++)
			{
				cout << b[i][j] << "\t";//cnt, iter, gen,time
				if (j == 0)
					fpw << b[i][j] << "\t";
				else
					fpw << setiosflags(ios::fixed) << setprecision(2) << b[i][j] / (float)b[i][0] << "\t";
			}
			cout << endl;
			fpw << endl;
		}
	}
	cout << endl;
	if (whe_curf_tm)
		return;
	rc += 1;
	for (int i = 0; i<rc; i++)
	{
		if (i == 0)
		{
			for (int j = 0; j<8; j++)
			{
				cout << a[i][j] << "\t";
				if (whe_of_sol_each_run)
					fpw << a[i][j] << "\t";
				a[N][j] = a[N + 1][j] = a[N + 2][j] = a[i][j];
			}
		}
		else
		{
			for (int j = 0; j<8; j++)
			{
				cout << a[i][j] << "\t";
				if (whe_of_sol_each_run)
					fpw << a[i][j] << "\t";
				if (a[N][j]>a[i][j])
					a[N][j] = a[i][j];
				if (a[N + 1][j]<a[i][j])
					a[N + 1][j] = a[i][j];
				a[N + 2][j] += a[i][j];
			}
		}
		cout << endl;
		if (whe_of_sol_each_run)
			fpw << endl;
	}
	cout << endl;
	if (whe_of_sol_each_run)
		fpw << endl;
	int min_f_cnt = 0, min_sum_tm = 0;
	for (int i = 0; i<rc; i++)
	{
		if (a[i][4] == a[N][4])
		{
			min_f_cnt += 1;
			min_sum_tm += a[i][7];
		}
	}
	for (int i = N; i<N + 3; i++)
	{
		if (i == N&&whe_of_sol_each_run == 0)
			fpw << "MIN\t";
		if (i == N + 1 && whe_of_sol_each_run == 0)
			fpw << "MAX\t";
		if (i == N + 2 && whe_of_sol_each_run == 0)
			fpw << "AVG\t";
		if (i != N + 2)//not the last row(the average)
		{
			for (int j = 0; j<8; j++)
			{
				cout << a[i][j] << "\t";
				fpw << a[i][j] << "\t";
			}
		}
		else
		{
			for (int j = 0; j<8; j++)//the average
			{
				printf("%.2f\t", a[i][j] / (float)rc);
				fpw << setiosflags(ios::fixed) << setprecision(2) << a[i][j] / (float)rc << "\t";
			}
			printf("%d\t%d\t%.2f\t", min_f_cnt, rc, min_sum_tm / (float)min_f_cnt);
			fpw << min_f_cnt << "\t" << rc << "\t" << a[N][4] << "\t"
				<< setiosflags(ios::fixed) << setprecision(2) << min_sum_tm / (float)min_f_cnt << "\t";
		}
		cout << endl;
		if (whe_of_sol_each_run)
			fpw << endl;
	}
	fpw << endl;

}
void analyze_result1(char *fnr, char *fnw, int whe_sol_seq, int whe_of_sol_each_run, int whe_curf_tm)
{	//for current f with respect to iteration, genteration, and time.
	FILE *fpr = fopen(fnr, "r");
	if (fpr == NULL)
	{
		cout << fpr << endl;	perror("fpr");
		exit(0);
	}
	ofstream fpw;
	if (whe_curf_tm == 0)
		fpw.open(fnw, ios::out);
	else
		fpw.open(fnw, ios::out);

	if (!fpw.is_open())
	{
		cout << " open failed" << endl;	perror("fpw");
		exit(0);
	}
	char strline1[2 * N];
	char strline2[2 * N];
	int a[N + 3][8];//run_cnt,A,B,C,  f,iter,gen,time(s)
	int b[N][5];
	int rc;
	memset(b, 0, sizeof(b));
	cout << fnr << endl;
	if (whe_curf_tm == 0)
	{
		if (whe_of_sol_each_run == 1)
			fpw << fnr << endl;
		else
			fpw << fnr << "\t";
	}
	int line_cnt = 0;
	while (fgets(strline1, 2 * N, fpr))
	{
		line_cnt += 1;
		if (whe_sol_seq == 1 && line_cnt % 2 == 0)
			continue;
		//fputs(strline1,stdout);
		sscanf(strline1, "%d", &rc);
		rc -= 1;
		sscanf(strline1, "%d %d %d %d %d %d %d %d", &a[rc][0], &a[rc][1]
			, &a[rc][2], &a[rc][3], &a[rc][4], &a[rc][5], &a[rc][6], &a[rc][7]);
		b[a[rc][4]][0] += 1;
		b[a[rc][4]][1] += a[rc][5];
		b[a[rc][4]][2] += a[rc][6];
		b[a[rc][4]][3] += a[rc][7];
		if (a[rc][4] == whe_curf_tm)
			fpw << strline1;
	}
}
LPCWSTR stringToLPCWSTR(std::string orig)
{
	size_t origsize = orig.length() + 1;
	const size_t newsize = 100;
	size_t convertedChars = 0;
	wchar_t *wcstring = (wchar_t *)malloc(sizeof(wchar_t)*(orig.length() - 1));
	mbstowcs_s(&convertedChars, wcstring, origsize, orig.c_str(), _TRUNCATE);

	return wcstring;
}
void add_best_found_solution_to_instance(string fnw,string fnr,string fn)
{
	fnr = fnr + "R" + fn; 
	ifstream ifs(fnr);
	if (!ifs.is_open())
	{
		cout << fnr << endl; perror("file_input.");
		exit(0);
	}
	fnw = fnw + fn;
	ofstream ofs(fnw, ios::app | ios::out);
	if (!ofs.is_open())
	{
		cout << fnw << endl; perror("file_output.");
		exit(0);
	}
	string strline;
	ifs >> strline;
	ifs >> strline;
	ifs >> strline;
	while (ifs.good())
	{
		ifs >> strline;
		if(strline!="")
		ofs << strline << endl;
		//cout << strline << endl;
		strline = "";
	}
	ifs.close();
	ofs.close();
}
void analyze_folder(string dir, string fw)
{
	//构造类对象
	CStatDir statdir;

	//设置要遍历的目录
	if (!statdir.SetInitDir(dir.c_str()))
	{
		puts("目录不存在。");
		return;
	}

	//开始遍历
	statdir.BeginBrowse("*.*");
		
	sort(statdir.vec_all_file.begin(), statdir.vec_all_file.end());
	int pic_cnt = 0;
	for (vector<string>::iterator iter = statdir.vec_all_file.begin();
	iter != statdir.vec_all_file.end(); iter++)
	{
		if ((*iter).find("barabasi_albert") != string::npos)//erdos_renyi,barabasi_albert
		{
			analyze_result(const_cast<char*>((*iter).data()), fw, 0, 0, 0);
			cout << *iter << endl;
		}
		add_best_found_solution_to_instance(dir, fw, (*iter).data());
	}
	printf("文件总数: %d\n子目录总数:%d\n", statdir.GetFileCount(), statdir.GetSubdirCount());
}

void friedman_test(string fnr, string fnw, int whe_sol_seq, int whe_of_sol_each_run, int run_time)
{	//for current f with respect to iteration, genteration, and time.
	FILE *fpr = fopen(fnr.c_str(), "r");
	if (fpr == NULL)
	{
		cout << fpr << endl;	perror("fpr");
		exit(0);
	}
	ofstream fpw;
	fpw.open(fnw, ios::out | ios::app);

	if (!fpw.is_open())
	{
		cout << " open failed" << endl;	perror("fpw");
		exit(0);
	}
	char strline1[2 * N];
	char strline2[2 * N];
	int a[N + 3][8];//run_cnt,A,B,C,  f,iter,gen,time(s)
	int b[N][5];
	int rc;
	memset(b, 0, sizeof(b));
	cout << fnr << endl;
	/*if (whe_of_sol_each_run == 1)
	fpw << fnr << endl;
	else
	fpw << fnr << "\t";*/
	int line_cnt = 0, run_cnt = 0;
	int flags[20] = { 0 };
	int temp[2][8];
	int min_f = 0;
	for (int i = 0; i < 8; i++)
		temp[1][i] = 0;
	while (fgets(strline1, 2 * N, fpr))
	{
		line_cnt += 1;
		if (whe_sol_seq == 1 && line_cnt % 2 == 0)
			continue;
		//fputs(strline1,stdout);
		sscanf(strline1, "%d", &rc);
		rc -= 1;
		if (rc>0 && flags[rc - 1] == 0)
		{
			if (whe_of_sol_each_run == 1)
				fpw << strline2;
			cout << strline2;
			flags[rc - 1] = 1;
			sscanf(strline2, "%d %d %d %d %d %d %d %d", &temp[0][0], &temp[0][1]
				, &temp[0][2], &temp[0][3], &temp[0][4], &temp[0][5], &temp[0][6], &temp[0][7]);
			temp[1][4] += temp[0][4];
			run_cnt += 1;
			if (run_cnt == 1)
				min_f = temp[0][4];
			else
			{
				if (min_f>temp[0][4])
					min_f = temp[0][4];
			}
		}
		sscanf(strline1, "%d %d %d %d %d %d %d %d", &a[rc][0], &a[rc][1]
			, &a[rc][2], &a[rc][3], &a[rc][4], &a[rc][5], &a[rc][6], &a[rc][7]);
		b[a[rc][4]][0] += 1;
		b[a[rc][4]][1] += a[rc][5];
		b[a[rc][4]][2] += a[rc][6];
		b[a[rc][4]][3] += a[rc][7];
		if (a[rc][7] > run_time&&flags[rc] == 0)
		{
			if (whe_of_sol_each_run == 1)
				fpw << strline2;
			cout << strline2;
			flags[rc] = 1;
			sscanf(strline2, "%d %d %d %d %d %d %d %d", &temp[0][0], &temp[0][1]
				, &temp[0][2], &temp[0][3], &temp[0][4], &temp[0][5], &temp[0][6], &temp[0][7]);
			temp[1][4] += temp[0][4];
			run_cnt += 1;
			if (run_cnt == 1)
				min_f = temp[0][4];
			else
			{
				if (min_f>temp[0][4])
					min_f = temp[0][4];
			}
		}
		strcpy(strline2, strline1);
	}
	if (flags[rc] == 0)
	{
		if (whe_of_sol_each_run == 1)
			fpw << strline2;
		cout << strline2;
		flags[rc] = 1;
		sscanf(strline2, "%d %d %d %d %d %d %d %d", &temp[0][0], &temp[0][1]
			, &temp[0][2], &temp[0][3], &temp[0][4], &temp[0][5], &temp[0][6], &temp[0][7]);
		temp[1][4] += temp[0][4];
		run_cnt += 1;
		if (run_cnt == 1)
			min_f = temp[0][4];
		else
		{
			if (min_f>temp[0][4])
				min_f = temp[0][4];
		}
	}
	fpw << fnr << "\t"
		<< run_time << "\t"
		<< min_f << "\t"
		<< setiosflags(ios::fixed) << setprecision(2)
		<< temp[1][4] / (float)run_cnt << endl;
}
int main(int argc, char **argv)
{
	char *rgv[]={"",//0
		"whe_if_sol_seq","0",//1,2
		"whe_of_sol_each_run","1",//3,4
		"_IFP","D:\\MyProjects\\PMS1\\PMS1\\instance\\BB_Problem_BestSolution\\",//5,6	
		"_OFP","D:\\MyProjects\\PMS1\\PMS1\\instance\\BB_BestSolution\\",//7,8
		"_OF","friedman_test_QD-HA_tm60",//9,10
		"_IF","G23_r1_r20_s3600_sc2_t10000_p13_alpha700_beta250_ws0_xa500_xb300_xc100_qa100_qb5000_t1a100_t1b200_t2a200_t2b200",//11,12
		"_whe_curf_tm","60"//13,14
	};
	argv=rgv;
	std::map<string, string> argv_map;
	for (int i = 1; i < sizeof(rgv) / sizeof(rgv[0]); i += 2)
		argv_map[string(argv[i])] = string(argv[i + 1]);
	string fr = argv_map.at("_IFP") + argv_map.at("_IF") + ".txt";
	string fw = argv_map.at("_OFP") + argv_map.at("_OF") + ".txt";
	/*strcat(fw,argv[12]);
	strcat(fw,"_f");
	strcat(fw,argv[14]);
	strcat(fw,"_analyze_boxf");*/

	analyze_folder(argv_map.at("_IFP"), argv_map.at("_OFP"));
	//analyze_result(fr,fw, atoi(argv[2]),atoi(argv[4]),atoi(argv[14]));	
	//friedman_test(fr, fw, atoi(argv[2]), atoi(argv[4]), atoi(argv[14]));
	cout << fr << " is analyzed, friedman test (min_f,avg_f at a specific time) the results are saved in:" << endl << fw;
}
#endif