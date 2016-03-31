#if 1
#include<iostream>
#include<fstream>
#include<string>
#include<time.h>
#include<iomanip>
#include<algorithm>
using namespace std;
typedef double proc_type;	// type of processing time
typedef double dete_type;	// type of deterioration effect
typedef double obj_type;	// type of objective value
typedef double perf_type;	// type of performance level

class PMS
{
	const int m, n;	// m and n are the number of machines and jobs
	const proc_type **p;	// processing time
	const dete_type **d;	// deterioration effect
	const int **s_opt;	// the optimal solution
	const obj_type sol_obj_opt;	// the optimal objective value
	const string file_input, file_output;
	const int sol_num;	// number of solution

	int ***s;	// solution
	obj_type *sol_obj;	// objective value
	
	clock_t start_tm, end_tm;
	class cmpSort;
public:
	PMS(int _m, int _n, proc_type ** _p, dete_type **_d,
		int **_s_opt, obj_type _sol_obj_opt, string _file_input,
		string _file_output, int _sol_num);
	void display_problem();
	void display_solution(int);
	void check_solution(int);
	void get_min_p(proc_type *);
	void init_solution(int);
	proc_type calculate_completion_time(int, int);
	void remove_job(int, int, int);
	void test();
};
void load_data(string file_input,string file_output,int sol_num, PMS *pms)
{
	int m, n;	
	
	ifstream ifs(file_input);
	if (!ifs.is_open())
	{
		cout << file_input << endl; perror("file_input.");
		exit(0);
	}
	ifs >> n >> m;
	proc_type **p = new proc_type *[m];
	dete_type **d = new dete_type *[m];

	obj_type sol_obj_opt;	// objective value
	int **s_opt = new int *[m + 1];
	for (int j = 0; j <= m; j++)
		s_opt[j] = new int[n + 1];

	for (int i = 0; i <= m; i++)
	{
		p[i] = new proc_type[n + 1];
		d[i] = new dete_type[n + 1];
	}
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			if (!ifs.good())
			{
				perror("input processing time error.");
				exit(0);
			}
			ifs >> p[i][j];
		}
	}
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			if (!ifs.good())
			{
				perror("input deterioration effect error.");
				exit(0);
			}
			ifs >> d[i][j];
			d[i][j] *= 0.01;
		}
	}
	for (int i = 1; i <= m; i++)
	{
		ifs >> s_opt[i][0];	// the optimal number of jobs to be assigned to machine i
		for (int j = 1; j <= s_opt[i][0]; j++)
		{
			ifs >> s_opt[i][j];	// the list of jobs to be assinged to machine i
		}
	}
	ifs >> sol_obj_opt;
	pms = new PMS(m, n, p, d, s_opt, sol_obj_opt, file_input, file_output, sol_num);
	//pms = new PMS(m, n, p, d, s_opt, sol_obj_opt, file_input, file_output, sol_num);
	ifs.close();
}

PMS::PMS( int _m, int _n, proc_type ** _p, dete_type **_d,
	int **_s_opt, obj_type _sol_obj_opt, string _file_input,
	string _file_output, int _sol_num) :m(_m), n(_n),
	p(_p), d(_d), s_opt(_s_opt), sol_obj_opt(_sol_obj_opt), file_input(_file_input),
	file_output(_file_output), sol_num(_sol_num) 
{
	s = new int **[sol_num];
	sol_obj = new obj_type[sol_num];
	for (int i = 0; i < sol_num; i++)
	{
		s[i] = new int *[m + 1];
		for (int j = 0; j <= m; j++)
			s[i][j] = new int[n + 1];
	}
	// copy the optimal solution to s[0]
	sol_obj[0] = sol_obj_opt;
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
			s[0][i][j] = s_opt[i][j];
	}
}

class PMS::cmpSort
{
public:
	cmpSort(double *as) :array_sort(as) {}
	bool operator()(const int &l, const int &r)
	{
		return array_sort[l] > array_sort[r];
	}
private:
	const double *array_sort;
};
void PMS::display_problem()
{
	cout << file_input << ":\nn: " << n << ", m: " << m << endl;
	cout << "processing time:" << endl;
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
			cout << p[i][j] << " ";
		cout << endl;
	}
	cout << "deterioration effect: " << endl;
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
			cout << d[i][j] << " ";
		cout << endl;
	}
}
void PMS::display_solution(int sol_index)
{
	cout << "*** sol_index:" << sol_index
		<< ", obj: " /*<< fixed << setprecision(4)*/
		<< sol_obj[sol_index] << " ***" << endl;
	for (int i = 1; i <= m; i++)
	{
		cout << s[sol_index][i][0] << endl;
		for (int j = 1; j <= s[sol_index][i][0]; j++)
			cout << s[sol_index][i][j] << " ";
		cout << endl;
	}
}
void PMS::check_solution(int sol_index)
{
	perf_type **q = new perf_type *[m+1];
	obj_type *c = new obj_type [m+1];
	memset(c, 0, (m + 1)*sizeof(obj_type));
	for (int i = 1; i <= m; i++)
	{
		q[i] = new perf_type[n + 1];
		q[i][0] = 0;
		q[i][1] = 1;
	}
	for (int i = 1; i <= m; i++)
	{
		c[i] += p[i][s[sol_index][i][1]];
		for (int j = 2; j <= s[sol_index][i][0]; j++)
		{
			q[i][j] = (1 - d[i][s[sol_index][i][j - 1]])*q[i][j - 1];
			c[i] += p[i][s[sol_index][i][j]] / q[i][j];
		}
		cout << "completion time of machine " << i << ": " <</*fixed<<setprecision(4)<<*/ c[i] << endl;
	}
}
void PMS::get_min_p(proc_type *min_p)
{
	for (int i = 1; i <= n; i++)
	{
		min_p[i] = p[1][i];
		for (int j = 2; j <= m; j++)
		{
			if (p[j][i]<min_p[i])
				min_p[i] = p[j][i];
		}
	}
}
proc_type PMS::calculate_completion_time(int sol_index, int mach_index)
{
	proc_type ret = 0;
	perf_type *q = new perf_type[n + 1];
	q[1] = 1;
	ret += p[mach_index][s[sol_index][mach_index][1]];
	for (int i = 2; i <= s[sol_index][mach_index][0]; i++)
	{
		q[i] = (1 - d[mach_index][s[sol_index][mach_index][i]]) / q[i - 1];
		ret += p[mach_index][s[sol_index][mach_index][i]] / q[i];
	}
	return ret;
}
void PMS::remove_job(int sol_index, int mach_index, int key)
{
	for (int i = 1; i <= s[sol_index][mach_index][0]; i++)
	{
		if (s[sol_index][mach_index][i] == key)
		{
			for (int j = i + 1; j <= s[sol_index][mach_index][0]; j++)
				s[sol_index][mach_index][j - 1] = s[sol_index][mach_index][j];
			break;
		}
	}
}
void PMS::init_solution(int sol_index)
{
	test();
	proc_type *min_p = new proc_type[n + 1];
	get_min_p(min_p);
	for (int i = 1; i <= n; i++)
	{
		cout << min_p[i] << "\t";
	}
	cout << endl;
	int *init_s = new int[n + 1];
	for (int i = 1; i <= n; i++)
		init_s[i] = i;
	sort(init_s + 1, init_s + n + 1, cmpSort(min_p));
	for (int i = 1; i <= n; i++)
	{
		cout << init_s[i] << " " << min_p[init_s[i]] << "\t";
	}
	cout << endl;
	int sol_index_try = 2;
	for (int i = 1; i <= m; i++)
	{
		s[sol_index][i][0] = 0;	// init the number of jobs assigned to machine i
		s[sol_index_try][i][0] = 0;
	}
	proc_type **r = new proc_type *[m + 1];
	for (int i = 1; i <= m; i++)
		r[i] = new proc_type[n + 1];
	for (int i = 1; i <= n; i++)
	{
		cout << "job " << i << ": " << init_s[i] << endl;
		proc_type min_completion_time = DBL_MAX;
		int min_ct_mach_index;
		for (int j = 1; j <= m; j++)
		{
			cout << "machine " << j << ": ";
			s[sol_index][j][0] += 1;
			s[sol_index][j][s[sol_index][j][0]] = init_s[i];
			r[j][s[sol_index][j][0]] = p[sol_index][init_s[i]] * (1 - d[sol_index][init_s[i]]) /
				d[sol_index][init_s[i]];
			for (int k = 1; k <= s[sol_index][j][0]; k++)
				cout << s[sol_index][j][k] << "," << r[j][k] << " ";
			sort(s[sol_index][j] + 1, s[sol_index][j] + s[sol_index][j][0] + 1, cmpSort(r[j]));
			cout << "|";
			for (int k = 1; k <= s[sol_index][j][0]; k++)
				cout << s[sol_index][j][k] << "," << r[j][k] << " ";
			proc_type temp_c = calculate_completion_time(sol_index, j);
			cout << ", ct: " << temp_c << endl;
			if (temp_c < min_completion_time)
			{
				min_completion_time = temp_c;
				min_ct_mach_index = j; 
				memcpy(s[sol_index_try][min_ct_mach_index], s[sol_index][min_ct_mach_index],
					(s[sol_index][min_ct_mach_index][0] + 1)*sizeof(int));				
				/*cout << "min save: ";
				for (int k = 1; k <= s[sol_index][min_ct_mach_index][0]; k++)
					cout << s[sol_index_try][min_ct_mach_index][k] << " ";
				cout << endl;*/
			}
			remove_job(sol_index, j, init_s[i]);
			s[sol_index][j][0] -= 1;
		}
		cout << "min_ct: " << min_completion_time << " min_ct_mach_index: " << min_ct_mach_index << endl;
		memcpy(s[sol_index][min_ct_mach_index], s[sol_index_try][min_ct_mach_index],
			(s[sol_index_try][min_ct_mach_index][0] + 1)*sizeof(int));
		sort(r[min_ct_mach_index] + 1, r[min_ct_mach_index] + s[sol_index_try][min_ct_mach_index][0]);
	}
}
void PMS::test()
{
	cout << endl << "This is test function." << endl;
	int cnt = 3;
	double si[] = { 0,1 ,2};
	double rd[] = { 7000,2000,4000 };

	sort(si, si + cnt, cmpSort(rd));
	for (int i = 0; i < cnt; i++)
		cout << si[i] << "," << rd[i] << " ";
}

int main(int argc, char **argv)
{
	cout << "This is the PMS problem solver." << endl;
	PMS *pms;
	load_data("instance\\OB_Problem_OptSolution\\Ni_8_2-1_1_1", "", 5, pms);
	//PMS pms("instance\\OB_Problem_OptSolution\\Ni_8_2-1_1_1", "",3);
	pms->display_problem();
	pms->display_solution(0);
	pms->check_solution(0);
	pms->init_solution(1);
	pms->display_solution(1);
	system("pause");
}
#endif
