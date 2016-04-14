#if 0
#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<fstream>
#include<string>
#include<time.h>
#include<iomanip>
#include<algorithm>
#include<math.h>
#include<map>
#include<vector>
using namespace std;
typedef double proc_type;	// type of processing time
typedef double dete_type;	// type of deterioration effect
typedef double obj_type;	// type of objective value
typedef double perf_type;	// type of performance level

class PMS
{
private:
	int m, n;	// m and n are the number of machines and jobs
	proc_type **p;	// processing time
	dete_type **d;	// deterioration effect
	int **s_opt;	// the optimal solution
	obj_type sol_obj_opt;	// the optimal objective value
	string file_input, file_output;
	int sol_num;	// number of solution

	int ***s;	// solution
	obj_type **c, obj_given;	// completion time
	int *mm;	// makespan machine

	perf_type **r;	// exact sort accordance
	perf_type **charact;	// characteristics sort accordance

	clock_t start_tm, end_tm;
	class cmpSort;
	const double MIN_EQUAL = 0.001;
	int rand_seed;
	ofstream ofs;
public:
	enum R_Mode { EXAC = 0, MINP, MAXP, MIND, MAXD, MINPDD, MAXPDD, MINPD, MAXPD, SIZE };
	enum Cmp_Mode { GREATER, LESS };
	enum NS_Mode { SWAP, INSERT };
	int whe_save_sol_seq;
	string ins_name;
	double control_para, temperature;
	PMS(string, string, int);
	~PMS();
	PMS(int _m, int _n, proc_type ** _p, dete_type **_d,
		int **_s_opt, obj_type _sol_obj_opt, string _file_input,
		string _file_output, int _sol_num);
	void display_problem();
	void display_solution(int);
	void check_solution(int);
	void init_solution(int, int, R_Mode);
	void calculate_completion_time(int, int);
	void calculate_obj(int);
	void remove_job(int, int, int);
	void local_search(int, int, NS_Mode);
	void local_search_hybrid(int, int, NS_Mode);
	void local_search_hybrid1(int, int, NS_Mode);
	void iterated_local_search(int, int, R_Mode, NS_Mode, int);
	void replace_solution(int, int);
	void perturb(int, int);
	void swap(int, int, int, int, int);
	void insert(int, int, int, int);
	void save_solution(int, int, int, int);
	void test();
};

class PMS::cmpSort
{
public:
	cmpSort(double *as, Cmp_Mode _cm = Cmp_Mode::GREATER) :
		array_sort(as), cm(_cm) {}
	bool operator()(int &l, int &r)
	{
		return (cm == Cmp_Mode::GREATER) ? (array_sort[l] > array_sort[r]) :
			(array_sort[l] < array_sort[r]);
	}
private:
	double *array_sort;
	Cmp_Mode cm;
};
PMS::PMS(string file_input, string file_output, int sol_num)
{
	ifstream ifs(file_input);
	if (!ifs.is_open())
	{
		cout << file_input << endl; perror("file_input.");
		exit(0);
	}
	ofs.open(file_output, ios::app | ios::out);
	if (!ofs.is_open())
	{
		cout << file_output << endl; perror("file_output.");
		exit(0);
	}
	ofs.setf(ios::fixed, ios::floatfield);  
	ofs.precision(6);  
	ifs >> n >> m;
	p = new proc_type *[m + 1];
	d = new dete_type *[m + 1];
	r = new perf_type*[m + 1];
	for (int i = 0; i <= m; i++)
	{
		p[i] = new proc_type[n + 1];
		d[i] = new dete_type[n + 1];
		r[i] = new perf_type[n + 1];
	}
	charact = new perf_type*[R_Mode::SIZE];
	for (int i = 1; i < R_Mode::SIZE; i++)
		charact[i] = new perf_type[n + 1];
	s = new int **[sol_num];
	c = new obj_type *[sol_num];
	mm = new int[sol_num];
	for (int i = 0; i < sol_num; i++)
	{
		s[i] = new int *[m + 1];
		c[i] = new obj_type[m + 1];
		for (int j = 0; j <= m; j++)
		{
			s[i][j] = new int[n + 1];
		}
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
	// calculate r
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			r[i][j] = p[i][j] * (1 - d[i][j]) / d[i][j];
		}
	}
	// calculate charact
	for (int rm = 1; rm < R_Mode::SIZE; rm++)
	{
		for (int j = 1; j <= n; j++)
		{
			if (rm == R_Mode::MINP || rm == R_Mode::MAXP)
				charact[rm][j] = p[1][j];
			else if (rm == R_Mode::MIND || rm == R_Mode::MAXD)
				charact[rm][j] = d[1][j];
			else if (rm == R_Mode::MINPDD || rm == R_Mode::MAXPDD)
				charact[rm][j] = p[1][j] * (1 - d[1][j]) / d[1][j];
			else if (rm == R_Mode::MINPD | rm == R_Mode::MAXPD)
				charact[rm][j] = p[1][j] * (1 - d[1][j]);
			for (int i = 2; i <= m; i++)
			{
				if (rm == R_Mode::MINP&&p[i][j] < charact[rm][j])
					charact[rm][j] = p[i][j];
				if (rm == R_Mode::MAXP&&p[i][j] > charact[rm][j])
					charact[rm][j] = p[i][j];
				if (rm == R_Mode::MIND&&d[i][j] < charact[rm][j])
					charact[rm][j] = d[i][j];
				if (rm == R_Mode::MAXD&&d[i][j] > charact[rm][j])
					charact[rm][j] = d[i][j];
				if (rm == R_Mode::MINPDD && (p[i][j] * (1 - d[i][j]) / d[i][j]) < charact[rm][j])
					charact[rm][j] = p[i][j] * (1 - d[i][j]) / d[i][j];
				if (rm == R_Mode::MAXPDD && (p[i][j] * (1 - d[i][j]) / d[i][j]) > charact[rm][j])
					charact[rm][j] = p[i][j] * (1 - d[i][j]) / d[i][j];
				if (rm == R_Mode::MINPD && (p[i][j] * (1 - d[i][j])) < charact[rm][j])
					charact[rm][j] = p[i][j] * (1 - d[i][j]);
				if (rm == R_Mode::MAXPD && (p[i][j] * (1 - d[i][j])) > charact[rm][j])
					charact[rm][j] = p[i][j] * (1 - d[i][j]);
			}
		}
	}
	int sol_index_opt = 0;
	for (int i = 1; i <= m; i++)
	{
		ifs >> s[sol_index_opt][i][0];	// the optimal number of jobs to be assigned to machine i
		for (int j = 1; j <= s[sol_index_opt][i][0]; j++)
		{
			ifs >> s[sol_index_opt][i][j];	// the list of jobs to be assinged to machine i
		}
	}
	ifs >> obj_given;
	for (int i = 1; i <= m; i++)
	{
		sort(s[sol_index_opt][i] + 1, s[sol_index_opt][i] + s[sol_index_opt][i][0] + 1, cmpSort(r[i]));
		calculate_completion_time(sol_index_opt, i);
	}
	calculate_obj(sol_index_opt);
	if (abs(obj_given - c[sol_index_opt][0]) > MIN_EQUAL)
	{
		cout << "the given optimal solution is wrong." 
			<< obj_given<<" "<< c[sol_index_opt][0]<<endl;
		//system("pause");
	}
	ifs.close();
}
PMS::~PMS()
{
	for (int i = 0; i < sol_num; i++)
	{
		for (int j = 0; j <= m; j++)
		{
			delete[]s[i][j];
		}
		delete[]s[i];
		delete[]c[i];
	}
	delete[]s;
	delete[]c;
	delete[]mm;

	for (int i = 1; i < R_Mode::SIZE; i++)
		delete[]charact[i];
	delete[]charact;

	for (int i = 0; i <= m; i++)
	{
		delete[]p[i];
		delete[]d[i];
		delete[]r[i];
	}
	delete[]p;
	delete[]d;
	delete[]r;
	ofs.close();
	//cout << "destruct function." << endl;
}
void PMS::save_solution(int sol_index, int sol_index_opt, int iterration, int run_cnt)
{
	//end_tm = clock();
	int result_improve,given_result_improve;
	if (abs(c[sol_index][0] - c[sol_index_opt][0]) <= MIN_EQUAL)
		result_improve = 1;	// equal
	else if (c[sol_index_opt][0] - c[sol_index][0] > MIN_EQUAL)
		result_improve = 2;	// improved
	else
		result_improve = 0;	// not improved
	if (abs(c[sol_index][0] - obj_given) <= MIN_EQUAL)
		given_result_improve = 1;
	else if (obj_given - c[sol_index][0] > MIN_EQUAL)
		given_result_improve = 2;
	else
		given_result_improve = 0;
	ofs << run_cnt << "\t"
		//<< c[sol_index_opt][0] << "\t"
		<< c[sol_index][0] << "\t"
		<< mm[sol_index] << "\t"
		<< iterration << "\t"
		<< (end_tm - start_tm) /*/ CLOCKS_PER_SEC*/ << "\t"
		<<given_result_improve<<"\t"
		<< result_improve
		<< endl;
	//cout << run_cnt << "\t"
	//	//<< c[sol_index_opt][0] << "\t"
	//	<< c[sol_index][0] << "\t"
	//	<< mm[sol_index] << "\t"
	//	<< iterration << "\t"
	//	<< (end_tm - start_tm) /*/ CLOCKS_PER_SEC*/ << "\t"
	//	<< given_result_improve << "\t"
	//	<< result_improve
	//	<< endl;
	if (whe_save_sol_seq)	// whe_save_sol_seq
	{
		for (int i = 1; i <= m; i++)
		{
			ofs << c[sol_index][i] << "\t"
				<< s[sol_index][i][0] << "\t";
			for (int j = 1; j <= s[sol_index][i][0]; j++)
				ofs << s[sol_index][i][j] << "\t";
			ofs << endl;
		}
	}
}
void PMS::display_problem()
{
	cout << file_input << "n: " << n << ", m: " << m
		<< ", opt: " << c[0][0] << endl;
	cout << "processing time:" << endl;
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
			cout << p[i][j] << "\t";
		cout << endl;
	}
	cout << "deterioration effect: " << endl;
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
			cout << d[i][j] << "\t";
		cout << endl;
	}
}
void PMS::display_solution(int sol_index)
{
	cout << "*** display sol_index:" << sol_index
		<< ", obj: " /*<< fixed << setprecision(4)*/
		<< c[sol_index][0] << ", " << mm[sol_index] << " ***" << endl;
	for (int i = 1; i <= m; i++)
	{
		obj_type ir = 0;
		cout << s[sol_index][i][0] << ", " << c[sol_index][i] << ": ";
		for (int j = 1; j <= s[sol_index][i][0]; j++)
		{
			cout << s[sol_index][i][j] /*<< ", " << r[i][s[sol_index][i][j]] */<< "\t";
			ir += (j*r[i][s[sol_index][i][j]]);
		}
		cout /*<<ir*/<< endl;
	}
}
void PMS::check_solution(int sol_index)
{
	cout << "*** check sol_index:" << sol_index
		<< ", obj: " << c[sol_index][0] << " ***" << endl;
	perf_type ql;
	obj_type cl, max_cl = 0;
	int max_cl_mach_index, sum_job_index = 0;
	for (int i = 1; i <= m; i++)
	{
		ql = 1;
		cl = 0;
		//cout << cl << " ";
		for (int j = 1; j <= s[sol_index][i][0]; j++)
		{
			cl += p[i][s[sol_index][i][j]] / ql;
			ql *= (1 - d[i][s[sol_index][i][j]]);
			sum_job_index += s[sol_index][i][j];
			//cout << cl << " ";
			if (j>1&&r[i][s[sol_index][i][j - 1]] - r[i][s[sol_index][i][j]] < MIN_EQUAL)
			{
				cout << "ERROR, sol_index: " << sol_index << " not sort by r"
					<< r[i][s[sol_index][i][j - 1]] << " " << r[i][s[sol_index][i][j]] << endl;
				display_solution(sol_index);
				//system("pause");
			}
		}
		if (max_cl < cl)
		{
			max_cl = cl;
			max_cl_mach_index = i;
		}
		//cout << "\ncompletion time of machine " << i << ": " <</*fixed<<setprecision(4)<<*/ cl << endl;
	}
	if (abs(max_cl - c[sol_index][0])>MIN_EQUAL/*||max_cl_mach_index!=mm[sol_index]*/)
	{
		cout << "ERROR, sol_index: " << sol_index
			<< ", real obj: " << max_cl << " " << c[sol_index][0] << endl;
		system("pause");
	}
	if (sum_job_index != (1 + n)*n / 2)
	{
		cout << "ERROR, sol_index: " << sol_index
			<< ", multiple job index" << endl;
		display_solution(sol_index);
		system("pause");
	}
}
void PMS::calculate_completion_time(int sol_index, int mach_index)
{
	/*if (s[sol_index][mach_index][0] == 0)
	{
		c[sol_index][mach_index] = 0;
		return;
	}
	perf_type q = 1;
	c[sol_index][mach_index] = p[mach_index][s[sol_index][mach_index][1]];
	for (int i = 2; i <= s[sol_index][mach_index][0]; i++)
	{
		q *= (1 - d[mach_index][s[sol_index][mach_index][i - 1]]);
		c[sol_index][mach_index] += p[mach_index][s[sol_index][mach_index][i]] / q;
	}*/
	perf_type q = 1;
	c[sol_index][mach_index] = 0;
	for (int i = 1; i <= s[sol_index][mach_index][0]; i++)
	{
		c[sol_index][mach_index] += p[mach_index][s[sol_index][mach_index][i]] / q;
		q *= (1 - d[mach_index][s[sol_index][mach_index][i]]);
	}
}
void PMS::calculate_obj(int sol_index)
{
	c[sol_index][0] = 0;
	for (int i = 1; i <= m; i++)
	{
		if ((c[sol_index][i] - c[sol_index][0]) > MIN_EQUAL)
			c[sol_index][0] = c[sol_index][i];
	}
	for (int i = 1; i <= m; i++)
	{
		if (abs(c[sol_index][i] - c[sol_index][0]) <= MIN_EQUAL)
			mm[sol_index] = i;
	}
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
	s[sol_index][mach_index][0] -= 1;
}
void PMS::init_solution(int sol_index, int sol_index_local, R_Mode rm)
{
	/*for (int i = 1; i <= n; i++)
	{
		cout << charact[rm][i] << "\t";
	}
	cout << endl;*/
	for (int i = 1; i <= n; i++)
		s[sol_index_local][0][i] = i;	// the machine 0 as the temp memory
	sort(s[sol_index_local][0] + 1, s[sol_index_local][0] + n + 1, cmpSort(charact[rm], Cmp_Mode::GREATER));
	/*for (int i = 1; i <= n; i++)
	{
		cout << s[sol_index_local][0][i] << " " << charact[rm][s[sol_index_local][0][i]] << "\t";
	}
	cout << endl;*/
	for (int i = 1; i <= m; i++)
	{
		s[sol_index][i][0] = 0;	// init the number of jobs assigned to machine i
		s[sol_index_local][i][0] = 0;
	}
	for (int j = 1; j <= n; j++)
	{
		/*cout << "job " << j << ": " << s[sol_index_local][0][j] << ", ";
		for (int i = 1; i <= m; i++)
			cout << r[i][s[sol_index_local][0][j]] << "\t";
		cout << endl;*/
		proc_type min_completion_time = DBL_MAX;
		int min_ct_mach_index;
		for (int i = 1; i <= m; i++)
		{
			//cout << "machine " << i << ": ";
			s[sol_index][i][0] += 1;
			s[sol_index][i][s[sol_index][i][0]] = s[sol_index_local][0][j];
			/*for (int k = 1; k <= s[sol_index][i][0]; k++)
				cout << s[sol_index][i][k] << "," << r[i][s[sol_index][i][k]] << " ";*/
			sort(s[sol_index][i] + 1, s[sol_index][i] + s[sol_index][i][0] + 1, cmpSort(r[i]));
			//cout << "| ";
			/*for (int k = 1; k <= s[sol_index][i][0]; k++)
				cout << s[sol_index][i][k] << "," << r[i][s[sol_index][i][k]] << " ";*/
			obj_type pre_obj = c[sol_index][i];
			calculate_completion_time(sol_index, i);
			//cout << ", c: " << c[sol_index][i] << endl;
			if (c[sol_index][i] < min_completion_time)
			{
				min_completion_time = c[sol_index][i];
				min_ct_mach_index = i;
				memcpy(s[sol_index_local][min_ct_mach_index], s[sol_index][min_ct_mach_index],
					(s[sol_index][min_ct_mach_index][0] + 1)*sizeof(int));
			}
			remove_job(sol_index, i, s[sol_index_local][0][j]);
			c[sol_index][i] = pre_obj;
		}
		//cout << "min_ct: " << min_completion_time << " min_ct_mach_index: " << min_ct_mach_index << endl;
		memcpy(s[sol_index][min_ct_mach_index], s[sol_index_local][min_ct_mach_index],
			(s[sol_index_local][min_ct_mach_index][0] + 1)*sizeof(int));
		c[sol_index][min_ct_mach_index] = min_completion_time;
	}
	calculate_obj(sol_index);
}
void PMS::replace_solution(int dest, int src)
{
	for (int i = 1; i <= m; i++)
	{
		memcpy(s[dest][i], s[src][i], (n + 1)*sizeof(int));
	}
	memcpy(c[dest], c[src], (m + 1)*sizeof(obj_type));
	mm[dest] = mm[src];
}
void PMS::swap(int sol_index, int mm, int mm_j, int mo, int mo_j)
{
	// swap neighborhood
	//cout << "m,i: " << pre_mm << ", " << i << " jm,j: " << s[sol_index_local][pre_mm][jm] << ", " << s[sol_index_local][i][j] << endl;
	int temp = s[sol_index][mm][mm_j];
	s[sol_index][mm][mm_j] = s[sol_index][mo][mo_j];
	s[sol_index][mo][mo_j] = temp;
	sort(s[sol_index][mm] + 1, s[sol_index][mm] + s[sol_index][mm][0] + 1, cmpSort(r[mm]));
	calculate_completion_time(sol_index, mm);
	sort(s[sol_index][mo] + 1, s[sol_index][mo] + s[sol_index][mo][0] + 1, cmpSort(r[mo]));
	calculate_completion_time(sol_index, mo);
	calculate_obj(sol_index);
}
void PMS::insert(int sol_index, int mm, int mm_j, int mo)
{
	// insert neighborhood
	//cout << "m,i: " << pre_mm << ", " << i << " jm,j: " << s[sol_index_local][pre_mm][jm] << endl;
	s[sol_index][mo][0] += 1;
	s[sol_index][mo][s[sol_index][mo][0]] = s[sol_index][mm][mm_j];
	remove_job(sol_index, mm, s[sol_index][mm][mm_j]);
	calculate_completion_time(sol_index, mm);
	sort(s[sol_index][mo] + 1, s[sol_index][mo] + s[sol_index][mo][0] + 1, cmpSort(r[mo]));
	calculate_completion_time(sol_index, mo);
	calculate_obj(sol_index);
}
void PMS::local_search(int sol_index, int sol_index_local, NS_Mode ns)
{
	replace_solution(sol_index_local, sol_index);
	bool is_still_improved = true;
	while (is_still_improved)
	{
		is_still_improved = false;
		int pre_mm = mm[sol_index_local];// s[sol_index_local][][0];
		for (int jm = 1; jm <= s[sol_index_local][pre_mm][0]; jm++)
		{
			for (int i = 1; i <= m; i++)
			{
				if (i == s[sol_index_local][pre_mm][0])
					continue;
				if (ns == NS_Mode::INSERT)
				{
					// insert neighborhood
					//cout << "m,i: " << pre_mm << ", " << i << " jm,j: " << s[sol_index_local][pre_mm][jm] << endl;
					s[sol_index_local][i][0] += 1;
					s[sol_index_local][i][s[sol_index_local][i][0]] = s[sol_index_local][pre_mm][jm];
					remove_job(sol_index_local, pre_mm, s[sol_index_local][pre_mm][jm]);
					/*sort(s[sol_index_local][mm[sol_index_local]] + 1, s[sol_index_local][mm[sol_index_local]] +
					s[sol_index_local][mm[sol_index_local]][0] + 1, cmpSort(r[mm[sol_index_local]]));*/
					calculate_completion_time(sol_index_local, pre_mm);
					sort(s[sol_index_local][i] + 1, s[sol_index_local][i] + s[sol_index_local][i][0] + 1, cmpSort(r[i]));
					calculate_completion_time(sol_index_local, i);
					calculate_obj(sol_index_local);
					//check_solution(sol_index_local);
					if (c[sol_index][0] - c[sol_index_local][0] > MIN_EQUAL)
					{
						replace_solution(sol_index, sol_index_local);
						//cout << c[sol_index][0] << " " << c[sol_index_local][0] << " improved @" << endl;
						is_still_improved = true;
					}
					else
					{
						replace_solution(sol_index_local, sol_index);
						//cout << c[sol_index][0] << " " << c[sol_index_local][0] << " not improved" << endl;
					}
				}
				else if (ns == NS_Mode::SWAP)
				{
					for (int j = 1; j <= s[sol_index_local][i][0]; j++)
					{
						// swap neighborhood
						//cout << "m,i: " << pre_mm << ", " << i << " jm,j: " << s[sol_index_local][pre_mm][jm] << ", " << s[sol_index_local][i][j] << endl;
						int temp = s[sol_index_local][pre_mm][jm];
						s[sol_index_local][pre_mm][jm] = s[sol_index_local][i][j];
						s[sol_index_local][i][j] = temp;
						sort(s[sol_index_local][pre_mm] + 1, s[sol_index_local][pre_mm] +
							s[sol_index_local][pre_mm][0] + 1, cmpSort(r[pre_mm]));
						calculate_completion_time(sol_index_local, pre_mm);
						sort(s[sol_index_local][i] + 1, s[sol_index_local][i] + s[sol_index_local][i][0] + 1, cmpSort(r[i]));
						calculate_completion_time(sol_index_local, i);
						calculate_obj(sol_index_local);
						//check_solution(sol_index_local);
						if (c[sol_index][0] - c[sol_index_local][0] > MIN_EQUAL)
						{
							replace_solution(sol_index, sol_index_local);
							//cout << c[sol_index][0] << " " << c[sol_index_local][0] << " improved @" << endl;
							is_still_improved = true;
						}
						else
						{
							replace_solution(sol_index_local, sol_index);
							//cout << c[sol_index][0] << " " << c[sol_index_local][0] << " not improved" << endl;
						}
					}
				}
			}
		}
	}
}
void PMS::local_search_hybrid(int sol_index, int sol_index_local, NS_Mode ns)
{
	replace_solution(sol_index_local, sol_index);
	bool is_still_improved = true;
	while (is_still_improved)
	{
		is_still_improved = false;
		bool is_mm_unchanged = true;
		int pre_mm = mm[sol_index_local];// s[sol_index_local][][0];
		for (int jm = 0 % s[sol_index_local][mm[sol_index_local]][0] + 1, jm_r = 1;
		jm_r <= s[sol_index_local][mm[sol_index_local]][0] && is_mm_unchanged;
			jm_r++, jm = jm%s[sol_index_local][mm[sol_index_local]][0] + 1)
		{
			for (int i = 0 % m + 1, i_r = 1; i_r <= m && is_mm_unchanged; i_r++, i = i%m + 1)
			{
				if (i == mm[sol_index_local])
					continue;
				insert(sol_index_local, mm[sol_index_local], jm, i);
				//check_solution(sol_index_local);
				if (c[sol_index][0] - c[sol_index_local][0] > MIN_EQUAL)
				{
					replace_solution(sol_index, sol_index_local);
					//cout << c[sol_index][0] << " " << c[sol_index_local][0] << " improved @" << endl;
					is_still_improved = true;
					break;
					if (pre_mm != mm[sol_index_local])
						is_mm_unchanged = false;
				}
				else
				{
					replace_solution(sol_index_local, sol_index);
					if (pre_mm != mm[sol_index_local])
						is_mm_unchanged = false;
					//cout << c[sol_index][0] << " " << c[sol_index_local][0] << " not improved" << endl;
				}
				if (jm > s[sol_index_local][mm[sol_index_local]][0])
				{
					is_mm_unchanged = false;
				}
				//check_solution(sol_index_local);
				/*cout << mm[sol_index_local] << " " << pre_mm << endl;*/

				for (int j = rand() % (s[sol_index_local][i][0]==0?1: s[sol_index_local][i][0]) + 1, j_r = 1;
				j_r <= s[sol_index_local][i][0] && is_mm_unchanged;
					j_r++, j = j%s[sol_index_local][i][0] + 1)
				{
					swap(sol_index_local, mm[sol_index_local], jm, i, j);
					//check_solution(sol_index_local);
					if (c[sol_index][0] - c[sol_index_local][0] > MIN_EQUAL)
					{
						replace_solution(sol_index, sol_index_local);
						//cout << c[sol_index][0] << " " << c[sol_index_local][0] << " improved @" << endl;
						is_still_improved = true;
						if (pre_mm != mm[sol_index_local])
							is_mm_unchanged = false;
						break;
					}
					else
					{
						replace_solution(sol_index_local, sol_index);
						if (pre_mm != mm[sol_index_local])
							is_mm_unchanged = false;
						//cout << c[sol_index][0] << " " << c[sol_index_local][0] << " not improved" << endl;
					}
					//check_solution(sol_index_local);
				}
			}
		}
	}
}
void PMS::local_search_hybrid1(int sol_index, int sol_index_local, NS_Mode ns)
{
	bool is_still_improve = true;
	while (is_still_improve)
	{
		is_still_improve = false;
		for (int jm = 0, jm_r = 1;
		jm_r <= s[sol_index][mm[sol_index]][0]; jm_r++)
		{
			jm = jm% s[sol_index][mm[sol_index]][0] + 1;
			obj_type min_delta_f = MIN_EQUAL;
			int min_delta_f_mach = 0;
			int min_delta_f_mach_j = 0;
			for (int i = 0, i_r = 1; i_r <= m; i_r++)
			{
				i = i%m + 1;
				if (i == mm[sol_index])
					continue;
				replace_solution(sol_index_local, sol_index);
				bool test_right = true;
				if (c[sol_index_local][i] + p[i][s[sol_index_local][mm[sol_index_local]][jm]] - c[sol_index_local][0] > MIN_EQUAL)
				{
					test_right = false;
					/*cout << c[sol_index_local][i] <<" "<< p[i][s[sol_index_local][mm[sol_index_local]][jm]] << " "
						<< c[sol_index_local][0] << endl;
					display_solution(sol_index_local);*/
				}
				//check_solution(sol_index_local);
				insert(sol_index_local, mm[sol_index_local], jm, i);
				//check_solution(sol_index_local);
				if (c[sol_index][0] - c[sol_index_local][0] > min_delta_f)
				{
					min_delta_f = c[sol_index][0] - c[sol_index_local][0];
					min_delta_f_mach = i;
				}
				if (c[sol_index][0] - c[sol_index_local][0] > MIN_EQUAL&&test_right == false)
				{
					cout << "********************************" << endl;
					cout << mm[sol_index_local] << " " << jm << " " << i << endl;
					cout << c[sol_index_local][0] << " " << c[sol_index][0] << endl;
					display_solution(sol_index_local);
					system("pause");
				}

			}
			for (int i = 1; i <= m; i++)
			{
				if (i == mm[sol_index])
					continue;
				for (int j = 1; j <= s[sol_index][i][0]; j++)
				{
					replace_solution(sol_index_local, sol_index);
					swap(sol_index_local, mm[sol_index_local], jm, i, j);
					//check_solution(sol_index_local);
					if (c[sol_index][0] - c[sol_index_local][0] > min_delta_f)
					{
						min_delta_f = c[sol_index][0] - c[sol_index_local][0];
						min_delta_f_mach = i;
						min_delta_f_mach_j = j;
					}
				}
			}
			if (min_delta_f_mach != 0)
			{
				if (min_delta_f_mach_j == 0)	// insert
					insert(sol_index, mm[sol_index], jm, min_delta_f_mach);
				else
					swap(sol_index, mm[sol_index], jm, min_delta_f_mach, min_delta_f_mach_j);
				is_still_improve = true;
			}
		}
	}
}
void PMS::perturb(int sol_index, int ptr_rate)
{
	for (int i = 0; i < s[sol_index][mm[sol_index]][0] * ptr_rate*0.01; i++)
	{
		int jm = rand() % s[sol_index][mm[sol_index]][0] + 1;
		int mach_r = rand() % m + 1;
		while (mach_r == mm[sol_index]|| s[sol_index][mach_r][0] ==0)
			mach_r = rand() % m + 1;
		int jo = rand() % s[sol_index][mach_r][0] + 1;
		int temp = s[sol_index][mm[sol_index]][jm];
		s[sol_index][mm[sol_index]][jm] = s[sol_index][mach_r][jo];
		s[sol_index][mach_r][jo] = temp;
	}
	for (int i = 1; i <= m; i++)
	{
		sort(s[sol_index][i] + 1, s[sol_index][i] + s[sol_index][i][0] + 1, cmpSort(r[i]));
		calculate_completion_time(sol_index, i);
	}
	calculate_obj(sol_index);
}
void PMS::test()
{
	/*cout << endl << "This is test function." << endl;*/
	int cnt = 3;
	int si[] = { 5,1 ,2 };
	double rd[] = { 7000,2000,4000,1400,900,3400 };

	sort(si, si + cnt, cmpSort(rd));
	/*for (int i = 0; i < cnt; i++)
	cout << si[i] << "," << rd[i] << " ";
	cout << endl << "end of test function." << endl;*/

	//cout << ri << endl;
	for (int j = 0; j < 10; j++)
	{
		for (int ri = rand() % 10 + 1, i = 1; i <= 10; i++, ri = ri % 10 + 1)
		{

			cout << ri << "\t";
		}
		cout << endl;
	}
}
void PMS::iterated_local_search(int iteration, int perturb_rate, R_Mode r_mode, NS_Mode ns, int run_cnt)
{
	int sol_index_opt = 0, sol_index_best = 1,
		sol_index_cur = 2, sol_index_local = 3, sol_index_ptr = 4;
	int opt_cnt = 0, imp_cnt = 0, non_imp_cnt = 0;
	obj_type sum_obj = 0;
	int rt = time(NULL);
	rt = 1460573277;
	srand(rt);
	ofs << ins_name << "\t" << n << "\t" << m << "\t"
		<< obj_given << "\t" << c[sol_index_opt][0] << "\t" << rt << endl;
	cout << ins_name << "\t" << n << "\t" << m << "\t"
		<< obj_given << "\t" << c[sol_index_opt][0] << "\t" << rt << endl;
	for (int rc = 1; rc <= run_cnt; rc++)
	{
		start_tm = clock();
		init_solution(sol_index_cur, sol_index_local, r_mode);	//PMS::MINPDD
		//display_solution(sol_index_cur);
		local_search_hybrid(sol_index_cur, sol_index_local, ns);
		replace_solution(sol_index_best, sol_index_cur);
		//save_solution(sol_index_best, sol_index_opt, 0, run_cnt);
		end_tm = clock();
		int min_obj_iter = 0;
		for (int i = 0; i < iteration; i++)
		{
			replace_solution(sol_index_ptr, sol_index_cur);
			perturb(sol_index_ptr, perturb_rate);
			//cout << c[sol_index_best][0] << ", " << c[sol_index_cur][0] << endl;
			local_search_hybrid(sol_index_ptr, sol_index_local, ns);
			if (c[sol_index_best][0] - c[sol_index_ptr][0]>MIN_EQUAL)
			{
				replace_solution(sol_index_best, sol_index_ptr);
				//save_solution(sol_index_best, sol_index_opt, i, rc);
				min_obj_iter = i;
				end_tm = clock();
				/*if (abs(c[sol_index_best][0] - c[sol_index_opt][0]) <= MIN_EQUAL)
					break;*/
			}
			if (c[sol_index_cur][0] - c[sol_index_ptr][0] > MIN_EQUAL ||
				rand() % 100 <= (100 * exp((c[sol_index_cur][0] - c[sol_index_ptr][0]) / temperature)))
				replace_solution(sol_index_cur, sol_index_ptr);
			temperature *= control_para;
		}
		//check_solution(sol_index_best); 
		save_solution(sol_index_best, sol_index_opt, min_obj_iter, rc);
		/*display_solution(sol_index_opt);*/
		//display_solution(sol_index_best);		
	}
}
void run_algorithm(int argc, char **argv,string input_file)
{
	char *rgv[] = { "",	//0
		"_fn","total_results10",	//1,2
		"_if","instance\\BB_Problem_BestSolution\\",	//3,4	
		"_of","results\\",//5,6	
		"_p","13",		//7,8
		"_r","20",		//9,10
		"_itr","2000",	//11,12
		"_ptr","50",	//13,14
		"_rm","1",	//15,16
		"_ns","0",	//17,18		
		"_r1","1",	//19,20
		"_r2","20",	//21,22
		"_ws","0",	//23,24
		"_t","2",	//25,26
		"_cp","90"	//27,28
	};
	/*argc = sizeof(rgv) / sizeof(rgv[0]);
	argv = rgv;*/
	std::map<string, string> argv_map;
	for (int i = 1; i < argc; i += 2)//sizeof(rgv) / sizeof(rgv[0])
	{
		argv_map[string(argv[i])] = string(argv[i + 1]);
		//cout << argv[i] << ", "<<argv[i + 1] << endl;
	}
	string fnr = argv_map.at("_if") + input_file;
	string fnw = argv_map.at("_of") + argv_map.at("_fn") +
		"_p" + argv_map.at("_p") +
		"_itr" + argv_map.at("_itr") +
		"_ptr" + argv_map.at("_ptr") +
		"_rm" + argv_map.at("_rm") +
		"_ns" + argv_map.at("_ns") +
		"_r" + argv_map.at("_r1") +
		"_r" + argv_map.at("_r2") + ".txt";
	
	/*cout << "\nThis program is operated by Junwen Ding, any question please contact Ding via 769172839@qq.com." << endl;
	cout << fnr << " from " << argv_map.at("_r1") << " to " << argv_map.at("_r2") << " ILS Starts at ";
	time_t tt = time(NULL);
	tm* t = localtime(&tt);
	printf("%d-%02d-%02d %02d:%02d:%02d.",
		t->tm_year + 1900,
		t->tm_mon + 1,
		t->tm_mday,
		t->tm_hour,
		t->tm_min,
		t->tm_sec);
	cout << "..." << endl;*/
	PMS *pms = new PMS(fnr, fnw, stoi(argv_map.at("_p")));//Ni_14_4-1_1_15
	pms->ins_name = input_file;
	pms->whe_save_sol_seq = stoi(argv_map.at("_ws"));
	pms->temperature = stoi(argv_map.at("_t"));
	pms->control_para = stoi(argv_map.at("_cp"))*0.01;

	//pms->display_problem();
	//pms->display_solution(0);
	//pms->check_solution(0);
	//enum R_Mode { EXAC = 0, MINP, MAXP, MIND, MAXD, MINPDD, MAXPDD, MINPD, MAXPD, SIZE };
	PMS::R_Mode r_mode;
	if (stoi(argv_map.at("_rm")) == 1)
		r_mode = PMS::MINP;
	else if (stoi(argv_map.at("_rm")) == 2)
		r_mode = PMS::MAXP;
	else if (stoi(argv_map.at("_rm")) == 3)
		r_mode = PMS::MIND;
	else if (stoi(argv_map.at("_rm")) == 4)
		r_mode = PMS::MAXD;
	else if (stoi(argv_map.at("_rm")) == 5)
		r_mode = PMS::MINPDD;
	else if (stoi(argv_map.at("_rm")) == 6)
		r_mode = PMS::MAXPDD;
	else if (stoi(argv_map.at("_rm")) == 7)
		r_mode = PMS::MINPD;
	else
		r_mode = PMS::MAXPD;
	PMS::NS_Mode ns_mode = stoi(argv_map.at("_ns")) == 0 ? PMS::SWAP : PMS::INSERT;
	pms->iterated_local_search(stoi(argv_map.at("_itr")), stoi(argv_map.at("_ptr")),
			r_mode, ns_mode, stoi(argv_map.at("_r")));
	delete pms;
	/*cout << "\nEnds at ";
	tt = time(NULL);
	t = localtime(&tt);
	printf("%d-%02d-%02d %02d:%02d:%02d.",
		t->tm_year + 1900,
		t->tm_mon + 1,
		t->tm_mday,
		t->tm_hour,
		t->tm_min,
		t->tm_sec);
	cout << "\n" << fnr << " from " << argv_map.at("_r1") << " to " << argv_map.at("_r2") << " ILS runing completed." << endl;*/
}
int main(int argc, char **argv)
{
	int rt = time(NULL);
	//rt = 1460296365;//1459751481
	//cout << rt << endl;
	//srand(rt);
	string ins_name="Ni_50_4-2_2_9";
	/*run_algorithm(argc,argv,ins_name); */
	/*ins_name = "Ni_50_10-1_1_13";
	run_algorithm(ins_name);*/
	vector<vector<int>> n_vec = { { 8, 11, 14},{20,35,50 } };
	vector<vector<int>> m_vec = { {2,3,4},{4,7,10} };
	for (int i = 0; i < 2; i++)
	{
		for (int ni = 0; ni < 3; ni++)
		{
			for (int j = 0; j < 2; j++)
			{
				if (i != j)
					continue;
				for (int mj = 0; mj < 3; mj++)
				{
					for (int pi = 1; pi <= 2; pi++)
					{
						for (int di = 1; di <= 2; di++)
						{
							for (int ins = 1; ins <= 25; ins++)
							{
								ins_name = "Ni_" + to_string(n_vec[i][ni]) + "_" + to_string(m_vec[j][mj]) + "-"
									+ to_string(pi) + "_" + to_string(di) + "_" + to_string(ins);
								run_algorithm(argc,argv,ins_name);
							}
						}
					}
				}
			}
		}			
	}	
	system("pause");
}
#endif
