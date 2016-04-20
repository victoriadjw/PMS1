#if 1
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
	enum R_Mode { EXAC = 0, MINP, MAXP, MIND, MAXD, MINPDD, MAXPDD, MINPD, MAXPD, RANDOM, SIZE };
	enum Cmp_Mode { GREATER, LESS };
	enum NS_Mode { SWAP, INSERT };
	int whe_save_sol_seq;
	string ins_name;
	double control_para, temperature;
	PMS(string, string, int);
	~PMS();
	/*PMS(int _m, int _n, proc_type ** _p, dete_type **_d,
		int **_s_opt, obj_type _sol_obj_opt, string _file_input,
		string _file_output, int _sol_num);*/
	void display_problem();
	void display_solution(int);
	void check_solution(int);
	void init_solution(int, int, R_Mode);
	void calculate_completion_time(int, int);
	void calculate_obj(int);
	void add_job(int, int, int);
	void remove_job(int, int, int);
	void local_search(int, int, NS_Mode);
	void local_search_hybrid(int, int, NS_Mode);
	void local_search_hybrid1(int, int, NS_Mode);
	void local_search_ejection_chain(int, int, NS_Mode);
	void iterated_local_search(int, int, R_Mode, NS_Mode, int, int);
	void replace_solution(int, int);
	void perturb(int, int);
	void swap(int, int, int, int, int);
	void insert(int, int, int, int);
	void trail_move(int, int, int, int, int);
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
	int si_opt = 0;
	for (int i = 1; i <= m; i++)
	{
		ifs >> s[si_opt][i][0];	// the optimal number of jobs to be assigned to machine i
		for (int j = 1; j <= s[si_opt][i][0]; j++)
		{
			ifs >> s[si_opt][i][j];	// the list of jobs to be assinged to machine i
		}
	}
	ifs >> obj_given;
	for (int i = 1; i <= m; i++)
	{
		sort(s[si_opt][i] + 1, s[si_opt][i] + s[si_opt][i][0] + 1, cmpSort(r[i]));
		calculate_completion_time(si_opt, i);
	}
	calculate_obj(si_opt);
	if (abs(obj_given - c[si_opt][0]) > MIN_EQUAL)
	{
		cout << "the given optimal solution is wrong."
			<< obj_given << " " << c[si_opt][0] << endl;
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
void PMS::save_solution(int si, int si_opt, int iterration, int run_cnt)
{
	//end_tm = clock();
	int result_improve, given_result_improve;
	if (abs(c[si][0] - c[si_opt][0]) <= MIN_EQUAL)
		result_improve = 1;	// equal
	else if (c[si_opt][0] - c[si][0] > MIN_EQUAL)
		result_improve = 2;	// improved
	else
		result_improve = 0;	// not improved
	if (abs(c[si][0] - obj_given) <= MIN_EQUAL)
		given_result_improve = 1;
	else if (obj_given - c[si][0] > MIN_EQUAL)
		given_result_improve = 2;
	else
		given_result_improve = 0;
	ofs << run_cnt << "\t"
		//<< c[si_opt][0] << "\t"
		<< c[si][0] << "\t"
		<< mm[si] << "\t"
		<< iterration << "\t"
		<< (end_tm - start_tm) /*/ CLOCKS_PER_SEC*/ << "\t"
		<< given_result_improve << "\t"
		<< result_improve
		<< endl;
	//cout << run_cnt << "\t"
	//	//<< c[si_opt][0] << "\t"
	//	<< c[si][0] << "\t"
	//	<< mm[si] << "\t"
	//	<< iterration << "\t"
	//	<< (end_tm - start_tm) /*/ CLOCKS_PER_SEC*/ << "\t"
	//	<< given_result_improve << "\t"
	//	<< result_improve
	//	<< endl;
	if (whe_save_sol_seq)	// whe_save_sol_seq
	{
		for (int i = 1; i <= m; i++)
		{
			ofs << c[si][i] << "\t"
				<< s[si][i][0] << "\t";
			for (int j = 1; j <= s[si][i][0]; j++)
				ofs << s[si][i][j] << "\t";
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
void PMS::display_solution(int si)
{
	cout << "*** display si:" << si
		<< ", obj: " /*<< fixed << setprecision(4)*/
		<< c[si][0] << ", " << mm[si] << " ***" << endl;
	for (int i = 1; i <= m; i++)
	{
		obj_type ir = 0;
		cout << s[si][i][0] << ", " << c[si][i] << ": ";
		for (int j = 1; j <= s[si][i][0]; j++)
		{
			cout << s[si][i][j] /*<< ", " << r[i][s[si][i][j]] */ << "\t";
			ir += (j*r[i][s[si][i][j]]);
		}
		cout /*<<ir*/ << endl;
	}
}
void PMS::check_solution(int si)
{
	cout << "*** check si:" << si
		<< ", obj: " << c[si][0] << " ***" << endl;
	perf_type ql;
	obj_type cl, max_cl = 0;
	int max_cl_mach_index, sum_job_index = 0;
	for (int i = 1; i <= m; i++)
	{
		ql = 1;
		cl = 0;
		//cout << cl << " ";
		for (int j = 1; j <= s[si][i][0]; j++)
		{
			cl += p[i][s[si][i][j]] / ql;
			ql *= (1 - d[i][s[si][i][j]]);
			sum_job_index += s[si][i][j];
			//cout << cl << " ";
			if (j > 1 && r[i][s[si][i][j - 1]] - r[i][s[si][i][j]] < MIN_EQUAL)
			{
				cout << "ERROR, si: " << si << " not sort by r"
					<< r[i][s[si][i][j - 1]] << " " << r[i][s[si][i][j]] << endl;
				display_solution(si);
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
	if (abs(max_cl - c[si][0])>MIN_EQUAL/*||max_cl_mach_index!=mm[si]*/)
	{
		cout << "ERROR, si: " << si
			<< ", real obj: " << max_cl << " " << c[si][0] << endl;
		system("pause");
	}
	if (sum_job_index != (1 + n)*n / 2)
	{
		cout << "ERROR, si: " << si
			<< ", multiple job index" << endl;
		display_solution(si);
		system("pause");
	}
}
void PMS::calculate_completion_time(int si, int mach_index)
{
	/*if (s[si][mach_index][0] == 0)
	{
		c[si][mach_index] = 0;
		return;
	}
	perf_type q = 1;
	c[si][mach_index] = p[mach_index][s[si][mach_index][1]];
	for (int i = 2; i <= s[si][mach_index][0]; i++)
	{
		q *= (1 - d[mach_index][s[si][mach_index][i - 1]]);
		c[si][mach_index] += p[mach_index][s[si][mach_index][i]] / q;
	}*/
	perf_type q = 1;
	c[si][mach_index] = 0;
	for (int i = 1; i <= s[si][mach_index][0]; i++)
	{
		c[si][mach_index] += p[mach_index][s[si][mach_index][i]] / q;
		q *= (1 - d[mach_index][s[si][mach_index][i]]);
	}
}
void PMS::calculate_obj(int si)
{
	c[si][0] = 0;
	for (int i = 1; i <= m; i++)
	{
		if ((c[si][i] - c[si][0]) > MIN_EQUAL)
			c[si][0] = c[si][i];
	}
	for (int i = 1; i <= m; i++)
	{
		if (abs(c[si][i] - c[si][0]) <= MIN_EQUAL)
			mm[si] = i;
	}
}
void PMS::add_job(int si, int mach_index, int job)
{
	s[si][mach_index][0] += 1;
	s[si][mach_index][s[si][mach_index][0]] = job;
}
void PMS::remove_job(int si, int mach_index, int job)
{
	for (int i = 1; i <= s[si][mach_index][0]; i++)
	{
		if (s[si][mach_index][i] == job)
		{
			for (int j = i + 1; j <= s[si][mach_index][0]; j++)
				s[si][mach_index][j - 1] = s[si][mach_index][j];
			break;
		}
	}
	s[si][mach_index][0] -= 1;
}
void PMS::init_solution(int si, int si_local, R_Mode rm)
{
	if (rm == R_Mode::RANDOM)
	{
		for (int j = 1; j <= n; j++)
			charact[rm][j] = rand() % 1000;
	}
	for (int i = 1; i <= n; i++)
		s[si_local][0][i] = i;	// the machine 0 as the temp memory
	sort(s[si_local][0] + 1, s[si_local][0] + n + 1, cmpSort(charact[rm], Cmp_Mode::GREATER));
	/*for (int i = 1; i <= n; i++)
	{
		cout << s[si_local][0][i] << " " << charact[rm][s[si_local][0][i]] << "\t";
	}
	cout << endl;*/
	for (int i = 1; i <= m; i++)
	{
		s[si][i][0] = 0;	// init the number of jobs assigned to machine i
		s[si_local][i][0] = 0;
	}
	for (int j = 1; j <= n; j++)
	{
		/*cout << "job " << j << ": " << s[si_local][0][j] << ", ";
		for (int i = 1; i <= m; i++)
			cout << r[i][s[si_local][0][j]] << "\t";
		cout << endl;*/
		proc_type min_completion_time = DBL_MAX;
		int min_ct_mach_index;
		for (int i = 1; i <= m; i++)
		{
			//cout << "machine " << i << ": ";
			s[si][i][0] += 1;
			s[si][i][s[si][i][0]] = s[si_local][0][j];
			/*for (int k = 1; k <= s[si][i][0]; k++)
				cout << s[si][i][k] << "," << r[i][s[si][i][k]] << " ";*/
			sort(s[si][i] + 1, s[si][i] + s[si][i][0] + 1, cmpSort(r[i]));
			//cout << "| ";
			/*for (int k = 1; k <= s[si][i][0]; k++)
				cout << s[si][i][k] << "," << r[i][s[si][i][k]] << " ";*/
			obj_type pre_obj = c[si][i];
			calculate_completion_time(si, i);
			//cout << ", c: " << c[si][i] << endl;
			if (c[si][i] < min_completion_time)
			{
				min_completion_time = c[si][i];
				min_ct_mach_index = i;
				memcpy(s[si_local][min_ct_mach_index], s[si][min_ct_mach_index],
					(s[si][min_ct_mach_index][0] + 1)*sizeof(int));
			}
			remove_job(si, i, s[si_local][0][j]);
			c[si][i] = pre_obj;
		}
		//cout << "min_ct: " << min_completion_time << " min_ct_mach_index: " << min_ct_mach_index << endl;
		memcpy(s[si][min_ct_mach_index], s[si_local][min_ct_mach_index],
			(s[si_local][min_ct_mach_index][0] + 1)*sizeof(int));
		c[si][min_ct_mach_index] = min_completion_time;
	}
	calculate_obj(si);
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
void PMS::swap(int si, int mm, int mm_j, int mo, int mo_j)
{
	// swap neighborhood
	//cout << "m,i: " << pre_mm << ", " << i << " jm,j: " << s[si_local][pre_mm][jm] << ", " << s[si_local][i][j] << endl;
	int temp = s[si][mm][mm_j];
	s[si][mm][mm_j] = s[si][mo][mo_j];
	s[si][mo][mo_j] = temp;
	sort(s[si][mm] + 1, s[si][mm] + s[si][mm][0] + 1, cmpSort(r[mm]));
	calculate_completion_time(si, mm);
	sort(s[si][mo] + 1, s[si][mo] + s[si][mo][0] + 1, cmpSort(r[mo]));
	calculate_completion_time(si, mo);
	calculate_obj(si);
}
void PMS::insert(int si, int mm, int mm_j, int mo)
{
	// insert neighborhood
	//cout << "m,i: " << pre_mm << ", " << i << " jm,j: " << s[si_local][pre_mm][jm] << endl;
	s[si][mo][0] += 1;
	s[si][mo][s[si][mo][0]] = s[si][mm][mm_j];
	remove_job(si, mm, s[si][mm][mm_j]);
	calculate_completion_time(si, mm);
	sort(s[si][mo] + 1, s[si][mo] + s[si][mo][0] + 1, cmpSort(r[mo]));
	calculate_completion_time(si, mo);
	calculate_obj(si);
}
void PMS::trail_move(int si, int m_mach, int ejected_job, int o_mach, int ejected_triggered_job)
{
	remove_job(si, o_mach, ejected_triggered_job);
	add_job(si, m_mach, ejected_triggered_job);
	add_job(si, o_mach, ejected_job);
	sort(s[si][m_mach] + 1, s[si][m_mach] + s[si][m_mach][0] + 1, cmpSort(r[m_mach]));
	calculate_completion_time(si, m_mach);
	sort(s[si][o_mach] + 1, s[si][o_mach] + s[si][o_mach][0] + 1, cmpSort(r[o_mach]));
	calculate_completion_time(si, o_mach);
	calculate_obj(si);
}
void PMS::local_search(int si, int si_local, NS_Mode ns)
{
	replace_solution(si_local, si);
	bool is_still_improved = true;
	while (is_still_improved)
	{
		is_still_improved = false;
		int pre_mm = mm[si_local];// s[si_local][][0];
		for (int jm = 1; jm <= s[si_local][pre_mm][0]; jm++)
		{
			for (int i = 1; i <= m; i++)
			{
				if (i == s[si_local][pre_mm][0])
					continue;
				if (ns == NS_Mode::INSERT)
				{
					// insert neighborhood
					//cout << "m,i: " << pre_mm << ", " << i << " jm,j: " << s[si_local][pre_mm][jm] << endl;
					s[si_local][i][0] += 1;
					s[si_local][i][s[si_local][i][0]] = s[si_local][pre_mm][jm];
					remove_job(si_local, pre_mm, s[si_local][pre_mm][jm]);
					/*sort(s[si_local][mm[si_local]] + 1, s[si_local][mm[si_local]] +
					s[si_local][mm[si_local]][0] + 1, cmpSort(r[mm[si_local]]));*/
					calculate_completion_time(si_local, pre_mm);
					sort(s[si_local][i] + 1, s[si_local][i] + s[si_local][i][0] + 1, cmpSort(r[i]));
					calculate_completion_time(si_local, i);
					calculate_obj(si_local);
					//check_solution(si_local);
					if (c[si][0] - c[si_local][0] > MIN_EQUAL)
					{
						replace_solution(si, si_local);
						//cout << c[si][0] << " " << c[si_local][0] << " improved @" << endl;
						is_still_improved = true;
					}
					else
					{
						replace_solution(si_local, si);
						//cout << c[si][0] << " " << c[si_local][0] << " not improved" << endl;
					}
				}
				else if (ns == NS_Mode::SWAP)
				{
					for (int j = 1; j <= s[si_local][i][0]; j++)
					{
						// swap neighborhood
						//cout << "m,i: " << pre_mm << ", " << i << " jm,j: " << s[si_local][pre_mm][jm] << ", " << s[si_local][i][j] << endl;
						int temp = s[si_local][pre_mm][jm];
						s[si_local][pre_mm][jm] = s[si_local][i][j];
						s[si_local][i][j] = temp;
						sort(s[si_local][pre_mm] + 1, s[si_local][pre_mm] +
							s[si_local][pre_mm][0] + 1, cmpSort(r[pre_mm]));
						calculate_completion_time(si_local, pre_mm);
						sort(s[si_local][i] + 1, s[si_local][i] + s[si_local][i][0] + 1, cmpSort(r[i]));
						calculate_completion_time(si_local, i);
						calculate_obj(si_local);
						//check_solution(si_local);
						if (c[si][0] - c[si_local][0] > MIN_EQUAL)
						{
							replace_solution(si, si_local);
							//cout << c[si][0] << " " << c[si_local][0] << " improved @" << endl;
							is_still_improved = true;
						}
						else
						{
							replace_solution(si_local, si);
							//cout << c[si][0] << " " << c[si_local][0] << " not improved" << endl;
						}
					}
				}
			}
		}
	}
}
void PMS::local_search_hybrid(int si, int si_local, NS_Mode ns)
{
	replace_solution(si_local, si);
	bool is_still_improved = true;
	while (is_still_improved)
	{
		is_still_improved = false;
		bool is_mm_unchanged = true;
		int pre_mm = mm[si_local];// s[si_local][][0];
		for (int jm = 0 % s[si_local][mm[si_local]][0] + 1, jm_r = 1;
		jm_r <= s[si_local][mm[si_local]][0] && is_mm_unchanged;
			jm_r++, jm = jm%s[si_local][mm[si_local]][0] + 1)
		{
			for (int i = 0 % m + 1, i_r = 1; i_r <= m && is_mm_unchanged; i_r++, i = i%m + 1)
			{
				if (i == mm[si_local])
					continue;
				insert(si_local, mm[si_local], jm, i);
				//check_solution(si_local);
				if (pre_mm != mm[si_local])
					is_mm_unchanged = false;
				if (c[si][0] - c[si_local][0] > MIN_EQUAL)
				{
					replace_solution(si, si_local);
					//cout << c[si][0] << " " << c[si_local][0] << " improved @" << endl;
					is_still_improved = true;
					break;
				}
				else
				{
					replace_solution(si_local, si);
					//cout << c[si][0] << " " << c[si_local][0] << " not improved" << endl;
				}
				//check_solution(si_local);
				/*cout << mm[si_local] << " " << pre_mm << endl;*/

				for (int j = rand() % (s[si_local][i][0] == 0 ? 1 : s[si_local][i][0]) + 1, j_r = 1;
				j_r <= s[si_local][i][0] && is_mm_unchanged;
					j_r++, j = j%s[si_local][i][0] + 1)
				{
					swap(si_local, mm[si_local], jm, i, j);
					if (pre_mm != mm[si_local])
						is_mm_unchanged = false;
					//check_solution(si_local);
					if (c[si][0] - c[si_local][0] > MIN_EQUAL)
					{
						replace_solution(si, si_local);
						//cout << c[si][0] << " " << c[si_local][0] << " improved @" << endl;
						is_still_improved = true;
						break;
					}
					else
					{
						replace_solution(si_local, si);
						//cout << c[si][0] << " " << c[si_local][0] << " not improved" << endl;
					}
					//check_solution(si_local);
				}
			}
		}
	}
}
void PMS::local_search_ejection_chain(int si, int si_local, NS_Mode ns)
{
	int si_ref = 7, si_trial = 8,si_cur=9;
	replace_solution(si_cur, si);
	bool is_still_improved = true;
	while (is_still_improved)
	{
		is_still_improved = false;
		bool is_mm_unchanged = true;
		int pre_mm = mm[si_cur];// s[si_local][][0];
		for (int jm = 0 % s[si_cur][mm[si_cur]][0] + 1, jm_r = 1;
		jm_r <= s[si_cur][mm[si_cur]][0]&&is_mm_unchanged;
			jm_r++, jm = jm%s[si_cur][mm[si_cur]][0] + 1)
		{
			// job s[si_local][mm[si_index]][jm] is ejected
			/*cout << "mm1: " << mm[si_local]
				<< "\t eje: " << s[si_local][mm[si_local]][jm]
				<< "\tf: " << prev_f << endl;
			for (int im = 1; im <= s[si_local][mm[si_local]][0]; im++)
				cout << s[si_local][mm[si_local]][im] << "\t";
			cout << endl;*/
			replace_solution(si_ref, si_cur);
			obj_type prev_f = c[si_ref][0];
			int ejected_job = s[si_ref][mm[si_ref]][jm];
			remove_job(si_ref, mm[si_ref], s[si_ref][mm[si_ref]][jm]);
			cout << prev_f << "\t" << mm[si_ref] << "\t" << ejected_job << endl;
			display_solution(si_ref);
			int len = 1;
			while (len <= 4)
			{
				int min_mach = 0, min_ejected_triggered_job = 0;
				obj_type min_df = DBL_MAX;
				if (rand() % 100 < 70)	// greedy strategy
				{
					for (int i = 0 % m + 1, i_r = 1; i_r <= m; i_r++, i = i%m + 1)
					{
						if (i == mm[si_ref] || s[si_ref][i][0] == 0)
							continue;
						for (int j = 0 % s[si_ref][i][0] + 1, j_r = 1;
						j_r <= s[si_ref][i][0]; j_r++, j = j%s[si_ref][i][0] + 1)
						{
							add_job(si_ref, mm[si_ref], s[si_ref][i][j]);
							sort(s[si_ref][mm[si_ref]] + 1,
								s[si_ref][mm[si_ref]] + s[si_ref][mm[si_ref]][0] + 1,
								cmpSort(r[mm[si_ref]]));
							calculate_completion_time(si_ref, mm[si_ref]);
							if (c[si_ref][mm[si_ref]] - prev_f < min_df&&c[si_ref][i]-p[i][s[si_ref][i][j]]+p[i][ejected_job]<prev_f)
							{
								min_df = c[si_ref][mm[si_ref]] - prev_f;
								min_mach = i;
								min_ejected_triggered_job = s[si_ref][i][j];
							}
							remove_job(si_ref, mm[si_ref], s[si_ref][i][j]);
						}
					}
				}
				else
				{
					min_mach = rand() % m + 1;
					while (min_mach == mm[si_ref] || s[si_ref][min_mach][0] == 0)
						min_mach = rand() % m + 1;
					min_ejected_triggered_job = s[si_ref][min_mach][rand() % s[si_ref][min_mach][0] + 1];
				}
				/*cout << s[si_local][mm[si_local]][jm] << endl;
				cout << "min_m1: " << min_mach << ", min_j: " <<min_ejected_triggered_job
					<< ", ct: " << c[si_local][min_mach] << endl;
				for (int im = 1; im <= s[si_local][min_mach][0]; im++)
					cout << s[si_local][min_mach][im] << "\t";
				cout << endl;*/

				cout << mm[si_ref] << "\t" << ejected_job << "\t"
					<< min_mach << "\t" << min_ejected_triggered_job << "\t";
				trail_move(si_ref, mm[si_ref], ejected_job, min_mach, min_ejected_triggered_job);
				
				cout<< c[si_cur][0] << "\t" << c[si_ref][0] << endl;
				if (prev_f - c[si_ref][0] > MIN_EQUAL)
				{
					replace_solution(si_cur, si_ref);
					is_still_improved = true;
					
					//check_solution(si_local);
					break;
				}
				else
				{
					remove_job(si_ref, min_mach, ejected_job);
				}
				display_solution(si_ref);
				len += 1;
			}
		}
	}
}
void PMS::local_search_hybrid1(int si, int si_local, NS_Mode ns)
{
	bool is_still_improve = true;
	while (is_still_improve)
	{
		is_still_improve = false;
		for (int jm = 0, jm_r = 1;
		jm_r <= s[si][mm[si]][0]; jm_r++)
		{
			jm = jm% s[si][mm[si]][0] + 1;
			obj_type min_delta_f = MIN_EQUAL;
			int min_delta_f_mach = 0;
			int min_delta_f_mach_j = 0;
			for (int i = 0, i_r = 1; i_r <= m; i_r++)
			{
				i = i%m + 1;
				if (i == mm[si])
					continue;
				replace_solution(si_local, si);
				bool test_right = true;
				if (c[si_local][i] + p[i][s[si_local][mm[si_local]][jm]] - c[si_local][0] > MIN_EQUAL)
				{
					test_right = false;
					/*cout << c[si_local][i] <<" "<< p[i][s[si_local][mm[si_local]][jm]] << " "
						<< c[si_local][0] << endl;
					display_solution(si_local);*/
				}
				//check_solution(si_local);
				insert(si_local, mm[si_local], jm, i);
				//check_solution(si_local);
				if (c[si][0] - c[si_local][0] > min_delta_f)
				{
					min_delta_f = c[si][0] - c[si_local][0];
					min_delta_f_mach = i;
				}
				if (c[si][0] - c[si_local][0] > MIN_EQUAL&&test_right == false)
				{
					cout << "********************************" << endl;
					cout << mm[si_local] << " " << jm << " " << i << endl;
					cout << c[si_local][0] << " " << c[si][0] << endl;
					display_solution(si_local);
					system("pause");
				}

			}
			for (int i = 1; i <= m; i++)
			{
				if (i == mm[si])
					continue;
				for (int j = 1; j <= s[si][i][0]; j++)
				{
					replace_solution(si_local, si);
					swap(si_local, mm[si_local], jm, i, j);
					//check_solution(si_local);
					if (c[si][0] - c[si_local][0] > min_delta_f)
					{
						min_delta_f = c[si][0] - c[si_local][0];
						min_delta_f_mach = i;
						min_delta_f_mach_j = j;
					}
				}
			}
			if (min_delta_f_mach != 0)
			{
				if (min_delta_f_mach_j == 0)	// insert
					insert(si, mm[si], jm, min_delta_f_mach);
				else
					swap(si, mm[si], jm, min_delta_f_mach, min_delta_f_mach_j);
				is_still_improve = true;
			}
		}
	}
}
void PMS::perturb(int si, int ptr_rate)
{
	for (int i = 0; i < s[si][mm[si]][0] * ptr_rate*0.01; i++)
	{
		int jm = rand() % s[si][mm[si]][0] + 1;
		int mach_r = rand() % m + 1;
		while (mach_r == mm[si] || s[si][mach_r][0] == 0)
			mach_r = rand() % m + 1;
		int jo = rand() % s[si][mach_r][0] + 1;
		int temp = s[si][mm[si]][jm];
		s[si][mm[si]][jm] = s[si][mach_r][jo];
		s[si][mach_r][jo] = temp;
	}
	for (int i = 1; i <= m; i++)
	{
		sort(s[si][i] + 1, s[si][i] + s[si][i][0] + 1, cmpSort(r[i]));
		calculate_completion_time(si, i);
	}
	calculate_obj(si);
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
void PMS::iterated_local_search(int iteration, int perturb_rate, R_Mode r_mode, NS_Mode ns, int run_cnt_from, int run_cnt_to)
{
	int si_opt = 0, si_best = 1,
		si_cur = 2, si_local = 3, si_ptr = 4;
	int opt_cnt = 0, imp_cnt = 0, non_imp_cnt = 0;
	obj_type sum_obj = 0;
	int rt = time(NULL);
	rt = 1461068961;
	srand(rt);
	ofs << ins_name << "\t" << n << "\t" << m << "\t"
		<< obj_given << "\t" << c[si_opt][0] << "\t" << rt << endl;
	cout << ins_name << "\t" << n << "\t" << m << "\t"
		<< obj_given << "\t" << c[si_opt][0] << "\t" << rt << endl;
	for (int rc = run_cnt_from; rc <= run_cnt_to; rc++)
	{
		start_tm = clock();
		init_solution(si_cur, si_local, r_mode);	//PMS::MINPDD
		//display_solution(si_cur);
		local_search_ejection_chain(si_cur, si_local, ns);
		replace_solution(si_best, si_cur);
		//save_solution(si_best, si_opt, 0, run_cnt);
		end_tm = clock();
		int min_obj_iter = 0;
		if (n <= 14)	// find the optimal solution for instances in OB set
		{
			if (abs(c[si_best][0] - c[si_opt][0]) <= MIN_EQUAL)
				iteration = 0;
			else
				iteration = INT16_MAX;
		}
		for (int i = 0; i < iteration; i++)
		{
			replace_solution(si_ptr, si_cur);
			perturb(si_ptr, perturb_rate);
			//cout << c[si_best][0] << ", " << c[si_cur][0] << endl;
			local_search_ejection_chain(si_ptr, si_local, ns);
			if (c[si_best][0] - c[si_ptr][0]>MIN_EQUAL)
			{
				replace_solution(si_best, si_ptr);
				//save_solution(si_best, si_opt, i, rc);
				min_obj_iter = i;
				end_tm = clock();
				if (n <= 14 && abs(c[si_best][0] - c[si_opt][0]) <= MIN_EQUAL)
					break;	// find the optimal solution for instances in OB set
			}
			if (c[si_cur][0] - c[si_ptr][0] > MIN_EQUAL ||
				rand() % 100 <= (100 * exp((c[si_cur][0] - c[si_ptr][0]) / temperature)))
				replace_solution(si_cur, si_ptr);
			temperature *= control_para;
		}
		//check_solution(si_best); 
		save_solution(si_best, si_opt, min_obj_iter, rc);
		/*display_solution(si_opt);*/
		//display_solution(si_best);		
	}
}
void run_algorithm(std::map<string, string> &argv_map, string ins_name)
{
	string fnr = argv_map.at("_if") + ins_name;
	string fnw = argv_map.at("_of") + argv_map.at("_px") +
		"_p" + argv_map.at("_p") +
		"_itr" + argv_map.at("_itr") +
		"_ptr" + argv_map.at("_ptr") +
		"_rm" + argv_map.at("_rm") +
		"_ns" + argv_map.at("_ns") +
		"_t" + argv_map.at("_t") +
		"_cp" + argv_map.at("_cp") +
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
	pms->ins_name = ins_name;
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
	else if (stoi(argv_map.at("_rm")) == 8)
		r_mode = PMS::MAXPD;
	else
		r_mode = PMS::RANDOM;
	PMS::NS_Mode ns_mode = stoi(argv_map.at("_ns")) == 0 ? PMS::SWAP : PMS::INSERT;
	pms->iterated_local_search(stoi(argv_map.at("_itr")), stoi(argv_map.at("_ptr")),
		r_mode, ns_mode, stoi(argv_map.at("_r1")), stoi(argv_map.at("_r2")));
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
	char *rgv_ins[] = { "",
		"_vi1","0",		"_vi2","2",
		"_ni1","0",		"_ni2","3",
		"_vj1","0",		"_vj2","2",
		"_mj1","0",		"_mj2","3",
		"_pi1","1",		"_pi2","2",
		"_di1","1",		"_di2","2",
		"_ins1","1",	"_ins2","25"
	};
	char *rgv[] = { "",	//0
		"_px","tr",	//prefix for output file name
		"_if","instance\\BB_Problem_BestSolution\\",	//input file directory
		"_of","results\\",// output file directory
		"_p","13",		// population size
		"_itr","2000",	// max iteration of ILS
		"_ptr","50",	// perturbation rate 
		"_rm","1",	// construction rules for initial solution
		"_ns","0",	// neighborhood search, 0:swap, 1:insert	
		"_r1","1",	// run cnt from
		"_r2","20",	// run cnt to
		"_ws","0",	// whe_save_sol_seq
		"_t","2",	// initial temperature T
		"_cp","70",	// control para cp, T=T*cp
		"_vi1","1",		"_vi2","2",
		"_ni1","2",		"_ni2","3",
		"_vj1","0",		"_vj2","2",
		"_mj1","1",		"_mj2","2",
		"_pi1","1",		"_pi2","1",
		"_di1","2",		"_di2","2",
		"_ins1","1",	"_ins2","1"
	};
	argc = sizeof(rgv) / sizeof(rgv[0]); argv = rgv;
	std::map<string, string> argv_map;
	for (int i = 1; i < sizeof(rgv_ins) / sizeof(rgv_ins[0]); i += 2)
		argv_map[string(rgv_ins[i])] = string(rgv_ins[i + 1]);
	for (int i = 1; i < argc; i += 2)
		argv_map[string(argv[i])] = string(argv[i + 1]);
	vector<vector<int>> n_vec = { { 8, 11, 14},{20,35,50 } };
	vector<vector<int>> m_vec = { {2,3,4},{4,7,10} };
	for (int vi = stoi(argv_map.at("_vi1")); vi < stoi(argv_map.at("_vi2")); vi++)	// 0, 2
	{
		for (int ni = stoi(argv_map.at("_ni1")); ni < stoi(argv_map.at("_ni2")); ni++)	// 0, 3
		{
			for (int vj = stoi(argv_map.at("_vj1")); vj < stoi(argv_map.at("_vj2")); vj++)	// 0, 2
			{
				if (vi != vj)
					continue;
				for (int mj = stoi(argv_map.at("_mj1")); mj < stoi(argv_map.at("_mj2")); mj++)	// 0, 3
				{
					for (int pi = stoi(argv_map.at("_pi1")); pi <= stoi(argv_map.at("_pi2")); pi++)	// 1, 2
					{
						for (int di = stoi(argv_map.at("_di1")); di <= stoi(argv_map.at("_di2")); di++)	// 1, 2
						{
							for (int ins = stoi(argv_map.at("_ins1")); ins <= stoi(argv_map.at("_ins2")); ins++)	// 1, 25
							{
								run_algorithm(argv_map, "Ni_" + to_string(n_vec[vi][ni]) + "_" + to_string(m_vec[vj][mj]) + "-"
									+ to_string(pi) + "_" + to_string(di) + "_" + to_string(ins));
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
