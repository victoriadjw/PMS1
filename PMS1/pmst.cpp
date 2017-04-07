#if 1
#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<fstream>
#include<string>
#include<time.h>
#include<iomanip>
#include<algorithm>
#include<math.h>
#include<float.h>
#include<map>
#include<vector>
#define DEBUG 
using namespace std;
typedef double proc_type;	// type of processing time
typedef double dete_type;	// type of deterioration effect
typedef double obj_type;	// type of objective value
typedef double perf_type;	// type of performance level
class PMS
{
private:
	int m, n;	// m and n are the number of machines and jobs
	vector<vector<proc_type>> p;	// processing time
	vector<vector<dete_type>> d;	// deterioration effect
	obj_type obj_given, makespan;	// the optimal objective value, makespan for total completion time
	string file_input, file_output;
	int sol_num;	// number of solution

	vector<vector<vector<int>>> s;	// solution, s[si][0][0] is the makespan machine
	vector<vector<vector<perf_type>>> q;	// perfermance level
	vector<vector<vector<obj_type>>> c;	// completion time, c[si][0][0] is the obj
	vector<vector<bool>> effe_mach;	// indicates the effectiveness of the machine
	vector<int> mm;	// indicates the makespan machine
	vector<obj_type> sol_obj;	// the objective value of the solution

	vector<vector<perf_type>> r;	// exact sort accordance
	vector<vector<perf_type>> charact;	// characteristics sort accordance
	vector<int>tabu_pool_update;	// tabu list in pool update
									//vector<vector<obj_type>> dis;	// distance between two parents in pool update rule
	clock_t start_tm, end_tm;
	class cmpSort;
	const double MIN_EQUAL = 0.001;
	int rand_seed;
	ofstream ofs;
	int ls_cnt1, ls_cnt2, rc, gen;
public:
	enum R_Mode { EXAC = 0, MINP, MAXP, MIND, MAXD, MINPDD, MAXPDD, MINPD, MAXPD, RANDOM, SIZE };
	enum Cmp_Mode { GREATER, LESS };
	enum NS_Mode { SWAP, INSERT };
	int whe_save_sol_seq, alpha_cx_ptr, non_popu, whe_dc, ls_method, pu_method, cx_method, cl, am;
	string ins_name;
	double control_para, temperature, pu_beta;
	PMS(string, string, int);
	~PMS();
	void display_problem();
	void display_solution(int);
	void check_solution(int);
	void check_solution_machine(int, int);
	void init_solution(int, R_Mode);
	void calculate_completion_time(int, int);
	void calculate_obj(int);
	void add_job(int, int, int&, int);
	void remove_job(int, int, int&, int);
	void exchange_job(int, int, int&, int&, int, int);
	void local_search(int);
	void local_search_hybrid(int);
	void local_search_hybrid_tri_insert(int);
	void local_search_hybrid_tri_swap(int);
	void local_search_hybrid_tri_insert_swap(int);
	void local_search_ejection_chain(int, int, NS_Mode);
	void iterated_local_search(int, int, R_Mode, NS_Mode, int, int);
	void ejection_chain_local_search(int, int, R_Mode, NS_Mode, int, int);
	void hma(int, int, R_Mode, NS_Mode, int, int);
	void gpx(int, int, int);
	void mpx(int, int, int);
	void ufx(int, int, int);
	void pool_update(int, int, int);
	void pool_update_qd(int, int, int);
	void divide_and_conquer(int si, int mach_num);
	void local_search1(int si, int *sub_mach, int sub_m, obj_type &sub_obj);
	void insert1(int si, int *sub_mach, int mm_j, int mo, obj_type &sub_obj);
	void swap1(int si, int *sub_mach, int mm_j, int mo, int mo_j, obj_type &sub_obj);
	void replace_solution(int, int);
	void perturb(int, int);
	void perturb1(int, int);
	void swap(int, int, int, int, int);
	void insert(int, int, int, int);
	void save_solution(int, int, int, int);
	void test();
};

class PMS::cmpSort
{
public:
	cmpSort(vector<proc_type> &as, Cmp_Mode _cm = Cmp_Mode::GREATER) :
		array_sort(as), cm(_cm) {}
	bool operator()(int &l, int &r)
	{
		return (cm == Cmp_Mode::GREATER) ? (array_sort[l] > array_sort[r]) :
			(array_sort[l] < array_sort[r]);
	}
private:
	vector<proc_type> &array_sort;
	Cmp_Mode cm;
};

PMS::PMS(string file_input, string file_output, int _sol_num) :sol_num(_sol_num)
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
	p.resize(m + 1, vector<proc_type>(n + 1, 0));
	d.resize(m + 1, vector<proc_type>(n + 1, 0));
	r.resize(m + 1, vector<proc_type>(n + 1, 0));
	charact.resize(R_Mode::SIZE, vector<perf_type>(n + 1));

	s.resize(sol_num, vector<vector<int>>(m + 1, vector<int>(1, 0)));
	effe_mach.resize(sol_num, vector<bool>(m + 1, false));
	mm.resize(sol_num, 0);
	q.resize(sol_num, vector<vector<perf_type>>(m + 1, vector<perf_type>(1, 1)));
	c.resize(sol_num, vector<vector<obj_type>>(m + 1, vector<obj_type>(1, 0)));
	sol_obj.resize(sol_num, 0);
	//dis.resize(sol_num - non_popu + 2, vector<obj_type>(sol_num - non_popu + 2, 0));
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
	int si_opt = 0, num_job;
	for (int i = 1; i <= m; i++)
	{
		ifs >> num_job;	// the optimal number of jobs to be assigned to machine i	
		s[si_opt][i].resize(num_job + 1);
		effe_mach[si_opt][i] = num_job > 0 ? true : false;	// the first element indicates the effectiveness of machine i
		q[si_opt][i].resize(num_job + 1);
		c[si_opt][i].resize(num_job + 1);
		for (int j = 1; j <= num_job; j++)
			ifs >> s[si_opt][i][j];	// the list of jobs to be assinged to machine i
	}
	ifs >> obj_given;
	for (int i = 1; i <= m; i++)
	{
		sort(s[si_opt][i].begin() + 1, s[si_opt][i].end(), cmpSort(r[i]));
		calculate_completion_time(si_opt, i);
	}
	calculate_obj(si_opt);
	if (fabs(obj_given - sol_obj[si_opt]) > MIN_EQUAL)
	{
		cout << "the given optimal solution is wrong."
			<< obj_given << " " << sol_obj[si_opt] << endl;
		//system("pause");
	}
	ifs.close();
}
PMS::~PMS()
{
	/*	for (int i = 0; i < sol_num; i++)
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
	delete[]sub_mach;

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
	delete[]r;*/
	ofs.close();
	//cout << "destruct function." << endl;
}
void PMS::display_problem()
{
	cout << ins_name << "n: " << n << ", m: " << m
		<< ", opt given: " << obj_given
		<< ", opt real: " << sol_obj[0] << endl;
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
	cout << "*** display solution: " << si
		<< ", obj: " /*<< fixed << setprecision(4)*/
		<< sol_obj[si] << ", " << mm[si] << " ***" << endl;
	for (int i = 1; i <= m; i++)
	{
		obj_type ir = 0;
		cout << i << ", " << effe_mach[si][i] << ", " << s[si][i].size() - 1 << ", " << c[si][i].back() << ": ";
		for (int j = 1; j < s[si][i].size(); j++)
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
		<< ", obj: " << sol_obj[si] << " ***" << endl;
	perf_type ql;
	obj_type cl, max_cl = 0, total_comp_t = 0;
	int max_cl_mach_index, sum_job_index = 0;
	for (int i = 1; i <= m; i++)
	{
		if (!effe_mach[si][i])
			continue;
		for (int j = 1; j < s[si][i].size(); j++)
			sum_job_index += s[si][i][j];
		ql = 1;
		cl = 0;
		for (int j = 1; j < s[si][i].size(); j++)
		{
			cl += p[i][s[si][i][j]] / ql;
			if (cl != c[si][i][j] || ql != q[si][i][j])
			{
				cout << "ERROR, si: " << si << " q or c is wrong "
					<< cl << ", " << c[si][i][j] << "\t"
					<< ql << ", " << q[si][i][j] << endl;
				system("pause");
			}
			ql *= (1 - d[i][s[si][i][j]]);
			if (j > 1 && r[i][s[si][i][j - 1]] - r[i][s[si][i][j]] < MIN_EQUAL)
			{
				cout << "ERROR, si: " << si << " not sort by r "
					<< si << "\t" << i << "\t" << j << "\t"
					<< r[i][s[si][i][j - 1]] << " " << r[i][s[si][i][j]] << endl;
				display_solution(si);
				system("pause");
			}
		}
		if (cl - max_cl > MIN_EQUAL)
		{
			max_cl = cl;
			max_cl_mach_index = i;
		}
		if (cl != c[si][i].back())
		{
			cout << "ERROR, completion time is wrong, si: " << si
				<< ", mach_index: " << i
				<< ", real completion time: " << cl << " " << c[si][i].back() << endl;
			system("pause");
		}
		total_comp_t += cl;
	}
	if (sum_job_index != (1 + n)*n / 2)
	{
		cout << "ERROR, si: " << si
			<< ", multiple job index" << endl;
		display_solution(si);
		system("pause");
	}
	if (total_comp_t != sol_obj[si])
	{
		cout << "ERROR, total completion time is wrong. " << si
			<< ", real total completion time: " << total_comp_t << " " << sol_obj[si] << endl;
	}
}
void PMS::check_solution_machine(int si, int mach_index)
{
	perf_type ql;
	obj_type cl, max_cl = 0;
	int max_cl_mach_index, sum_job_index = 0;
	int i = mach_index;
	for (int j = 1; j < s[si][i].size(); j++)
		sum_job_index += s[si][i][j];
	if (!effe_mach[si][i])
		return;
	ql = 1;
	cl = 0;
	for (int j = 1; j < s[si][i].size(); j++)
	{
		cl += p[i][s[si][i][j]] / ql;
		if (cl != c[si][i][j] || ql != q[si][i][j])
		{
			cout << "ERROR, si: " << si << " q or c is wrong "
				<< cl << ", " << c[si][i][j] << "\t"
				<< ql << ", " << q[si][i][j] << endl;
			system("pause");
		}
		ql *= (1 - d[i][s[si][i][j]]);
		if (j > 1 && r[i][s[si][i][j - 1]] - r[i][s[si][i][j]] < MIN_EQUAL)
		{
			cout << "ERROR, si: " << si << " not sort by r "
				<< si << "\t" << i << "\t" << j << "\t"
				<< r[i][s[si][i][j - 1]] << " " << r[i][s[si][i][j]] << endl;
			display_solution(si);
			system("pause");
		}
	}
	if (cl != c[si][i].back())
	{
		cout << "ERROR, completion time is wrong, si: " << si
			<< ", mach_index: " << i
			<< ", real completion time: " << cl << " " << c[si][i].back() << endl;
		system("pause");
	}

}
void PMS::calculate_completion_time(int si, int mach_index)
{
	for (int i = 1; i < s[si][mach_index].size(); i++)
	{
		q[si][mach_index][i] = (1 - d[mach_index][s[si][mach_index][i - 1]])*q[si][mach_index][i - 1];
		c[si][mach_index][i] = c[si][mach_index][i - 1] + p[mach_index][s[si][mach_index][i]] / q[si][mach_index][i];
	}
}
void PMS::calculate_obj(int si)
{
	sol_obj[si] = 0;
	for (int i = 1; i <= m; i++)
	{
		if (effe_mach[si][i] == false)
			continue;
		sol_obj[si] += c[si][i].back();
	}
}
void PMS::add_job(int si, int mach_index, int &position, int job)
{
	if (s[si][mach_index].size() == 1)	// empty machine
	{
		effe_mach[si][mach_index] = true;	// indicate it is effective
		s[si][mach_index].push_back(job);
		q[si][mach_index].push_back(1);
		c[si][mach_index].push_back(p[mach_index][job]);
	}
	else
	{
		if (position == 0)
		{
			position = 1;
			while (position < s[si][mach_index].size()
				&& r[mach_index][job] < r[mach_index][s[si][mach_index][position]])
				position += 1;
		}
		s[si][mach_index].insert(s[si][mach_index].begin() + position, job);
		q[si][mach_index].insert(q[si][mach_index].begin() + position,
			(1 - d[mach_index][s[si][mach_index][position - 1]])*q[si][mach_index][position - 1]);
		c[si][mach_index].insert(c[si][mach_index].begin() + position, c[si][mach_index][position - 1] +
			p[mach_index][s[si][mach_index][position]] / q[si][mach_index][position]);
		for (int i = position + 1; i < s[si][mach_index].size(); i++)
		{
			q[si][mach_index][i] = (1 - d[mach_index][s[si][mach_index][i - 1]])*q[si][mach_index][i - 1];
			c[si][mach_index][i] = c[si][mach_index][i - 1] + p[mach_index][s[si][mach_index][i]] / q[si][mach_index][i];
		}
	}
}
void PMS::remove_job(int si, int mach_index, int &position, int job)
{
	if (s[si][mach_index].empty())
	{
		cout << "ERROR, remove job from empty machine" << endl;
		system("pause");
	}
	if (position == 0)
	{
		position = 1;
		while (position < s[si][mach_index].size()
			&& s[si][mach_index][position] != job)
			position += 1;
	}
	else
	{
		if (s[si][mach_index][position] != job)
		{
			cout << "ERROR, remove the wrong job" << endl;
			system("pause");
		}
	}
	s[si][mach_index].erase(s[si][mach_index].begin() + position);
	q[si][mach_index].erase(q[si][mach_index].begin() + position);
	c[si][mach_index].erase(c[si][mach_index].begin() + position);
	for (int i = position; i < s[si][mach_index].size(); i++)
	{
		q[si][mach_index][i] = (1 - d[mach_index][s[si][mach_index][i - 1]])*q[si][mach_index][i - 1];
		c[si][mach_index][i] = c[si][mach_index][i - 1] + p[mach_index][s[si][mach_index][i]] / q[si][mach_index][i];
	}
	if (s[si][mach_index].size() == 1)
		effe_mach[si][mach_index] = false;	// indicate it is not effective
}
void PMS::exchange_job(int si, int mach_index, int &pos_add, int &pos_remove, int job_add, int job_remove)
{
	if (pos_remove == 0)
	{
		pos_remove = 1;
		while (pos_remove < s[si][mach_index].size()
			&& s[si][mach_index][pos_remove] != job_remove)
			pos_remove += 1;
	}
	else
	{
		if (s[si][mach_index][pos_remove] != job_remove)
		{
			cout << "ERROR, remove the wrong job in exchange job" << endl;
			system("pause");
		}
	}
	s[si][mach_index].erase(s[si][mach_index].begin() + pos_remove);
	if (pos_add == 0)
	{
		pos_add = 1;
		while (pos_add < s[si][mach_index].size()
			&& r[mach_index][job_add] < r[mach_index][s[si][mach_index][pos_add]])
			pos_add += 1;
	}
	s[si][mach_index].insert(s[si][mach_index].begin() + pos_add, job_add);
	for (int i = pos_add < pos_remove ? pos_add : pos_remove; i < s[si][mach_index].size(); i++)
	{
		q[si][mach_index][i] = (1 - d[mach_index][s[si][mach_index][i - 1]])*q[si][mach_index][i - 1];
		c[si][mach_index][i] = c[si][mach_index][i - 1] + p[mach_index][s[si][mach_index][i]] / q[si][mach_index][i];
	}
}
void PMS::replace_solution(int dest, int src)
{
	for (int i = 1; i <= m; i++)
	{
		s[dest][i].assign(s[src][i].begin(), s[src][i].end());
		q[dest][i].assign(q[src][i].begin(), q[src][i].end());
		c[dest][i].assign(c[src][i].begin(), c[src][i].end());
	}
	effe_mach[dest].assign(effe_mach[src].begin(), effe_mach[src].end());
	mm[dest] = mm[src];
	sol_obj[dest] = sol_obj[src];
}
void PMS::save_solution(int si, int si_opt, int iterration, int run_cnt)
{
	int result_improve, given_result_improve;
	if (fabs(makespan - sol_obj[si_opt]) <= MIN_EQUAL)
		result_improve = 1;	// equal
	else if (makespan - sol_obj[si] > MIN_EQUAL)
		result_improve = 2;	// improved
	else
		result_improve = 0;	// not improved
	if (fabs(makespan - obj_given) <= MIN_EQUAL)
		given_result_improve = 1;
	else if (obj_given - makespan > MIN_EQUAL)
		given_result_improve = 2;
	else
		given_result_improve = 0;
	ofs << run_cnt << "\t"
		//<< c[si_opt][0] << "\t"
		<< sol_obj[si] << "\t"
		<< makespan << "\t"
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
			ofs << c[si][i].back() << "\t"
				<< s[si][i].size() - 1 << "\t";
			for (int j = 1; j < s[si][i].size(); j++)
				ofs << s[si][i][j] << "\t";
			ofs << endl;
		}
	}
}
void PMS::init_solution(int si, R_Mode rm)
{
	if (rm == R_Mode::RANDOM)
	{
		for (int j = 1; j <= n; j++)
			charact[rm][j] = rand() % 1000;
	}
	vector<int> job_sort_charact_vec(n + 1);
	for (int i = 1; i <= n; i++)
		job_sort_charact_vec[i] = i;	// the machine 0 as the temp memory
	sort(job_sort_charact_vec.begin() + 1, job_sort_charact_vec.end(), cmpSort(charact[rm], Cmp_Mode::GREATER));
	/*for (int i = 1; i <= n; i++)
	cout << job_sort_charact_vec[i] << " " << charact[rm][job_sort_charact_vec[i]] << "\t";
	cout << endl;*/
	for (int i = 1; i <= m; i++)
	{
		s[si][i].resize(1, 0);	// indicate it is effective
		q[si][i].resize(1, 1);	// out of essence
		c[si][i].resize(1, 0);	// out of essence
	}
	effe_mach[si].assign(m + 1, true);
	mm[si] = 0;
	for (int j = 1; j <= n; j++)
	{
		proc_type min_completion_time = DBL_MAX;
		int min_ct_mach_index, min_position;
		for (int i = 1; i <= m; i++)
		{
			int position = 0;
			add_job(si, i, position, job_sort_charact_vec[j]);
			if (c[si][i].back() < min_completion_time)
			{
				min_completion_time = c[si][i].back();
				min_ct_mach_index = i;
				min_position = position;
			}
			remove_job(si, i, position, job_sort_charact_vec[j]);
		}
		add_job(si, min_ct_mach_index, min_position, job_sort_charact_vec[j]);
	}
	calculate_obj(si);
	//check_solution(si);
}
void PMS::perturb(int si, int ptr_rate)
{
	for (int i = 0; i < (s[si][mm[si]].size() - 1) * ptr_rate*0.01; i++)
	{
		int pos_remove_mm = rand() % (s[si][mm[si]].size() - 1) + 1, pos_add_mm = 0;
		int job_add = s[si][mm[si]][pos_remove_mm];
		int mach_i = rand() % m + 1;
		while (mach_i == mm[si] || s[si][mach_i].size() == 1)
			mach_i = rand() % m + 1;
		int pos_remove_mi = rand() % (s[si][mach_i].size() - 1) + 1, pos_add_mi = 0;
		int job_remove = s[si][mach_i][pos_remove_mi];
		exchange_job(si, mach_i, pos_add_mi, pos_remove_mi, job_add, job_remove);
		exchange_job(si, mm[si], pos_add_mm, pos_remove_mm, job_remove, job_add);
	}
	calculate_obj(si);
}
void PMS::perturb1(int si, int ptr_rate)
{
	int mach1, mach2, job1, job2, effe_num = 0;
	bool whe_swap = true;
	for (int i = 1; i <= m; i++)
	{
		if (effe_mach[si][i])
			effe_num += 1;
	}
	if (effe_num <= 1)
		whe_swap = false;
	for (int i = 0; i < n*ptr_rate*0.01&&whe_swap; i++)
	{
		mach1 = rand() % m + 1;
		while (s[si][mach1].size() == 1)
			mach1 = rand() % m + 1;
		mach2 = rand() % m + 1;
		while (mach1 == mach2 || s[si][mach2].size() == 1)
			mach2 = rand() % m + 1;
		job1 = rand() % (s[si][mach1].size() - 1) + 1;
		job2 = rand() % (s[si][mach2].size() - 1) + 1;
		int temp = s[si][mach1][job1];
		s[si][mach1][job1] = s[si][mach2][job2];
		s[si][mach2][job2] = temp;
	}
	for (int k = 1; k <= m + m*ptr_rate*0.01; k++)
	{
		int i = rand() % m + 1;
		while (s[si][i].size() <= 2)
			i = rand() % m + 1;
		int pos_remove = rand() % (s[si][i].size() - 1) + 1;
		int mach = rand() % m + 1;
		while (mach == i)
			mach = rand() % m + 1;
		s[si][mach].push_back(s[si][i][pos_remove]);
		if (s[si][mach].size() > 1)
			effe_mach[si][mach] = true;
		s[si][i].erase(s[si][i].begin() + pos_remove);
	}
	for (int i = 1; i <= m; i++)
	{
		q[si][i].resize(s[si][i].size());
		c[si][i].resize(s[si][i].size());
	}
	/*for (int i = 0; i < (s[si][mm[si]].size() - 1) * ptr_rate*0.01; i++)
	{
		int jm = rand() % (s[si][mm[si]].size() - 1) + 1;
		int mach_r = rand() % m + 1;
		while (mach_r == mm[si] || s[si][mach_r].size() == 1)
			mach_r = rand() % m + 1;
		int jo = rand() % (s[si][mach_r].size() - 1) + 1;
		int temp = s[si][mm[si]][jm];
		s[si][mm[si]][jm] = s[si][mach_r][jo];
		s[si][mach_r][jo] = temp;
	}*/
	for (int i = 1; i <= m; i++)
	{
		sort(s[si][i].begin() + 1, s[si][i].end(), cmpSort(r[i]));
		calculate_completion_time(si, i);
	}
	calculate_obj(si);
}
void PMS::local_search_hybrid(int si)
{
	bool is_still_improved = true;
	while (is_still_improved)
	{
		is_still_improved = false;
		int pre_mm = mm[si];
		for (int jm = 0 % (s[si][mm[si]].size() - 1) + 1, jm_r = 1;
			jm_r < s[si][mm[si]].size();
			jm_r++, jm = jm % (s[si][mm[si]].size() - 1) + 1)
		{
			for (int i = 0 % m + 1, i_r = 1; i_r <= m; i_r++, i = i%m + 1)
			{
				if (i == mm[si] || !effe_mach[si][i])
					continue;
				//display_solution(si);
				int pos_remove = 0, pos_add = 0, cur_job = s[si][mm[si]][jm];
				obj_type pre_obj = sol_obj[si];
				remove_job(si, mm[si], pos_remove, cur_job);
				add_job(si, i, pos_add, cur_job);
				calculate_obj(si);
				//check_solution(si_local);
				if (pre_obj - sol_obj[si] > MIN_EQUAL)
				{
					is_still_improved = true;
					break;
				}
				else
				{
					remove_job(si, i, pos_add, cur_job);
					add_job(si, pre_mm, pos_remove, cur_job);
					//calculate_obj(si);
					sol_obj[si] = pre_obj;
				}
				/*display_solution(si);
				check_solution(si);*/
				for (int j = rand() % ((s[si][i].size() - 1) == 0 ? 1 : (s[si][i].size() - 1)) + 1, j_r = 1;
					j_r < s[si][i].size(); j_r++, j = j % (s[si][i].size() - 1) + 1)
				{
					int pos_remove_mi = j, pos_add_mi = 0, pos_remove_mm = jm, pos_add_mm = 0,
						job_remove = s[si][i][j], job_add = s[si][mm[si]][jm];
					obj_type pre_obj = sol_obj[si];
					/*if (i == mm[si])
					{
					system("pause");
					}*/
					exchange_job(si, i, pos_add_mi, pos_remove_mi, job_add, job_remove);
					exchange_job(si, mm[si], pos_add_mm, pos_remove_mm, job_remove, job_add);
					calculate_obj(si);
					//check_solution(si);
					/*display_solution(si);*/
					if (pre_obj - sol_obj[si] > MIN_EQUAL)
					{
						is_still_improved = true;
						break;
					}
					else
					{
						exchange_job(si, i, pos_remove_mi, pos_add_mi, job_remove, job_add);
						exchange_job(si, pre_mm, pos_remove_mm, pos_add_mm, job_add, job_remove);
						//calculate_obj(si);
						sol_obj[si] = pre_obj;
						mm[si] = pre_mm;
					}
					/*check_solution(si);
					display_solution(si);*/
				}
			}
		}
	}
}

void PMS::local_search(int si)
{
	int num_imp = 0, num_ins_swap = 0;
	bool is_still_improve = true;
	while (is_still_improve)
	{
		num_imp += 1;
		is_still_improve = false;
		for (int i1 = 1; i1 <= m; i1++)
		{
			if (!effe_mach[si][i1])
				continue;
			for (int j1 = 1; j1 < s[si][i1].size(); j1++)
			{
				int pos_remove_j1 = j1, job_i1 = s[si][i1][j1];
				bool is_improve_swap = true;
				for (int i2 = 1; i2 <= m&&is_improve_swap; i2++)
				{
					if (i1 == i2 || !effe_mach[si][i2])
						continue;
					int pos_add_j2 = 0;
					obj_type pre_obj = sol_obj[si];
					remove_job(si, i1, pos_remove_j1, job_i1);
					add_job(si, i2, pos_add_j2, job_i1);
					calculate_obj(si);
					//check_solution(si);
					if (pre_obj - sol_obj[si] > MIN_EQUAL)
					{
						is_still_improve = true;
						num_ins_swap += 1;
						break;
					}
					else
					{
						remove_job(si, i2, pos_add_j2, job_i1);
						add_job(si, i1, pos_remove_j1, job_i1);
						sol_obj[si] = pre_obj;
					}
					for (int j2 = 1; j2 < s[si][i2].size(); j2++)
					{
						int pos_add_j1 = 0, pos_remove_j2 = j2, job_i2 = s[si][i2][j2];
						pos_add_j2 = 0;	// consider to remove this line
						pre_obj = sol_obj[si];
						exchange_job(si, i1, pos_add_j1, pos_remove_j1, job_i2, job_i1);
						exchange_job(si, i2, pos_add_j2, pos_remove_j2, job_i1, job_i2);
						calculate_obj(si);
						if (pre_obj - sol_obj[si] > MIN_EQUAL)
						{
							is_still_improve = true;
							is_improve_swap = false;
							num_ins_swap += 1;
							break;
						}
						else
						{
							exchange_job(si, i1, pos_remove_j1, pos_add_j1, job_i1, job_i2);
							exchange_job(si, i2, pos_remove_j2, pos_add_j2, job_i2, job_i1);
							//calculate_obj(si);
							sol_obj[si] = pre_obj;
						}
					}
				}
			}
		}
	}
	/*cout << "num_cmp: " << num_imp
		<< ", num_ins_swap: " << num_ins_swap << endl;*/
		//check_solution(si);
}
void PMS::local_search_hybrid_tri_insert(int si)
{
	bool is_still_improved = true;
	while (is_still_improved)
	{
		is_still_improved = false;
		bool is_mm_unchanged = true;
		int pre_mm = mm[si];
		for (int jm = 0 % (s[si][mm[si]].size() - 1) + 1, jm_r = 1;
			jm_r < s[si][mm[si]].size() && is_mm_unchanged;
			jm_r++, jm = jm % (s[si][mm[si]].size() - 1) + 1)
		{
			bool is_still_improved_swap = true;
			for (int i = 0 % m + 1, i_r = 1; i_r <= m && is_mm_unchanged&&is_still_improved_swap; i_r++, i = i%m + 1)
			{
				//is_still_improved_swap = false;
				if (i == mm[si] || !effe_mach[si][i])
					continue;
				//display_solution(si);
				int pos_remove = 0, pos_add = 0, cur_job = s[si][mm[si]][jm];
				obj_type pre_obj = sol_obj[si];
				remove_job(si, mm[si], pos_remove, cur_job);
				add_job(si, i, pos_add, cur_job);
				calculate_obj(si);
				//check_solution(si_local);
				if (pre_obj - sol_obj[si] > MIN_EQUAL)
				{
					if (pre_mm != mm[si])
						is_mm_unchanged = false;
					is_still_improved = true;
					break;
				}
				else
				{
					remove_job(si, i, pos_add, cur_job);
					add_job(si, pre_mm, pos_remove, cur_job);
					//calculate_obj(si);
					sol_obj[si] = pre_obj;
					mm[si] = pre_mm;
				}
				/*display_solution(si);
				check_solution(si);*/
				for (int j = rand() % ((s[si][i].size() - 1) == 0 ? 1 : (s[si][i].size() - 1)) + 1, j_r = 1;
					j_r < s[si][i].size() && is_mm_unchanged; j_r++, j = j % (s[si][i].size() - 1) + 1)
				{
					int pos_remove_mi = j, pos_add_mi = 0, pos_remove_mm = jm, pos_add_mm = 0,
						job_remove = s[si][i][j], job_add = s[si][mm[si]][jm];
					obj_type pre_obj = sol_obj[si];
					/*if (i == mm[si])
					{
					system("pause");
					}*/
					exchange_job(si, i, pos_add_mi, pos_remove_mi, job_add, job_remove);
					exchange_job(si, mm[si], pos_add_mm, pos_remove_mm, job_remove, job_add);
					calculate_obj(si);
					//check_solution(si);
					/*display_solution(si);*/
					if (pre_obj - sol_obj[si] > MIN_EQUAL)
					{
						if (pre_mm != mm[si])
							is_mm_unchanged = false;
						is_still_improved_swap = false;
						is_still_improved = true;
						break;
					}
					else
					{
						exchange_job(si, i, pos_remove_mi, pos_add_mi, job_remove, job_add);
						exchange_job(si, pre_mm, pos_remove_mm, pos_add_mm, job_add, job_remove);
						//calculate_obj(si);
						sol_obj[si] = pre_obj;
						mm[si] = pre_mm;
					}
					/*check_solution(si);
					display_solution(si);*/
				}
			}
		}
		if (is_still_improved == false)
		{
			bool is_triple_insert = true;
			for (int i1 = 1; i1 <= m&&is_triple_insert; i1++)
			{
				if (i1 == mm[si] || !effe_mach[si][i1])
					continue;
				for (int j1 = 1; j1 < s[si][i1].size() && is_triple_insert; j1++)
				{
					for (int i2 = i1 + 1; i2 <= m&&is_triple_insert; i2++)
					{
						if (i2 == mm[si] || !effe_mach[si][i2])
							continue;
						int pos_remove = 0, pos_add_i2 = 0, cur_job1 = s[si][i1][j1];
						add_job(si, i2, pos_add_i2, cur_job1);	// try to add cur_job1 from i1 to i2
																//check_solution_machine(si, i2);
						if (sol_obj[si] - c[si][i2].back() > MIN_EQUAL)
						{
							remove_job(si, i1, pos_remove, cur_job1);	// actually remove cur_job1 from i1 
							for (int jm = 1; jm < s[si][mm[si]].size() && is_triple_insert; jm++)
							{
								int pos_add_i1 = 0;
								int job_mm = s[si][mm[si]][jm];
								add_job(si, i1, pos_add_i1, job_mm);	// add job mm to i1
																		//check_solution_machine(si, i1);
								if (sol_obj[si] - c[si][i1].back() > MIN_EQUAL)
								{
									int pos_remove_mm = 0;
									remove_job(si, mm[si], pos_remove_mm, job_mm);	// remove job mm from mm
									calculate_obj(si);
									//check_solution(si);
									is_still_improved = true;
									is_triple_insert = false;
									break;
								}
								remove_job(si, i1, pos_add_i1, job_mm);	// remove job mm from i1
							}
							if (is_still_improved == false)
							{
								add_job(si, i1, pos_remove, cur_job1);	// repair i1 by adding cur_job1 to i1
																		//check_solution_machine(si, i1);
							}
						}
						if (is_still_improved == false)
							remove_job(si, i2, pos_add_i2, cur_job1);
					}
				}
			}
		}
	}
}
void PMS::local_search_hybrid_tri_swap(int si)
{
	bool is_still_improved = true;
	while (is_still_improved)
	{
		is_still_improved = false;
		bool is_mm_unchanged = true;
		int pre_mm = mm[si];
		for (int jm = 0 % (s[si][mm[si]].size() - 1) + 1, jm_r = 1;
			jm_r < s[si][mm[si]].size() && is_mm_unchanged;
			jm_r++, jm = jm % (s[si][mm[si]].size() - 1) + 1)
		{
			bool is_still_improved_swap = true;
			for (int i = 0 % m + 1, i_r = 1; i_r <= m && is_mm_unchanged&&is_still_improved_swap; i_r++, i = i%m + 1)
			{
				//is_still_improved_swap = false;
				if (i == mm[si] || !effe_mach[si][i])
					continue;
				//display_solution(si);
				int pos_remove = 0, pos_add = 0, cur_job = s[si][mm[si]][jm];
				obj_type pre_obj = sol_obj[si];
				remove_job(si, mm[si], pos_remove, cur_job);
				add_job(si, i, pos_add, cur_job);
				calculate_obj(si);
				//check_solution(si_local);
				if (pre_obj - sol_obj[si] > MIN_EQUAL)
				{
					if (pre_mm != mm[si])
						is_mm_unchanged = false;
					is_still_improved = true;
					break;
				}
				else
				{
					remove_job(si, i, pos_add, cur_job);
					add_job(si, pre_mm, pos_remove, cur_job);
					//calculate_obj(si);
					sol_obj[si] = pre_obj;
					mm[si] = pre_mm;
				}
				/*display_solution(si);
				check_solution(si);*/
				for (int j = rand() % ((s[si][i].size() - 1) == 0 ? 1 : (s[si][i].size() - 1)) + 1, j_r = 1;
					j_r < s[si][i].size() && is_mm_unchanged; j_r++, j = j % (s[si][i].size() - 1) + 1)
				{
					int pos_remove_mi = j, pos_add_mi = 0, pos_remove_mm = jm, pos_add_mm = 0,
						job_remove = s[si][i][j], job_add = s[si][mm[si]][jm];
					obj_type pre_obj = sol_obj[si];
					/*if (i == mm[si])
					{
					system("pause");
					}*/
					exchange_job(si, i, pos_add_mi, pos_remove_mi, job_add, job_remove);
					exchange_job(si, mm[si], pos_add_mm, pos_remove_mm, job_remove, job_add);
					calculate_obj(si);
					//check_solution(si);
					/*display_solution(si);*/
					if (pre_obj - sol_obj[si] > MIN_EQUAL)
					{
						if (pre_mm != mm[si])
							is_mm_unchanged = false;
						is_still_improved_swap = false;
						is_still_improved = true;
						break;
					}
					else
					{
						exchange_job(si, i, pos_remove_mi, pos_add_mi, job_remove, job_add);
						exchange_job(si, pre_mm, pos_remove_mm, pos_add_mm, job_add, job_remove);
						//calculate_obj(si);
						sol_obj[si] = pre_obj;
						mm[si] = pre_mm;
					}
					/*check_solution(si);
					display_solution(si);*/
				}
			}
		}
		if (is_still_improved == false)
		{
			bool is_triple_swap = true;
			for (int i1 = 1; i1 <= m&&is_triple_swap; i1++)
			{
				if (i1 == mm[si] || !effe_mach[si][i1])
					continue;
				for (int j1 = 1; j1 < s[si][i1].size() && is_triple_swap; j1++)
				{
					for (int i2 = i1 + 1; i2 <= m&&is_triple_swap; i2++)
					{
						if (i2 == mm[si] || !effe_mach[si][i2])
							continue;
						for (int j2 = 1; j2 < s[si][i2].size() && is_triple_swap; j2++)
						{
							int job_i1 = s[si][i1][j1], job_i2 = s[si][i2][j2];
							int pos_remove_i1 = j1, pos_remove_i2 = j2, pos_add_i1 = 0, pos_add_i2 = 0;
							exchange_job(si, i1, pos_add_i1, pos_remove_i1, job_i2, job_i1);
							if (sol_obj[si] - c[si][i1].back() > MIN_EQUAL)
							{
								exchange_job(si, i2, pos_add_i2, pos_remove_i2, job_i1, job_i2);
								if (sol_obj[si] - c[si][i2].back() > MIN_EQUAL)
								{
									for (int jm = 1; jm < s[si][mm[si]].size() && is_triple_swap; jm++)
									{
										int job_mm = s[si][mm[si]][jm], pos_add_i = 0;
										add_job(si, i1, pos_add_i, job_mm);	// try to add job_mm to i1
										if (sol_obj[si] - c[si][i1].back() > MIN_EQUAL)
										{
											int pos_remove_mm = 0;
											remove_job(si, mm[si], pos_remove_mm, job_mm);	// actually remove job_mm from mm
											calculate_obj(si);
											is_still_improved = true;
											is_triple_swap = false;
											break;
										}
										remove_job(si, i1, pos_add_i, job_mm);	// can not improve, remove job_mm from i1
										pos_add_i = 0;
										add_job(si, i2, pos_add_i, job_mm);	// try to add job_mm to i2
										if (sol_obj[si] - c[si][i2].back() > MIN_EQUAL)
										{
											int pos_remove_mm = 0;
											remove_job(si, mm[si], pos_remove_mm, job_mm);	// actually remove job_mm from mm
											calculate_obj(si);
											is_still_improved = true;
											is_triple_swap = false;
											break;
										}
										remove_job(si, i2, pos_add_i, job_mm);	// can not improve, remove job_mm from i2
									}
								}
								if (is_still_improved == false)
									exchange_job(si, i2, pos_remove_i2, pos_add_i2, job_i2, job_i1);
							}
							if (is_still_improved == false)
								exchange_job(si, i1, pos_remove_i1, pos_add_i1, job_i1, job_i2);
						}
					}
				}
			}
		}
	}
}
void PMS::local_search_hybrid_tri_insert_swap(int si)
{
	bool is_still_improved = true;
	while (is_still_improved)
	{
		is_still_improved = false;
		bool is_mm_unchanged = true;
		int pre_mm = mm[si];
		for (int jm = 0 % (s[si][mm[si]].size() - 1) + 1, jm_r = 1;
			jm_r < s[si][mm[si]].size() && is_mm_unchanged;
			jm_r++, jm = jm % (s[si][mm[si]].size() - 1) + 1)
		{
			bool is_still_improved_swap = true;
			for (int i = 0 % m + 1, i_r = 1; i_r <= m && is_mm_unchanged&&is_still_improved_swap; i_r++, i = i%m + 1)
			{
				//is_still_improved_swap = false;
				if (i == mm[si] || !effe_mach[si][i])
					continue;
				//display_solution(si);
				int pos_remove = 0, pos_add = 0, cur_job = s[si][mm[si]][jm];
				obj_type pre_obj = sol_obj[si];
				remove_job(si, mm[si], pos_remove, cur_job);
				add_job(si, i, pos_add, cur_job);
				calculate_obj(si);
				//check_solution(si_local);
				if (pre_obj - sol_obj[si] > MIN_EQUAL)
				{
					if (pre_mm != mm[si])
						is_mm_unchanged = false;
					is_still_improved = true;
					break;
				}
				else
				{
					remove_job(si, i, pos_add, cur_job);
					add_job(si, pre_mm, pos_remove, cur_job);
					//calculate_obj(si);
					sol_obj[si] = pre_obj;
					mm[si] = pre_mm;
				}
				/*display_solution(si);
				check_solution(si);*/
				for (int j = rand() % ((s[si][i].size() - 1) == 0 ? 1 : (s[si][i].size() - 1)) + 1, j_r = 1;
					j_r < s[si][i].size() && is_mm_unchanged; j_r++, j = j % (s[si][i].size() - 1) + 1)
				{
					int pos_remove_mi = j, pos_add_mi = 0, pos_remove_mm = jm, pos_add_mm = 0,
						job_remove = s[si][i][j], job_add = s[si][mm[si]][jm];
					obj_type pre_obj = sol_obj[si];
					/*if (i == mm[si])
					{
					system("pause");
					}*/
					exchange_job(si, i, pos_add_mi, pos_remove_mi, job_add, job_remove);
					exchange_job(si, mm[si], pos_add_mm, pos_remove_mm, job_remove, job_add);
					calculate_obj(si);
					//check_solution(si);
					/*display_solution(si);*/
					if (pre_obj - sol_obj[si] > MIN_EQUAL)
					{
						if (pre_mm != mm[si])
							is_mm_unchanged = false;
						is_still_improved_swap = false;
						is_still_improved = true;
						break;
					}
					else
					{
						exchange_job(si, i, pos_remove_mi, pos_add_mi, job_remove, job_add);
						exchange_job(si, pre_mm, pos_remove_mm, pos_add_mm, job_add, job_remove);
						//calculate_obj(si);
						sol_obj[si] = pre_obj;
						mm[si] = pre_mm;
					}
					/*check_solution(si);
					display_solution(si);*/
				}
			}
		}
		if (is_still_improved == false)
		{
			bool is_triple_insert = true;
			for (int i1 = 1; i1 <= m&&is_triple_insert; i1++)
			{
				if (i1 == mm[si] || !effe_mach[si][i1])
					continue;
				for (int j1 = 1; j1 < s[si][i1].size() && is_triple_insert; j1++)
				{
					for (int i2 = i1 + 1; i2 <= m&&is_triple_insert; i2++)
					{
						if (i2 == mm[si] || !effe_mach[si][i2])
							continue;
						int pos_remove = 0, pos_add_i2 = 0, cur_job1 = s[si][i1][j1];
						add_job(si, i2, pos_add_i2, cur_job1);	// try to add cur_job1 from i1 to i2
																//check_solution_machine(si, i2);
						if (sol_obj[si] - c[si][i2].back() > MIN_EQUAL)
						{
							remove_job(si, i1, pos_remove, cur_job1);	// actually remove cur_job1 from i1 
							for (int jm = 1; jm < s[si][mm[si]].size() && is_triple_insert; jm++)
							{
								int pos_add_i1 = 0;
								int job_mm = s[si][mm[si]][jm];
								add_job(si, i1, pos_add_i1, job_mm);	// add job mm to i1
																		//check_solution_machine(si, i1);
								if (sol_obj[si] - c[si][i1].back() > MIN_EQUAL)
								{
									int pos_remove_mm = 0;
									remove_job(si, mm[si], pos_remove_mm, job_mm);	// remove job mm from mm
									calculate_obj(si);
									//check_solution(si);
									is_still_improved = true;
									is_triple_insert = false;
									break;
								}
								remove_job(si, i1, pos_add_i1, job_mm);	// remove job mm from i1
							}
							if (is_still_improved == false)
							{
								add_job(si, i1, pos_remove, cur_job1);	// repair i1 by adding cur_job1 to i1
																		//check_solution_machine(si, i1);
							}
						}
						if (is_still_improved == false)
							remove_job(si, i2, pos_add_i2, cur_job1);
					}
				}
			}
		}
		if (is_still_improved == false)
		{
			bool is_triple_swap = true;
			for (int i1 = 1; i1 <= m&&is_triple_swap; i1++)
			{
				if (i1 == mm[si] || !effe_mach[si][i1])
					continue;
				for (int j1 = 1; j1 < s[si][i1].size() && is_triple_swap; j1++)
				{
					for (int i2 = i1 + 1; i2 <= m&&is_triple_swap; i2++)
					{
						if (i2 == mm[si] || !effe_mach[si][i2])
							continue;
						for (int j2 = 1; j2 < s[si][i2].size() && is_triple_swap; j2++)
						{
							int job_i1 = s[si][i1][j1], job_i2 = s[si][i2][j2];
							int pos_remove_i1 = j1, pos_remove_i2 = j2, pos_add_i1 = 0, pos_add_i2 = 0;
							exchange_job(si, i1, pos_add_i1, pos_remove_i1, job_i2, job_i1);
							if (sol_obj[si] - c[si][i1].back() > MIN_EQUAL)
							{
								exchange_job(si, i2, pos_add_i2, pos_remove_i2, job_i1, job_i2);
								if (sol_obj[si] - c[si][i2].back() > MIN_EQUAL)
								{
									for (int jm = 1; jm < s[si][mm[si]].size() && is_triple_swap; jm++)
									{
										int job_mm = s[si][mm[si]][jm], pos_add_i = 0;
										add_job(si, i1, pos_add_i, job_mm);	// try to add job_mm to i1
										if (sol_obj[si] - c[si][i1].back() > MIN_EQUAL)
										{
											int pos_remove_mm = 0;
											remove_job(si, mm[si], pos_remove_mm, job_mm);	// actually remove job_mm from mm
											calculate_obj(si);
											is_still_improved = true;
											is_triple_swap = false;
											break;
										}
										remove_job(si, i1, pos_add_i, job_mm);	// can not improve, remove job_mm from i1
										pos_add_i = 0;
										add_job(si, i2, pos_add_i, job_mm);	// try to add job_mm to i2
										if (sol_obj[si] - c[si][i2].back() > MIN_EQUAL)
										{
											int pos_remove_mm = 0;
											remove_job(si, mm[si], pos_remove_mm, job_mm);	// actually remove job_mm from mm
											calculate_obj(si);
											is_still_improved = true;
											is_triple_swap = false;
											break;
										}
										remove_job(si, i2, pos_add_i, job_mm);	// can not improve, remove job_mm from i2
									}
								}
								if (is_still_improved == false)
									exchange_job(si, i2, pos_remove_i2, pos_add_i2, job_i2, job_i1);
							}
							if (is_still_improved == false)
								exchange_job(si, i1, pos_remove_i1, pos_add_i1, job_i1, job_i2);
						}
					}
				}
			}
		}
	}
}

void PMS::divide_and_conquer(int si, int mach_num)
{
	if (mach_num == 2)
	{
		if (ls_method == 0)
			local_search_hybrid(si);
		else if (ls_method == 1)
			local_search_hybrid_tri_insert(si);
		else if (ls_method == 2)
			local_search_hybrid_tri_swap(si);
		else if (ls_method == 3)
			//local_search_hybrid_tri_insert_swap(si);
			local_search(si);
		//ls_cnt1 += 1;
		/*cout << "exe 2: ";
		for (int i = 1; i <= m; i++)
		{
		if (effe_mach[si][i])
		cout << i << "\t";
		}
		cout << endl;*/
		return;
	}
	else if (mach_num == 1)
	{
		//ls_cnt1 += 1;
		/*cout << "exe 1: ";
		for (int i = 1; i <= m; i++)
		{
		if (effe_mach[si][i])
		cout << i << "\t";
		}
		cout << endl;*/
		return;
	}
	vector<int> sub_mach1, sub_mach2;
	for (int i = 1; i <= m; i++)
	{
		if (effe_mach[si][i])
			sub_mach1.push_back(i);
	}
	int half = sub_mach1.size() / 2;
	for (int i = 0; i < half; i++)
	{
		/*sub_mach2.push_back(sub_mach1.back());
		sub_mach1.pop_back();*/
		int rd = rand() % sub_mach1.size();
		sub_mach2.push_back(sub_mach1[rd]);
		sub_mach1.erase(sub_mach1.begin() + rd);
	}
	sort(sub_mach2.begin(), sub_mach2.end());
	/*cout << "before sub1: ";
	for (vector<int>::iterator iter = sub_mach1.begin();
	iter != sub_mach1.end(); iter++)
	{
	cout << *iter << "\t";
	}
	cout << ", sub2: ";
	for (vector<int>::iterator iter = sub_mach2.begin();
	iter != sub_mach2.end(); iter++)
	{
	cout << *iter << "\t";
	}
	cout << endl;*/

	sol_obj[si] = 0;
	effe_mach[si].assign(m + 1, false);
	for (int i = 0; i < sub_mach1.size(); i++)
	{
		sol_obj[si] += c[si][sub_mach1[i]].back();
	}
	divide_and_conquer(si, sub_mach1.size());

	sol_obj[si] = 0;
	effe_mach[si].assign(m + 1, false);
	for (int i = 0; i < sub_mach2.size(); i++)
	{
		sol_obj[si] += c[si][sub_mach2[i]].back();
	}
	divide_and_conquer(si, sub_mach2.size());

	effe_mach[si].assign(m + 1, false);
	//cout << "after sub1: ";
	for (vector<int>::iterator iter = sub_mach1.begin();
		iter != sub_mach1.end(); iter++)
	{
		effe_mach[si][*iter] = true;
		//cout << *iter << "\t";
	}
	//cout << ", sub2: ";
	for (vector<int>::iterator iter = sub_mach2.begin();
		iter != sub_mach2.end(); iter++)
	{
		effe_mach[si][*iter] = true;
		//cout << *iter << "\t";
	}
	//cout << endl;
	sol_obj[si] = 0;
	for (int i = 1; i <= m; i++)
	{
		if (effe_mach[si][i] == false)
			continue;
		sol_obj[si] += c[si][i].back();
	}
	if (ls_method == 0)
		local_search_hybrid(si);
	else if (ls_method == 1)
		local_search_hybrid_tri_insert(si);
	else if (ls_method == 2)
		local_search_hybrid_tri_swap(si);
	else if (ls_method == 3)
		//local_search_hybrid_tri_insert_swap(si);
		local_search(si);
	ls_cnt2 += 1;
	//check_solution(si);
}
void PMS::iterated_local_search(int iteration, int perturb_rate, R_Mode r_mode, NS_Mode ns, int run_cnt_from, int run_cnt_to)
{
	int si_opt = 0, si_best = 1,
		si_cur = 2, si_local = 3, si_ptr = 4;
	int opt_cnt = 0, imp_cnt = 0, non_imp_cnt = 0;
	obj_type sum_delta_obj = 0;
	int rt = time(NULL);
	//rt = 1461506514;
	srand(rt);
	ofs << ins_name << "\t" << n << "\t" << m << "\t"
		<< obj_given << "\t" << sol_obj[si_opt] << "\t" << rt << endl;
	cout << ins_name << "\t" << n << "\t" << m << "\t"
		<< obj_given << "\t" << sol_obj[si_opt] << "\t" << rt << endl;
	double temperature0 = temperature;
	for (int rc = run_cnt_from; rc <= run_cnt_to; rc++)
	{
		start_tm = clock();
		init_solution(si_cur, r_mode);	//PMS::MINPDD
		//display_solution(si_cur);
		local_search(si_cur);
		//local_search_hybrid(si_cur);
		//local_search_ejection_chain(si_cur, si_local, ns);
		replace_solution(si_best, si_cur);
		end_tm = clock();
		//return;
		//cout << rc << "\t0\t" << sol_obj[si_best] << endl;
		//save_solution(si_best, si_opt, 0, run_cnt);
		int min_obj_iter = 0;
		//if (n <= 14)	// find the optimal solution for instances in OB set
		//{
		//	if (fabs(sol_obj[si_best] - sol_obj[si_opt]) <= MIN_EQUAL)
		//		iteration = 0;
		//	else
		//		iteration = INT16_MAX;
		//}
		for (int i = 0; i < iteration; i++)
		{
			replace_solution(si_ptr, si_cur);
			//display_solution(si_ptr);
			//check_solution(si_ptr);
			perturb1(si_ptr, perturb_rate);
			//display_solution(si_ptr);
			//cout << sol_obj[si_best]<< ", " << sol_obj[si_cur] << endl;
			//local_search_ejection_chain(si_ptr, si_local, ns);
			//check_solution(si_ptr);
			//local_search_hybrid(si_ptr);
			local_search(si_ptr);
			//check_solution(si_ptr);
			if (sol_obj[si_best] - sol_obj[si_ptr] > MIN_EQUAL)
			{
				replace_solution(si_best, si_ptr);
				end_tm = clock();
				min_obj_iter = i;
				//cout << rc << "\t" << i<<"\t"<<sol_obj[si_best] << endl;
				if (n <= 14 && fabs(sol_obj[si_best] - sol_obj[si_opt]) <= MIN_EQUAL)
					break;
			}
			int accept = 0;
			int rand_num = rand() % 100;
			obj_type obj_si_cur = sol_obj[si_cur];
			//if (sol_obj[si_cur] - sol_obj[si_ptr] > MIN_EQUAL ||
			//	rand() % 100 <= (control_para*100))
			//{
			//	replace_solution(si_cur, si_ptr);
			//	//accept = 1;
			//}
			if (sol_obj[si_cur] - sol_obj[si_ptr] > MIN_EQUAL ||
				rand_num <= (100 * exp((sol_obj[si_cur] - sol_obj[si_ptr]) / temperature)))
			{
				replace_solution(si_cur, si_ptr);
				accept = 1;
			}
			//replace_solution(si_cur, si_best);
			/*cout << temperature << "\t"
			<< obj_si_cur <<"\t"
			<<sol_obj[si_ptr]<<"\t"
			<<rand_num<<"\t"
			<< (100 * exp((obj_si_cur - sol_obj[si_ptr]) / (1*temperature)))<<"\t"
			<<accept
			<< endl;*/
			temperature *= control_para;
			if (temperature < 1)
				temperature = temperature0;
		}
		//check_solution(si_best); 
		//display_solution(si_best);
		makespan = 0;
		for (int i = 1; i <= m; i++)
		{
			if (!effe_mach[si_best][i])
				continue;
			if (c[si_best][i].back() - makespan > MIN_EQUAL)
			{
				makespan = c[si_best][i].back();
				mm[si_best] = i;
			}
		}
		save_solution(si_best, si_opt, min_obj_iter, rc);
		cout << sol_obj[si_best] << "\t" << sol_obj[si_opt] - sol_obj[si_best] << "\t"
			<< makespan << "\t" << mm[si_best]
			<< endl;
		sum_delta_obj += (sol_obj[si_opt] - sol_obj[si_best]);
		/*display_solution(si_opt);*/
		//display_solution(si_best);		
	}
	cout << sum_delta_obj << endl;
}
void PMS::ejection_chain_local_search(int iteration, int perturb_rate, R_Mode r_mode, NS_Mode ns, int run_cnt_from, int run_cnt_to)
{
	int si_opt = 0, si_best = 1,
		si_cur = 2, si_local = 3, si_ptr = 4, si_ref = 5, si_trial = 6;
	int opt_cnt = 0, imp_cnt = 0, non_imp_cnt = 0;
	int rt = time(NULL);
	//rt = 1491483583;
	srand(rt);
	ofs << ins_name << "\t" << n << "\t" << m << "\t"
		<< obj_given << "\t" << sol_obj[si_opt] << "\t" << rt << endl;
	cout << ins_name << "\t" << n << "\t" << m << "\t"
		<< obj_given << "\t" << sol_obj[si_opt] << "\t" << rt << endl;
	for (int rc = run_cnt_from; rc <= run_cnt_to; rc++)
	{
		start_tm = clock();
		init_solution(si_cur, r_mode);
		//display_solution(si_cur);
		replace_solution(si_best, si_cur);
		bool is_still_improve = true;
		int min_obj_iter = 0;
		while (is_still_improve)
		{
			is_still_improve = false;
			bool is_trial_improve = true;

			for (int eject_job1 = rand()%n+1,eje=1; eje <= n&&is_trial_improve;eje++, eject_job1=eject_job1%n+1)
			{
				bool is_located = true;
				int i1, j1;
				for (int i3 = 1; i3 <= m&&is_located; i3++)
				{
					if (!effe_mach[si_cur][i3])
						continue;
					for (int j3 = 1; j3 < s[si_cur][i3].size(); j3++)
					{
						if (s[si_cur][i3][j3] == eject_job1)
						{
							i1 = i3;
							j1 = j3;
							is_located = false;
							break;
						}
					}
				}
				int len = 1, just_eject_mach = i1;
				//display_solution(si_cur);
				replace_solution(si_ref, si_cur);
				remove_job(si_ref, i1, j1, s[si_ref][i1][j1]);
				//display_solution(si_ref);
				while (len <= cl)
				{
					//trial_move(si_trial, i1, eject_job1, si_ref);
					replace_solution(si_trial, si_ref);
					obj_type min_c = DBL_MAX;
					int min_c_mach, min_c_pos;
					for (int i4 = 1; i4 <= m; i4++)
					{
						if (i4 == i1)
							continue;
						int probe_pos = 0;
						obj_type min_delta_c = -c[si_trial][i4].back();
						add_job(si_trial, i4, probe_pos, eject_job1);
						min_delta_c += c[si_trial][i4].back();
						if (min_c - min_delta_c > MIN_EQUAL)
						{
							min_c = min_delta_c;
							min_c_mach = i4;
							min_c_pos = probe_pos;
						}
						remove_job(si_trial, i4, probe_pos, eject_job1);
					}
					add_job(si_trial, min_c_mach, min_c_pos, eject_job1);
					calculate_obj(si_trial);
					min_obj_iter += 1;
					//display_solution(si_trial);
					//check_solution(si_trial);
					if (sol_obj[si_cur] - sol_obj[si_trial] > MIN_EQUAL)
					{
						replace_solution(si_cur, si_trial);
						is_trial_improve = false;
						break;
					}
					if (rand() % 100 < perturb_rate)
					{
						int min_mach, min_job, pos_eject2;
						obj_type min_delta_obj = DBL_MAX;
						if (m == 2)
						{
							if (effe_mach[si_ref][1] == false)
								just_eject_mach = 1;
							if (effe_mach[si_ref][2] == false)
								just_eject_mach = 2;
						}
						for (int i2 = 1; i2 <= m; i2++)
						{
							if (i2 == just_eject_mach || !effe_mach[si_ref][i2])
								continue;
							for (int j2 = 1; j2 < s[si_ref][i2].size(); j2++)
							{
								obj_type delta_obj = -c[si_ref][just_eject_mach].back() - c[si_ref][i2].back();
								int eject_job2 = s[si_ref][i2][j2];
								pos_eject2 = 0;
								remove_job(si_ref, i2, j2, eject_job2);
								add_job(si_ref, just_eject_mach, pos_eject2, eject_job2);
								delta_obj += (c[si_ref][just_eject_mach].back() + c[si_ref][i2].back());
								if (min_delta_obj - delta_obj > MIN_EQUAL)
								{
									min_delta_obj = delta_obj;
									min_mach = i2;
									min_job = j2;
								}
								add_job(si_ref, i2, j2, eject_job2);
								remove_job(si_ref, just_eject_mach, pos_eject2, eject_job2);
							}
						}
						pos_eject2 = 0;
						add_job(si_ref, just_eject_mach, pos_eject2, s[si_ref][min_mach][min_job]);
						remove_job(si_ref, min_mach, min_job, s[si_ref][min_mach][min_job]);
						just_eject_mach = min_mach;
					}
					else
					{
						int rand_m = rand() % m + 1;
						if (m == 2)
						{
							if (effe_mach[si_ref][1] == false)
								just_eject_mach = 1;
							if (effe_mach[si_ref][2] == false)
								just_eject_mach = 2;
						}
						while (rand_m == just_eject_mach || !effe_mach[si_ref][rand_m])
							rand_m = rand() % m + 1;
						int rand_j = rand() % (s[si_ref][rand_m].size() - 1) + 1;
						int eject_job2 = s[si_ref][rand_m][rand_j];
						int pos_eject2 = 0;
						//display_solution(si_ref);
						add_job(si_ref, just_eject_mach, pos_eject2, eject_job2);
						remove_job(si_ref, rand_m, rand_j, eject_job2);
						just_eject_mach = rand_m;
						//display_solution(si_ref);
					}
					len += 1;
				}
			}
			if (sol_obj[si_best] - sol_obj[si_cur] > MIN_EQUAL)
			{
				end_tm = clock();
				replace_solution(si_best, si_cur);
				//check_solution(si_best);
				is_still_improve = true;
				//cout << sol_obj[si_best] << endl;
			}
		}
		makespan = 0;
		for (int i3 = 1; i3 <= m; i3++)
		{
			if (!effe_mach[si_best][i3])
				continue;
			if (c[si_best][i3].back() - makespan > MIN_EQUAL)
			{
				makespan = c[si_best][i3].back();
				mm[si_best] = i3;
			}
		}
		save_solution(si_best, si_opt, min_obj_iter, rc);
		//display_solution(si_best);
		cout << sol_obj[si_best] << "\t" << sol_obj[si_opt] - sol_obj[si_best] << "\t"
			<< makespan << "\t" << mm[si_best]
			<< endl;
	}
}
void PMS::hma(int iteration, int perturb_rate, R_Mode r_mode, NS_Mode ns, int run_cnt_from, int run_cnt_to)
{
	int si_opt = 0, si_best = sol_num - non_popu + 2,
		si_offspring = sol_num - non_popu + 1;
	int opt_cnt = 0, imp_cnt = 0, non_imp_cnt = 0;
	obj_type sum_delta_obj = 0;
	int rt = time(NULL);
#ifdef DEBUG
	//rt = 1462129579;
#endif
	srand(rt);
	ofs << ins_name << "\t" << n << "\t" << m << "\t"
		<< obj_given << "\t" << sol_obj[si_opt] << "\t" << rt << endl;
	cout << ins_name << "\t" << n << "\t" << m << "\t"
		<< obj_given << "\t" << sol_obj[si_opt] << "\t" << rt << endl;
	for (rc = run_cnt_from; rc <= run_cnt_to; rc++)
	{
		tabu_pool_update.assign(sol_num, 0);
		start_tm = clock();
		init_solution(si_best, r_mode);
		end_tm = clock();
		for (int p = 1; p <= sol_num - non_popu; p++)
		{
			init_solution(p, r_mode);
			//display_solution(p);
			if (whe_dc == 1)
				divide_and_conquer(p, m);
			else
			{
				if (ls_method == 0)
					local_search_hybrid(p);
				else if (ls_method == 1)
					local_search_hybrid_tri_insert(p);
				else if (ls_method == 2)
					local_search_hybrid_tri_swap(p);
				else if (ls_method == 3)
					//local_search_hybrid_tri_insert_swap(p);
					local_search(p);
			}
			if (sol_obj[p] < sol_obj[si_best])
			{
				end_tm = clock();
				replace_solution(si_best, p);
				if (n <= 14 && fabs(sol_obj[si_best] - sol_obj[si_opt]) <= MIN_EQUAL)
					break;
			}
		}
		int min_obj_iter = 0;
		if (n <= 14)	// find the optimal solution for instances in OB set
		{
			if (fabs(sol_obj[si_best] - sol_obj[si_opt]) <= MIN_EQUAL)
				iteration = 0;
			else
				iteration = INT16_MAX;
		}
		for (gen = 0; gen < iteration; gen++)
		{
			int p1 = rand() % (sol_num - non_popu) + 1;
			int p2 = rand() % (sol_num - non_popu) + 1;
			while (p1 == p2)
				p2 = rand() % (sol_num - non_popu) + 1;
			if (rand() % 100 < alpha_cx_ptr)
			{
				if (cx_method == 0)
					ufx(si_offspring, p1, p2);
				else if (cx_method == 1)
					mpx(si_offspring, p1, p2);
				else if (cx_method == 2)
					gpx(si_offspring, p1, p2);
			}
			else
			{
				replace_solution(si_offspring, p1);
				perturb1(si_offspring, perturb_rate);
			}
			if (whe_dc == 1)
				divide_and_conquer(si_offspring, m);
			else
			{
				if (ls_method == 0)
					local_search_hybrid(si_offspring);
				else if (ls_method == 1)
					local_search_hybrid_tri_insert(si_offspring);
				else if (ls_method == 2)
					local_search_hybrid_tri_swap(si_offspring);
				else if (ls_method == 3)
					//local_search_hybrid_tri_insert_swap(si_offspring);
					local_search(si_offspring);
			}
			if (sol_obj[si_best] - sol_obj[si_offspring] > MIN_EQUAL)
			{
				end_tm = clock();
				replace_solution(si_best, si_offspring);
				min_obj_iter = gen;
				//cout << rc << "\t" << i<<"\t"<<sol_obj[si_best] << endl;
				if (n <= 14 && fabs(sol_obj[si_best] - sol_obj[si_opt]) <= MIN_EQUAL)
					break;
			}
			/*cout << rc << "\t" << gen << "\t";
			for (int i = 1; i <= sol_num - non_popu; i++)
			cout << sol_obj[i] << "\t";
			cout << sol_obj[si_offspring] << "\t" << sol_obj[si_best] << "\t";*/
			/*if (rc == 6&&(gen==32||gen==33))
			check_solution(si_offspring);*/
			if (pu_method == 0)
				pool_update(si_offspring, gen, 0);
			else if (pu_method == 1)
				pool_update_qd(si_offspring, gen, 0);
			//cout << endl;
			/*cout << temperature << "\t"
			<< obj_si_cur <<"\t"
			<<sol_obj[si_ptr]<<"\t"
			<<rand_num<<"\t"
			<< (100 * exp((obj_si_cur - sol_obj[si_ptr]) / (1*temperature)))<<"\t"
			<<accept
			<< endl;*/
			/*temperature *= control_para;
			if (temperature < 1)
			temperature = temperature0;*/
		}
		//check_solution(si_best); 
		save_solution(si_best, si_opt, min_obj_iter, rc);
#ifdef DEBUG
		cout << sol_obj[si_best] << "\t" << sol_obj[si_opt] - sol_obj[si_best] << endl;
		sum_delta_obj += (sol_obj[si_opt] - sol_obj[si_best]);
#endif
		/*display_solution(si_opt);*/
		//display_solution(si_best);	
	}
#ifdef DEBUG
	cout << sum_delta_obj << endl;
#endif
}
void PMS::gpx(int offspring, int pa1, int pa2)
{
	int p1 = sol_num - non_popu + 3, p2 = sol_num - non_popu + 4;
	replace_solution(p1, pa1);
	replace_solution(p2, pa2);
	vector<bool> derived_vec(m + 1, true);
	if (sol_obj[p2] - sol_obj[p1] > MIN_EQUAL)
	{
		int temp = p1;
		p1 = p2;
		p2 = temp;
	}

	/*if (rc == 6 && gen == 33)
	{
	for (int i = 1; i <= m; i++)
	{
	s[offspring][i].resize(1, 0);
	q[offspring][i].resize(1, 1);
	c[offspring][i].resize(1, 0);
	}
	display_solution(p1);
	display_solution(p2);
	check_solution(p1);
	check_solution(p2);
	}*/
	for (int i = 1; i <= m; i++)
	{
		int p = (i % 2 == 1) ? p1 : p2;
		int po = (i % 2 == 1) ? p2 : p1;
		obj_type max_obj = -1;
		int max_obj_mach;
		for (int j = 1; j <= m; j++)
		{
			if (!derived_vec[j])
				continue;
			if (c[p][j].back() - max_obj > MIN_EQUAL)
			{
				max_obj = c[p][j].back();
				max_obj_mach = j;
			}
		}
		s[offspring][max_obj_mach].assign(s[p][max_obj_mach].begin(), s[p][max_obj_mach].end());
		q[offspring][max_obj_mach].assign(q[p][max_obj_mach].begin(), q[p][max_obj_mach].end());
		c[offspring][max_obj_mach].assign(c[p][max_obj_mach].begin(), c[p][max_obj_mach].end());
		effe_mach[offspring][max_obj_mach] = true;
		for (int j = 1; j < s[offspring][max_obj_mach].size(); j++)
		{
			for (int k = 1; k <= m; k++)
			{
				for (int j1 = 1; j1 < s[po][k].size(); j1++)
				{
					if (s[po][k][j1] == s[offspring][max_obj_mach][j])
					{
						remove_job(po, k, j1, s[po][k][j1]);
						break;
					}
				}
			}
		}
		s[p][max_obj_mach].resize(1, 0);
		q[p][max_obj_mach].resize(1, 0);
		c[p][max_obj_mach].resize(1, 0);
		derived_vec[max_obj_mach] = false;
		/*if (rc == 6 && gen == 33)
		{
		display_solution(p1);
		display_solution(p2);
		cout << p << "\t" << max_obj_mach << endl;
		display_solution(offspring);
		}*/
	}
	calculate_obj(offspring);
	/*if (rc == 6 && gen == 33)
	{
	display_solution(offspring);
	check_solution(offspring);
	}*/
	for (int i = 1; i <= m; i++)
	{
		if (s[p1][i].size() == 1)
			continue;
		for (int j = 1; j < s[p1][i].size(); j++)
		{
			proc_type min_completion_time = DBL_MAX;
			int min_ct_mach_index, min_position;
			for (int k = 1; k <= m; k++)
			{
				int position = 0;
				add_job(offspring, k, position, s[p1][i][j]);
				if (c[offspring][k].back() < min_completion_time)
				{
					min_completion_time = c[offspring][k].back();
					min_ct_mach_index = k;
					min_position = position;
				}
				remove_job(offspring, k, position, s[p1][i][j]);
			}
			add_job(offspring, min_ct_mach_index, min_position, s[p1][i][j]);
		}
	}
	calculate_obj(offspring);
}
void PMS::mpx(int offspring, int pa1, int pa2)
{
	int p1 = sol_num - non_popu + 3, p2 = sol_num - non_popu + 4,
		p_flag = sol_num - non_popu + 5;
	replace_solution(p1, pa1);
	replace_solution(p2, pa2);
	replace_solution(p_flag, pa1);
	for (int i = 1; i <= m; i++)
	{
		s[offspring][i].resize(1, 0);
		q[offspring][i].resize(1, 1);
		c[offspring][i].resize(1, 0);
		for (int j = 1; j < s[p_flag][i].size(); j++)
			s[p_flag][i][j] = 1;
	}
	//display_solution(p_flag);
	for (int i = 1; i <= m; i++)
	{
		for (int j1 = 1; j1 < s[p1][i].size(); j1++)
		{
			for (int j2 = 1; j2 < s[p2][i].size(); j2++)
			{
				if (s[p1][i][j1] == s[p2][i][j2])
				{
					int position = 0;
					add_job(offspring, i, position, s[p1][i][j1]);
					s[p_flag][i][j1] = 0;
					break;
				}
			}
		}
	}
	//calculate_obj(offspring);
	/*display_solution(p1);
	display_solution(p2);
	display_solution(p_flag);
	display_solution(offspring);*/
	//check_solution(offspring);
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j < s[p1][i].size(); j++)
		{
			if (s[p_flag][i][j] == 0)	// common jobs derived
				continue;
			proc_type min_completion_time = DBL_MAX;
			int min_ct_mach_index, min_position;
			for (int k = 1; k <= m; k++)
			{
				int position = 0;
				add_job(offspring, k, position, s[p1][i][j]);
				if (c[offspring][k].back() < min_completion_time)
				{
					min_completion_time = c[offspring][k].back();
					min_ct_mach_index = k;
					min_position = position;
				}
				remove_job(offspring, k, position, s[p1][i][j]);
			}
			add_job(offspring, min_ct_mach_index, min_position, s[p1][i][j]);
		}
	}
	calculate_obj(offspring);
	//display_solution(offspring);
	//check_solution(offspring);
}
void PMS::ufx(int offspring, int pa1, int pa2)
{
	int p1 = sol_num - non_popu + 3;
	replace_solution(p1, pa1);
	for (int i = 1; i <= m; i++)
	{
		s[offspring][i].resize(1, 0);
		q[offspring][i].resize(1, 1);
		c[offspring][i].resize(1, 0);
	}
	for (int i = 1; i <= m; i++)
	{
		for (int j1 = 1; j1 < s[p1][i].size(); j1++)
		{
			if (rand() % 2 == 1)
				continue;
			int position_add = 0, position_remove = j1;
			add_job(offspring, i, position_add, s[p1][i][j1]);
			remove_job(p1, i, position_remove, s[p1][i][j1]);
		}
	}
	//calculate_obj(offspring);
	/*display_solution(pa1);
	display_solution(p1);
	display_solution(p2);
	display_solution(offspring);*/
	//check_solution(offspring);
	for (int i1 = 1; i1 <= m; i1++)
	{
		for (int j1 = 1; j1 < s[p1][i1].size(); j1++)
		{
			bool is_not_derived = true;
			for (int i2 = 1; i2 <= m&&is_not_derived; i2++)
			{
				for (int j2 = 1; j2 < s[pa2][i2].size(); j2++)
				{
					if (s[p1][i1][j1] == s[pa2][i2][j2])
					{
						int position_add = 0;
						add_job(offspring, i2, position_add, s[pa2][i2][j2]);
						is_not_derived = false;
						break;
					}
				}
			}
		}
	}
	calculate_obj(offspring);
	/*display_solution(offspring);
	check_solution(offspring);*/
}
void PMS::pool_update(int offspring, int gen, int p2)
{
	obj_type obj_worst = 0;
	int p_worst, equ_cnt;
	int rand_tt = sol_num / 4 + rand() % (sol_num / 2);
	for (int i = 1; i <= sol_num - non_popu; i++)
	{
		if (sol_obj[i] - obj_worst > MIN_EQUAL && (gen >= tabu_pool_update[i] + rand_tt))
		{
			obj_worst = sol_obj[i];
			equ_cnt = 1;
			p_worst = i;
		}
		else if (fabs(sol_obj[i] - obj_worst) <= MIN_EQUAL)
		{
			equ_cnt += 1;
			if (rand() % equ_cnt == 0)
				p_worst = i;
		}
	}
	if (obj_worst == 0)
	{
		p_worst = rand() % (sol_num - non_popu) + 1;
		//cout<<" at ";
	}
	if (sol_obj[p_worst] - sol_obj[offspring] > MIN_EQUAL ||
		rand() % 100 < 30)
	{
		//cout << p_worst;
		replace_solution(p_worst, offspring);
		tabu_pool_update[p_worst] = gen;
	}
}
// quality and distance based pool updating rule
void PMS::pool_update_qd(int offspring, int gen, int p2)
{
	obj_type obj_worst = 0;
	int p_worst, equ_cnt;
	vector<vector<obj_type>> dis(sol_num - non_popu + 2, vector<obj_type>(sol_num - non_popu + 2, 0));
	int rand_tt = sol_num / 4 + rand() % (sol_num / 2);
	for (int ins1 = 1; ins1 < sol_num - non_popu + 1; ins1++)
	{
		for (int ins2 = ins1 + 1; ins2 <= sol_num - non_popu + 1; ins2++)
		{
			dis[ins1][ins2] = 0;
			for (int i = 1; i <= m; i++)
			{
				for (int j1 = 1; j1 < s[ins1][i].size(); j1++)
				{
					for (int j2 = 1; j2 < s[ins2][i].size(); j2++)
					{
						if (s[ins1][i][j1] == s[ins2][i][j2])
						{
							dis[ins1][ins2] += 1;
							break;
						}
					}
				}
			}
			dis[ins1][ins2] /= n;
			dis[ins1][ins2] = 1 - dis[ins1][ins2];
			dis[ins2][ins1] = dis[ins1][ins2];
		}
	}
	vector<obj_type> distance(sol_num - non_popu + 2, 1);
	vector<obj_type> dis_port(2, 0), obj_port(2, 0);
	int min = 0, max = 1;
	for (int ins1 = 1; ins1 <= sol_num - non_popu + 1; ins1++)
	{
		for (int ins2 = 1; ins2 <= sol_num - non_popu + 1; ins2++)
		{
			if (ins1 == ins2)
				continue;
			if (distance[ins1] > dis[ins1][ins2])
				distance[ins1] = dis[ins1][ins2];
		}
		if (ins1 == 1)
		{
			dis_port[min] = dis_port[max] = distance[ins1];
			obj_port[min] = obj_port[max] = sol_obj[ins1];
		}
		else
		{
			if (dis_port[min] > distance[ins1])
				dis_port[min] = distance[ins1];
			if (dis_port[max] < distance[ins1])
				dis_port[max] = distance[ins1];
			if (obj_port[min] > sol_obj[ins1])
				obj_port[min] = sol_obj[ins1];
			if (obj_port[max] < sol_obj[ins1])
				obj_port[max] = sol_obj[ins1];
		}
	}
	vector<obj_type> goodscore(sol_num - non_popu + 2, 1);
	int worst_p;
	obj_type min_gs = 1;
	for (int ins = 1; ins <= sol_num - non_popu + 1; ins++)
	{
		goodscore[ins] = pu_beta*(distance[ins] - dis_port[min]) / (dis_port[max] - dis_port[min] + 1) +
			(1 - pu_beta)*(obj_port[max] - sol_obj[ins]) / (obj_port[max] - obj_port[min] + 1);
		if (min_gs > goodscore[ins] && ins < sol_num - non_popu + 1)
		{
			min_gs = goodscore[ins];
			worst_p = ins;
		}
		//cout << goodscore[ins] << "\t";
	}
	//cout << worst_p << "\t"<<min_gs<<"\t";
	if (goodscore[offspring] > goodscore[worst_p] || rand() % 100 < 30)
	{
		replace_solution(worst_p, offspring);
		//cout << "r";
	}
	//cout << endl;
}
void PMS::test()
{
	/*cout << endl << "This is test function." << endl;*/
	int si[] = { 5,1 ,2 ,3,8,9 };
	proc_type rd[] = { 7000,2000,4000,1400,900,3400,1239,987,2390,5542 };
	vector<proc_type> rd_vec(rd, rd + sizeof(rd) / sizeof(proc_type));
	vector<int> si_vec(si, si + sizeof(si) / sizeof(int));

	sort(si_vec.begin() + 2, si_vec.end(), cmpSort(rd_vec));
	for (int i = 0; i < si_vec.size(); i++)
		cout << si_vec[i] << "," << rd_vec[si_vec[i]] << "\t";
	cout << endl << "end of test function." << endl;

	std::vector<int> second_vec(5, 100);
	std::vector<int> myvector;
	myvector.assign(second_vec.begin(), second_vec.end());
	std::copy(second_vec.begin(), second_vec.end(), myvector.begin());
	std::cout << "myvector contains:";
	for (std::vector<int>::iterator it = myvector.begin(); it != myvector.end(); ++it)
		std::cout << ' ' << *it;
	std::cout << '\n';
}
void run_algorithm(std::map<string, string> &argv_map, string ins_name)
{
	string fnr = argv_map.at("_if") + ins_name;
	string fnw = argv_map.at("_of") /*+ argv_map.at("_exe_name")*/ + argv_map.at("_px") +
		"_am" + argv_map.at("_am") +
		"_dc" + argv_map.at("_dc") +
		"_ls" + argv_map.at("_ls") +
		"_ox" + argv_map.at("_xo") +
		"_pu" + argv_map.at("_pu") +
		"_p" + argv_map.at("_p") +
		"_np" + argv_map.at("_non_popu") +
		"_itr" + argv_map.at("_itr") +
		"_ptr" + argv_map.at("_ptr") +
		"_cl" + argv_map.at("_cl") +
		"_rm" + argv_map.at("_rm") +
		//"_ns" + argv_map.at("_ns") +
		"_alpha" + argv_map.at("_cx") +
		"_beta" + argv_map.at("_pu_beta") +
		"_r" + argv_map.at("_r1") +
		"_r" + argv_map.at("_r2") + ".txt";
	PMS *pms = new PMS(fnr, fnw, stoi(argv_map.at("_p")));//Ni_14_4-1_1_15
	pms->ins_name = ins_name;
	pms->whe_save_sol_seq = stoi(argv_map.at("_ws"));
	pms->temperature = stoi(argv_map.at("_t"));
	pms->control_para = stoi(argv_map.at("_cp"))*0.01;
	pms->alpha_cx_ptr = stoi(argv_map.at("_cx"));
	pms->non_popu = stoi(argv_map.at("_non_popu"));
	pms->pu_beta = stoi(argv_map.at("_pu_beta"))*0.01;
	pms->whe_dc = stoi(argv_map.at("_dc"));
	pms->ls_method = stoi(argv_map.at("_ls"));
	pms->pu_method = stoi(argv_map.at("_pu"));
	pms->cx_method = stoi(argv_map.at("_xo"));
	pms->cl = stoi(argv_map.at("_cl"));
	pms->am = stoi(argv_map.at("_am"));
	//pms->test();
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
	/*pms->hma(stoi(argv_map.at("_itr")), stoi(argv_map.at("_ptr")),
		r_mode, ns_mode, stoi(argv_map.at("_r1")), stoi(argv_map.at("_r2")));*/
	if (pms->am == 0)
		pms->iterated_local_search(stoi(argv_map.at("_itr")), stoi(argv_map.at("_ptr")),
			r_mode, ns_mode, stoi(argv_map.at("_r1")), stoi(argv_map.at("_r2")));
	else if (pms->am == 1)
		pms->ejection_chain_local_search(stoi(argv_map.at("_itr")), stoi(argv_map.at("_ptr")),
			r_mode, ns_mode, stoi(argv_map.at("_r1")), stoi(argv_map.at("_r2")));
	delete pms;
}
int main(int argc, char **argv)
{
	char *rgv_ins[] = { "",
		"_p","20", "_itr","2000",
		"_rm","9", "_ptr","70",
		"_r1","1", "_r2","20",
		"_ws","0", "_t","50",
		"_cp","95", "_ns","0",
		"_non_popu","10","_cx","30",
		"_pu_beta","60", "_px","hma",
		"_dc","0", "_ls","3",
		"_pu","1",	"_xo","2",
		"_cl", "5",	"_am", "1",

		"_vi1","0",		"_vi2","2",
		"_ni1","0",		"_ni2","3",
		"_vj1","0",		"_vj2","2",
		"_mj1","0",		"_mj2","3",
		"_pi1","1",		"_pi2","2",
		"_di1","1",		"_di2","2",
		"_ins1","1",	"_ins2","25"
	};
	char *rgv[] = { "",	//0
		"_px","hma_test",	//prefix for output file name
		"_if","instance\\BB_Problem_BestSolution\\",	//input file directory
		"_of","results\\",// output file directory
		"_p","20",		// population size
		"_itr","2000",	// max iteration of ILS
		"_ptr","70",	// perturbation rate 
		"_rm","9",	// construction rules for initial solution
		"_ns","0",	// neighborhood search, 0:swap, 1:insert	
		"_r1","1",	// run cnt from
		"_r2","20",	// run cnt to
		"_ws","0",	// whe_save_sol_seq
		"_t","50",	// initial temperature T
		"_cp","95",	// control para cp, T=T*cp
		"_cx","30",	// crossover vs perturbation rate
		"_pu_beta","60",	//pool update coefficient
		"_dc", "0",	// whether use divide and conquer
		"_ls", "3",	// local search method
		"_pu", "1",	// pool update method
		"_xo", "2",	// crossover method
		"_cl", "15", // length of the chain
		"_am", "1", // algorithm method
		"_vi1","1",		"_vi2","2",
		"_ni1","2",		"_ni2","3",
		"_vj1","1",		"_vj2","2",
		"_mj1","2",		"_mj2","3",
		"_pi1","2",		"_pi2","2",
		"_di1","2",		"_di2","2",
		"_ins1","25",	"_ins2","25"
	};
#ifdef DEBUG
	argc = sizeof(rgv) / sizeof(rgv[0]); argv = rgv;
#endif
	std::map<string, string> argv_map;
	argv_map["_exe_name"] = argv[0];	// add the exe file name to argv map, to append to the output file name
	for (int i = 1; i < sizeof(rgv_ins) / sizeof(rgv_ins[0]); i += 2)
		argv_map[string(rgv_ins[i])] = string(rgv_ins[i + 1]);
	for (int i = 1; i < argc; i += 2)
		argv_map[string(argv[i])] = string(argv[i + 1]);
	vector<vector<int>> n_vec{ { 8, 11, 14 },{ 20, 35, 50 } };
	vector<vector<int>> m_vec{ { 2, 3, 4 },{ 4, 7, 10 } };
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
#ifdef DEBUG
	system("pause");
#endif	
}
#endif