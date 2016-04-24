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
	vector<vector<proc_type>> p;	// processing time
	vector<vector<dete_type>> d;	// deterioration effect
	obj_type obj_given;	// the optimal objective value
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

	clock_t start_tm, end_tm;
	class cmpSort;
	const double MIN_EQUAL = 0.001;
	int rand_seed;
	ofstream ofs;
	int ls_cnt1, ls_cnt2;
public:
	enum R_Mode { EXAC = 0, MINP, MAXP, MIND, MAXD, MINPDD, MAXPDD, MINPD, MAXPD, RANDOM, SIZE };
	enum Cmp_Mode { GREATER, LESS };
	enum NS_Mode { SWAP, INSERT };
	int whe_save_sol_seq;
	string ins_name;
	double control_para, temperature;
	PMS(string, string, int);
	~PMS();
	void display_problem();
	void display_solution(int);
	void check_solution(int);
	void init_solution(int, R_Mode);
	void calculate_completion_time(int, int);
	void calculate_obj(int);
	void add_job(int, int, int&, int);
	void remove_job(int, int,int&, int);
	void exchange_job(int, int, int&, int&, int, int);
	void local_search(int, int, NS_Mode);
	void local_search_hybrid(int);
	void local_search_hybrid1(int, int, NS_Mode);
	void local_search_ejection_chain(int, int, NS_Mode);
	void iterated_local_search(int, int, R_Mode, NS_Mode, int, int);
	void local_search_recursion(int);
	void divide_and_conquer(int si,int mach_num);
	void local_search1(int si, int *sub_mach, int sub_m, obj_type &sub_obj);
	void insert1(int si, int *sub_mach, int mm_j, int mo, obj_type &sub_obj);
	void swap1(int si, int *sub_mach, int mm_j, int mo, int mo_j, obj_type &sub_obj);
	void replace_solution(int, int);
	void perturb(int, int);
	void perturb1(int, int);
	void swap(int, int, int, int, int);
	void insert(int, int, int, int);
	void trail_move(int, int, int, int, int);
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
	int si_opt = 0,num_job;
	for (int i = 1; i <= m; i++)
	{
		ifs >> num_job;	// the optimal number of jobs to be assigned to machine i	
		s[si_opt][i].resize(num_job+1);
		effe_mach[si_opt][i] = num_job > 0 ? true : false;	// the first element indicates the effectiveness of machine i
		q[si_opt][i].resize(num_job + 1);
		c[si_opt][i].resize(num_job + 1);
		for (int j = 1; j <= num_job; j++)
			ifs >> s[si_opt][i][j];	// the list of jobs to be assinged to machine i
	}
	ifs >> obj_given;
	for (int i = 1; i <= m; i++)
	{
		sort(s[si_opt][i].begin()+1, s[si_opt][i].end(), cmpSort(r[i]));
		calculate_completion_time(si_opt, i);
	}
	calculate_obj(si_opt);
	if (abs(obj_given - sol_obj[si_opt]) > MIN_EQUAL)
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
		cout <<effe_mach[si][i]<<"\t"<< s[si][i].size() - 1 << ", " << c[si][i].back() << ": ";
		for (int j = 1; j <s[si][i].size(); j++)
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
	obj_type cl, max_cl = 0;
	int max_cl_mach_index, sum_job_index = 0;
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j < s[si][i].size(); j++)
			sum_job_index += s[si][i][j];
		if (!effe_mach[si][i])
			continue;
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
				cout << "ERROR, si: " << si << " not sort by r"
					<< r[i][s[si][i][j - 1]] << " " << r[i][s[si][i][j]] << endl;
				display_solution(si);
				//system("pause");
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
				<<", mach_index: "<<i
				<< ", real completion time: " << cl << " " << c[si][i].back() << endl;
			system("pause");
		}
	}
	if (abs(max_cl - sol_obj[si])>MIN_EQUAL||max_cl_mach_index!=mm[si])
	{
		cout << "ERROR, or mm is wrong, si: " << si
			<< ", real obj: " << max_cl << " " << sol_obj[si] << endl;
		display_solution(si);
		system("pause");
	}
	if (sum_job_index != (1 + n)*n / 2)
	{
		cout << "ERROR, si: " << si
			<< ", multiple job index" << endl;
		//display_solution(si);
		system("pause");
	}
}
void PMS::calculate_completion_time(int si, int mach_index)
{
	//q[si][mach_index][1] = 1;
	//c[si][mach_index][1] = p[mach_index][s[si][mach_index][1]];
	for (int i = 1; i < s[si][mach_index].size(); i++)
	{
		q[si][mach_index][i] = (1 - d[mach_index][s[si][mach_index][i - 1]])*q[si][mach_index][i - 1];
		c[si][mach_index][i] = c[si][mach_index][i - 1] + p[mach_index][s[si][mach_index][i]] / q[si][mach_index][i];
	}
}
void PMS::calculate_obj(int si)
{
	sol_obj[si] = -1;
	for (int i = 1; i <= m; i++)
	{
		if (effe_mach[si][i] == true && (c[si][i].back() - sol_obj[si]) > MIN_EQUAL)
		{
			sol_obj[si] = c[si][i].back();
			mm[si] = i;
		}
	}
}
void PMS::add_job(int si, int mach_index, int &position, int job)
{
	if (s[si][mach_index].size()==1)	// empty machine
	{
		effe_mach[si][mach_index] = true;	// indicate it is effective
		s[si][mach_index].push_back(job);
		q[si][mach_index].push_back(1);
		c[si][mach_index].push_back(p[mach_index][job]);
	}else
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
		c[si][mach_index].insert(c[si][mach_index].begin() + position, c[si][mach_index][position - 1]+
			p[mach_index][s[si][mach_index][position]] / q[si][mach_index][position]);
		for (int i = position + 1; i < s[si][mach_index].size(); i++)
		{
			q[si][mach_index][i] = (1 - d[mach_index][s[si][mach_index][i - 1]])*q[si][mach_index][i - 1];
			c[si][mach_index][i] = c[si][mach_index][i - 1] + p[mach_index][s[si][mach_index][i]] / q[si][mach_index][i];
		}
	}
}void PMS::remove_job(int si, int mach_index, int &position, int job)
{
	if (s[si][mach_index].empty())
	{
		cout << "ERROR, remove job from empty machine" << endl;
		system("pause");
	}
	if( position ==0)
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
}void PMS::save_solution(int si, int si_opt, int iterration, int run_cnt)
{
	//end_tm = clock();
	int result_improve, given_result_improve;
	if (abs(sol_obj[si] - sol_obj[si_opt]) <= MIN_EQUAL)
		result_improve = 1;	// equal
	else if (sol_obj[si_opt] - sol_obj[si] > MIN_EQUAL)
		result_improve = 2;	// improved
	else
		result_improve = 0;	// not improved
	if (abs(sol_obj[si] - obj_given) <= MIN_EQUAL)
		given_result_improve = 1;
	else if (obj_given - sol_obj[si] > MIN_EQUAL)
		given_result_improve = 2;
	else
		given_result_improve = 0;
	ofs << run_cnt << "\t"
		//<< c[si_opt][0] << "\t"
		<< sol_obj[si] << "\t"
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
				<< s[si][i].size()-1 << "\t";
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
		int min_ct_mach_index,min_position;
		for (int i = 1; i <= m; i++)
		{
			int position = 0;
			add_job(si, i, position,job_sort_charact_vec[j]);	
			if (c[si][i].back() < min_completion_time)
			{
				min_completion_time = c[si][i].back();
				min_ct_mach_index = i;
				min_position = position;
			}
			remove_job(si, i, position, job_sort_charact_vec[j]);
		}
		add_job(si, min_ct_mach_index,min_position, job_sort_charact_vec[j]);
	}
	calculate_obj(si);
	//check_solution(si);
}
void PMS::perturb(int si, int ptr_rate)
{
	for (int i = 0; i < (s[si][mm[si]].size()-1) * ptr_rate*0.01; i++)
	{
		int pos_remove_mm = rand() % (s[si][mm[si]].size() - 1) + 1, pos_add_mm = 0;
		int job_add = s[si][mm[si]][pos_remove_mm];
		int mach_i = rand() % m + 1;
		while (mach_i == mm[si] || s[si][mach_i].size()== 1)
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
	for (int i = 0; i < (s[si][mm[si]].size() - 1) * ptr_rate*0.01; i++)
	{
		int jm = rand() % (s[si][mm[si]].size() - 1) + 1;
		int mach_r = rand() % m + 1;
		while (mach_r == mm[si] || s[si][mach_r].size() == 1)
			mach_r = rand() % m + 1;
		int jo = rand() % (s[si][mach_r].size() - 1) + 1;
		int temp = s[si][mm[si]][jm];
		s[si][mm[si]][jm] = s[si][mach_r][jo];
		s[si][mach_r][jo] = temp;
	}
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
		bool is_mm_unchanged = true;
		int pre_mm = mm[si];
		for (int jm = 0 % (s[si][mm[si]].size()-1) + 1, jm_r = 1;
		jm_r < s[si][mm[si]].size() && is_mm_unchanged;
			jm_r++, jm = jm % (s[si][mm[si]].size() - 1) + 1)
		{
			for (int i = 0 % m + 1, i_r = 1; i_r <= m && is_mm_unchanged; i_r++, i = i%m + 1)
			{
				if (i == mm[si]||!effe_mach[si][i])
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
					remove_job(si,i, pos_add, cur_job);
					add_job(si, pre_mm, pos_remove, cur_job);
					//calculate_obj(si);
					sol_obj[si] = pre_obj;
					mm[si] = pre_mm;
				}
				/*display_solution(si);
				check_solution(si);*/
				for (int j = rand() % ((s[si][i].size()-1) == 0 ? 1 : (s[si][i].size()-1)) + 1, j_r = 1;
				j_r < s[si][i].size() && is_mm_unchanged; j_r++, j = j%(s[si][i].size()-1) + 1)
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
					if (pre_obj-sol_obj[si] > MIN_EQUAL)
					{
						if (pre_mm != mm[si])
							is_mm_unchanged = false;
						is_still_improved = true;
						break;
					}
					else
					{
						exchange_job(si, i, pos_remove_mi,  pos_add_mi, job_remove,job_add);
						exchange_job(si, pre_mm,  pos_remove_mm, pos_add_mm, job_add,job_remove);
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
void PMS::local_search_recursion(int si)
{
	
	int sub_m = m;
	obj_type sub_obj = sol_obj[si];
	//ls_cnt1 = ls_cnt2 = 0;
	divide_and_conquer(si,m);
	//cout << sub_obj << endl;
	/*mm[si] = sub_mach[0];
	c[si][0] = sub_obj;*/
	//check_solution(si);
	//cout << ls_cnt1 << "\t" << ls_cnt2 << endl;
}

void PMS::divide_and_conquer(int si,int mach_num)
{
	if (mach_num == 2)
	{
		local_search_hybrid(si);
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
		effe_mach[si][sub_mach1[i]] = true;
		if (c[si][sub_mach1[i]].back() - sol_obj[si]>MIN_EQUAL)
		{
			sol_obj[si] = c[si][sub_mach1[i]].back();
			mm[si] = sub_mach1[i];
		}
	}
	divide_and_conquer(si, sub_mach1.size());	

	sol_obj[si] = 0;
	effe_mach[si].assign(m + 1, false);
	for (int i = 0; i < sub_mach2.size(); i++)
	{
		effe_mach[si][sub_mach2[i]] = true;
		if (c[si][sub_mach2[i]].back() - sol_obj[si]>MIN_EQUAL)
		{
			sol_obj[si] = c[si][sub_mach2[i]].back();
			mm[si] = sub_mach2[i];
		}
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
		if (c[si][i].back() - sol_obj[si] > MIN_EQUAL)
		{
			sol_obj[si] = c[si][i].back();
			mm[si] = i;
		}
	}
	local_search_hybrid(si);
	ls_cnt2 += 1;	
	//check_solution(si);
}
void PMS::iterated_local_search(int iteration, int perturb_rate, R_Mode r_mode, NS_Mode ns, int run_cnt_from, int run_cnt_to)
{
	int si_opt = 0, si_best = 1,
		si_cur = 2, si_local = 3, si_ptr = 4;
	int opt_cnt = 0, imp_cnt = 0, non_imp_cnt = 0;
	obj_type sum_obj = 0;
	int rt = time(NULL);
	//rt = 1461506514;
	srand(rt);
	ofs << ins_name << "\t" << n << "\t" << m << "\t"
		<< obj_given << "\t" << c[si_opt][0][0] << "\t" << rt << endl;
	cout << ins_name << "\t" << n << "\t" << m << "\t"
		<< obj_given << "\t" << c[si_opt][0][0] << "\t" << rt << endl;
	//double temperature0 = temperature;
	//temperature = temperature0;
	for (int rc = run_cnt_from; rc <= run_cnt_to; rc++)
	{
		start_tm = clock();
		init_solution(si_cur, r_mode);	//PMS::MINPDD
		//display_solution(si_cur);
		local_search_recursion(si_cur);
		//continue;
		//local_search_hybrid(si_cur);
		//local_search_ejection_chain(si_cur, si_local, ns);
		replace_solution(si_best, si_cur);
		//return;
		cout << rc << "\t0\t" << sol_obj[si_best] << endl;
		//save_solution(si_best, si_opt, 0, run_cnt);
		end_tm = clock();
		int min_obj_iter = 0;
		if (n <= 14)	// find the optimal solution for instances in OB set
		{
			if (abs(sol_obj[si_best] - sol_obj[si_opt]) <= MIN_EQUAL)
				iteration = 0;
			else
				iteration = INT16_MAX;
		}
		for (int i = 0; i < iteration; i++)
		{
			replace_solution(si_ptr, si_cur);
			perturb1(si_ptr, perturb_rate);
			//cout << sol_obj[si_best]<< ", " << sol_obj[si_cur] << endl;
			//local_search_ejection_chain(si_ptr, si_local, ns);
			//check_solution(si_ptr);
			//local_search_hybrid(si_ptr);
			local_search_recursion(si_ptr);
			if (sol_obj[si_best] - sol_obj[si_ptr]>MIN_EQUAL)
			{
				replace_solution(si_best, si_ptr);
				min_obj_iter = i;
				end_tm = clock();
				cout << rc << "\t" << i<<"\t"<<sol_obj[si_best] << endl;
				if (n <= 14 && abs(sol_obj[si_best] - sol_obj[si_opt]) <= MIN_EQUAL)
					break;	
			}
			/*int accept = 0;
			int rand_num = rand() % 100;
			obj_type obj_si_cur = sol_obj[si_cur];*/
			if (sol_obj[si_cur] - sol_obj[si_ptr] > MIN_EQUAL ||
				rand() % 100 <= (control_para*100))
			{
				replace_solution(si_cur, si_ptr);
				//accept = 1;
			}
			//cout << temperature << "\t"
			//	<< obj_si_cur <<"\t"
			//	<<sol_obj[si_ptr]<<"\t"
			//	<<rand_num<<"\t"
			//	//<< (100 * exp((obj_si_cur - sol_obj[si_ptr]) / (1*temperature)))<<"\t"
			//	<<accept
			//	<< endl;
			/*temperature *= control_para;
			if (temperature < 0.1)
				temperature = temperature0;*/
		}
		//check_solution(si_best); 
		save_solution(si_best, si_opt, min_obj_iter, rc);
		/*display_solution(si_opt);*/
		//display_solution(si_best);		
	}
}
void PMS::test()
{
	/*cout << endl << "This is test function." << endl;*/
	int si[] = { 5,1 ,2 ,3,8,9};
	proc_type rd[] = { 7000,2000,4000,1400,900,3400,1239,987,2390,5542 };
	vector<proc_type> rd_vec(rd, rd + sizeof(rd) / sizeof(proc_type));
	vector<int> si_vec(si, si + sizeof(si) / sizeof(int));

	sort(si_vec.begin()+2, si_vec.end(), cmpSort(rd_vec));
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
	PMS *pms = new PMS(fnr, fnw, stoi(argv_map.at("_p")));//Ni_14_4-1_1_15
	pms->ins_name = ins_name;
	pms->whe_save_sol_seq = stoi(argv_map.at("_ws"));
	pms->temperature = stoi(argv_map.at("_t"));
	pms->control_para = stoi(argv_map.at("_cp"))*0.01;
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
	pms->iterated_local_search(stoi(argv_map.at("_itr")), stoi(argv_map.at("_ptr")),
		r_mode, ns_mode, stoi(argv_map.at("_r1")), stoi(argv_map.at("_r2")));
	delete pms;
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
		"_px","tr_test",	//prefix for output file name
		"_if","instance\\BB_Problem_BestSolution\\",	//input file directory
		"_of","results\\",// output file directory
		"_p","13",		// population size
		"_itr","2000",	// max iteration of ILS
		"_ptr","50",	// perturbation rate 
		"_rm","9",	// construction rules for initial solution
		"_ns","0",	// neighborhood search, 0:swap, 1:insert	
		"_r1","1",	// run cnt from
		"_r2","20",	// run cnt to
		"_ws","0",	// whe_save_sol_seq
		"_t","2",	// initial temperature T
		"_cp","10",	// control para cp, T=T*cp
		"_vi1","1",		"_vi2","2",
		"_ni1","2",		"_ni2","3",
		"_vj1","0",		"_vj2","2",
		"_mj1","2",		"_mj2","3",
		"_pi1","1",		"_pi2","1",
		"_di1","1",		"_di2","1",
		"_ins1","1",	"_ins2","1"
	};
	argc = sizeof(rgv) / sizeof(rgv[0]); argv = rgv;
	std::map<string, string> argv_map;
	for (int i = 1; i < sizeof(rgv_ins) / sizeof(rgv_ins[0]); i += 2)
		argv_map[string(rgv_ins[i])] = string(rgv_ins[i + 1]);
	for (int i = 1; i < argc; i += 2)
		argv_map[string(argv[i])] = string(argv[i + 1]);
	vector<vector<int>> n_vec = { { 8, 11, 14 },{ 20,35,50 } };
	vector<vector<int>> m_vec = { { 2,3,4 },{ 4,7,10 } };
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