#if 0
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<map>
#include<math.h>
using namespace std;
const double MIN_EQUAL = 0.001;

template<typename T>
void split_generic(vector<T> &v, const T & str, const T & delimiters) {
	v.clear();	// clear v to be empty
	typename T::size_type start = 0;
	auto pos = str.find_first_of(delimiters, start);
	while (pos != T::npos) {
		if (pos != start) // ignore empty tokens
			v.emplace_back(str, start, pos - start);
		start = pos + 1;
		pos = str.find_first_of(delimiters, start);
	}
	if (start < str.length()) // ignore trailing delimiter
		v.emplace_back(str, start, str.length() - start); // add what's left of the string
}
// replace all occurance of t in s to w  
void replace_all(std::string & s, std::string const & t, std::string const & w)
{
	string::size_type pos = s.find(t), t_size = t.size(), w_size = w.size();
	while (pos != std::string::npos) { // found   
		s.replace(pos, t_size, w);
		pos = s.find(t, pos + w_size);
	}
}
class SolutionInfo
{
	/*ofs << run_cnt << "\t"
		<< sol_obj[si] << "\t"
		<< eca_imp_cnt << "\t"
		<< ls_iter << "\t"
		<< (double)(end_tm - start_tm) / CLOCKS_PER_SEC << "\t"
		<< endl;*/
public:
	int run_cnt, eca_imp_cnt, ls_iter;
	double obj, tm;
	SolutionInfo(int _rc, double _obj, int _eic, int _iter, double _tm)
		:run_cnt(_rc), obj(_obj), eca_imp_cnt(_eic), ls_iter(_iter), tm(_tm)
	{}
};
class InstanceInfo
{
public:
	string filename;
	int m, n, p, d, ins, rand_seed;
	vector<SolutionInfo*> &sol_info_vec;
	InstanceInfo(string _fn, int _n, int _m, int _p, int _d, int _ins, int _rd, vector<SolutionInfo*> &_siv)
		:filename(_fn), m(_m), n(_n), p(_p), d(_d), ins(_ins), rand_seed(_rd), sol_info_vec(_siv)
	{}
};
class Analyze {
public:
	class BestSolution { 
	public: 
		string ins_name,file_name;
		double obj;
		BestSolution(string _ins, string _fn, double _obj) :
			ins_name(_ins), file_name(_fn), obj(_obj) {}
	};
	Analyze(string, string, string, string, int);
	~Analyze();
	void total_info();
	void para_setting();
	double f, d, h, t, opt_best, set_ins_cnt;
	enum BM { MIN, AVG, MAX };
	int whe_save_each_file_result;
private:
	string fnr, fnw, fnwt, fnr_best;
	bool whe_input_table_head_fnw, whe_input_table_head_fnwt, best_obj_need_update;
	ofstream ofs, ofst,ofs_best_obj;
	vector<InstanceInfo*>insinfo_vec;
	map<string, BestSolution*> opt_best_map;
};
Analyze::Analyze(string _fnr, string _fnw, string _fnwt, string _fnr_best, int _ws) :
	fnr(_fnr), fnw(_fnw), fnwt(_fnwt), fnr_best(_fnr_best), whe_save_each_file_result(_ws),best_obj_need_update(false)
{
	ifstream ifs(fnr);
	if (!ifs.is_open())
	{
		cout << fnr << endl; perror("file_input fnr.");
		exit(0);
	}
	ifstream ifs_test(fnw);
	if (!ifs_test.is_open())
		whe_input_table_head_fnw = true;
	else
		whe_input_table_head_fnw = false;
	ifs_test.close();
	ifs_test.open(fnwt);
	if (!ifs_test.is_open())
		whe_input_table_head_fnwt = true;
	else
		whe_input_table_head_fnwt = false;
	ifs_test.close();
	if (whe_save_each_file_result)
	{
		ofs.open(fnw, ios::trunc | ios::out);
		ofs.setf(ios::fixed, ios::floatfield);
		ofs.precision(6);
	}
	ofst.open(fnwt, ios::app | ios::out);
	ofst.setf(ios::fixed, ios::floatfield);
	ofst.precision(6);

	string strline;
	vector<string> fields_vec, ins_name_vec;
	vector<SolutionInfo*>  *sol_info_vec = new vector<SolutionInfo*>;
	while (getline(ifs, strline))
	{
		InstanceInfo *insinfo;
		string filename;
		int m, n, rand_seed;
		double opt_obj_given, opt_obj_real;
		split_generic<string>(fields_vec, strline, "\t");

		if (fields_vec.front().find("N_") != string::npos)
		{
			split_generic<string>(ins_name_vec, fields_vec[0], "_.");
			replace_all(fields_vec[0], ".txt", "");
			sol_info_vec = new vector<SolutionInfo*>;
			insinfo = new InstanceInfo(fields_vec[0],	// file name
				stoi(fields_vec[1]), stoi(fields_vec[2]),	// n, m
				stoi(ins_name_vec[3]), stoi(ins_name_vec[4]), stoi(ins_name_vec[5]),	// p, d, ins
				stod(fields_vec[3]),	// rand_seed
				*sol_info_vec);
			insinfo_vec.push_back(insinfo);
		}
		else
		{
			if (sol_info_vec->size() == 0 || sol_info_vec->back()->run_cnt != stoi(fields_vec[0]))
			{
				sol_info_vec->push_back(new SolutionInfo(stoi(fields_vec[0]),//run_cnt
					stod(fields_vec[1]),	// obj
					stoi(fields_vec[2]),	// eca_imp_cnt
					stoi(fields_vec[3]),	// ls_iter
					stod(fields_vec[4]))	// tm
					);
			}
			else
			{
				sol_info_vec->back()->run_cnt = stoi(fields_vec[0]);
				sol_info_vec->back()->obj = stod(fields_vec[1]);
				sol_info_vec->back()->eca_imp_cnt = stoi(fields_vec[2]);
				sol_info_vec->back()->ls_iter = stoi(fields_vec[3]);
				sol_info_vec->back()->tm = stod(fields_vec[4]);
			}
		}
	}
	ifs.close();
	ifstream ifs_best(fnr_best);
	if (ifs_best.is_open())
	{
		while (getline(ifs_best, strline))
		{
			split_generic<string>(fields_vec, strline, "\t");
			opt_best_map[fields_vec[0]] = new BestSolution(fields_vec[0], fields_vec[2], stod(fields_vec[1]));
		}
	}
	ifs_best.close();
}
Analyze::~Analyze()
{
	if(best_obj_need_update)
	{
		ofs_best_obj.open(fnr_best, ios::trunc | ios::out);
		ofs_best_obj.setf(ios::fixed, ios::floatfield);
		ofs_best_obj.precision(6);
		if (!ofs_best_obj.is_open())
		{
			cout << fnw << endl; perror("file_output update best obj.");
			exit(0);
		}
		for (map<string, BestSolution*>::iterator iter = opt_best_map.begin();
		iter != opt_best_map.end(); iter++)
		{
			ofs_best_obj << iter->second->ins_name << "\t"
				<< iter->second->obj << "\t"
				<< iter->second->file_name << "\t"
				<< endl;
		}
		ofs_best_obj.close();
	}
	ofs.close();
	ofst.close();
}
void Analyze::para_setting()
{
	if (whe_input_table_head_fnw)
		ofs << "(*iter)->filename  \t(*iter)->n  \t(*iter)->m  \t"
		<< "(*iter)->rand_seed  \t(*iter)->sol_info_vec.size()  \t"
		<< "obj[MIN]  \t  obj[AVG]  \t  obj[MAX]  \t"
		<< "eca_imp_cnt[MIN]  \t  eca_imp_cnt[AVG]  \t  eca_imp_cnt[MAX]  \t"
		<< "ls_iter[MIN]  \t  ls_iter[AVG]  \t  ls_iter[MAX]  \t"
		<< "tm[MIN]  \t  tm[AVG]  \t  tm[MAX]  \t"
		<< "hit_cnt"
		<< endl;
	if (whe_input_table_head_fnwt)
		ofst << " fnr  \t f  \t  d  \t  h  \t  t  \t"
		<< endl;
	int num_best_ins = 0;
	for (vector<InstanceInfo*>::iterator iter = insinfo_vec.begin();
	iter != insinfo_vec.end(); iter++)
	{
		cout << (*iter)->filename << ", " << (*iter)->m << ", " << (*iter)->n << ", "
			<< (*iter)->rand_seed << ", " << (*iter)->sol_info_vec.size() << endl;
		double obj[3], tm[3];
		int eca_imp_cnt[3], ls_iter[3];
		int hit_cnt = 0, min_obj_cnt = 0, best_obj_cnt = 0;
		for (vector<SolutionInfo*>::iterator iter_sol = (*iter)->sol_info_vec.begin();
		iter_sol != (*iter)->sol_info_vec.end(); iter_sol++)
		{
			if (opt_best_map.find((*iter)->filename) != opt_best_map.end())	// the instance exists
			{
				if (opt_best_map.at((*iter)->filename)->obj - (*iter_sol)->obj > MIN_EQUAL)
				{
					opt_best_map.at((*iter)->filename)->obj = (*iter_sol)->obj;
					opt_best_map.at((*iter)->filename)->file_name = fnr;
					best_obj_cnt = 1;
					best_obj_need_update = true;	// need update
				}
				else if (fabs(opt_best_map.at((*iter)->filename)->obj - (*iter_sol)->obj) <= MIN_EQUAL)
					best_obj_cnt += 1;
			}
			else // the instance does not exist
			{
				opt_best_map[(*iter)->filename] = new BestSolution((*iter)->filename, fnr, (*iter_sol)->obj);
				best_obj_need_update = true;	// need update
				best_obj_cnt = 1;
			}
			if (iter_sol == (*iter)->sol_info_vec.begin())
			{
				obj[MIN] = obj[AVG] = obj[MAX] = (*iter_sol)->obj;
				eca_imp_cnt[MIN] = eca_imp_cnt[AVG] = eca_imp_cnt[MAX] = (*iter_sol)->eca_imp_cnt;
				ls_iter[MIN] = ls_iter[AVG] = ls_iter[MAX] = (*iter_sol)->ls_iter;
				tm[MIN] = tm[AVG] = tm[MAX] = (*iter_sol)->tm;
				min_obj_cnt = 1;
			}
			else
			{
				obj[AVG] += (*iter_sol)->obj;
				eca_imp_cnt[AVG] += (*iter_sol)->eca_imp_cnt;
				ls_iter[AVG] += (*iter_sol)->ls_iter;
				tm[AVG] += (*iter_sol)->tm;

				if (obj[MIN] - (*iter_sol)->obj > MIN_EQUAL)
				{
					obj[MIN] = (*iter_sol)->obj;
					min_obj_cnt = 1;
				}
				else if (fabs(obj[MIN] - (*iter_sol)->obj )<= MIN_EQUAL)
					min_obj_cnt += 1;
				if ((*iter_sol)->obj - obj[MAX] > MIN_EQUAL)
					obj[MAX] = (*iter_sol)->obj;

				if (eca_imp_cnt[MIN] - (*iter_sol)->eca_imp_cnt > MIN_EQUAL)
					eca_imp_cnt[MIN] = (*iter_sol)->eca_imp_cnt;
				if ((*iter_sol)->eca_imp_cnt - eca_imp_cnt[MAX] > MIN_EQUAL)
					eca_imp_cnt[MAX] = (*iter_sol)->eca_imp_cnt;

				if (ls_iter[MIN] - (*iter_sol)->ls_iter > MIN_EQUAL)
					ls_iter[MIN] = (*iter_sol)->ls_iter;
				if ((*iter_sol)->ls_iter - ls_iter[MAX] > MIN_EQUAL)
					ls_iter[MAX] = (*iter_sol)->ls_iter;

				if (tm[MIN] - (*iter_sol)->tm > MIN_EQUAL)
					tm[MIN] = (*iter_sol)->tm;
				if ((*iter_sol)->tm - tm[MAX] > MIN_EQUAL)
					tm[MAX] = (*iter_sol)->tm;
			}
			/*if (fabs((*iter_sol)->obj - opt_best_map.at((*iter)->filename)) <= MIN_EQUAL)
				hit_cnt += 1;*/
			/*if (opt_best_map.at((*iter)->filename) - (*iter_sol)->obj > MIN_EQUAL)
			{
				system("pause");
			}*/
		}
		hit_cnt = best_obj_cnt;
		if (best_obj_cnt > 0)
			opt_best += 1;
		obj[AVG] /= (*iter)->sol_info_vec.size();
		eca_imp_cnt[AVG] /= (*iter)->sol_info_vec.size();
		ls_iter[AVG] /= (*iter)->sol_info_vec.size();
		tm[AVG] /= (*iter)->sol_info_vec.size();
		if (whe_save_each_file_result)
			ofs << (*iter)->filename << "\t" << (*iter)->n << "\t" << (*iter)->m << "\t"
			<< (*iter)->rand_seed << "\t" << (*iter)->sol_info_vec.size() << "\t"
			<< obj[MIN] << "\t" << obj[AVG] << "\t" << obj[MAX] << "\t"
			<< eca_imp_cnt[MIN] << "\t" << eca_imp_cnt[AVG] << "\t" << eca_imp_cnt[MAX] << "\t"
			<< ls_iter[MIN] << "\t" << ls_iter[AVG] << "\t" << ls_iter[MAX] << "\t"
			<< tm[MIN] << "\t" << tm[AVG] << "\t" << tm[MAX] << "\t"
			<< hit_cnt << endl;
		f += obj[AVG];
		d += (obj[AVG] - obj[MIN]) / obj[MIN];
		h += ((double)hit_cnt / (*iter)->sol_info_vec.size());
		t += tm[AVG];

		set_ins_cnt += 1;
	}
	int sm = 1;
	f /= set_ins_cnt;
	d /= set_ins_cnt, d *= 100;
	h /= set_ins_cnt, h *= 100;
	t /= set_ins_cnt;

	ofst << fnr << "\t"
		<< f << "\t" << d << "\t" << h << "\t" << t << "\t"
		<< opt_best << "\t" << set_ins_cnt << "\t"
		<< endl;
	cout << fnr << "\t"
		<< f << "\t" << d << "\t" << h << "\t" << t << "\t"
		<< opt_best << "\t" << set_ins_cnt << "\t"
		<< endl;

}
#if 0
void Analyze::total_info()
{
	if (whe_input_table_head_fnwt)
		ofst << "fnr \t"
		<< "f[OB] \t d[OB] \t h[OB] \t t[OB] \t opt_best[OB] \t"
		<< "cmp_give[OB][IMPROVED] \t"
		<< "cmp_give[OB][EQUAL] \t"
		<< "cmp_give[OB][NONIMPROVED] \t"
		<< "f[BB] \t d[BB] \t h[BB] \t t[BB] \t opt_best[BB] \t"
		<< "cmp_give[BB][IMPROVED] \t"
		<< "cmp_give[BB][EQUAL] \t"
		<< "cmp_give[BB][NONIMPROVED]"
		<< endl;
	if (whe_save_each_file_result)
		ofs << "(*iter)->filename  \t  (*iter)->n  \t  (*iter)->m  \t"
		<< "(*iter)->opt_obj_given  \t  (*iter)->opt_obj_real  \t  opt_best_map.at((*iter)->filename)  \t"
		<< "(*iter)->rand_seed  \t  (*iter)->sol_info_vec.size()  \t"
		<< "tct[0]  \t  tct[1]  \t  tct[2]  \t"
		<< "obj[0]  \t  obj[1]  \t  obj[2]  \t"
		<< "iteration[0]  \t  iteration[1]  \t  iteration[2]  \t"
		<< "tm[0]  \t  tm[1]  \t  tm[2]  \t"
		<< "cmp_give_result_cnt[0]  \t  cmp_give_result_cnt[1]  \t  cmp_give_result_cnt[2]  \t"
		<< "cmp_real_result_cnt[0]  \t  cmp_real_result_cnt[1]  \t  cmp_real_result_cnt[2]  \t"
		<< "hit_cnt"
		<< endl;
	for (vector<InstanceInfo*>::iterator iter = insinfo_vec.begin();
	iter != insinfo_vec.end(); iter++)
	{
		cout << (*iter)->filename << ", " << (*iter)->m << ", " << (*iter)->n << ", "
			<< (*iter)->opt_obj_given << ", " << (*iter)->opt_obj_real << ", "
			<< (*iter)->rand_seed << ", " << (*iter)->sol_info_vec.size() << endl;
		double obj[3];
		double tct[3];
		double iteration[3];
		double tm[3];
		int hit_cnt = 0;
		int cmp_give_result_cnt[3] = { 0,0,0 }, cmp_real_result_cnt[3] = { 0,0,0 };
		for (vector<SolutionInfo*>::iterator iter_sol = (*iter)->sol_info_vec.begin();
		iter_sol != (*iter)->sol_info_vec.end(); iter_sol++)
		{
			if (iter_sol == (*iter)->sol_info_vec.begin())
			{
				tct[MIN] = tct[AVG] = tct[MAX] = (*iter_sol)->tct;
				obj[MIN] = obj[AVG] = obj[MAX] = (*iter_sol)->obj;
				iteration[MIN] = iteration[AVG] = iteration[MAX] = (*iter_sol)->iteration;
				tm[MIN] = tm[AVG] = tm[MAX] = (*iter_sol)->tm;
			}
			else
			{
				tct[AVG] += (*iter_sol)->tct;
				obj[AVG] += (*iter_sol)->obj;
				iteration[AVG] += (*iter_sol)->iteration;
				tm[AVG] += (*iter_sol)->tm;

				if (tct[MIN] - (*iter_sol)->tct > MIN_EQUAL)
					tct[MIN] = (*iter_sol)->tct;
				if ((*iter_sol)->tct - tct[MAX] > MIN_EQUAL)
					tct[MAX] = (*iter_sol)->tct;

				if (obj[MIN] - (*iter_sol)->obj > MIN_EQUAL)
					obj[MIN] = (*iter_sol)->obj;
				if ((*iter_sol)->obj - obj[MAX] > MIN_EQUAL)
					obj[MAX] = (*iter_sol)->obj;

				if (iteration[MIN] - (*iter_sol)->iteration > MIN_EQUAL)
					iteration[MIN] = (*iter_sol)->iteration;
				if ((*iter_sol)->iteration - iteration[MAX] > MIN_EQUAL)
					iteration[MAX] = (*iter_sol)->iteration;

				if (tm[MIN] - (*iter_sol)->tm > MIN_EQUAL)
					tm[MIN] = (*iter_sol)->tm;
				if ((*iter_sol)->tm - tm[MAX] > MIN_EQUAL)
					tm[MAX] = (*iter_sol)->tm;
			}
			if (fabs((*iter_sol)->obj - opt_best_map.at((*iter)->filename)) <= MIN_EQUAL)
				hit_cnt += 1;
			if (opt_best_map.at((*iter)->filename) - (*iter_sol)->obj > MIN_EQUAL)
			{
				system("pause");
			}
			if ((*iter_sol)->result_improve == 1)
			{
				cmp_real_result_cnt[EQUAL] += 1;// equal
			}
			else if ((*iter_sol)->result_improve == 2)
			{
				cmp_real_result_cnt[IMPROVED] += 1;// improved
			}
			else
				cmp_real_result_cnt[NONIMPROVED] += 1;// not improved

			if ((*iter_sol)->given_result_improve == 1)
			{
				cmp_give_result_cnt[EQUAL] += 1;// equal
			}
			else if ((*iter_sol)->given_result_improve == 2)
			{
				cmp_give_result_cnt[IMPROVED] += 1;// improved
			}
			else
				cmp_give_result_cnt[NONIMPROVED] += 1;// not improved		
		}
		tct[AVG] /= (*iter)->sol_info_vec.size();
		obj[AVG] /= (*iter)->sol_info_vec.size();
		iteration[AVG] /= (*iter)->sol_info_vec.size();
		tm[AVG] /= (*iter)->sol_info_vec.size();
		if (whe_save_each_file_result)
			ofs << (*iter)->filename << "\t" << (*iter)->n << "\t" << (*iter)->m << "\t"
			<< (*iter)->opt_obj_given << "\t" << (*iter)->opt_obj_real << "\t" << opt_best_map.at((*iter)->filename) << "\t"
			<< (*iter)->rand_seed << "\t" << (*iter)->sol_info_vec.size() << "\t"
			<< tct[MIN] << "\t" << tct[AVG] << "\t" << tct[MAX] << "\t"
			<< obj[MIN] << "\t" << obj[AVG] << "\t" << obj[MAX] << "\t"
			<< iteration[MIN] << "\t" << iteration[AVG] << "\t" << iteration[MAX] << "\t"
			<< tm[MIN] << "\t" << tm[AVG] << "\t" << tm[MAX] << "\t"
			<< cmp_give_result_cnt[IMPROVED] << "\t" << cmp_give_result_cnt[EQUAL] << "\t" << cmp_give_result_cnt[NONIMPROVED] << "\t"
			<< cmp_real_result_cnt[IMPROVED] << "\t" << cmp_real_result_cnt[EQUAL] << "\t" << cmp_real_result_cnt[NONIMPROVED] << "\t"
			<< hit_cnt << endl;
		int set_index = OB;
		if ((*iter)->n > 14)	// OB set
			set_index = BB;
		f[set_index] += obj[AVG];
		d[set_index] += ((obj[AVG] - obj[MIN]) / obj[MIN]);
		double cur_instance_hr = (double)cmp_give_result_cnt[set_index == OB ? EQUAL : IMPROVED] / (*iter)->sol_info_vec.size();
		//h[set_index] += cur_instance_hr;
		h[set_index] += ((double)hit_cnt / (*iter)->sol_info_vec.size());
		t[set_index] += tm[AVG];
		opt_best[set_index] += opt_best_map.at((*iter)->filename);
		if (cmp_give_result_cnt[IMPROVED] > 0)
			cmp_give_result_total_cnt[set_index][IMPROVED] += 1;
		if (cmp_give_result_cnt[IMPROVED] == 0 && cmp_give_result_cnt[EQUAL] > 0)
			cmp_give_result_total_cnt[set_index][EQUAL] += 1;
		if (cmp_give_result_cnt[IMPROVED] == 0 && cmp_give_result_cnt[EQUAL] == 0)
			cmp_give_result_total_cnt[set_index][NONIMPROVED] += 1;

		int p_para, d_para, m_para, n_para;
		p_para = (*iter)->p == 1 ? P1 : P2;
		d_para = (*iter)->d == 1 ? D1 : D2;
		if ((set_index == OB && (*iter)->m == 2) || (set_index == BB && (*iter)->m == 4))
			m_para = M1;
		else if ((set_index == OB && (*iter)->m == 3) || (set_index == BB && (*iter)->m == 7))
			m_para = M2;
		else
			m_para = M3;
		if ((set_index == OB && (*iter)->n == 8) || (set_index == BB && (*iter)->n == 20))
			n_para = N1;
		else if ((set_index == OB && (*iter)->n == 11) || (set_index == BB && (*iter)->n == 35))
			n_para = N2;
		else
			n_para = N3;
		optfound_pct[set_index][p_para][HR_JUDGE] += cur_instance_hr;
		optfound_pct[set_index][d_para][HR_JUDGE] += cur_instance_hr;
		optfound_pct[set_index][m_para][HR_JUDGE] += cur_instance_hr;
		optfound_pct[set_index][n_para][HR_JUDGE] += cur_instance_hr;
		optfound_pct[set_index][p_para][DEV_JUDGE] += (obj[AVG] - obj[MIN]) / obj[MIN];
		optfound_pct[set_index][d_para][DEV_JUDGE] += (obj[AVG] - obj[MIN]) / obj[MIN];
		optfound_pct[set_index][m_para][DEV_JUDGE] += (obj[AVG] - obj[MIN]) / obj[MIN];
		optfound_pct[set_index][n_para][DEV_JUDGE] += (obj[AVG] - obj[MIN]) / obj[MIN];
		optfound_pct[set_index][p_para][TIME_JUDGE] += tm[AVG];
		optfound_pct[set_index][d_para][TIME_JUDGE] += tm[AVG];
		optfound_pct[set_index][m_para][TIME_JUDGE] += tm[AVG];
		optfound_pct[set_index][n_para][TIME_JUDGE] += tm[AVG];
		optfound_pct[set_index][p_para][OBJ_JUDGE] += obj[AVG];
		optfound_pct[set_index][d_para][OBJ_JUDGE] += obj[AVG];
		optfound_pct[set_index][m_para][OBJ_JUDGE] += obj[AVG];
		optfound_pct[set_index][n_para][OBJ_JUDGE] += obj[AVG];

		set_ins_cnt[set_index] += 1;
	}
	for (int sm = SM::OB; sm <= SM::BB; sm++)
	{
		f[sm] /= set_ins_cnt[sm];
		d[sm] /= set_ins_cnt[sm], d[sm] *= 100;
		h[sm] /= set_ins_cnt[sm], h[sm] *= 100;
		t[sm] /= set_ins_cnt[sm];
		opt_best[sm] /= set_ins_cnt[sm];
		for (int para = PM::P1; para <= PM::N3; para++)
		{
			for (int judm = JM::HR_JUDGE; judm <= JM::OBJ_JUDGE; judm++)
			{
				if (para <= PM::D2)
					optfound_pct[sm][para][judm] /= (set_ins_cnt[sm] / 2);
				else
					optfound_pct[sm][para][judm] /= (set_ins_cnt[sm] / 3);
				if (judm == JM::DEV_JUDGE || judm == JM::HR_JUDGE)
					optfound_pct[sm][para][judm] *= 100;
			}
		}
	}
	/*ofst << "fnr \t"
	<< "optfound_pct[OB][P1] \t optfound_pct[OB][P2] \t"
	<< "optfound_pct[OB][D1] \t optfound_pct[OB][D2] \t"
	<< "optfound_pct[OB][M1] \t optfound_pct[OB][M2] \t optfound_pct[OB][M3] \t"
	<< "optfound_pct[OB][N1] \t optfound_pct[OB][N2] \t optfound_pct[OB][N3] \t"
	<< "optfound_pct[BB][P1] \t optfound_pct[BB][P2] \t"
	<< "optfound_pct[BB][D1] \t optfound_pct[BB][D2] \t"
	<< "optfound_pct[BB][M1] \t optfound_pct[BB][M2] \t optfound_pct[BB][M3] \t"
	<< "optfound_pct[BB][N1] \t optfound_pct[BB][N2] \t optfound_pct[BB][N3] \t"
	<< endl;*/
	for (int judm = JM::HR_JUDGE; judm <= JM::OBJ_JUDGE; judm++)
	{
		string str_jm = "HR_JUDGE";
		if (judm == JM::DEV_JUDGE)
			str_jm = "DEV_JUDGE";
		if (judm == JM::OBJ_JUDGE)
			str_jm = "OBJ_JUDGE";
		if (judm == JM::TIME_JUDGE)
			str_jm = "TIME_JUDGE";
		/*ofst << fnr << "\t" << str_jm << "\t"
		<< optfound_pct[OB][P1][judm] << "\t" << optfound_pct[OB][P2][judm] << "\t"
		<< optfound_pct[OB][D1][judm] << "\t" << optfound_pct[OB][D2][judm] << "\t"
		<< optfound_pct[OB][M1][judm] << "\t" << optfound_pct[OB][M2][judm] << "\t" << optfound_pct[OB][M3][judm] << "\t"
		<< optfound_pct[OB][M1][judm] << "\t" << optfound_pct[OB][M2][judm] << "\t" << optfound_pct[OB][M3][judm] << "\t"
		<< optfound_pct[BB][P1][judm] << "\t" << optfound_pct[BB][P2][judm] << "\t"
		<< optfound_pct[BB][D1][judm] << "\t" << optfound_pct[BB][D2][judm] << "\t"
		<< optfound_pct[BB][M1][judm] << "\t" << optfound_pct[BB][M2][judm] << "\t" << optfound_pct[BB][M3][judm] << "\t"
		<< optfound_pct[BB][M1][judm] << "\t" << optfound_pct[BB][M2][judm] << "\t" << optfound_pct[BB][M3][judm] << "\t"
		<< endl;*/
	}

	ofst << fnr << "\t"
		<< f[OB] << "\t" << d[OB] << "\t" << h[OB] << "\t" << t[OB] << "\t" << opt_best[OB] << "\t"
		<< cmp_give_result_total_cnt[OB][IMPROVED] << "\t"
		<< cmp_give_result_total_cnt[OB][EQUAL] << "\t"
		<< cmp_give_result_total_cnt[OB][NONIMPROVED] << "\t"
		<< f[BB] << "\t" << d[BB] << "\t" << h[BB] << "\t" << t[BB] << "\t" << opt_best[BB] << "\t"
		<< cmp_give_result_total_cnt[BB][IMPROVED] << "\t"
		<< cmp_give_result_total_cnt[BB][EQUAL] << "\t"
		<< cmp_give_result_total_cnt[BB][NONIMPROVED]
		<< endl;
	cout << fnr << "\t"
		<< f[OB] << "\t" << d[OB] << "\t" << h[OB] << "\t" << t[OB] << "\t" << opt_best[OB] << "\t"
		<< cmp_give_result_total_cnt[OB][IMPROVED] << "\t"
		<< cmp_give_result_total_cnt[OB][EQUAL] << "\t"
		<< cmp_give_result_total_cnt[OB][NONIMPROVED] << "\t"
		<< f[BB] << "\t" << d[BB] << "\t" << h[BB] << "\t" << t[BB] << "\t" << opt_best[BB] << "\t"
		<< cmp_give_result_total_cnt[BB][IMPROVED] << "\t"
		<< cmp_give_result_total_cnt[BB][EQUAL] << "\t"
		<< cmp_give_result_total_cnt[BB][NONIMPROVED]
		<< endl;
}
void analyze_total_file(string fnr, string fnw)
{
	ifstream ifs(fnr);
	if (!ifs.is_open())
	{
		cout << fnr << endl; perror("file_input.");
		exit(0);
	}
	ofstream ofs(fnw, ios::trunc | ios::out);
	if (!ofs.is_open())
	{
		cout << fnw << endl; perror("file_output.");
		exit(0);
	}
	string strline;
	vector<string> fields_vec, ins_name_vec;
	vector<InstanceInfo*>insinfo_vec;
	vector<SolutionInfo*>  *sol_info_vec = new vector<SolutionInfo*>;
	while (getline(ifs, strline))
	{
		//	cout << strline << endl;
		InstanceInfo *insinfo;

		string filename;
		int m, n, rand_seed;
		double opt_obj_given, opt_obj_real;

		//split(fields_vec, strline, is_any_of("\t"));
		split_generic<string>(fields_vec, strline, "\t");
		/*for (vector<string>::iterator iter = fields_vec.begin();
		iter != fields_vec.end(); iter++)
		{
		cout << *iter << "\t";

		}*/
		if (fields_vec.front().find("Ni_") != string::npos)
		{
			//split(ins_name_vec, fields_vec[0], is_any_of("_-"));
			split_generic<string>(ins_name_vec, fields_vec[0], "_-");
			sol_info_vec = new vector<SolutionInfo*>;
			insinfo = new InstanceInfo(fields_vec[0],
				stoi(fields_vec[1]), stoi(fields_vec[2]),
				stoi(ins_name_vec[3]), stoi(ins_name_vec[4]), stoi(ins_name_vec[5]),
				stod(fields_vec[3]),
				stod(fields_vec[4]), stoi(fields_vec[5]),
				*sol_info_vec);
			insinfo_vec.push_back(insinfo);
		}
		else
		{
			sol_info_vec->push_back(new SolutionInfo(stoi(fields_vec[0]),//run_cnt
				stod(fields_vec[1]), stod(fields_vec[2]), stoi(fields_vec[3]),//obj,mm
				stoi(fields_vec[4]), stoi(fields_vec[5]),//iter, tm
				stoi(fields_vec[6]), stoi(fields_vec[7])//given_result_improve,result_improve
				));
			//cout << stoi(fields_vec[0]) << "\t"//run_cnt
			//	<< boost::lexical_cast<double>(fields_vec[1]) << "\t" << stoi(fields_vec[2]) << "\t"//obj<<"\t"mm
			//	<< stoi(fields_vec[3]) << "\t" << stoi(fields_vec[4]) << "\t"//iter<<"\t" tm
			//	<< stoi(fields_vec[5]) << "\t" << stoi(fields_vec[6]) << endl;
		}
		//cout << endl;
	}
	ofs << "(*iter)->filename  \t  (*iter)->n  \t  (*iter)->m  \t"
		<< "(*iter)->opt_obj_given  \t  (*iter)->opt_obj_real  \t"
		<< "(*iter)->rand_seed  \t  (*iter)->sol_info_vec.size()  \t"
		<< "obj[0]  \t  obj[1]  \t  obj[2]  \t"
		<< "iteration[0]  \t  iteration[1]  \t  iteration[2]  \t"
		<< "tm[0]  \t  tm[1]  \t  tm[2]  \t"
		<< "cmp_give_result_cnt[0]  \t  cmp_give_result_cnt[1]  \t  cmp_give_result_cnt[2]  \t"
		<< "cmp_give_result[0]  \t  cmp_give_result[1]  \t  cmp_give_result[2]  \t"
		<< "cmp_real_result_cnt[0]  \t  cmp_real_result_cnt[1]  \t  cmp_real_result_cnt[2]  \t"
		<< "cmp_real_result[0]  \t  cmp_real_result[1]  \t  cmp_real_result[2]"
		<< endl;

	for (vector<InstanceInfo*>::iterator iter = insinfo_vec.begin();
	iter != insinfo_vec.end(); iter++)
	{
		cout << (*iter)->filename << ", " << (*iter)->m << ", " << (*iter)->n << ", "
			<< (*iter)->opt_obj_given << ", " << (*iter)->opt_obj_real << ", "
			<< (*iter)->rand_seed << ", " << (*iter)->sol_info_vec.size() << endl;
		double obj[3];
		int iteration[3];
		int tm[3];
		int //cmp_give_result[3] = { 0,0,0 }, cmp_real_result[3] = { 0,0,0 }, 
			cmp_give_result_cnt[3] = { 0,0,0 }, cmp_real_result_cnt[3] = { 0,0,0 };
		//cmp_give_result_t,com_real_result_t;
		for (vector<SolutionInfo*>::iterator iter_sol = (*iter)->sol_info_vec.begin();
		iter_sol != (*iter)->sol_info_vec.end(); iter_sol++)
		{
			//ofs << run_cnt << "\t"
			//	//<< c[si_opt][0] << "\t"
			//	<< c[si][0] << "\t"
			//	<< mm[si] << "\t"
			//	<< iterration << "\t"
			//	<< (end_tm - start_tm) /*/ CLOCKS_PER_SEC*/ << "\t"
			//	<< given_result_improve << "\t"
			//	<< result_improve
			//	<< endl;
			if (iter_sol == (*iter)->sol_info_vec.begin())
			{
				obj[0] = obj[1] = obj[2] = (*iter_sol)->obj;
				iteration[0] = iteration[1] = iteration[2] = (*iter_sol)->iteration;
				tm[0] = tm[1] = tm[2] = (*iter_sol)->tm;
			}
			else
			{
				obj[1] += (*iter_sol)->obj;
				iteration[1] += (*iter_sol)->iteration;
				tm[1] += (*iter_sol)->tm;

				if (obj[0] - (*iter_sol)->obj > MIN_EQUAL)
					obj[0] = (*iter_sol)->obj;
				if ((*iter_sol)->obj - obj[2] > MIN_EQUAL)
					obj[2] = (*iter_sol)->obj;

				if (iteration[0] - (*iter_sol)->iteration > MIN_EQUAL)
					iteration[0] = (*iter_sol)->iteration;
				if ((*iter_sol)->iteration - iteration[2] > MIN_EQUAL)
					iteration[2] = (*iter_sol)->iteration;

				if (tm[0] - (*iter_sol)->tm > MIN_EQUAL)
					tm[0] = (*iter_sol)->tm;
				if ((*iter_sol)->tm - tm[2] > MIN_EQUAL)
					tm[2] = (*iter_sol)->tm;
			}

			if ((*iter_sol)->result_improve == 1)
			{
				cmp_real_result_cnt[1] += 1;// equal
			}
			else if ((*iter_sol)->result_improve == 2)
			{
				cmp_real_result_cnt[0] += 1;// improved
			}
			else
				cmp_real_result_cnt[2] += 1;// not improved

			if ((*iter_sol)->given_result_improve == 1)
			{
				cmp_give_result_cnt[1] += 1;// equal
			}
			else if ((*iter_sol)->given_result_improve == 2)
			{
				cmp_give_result_cnt[0] += 1;// improved
			}
			else
				cmp_give_result_cnt[2] += 1;// not improved		
		}
		obj[1] /= (*iter)->sol_info_vec.size();
		iteration[1] /= (*iter)->sol_info_vec.size();
		tm[1] /= (*iter)->sol_info_vec.size();
		ofs << (*iter)->filename << "\t" << (*iter)->n << "\t" << (*iter)->m << "\t"
			<< (*iter)->opt_obj_given << "\t" << (*iter)->opt_obj_real << "\t"
			<< (*iter)->rand_seed << "\t" << (*iter)->sol_info_vec.size() << "\t"
			<< obj[0] << "\t" << obj[1] << "\t" << obj[2] << "\t"
			<< iteration[0] << "\t" << iteration[1] << "\t" << iteration[2] << "\t"
			<< tm[0] << "\t" << tm[1] << "\t" << tm[2] << "\t"
			<< cmp_give_result_cnt[0] << "\t" << cmp_give_result_cnt[1] << "\t" << cmp_give_result_cnt[2] << "\t"
			//<< cmp_give_result[0] << "\t" << cmp_give_result[1] << "\t" << cmp_give_result[2] << "\t"
			<< cmp_real_result_cnt[0] << "\t" << cmp_real_result_cnt[1] << "\t" << cmp_real_result_cnt[2] << endl;
		//<< cmp_real_result[0] << "\t" << cmp_real_result[1] << "\t" << cmp_real_result[2] << endl;

}
	ifs.close();
}
#endif
int main(int argc, char **argv)
{
	char *rgv_win[] = { "",	//0
		"_fn", "ECA_am1_p20_itr2000_ptr10_ecag70_cl5_r1_r10",	//1,2
		"_if", "instance\\BB_Problem_BestSolution\\",	//3,4	
		"_of", "results\\",	// output file path
		"_efn", "_analyze0",	// each output file name, append to the source file, trunc mode
		"_tfn", "ECA_T0", // total information of all the files, append mode
		"_optfn", "pms_twc_best",	// best found solution values file name
		"_whe_save_each","1",		// whether save the results of each file	
		"_total_info_per_para", "0",	// invoke total_information function or per parameter function
		"_sfix",".txt"	// suffix of input file
	};
	std::map<string, string> argv_map;
	argv_map["_exe_name"] = argv[0];
	for (int i = 1; i < sizeof(rgv_win) / sizeof(rgv_win[0]); i += 2)
		argv_map[string(rgv_win[i])] = string(rgv_win[i + 1]);
	for (int i = 1; i < argc; i += 2)
		argv_map[string(argv[i])] = string(argv[i + 1]);
#ifdef __linux__
	replace_all(argv_map["_if"], "\\", "//");
	replace_all(argv_map["_of"], "\\", "//");
#endif
	string fnr, fnw, fnwt, fnr_best;
	fnr = argv_map.at("_of") + argv_map.at("_fn") + argv_map.at("_sfix");
	fnw = argv_map.at("_of") + argv_map.at("_fn") + argv_map.at("_efn") + argv_map.at("_sfix");
	fnwt = argv_map.at("_of") + argv_map.at("_tfn") + argv_map.at("_sfix");
	fnr_best = argv_map.at("_of") + argv_map.at("_optfn") + argv_map.at("_sfix");
	Analyze *an = new Analyze(fnr, fnw, fnwt, fnr_best, stoi(argv_map.at("_whe_save_each")));
	if (stoi(argv_map.at("_total_info_per_para")) == 0)
		an->para_setting();
	/*else if (stoi(argv_map.at("_total_info_per_para")) == 1)
		an->total_info();*/
		//analyze_total_file(fnr,fnw);
	delete an;
#ifdef _WIN32
	system("pause");
#endif	
}
#endif