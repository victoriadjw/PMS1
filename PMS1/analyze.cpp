#if 0
#include <boost/lambda/lambda.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<map>
using namespace std;
using namespace boost;
const double MIN_EQUAL = 0.001;
class SolutionInfo
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
public:
	int run_cnt, mm, iteration, tm, given_result_improve, result_improve;
	double obj;
	SolutionInfo(int _rc, double _obj, int _mm, int _iter, int _tm, int _gri, int _ri)
		:run_cnt(_rc), obj(_obj), mm(_mm), iteration(_iter), tm(_tm), given_result_improve(_gri), result_improve(_ri)
	{}

};
class InstanceInfo
{
public:
	string filename;
	int m, n, p,d,ins,rand_seed;
	double opt_obj_given, opt_obj_real;
	vector<SolutionInfo*> &sol_info_vec;
	InstanceInfo(string _fn,int _n,int _m,int _p, int _d, int _ins,double _oog,double _oor,int _rd,vector<SolutionInfo*> &_siv)
		:filename(_fn),m(_m),n(_n),p(_p),d(_d),ins(_ins),opt_obj_given(_oog),opt_obj_real(_oor),rand_seed(_rd),sol_info_vec(_siv)
	{}
	
};
class Analyze {
public:
	Analyze(string, string,string);
	void total_info();
	void para_setting();
	double f[2], d[2], h[2], t[2];
	int cmp_give_result_total_cnt[2][3], set_ins_cnt[2];
	double optfound_pct[2][10][4];
	const enum SM {OB,BB};
	const enum CM {IMPROVED,EQUAL,NONIMPROVED};
	const enum BM {MIN,AVG,MAX};
	const enum PM {P1,P2,D1,D2,M1,M2,M3,N1,N2,N3};
	const enum JM { HR_JUDGE, DEV_JUDGE, TIME_JUDGE, OBJ_JUDGE };
private:
	string fnr,fnw,fnwt;
	bool whe_input_table_head;
	ofstream ofs,ofst;
	vector<InstanceInfo*>insinfo_vec; 
};
Analyze::Analyze(string _fnr, string _fnw,string _fnwt):fnr(_fnr),fnw(_fnw),fnwt(_fnwt)
{
	ifstream ifs(fnwt);
	if (!ifs.is_open())
		whe_input_table_head = true;
	else
		whe_input_table_head = false;
	ifs.close();
	ifs.open(fnr);
	if (!ifs.is_open())
	{
		cout << fnr << endl; perror("file_input.");
		exit(0);
	}
	ofs.open(fnw, ios::trunc | ios::out);
	if (!ofs.is_open())
	{
		cout << fnw << endl; perror("file_output.");
		exit(0);
	}
	ofs.setf(ios::fixed, ios::floatfield);
	ofs.precision(6);
	ofst.open(fnwt, ios::app | ios::out);
	if (!ofst.is_open())
	{
		cout << fnwt << endl; perror("file_output.");
		exit(0);
	}
	ofst.setf(ios::fixed, ios::floatfield);
	ofst.precision(6);
	string strline;
	vector<string> fields_vec,ins_name_vec;
	//vector<InstanceInfo*>insinfo_vec;
	vector<SolutionInfo*>  *sol_info_vec = new vector<SolutionInfo*>;
	while (getline(ifs, strline))
	{
		//	cout << strline << endl;
		InstanceInfo *insinfo;

		string filename;
		int m, n, rand_seed;
		double opt_obj_given, opt_obj_real;

		split(fields_vec, strline, is_any_of("\t"));
		/*for (vector<string>::iterator iter = fields_vec.begin();
		iter != fields_vec.end(); iter++)
		{
		cout << *iter << "\t";

		}*/
		if (fields_vec.front().find("Ni_") != string::npos)
		{
			split(ins_name_vec, fields_vec[0], is_any_of("_-"));
			sol_info_vec = new vector<SolutionInfo*>;
			insinfo = new InstanceInfo(fields_vec[0],
				stoi(fields_vec[1]), stoi(fields_vec[2]),
				stoi(ins_name_vec[3]), stoi(ins_name_vec[4]), stoi(ins_name_vec[5]),
				boost::lexical_cast<double>(fields_vec[3]),
				boost::lexical_cast<double>(fields_vec[4]), stoi(fields_vec[5]),
				*sol_info_vec);
			insinfo_vec.push_back(insinfo);
		}
		else
		{
			sol_info_vec->push_back(new SolutionInfo(stoi(fields_vec[0]),//run_cnt
				boost::lexical_cast<double>(fields_vec[1]), stoi(fields_vec[2]),//obj,mm
				stoi(fields_vec[3]), stoi(fields_vec[4]),//iter, tm
				stoi(fields_vec[5]), stoi(fields_vec[6])//given_result_improve,result_improve
				));
			//cout << stoi(fields_vec[0]) << "\t"//run_cnt
			//	<< boost::lexical_cast<double>(fields_vec[1]) << "\t" << stoi(fields_vec[2]) << "\t"//obj<<"\t"mm
			//	<< stoi(fields_vec[3]) << "\t" << stoi(fields_vec[4]) << "\t"//iter<<"\t" tm
			//	<< stoi(fields_vec[5]) << "\t" << stoi(fields_vec[6]) << endl;
		}
		//cout << endl;
	}
	ifs.close();
}
void Analyze::para_setting()
{
	if (whe_input_table_head)
		ofst << "fnr \t"
		<< "(*iter)->filename  \t  (*iter)->n  \t  (*iter)->m  \t"
		<< "(*iter)->opt_obj_given  \t  (*iter)->opt_obj_real  \t"
		<< "(*iter)->rand_seed  \t  (*iter)->sol_info_vec.size()  \t"
		<< "obj[0]  \t  obj[1]  \t  obj[2]  \t"
		<< "iteration[0]  \t  iteration[1]  \t  iteration[2]  \t"
		<< "tm[0]  \t  tm[1]  \t  tm[2]  \t"
		<< "cmp_give_result_cnt[0]  \t  cmp_give_result_cnt[1]  \t  cmp_give_result_cnt[2]  \t"
		<< "cmp_real_result_cnt[0]  \t  cmp_real_result_cnt[1]  \t  cmp_real_result_cnt[2]  \t"
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
		int cmp_give_result_cnt[3] = { 0,0,0 }, cmp_real_result_cnt[3] = { 0,0,0 };
		for (vector<SolutionInfo*>::iterator iter_sol = (*iter)->sol_info_vec.begin();
		iter_sol != (*iter)->sol_info_vec.end(); iter_sol++)
		{
			if (iter_sol == (*iter)->sol_info_vec.begin())
			{
				obj[MIN] = obj[AVG] = obj[MAX] = (*iter_sol)->obj;
				iteration[MIN] = iteration[AVG] = iteration[MAX] = (*iter_sol)->iteration;
				tm[MIN] = tm[AVG] = tm[MAX] = (*iter_sol)->tm;
			}
			else
			{
				obj[AVG] += (*iter_sol)->obj;
				iteration[AVG] += (*iter_sol)->iteration;
				tm[AVG] += (*iter_sol)->tm;

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
		obj[AVG] /= (*iter)->sol_info_vec.size();
		iteration[AVG] /= (*iter)->sol_info_vec.size();
		tm[AVG] /= (*iter)->sol_info_vec.size();
		ofst << fnr << "\t" << (*iter)->filename << "\t" << (*iter)->n << "\t" << (*iter)->m << "\t"
			<< (*iter)->opt_obj_given << "\t" << (*iter)->opt_obj_real << "\t"
			<< (*iter)->rand_seed << "\t" << (*iter)->sol_info_vec.size() << "\t"
			<< obj[MIN] << "\t" << obj[AVG] << "\t" << obj[MAX] << "\t"
			<< iteration[MIN] << "\t" << iteration[AVG] << "\t" << iteration[MAX] << "\t"
			<< tm[MIN] << "\t" << tm[AVG] << "\t" << tm[MAX] << "\t"
			<< cmp_give_result_cnt[IMPROVED] << "\t" << cmp_give_result_cnt[EQUAL] << "\t" << cmp_give_result_cnt[NONIMPROVED] << "\t"
			<< cmp_real_result_cnt[IMPROVED] << "\t" << cmp_real_result_cnt[EQUAL] << "\t" << cmp_real_result_cnt[NONIMPROVED] << endl;
		if (cmp_give_result_cnt[IMPROVED] > 0)
		{
			/*ofst << (*iter)->filename << "\t" << (*iter)->n << "\t" << (*iter)->m << "\t"
			<< (*iter)->opt_obj_given << "\t" << (*iter)->opt_obj_real << "\t"
			<< (*iter)->rand_seed << "\t" << (*iter)->sol_info_vec.size() << "\t"
			<< obj[MIN] << "\t" << obj[AVG] << "\t" << obj[MAX] << "\t"
			<< iteration[MIN] << "\t" << iteration[AVG] << "\t" << iteration[MAX] << "\t"
			<< tm[MIN] << "\t" << tm[AVG] << "\t" << tm[MAX] << "\t"
			<< cmp_give_result_cnt[IMPROVED] << "\t" << cmp_give_result_cnt[EQUAL] << "\t" << cmp_give_result_cnt[NONIMPROVED] << "\t"
			<< cmp_real_result_cnt[IMPROVED] << "\t" << cmp_real_result_cnt[EQUAL] << "\t" << cmp_real_result_cnt[NONIMPROVED] << endl;*/
		}
		int set_index = OB;
		if ((*iter)->n > 14)	// OB set
			set_index = BB;
		f[set_index] += obj[AVG];
		d[set_index] += (obj[AVG] - obj[MIN]) / obj[MIN];
		double cur_instance_hr = (double)cmp_give_result_cnt[set_index == OB ? EQUAL : IMPROVED] / (*iter)->sol_info_vec.size();
		h[set_index] += cur_instance_hr;
		t[set_index] += tm[AVG];
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
	/*ofst << "fnr \t"
	<< "f[OB] \t d[OB] \t h[OB] \t t[OB] \t"
	<< "cmp_give[OB][IMPROVED] \t"
	<< "cmp_give[OB][EQUAL] \t"
	<< "cmp_give[OB][NONIMPROVED] \t"
	<< "f[BB] \t d[BB] \t h[BB] \t t[BB] \t"
	<< "cmp_give[BB][IMPROVED] \t"
	<< "cmp_give[BB][EQUAL] \t"
	<< "cmp_give[BB][NONIMPROVED]"
	<< endl;*/
	cout << fnr << "\t"
		<< f[OB] << "\t" << d[OB] << "\t" << h[OB] << "\t" << t[OB] << "\t"
		<< cmp_give_result_total_cnt[OB][IMPROVED] << "\t"
		<< cmp_give_result_total_cnt[OB][EQUAL] << "\t"
		<< cmp_give_result_total_cnt[OB][NONIMPROVED] << "\t"
		<< f[BB] << "\t" << d[BB] << "\t" << h[BB] << "\t" << t[BB] << "\t"
		<< cmp_give_result_total_cnt[BB][IMPROVED] << "\t"
		<< cmp_give_result_total_cnt[BB][EQUAL] << "\t"
		<< cmp_give_result_total_cnt[BB][NONIMPROVED]
		<< endl;
	cout << fnr << "\t"
		<< f[OB] << "\t" << d[OB] << "\t" << h[OB] << "\t" << t[OB] << "\t"
		<< cmp_give_result_total_cnt[OB][IMPROVED] << "\t"
		<< cmp_give_result_total_cnt[OB][EQUAL] << "\t"
		<< cmp_give_result_total_cnt[OB][NONIMPROVED] << "\t"
		<< f[BB] << "\t" << d[BB] << "\t" << h[BB] << "\t" << t[BB] << "\t"
		<< cmp_give_result_total_cnt[BB][IMPROVED] << "\t"
		<< cmp_give_result_total_cnt[BB][EQUAL] << "\t"
		<< cmp_give_result_total_cnt[BB][NONIMPROVED]
		<< endl;

}
void Analyze::total_info()
{
	if (whe_input_table_head)
		ofst << "fnr \t"
		<< "f[OB] \t d[OB] \t h[OB] \t t[OB] \t"
		<< "cmp_give[OB][IMPROVED] \t"
		<< "cmp_give[OB][EQUAL] \t"
		<< "cmp_give[OB][NONIMPROVED] \t"
		<< "f[BB] \t d[BB] \t h[BB] \t t[BB] \t"
		<< "cmp_give[BB][IMPROVED] \t"
		<< "cmp_give[BB][EQUAL] \t"
		<< "cmp_give[BB][NONIMPROVED]"
		<< endl;
	ofs << "(*iter)->filename  \t  (*iter)->n  \t  (*iter)->m  \t"
		<< "(*iter)->opt_obj_given  \t  (*iter)->opt_obj_real  \t"
		<< "(*iter)->rand_seed  \t  (*iter)->sol_info_vec.size()  \t"
		<< "obj[0]  \t  obj[1]  \t  obj[2]  \t"
		<< "iteration[0]  \t  iteration[1]  \t  iteration[2]  \t"
		<< "tm[0]  \t  tm[1]  \t  tm[2]  \t"
		<< "cmp_give_result_cnt[0]  \t  cmp_give_result_cnt[1]  \t  cmp_give_result_cnt[2]  \t"
		<< "cmp_real_result_cnt[0]  \t  cmp_real_result_cnt[1]  \t  cmp_real_result_cnt[2]  \t"
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
		int cmp_give_result_cnt[3] = { 0,0,0 }, cmp_real_result_cnt[3] = { 0,0,0 };
		for (vector<SolutionInfo*>::iterator iter_sol = (*iter)->sol_info_vec.begin();
		iter_sol != (*iter)->sol_info_vec.end(); iter_sol++)
		{
			if (iter_sol == (*iter)->sol_info_vec.begin())
			{
				obj[MIN] = obj[AVG] = obj[MAX] = (*iter_sol)->obj;
				iteration[MIN] = iteration[AVG] = iteration[MAX] = (*iter_sol)->iteration;
				tm[MIN] = tm[AVG] = tm[MAX] = (*iter_sol)->tm;
			}
			else
			{
				obj[AVG] += (*iter_sol)->obj;
				iteration[AVG] += (*iter_sol)->iteration;
				tm[AVG] += (*iter_sol)->tm;

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
		obj[AVG] /= (*iter)->sol_info_vec.size();
		iteration[AVG] /= (*iter)->sol_info_vec.size();
		tm[AVG] /= (*iter)->sol_info_vec.size();
		ofs << (*iter)->filename << "\t" << (*iter)->n << "\t" << (*iter)->m << "\t"
			<< (*iter)->opt_obj_given << "\t" << (*iter)->opt_obj_real << "\t"
			<< (*iter)->rand_seed << "\t" << (*iter)->sol_info_vec.size() << "\t"
			<< obj[MIN] << "\t" << obj[AVG] << "\t" << obj[MAX] << "\t"
			<< iteration[MIN] << "\t" << iteration[AVG] << "\t" << iteration[MAX] << "\t"
			<< tm[MIN] << "\t" << tm[AVG] << "\t" << tm[MAX] << "\t"
			<< cmp_give_result_cnt[IMPROVED] << "\t" << cmp_give_result_cnt[EQUAL] << "\t" << cmp_give_result_cnt[NONIMPROVED] << "\t"
			<< cmp_real_result_cnt[IMPROVED] << "\t" << cmp_real_result_cnt[EQUAL] << "\t" << cmp_real_result_cnt[NONIMPROVED] << endl;
		if (cmp_give_result_cnt[IMPROVED] > 0)
		{
			/*ofst << (*iter)->filename << "\t" << (*iter)->n << "\t" << (*iter)->m << "\t"
				<< (*iter)->opt_obj_given << "\t" << (*iter)->opt_obj_real << "\t"
				<< (*iter)->rand_seed << "\t" << (*iter)->sol_info_vec.size() << "\t"
				<< obj[MIN] << "\t" << obj[AVG] << "\t" << obj[MAX] << "\t"
				<< iteration[MIN] << "\t" << iteration[AVG] << "\t" << iteration[MAX] << "\t"
				<< tm[MIN] << "\t" << tm[AVG] << "\t" << tm[MAX] << "\t"
				<< cmp_give_result_cnt[IMPROVED] << "\t" << cmp_give_result_cnt[EQUAL] << "\t" << cmp_give_result_cnt[NONIMPROVED] << "\t"
				<< cmp_real_result_cnt[IMPROVED] << "\t" << cmp_real_result_cnt[EQUAL] << "\t" << cmp_real_result_cnt[NONIMPROVED] << endl;*/
		}
		int set_index = OB;
		if ((*iter)->n > 14)	// OB set
			set_index = BB;
		f[set_index] += obj[AVG];
		d[set_index] += (obj[AVG] - obj[MIN]) / obj[MIN];
		double cur_instance_hr = (double)cmp_give_result_cnt[set_index == OB ? EQUAL : IMPROVED] / (*iter)->sol_info_vec.size();
		h[set_index] += cur_instance_hr;
		t[set_index] += tm[AVG];
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
		<< f[OB] << "\t" << d[OB] << "\t" << h[OB] << "\t" << t[OB] << "\t"
		<< cmp_give_result_total_cnt[OB][IMPROVED] << "\t"
		<< cmp_give_result_total_cnt[OB][EQUAL] << "\t"
		<< cmp_give_result_total_cnt[OB][NONIMPROVED] << "\t"
		<< f[BB] << "\t" << d[BB] << "\t" << h[BB] << "\t" << t[BB] << "\t"
		<< cmp_give_result_total_cnt[BB][IMPROVED] << "\t"
		<< cmp_give_result_total_cnt[BB][EQUAL] << "\t"
		<< cmp_give_result_total_cnt[BB][NONIMPROVED]
		<< endl;
	cout << fnr << "\t"
		<< f[OB] << "\t" << d[OB] << "\t" << h[OB] << "\t" << t[OB] << "\t"
		<< cmp_give_result_total_cnt[OB][IMPROVED] << "\t"
		<< cmp_give_result_total_cnt[OB][EQUAL] << "\t"
		<< cmp_give_result_total_cnt[OB][NONIMPROVED] << "\t"
		<< f[BB] << "\t" << d[BB] << "\t" << h[BB] << "\t" << t[BB] << "\t"
		<< cmp_give_result_total_cnt[BB][IMPROVED] << "\t"
		<< cmp_give_result_total_cnt[BB][EQUAL] << "\t"
		<< cmp_give_result_total_cnt[BB][NONIMPROVED]
		<< endl;

}
void analyze_total_file(string fnr,string fnw)
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

		split(fields_vec, strline, is_any_of("\t"));
		/*for (vector<string>::iterator iter = fields_vec.begin();
		iter != fields_vec.end(); iter++)
		{
			cout << *iter << "\t";
			
		}*/
		if (fields_vec.front().find("Ni_") != string::npos)
		{
			split(ins_name_vec, fields_vec[0], is_any_of("_-"));
			sol_info_vec = new vector<SolutionInfo*>;
			insinfo = new InstanceInfo(fields_vec[0],
				stoi(fields_vec[1]), stoi(fields_vec[2]), 
				stoi(ins_name_vec[3]), stoi(ins_name_vec[4]), stoi(ins_name_vec[5]),
				boost::lexical_cast<double>(fields_vec[3]),
				boost::lexical_cast<double>(fields_vec[4]),stoi(fields_vec[5]),
				*sol_info_vec);
			insinfo_vec.push_back(insinfo);
		}
		else
		{
			sol_info_vec->push_back(new SolutionInfo(stoi(fields_vec[0]),//run_cnt
				boost::lexical_cast<double>(fields_vec[1]), stoi(fields_vec[2]),//obj,mm
				stoi(fields_vec[3]), stoi(fields_vec[4]),//iter, tm
				stoi(fields_vec[5]), stoi(fields_vec[6])//given_result_improve,result_improve
				));
			//cout << stoi(fields_vec[0]) << "\t"//run_cnt
			//	<< boost::lexical_cast<double>(fields_vec[1]) << "\t" << stoi(fields_vec[2]) << "\t"//obj<<"\t"mm
			//	<< stoi(fields_vec[3]) << "\t" << stoi(fields_vec[4]) << "\t"//iter<<"\t" tm
			//	<< stoi(fields_vec[5]) << "\t" << stoi(fields_vec[6]) << endl;
		}			
		//cout << endl;
	}
	ofs<< "(*iter)->filename  \t  (*iter)->n  \t  (*iter)->m  \t"
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
	iter!=insinfo_vec.end();iter++)
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
int main(int argc, char **argv)
{
	char *rgv[] = { "",	//0
		"_fn","tr_p13_itr2000_ptr50_rm2_ns0_t2_cp70_r1_r20",	//1,2
		"_if","instance\\BB_Problem_BestSolution\\",	//3,4	
		"_of","results\\",//5,6	
		"_p","13",		//7,8
		"_r","20",		//9,10
		"_itr","200",	//11,12
		"_ptr","30",	//13,14
		"_rm","1",	//15,16
		"_ns","0",	//17,18		
		"_r1","1",	//19,20
		"_r2","20",	//21,22
		"_ws","0",	//23,24
		"_t","2",	//25,26
		"_cp","90"	//27,28
	};
	argv = rgv;
	std::map<string, string> argv_map;
	for (int i = 1; i < sizeof(rgv) / sizeof(rgv[0]); i += 2)
		argv_map[string(argv[i])] = string(argv[i + 1]);
	string fnr, fnw,fnwt;
	fnr = argv_map.at("_of") + argv_map.at("_fn")+".txt";
	fnw = argv_map.at("_of") + argv_map.at("_fn") + "_analyze3.txt";
	fnwt = argv_map.at("_of") + "total_information13.txt";
	Analyze *an = new Analyze(fnr, fnw,fnwt);	
	//an->para_setting();
	an->total_info();
	//analyze_total_file(fnr,fnw);
	system("pause");
}
#endif