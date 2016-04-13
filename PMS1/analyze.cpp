#if 1
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
	//	//<< c[sol_index_opt][0] << "\t"
	//	<< c[sol_index][0] << "\t"
	//	<< mm[sol_index] << "\t"
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
	int m, n, rand_seed;
	double opt_obj_given, opt_obj_real;
	vector<SolutionInfo*> &sol_info_vec;
	InstanceInfo(string _fn,int _m,int _n,double _oog,double _oor,int _rd,vector<SolutionInfo*> &_siv)
		:filename(_fn),m(_m),n(_n),opt_obj_given(_oog),opt_obj_real(_oor),rand_seed(_rd),sol_info_vec(_siv)
	{}
	
};
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
	vector<string> fields_vec;
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
		for (vector<string>::iterator iter = fields_vec.begin();
		iter != fields_vec.end(); iter++)
		{
		//	cout << *iter << "\t";
			
		}
		if (fields_vec.front().find("Ni_") != string::npos)
		{
			sol_info_vec = new vector<SolutionInfo*>;
			insinfo = new InstanceInfo(fields_vec[0],
				stoi(fields_vec[1]), stoi(fields_vec[2]), 
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
		int cmp_give_result[3] = { 0,0,0 }, cmp_real_result[3] = { 0,0,0 }, 
			cmp_give_result_cnt[3] = {0,0,0}, cmp_real_result_cnt[3] = { 0,0,0 },
			cmp_give_result_t,com_real_result_t;
		for (vector<SolutionInfo*>::iterator iter_sol = (*iter)->sol_info_vec.begin();
		iter_sol != (*iter)->sol_info_vec.end(); iter_sol++)
		{
			//ofs << run_cnt << "\t"
			//	//<< c[sol_index_opt][0] << "\t"
			//	<< c[sol_index][0] << "\t"
			//	<< mm[sol_index] << "\t"
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
			if (abs((*iter_sol)->obj - (*iter)->opt_obj_real) <= MIN_EQUAL)
			{
				cmp_real_result[1] += 1;	// equal
				com_real_result_t = 1;
			}
			else if ((*iter)->opt_obj_real - (*iter_sol)->obj > MIN_EQUAL)
			{
				cmp_real_result[0] += 1;	// improved
				com_real_result_t = 2;
			}
			else
			{
				cmp_real_result[2] += 1;	// not improved
				com_real_result_t = 0;
			}
			if (abs((*iter_sol)->obj - (*iter)->opt_obj_given) <= MIN_EQUAL)
			{
				cmp_give_result[1] += 1;	// equal
				cmp_give_result_t = 1;
			}
			else if ((*iter)->opt_obj_given - (*iter_sol)->obj > MIN_EQUAL)
			{
				cmp_give_result[0] += 1;	// improved
				cmp_give_result_t = 2;
			}
			else
			{
				cmp_give_result[2] += 1;	// not improved
				cmp_give_result_t = 0;
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
			if (com_real_result_t != (*iter_sol)->result_improve ||
				cmp_give_result_t != (*iter_sol)->given_result_improve ||
				cmp_give_result_cnt[0] != cmp_give_result[0] ||
				cmp_give_result_cnt[1] != cmp_give_result[1] ||
				cmp_give_result_cnt[2] != cmp_give_result[2] ||
				cmp_real_result_cnt[0] != cmp_real_result[0]||
				cmp_real_result_cnt[1] != cmp_real_result[1] || 
				cmp_real_result_cnt[2] != cmp_real_result[2] )
			{
				cout << "ERROR, cmp_result wrong." << endl;
				cout << (*iter)->filename << "\t"
					<< (*iter)->opt_obj_given << "\t"
					<< (*iter)->opt_obj_real << "\t"
					<< (*iter_sol)->run_cnt << "\t"
					<< (*iter_sol)->obj << endl;
				system("pause");
			}
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
			<< cmp_give_result[0] << "\t" << cmp_give_result[1] << "\t" << cmp_give_result[2] << "\t"
			<< cmp_real_result_cnt[0] << "\t" << cmp_real_result_cnt[1] << "\t" << cmp_real_result_cnt[2] << "\t"
			<< cmp_real_result[0] << "\t" << cmp_real_result[1] << "\t" << cmp_real_result[2] << endl;
			
	}
	ifs.close();
}
int main(int argc, char **argv)
{
	char *rgv[] = { "",	//0
		"_fn","total_results9_p13_itr2000_ptr50_rm1_ns0_r1_r20",	//1,2
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
	string fnr, fnw;
	fnr = argv_map.at("_of") + argv_map.at("_fn")+".txt";
	fnw = argv_map.at("_of") + argv_map.at("_fn") + "_analyze.txt";
	analyze_total_file(fnr,fnw);
	system("pause");
}
#endif