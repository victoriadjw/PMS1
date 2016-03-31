#if 0
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

class TestClass
{
private:
	const proc_type **p;	// processing time
public:
	TestClass(proc_type**);
};
TestClass::TestClass(proc_type **_p) :p(_p)
{
	cout << **p << endl;
}
void main()
{
	proc_type a = 123;
	proc_type **p = new proc_type*[a];
	for (int i = 0; i < a; i++)
		p[i] = new proc_type[10];
	TestClass tc(p);
	system("pause");
}
#endif