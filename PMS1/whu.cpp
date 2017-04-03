#if 1
#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<fstream>
#include<sys/time.h>
#include<limits.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
using namespace std;
		
int main(int argc, char **argv)
{
	int cur_t=time(NULL);
	srand(cur_t);
	cout<<"this is the program to test the server of WHU."<<endl;
	cout<<cur_t<<endl;
	return 0;
}
#endif 