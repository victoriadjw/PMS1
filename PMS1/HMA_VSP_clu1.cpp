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
typedef float weight_type;
class Graph
{
public:
	static const int MAX_VERTEX_NUM=1000;
	static const int MAX_DELTA_OBJ=9999999;
	long long iteration,global_iteration;
	struct timeval start_tm,end_tm;
	//clock_t start_tm,end_tm;
	int whether_save_sol_seq,cx_A,cx_B,cx_C,alpha_cx_qd,qd_tabu_iter_a,qd_tabu_iter_b;
	float beta_qd,t1a,t1b,t2a,t2b;
	struct Edge
	{
		int vertex;
		Edge *next;
	};
	struct separator
	{
		int next;
		int pre;
	};
	//typedef separator* type_list;
	separator **sep;// the separator list
	int *first;//the first element in the separator list
	int ***constraint_matrix,***obj_gain_matrix,***dis_gain_matrix;
	//used in pool update
	int *tabu_pool_update,**dis;
	double *g;
	int result_best[55];	

	Edge **head;
	weight_type *vw;
	weight_type *sol_obj_phase1,*sol_obj_phase2;
	int *tabu[3];
	int **sol_partition;
	int b,vertex_num,edge_num,sol_num,ins_num;
	Graph(char *,char *,int,int);
	void add_edge(int,int);
	void init_solution(int,int);
	void display_graph();
	void display_solution(int);
	void check_solution(int,int,int);
	void move_gain1(int,int,int,weight_type&,weight_type&,weight_type&);
	void move_gain2(int,int,int,weight_type&,weight_type&,weight_type&);
	void perform_move1(int,int,int,int,weight_type,weight_type);
	void perform_move2(int,int,int,weight_type,weight_type);
	void replace_solution(int,int);
	void check_single_move(int,int,int);
	void check_improve(int,int);
	void tabu_search1(int,int,int,int,char*,int,int);
	void tabu_search2(int,int,int,int,char*,int,int);
	void save_solution(int,char *,int,int);
	void gpx(int,int,int);
	void gpx1(int,int,int);
	void couple_search(int,char *,int,int,int);
	void pool_update(int,int*,int,int);
	void calculate_f1_history_best(int,int);
	void add_to_separator_list(int,separator *list, int &first,int num);
	void delete_from_separator(int,separator *list, int &first,int num);
	void check_gain_matrix1(int,int);
	void check_gain_matrix2(int);
private:
	Edge *edge_set;
	//int cmp(const void *,const void *);
	int top;
};
Graph::Graph(char *filename, char *file_best_result, int _sol_num,int _ins_num):sol_num(_sol_num),ins_num(_ins_num)
{
	FILE *fp=fopen(filename, "r");
	FILE *fprr=fopen(file_best_result, "r");
	int ch;
	char buff[2048];
	fscanf(fp,"%d%d%d",&b,&vertex_num,&edge_num);
	head=new Edge *[vertex_num];
	edge_set=new Edge[edge_num*2];
	for(int i=0;i<edge_num*2;i++)
		edge_set[i].next=NULL;
	for(int i=0;i<vertex_num;i++)
		head[i]=NULL;
	for(int i=0;i<3;i++)
	{
		tabu[i]=new int[vertex_num];
	}
	sol_partition=new int *[sol_num];
	for(int i=0;i<sol_num;i++)
		sol_partition[i]=new int[vertex_num+3];//最后一个存储obj
	vw=new weight_type[vertex_num];
	sol_obj_phase1=new weight_type[sol_num];
	sol_obj_phase2=new weight_type[sol_num];
	top=0;
	for(int i=0;i<vertex_num;i++)
		fscanf(fp,"%*s%d%f",&ch,&vw[i]);
	int v1,v2;
	float weight;
	for(int i=0;i<edge_num;i++)
	{
		fscanf(fp,"%*s%d%d%f",&v1,&v2,&weight);
		//printf("%d,%d,%f\n",v1,v2,weight);
		add_edge(v1,v2);
		add_edge(v2,v1);
	}
	sep=new separator *[sol_num];
	first=new int[sol_num];
	for(int i=0;i<sol_num;i++)
		sep[i]=new separator[vertex_num+1];
	obj_gain_matrix=new int**[sol_num];
	dis_gain_matrix=new int**[sol_num];
	constraint_matrix=new int**[sol_num];

	tabu_pool_update=new int[sol_num];
	dis=new int *[sol_num];
	g=new double[sol_num];
	for(int i=0;i<sol_num;i++)
	{
		obj_gain_matrix[i]=new int *[3];
		dis_gain_matrix[i]=new int *[3];
		constraint_matrix[i]=new int *[3];
		dis[i]=new int[sol_num];
		for(int j=0;j<3;j++)
		{
			obj_gain_matrix[i][j]=new int[vertex_num];
			dis_gain_matrix[i][j]=new int[vertex_num];
			constraint_matrix[i][j]=new int[vertex_num];
			memset(obj_gain_matrix[i][j],0,vertex_num*sizeof(int));
			memset(dis_gain_matrix[i][j],0,vertex_num*sizeof(int));
			memset(constraint_matrix[i][j],0,vertex_num*sizeof(int));
		}
	}
	for(int i=1;i<55;i++)
	{
		fscanf(fprr,"%d",&result_best[i]);
	}
	////for g1.rud////
	/*fscanf(fp,"%d%d",&vertex_num,&edge_num);
	head=new Edge *[vertex_num];
	vw=new weight_type [vertex_num];
	edge_set=new Edge[edge_num*2];
	for(int i=0;i<edge_num*2;i++)
		edge_set[i].next=NULL;
	for(int i=0;i<vertex_num;i++)
		head[i]=NULL;
	sol_partition=new int *[sol_num];
	sol_obj_phase2=new weight_type[sol_num];
	for(int i=0;i<sol_num;i++)
		sol_partition[i]=new int[vertex_num+4];//最后一个存储obj
	for(int i=0;i<vertex_num;i++)
		vw[i]=1;
	top=0;
	int v1,v2;
	float weight;
	for(int i=0;i<edge_num;i++)
	{
		fscanf(fp,"%d%d%f",&v1,&v2,&weight);
		//printf("%d,%d,%f\n",v1,v2,weight);
		add_edge(v1,v2);
		add_edge(v2,v1);
	}
	b=0.525*vertex_num;*/
}
void Graph::add_edge(int v1,int v2)
{
	edge_set[top].vertex=v2;
	edge_set[top].next=head[v1];
	head[v1]=edge_set+top++;
}
void Graph::display_graph()
{
	for(int i=0;i<vertex_num;i++)
	{
		cout<<i<<":";
		for(Edge *e=head[i];e;e=e->next)
			cout<<e->vertex<<"\t";
		cout<<endl;
	}
	cout<<"vertex_num:"<<vertex_num
		<<"edge_num:"<<edge_num
		<<"b:"<<b<<endl;
}
void Graph::check_solution(int sol_index,int sol_index_history_best,int phase)
{
	int cnt_a,cnt_b,cnt_c;
	weight_type real_obj=0;
	cnt_a=cnt_b=cnt_c=0;
	int dis=0;
	for(int i=0;i<vertex_num;i++)
	{
		if(sol_partition[sol_index][i]==1)
			cnt_a+=1;
		else if(sol_partition[sol_index][i]==2)
			cnt_b+=1;
		else
		{
			cnt_c+=1;
			real_obj+=vw[i];
			if(phase==1&&sol_partition[sol_index_history_best][i]==0)
				dis+=1;
		}
		if(sol_partition[sol_index][i]!=1)
			continue;
		for(Edge *e=head[i];e;e=e->next)
		{
			if(sol_partition[sol_index][e->vertex]!=2)
				continue;
			cout<<"invaild solution: "<<sol_partition[sol_index][i]<<","
			<<sol_partition[sol_index][e->vertex]<<endl;
			break;
		}
	}
	if(phase==1)
	{
		if(cnt_a==sol_partition[sol_index][vertex_num+1]&&cnt_a>=1&&cnt_a<=b&&
			cnt_b==sol_partition[sol_index][vertex_num+2]&&cnt_b>=1&&cnt_b<=b&&
			cnt_c==vertex_num-sol_partition[sol_index][vertex_num+1]-sol_partition[sol_index][vertex_num+2]&&
			real_obj+dis*beta_qd==sol_obj_phase1[sol_index])
			;//cout<<"chech:A,B,C size,const,obj are right,phase1."<<endl;
		else 
			cout<<"check error, size wrong,phase1."
			<<cnt_a<<"\t"
			<<cnt_b<<"\t"
			<<cnt_c<<"\n";
	}else
	{
		if(cnt_a==sol_partition[sol_index][vertex_num+1]&&cnt_a>=1&&cnt_a<=b&&
			cnt_b==sol_partition[sol_index][vertex_num+2]&&cnt_b>=1&&cnt_b<=b&&
			cnt_c==vertex_num-sol_partition[sol_index][vertex_num+1]-sol_partition[sol_index][vertex_num+2]&&
			real_obj==sol_obj_phase2[sol_index])
			;//cout<<"chech:A,B,C size,const,obj are right,phase2."<<endl;
		else 
			cout<<"check error, size wrong,phase2."
			<<cnt_a<<"\t"
			<<cnt_b<<"\t"
			<<cnt_c<<"\n";
	}
}
void Graph::check_improve(int sol_index,int phase)
{
	for(int i=0;i<vertex_num;i++)
	{
		if(sol_partition[sol_index][i]!=0)
			continue;
		int sol_index_history_best1,sol_index_history_best2;
		weight_type dis_gain1,obj_gain1,total_gain1,dis_gain2,obj_gain2,total_gain2;
		move_gain1(sol_index,i,1,obj_gain1,dis_gain1,total_gain2);
		move_gain1(sol_index,i,2,obj_gain2,dis_gain2,total_gain2);
		if(total_gain1<0||total_gain2<0)
		{
			cout<<"still can be improved"<<endl;
		}
	}
}
void Graph::check_gain_matrix1(int sol_index,int sol_index_history_best)
{
	for(int iter=first[sol_index];iter!=-1;iter=sep[sol_index][iter].next)
	{
		int delta_obj1,delta_obj2,constraint1,constraint2,dis1,dis2;
		delta_obj1=delta_obj2=-vw[iter];
		dis1=dis2=0;
		if(sol_partition[sol_index_history_best][iter]==0)
			dis1=dis2=-1;
		constraint1=constraint2=0;
		for(Edge *e=head[iter];e;e=e->next)
		{
			if(sol_partition[sol_index][e->vertex]==1)
			{
				delta_obj2+=vw[e->vertex];
				constraint2-=1;
				if(sol_partition[sol_index_history_best][e->vertex]==0)
					dis2+=1;
			}
			else if(sol_partition[sol_index][e->vertex]==2)
			{
				delta_obj1+=vw[e->vertex];
				constraint1-=1;
				if(sol_partition[sol_index_history_best][e->vertex]==0)
					dis1+=1;
			}
		}
		if(delta_obj1!=obj_gain_matrix[sol_index][1][iter]||delta_obj2!=obj_gain_matrix[sol_index][2][iter])
		{
			system("pause");
		}
		if(constraint1!=constraint_matrix[sol_index][1][iter]||constraint2!=constraint_matrix[sol_index][2][iter])
		{
			system("pause");
		}
		if(dis1!=dis_gain_matrix[sol_index][1][iter]||dis2!=dis_gain_matrix[sol_index][2][iter])
		{
			system("pause");
		}
	}
}
void Graph::check_gain_matrix2(int sol_index)
{
	for(int iter=first[sol_index];iter!=-1;iter=sep[sol_index][iter].next)
	{
		int delta_obj1,delta_obj2,constraint1,constraint2,dis1,dis2;
		delta_obj1=delta_obj2=-vw[iter];
		constraint1=constraint2=0;
		for(Edge *e=head[iter];e;e=e->next)
		{
			if(sol_partition[sol_index][e->vertex]==1)
			{
				delta_obj2+=vw[e->vertex];
				constraint2-=1;
			}
			else if(sol_partition[sol_index][e->vertex]==2)
			{
				delta_obj1+=vw[e->vertex];
				constraint1-=1;
			}
		}
		if(delta_obj1!=obj_gain_matrix[sol_index][1][iter]||delta_obj2!=obj_gain_matrix[sol_index][2][iter])
		{
			system("pause");
		}
		if(constraint1!=constraint_matrix[sol_index][1][iter]||constraint2!=constraint_matrix[sol_index][2][iter])
		{
			system("pause");
		}
	}
}
void Graph::init_solution(int sol_index, int ini)
{
	for(int i=1;i<=2;i++)
		memset(tabu[i],0,vertex_num*sizeof(int));
	iteration=0;

	if(ini)
	{
		sol_partition[sol_index][vertex_num+1]=0;
		sol_partition[sol_index][vertex_num+2]=0;
		sol_obj_phase1[sol_index]=0;
		sol_obj_phase2[sol_index]=0;
		
		first[sol_index]=-1;
		for(int i=0;i<vertex_num;i++)
		{
			sep[sol_index][i].next=sep[sol_index][i].pre=-1;
			int flag=0;
			if(rand()%2)
			{
				if(sol_partition[sol_index][vertex_num+1]<b)
				{
					sol_partition[sol_index][i]=1;
					sol_partition[sol_index][vertex_num+1]+=1;
					flag=1;
				}
			}else
			{
				if(sol_partition[sol_index][vertex_num+2]<b)
				{
					sol_partition[sol_index][i]=2;
					sol_partition[sol_index][vertex_num+2]+=1;
					flag=1;
				}
			}
		}
		for(int i=0;i<vertex_num;i++)
		{
			if(sol_partition[sol_index][i]!=1)
				continue;
			for(Edge *e=head[i];e;e=e->next)
			{
				if(sol_partition[sol_index][e->vertex]!=2)
					continue;
				if(rand()%2)
				{
					sol_partition[sol_index][i]=0;
					sol_partition[sol_index][vertex_num+1]-=1;
					sol_obj_phase2[sol_index]+=vw[i];
					add_to_separator_list(sol_index,sep[sol_index],first[sol_index],i);
					break;
				}
				else
				{
					sol_partition[sol_index][e->vertex]=0;
					sol_partition[sol_index][vertex_num+2]-=1;	
					sol_obj_phase2[sol_index]+=vw[e->vertex];
					add_to_separator_list(sol_index,sep[sol_index],first[sol_index],e->vertex);
				}
			}
		}
		if(sol_partition[sol_index][vertex_num+1]<=0||sol_partition[sol_index][vertex_num+2]<=0)
		{
			cout<<"initialize error"<<endl;
			system("pause");
		}
		//////initialize the gain and constraint matrix		
		for(int iter=first[sol_index];iter!=-1;iter=sep[sol_index][iter].next)
		{
			obj_gain_matrix[sol_index][1][iter]=obj_gain_matrix[sol_index][2][iter]=-vw[iter];
			constraint_matrix[sol_index][1][iter]=constraint_matrix[sol_index][2][iter]=0;
			for(Edge *e=head[iter];e;e=e->next)
			{				
				if(sol_partition[sol_index][e->vertex]==0)
					continue;
				obj_gain_matrix[sol_index][3-sol_partition[sol_index][e->vertex]][iter]+=vw[e->vertex];
				constraint_matrix[sol_index][3-sol_partition[sol_index][e->vertex]][iter]-=1;
			}
		}
		
	}
	/*cout<<"init:"<<sol_partition[sol_index][vertex_num+1]<<"\t"
	<<sol_partition[sol_index][vertex_num+2]<<"\t"
	<<sol_partition[sol_index][vertex_num+3]<<"\t"
	<<sol_obj_phase2[sol_index]<<"\n";*/
	//display_solution(sol_index);
}
void Graph::display_solution(int sol_index)
{
	int part_cnt=0,sep_cnt=0;
	for(int i=0;i<vertex_num;i++)
	{
		if(sol_partition[sol_index][i]==0)
		{
			cout<<i<<" ";
			part_cnt+=1;
		}
	}
	cout<<"*"<<part_cnt<<endl;
	int iter=first[sol_index];
	while(iter!=-1)
	{
		cout<<iter<<" ";
		iter=sep[sol_index][iter].next;
		sep_cnt+=1;
	}
	cout<<"*"<<sep_cnt<<endl;
}
void Graph::add_to_separator_list(int sol_index,separator *list, int &first,int num)
{
	if(first==-1)
	{
		first=num;
		list[num].pre=-1;
		list[num].next=-1;
	}else
	{
		list[first].pre=num;
		list[num].next=first;
		list[num].pre=-1;
		first=num;
	}
}
void Graph::delete_from_separator(int sol_index,separator *list, int &first, int num)
{
	if(list[num].pre==-1)
		first=list[num].next;
	else
		list[list[num].pre].next=list[num].next;
	if(list[num].next==-1)
		list[list[num].pre].next=-1;
	else
		list[list[num].next].pre=list[num].pre;
	list[num].next=list[num].pre=-1;
}
void Graph::move_gain1(int sol_index, int sel_v, int sel_s,weight_type &obj_gain,weight_type &dis_gain, weight_type &total_gain)
{
	if(sol_partition[sol_index][vertex_num+sel_s]==b)// violate |A|<=b
	{
		//cout<<"violate |A|<=b"<<endl;
		total_gain= MAX_DELTA_OBJ+1;
		return;
	}
	if(sol_partition[sol_index][vertex_num+3-sel_s]+constraint_matrix[sol_index][sel_s][sel_v]==0)// violate 1<=|A|
	{
		//cout<<"violate 1<=|A|"<<endl;
		total_gain= MAX_DELTA_OBJ+2;
		return;
	}

	/*dis_gain=0;
	if(sol_partition[sol_index_history_best][sel_v]==0)
		dis_gain=-1;
	for(Edge *e=head[sel_v];e;e=e->next)
	{
		if(sol_partition[sol_index][e->vertex]==3-sel_s)
		{
			if(sol_partition[sol_index_history_best][e->vertex]==0)
				dis_gain+=1;
		}
	}*/
	/*if(dis_gain!=dis_gain_matrix[sol_index][sel_s][sel_v])
	{
		system("pause");
	}
	if(obj_gain!=obj_gain_matrix[sol_index][sel_s][sel_v])
	{
		system("pause");
	}
	if(op_remove_cnt+constraint_matrix[sol_index][sel_s][sel_v]!=0)
	{
		system("pause");
	}*/
	obj_gain=obj_gain_matrix[sol_index][sel_s][sel_v];
	dis_gain=dis_gain_matrix[sol_index][sel_s][sel_v];
	total_gain=obj_gain+dis_gain*beta_qd;
	
}
void Graph::move_gain2(int sol_index, int sel_v, int sel_s,weight_type &obj_gain,weight_type &dis_gain, weight_type &total_gain)
{
	if(sol_partition[sol_index][vertex_num+sel_s]==b)// violate |A|<=b
	{
		//cout<<"violate |A|<=b"<<endl;
		total_gain= MAX_DELTA_OBJ+1;
		return;
	}
	if(sol_partition[sol_index][vertex_num+3-sel_s]+constraint_matrix[sol_index][sel_s][sel_v]==0)// violate 1<=|A|
	{
		//cout<<"violate 1<=|A|"<<endl;
		total_gain= MAX_DELTA_OBJ+2;
		return;
	}
	/*if(obj_gain!=obj_gain_matrix[sol_index][sel_s][sel_v])
	{
		system("pause");
	}
	if(op_remove_cnt+constraint_matrix[sol_index][sel_s][sel_v]!=0)
	{
		system("pause");
	}*/
	total_gain=obj_gain=obj_gain_matrix[sol_index][sel_s][sel_v];
}
void Graph::perform_move1(int sol_index,int sol_index_history_best,int sel_v,int sel_s,weight_type delta_obj,weight_type obj_gain)
{
	sol_partition[sol_index][sel_v]=sel_s;
	sol_partition[sol_index][vertex_num+sel_s]+=1;
	delete_from_separator(sol_index,sep[sol_index],first[sol_index],sel_v);
	obj_gain_matrix[sol_index][1][sel_v]=obj_gain_matrix[sol_index][2][sel_v]=0;
	dis_gain_matrix[sol_index][1][sel_v]=dis_gain_matrix[sol_index][2][sel_v]=0;
	constraint_matrix[sol_index][1][sel_v]=constraint_matrix[sol_index][2][sel_v]=0;
	sol_obj_phase1[sol_index]+=delta_obj;
	sol_obj_phase2[sol_index]+=obj_gain;
	for(Edge *e=head[sel_v];e;e=e->next)
	{
		if(sol_partition[sol_index][e->vertex]==3-sel_s)
		{
			//int prev_set=sol_partition[sol_index][e->vertex];
			sol_partition[sol_index][e->vertex]=0;
			sol_partition[sol_index][vertex_num+3-sel_s]-=1;
			add_to_separator_list(sol_index,sep[sol_index],first[sol_index],e->vertex);
				
			obj_gain_matrix[sol_index][1][e->vertex]=obj_gain_matrix[sol_index][2][e->vertex]=-vw[e->vertex];
			dis_gain_matrix[sol_index][1][e->vertex]=dis_gain_matrix[sol_index][2][e->vertex]=0;
			if(sol_partition[sol_index_history_best][e->vertex]==0)
				dis_gain_matrix[sol_index][1][e->vertex]=dis_gain_matrix[sol_index][2][e->vertex]=-1;
			constraint_matrix[sol_index][1][e->vertex]=constraint_matrix[sol_index][2][e->vertex]=0;
			for(Edge *g=head[e->vertex];g;g=g->next)
			{				
				if(sol_partition[sol_index][g->vertex]==0)
				{
					obj_gain_matrix[sol_index][sel_s][g->vertex]-=vw[e->vertex];
					if(sol_partition[sol_index_history_best][e->vertex]==0)
						dis_gain_matrix[sol_index][sel_s][g->vertex]-=1;
					constraint_matrix[sol_index][sel_s][g->vertex]+=1;
				}else
				{
					obj_gain_matrix[sol_index][3-sol_partition[sol_index][g->vertex]][e->vertex]+=vw[g->vertex];
					if(sol_partition[sol_index_history_best][g->vertex]==0)
						dis_gain_matrix[sol_index][3-sol_partition[sol_index][g->vertex]][e->vertex]+=1;
					constraint_matrix[sol_index][3-sol_partition[sol_index][g->vertex]][e->vertex]-=1;
				}
			}
			
		}else if(sol_partition[sol_index][e->vertex]==0)
		{
			obj_gain_matrix[sol_index][3-sel_s][e->vertex]+=vw[sel_v];
			if(sol_partition[sol_index_history_best][sel_v]==0)
				dis_gain_matrix[sol_index][3-sel_s][e->vertex]+=1;
			constraint_matrix[sol_index][3-sel_s][e->vertex]-=1;
		}
	}
}
void Graph::perform_move2(int sol_index,int sel_v,int sel_s,weight_type delta_obj,weight_type obj_gain)
{
	sol_partition[sol_index][sel_v]=sel_s;
	sol_partition[sol_index][vertex_num+sel_s]+=1;
	delete_from_separator(sol_index,sep[sol_index],first[sol_index],sel_v);
	obj_gain_matrix[sol_index][1][sel_v]=obj_gain_matrix[sol_index][2][sel_v]=0;
	constraint_matrix[sol_index][1][sel_v]=constraint_matrix[sol_index][2][sel_v]=0;
	sol_obj_phase2[sol_index]+=obj_gain;
	for(Edge *e=head[sel_v];e;e=e->next)
	{
		if(sol_partition[sol_index][e->vertex]==3-sel_s)
		{
			//int prev_set=sol_partition[sol_index][e->vertex];
			sol_partition[sol_index][e->vertex]=0;
			sol_partition[sol_index][vertex_num+3-sel_s]-=1;
			add_to_separator_list(sol_index,sep[sol_index],first[sol_index],e->vertex);
				
			obj_gain_matrix[sol_index][1][e->vertex]=obj_gain_matrix[sol_index][2][e->vertex]=-vw[e->vertex];
			constraint_matrix[sol_index][1][e->vertex]=constraint_matrix[sol_index][2][e->vertex]=0;
			for(Edge *g=head[e->vertex];g;g=g->next)
			{				
				if(sol_partition[sol_index][g->vertex]==0)
				{
					obj_gain_matrix[sol_index][sel_s][g->vertex]-=vw[e->vertex];
					constraint_matrix[sol_index][sel_s][g->vertex]+=1;
				}else
				{
					obj_gain_matrix[sol_index][3-sol_partition[sol_index][g->vertex]][e->vertex]+=vw[g->vertex];
					constraint_matrix[sol_index][3-sol_partition[sol_index][g->vertex]][e->vertex]-=1;
				}
			}
			
		}else if(sol_partition[sol_index][e->vertex]==0)
		{
			obj_gain_matrix[sol_index][3-sel_s][e->vertex]+=vw[sel_v];
			constraint_matrix[sol_index][3-sel_s][e->vertex]-=1;
		}
	}
}
void Graph::replace_solution(int dest,int src)
{
	memcpy(sol_partition[dest],sol_partition[src],(vertex_num+3)*sizeof(int));
	memcpy(sep[dest],sep[src],(vertex_num+1)*sizeof(separator));
	first[dest]=first[src];
	sol_obj_phase2[dest]=sol_obj_phase2[src];
	for(int i=1;i<=2;i++)
	{
		memcpy(obj_gain_matrix[dest][i],obj_gain_matrix[src][i],vertex_num*sizeof(int));
		memcpy(constraint_matrix[dest][i],constraint_matrix[src][i],vertex_num*sizeof(int));
	}
}
void Graph::calculate_f1_history_best(int sol_index,int sol_index_history_best)
{
	sol_obj_phase1[sol_index]=sol_obj_phase2[sol_index];
	int same_cnt=0;
	for(int i=first[sol_index];i!=-1;i=sep[sol_index][i].next)
	{
		dis_gain_matrix[sol_index][1][i]=dis_gain_matrix[sol_index][2][i]=0;
		if(sol_partition[sol_index_history_best][i]==0)
		{
			same_cnt+=1;
			dis_gain_matrix[sol_index][1][i]=dis_gain_matrix[sol_index][2][i]=-1;
		}
		for(Edge *e=head[i];e;e=e->next)
		{
			if(sol_partition[sol_index][e->vertex]==0)
				continue;
			if(sol_partition[sol_index_history_best][e->vertex]==0)
				dis_gain_matrix[sol_index][3-sol_partition[sol_index][e->vertex]][i]+=1;
		}
	}
	sol_obj_phase1[sol_index]+=beta_qd*same_cnt;
}
void Graph::tabu_search1(int sol_index,int itr_num,int is_display,int sol_index_history_best,char *fnw, int gen,int run_num)
{
	int *origin_index=new int[vertex_num];
	int *random_sel_index=new int[vertex_num];
	calculate_f1_history_best(sol_index,sol_index_history_best);
	//check_gain_matrix1(sol_index,sol_index_history_best);
	replace_solution(0,sol_index);
	sol_obj_phase1[0]=sol_obj_phase1[sol_index];
	//check_solution(sol_index,sol_index_history_best,phase);
	int equ_cnt_phase1,equ_cnt_phase2;
	equ_cnt_phase1=equ_cnt_phase2=1;
	for(int iter=0;iter<itr_num;iter++)
	{
		weight_type max_gain=2*vertex_num;
		int equ_cnt;
		int sel_v,sel_s;
		weight_type obj_gain,dis_gain;
		//int n_s3=vertex_num-sol_partition[sol_index][vertex_num+1]-sol_partition[sol_index][vertex_num+2];
		//int rand_tt=(int)(ceil(0.3*n_s3+rand()%(int)ceil(0.2*n_s3)));		
		int rand_tt=(int)(t1a*sol_obj_phase2[sol_index])+rand()%(int)(t1b*sol_obj_phase2[sol_index]);		
		for(int sep_v=first[sol_index];sep_v!=-1;sep_v=sep[sol_index][sep_v].next)
		{
			int rd_s=rand()%2+1;
			weight_type obj_gain1,obj_gain2,dis_gain1,dis_gain2,total_gain1,total_gain2;
			move_gain1(sol_index,sep_v,rd_s,obj_gain1,dis_gain1,total_gain1);
			move_gain1(sol_index,sep_v,3-rd_s,obj_gain2,dis_gain2,total_gain2);
			//cout<<i<<"\t"<<rd_s<<"\t"<<obj_gain1<<"\t"<<dis_gain1<<"\t"<<gain1<<endl;
			//cout<<i<<"\t"<<3-rd_s<<"\t"<<obj_gain2<<"\t"<<dis_gain2<<"\t"<<gain2<<endl;
			if(total_gain1<=max_gain&&(total_gain1+sol_obj_phase1[sol_index]<sol_obj_phase1[0]||iteration>=tabu[rd_s][sep_v]+rand_tt))
			{
				if(total_gain1<max_gain)
				{
					equ_cnt=1;
					max_gain=total_gain1;
					sel_v=sep_v;
					sel_s=rd_s;
					obj_gain=obj_gain1;
					/*if(obj_gain==-1)
						break;*/
				}else if(total_gain1==max_gain)
				{
					equ_cnt+=1;
					if(rand()%equ_cnt==0)
					{
						sel_v=sep_v;
						sel_s=rd_s;
						obj_gain=obj_gain1;
					/*if(obj_gain==-1)
						break;*/
					}
				}
			}
			if(total_gain2<=max_gain&&(total_gain2+sol_obj_phase1[sol_index]<sol_obj_phase1[0]||iteration>=tabu[3-rd_s][sep_v]+rand_tt))
			{
				if(total_gain2<max_gain)
				{
					equ_cnt=1;
					max_gain=total_gain2;
					sel_v=sep_v;
					sel_s=3-rd_s;
					obj_gain=obj_gain2;
					/*if(obj_gain==-1)
						break;*/
				}else if(total_gain2==max_gain)
				{
					equ_cnt+=1;
					if(rand()&equ_cnt==0)
					{
						sel_v=sep_v;
						sel_s=3-rd_s;
						obj_gain=obj_gain2;
					/*if(obj_gain==-1)
						break;*/
					}
				}
			}
		}
		global_iteration++;
		if(max_gain==2*vertex_num)
		{
			iteration++;
			continue;
		}
		perform_move1(sol_index,sol_index_history_best,sel_v,sel_s,max_gain,obj_gain);
		//check_gain_matrix1(sol_index,sol_index_history_best);
		/*if(obj_gain==0)
			system("pause");*/
		tabu[sel_s][sel_v]=iteration;
		iteration++;
		if(is_display&&iter%1000==0)
		{
			cout<<iter<<"\t"<<sol_obj_phase1[sol_index]<<"\t"<<sol_obj_phase1[0]<<"\t"
				<<vertex_num-sol_partition[sol_index][vertex_num+1]-sol_partition[sol_index][vertex_num+2]<<"\tb "
				<<sol_obj_phase2[sol_index_history_best]<<"...";
			cout<<sel_v<<"\t"<<sel_s<<"\t"<<obj_gain<<"\t"<<dis_gain<<endl;
		}

		if(sol_obj_phase1[sol_index]<sol_obj_phase1[0])
		{
			replace_solution(0,sol_index);
			sol_obj_phase1[0]=sol_obj_phase1[sol_index];
			equ_cnt_phase1=1;
		}else if(sol_obj_phase1[sol_index]==sol_obj_phase1[0])
		{
			equ_cnt_phase1+=1;
			if(rand()%equ_cnt_phase1==0)
			{
				replace_solution(0,sol_index);
				sol_obj_phase1[0]=sol_obj_phase1[sol_index];
			}
		}
		if(sol_obj_phase2[sol_index]<sol_obj_phase2[sol_index_history_best])
		{
			replace_solution(sol_index_history_best,sol_index);
			memset(dis_gain_matrix[sol_index][1],-1,vertex_num*sizeof(int));
			memset(dis_gain_matrix[sol_index][2],-1,vertex_num*sizeof(int));
			sol_obj_phase1[sol_index]=(1+beta_qd)*sol_obj_phase2[sol_index];
			replace_solution(0,sol_index);
			sol_obj_phase1[0]=sol_obj_phase1[sol_index];
			equ_cnt_phase2=1;
			save_solution(sol_index_history_best,fnw,gen,run_num);
		}else if(sol_obj_phase2[sol_index]==sol_obj_phase2[sol_index_history_best])
		{
			equ_cnt_phase2+=1;
			if(rand()%equ_cnt_phase2==0)
			{
				replace_solution(sol_index_history_best,sol_index);
				memset(dis_gain_matrix[sol_index][1],-1,vertex_num*sizeof(int));
				memset(dis_gain_matrix[sol_index][2],-1,vertex_num*sizeof(int));
				sol_obj_phase1[sol_index]=(1+beta_qd)*sol_obj_phase2[sol_index];
				replace_solution(0,sol_index);
				sol_obj_phase1[0]=sol_obj_phase1[sol_index];
			}
		}
		//check_solution(sol_index,sol_index_history_best,phase);
	}
	replace_solution(sol_index,0);
	//display_solution(sol_index);
}
void Graph::tabu_search2(int sol_index,int itr_num,int is_display,int sol_index_history_best,char *fnw, int gen,int run_num)
{
	int *origin_index=new int[vertex_num];
	int *random_sel_index=new int[vertex_num];
	replace_solution(0,sol_index);
	//check_gain_matrix(sol_index,sol_index_history_best);
	//check_solution(sol_index,sol_index_history_best,phase);
	int equ_cnt_phase1,equ_cnt_phase2;
	equ_cnt_phase1=equ_cnt_phase2=1;
	for(int iter=0;iter<itr_num;iter++)
	{
		weight_type max_gain=2*vertex_num;
		int equ_cnt;
		int sel_v,sel_s;
		weight_type obj_gain,dis_gain;
		int rand_tt=(int)(t2a*sol_obj_phase2[sol_index])+rand()%(int)(t2b*sol_obj_phase2[sol_index]);		
		for(int sep_v=first[sol_index];sep_v!=-1;sep_v=sep[sol_index][sep_v].next)
		{
			int rd_s=rand()%2+1;
			weight_type obj_gain1,obj_gain2,dis_gain1,dis_gain2,total_gain1,total_gain2;
			move_gain2(sol_index,sep_v,rd_s,obj_gain1,dis_gain1,total_gain1);
			move_gain2(sol_index,sep_v,3-rd_s,obj_gain2,dis_gain2,total_gain2);
			if(total_gain1<=max_gain&&(total_gain1+sol_obj_phase2[sol_index]<sol_obj_phase2[0]||iteration>=tabu[rd_s][sep_v]+rand_tt))
			{
				if(total_gain1<max_gain)
				{
					equ_cnt=1;
					max_gain=total_gain1;
					sel_v=sep_v;
					sel_s=rd_s;
					obj_gain=obj_gain1;
					/*if(obj_gain==-1)
						break;*/
				}else if(total_gain1==max_gain)
				{
					equ_cnt+=1;
					if(rand()%equ_cnt==0)
					{
						sel_v=sep_v;
						sel_s=rd_s;
						obj_gain=obj_gain1;
					/*if(obj_gain==-1)
						break;*/
					}
				}
			}
			if(total_gain2<=max_gain&&(total_gain2+sol_obj_phase2[sol_index]<sol_obj_phase2[0]||iteration>=tabu[3-rd_s][sep_v]+rand_tt))
			{
				if(total_gain2<max_gain)
				{
					equ_cnt=1;
					max_gain=total_gain2;
					sel_v=sep_v;
					sel_s=3-rd_s;
					obj_gain=obj_gain2;
					/*if(obj_gain==-1)
						break;*/
				}else if(total_gain2==max_gain)
				{
					equ_cnt+=1;
					if(rand()&equ_cnt==0)
					{
						sel_v=sep_v;
						sel_s=3-rd_s;
						obj_gain=obj_gain2;
					/*if(obj_gain==-1)
						break;*/
					}
				}
			}
		}
		global_iteration++;
		if(max_gain==2*vertex_num)
		{
			iteration++;
			continue;
		}
		perform_move2(sol_index,sel_v,sel_s,max_gain,obj_gain);
		//check_gain_matrix(sol_index,sol_index_history_best);
		/*if(obj_gain==0)
			system("pause");*/
		tabu[sel_s][sel_v]=iteration;
		iteration++;
		if(is_display&&iter%1000==0)
			cout<<iter<<"\t"<<sol_obj_phase2[sol_index]<<"\t"<<sol_obj_phase2[0]<<endl;
		if(sol_obj_phase2[sol_index]<sol_obj_phase2[0])
		{
			replace_solution(0,sol_index);
			equ_cnt_phase1=1;
			if(sol_obj_phase2[sol_index]<sol_obj_phase2[sol_index_history_best])
				save_solution(sol_index,fnw,gen,run_num);
		}else if(sol_obj_phase2[sol_index]==sol_obj_phase2[0])
		{
			equ_cnt_phase1+=1;
			if(rand()%equ_cnt_phase1==0)
				replace_solution(0,sol_index);
		}
		//check_solution(sol_index,sol_index_history_best,phase);
	}
	replace_solution(sol_index,0);
	//display_solution(sol_index);
}
void Graph::save_solution(int sol_index,char *filename,int gen,int run_cnt)
{
	//end_tm=clock();
	gettimeofday(&end_tm, NULL);
	ofstream fp(filename,ios::app|ios::out);
	fp<<run_cnt<<"\t"
		<<sol_partition[sol_index][vertex_num+1]<<"\t"
		<<sol_partition[sol_index][vertex_num+2]<<"\t"
		<<vertex_num-sol_partition[sol_index][vertex_num+1]-sol_partition[sol_index][vertex_num+2]<<"\t"
		<<sol_obj_phase2[sol_index]<<"\t"
		<<global_iteration<<"\t"
		<<gen<<"\t"
		<<end_tm.tv_sec-start_tm.tv_sec
		//<<(end_tm-start_tm)/CLOCKS_PER_SEC
		<<endl;
	/*cout<<run_cnt<<"\t"
		<<sol_partition[sol_index][vertex_num+1]<<"\t"
		<<sol_partition[sol_index][vertex_num+2]<<"\t"
		<<vertex_num-sol_partition[sol_index][vertex_num+1]-sol_partition[sol_index][vertex_num+2]<<"\t"
		<<sol_obj_phase2[sol_index]<<"\t"
		<<global_iteration<<"\t"
		<<gen<<"\t"
		<<(end_tm-start_tm)/CLOCKS_PER_SEC
		<<endl;*/
	if(whether_save_sol_seq)
	{
		for(int i=0;i<vertex_num;i++)
			fp<<sol_partition[sol_index][i]<<" ";
		fp<<endl;
	}
}
void Graph::gpx(int dest,int p1,int p2)
{
	int same_cnt=0,same_p1=0,same_p2=0;
	replace_solution(dest,p1);
	for(int i=0;i<vertex_num;i++)
	{
		if(sol_partition[p1][i]==0||sol_partition[p2][i]==0)
		{
			if((sol_partition[p1][i]==0&&sol_partition[p2][i]!=0&&rand()%1000<cx_A)||
				(sol_partition[p1][i]==0&&sol_partition[p2][i]==0&&rand()%1000<cx_C))
			{
				int rd_s=rand()%2+1;
				weight_type obj_gain1,dis_gain1,total_gain1,obj_gain2,dis_gain2,total_gain2;
				move_gain2(dest,i,rd_s,obj_gain1,dis_gain1,total_gain1);
				move_gain2(dest,i,3-rd_s,obj_gain2,dis_gain2,total_gain2);
				if(total_gain1<=MAX_DELTA_OBJ&&total_gain1<total_gain2)
					perform_move2(dest,i,rd_s,total_gain1,obj_gain1);
				else if(total_gain2<=MAX_DELTA_OBJ&&total_gain2<total_gain1)
					perform_move2(dest,i,3-rd_s,total_gain2,obj_gain2);
			}else if(sol_partition[p1][i]!=0&&sol_partition[p2][i]==0&&rand()%1000<cx_B)
			{
				if(sol_partition[dest][vertex_num+sol_partition[dest][i]]>1)
				{
					if(sol_partition[dest][i]!=0)
					{
						int prev_set=sol_partition[dest][i];
						sol_partition[dest][vertex_num+sol_partition[dest][i]]-=1;
						sol_obj_phase2[dest]+=vw[i];
						sol_partition[dest][i]=0;
						add_to_separator_list(dest,sep[dest],first[dest],i);

						obj_gain_matrix[dest][1][i]=obj_gain_matrix[dest][2][i]=-vw[i];
						constraint_matrix[dest][1][i]=constraint_matrix[dest][2][i]=0;
						for(Edge *e=head[i];e;e=e->next)
						{
							if(sol_partition[dest][e->vertex]==0)
							{
								obj_gain_matrix[dest][3-prev_set][e->vertex]-=vw[i];
								constraint_matrix[dest][3-prev_set][e->vertex]+=1;
							}else
							{
								obj_gain_matrix[dest][3-sol_partition[dest][e->vertex]][i]+=vw[e->vertex];
								constraint_matrix[dest][3-sol_partition[dest][e->vertex]][i]-=1;
							}
						}
					}
				}
			}
		}
	}
}
void Graph::pool_update(int offspring,int *tabu_pool_update,int gen,int phase)
{
	double beta_pool_update=0.6;
	for(int i=1;i<sol_num-1;i++)
	{
		for(int j=i+1;j<sol_num-1;j++)
		{
			dis[i][j]=0;
			for(int k=0;k<vertex_num;k++)
			{
				if(sol_partition[i][k]==0&&sol_partition[j][k]!=0)
				{
					dis[i][j]+=1;
				}else if(sol_partition[i][k]!=0&&sol_partition[j][k]==0)
				{
					dis[i][j]+=1;
				}
			}
		}
	}
	int max_dis,min_dis;
	weight_type max_obj,min_obj;
	for(int i=1;i<sol_num-1;i++)//from 1 to sol_num-2;0 is used in TS,
	{// sol_num-1 is the best solution, sol_num-2 is the offspring solution.
		dis[sol_num-1][i]=0;//distance between one solution and a population.
		for(int j=1;j<sol_num-1;j++)
		{
			if(i<j)
			{
				dis[sol_num-1][i]+=dis[i][j];
			}else if(i>j)
			{
				dis[sol_num-1][i]+=dis[j][i];
			}
		}
		if(i==1)
		{
			max_dis=min_dis=dis[sol_num-1][i];
			max_obj=min_obj=sol_obj_phase2[i];
		}else
		{
			if(max_dis<dis[sol_num-1][i])
				max_dis=dis[sol_num-1][i];
			if(min_dis>dis[sol_num-1][i])
				min_dis=dis[sol_num-1][i];
			if(max_obj<sol_obj_phase2[i])
				max_obj=sol_obj_phase2[i];
			if(min_obj>sol_obj_phase2[i])
				min_obj=sol_obj_phase2[i];
		}
	}
	int equ_cnt,gw_index;
	for(int i=1;i<=sol_num-2;i++)
	{
		g[i]=beta_pool_update*(max_obj-sol_obj_phase2[i])/(double)(max_obj-min_obj+1);
		g[i]+=(1-beta_pool_update)*(dis[sol_num-1][i]-min_dis)/(double)(max_dis-min_dis+1);
	}
	double gw=2;
	int rand_tt=sol_num/3+rand()%(sol_num/2);
	for(int i=1;i<sol_num-2;i++)
	{
		if(gw>g[i]&&(gen>=tabu_pool_update[i]+rand_tt))
		{
			gw=g[i];
			gw_index=i;
			equ_cnt=1;
		}else if(gw==g[i])
		{
			equ_cnt+=1;
			if(rand()%equ_cnt==0)
				gw_index=i;
		}
	}
	if(gw==2)
	{
		gw_index=rand()%(sol_num-3)+1;
		//cout<<" at ";
	}
	if(g[offspring]>=g[gw_index]||rand()%1000<300)
	{
		replace_solution(gw_index,offspring);
		tabu_pool_update[gw_index]=gen;
		//if(g[offspring]>=g[gw_index])
		//	cout<<gw_index<<"$ ";
		//else
		//	cout<<gw_index<<" ";
	}
}
void Graph::couple_search(int tabu_iter,char *fnw,int run_num, int stop_non_imp_cnt,int stop_gen_condition)
{
	int best=sol_num-1,offspring=sol_num-2;
	int tabu_perturb_iter=tabu_iter;
	memset(tabu_pool_update,0,sol_num*sizeof(int));
	gettimeofday(&start_tm, NULL);
	//start_tm=clock();
	global_iteration=0;
	init_solution(best,1);
	int equ_best_cnt=1;
	for(int i=1;i<=sol_num-2;i++)
	{
		init_solution(i,1);
		tabu_search2(i,tabu_iter,0,best,fnw,0,run_num);
		if(sol_obj_phase2[i]<sol_obj_phase2[best])
		{
			replace_solution(best,i);
			equ_best_cnt=1;
			if(sol_obj_phase2[best]==result_best[ins_num])
				return;
		}
		else if(sol_obj_phase2[i]==sol_obj_phase2[best])
		{
			equ_best_cnt+=1;
			if(rand()%equ_best_cnt==0)
				replace_solution(best,i);
		}
		//display_solution(i);
		//cout<<"\n"<<sep_size[i]<<"\t"<<sol_obj_phase2[i]<<endl;
	}
	int non_con_imp_cnt=0;
	int generation=stop_gen_condition==1?stop_non_imp_cnt:INT_MAX;
	for(int gen=0;gen<generation;gen++)
	{
		int p1=rand()%(sol_num-3)+1;
		int p2=rand()%(sol_num-3)+1;
		while(p1==p2)
			p2=rand()%(sol_num-3)+1;
		
		if(rand()%1000<alpha_cx_qd)
		{
			gpx(offspring,p1,p2);
			
		}else
		{
			int iter_tabu_phase1=qd_tabu_iter_a+rand()%qd_tabu_iter_b;
			replace_solution(offspring,p1);
			init_solution(offspring,0);
			tabu_search1(offspring,iter_tabu_phase1,0,best,fnw,gen+1,run_num);
			//check_solution(offspring,best,1);
		}
		//check_solution(offspring,best,1);
		init_solution(offspring,0);
		tabu_search2(offspring,tabu_iter,0,best,fnw,gen+1,run_num);
		//check_solution(offspring,best,2);

		//display_solution(offspring);
		if(sol_obj_phase2[offspring]<sol_obj_phase2[best])
		{
			replace_solution(best,offspring);			
			if(sol_obj_phase2[best]==result_best[ins_num])
				return;
			equ_best_cnt=1;
			non_con_imp_cnt=0;
		}
		else if(sol_obj_phase2[offspring]==sol_obj_phase2[best])
		{
			equ_best_cnt+=1;
			if(rand()%equ_best_cnt==0)
				replace_solution(best,offspring);
			non_con_imp_cnt+=1;
		}else
			non_con_imp_cnt+=1;
		if(stop_gen_condition==0&&non_con_imp_cnt>stop_non_imp_cnt)
			return;
		//cout<<gen<<":";
		//for(int i=0;i<sol_num;i++)
			//cout<<sol_obj_phase2[i]<<" ";
		pool_update(offspring,tabu_pool_update,gen,2);
		//cout<<endl;		
	}
}
int main(int argc, char **argv)
{
	srand(time(NULL));
	int i;
	for(i=0;argv[36][i]!='\r'&&i<10;i++)
	;argv[36][i]='\0';
	/*const char *rgv[]={"",//0
		"_p","13",//1,2
		"_ws","0",//3,4
		"_xa","500",//5,6
		"_xb","300",//7,8
		"_xc","100",//9,10
		"_alpha","700",//11,12
		"_beta","250",//13,14
		"_qa","100",//15,16
		"_qb","5000",//17,18

		"_fn","G14",//19,20
		"_r","20",//21,22
		"_s","3000",//23,24
		"_t","10000",//25,26
		"_sc","1",//27,28
		"_t1a","100",//29,30
		"_t1b","200",//31,32
		"_t2a","200",//33,34
		"_t2b","200",//35,36	
		"_rb","vsp_best54"//37,38	
	};
	const char **argv=rgv;*/
	char fnr[1000]="instance//una_instance//";
	char fnrr[1000]="";
	char fnw[1000]="results//";
	strcpy(fnrr,fnr);
	strcat(fnr,argv[20]);//filename
	strcat(fnr,".txt");
	strcat(fnrr,argv[38]);//filename of best result
	strcat(fnrr,".txt");

	strcat(fnw,argv[20]);//filename
	strcat(fnw,argv[21]);//"r"
	strcat(fnw,argv[22]);//run_cnt
	strcat(fnw,argv[23]);//"s"
	strcat(fnw,argv[24]);//stop_non_imp_cnt or max gen.
	strcat(fnw,argv[27]);//"sc"
	strcat(fnw,argv[28]);//stop_gen_condition
	strcat(fnw,argv[25]);//"t"
	strcat(fnw,argv[26]);//tabu_iter
	strcat(fnw,argv[1]);//"p"
	strcat(fnw,argv[2]);//population size

	strcat(fnw,argv[11]);//"alpha"
	strcat(fnw,argv[12]);//alpha_cx_qd
	strcat(fnw,argv[13]);//"beta"
	strcat(fnw,argv[14]);//beta_qd
	strcat(fnw,argv[3]);//"ws"
	strcat(fnw,argv[4]);//whether_sav_sol_seq
	strcat(fnw,argv[5]);//"xa"
	strcat(fnw,argv[6]);//cx_A
	strcat(fnw,argv[7]);//"xb"
	strcat(fnw,argv[8]);//cx_B
	strcat(fnw,argv[9]);//"xc"
	strcat(fnw,argv[10]);//cx_C
	strcat(fnw,argv[15]);//"qa"
	strcat(fnw,argv[16]);//qd_tabu_iter_a
	strcat(fnw,argv[17]);//"qb"
	strcat(fnw,argv[18]);//qd_tabu_iter_b
	strcat(fnw,argv[29]);//"t1a"
	strcat(fnw,argv[30]);//t1a
	strcat(fnw,argv[31]);//"t1b"
	strcat(fnw,argv[32]);//t1b
	strcat(fnw,argv[33]);//"t2a"
	strcat(fnw,argv[34]);//t2a
	strcat(fnw,argv[35]);//"t2b"
	strcat(fnw,argv[36]);//t2b
	strcat(fnw,".txt");

	char ins_n[10]="";
	for(i=0;argv[20][i]!='\0';i++)	
	{
		if(i!=0)
		ins_n[i-1]=argv[20][i];
	};
	for(i=0;i<atoi(argv[22]);i++)
	{		
		Graph *g=new Graph(fnr,fnrr,atoi(argv[2]),atoi(ins_n));	
		g->whether_save_sol_seq=atoi(argv[4]);
		g->cx_A=atoi(argv[6]);
		g->cx_B=atoi(argv[8]);
		g->cx_C=atoi(argv[10]);
		g->alpha_cx_qd=atoi(argv[12]);
		g->beta_qd=atoi(argv[14])*0.01;
		g->qd_tabu_iter_a=atoi(argv[16]);
		g->qd_tabu_iter_b=atoi(argv[18]);
		g->t1a=atoi(argv[30])*0.001;
		g->t1b=atoi(argv[32])*0.001;
		g->t2a=atoi(argv[34])*0.001; 
		g->t2b=atoi(argv[36])*0.001;
		g->couple_search(atoi(argv[26]),fnw,i+1,atoi(argv[24]),atoi(argv[28]));
		delete g;
	}
	return 0;
}
#endif 