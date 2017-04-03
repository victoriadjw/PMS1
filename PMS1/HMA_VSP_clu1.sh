#!/bin/bash
#const char *rgv[]={"",
#		"_p","13",//1,2			population size; DEFAULT VALUE:13;
#		"_ws","0",//3,4			whether_sav_sol_seq; DEFAULT VALUE:0;{0,1},0:not save,1:save;
#		"_xa","500",//5,6		cx_A; DEFAULT VALUE:500; [0,1000]; the proportion of p1=0,p2!=0;
#		"_xb","300",//7,8		cx_B; DEFAULT VALUE:300; [0,1000]; the proportion of p1!=0,p2=0;
#		"_xc","100",//9,10		cx_C; DEFAULT VALUE:100; [0,1000]; the proportion of p1=0,p2=0;
#		"_alpha","700",//11,12	alpha_cx_qd; DEFAULT VALUE:700; [0,1000]; the proportion of applying cx versus qd;
#		"_beta","700",//13,14	beta_qd; DEFAULT VALUE:1; [0,+INF]; the proportion of quality and distance, F=f+beta*d; NOTE: beta=beta*0.001;
#		"_qa","100",//15,16		qd_tabu_iter_a; DEFAULT VALUE:100; the iteration used in qd_TS = qd_tabu_iter_a+rand()%qd_tabu_iter_b; 
#		"_qb","5000",//17,18	qd_tabu_iter_b; DEFAULT VALUE:5000; the iteration used in qd_TS = qd_tabu_iter_a+rand()%qd_tabu_iter_b; 
#		"_fn","G1",//19,20		filename; without extensions;
#		"_r","20",//21,22		run_cnt; DEFAULT VALUE:20;
#		"_s","100",//23,24		stop_non_imp_cnt (maximum generations); consecutive non-improvement occurs in generations, stop condition;
#		"_t","1099",//25,26		tabu_iter; DEFAULT VALUE:10000; the iteration used in basic TS;
#		"_sc","1",//27,28		stop_gen_condition; DEFAULT VALUE:1; {0,1}; 1: program stop when it reaches the maximum generations(stop_non_imp_cnt); 0: stop not depend on generations, stop when consecutive non-improvement generation number reaches stop_non_imp_cnt;
#		"_t1a","100",//29,30	in tabu_search1, tt=t1a*obj2+rand()%(t1b*obj2); quality and distance guided tabu search.
#		"_t1b","200",//31,32	in tabu_search1, tt=t1a*obj2+rand()%(t1b*obj2); quality and distance guided tabu search.
#		"_t2a","200",//33,34	in tabu_search1, tt=t2a*obj2+rand()%(t2b*obj2); quality only guided tabu search.
#		"_t2b","200",//35,36	in tabu_search1, tt=t2a*obj2+rand()%(t2b*obj2); quality only guided tabu search.
#		"_rb","vsp_best54"//37,38	filename of the best results of the 54 large instances. The program stop when it reaches the history best found solution.
#		"_if","\\instance\\una_instance\\"//39,40	filename of the best results of the 54 large instances. The program stop when it reaches the history best found solution.
#		"_of","F:\\MyProjects\\VSP\\results\\"//41,42	filename of the best results of the 54 large instances. The program stop when it reaches the history best found solution.

./HMA_VSP_clu1.out _p 13 _ws 0 _xa 500 _xb 300 _xc 100 _alpha 700  _beta 250 _qa 100 _qb 5000 _fn G14 _r 5 _s 300 _t 100 _sc 1 _t1a 100 _t1b 200 _t2a 200 _t2b 200 _rb vsp_best54
