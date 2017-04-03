#!/bin/bash
# ****the setting and meaning of parameters****
#	_p 	20						[sol_num]; population size; default value: 13;
#	_non_p	10					[non_popu]; non population size; default value: 10; the real population size in hma is _p-_non_p;
#	_rm	1						[r_mode] in hma; for initial solution, jobs are sorted by non-decreasing order of [charact]; default value: 1; [1,9]; corresponding to enum R_Mode { EXAC = 0, MINP, MAXP, MIND, MAXD, MINPDD, MAXPDD, MINPD, MAXPD, RANDOM, SIZE };
#	_ns	0						[ns_mode] in hma; neighborhood search; default value: 0; corresponding to enum NS_Mode { SWAP, INSERT }; In hma,	compound neighborhood is used, so this parameter is out of use;
#	_itr	2000				[iteration] in hma; the maximum number of generation in hma; for OB set, hma stops when it finds the optimal 		solutions; for BB set, hma stops when it reaches _itr generations;
#	_ptr	500					[perturb_rate] in hma; _ptr*0.001; the proportion of jobs in the makespan machine needed to be swapped with jobs in other machine; default value: 500; [0,1000];
#	_cx	300						[alpha_cx_ptr]; rand()%1000<alpha_cx_ptr; the proportion of applying cx versus perturb; default value: 300; [0,1000];
#	_pu_beta	600				[pu_beta]; pu_beta=_pu_beta*0.001; pool update coefficient; default value: 600; [0,1000];
#	_ws	0						[whe_save_sol_seq]; whether save the solution sequence or not, used in save_solution; default value: 0; {0,1};
#	_r1	1						[run_cnt_from] in hma; default value: 1;
#	_r2 20						[run_cnt_to] in hma; default value: 20;
#	_px	null					prefix for output file name, default value: null;
#	_if 						input file path, default value: instance//BB_Problem_BestSolution// 
#	_of 						output file path, default value: results//
#	n_vec{ { 8, 11, 14 }, { 20, 35, 50 } };
#	m_vec{ { 2, 3, 4 }, { 4, 7, 10 } };
#	_vi1 0; _vi2 2,				[vi]; the domain of vi; [0,2];
#	_ni1 0; _ni2 3,				[ni]; the domain of ni; [0,3];
#	_vj1 0; _vj2 2,				[vj]; the domain of vj; [0,2];
#	_mj1 0; _mj2 3,				[mj]; the domain of mj; [0,3]; n_vec[vi][ni]; m_vec[vj][mj]; vi != vj;
#	_pi1 1; _pi2 2,				[pi]; the domain of pi; [1,2];
#	_di1 1; _di2 2,				[di]; the domain of di; [1,2];
#	_ins1 1, _ins2 25			[ins]; the domain of ins; [1,25];
#	_dc	1						[whe_dc]; whether use divide and conquer; default value: 1; {0,1};
#	_ls 3						[ls_method]; local search method, 0: basic compound (swap and insert); 1: add triple insert; 2: add triple swap; 3: add trip insert and trip swap; default value: 3; {0,1,2,3};
#	_pu 1						[pu_method]; pool updating method; default value: 1; 0: replace the worst solution; 1: based on quality and distance; default value: 1; {0,1};
#	_xo	1						[cx_method]; crossover method; 0: uniform crossover; 1: maximum common jobs derived; 2: objective oriented derived, resembles the crossover operator in graph coloring; default value: 1; {0,1,2};

#pms_s _px hma12_pms_s _if ..\\PMS1\\instance\\BB_Problem_BestSolution\\ _of ..\\PMS1\\results\\ 
#pms_i _px hma13_pms_i _if ..\\PMS1\\instance\\BB_Problem_BestSolution\\ _of ..\\PMS1\\results\\ 
#pms_non _px hma14_pms_non _if ..\\PMS1\\instance\\BB_Problem_BestSolution\\ _of ..\\PMS1\\results\\ 
#pms_ufx _px hma15_pms_ufx _if ..\\PMS1\\instance\\BB_Problem_BestSolution\\ _of ..\\PMS1\\results\\ 
#pms_mpx _px hma16_pms_mpx1 _if ..\\PMS1\\instance\\BB_Problem_BestSolution\\ _of ..\\PMS1\\results\\  _vi1 1 _ni1 2 _mj1 1
#pms_puqd_IS _px hma17_pms_puqd_is _if ..\\PMS1\\instance\\BB_Problem_BestSolution\\ _of ..\\PMS1\\results\\ 
#PMS_puqd_mpx _px hma17_pms_puqd_mpx _if ..\\PMS1\\instance\\BB_Problem_BestSolution\\ _of ..\\PMS1\\results\\ 
#PMS_nodc _px hma18_pms_nodc_mpx _if ..\\PMS1\\instance\\BB_Problem_BestSolution\\ _of ..\\PMS1\\results\\ 
#pms_puqd_i _px hma19_pms_puqd_i _if ..\\PMS1\\instance\\BB_Problem_BestSolution\\ _of ..\\PMS1\\results\\ 
#pms_puqd_s _px hma19_pms_puqd_s _if ..\\PMS1\\instance\\BB_Problem_BestSolution\\ _of ..\\PMS1\\results\\ 
#PMS1_para_setting _fn total_results14 _if ..\\PMS1\\instance\\BB_Problem_BestSolution\\ _of ..\\PMS1\\results\\ _p 13 _r 100 _itr 2000 _ptr 50 _rm 1 _ns 0 _r1 1 _r2 20 _ws 0 _t 2 _cp 70

#	g++ -O3 -std=c++0x(c++11) -o pms_dc.out pms_dc.cpp -lm
#	qsub -d $PWD -o $OEPATH($PWD/oe) -e $OEPATH pms_dc.sh
#	parameter settings _dc 1 _ls 3 _xo 1 _pu 1 in hma for pms_dc
	./pms_dc.out  _if instance//BB_Problem_BestSolution// _of results// _dc 1 _ls 3 _xo 1 _pu 1 _pu_beta 100
