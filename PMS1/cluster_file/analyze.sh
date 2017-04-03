#!/bin/bash
# ****the setting and meaning of parameters****
#	_fn							[ifs], fnr; the filename to be analyzed;
#	_if 						input file path, default value: instance//BB_Problem_BestSolution//;
#	_of 						output file path, default value： results//;
#	_efn 						[ofs], fnw, each output file name, append to the source file, trunc mode; default value: _analyze0;
#	_tfn 						[ofst], fnwt; total information of all the files, append mode; default value: total_information0;
#	_optfn  					[ifs_best], best found solution values file name; default value: opt_best_results;
#	_whe_save_each				[whe_save_each_file_result]; whether save the results of each file;
#	_total_info_per_para		[which function to run], invoke total_information function or per parameter function; 1: tota_information; 0: per parameter; default value 1;

#	g++ -O3 -std=c++0x(c++11) -o pms_dc.out pms_dc.cpp -lm
#	qsub -d $PWD -o $OEPATH($PWD/oe) -e $OEPATH pms_dc.sh
#	change the \\ to // for _if and _of, re-compile the analyze.cpp
#	be aware that _ arguments are pairwise appeared

#./analyze.out  _tfn per_para _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr50_rm1_alpha30_beta60_r1_r20.txt	_total_info_per_para 0

./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha0_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha25_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha50_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha75_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha100_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha125_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha150_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha175_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha200_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha225_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha250_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha275_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha300_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha325_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha350_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha375_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha400_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha425_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha450_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha475_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha500_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha525_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha550_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha575_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha600_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha625_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha650_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha675_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha700_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha725_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha750_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha775_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha800_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha825_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha850_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha875_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha900_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha925_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha950_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha975_beta600_r1_r20	
./analyze.out  _tfn total_alpha1 _fn	hma_dc1_ls3_ox1_pu1_p20_np10_itr2000_ptr500_rm1_alpha1000_beta600_r1_r20	

	
		
