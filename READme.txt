Welcome to k-normCode! 
Authors: Can Wu, Ying Cui, Donghui Li, Defeng Sun. 

A MATLAB software package for solving the convex and nonconvex risk-based 
linear regression at scale

===============================================================================================================================
SETTINGS:

(1) Unpack the software.
(2) Run Matlab in the directory knormCode.
(3) In the Matlab command window, type: 
    >> startup 
(4) generate the UCI dataset as follows:
   - Download the original datasets tested in the paper from UCI data repository (or download the UCI dataset in LIBSVM
     format at https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression.html);
   - Put them into the folder [\genUCIdatafun\UCIdataorg];
   - Run genUCIdata.m (For demonstration purpose, "mpg_scale.txt" and "housing_scale.txt" are left in the folder
     [\genUCIdatafun\UCIdataorg]).
(5) Download and install the version 9.5.2 of the Gurobi optimizer from https://www.gurobi.com/downloads/gurobi-software/.
(6) Download and install the version 3.2 of the PSG program from http://www.aorda.com/index.php/downloading/.

By now, k-normCode is ready for you to use.
Notice: Run all the files are provided for only demonstration purpose.
===============================================================================================================================

===============================================================================================================================
CONTENTS:

(1) Test1_fixed_lambda: Solvers for the convex CVaR-based sparse linear regression with a fixed lambda: N-ALM, ADMM, Gurobi, S-IRPN
PSG solvers
    (a1) Test_NALM_UCI: N-ALM for the convex CVaR-based models with UCI data; 
    (a2) Test_NALM_random: N-ALM for the convex CVaR-based models with random data; 
    (b1) Test_ADMM_UCI: ADMM for the convex CVaR-based models with UCI data; 
    (b2) Test_ADMM_random: ADMM for the convex CVaR-based models with random data; 
    (c1) Test_Barrier_Gurobi_UCI: Barrier method in Gurobi for the convex CVaR-based models with UCI data; 
    (c2) Test_Barrier_Gurobi_random: Barrier method in Gurobi for the convex CVaR-based models with random data; 
    (d1) Test_SIRPN_UCI: S-IRPN for the convex CVaR-based models with UCI data; 
    (d2) Test_SIRPN_random: S-IRPN for the convex CVaR-based models with random data;
    (e1) Test_PSG_UCI: PSG solvers for the convex CVaR-based models with UCI data;
    (e2) Test_PSG_random: PSG solvers for the convex CVaR-based models with random data.

(2) Test2_Solution_path: Solvers to generate the whole solution path of the convex CVaR-based models: AS+N-ALM, Warm+N-ALM, N-ALM
    (a1) Test_AS_NALM_path_UCI: AS+N-ALM for the computation of a solution path with UCI data; 
    (a2) Test_AS_NALM_path_random: AS+N-ALM for the computation of a solution path with random data; 
    (b1) Test_warm_NALM_path_UCI: Warm+N-ALM for the computation of a solution path with UCI data; 
    (b2) Test_warm_NALM_path_random: Warm+N-ALM for the computation of a solution path with random data; 
    (c1) Test_NALM_path_UCI: N-ALM for the computation of a solution path with UCI data; 
    (c2) Test_NALM_path_random: N-ALM for the computation of a solution path with random data;
    (d1) Test_path_figure_random: Generate figures for solution paths on random data;
    (d2) Test_path_figure_UCI: Generate figures for solution paths on UCI data.

(3) Test3_truncated_CVaR_MM: Solvers for the nonconvex truncated CVaR-based sparse linear regression: MM+N-PPA, MM+Gurobi
    (a1) Test_MM_NPPA_UCI: MM+N-PPA for the nonconvex truncated CVaR-based model with UCI data; 
    (a2) Test_MM_NPPA_random: MM+N-PPA for the nonconvex truncated CVaR-based model with random data; 
    (b1) Test_MM_Gurobi_UCI: MM+Gurobi for the nonconvex truncated CVaR-based model with UCI data; 
    (b2) Test_MM_Gurobi_random: MM+Gurobi for the nonconvex truncated CVaR-based model with random data. 

(4) Test4_Risk-averse_model: Comparison of the risk-sensitive and risk-neutral regression models
    (a) Test_convex_CVaR_model_random: Comparison of the convex CVaR-based linear regression model, Average model and Maximum model; 
    (b) Test_truncated_CVaR_model_random: Comparison of the truncated CVaR-based linear regression model and Average model.

(5) genUCIdatafun: Convert the UCI data from LIBSVM format to MAT format;

(6) mainfun: All main functions for the above solvers;

(7) mexfun: two funtions in MEX format;

(8) solver: all the functions called by the above solvers;

(9) UCIdata: all availiable data for the UCI data and the high-accurate objective values.
===============================================================================================================================

===============================================================================================================================
EXAMPLES:

Example 1: If you want to test the performance of the N-ALM for solving the convex CVaR-based models under a fixed lambda with 
UCI data, please perform the following steps:

(1) Open the file "Test_NALM_UCI.m" in the folder "Test1_fixed_lambda";
(2) Input the values of prob and flag_tol, and then run the file "Test_NALM_UCI.m" in its directory;
(3) The final result can be found in the workspace.

Example 2: If you want to generate a figure of the cumulative time comparison by the AS+N-ALM, the Warm+N-ALM and the N-ALM 
with UCI data, please perform the following steps:

(1) In th folder "Test2_Solution_path", open and run the files "Test_AS_NALM_path_UCI.m" "Test_warm_NALM_path_UCI.m" and 
    "Test_NALM_path_UCI.m" in sequence;
(2) Open and run the file "Test_path_figure_UCI.m" in the folder "Result_solution_path_figure_UCI";
(3) The figure "Time_comparison_9_4_2.fig" will be generated.
===============================================================================================================================
