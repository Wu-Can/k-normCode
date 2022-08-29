Welcome to k-normCode! 
Authors: Can Wu, Ying Cui, Donghui Li, Defeng Sun. 

A MATLAB software package for solving the convex and nonconvex risk-based 
linear regression at scale

Set up
(1) Unpack the software
(2) Run Matlab in the directory knormCode
(3) In the Matlab command window, type: 
    >> startup 
(4) generate the UCI dataset as follows:

   - Download the original datasets tested in the paper from UCI data repository (or download the UCI dataset in LIBSVM
     format at https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression.html);
   - Put them into the folder [\genUCIdatafun\UCIdataorg];
   - Run genUCIdata.m (For demonstration purpose, "mpg_scale.txt" is left in the folder
     [\genUCIdatafun\UCIdataorg]).


By now, k-normCode is ready for you to use.
Notice: Run all the files are provided for only demonstration purpose.
  
    
Contents:
===============================================================================================================================
Test1_Risk-averse_model: Comparison of the risk-sensitive and risk-neutral regression models

(a) Test_convex_CVaR_model_random: Comparison of the convex CVaR-based linear regression model, Average model and Maximum model; 
(b) Test_truncated_CVaR_model_random: Comparison of the truncated CVaR-based linear regression model and Average model.
===============================================================================================================================

===============================================================================================================================
Test2_fixed_lambda: Solvers for the convex CVaR-based sparse linear regression with a fixed lambda: N-ALM, ADMM, Gurobi, S-IRPN

(a1) Test_NALM_UCI: N-ALM for the convex CVaR-based models with UCI data; 
(a2) Test_NALM_random: N-ALM for the convex CVaR-based models with random data; 
(b1) Test_ADMM_UCI: ADMM for the convex CVaR-based models with UCI data; 
(b2) Test_ADMM_random: ADMM for the convex CVaR-based models with random data; 
(c1) Test_Barrier_Gurobi_UCI: Barrier method in Gurobi for the convex CVaR-based models with UCI data; 
(c2) Test_Barrier_Gurobi_random: Barrier method in Gurobi for the convex CVaR-based models with random data; 
(d1) Test_SIRPN_UCI: S-IRPN for the convex CVaR-based models with UCI data; 
(d2) Test_SIRPN_random: S-IRPN for the convex CVaR-based models with random data.
===============================================================================================================================

===============================================================================================================================
Test3_Solution_path: Solvers to generate the whole solution path of the convex CVaR-based models: AS+N-ALM, Warm+N-ALM, N-ALM

(a1) Test_AS_NALM_path_UCI: AS+N-ALM for the computation of a solution path with UCI data; 
(a2) Test_AS_NALM_path_random: AS+N-ALM for the computation of a solution path with random data; 
(b1) Test_warm_NALM_path_UCI: Warm+N-ALM for the computation of a solution path with UCI data; 
(b2) Test_warm_NALM_path_random: Warm+N-ALM for the computation of a solution path with random data; 
(c1) Test_NALM_path_UCI: N-ALM for the computation of a solution path with UCI data; 
(c2) Test_NALM_path_random: N-ALM for the computation of a solution path with random data;
(d1) Test_path_figure_random: Generate figures for solution paths on random data;
(d2) Test_path_figure_UCI: Generate figures for solution paths on UCI data;
===============================================================================================================================

===============================================================================================================================
Test4_truncated_CVaR_MM: Solvers for the nonconvex truncated CVaR-based sparse linear regression: MM+N-PPA, MM+Gurobi

(a1) Test_MM_NPPA_UCI: MM+N-PPA for the nonconvex truncated CVaR-based model with UCI data; 
(a2) Test_MM_NPPA_random: MM+N-PPA for the nonconvex truncated CVaR-based model with random data; 
(b1) Test_MM_Gurobi_UCI: MM+Gurobi for the nonconvex truncated CVaR-based model with UCI data; 
(b2) Test_MM_Gurobi_random: MM+Gurobi for the nonconvex truncated CVaR-based model with random data; 
===============================================================================================================================

===============================================================================================================================
genUCIdatafun: Convert the UCI data from LIBSVM format to MAT format;
mainfun: All main functions for the above solvers;
mexfun: two funtions in MEX format;
solver: all the functions called by the above solvers;
UCIdata: all availiable data for the UCI data and the high-accurate objective values.
===============================================================================================================================





    
