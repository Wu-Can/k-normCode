# k-normCode
A MATLAB software package for solving the convex and nonconvex risk-based  linear regression at scale

Authors: Can Wu, Ying Cui, Donghui Li, Defeng Sun. 

Set up

(1) Unpack the software

(2) Run Matlab in the directory k-normCode

(3) In the Matlab command window, type: 
    >> startup 
    
(4) Generate the UCI dataset as follows:

-  Download the original datasets tested in the paper from UCI data repository (or download the UCI dataset in LIBSVM
    format at https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression.html);   
- Put them into the folder [\genUCIdatafun\UCIdataorg];
- Run genUCIdata.m (For demonstration purpose, "mpg_scale.txt" is left in the folder
    [\genUCIdatafun\UCIdataorg]).
    
(5) Download the latest version 9.5.2 of the Gurobi optimizer from https://www.gurobi.com/downloads/gurobi-software/.

By now, k-normCode is ready for you to use. Notice: Run all the files are provided for only demonstration purpose.

Contents:

(1) Test1_Risk-averse_model: Comparison of the risk-sensitive and risk-neutral regression models

- Test_convex_CVaR_model_random: Comparison of the convex CVaR-based linear regression model, Average model and Maximum model; 
- Test_truncated_CVaR_model_random: Comparison of the truncated CVaR-based linear regression model and Average model.

(2) Test2_fixed_lambda: Solvers for the convex CVaR-based sparse linear regression with a fixed lambda: N-ALM, ADMM, Gurobi, S-IRPN

- Test_NALM_UCI: N-ALM for the convex CVaR-based models with UCI data; 
- Test_NALM_random: N-ALM for the convex CVaR-based models with random data; 
- Test_ADMM_UCI: ADMM for the convex CVaR-based models with UCI data; 
- Test_ADMM_random: ADMM for the convex CVaR-based models with random data; 
- Test_Barrier_Gurobi_UCI: Barrier method in Gurobi for the convex CVaR-based models with UCI data; 
- Test_Barrier_Gurobi_random: Barrier method in Gurobi for the convex CVaR-based models with random data; 
- Test_SIRPN_UCI: S-IRPN for the convex CVaR-based models with UCI data; 
- Test_SIRPN_random: S-IRPN for the convex CVaR-based models with random data.

(3) Test3_Solution_path: Solvers to generate the whole solution path of the convex CVaR-based models: AS+N-ALM, Warm+N-ALM, N-ALM

- Test_AS_NALM_path_UCI: AS+N-ALM for the computation of a solution path with UCI data; 
- Test_AS_NALM_path_random: AS+N-ALM for the computation of a solution path with random data; 
- Test_warm_NALM_path_UCI: Warm+N-ALM for the computation of a solution path with UCI data; 
- Test_warm_NALM_path_random: Warm+N-ALM for the computation of a solution path with random data; 
- Test_NALM_path_UCI: N-ALM for the computation of a solution path with UCI data; 
- Test_NALM_path_random: N-ALM for the computation of a solution path with random data.
- Test_path_figure_random: Generate figures for solution paths on random data;
- Test_path_figure_UCI: Generate figures for solution paths on UCI data.

(4) Test4_truncated_CVaR_MM: Solvers for the nonconvex truncated CVaR-based sparse linear regression: MM+N-PPA, MM+Gurobi

- Test_MM_NPPA_UCI: MM+N-PPA for the nonconvex truncated CVaR-based model with UCI data; 
- Test_MM_NPPA_random: MM+N-PPA for the nonconvex truncated CVaR-based model with random data; 
- Test_MM_Gurobi_UCI: MM+Gurobi for the nonconvex truncated CVaR-based model with UCI data; 
- Test_MM_Gurobi_random: MM+Gurobi for the nonconvex truncated CVaR-based model with random data.

(5) genUCIdatafun: Convert the UCI data from LIBSVM format to MAT format.

(6) mainfun: All main functions for the above solvers.

(7) mexfun: two funtions in MEX format.

(8) solver: all the functions called by the above solvers.

(9) UCIdata: all availiable data for the UCI data and the high-accurate objective values.




