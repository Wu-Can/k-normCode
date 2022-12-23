# mainfun
This folder contains twelve main functions and their main subfunctions.
- `mainfun/NALM.m`: the main function of the solver ***N-ALM***.
- `mainfun/SIRPN_CD.m`: the main function of the solver ***S-IRPN***.
- `mainfun/ADMM_knorm.m`: the main function of the solver ***ADMM***.
- `mainfun/PSG_solver.m`: the main function of the solver ***PSG solvers***.
- `mainfun/Gurobi_knorm_lp.m`: the main function of the solver ***Gurobi***.
- `mainfun/AS_NALM_path.m`:  the main function of the solver ***AS+N-ALM***.
- `mainfun/MM_NPPA.m`:  the main function of the solver ***MM+N-PPA***.
- `mainfun/dPPA_SSN.m`: the main subfunction of *MM_NPPA.m*.
- `mainfun/MM_Gurobi.m`: the main function of the solver ***MM+Gurobi***.
- `mainfun/qp_Gurobi.m`: the main subfunction of *MM_Gurobi.m*.
- `mainfun/Cross_validation_CVaR.m`: the main function called by *Test_convex_CVaR_model_random.m* in `Test4_Risk-averse_model`.
- `mainfun/Cross_validation_truncated_CVaR.m`: the main function called by *Test_truncated_CVaR_model_random.m* in `Test4_Risk-averse_model`.

All main functions are called directly by the scripts in the folders `Test1_fixed_lambda`, `Test2_Solution_path`, `Test3_truncated_CVaR_MM` and `Test4_Risk-averse_model`, respectively.
