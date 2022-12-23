# mexfun
This folder contains two C source files and one M file.

- `mexfun/mexAx.c`: compute matrix-vector products called by the main functions *ADMM_knorm.m*, *AS_NALM_path.m*, *dPPA_SSN.m*, *MM_Gurobi.m*, *MM_NPPA.m* and *NALM.m* in `mainfun`.
- `mexfun/OnceCD.c`: one loop for the coordinate descent algorithm called by the subfunction *CD.m* in `subfun`.
- `mexfun/BuildMex.m`: compile *mexAx.c* and *OnceCD.c* into MEX files. 
