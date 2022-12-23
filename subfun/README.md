# subfun
This folder contains 22 subfunctions of the corresponding main functions in `mainfun`, 5 functions for generating random data as well as another subfunction.

### 22 subfunctions of the corresponding main functions in `mainfun`

#### Subfunctions called by the main function *NALM.m* include: 
  - *CholHess.m*:  Cholesky factorization of Hessian
  - *Cumpute_matrix_C.m*: calculate a matrix $C$ as a part of Hessian
  - *findstep*: find a step length by the line search
  - *Generatedash_Jacobi*: obtain some information of the generalized Jacobian matrix
  - *invDvec.m*: calculate $D^{-1}\mbox{rhs}$ and $T$ defined in Section 3.3 of the paper
  - *linsysolvefun.m*: solve a linear system
  - *matvec_N_ALM.m*: compute Hessian-vector products
  - *Proj_dknorm.m*: compute the projection on the ball with center 0 and radius $r$ in the sense of the dual norm of k-norm
  - *Proj_inf.m*: compute the projection on the ball with center 0 and radius $r$ in the sense of $\ell_{\infty}$ norm
  - *psqmr_knorm_N_ALM.m*: solve a Newton linear system by the preconditioned symmetric QMR
#### Subfunctions called by the main function *SIRPN_CD.m* include: 
  - *CD.m*: coordinate descent algorithm for solving subproblems of ***S-IRPN***.
  - *phi_eps_fun.m*: calculate the function value, first-order derivative, second-order derivative of $\phi_{\epsilon}$ defined in Appendix E of the supplementary materials
  - *psi_eps_fun.m*: calculate the function value, first-order derivative, second-order derivative of $\psi_{\epsilon}$ defined in Appendix E of the supplementary materials
  - *Proj_inf.m*: compute the projection on the ball with center 0 and radius $r$ in the sense of $\ell_{\infty}$ norm
#### Subfunctions called by the main function *ADMM_knorm.m* include: 
  - *CholHess.m*:  Cholesky factorization of Hessian
  - *linsysolvefun.m*: solve a linear system 
  - *Matvecmu.m*: compute a matrix times a vector when the sample size is not greater that the feature size
  - *Matvecnx.m*: compute a matrix times a vector when the sample size is greater that the feature size
  - *Proj_dknorm.m*: compute the projection on the ball with center 0 and radius $r$ in the sense of the dual norm of k-norm
  - *Proj_inf.m*: compute the projection on the ball with center 0 and radius $r$ in the sense of $\ell_{\infty}$ norm
  - *psqmr_ADMM_new.m*: solve a linear system of the ADMM subproblems by the preconditioned symmetric QMR
#### Subfunctions called by the main function *Gurobi_knorm_lp.m* include: 
  - *linprog.m*: solve LP using the Gurobi MATLAB interface
#### Subfunctions called by the main function *AS_NALM_path.m* include: 
  - *Proj_dknorm.m*: compute the projection on the ball with center 0 and radius $r$ in the sense of the dual norm of k-norm
  - *Proj_inf.m*: compute the projection on the ball with center 0 and radius $r$ in the sense of $\ell_{\infty}$ norm
#### Subfunctions called by the main subfunction *dPPA_SSN.m* of *MM_NPPA.m* include: 
  - *CholHess.m*:  Cholesky factorization of Hessian 
  - *Cumpute_matrix_W_MM.m*: calculate an element $W$ in the Clarke generalized Jacobian of the proximal mapping with respect to k-norm
  - *findstep_new_MM.m*: find a step length by the line search
  - Generatedash_Jacobi_MM.m: obtain some information of the generalized Jacobian matrix
  - *invDvec_MM.m*: calculate $D^{-1}\mbox{rhs}$ and $T$ 
  - *linsysolvefun.m*: solve a linear system
  - *Proj_dknorm.m*: compute the projection on the ball with center 0 and radius $r$ in the sense of the dual norm of k-norm
  - *Proj_inf.m*: compute the projection on the ball with center 0 and radius $r$ in the sense of $\ell_{\infty}$ norm
  
### 5 functions for generating random data based on Table 1 or Table 3 in the supplementary materials
  - *Generate_A_b_t.m*: generate the small-scale random data based on Table 1 
  - *Generate_A_b_t_large.m*: generate the large-scale random data based on Table 1 
  - *Generate_A_b_Configuration.m*: generate the small-scale random data based on Table 3
  - *Generate_A_b_Conf_large.m*: generate the large-scale random data based on Table 3 
  - *Generate_Toeplitz_matrix.m*: generate $n$-dimensional Toeplitz covariance matrix
  
### Another subfunction:
  - *HKLeq.m*: subfunction of *Proj_dknorm.m*.
 
