%%*************************************************************************
%% PSG_solver: run the mpsg_solver for solving the convex CVaR-based
%% sparse linear regression:
%%
%% (P)   minimize_{x in R^n} {||A*x - b||_(k) + lambda ||x||_1}
%%
%% [solution_str,outargstruc_arr] = PSG_solver(b,A,alpha,abs_data,...
%%        para_kk,para_lambda,ABCDEF)
%%
%% Input:
%% A, b = the design matrix A and the response vector b in (P)
%% alpha = the value of alpha satisfying k=(1-alpha)*m with the 
%%        sample size m
%% abs_data = the n-dimensional vector of all ones
%% para_kk = the value of parameter k in (P)
%% para_lambda = the value of parameter lambda in (P)
%% ABCDEF = one of solvers {'VAN','TANK','CAR','BULDOZER','VANGRB',
%%         'CARGRB','HELI'}
%% Output: 
%% output_structure = tbpsg_solution_struct(solution_str,
%%                    outargstruc_arr);
%% point_data = tbpsg_optimal_point_data(solution_str, 
%%              outargstruc_arr);
%% output_structure.objective = the output objective value of (P)
%% output_structure.gap = difference between objective value in 
%%                        obtained point and lower estimate of 
%%                        optimal value
%% output_structure.time(1) = data loading time
%% output_structure.time(2) = preprocessing time
%% output_structure.time(3) = solving time
%% output_structure.status = status of solution of (P) 'optimal' 
%%                           or 'feasible' or 'infeasible' or 
%%                           'unbounded' or 'calculated'
%% point_data = the output solution of (P)
%% For more details, please see the online manual of PSG at 
%% http://www.aorda.com/html/PSG_Help_HTML/index.html?why_psg.htm
%%*************************************************************************
function [solution_str,outargstruc_arr] = PSG_solver(b,A,alpha,abs_data,para_kk,para_lambda,ABCDEF)
%Pack data to PSG structure:
iargstruc_arr(1) = matrix_pack('matrix_CVaR_abs',A,[],b,[]);
iargstruc_arr(2) = matrix_pack('matrix_polynom_abs',abs_data,[],[]);

%Define problem statement:
problem_statement = sprintf('%s\n',...
'minimize',...
'Objective: linearize = 0',...
'para_kk*cvar_risk(alpha, abs(matrix_CVaR_abs))',...
'+ para_lambda*polynom_abs(matrix_polynom_abs)',...
'Solver: ABCDEF, precision = 6, stages = 6, timelimit = 7200');

%Solvers: VAN,TANK,CAR,BULDOZER,VANGRB,CARGRB,HELI

%Change parameters with their values in problem statement:
problem_statement = strrep(problem_statement,'alpha',num2str(alpha));
problem_statement = strrep(problem_statement,'para_kk',num2str(para_kk));
problem_statement = strrep(problem_statement,'para_lambda',num2str(para_lambda));
problem_statement = strrep(problem_statement,'ABCDEF',ABCDEF);

%Optimize problem using mpsg_solver function:
[solution_str, outargstruc_arr] = mpsg_solver(problem_statement, iargstruc_arr);
results = psg_convert_outargstruc(solution_str,outargstruc_arr);
outargstruc_arr = results.data;
solution_str = results.string;

%=======================================================================|
%American Optimal Decisions, Inc. Copyright                             |
%Copyright 〢merican Optimal Decisions, Inc. 2007-2016.                 |
%American Optimal Decisions (AOD) retains copyrights to this material.  |
%                                                                       |
%Permission to reproduce this document and to prepare derivative works  |
%from this document for internal use is granted, provided the copyright |
%and 揘o Warranty�  statements are included with all reproductions       |
%and derivative works.                                                  |
%                                                                       |
%For information regarding external or commercial use of copyrighted    |
%materials owned by AOD, contact AOD at support@aorda.com.              |
%=======================================================================|
 

