function [solution_str,outargstruc_arr] = PSG_solver(b,A,alpha,abs_data,para_kk,para_lambda,ABCDEF)
%% This function is used to run the mpsg_solver for solving the the convex CVaR-based models

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
 

