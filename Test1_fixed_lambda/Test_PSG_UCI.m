%%=========================================================================
%% Test the PSG solvers for the convex CVaR-based models under a 
%% fixed lambda with UCI data
%% INPUT:
%% prob = the order number of problems;
%% SOLVER = one of seven solvers in the PSG package.                       
%%=========================================================================
clear all; clc;
rng('default');
%%
%================================== INPUT =================================
prob = [7:9];
SOLVER = 'VAN';%{'VAN','TANK','CAR','BULDOZER','VANGRB','CARGRB','HELI'};
%==========================================================================

HOME = pwd;
addpath(genpath(HOME));
filepath = fileparts(HOME);
datadir = [filepath,filesep,'UCIdata']; 

fname{1} = 'triazines_scale_expanded4'; 
fname{2} = 'pyrim_scale_expanded5'; 
fname{3} = 'log1p.E2006.test';
fname{4} = 'bodyfat_scale_expanded7'; 
fname{5} = 'housing_scale_expanded7'; 
fname{6} = 'log1p.E2006.train';
fname{7} = 'mpg_scale_expanded7'; 
fname{8} = 'abalone_scale_expanded7'; 
fname{9} = 'space_ga_scale_expanded9'; 
fname{10} = 'E2006.test';
fname{11} = 'E2006.train';

alpha_vec = [0.9 0.5 0.1];
lenprob = length(prob);
lenalp = length(alpha_vec);
result = zeros(lenprob*lenalp,8);
result_status = cell(lenprob*lenalp,1);
for ii = 1:lenprob
    i = prob(ii);
    probname = [datadir,filesep,fname{i}];
    fprintf('\n Problem name: %s', fname{i});
    if exist([probname,'.mat'])
        load([probname,'.mat'])
    else
        fprintf('\n Can not find the file in UCIdata');
        fprintf('\n ');
        return;
    end
    [m,n] = size(A);
    if i == 1
        lamc = 1e-4;
    elseif i == 2
        lamc = 1e-5;
    elseif i == 3
        lamc = 1e-6;
    elseif i <= 7
        lamc = 1e-7;
    elseif i <= 9
        lamc = 1e-8;
    elseif i == 10
        lamc = 1e-9;
    elseif i == 11
        lamc = 1e-10;
    end
    lamcmaxbTA = lamc*max(abs(b'*A));
    abs_data = ones(1,n);
    for jj = 1:lenalp
        alpha = alpha_vec(jj);
        kk = ceil((1-alpha)*m);
        lambda = kk*lamcmaxbTA;
        alpha = 1 - kk/m;
        
        [solution_str,outargstruc_arr] = PSG_solver(b,A,alpha,abs_data,kk,lambda,SOLVER);
        
        output_structure = tbpsg_solution_struct(solution_str, outargstruc_arr);
        point_data = tbpsg_optimal_point_data(solution_str, outargstruc_arr);
        %----------------- compute nnz of x --------------
        [sortx,index_x] = sort(abs(point_data),'descend');
        normx1 = 0.999*norm(point_data,1);
        tmpidex = find(cumsum(sortx) > normx1);
        if isempty(tmpidex)
            nnzx = 0;
        else
            nnzx = tmpidex(1);
        end
        %------------------------------------------------
        result(jj + (ii-1)*lenalp,1) = kk;
        result(jj + (ii-1)*lenalp,2) = lamc;
        result(jj + (ii-1)*lenalp,3) = output_structure.objective;
        result(jj + (ii-1)*lenalp,4) = output_structure.gap;
        result(jj + (ii-1)*lenalp,5) = output_structure.time(1); % data loading time
        result(jj + (ii-1)*lenalp,6) = output_structure.time(2); % preprocessing time
        result(jj + (ii-1)*lenalp,7) = output_structure.time(3); % solving time
        result(jj + (ii-1)*lenalp,8) = nnzx;
        
        result_status(jj + (ii-1)*lenalp) = output_structure.status;
    end
end

save result_PSG_UCI.mat result



