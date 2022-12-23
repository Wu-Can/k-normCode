%%=========================================================================
%% Test the PSG solvers for the convex CVaR-based models under a 
%% fixed lambda with random data
%% INPUT:
%% prob = the order number of problems;
%% SOLVER = solvers in the PSG package.                        
%%=========================================================================
clear all; clc;
rng('default');
HOME = pwd;
addpath(genpath(HOME));
%%
%profile on
%================================== INPUT =================================
prob = 1;%[1:7];
solvers = {'VAN'};%{ 'VAN','TANK','CAR','BULDOZER','VANGRB','CARGRB','HELI'};
%==========================================================================

HOME = pwd;
addpath(genpath(HOME));
datadir = 'D:\k_norm_code\data\AS_SNIPAL_random_data';

alpha_vec = [0.9, 0.5, 0.1];
lenprob = length(prob);
lensol = length(solvers);
lenalp = length(alpha_vec);
result = zeros(lenprob*lensol*lenalp,10);
result_status = cell(lenprob*lensol*lenalp,1);

for ii = 1:lenprob
    i = prob(ii);
    switch i
        case 1
            m = 300; n = 1000;
        case 2
            m = 500; n = 3000;
        case 3
            m = 800; n = 8000;
        case 4
            m = 1000; n = 10000;
        case 5
            m = 300; n = 50000;
        case 6
            m = 1000; n = 100000;
        case 7
            m = 3000; n = 500000;
    end
    
    %% Generate A and b
    [A,b] = Generate_A_b_t_large(m,n,1);
    
    lamb = 0.12;
    abs_data = ones(1,n);
    kk_vec = ceil((1-alpha_vec)*m);
    alpha_tmp = 1 - kk_vec/m;
    for pp = 1:lensol
        SOLVER = solvers{pp};
        for jj = 1:lenalp
            alpha = alpha_tmp(jj);
            kk = kk_vec(jj);
            lambda = lamb*kk;
            
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
            
            result(jj + (pp-1)*lenalp + (ii-1)*lenalp*lensol,1) = m;
            result(jj + (pp-1)*lenalp + (ii-1)*lenalp*lensol,2) = n;            
            result(jj + (pp-1)*lenalp + (ii-1)*lenalp*lensol,3) = lambda;
            result(jj + (pp-1)*lenalp + (ii-1)*lenalp*lensol,4) = kk;
            result(jj + (pp-1)*lenalp + (ii-1)*lenalp*lensol,5) = output_structure.objective;
            result(jj + (pp-1)*lenalp + (ii-1)*lenalp*lensol,6) = output_structure.gap;
            result(jj + (pp-1)*lenalp + (ii-1)*lenalp*lensol,7) = output_structure.time(1);
            result(jj + (pp-1)*lenalp + (ii-1)*lenalp*lensol,8) = output_structure.time(2);
            result(jj + (pp-1)*lenalp + (ii-1)*lenalp*lensol,9) = output_structure.time(3);
            result(jj + (pp-1)*lenalp + (ii-1)*lenalp*lensol,10) = nnzx;
            result_status(jj + (pp-1)*lenalp + (ii-1)*lenalp*lensol) = output_structure.status;
        end
    end
end

save result_PSG_random.mat result
% save result_PSG_status_random.mat result_status


