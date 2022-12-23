%%================================================================================================
%% Test the barrier method in Gurobi for the convex CVaR-based models with random data
%% Input:
%% m_vec = sample size vector;
%% n_vec = feature size vector;
%%================================================================================================
clear all; clc;
rng('default');
HOME = pwd;
addpath(genpath(HOME));
%%
%profile on
%================================== INPUT =================================
m_vec = 300;%[300,500,800,1000,300,1000,3000];
n_vec = 1000;%[1000,3000,8000,10000,50000,100000,500000];
%==========================================================================
%% -------------------------Generate A and b------------------------
len_m = length(m_vec);
len_n = length(n_vec);
if len_m ~= len_n
   fprintf('m_vec and n_vec must be the same length !'); 
   return;
end
alpha_vec = [0.9;0.5;0.1];
lenalp = length(alpha_vec);
result = zeros(lenalp*len_m,10);
for pp = 1:len_m
    m = m_vec(pp);
    n = n_vec(pp);
   [A,b] = Generate_A_b_t_large(m,n,1);
    
    OPTIONS.m = m;
    OPTIONS.n = n;
    lamb = 0.12;
    
    for jj = 1:lenalp
        alpha = alpha_vec(jj);
        OPTIONS.kk = ceil((1-alpha)*m);
        OPTIONS.lambda = lamb*OPTIONS.kk;
        
        [sol,fval_bar,exitflag_bar,output_bar,~] = Gurobi_knorm_lp(A,b,OPTIONS);
        [sortx,~] = sort(abs(sol.x),'descend');
        normx1 = 0.999*norm(sol.x,1);
        tmpidex = find(cumsum(sortx) > normx1);
        if isempty(tmpidex)
            nnzeros = 0;
        else
            nnzeros = tmpidex(1);
        end
        
        result(jj+(pp-1)*lenalp,1) = OPTIONS.m;
        result(jj+(pp-1)*lenalp,2) = OPTIONS.n;
        result(jj+(pp-1)*lenalp,3) = OPTIONS.lambda;
        result(jj+(pp-1)*lenalp,4) = OPTIONS.kk;
        result(jj+(pp-1)*lenalp,5) = fval_bar;
        result(jj+(pp-1)*lenalp,6) = output_bar.baritercount;
        result(jj+(pp-1)*lenalp,7) = output_bar.time;
        result(jj+(pp-1)*lenalp,8) = exitflag_bar;
        result(jj+(pp-1)*lenalp,9) = nnzeros;
        result(jj+(pp-1)*lenalp,10) = output_bar.constrviolation;
    end
end

save Result_Gorubi_random.mat result
% profile viewer

%%*********************************************************************
