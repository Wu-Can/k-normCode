%%=========================================================================
%% Test the ADMM for the convex CVaR-based models with random data
%% Input:
%% prob = the order number of problems;
%% flag_tol = 0, adopt eta_res < tol as the stopping criterion;                      
%%          = 1, adopt relobj < tol as the stopping criterion.
%%=========================================================================
clear all; clc;
rng('default');
HOME = pwd;
addpath(genpath(HOME));
%%
%profile on
%========================= INPUT =========================
prob = 1;%[1:7];
flag_tol = 1;
%=========================================================
if (flag_tol ~= 0) && (flag_tol ~= 1) 
    fprintf('The value of flag_tol must be 0 or 1!');
    return;
end

%% Input the objective values with the tolerance 1e-9 if flag_tol = 1
if flag_tol == 1
    filepath = fileparts(HOME);
    datadir = [filepath,filesep,'UCIdata'];
    addpath(genpath(datadir));
    fname{1} = 'Result_Gorubi_random_123456_1e-9';
    fname{2} = 'Resulthist_SNIPAL_random_3000_500000';
    for rr = 1:2
        resultname = [datadir,filesep,fname{rr}];
        if exist([resultname,'.mat'])
            load([resultname,'.mat'])
            if rr == 1
                obj_opt_vec1 = result(:,5); clear result;
            else
                obj_opt_vec2 = result(:,end); clear result;
            end
        else
            fprintf('\n Can not find the file');
            fprintf('\n ');
            return;
        end
    end
end

alpha_vec = [0.9;0.5;0.1];
lenprob = length(prob);
lenalp = length(alpha_vec);
result = zeros(lenalp*lenprob,11);
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
    
    OPTIONS.tol = 1e-3;
    OPTIONS.tau_ADMM = 1.618;
    OPTIONS.gamma = 1;
    OPTIONS.m = m;
    OPTIONS.n = n;
    OPTIONS.flag_tol = flag_tol;
    lamb = 0.12;
    for jj = 1:lenalp
        alpha = alpha_vec(jj);
        if flag_tol == 1
            if alpha == 0.9
                j = 1;
            elseif alpha == 0.5
                j = 2;
            else
                j = 3;
            end
            if i == 7
                OPTIONS.obj_opt = obj_opt_vec2(j);
            else
                OPTIONS.obj_opt = obj_opt_vec1(j+(i-1)*3);
            end
        end
         
        OPTIONS.kk = ceil((1-alpha)*m);
        OPTIONS.lambda = lamb*OPTIONS.kk;
        [obj,~,~,runhist,info] = ADMM_knorm(A,b,OPTIONS);
        
        result(jj+(ii-1)*lenalp,1) =  OPTIONS.m;
        result(jj+(ii-1)*lenalp,2) =  OPTIONS.n; 
        result(jj+(ii-1)*lenalp,3) =  OPTIONS.lambda;
        result(jj+(ii-1)*lenalp,4) =  OPTIONS.kk;
        result(jj+(ii-1)*lenalp,5) = info.xnnz;
        result(jj+(ii-1)*lenalp,6) = info.eta_res;
        result(jj+(ii-1)*lenalp,7) = info.res_kkt;
        result(jj+(ii-1)*lenalp,8) = info.iter;
        result(jj+(ii-1)*lenalp,9) = info.time;
        result(jj+(ii-1)*lenalp,10) = obj(1);
        if OPTIONS.flag_tol == 1
            result(jj+(ii-1)*lenalp,11) = info.relobj;
        end
    end
end

save Result_ADMM_knorm_random.mat result
%profile viewer




%%*********************************************************************
