%%============================================================================
%% Test the S-IRPN for the convex CVaR-based models with random data
%%============================================================================
clear all; clc;
rng('default');
HOME = pwd;
addpath(genpath(HOME));
filepath = fileparts(HOME);
datadir = [filepath,filesep,'UCIdata']; 
addpath(genpath(datadir));
%%
%profile on
%================== INPUT ================== 
Prob = 1;%[1:7];
%===========================================

%% Input the objective values with the tolerance 1e-9
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

alpha_vec = [0.9;0.5;0.1];
lenprob = length(Prob);
lenalp = length(alpha_vec);
result = zeros(lenalp*lenprob,10);

for pp = 1:lenprob
    i = Prob(pp);
    switch i
        case 1
            n = 300; d = 1000;
        case 2
            n = 500; d = 3000;
        case 3
            n = 800; d = 8000;
        case 4
            n = 1000; d = 10000;
        case 5
            n = 300; d = 50000;
        case 6
            n = 1000; d = 100000;
        case 7
            n = 3000; d = 500000;
    end
    
    %% generate A and b
    [A,b] = Generate_A_b_t_large(n,d,1);
    
    lamb = 0.12;
    OPTIONS.tol = 1e-3;
    OPTIONS.n = n;
    OPTIONS.d = d;
    OPTIONS.epsilon = 1;
    OPTIONS.tolsub = 1e-2;
    OPTIONS.mu = 1e-5;
    OPTIONS.maxiter = 9;
    OPTIONS.maxitersub = 200;
    OPTIONS.maxtime = 7200;
    
    for jj = 1:lenalp
        alpha = alpha_vec(jj);
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
        OPTIONS.kk = ceil((1-alpha)*n);
        OPTIONS.lambda = lamb;
        
        [obj,x,t,runhist,info] = SIRPN_CD(A,b,OPTIONS);
        
        result(jj+(pp-1)*lenalp,1) = OPTIONS.n;
        result(jj+(pp-1)*lenalp,2) = OPTIONS.d;
        result(jj+(pp-1)*lenalp,3) = OPTIONS.lambda;
        result(jj+(pp-1)*lenalp,4) = OPTIONS.kk;
        result(jj+(pp-1)*lenalp,5) = info.xnnz;
        result(jj+(pp-1)*lenalp,6) = info.relobj;
        result(jj+(pp-1)*lenalp,7) = info.iter;
        result(jj+(pp-1)*lenalp,8) = info.itersub;
        result(jj+(pp-1)*lenalp,9) = info.numiterCD;
        result(jj+(pp-1)*lenalp,10) = info.time;
        result(jj+(pp-1)*lenalp,11) = obj;
    end
end

save Result_SIRPN_random.mat result
%profile viewer




%%*********************************************************************
