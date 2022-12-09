%%=================================================================================================
%% Test the N-ALM for the convex CVaR-based models under a fixed lambda with UCI data
%% INPUT:
%% prob = the order number of problems;
%% flag_tol = 0, adopt eta_res < tol as the stopping criterion of the N-ALM;                         
%%          = 1, adopt  relobj < tol as the stopping criterion of the N-ALM;
%%          = 2, adopt  relkkt < tol as the stopping criterion of the N-ALM;
%% tol = the tolerance of the N-ALM.
%% alpha_vec = the vector of the confident level alpha satisfying k:=ceil((1-alpha)*m)
%%=================================================================================================
clear all; clc;
rng('default');
%%
%profile on
%================================== INPUT =================================
prob = [7 8 9];%[1:11];
flag_tol = 0;
tol = 1e-6;
alpha_vec = [0.9 0.5 0.1];
%==========================================================================

%eval(['diary ''Diary_NALM_UCI_flagtol_',num2str(flag_tol),'_tol_',num2str(tol),'.txt''']);
if (flag_tol ~= 0) && (flag_tol ~= 1) && (flag_tol ~= 2)
    fprintf('The value of flag_tol must be 0, 1 or 2 !');
    return;
end
if flag_tol == 1
    load Result_Gorubi_UCI_1_11_1e-9.mat
    obj_opt_vec = result(:,3);
end

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

%%
numProb =  length(prob);
lenalp = length(alpha_vec);
result = zeros(numProb*lenalp,9);
for ii = 1:numProb
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
    
    OPTIONS.tol = tol;
    OPTIONS.maxiter = 200;
    OPTIONS.maxitersub = 200; 
    OPTIONS.maxtime = 7200;
    OPTIONS.m = m;
    OPTIONS.n = n;
    OPTIONS.sigmascale = 1.3;
    OPTIONS.tauscale = 1.2;
    OPTIONS.flag_tol = flag_tol;
    maxbTA = lamc*max(abs(b'*A));
    
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
            OPTIONS.obj_opt = obj_opt_vec(j+(i-1)*3);
        end
        OPTIONS.kk = ceil((1-alpha)*m);
        OPTIONS.lambda = maxbTA*OPTIONS.kk;
        OPTIONS.sigma = 1e-4/(lamc*OPTIONS.kk);
        OPTIONS.tau = 1e-4/(lamc*OPTIONS.kk);
        
        [obj,~,~,~,runhist,info]= NALM(A,b,OPTIONS);
        
        result(jj+(ii-1)*lenalp,1) =  lamc;
        result(jj+(ii-1)*lenalp,2) =  OPTIONS.kk;
        result(jj+(ii-1)*lenalp,3) = info.xnnz;
        result(jj+(ii-1)*lenalp,4) = info.eta_res;
        result(jj+(ii-1)*lenalp,5) = info.res_kkt_final;
        result(jj+(ii-1)*lenalp,6) = info.iter;
        result(jj+(ii-1)*lenalp,7) = info.numSSN;
        result(jj+(ii-1)*lenalp,8) = info.totletime;
        result(jj+(ii-1)*lenalp,9) = obj(1);
        if flag_tol == 1
            result(jj+(ii-1)*lenalp,10) = info.relobj;
        end
    end
end

eval(['save Result_NALM_UCI_flagtol_',num2str(flag_tol),'_tol_',num2str(tol),'.mat result']);

%profile viewer
%diary off




%%*********************************************************************