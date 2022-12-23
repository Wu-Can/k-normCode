%%===============================================================
%% Test a solution path for the AS + N-ALM on UCI data
%%===============================================================

clear all; clc;
rng('default');
%%
%profile on
%========== INPUT ==========
prob = [9 4 2];
%===========================
%% -------------------------input A and b------------------
HOME = pwd;
addpath(genpath(HOME));
filepath = fileparts(HOME);
datadir = [filepath,filesep,'UCIdata'];
pathname = [HOME,filesep,'Result_solution_path_figure_UCI',...
    filesep,'Result_AS_NALM_path_UCI',filesep];
if ~exist(pathname,'dir'),mkdir(pathname); end
addpath(genpath(pathname));

fname{1} = 'E2006.train';
fname{2} = 'log1p.E2006.train';
fname{3} = 'E2006.test';
fname{4} = 'log1p.E2006.test';
fname{5} = 'pyrim_scale_expanded5';
fname{6} = 'triazines_scale_expanded4';
fname{7} = 'abalone_scale_expanded7';
fname{8} = 'bodyfat_scale_expanded7';
fname{9} = 'housing_scale_expanded7';
fname{10} = 'mpg_scale_expanded7';
fname{11} = 'space_ga_scale_expanded9';

alpha_vec = [0.9 0.5 0.1];
lenprob = length(prob);
lenalp = length(alpha_vec);
resulthist = zeros(lenalp*lenprob,10);
for pp = 1:lenprob
    i = prob(pp);
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
    norm_bTA_inf = max(abs(b'*A));
    OPTIONS.tol = 1e-6;
    OPTIONS.m = m;
    OPTIONS.n = n;
    OPTIONS.sigmascale = 1.3;
    OPTIONS.tauscale = 1.2;
    OPTIONS.maxitersub = 100;
    if (i == 9)
        edgp_vec = [20 30 40];
    elseif (i == 4)
        edgp_vec = [30 60 90];
    elseif (i == 2)
        edgp_vec = [100 200 300];
    end
    
    lengp = length(edgp_vec);
    for qq = 1:lengp
        edgp = edgp_vec(qq);
        if i == 9
            lambda_vec = [4.5:(1.5-4.5)/edgp:1.5;4.5:(1.5-4.5)/edgp:1.5;4.5:(1.5-4.5)/edgp:1.5];
        elseif i == 4
            lambda_vec = [45:(20-45)/edgp:20;75:(45-75)/edgp:45;85:(60-85)/edgp:60];
        elseif i == 2
            lambda_vec = [140:-100/edgp:40;300:-200/edgp:100;330:-200/edgp:130];
        end
        
        lenlam = size(lambda_vec,2);
        resulthist = zeros(10*lenalp,lenlam);
        result_time_path = zeros(lenalp,lenlam);
        result_xnnz_path = zeros(lenalp,lenlam);
        result_n_mean = zeros(lenalp,lenlam);
        
        for jj = 1:lenalp
            OPTIONS.kk = ceil((1-alpha_vec(jj))*m);
            OPTIONS.lambda = lambda_vec(jj,:);
            lamc_vec = OPTIONS.lambda./norm_bTA_inf;
            OPTIONS.sigma = 1e-4./lamc_vec;
            OPTIONS.tau = 1e-4./lamc_vec;
            
            [~,~,~,runhistAS,~] = AS_NALM_path(A,b,OPTIONS);
            result_time_path(jj,:) = runhistAS.ttime_path;
            result_xnnz_path(jj,:) = runhistAS.xnnz_path;
            result_n_mean(jj,:) = runhistAS.n_mean;
            
            resulthist(1+(jj-1)*10,:) =  OPTIONS.lambda;
            resulthist(2+(jj-1)*10,:) =  runhistAS.relkkt;
            resulthist(3+(jj-1)*10,:) =  runhistAS.relgap;
            resulthist(4+(jj-1)*10,:) =  runhistAS.iterAS;
            resulthist(5+(jj-1)*10,:) =  runhistAS.iterPAL;
            resulthist(6+(jj-1)*10,:) =  runhistAS.iterSSN;
            resulthist(7+(jj-1)*10,:) =  runhistAS.time;
            resulthist(8+(jj-1)*10,:) = runhistAS.xnnz_path;
            resulthist(9+(jj-1)*10,:) = runhistAS.n_mean;
            resulthist(10+(jj-1)*10,:) = runhistAS.ttime_path;
        end
        
        eval(['filename_time = ''new_Result_AS_NALM_path_time_',num2str(m),'_',num2str(n),'_',num2str(edgp),''';']);
        eval(['filename_nnz = ''new_Result_AS_NALM_path_nnz_',num2str(m),'_',num2str(n),'_',num2str(edgp),''';']);
        eval(['filename_mean = ''new_Result_AS_NALM_path_n_mean_',num2str(m),'_',num2str(n),'_',num2str(edgp),''';']);
        eval(['filename_result = ''new_Resulthist_AS_NALM_path_',num2str(m),'_',num2str(n),'_',num2str(edgp),''';']);
        save([pathname,filename_time],'result_time_path');
        save([pathname,filename_nnz],'result_xnnz_path');
        save([pathname,filename_mean],'result_n_mean');
        save([pathname,filename_result],'resulthist');
    end
end

% profile viewer




%%*********************************************************************
