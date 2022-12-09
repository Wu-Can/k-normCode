%%=========================================================================
%% Test a solution path for the AS + N-ALM on random data
%% Input:
%% m_vec = sample size vector;
%% n_vec = feature size vector;
%% edgp_mat = the matrix whose each row is the number of the 
%%            spaced grid points of lambda for a given k
%%=========================================================================
clear all; clc;
rng('default');
HOME = pwd;
addpath(genpath(HOME));
pathname = [HOME,filesep,'Result_solution_path_figure_random',...
    filesep,'Result_AS_NALM_path_random',filesep]; 
if ~exist(pathname,'dir'),mkdir(pathname); end
addpath(genpath(pathname));
%%
%profile on
%================================== INPUT =================================
m_vec = [500 1000 3000];
n_vec = [100000 200000 500000];
edgp_mat = [15 30 45;15 30 45;30 60 90]; 
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
lengp = size(edgp_mat,2);
result = zeros(lenalp*len_m,10);
for pp = 1:len_m
    m = m_vec(pp);
    n = n_vec(pp);
    [A,b] = Generate_A_b_t_large(m,n,1);
    
    OPTIONS.tol = 1e-6;
    OPTIONS.m = m;
    OPTIONS.n = n;
    OPTIONS.sigmascale = 1.3;
    OPTIONS.tauscale = 1.2;
    for qq = 1:lengp
        edgp = edgp_mat(pp,qq);
        if pp == 1
            lambda_vec = [20:(8-20)/edgp:8;60:(45-60)/edgp:45;70:(55-70)/edgp:55];
        elseif pp == 2
            lambda_vec = [40:(25-40)/edgp:25;80:(65-80)/edgp:65;110:(95-110)/edgp:95];
        elseif pp == 3
            lambda_vec = [80:-30/edgp:50;160:-30/edgp:130;200:-30/edgp:170];
        end
        
        lenlam = size(lambda_vec,2);
        resulthist = zeros(10*lenalp,lenlam);
        result_time_path = zeros(lenalp,lenlam);
        result_xnnz_path = zeros(lenalp,lenlam);
        result_n_mean = zeros(lenalp,lenlam);
        
        for jj = 1:lenalp
            OPTIONS.lambda = lambda_vec(jj,:);
            OPTIONS.sigma = 20./OPTIONS.lambda;
            OPTIONS.tau = 20./OPTIONS.lambda;
            OPTIONS.kk = ceil((1-alpha_vec(jj))*m);
            
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
        eval(['filename_result = ''new_Resulthist_AS_SNIPAL_path_',num2str(m),'_',num2str(n),'_',num2str(edgp),''';']);
        
        save([pathname,filename_time],'result_time_path');
        save([pathname,filename_nnz],'result_xnnz_path');
        save([pathname,filename_mean],'result_n_mean');
        save([pathname,filename_result],'resulthist');
    end
end

% profile viewer





%%*********************************************************************
