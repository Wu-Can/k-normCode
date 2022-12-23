%%=========================================================================
%% Test a solution path for the warm-started N-ALM on random data
%% Input:
%% m_vec = sample size vector;
%% n_vec = feature size vector;
%% edgp_mat = the matrix whose each row is the number of the 
%%            spaced grid points lambda for a given k.
%%=========================================================================
clear all; clc;
rng('default');
HOME = pwd;
addpath(genpath(HOME));
pathname = [HOME,filesep,'Result_solution_path_figure_random',...
    filesep,'Result_warm_NALM_path_random',filesep];
if ~exist(pathname,'dir'),mkdir(pathname); end
addpath(genpath(pathname));
%%
%profile on
%================================== INPUT =================================
m_vec = [500 1000 3000];
n_vec = [100000 200000 500000];
edgp_vec = [15 30 45;15 30 45;30 60 90];
%==========================================================================
%% -------------------------Generate A and b--------------------------
len_m = length(m_vec);
len_n = length(n_vec);
if len_m ~= len_n
    fprintf('m_vec and n_vec must be the same length !');
    return;
end
alpha_vec = [0.9;0.5;0.1];
lenalp = length(alpha_vec);
lengp = size(edgp_vec,2);
for pp = 1:len_m
    m = m_vec(pp);
    n = n_vec(pp);
    [A,b] = Generate_A_b_t_large(m,n,1);
    OPTIONS.m = m;
    OPTIONS.n = n;
    OPTIONS.sigmascale = 1.3;
    OPTIONS.tauscale = 1.2;
    for qq  = 1:lengp
        edgp  = edgp_vec(pp,qq);
        if pp == 1
            lambda_vec = [20:(8-20)/edgp:8;60:(45-60)/edgp:45;70:(55-70)/edgp:55];
        elseif pp == 2
            lambda_vec = [40:(25-40)/edgp:25;80:(65-80)/edgp:65;110:(95-110)/edgp:95];
        elseif pp == 3
            lambda_vec = [80:-30/edgp:50;160:-30/edgp:130;200:-30/edgp:170];
        end
        
        lenlam = size(lambda_vec,2);
        resulthist = zeros(lenalp*8,lenlam);
        result_ttime_path = zeros(lenalp,lenlam);
        
        for jj = 1:lenalp
            OPTIONS.kk = ceil((1-alpha_vec(jj))*m);
            totletime = 0;
            for ii = 1:lenlam
                OPTIONS.lambda = lambda_vec(jj,ii);
                OPTIONS.sigma = 20/OPTIONS.lambda;
                OPTIONS.tau = 20/OPTIONS.lambda;
                if (ii > 1)
                    OPTIONS.x0 = solution_x;
                    OPTIONS.z0 = solution_z;
                    OPTIONS.u0 = solution_u;
                elseif (pp > 1)
                    field = {'x0','z0','u0'};
                    OPTIONS = rmfield(OPTIONS,field);
                end
                fprintf('\n ============* pp = %2d, qq = %2d, jj = %2d, ii = %2d *============',...
                    pp, qq, jj, ii);
                
                [~,~,~,~,~,info] = NALM(A,b,OPTIONS);
                
                totletime = totletime + info.totletime;
                result_ttime_path(jj,ii) = totletime;
                
                resulthist(1+(jj-1)*8,ii) = OPTIONS.lambda;
                resulthist(2+(jj-1)*8,ii) = info.xnnz;
                resulthist(3+(jj-1)*8,ii) =  info.res_kkt_final;
                resulthist(4+(jj-1)*8,ii) =  info.res_gap_final;
                resulthist(5+(jj-1)*8,ii) =  info.iter;
                resulthist(6+(jj-1)*8,ii) =  info.numSSN;
                resulthist(7+(jj-1)*8,ii) =  info.totletime;
                resulthist(8+(jj-1)*8,ii) =  totletime;
                
                solution_x = info.x;
                solution_z = info.z;
                solution_u = info.u;
            end
        end
        
        eval(['filename_time = ''new_Result_warm_NALM_time_',num2str(m),'_',num2str(n),'_',num2str(edgp),''';']);
        eval(['filename_result = ''new_Resulthist_warm_NALM_',num2str(m),'_',num2str(n),'_',num2str(edgp),''';']);
        save([pathname,filename_time],'result_ttime_path');
        save([pathname,filename_result],'resulthist');
    end
end

%profile viewer



%%*********************************************************************
