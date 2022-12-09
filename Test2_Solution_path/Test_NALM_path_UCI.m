%%================================================================
%% Test a solution path for the pure N-ALM on UCI data
%%================================================================
clear all; clc;
rng('default');
HOME = pwd;
addpath(genpath(HOME));
filepath = fileparts(HOME);
datadir = [filepath,filesep,'UCIdata'];
pathname = [HOME,filesep,'Result_solution_path_figure_UCI',...
    filesep,'Result_NALM_path_UCI',filesep];
if ~exist(pathname,'dir'),mkdir(pathname); end
addpath(genpath(pathname));
%%
%profile on
%========= INPUT =========
prob = [9 4 2];
%=========================
%% -------------------------input A and b------------------------
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
    OPTIONS.m = m;
    OPTIONS.n = n;
    OPTIONS.sigmascale = 1.3;
    OPTIONS.tauscale = 1.2;
    norm_bTA_inf = max(abs(b'*A));
    if (i == 9)
        edgp_vec = [20 30 40];
    elseif (i == 4)
        edgp_vec = [30 60 90];
    elseif (i == 2)
        edgp_vec = [100 200 300];
    end
    %%
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
        resulthist = zeros(lenalp*8,lenlam);
        result_ttime_path = zeros(lenalp,lenlam);
        
        for jj = 1:lenalp
            OPTIONS.kk = ceil((1-alpha_vec(jj))*m);
            totletime = 0;
            for ii = 1:lenlam
                OPTIONS.lambda = lambda_vec(jj,ii);
                lamc = OPTIONS.lambda/norm_bTA_inf;
                OPTIONS.sigma = 1e-4/lamc;
                OPTIONS.tau = 1e-4/lamc;
                fprintf('\n ============* pp = %2d, qq = %2d, jj = %2d, ii = %2d *============',...
                    pp, qq, jj, ii);
                
                [obj,~,~,~,runhist,info] = NALM(A,b,OPTIONS);
                
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
            end
        end
        eval(['filename_time = ''new_Result_NALM_path_time_',num2str(m),'_',num2str(n),'_',num2str(edgp),''';']);
        eval(['filename_result = ''new_Resulthist_NALM_path_',num2str(m),'_',num2str(n),'_',num2str(edgp),''';']);
        save([pathname,filename_time],'result_ttime_path');
        save([pathname,filename_result],'resulthist');
    end
end

%profile viewer



%%*********************************************************************
