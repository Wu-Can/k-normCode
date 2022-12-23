%%==============================================================================================
%% Test the N-ALM for the convex CVaR-based models with random data
%% Input:
%% prob = the order number of problems;
%% flag_tol = 0, adopt eta_res < tol as the stopping criterion of the N-ALM;                         
%%          = 1, adopt  relobj < tol as the stopping criterion of the N-ALM;
%%          = 2, adopt  relkkt < tol as the stopping criterion of the N-ALM.
%% tol = the tolerance of the N-ALM. 
%% alpha_vec = the vector of the confidence level alpha such that k:=ceil((1-alpha)*m)
%% flag_J = 1, test the influence of the parameter k on the cardinality s of
%%             the index set J 
%%==============================================================================================
%clear all; clc;
rng('default');
HOME = pwd;
addpath(genpath(HOME));
%%
%profile on
%================================== Input =================================
Prob = 7;%[5:7];
flag_tol = 0;
tol = 1e-6; 
alpha_vec = [0.4:-0.1:0.1];%[0.9;0.5;0.1]; 
flag_J = 1;
%==========================================================================
if (flag_tol ~= 0) && (flag_tol ~= 1) && (flag_tol ~= 2)
    fprintf('The value of flag_tol must be 0, 1 or 2 !');
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
lenalp = length(alpha_vec);
lenprob = length(Prob);
result = zeros(lenalp*lenprob,11);
for pp = 1:lenprob
    i = Prob(pp);
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
    
    %% Input parameters
    OPTIONS.tol = tol;
    OPTIONS.m = m;
    OPTIONS.n = n;
    OPTIONS.maxiter = 200;
    OPTIONS.sigmascale = 1.3;
    OPTIONS.tauscale = 1.2;
    OPTIONS.flag_tol = flag_tol;
    lamb = 0.12;
    
    eval(['result_',num2str(m),'_',num2str(n),'_','J = zeros(',num2str(lenalp),...
        ',',num2str(OPTIONS.maxiter),');']);
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
        OPTIONS.sigma = 20/OPTIONS.lambda;
        OPTIONS.tau = 20/OPTIONS.lambda;
        
        [obj,~,~,~,runhist,info] = NALM(A,b,OPTIONS);
        
        result(jj+(pp-1)*lenalp,1) = OPTIONS.m;
        result(jj+(pp-1)*lenalp,2) = OPTIONS.n;
        result(jj+(pp-1)*lenalp,3) = OPTIONS.lambda;
        result(jj+(pp-1)*lenalp,4) = OPTIONS.kk;
        result(jj+(pp-1)*lenalp,5) = info.xnnz;
        result(jj+(pp-1)*lenalp,6) = info.eta_res;
        result(jj+(pp-1)*lenalp,7) = info.res_kkt_final;
        result(jj+(pp-1)*lenalp,8) = info.iter;
        result(jj+(pp-1)*lenalp,9) = info.numSSN;
        result(jj+(pp-1)*lenalp,10) = info.totletime;
        result(jj+(pp-1)*lenalp,11) = obj(1);
        if OPTIONS.flag_tol == 1
            result(jj+(pp-1)*lenalp,12) = info.relobj;
        end
        
        if flag_J == 1
            iter = info.iter; r_indexJ_vec = runhist.r_indexJ;
            eval(['result_',num2str(m),'_',num2str(n),'_J(',num2str(jj),',1:',...
                num2str(iter),') = runhist.r_indexJ;']);
        end
    end
    eval(['save result_',num2str(m),'_',num2str(n),'_J.mat result_',num2str(m),'_',num2str(n),'_J']);
end

eval(['save Result_NALM_random_flagtol_',num2str(flag_tol),'_tol_',num2str(tol),...
    '_flagJ_',num2str(flag_J),'.mat result']); 

%profile viewer





%%*********************************************************************