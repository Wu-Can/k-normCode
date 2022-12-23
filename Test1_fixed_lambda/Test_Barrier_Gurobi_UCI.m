%%==============================================================================================
%% Test the barrier method in Gurobi for the convex CVaR-based models with UCI data
%%==============================================================================================
clear all; clc;
rng('default');
%%
%profile on
%======== INPUT ========
prob = [1:11];
%=======================

%diary 'Diary_Barrier_Gurobi_UCI.txt'
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
lenprob = length(prob);
result = zeros(lenprob*2,8);
for ii = 1:lenprob
    i = prob(ii);
    probname = [datadir,filesep,fname{i}];
    fprintf('\n Problem name: %s', fname{i});
    if exist([probname,'.mat'])
        load([probname,'.mat'])
    else
        fprintf('\n Can not find the file in UCIdata');
        fprintf('\n ');
        return
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

    OPTIONS.m = m;
    OPTIONS.n = n;
    lamc_maxbTA = lamc*max(abs(b'*A));
    
    alpha = [0.9,0.5,0.1];
    lenalp = length(alpha);
    for jj = 1:lenalp
        OPTIONS.kk = ceil((1-alpha(jj))*m);
        OPTIONS.lambda = OPTIONS.kk*lamc_maxbTA;

        [sol, fval_bar,exitflag_bar,output_bar,~] = Gurobi_knorm_lp(A,b,OPTIONS);
        [sortx,~] = sort(abs(sol.x),'descend');
        normx1 = 0.999*norm(sol.x,1);
        tmpidex = find(cumsum(sortx) > normx1);
        if isempty(tmpidex)
            nnzeros = 0;
        else
            nnzeros = tmpidex(1);
        end
        
        if output_bar.time > 7200
            fval_bar = -1;
            output_bar.constrviolation = 0;
        end
        
        result(jj+(ii-1)*lenalp,1) = lamc;
        result(jj+(ii-1)*lenalp,2) = OPTIONS.kk;
        result(jj+(ii-1)*lenalp,3) = fval_bar;
        result(jj+(ii-1)*lenalp,4) = output_bar.baritercount;
        result(jj+(ii-1)*lenalp,5) = output_bar.time;
        result(jj+(ii-1)*lenalp,6) = exitflag_bar;
        result(jj+(ii-1)*lenalp,7) = nnzeros;
        result(jj+(ii-1)*lenalp,8) = output_bar.constrviolation;
    end
end

save Result_Gorubi_UCI3.mat result
%profile viewer
%diary off
