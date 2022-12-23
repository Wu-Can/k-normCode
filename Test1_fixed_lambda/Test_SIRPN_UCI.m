%%=========================================================================
%% Test the S-IRPN for the convex CVaR-based models with UCI data
%%=========================================================================
clear all; clc;
rng('default');
HOME = pwd;
addpath(genpath(HOME));
filepath = fileparts(HOME);
datadir = [filepath,filesep,'UCIdata']; 
%%
%profile on
%=================== INPUT ===================
prob = [1:11];
%=============================================

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

load Result_Gorubi_UCI_1_11_1e-9.mat
obj_opt_vec = result(:,3);
clear result;
%%
alpha_vec = [0.9 0.5 0.1];
lenalp = length(alpha_vec);
numProb =  length(prob);
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
    [n,d] = size(A);
    
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

    OPTIONS.tol = 1e-3;
    OPTIONS.n = n;
    OPTIONS.d = d;
    OPTIONS.lambda = lamc*max(abs(b'*A));
    OPTIONS.epsilon = 1;
    OPTIONS.tolsub = 1e-2;
    OPTIONS.mu = 1e-5;
    OPTIONS.maxiter = 9;
    OPTIONS.maxitersub = 200;
    OPTIONS.maxtime = 7200;
    OPTIONS.UCI = 1;
     
    for jj = 1:lenalp
        alpha = alpha_vec(jj);
        if alpha == 0.9
           j = 1;
        elseif alpha == 0.5
            j = 2;
        else
            j = 3;
        end
        OPTIONS.kk = ceil((1-alpha)*n);
        OPTIONS.obj_opt = obj_opt_vec(j+(i-1)*3);
                
        [obj,x,t,runhist,info] = SIRPN_CD(A,b,OPTIONS);

        result(jj+(ii-1)*lenalp,1) = OPTIONS.lambda;
        result(jj+(ii-1)*lenalp,2) = OPTIONS.kk;
        result(jj+(ii-1)*lenalp,3) = info.xnnz;
        result(jj+(ii-1)*lenalp,4) = info.relobj;
        result(jj+(ii-1)*lenalp,5) = info.iter;
        result(jj+(ii-1)*lenalp,6) = info.itersub;
        result(jj+(ii-1)*lenalp,7) = info.numiterCD;
        result(jj+(ii-1)*lenalp,8) = info.time;
        result(jj+(ii-1)*lenalp,9) = obj;
    end
end

save Result_SIRPN_UCI.mat result
%profile viewer
