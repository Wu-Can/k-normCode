%%==========================================================================================
%% Test the MM+Gurobi for the nonconvex truncated CVaR-based model with UCI data
%%==========================================================================================
clear all; clc;
rng('default');
%%
%profile on
%========================== INPUT =========================
prob = [10 11 5 8 9];
%==========================================================
HOME = pwd;
addpath(genpath(HOME));
filepath = fileparts(HOME);
datadir = [filepath,filesep,'UCIdata']; 

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

lamcvec = [1e-5 1e-6];
lenprob = length(prob);
lenlam = length(lamcvec);
result = zeros(lenprob*6,10);
for ii = 1:lenprob
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

    OPTIONS.tol = 1e-6;
    OPTIONS.m = m;
    OPTIONS.n = n;
    OPTIONS.rho0 = 15;
    OPTIONS.rho_iter = 4;
    
    alpha = [0.1,0.9];
    veckk1 = [ceil((1-alpha(1))*m),ceil((1-alpha(2))*m),ceil((1-alpha(1))*m)];
    veckk2 = [veckk1(1)-1,veckk1(2)-1,veckk1(2)];
    lenveckk = length(veckk1);
    for jj =  1:lenveckk
        kk1 = veckk1(jj);
        kk2 = veckk2(jj);
        if (abs(kk1-kk2) > 0.5*m) || (( m < 300) && (n/m) > 400)
            OPTIONS.rhoscale = 0.2;
            OPTIONS.sigmascale = 0.2;
        else
            OPTIONS.rhoscale = 0.7;
            OPTIONS.sigmascale = 0.6;
        end
        for qq = 1:lenlam
            lamc = lamcvec(qq);
            lambda = (kk1-kk2)*lamc*max(abs(b'*A));
            OPTIONS.lamc = (kk1-kk2)*lamc;
            if qq == 1
                OPTIONS.Warm_starting = 0;
            else
                OPTIONS.Warm_starting = 1;
            end
            if abs(kk1-kk2) == 1
                OPTIONS.sigma0 = -10*log10(lamc)-45;
            else
                OPTIONS.sigma0 = -5*log10(lamc)-10;
            end
            
            [~,obj,info,runhist] = MM_Gurobi(A,b,kk1,kk2,lambda,OPTIONS);
            
            result(qq+(jj-1)*lenlam+(ii-1)*lenveckk*lenlam,1) = OPTIONS.m;
            result(qq+(jj-1)*lenlam+(ii-1)*lenveckk*lenlam,2) = OPTIONS.n;            
            result(qq+(jj-1)*lenlam+(ii-1)*lenveckk*lenlam,3) = kk1;
            result(qq+(jj-1)*lenlam+(ii-1)*lenveckk*lenlam,4) = kk2;
            result(qq+(jj-1)*lenlam+(ii-1)*lenveckk*lenlam,5) = lamc;
            result(qq+(jj-1)*lenlam+(ii-1)*lenveckk*lenlam,6) = info.nnzeros_x;
            result(qq+(jj-1)*lenlam+(ii-1)*lenveckk*lenlam,7) = info.obj_gap;
            result(qq+(jj-1)*lenlam+(ii-1)*lenveckk*lenlam,8) = info.iter;
            result(qq+(jj-1)*lenlam+(ii-1)*lenveckk*lenlam,9) = info.iterBar;
            result(qq+(jj-1)*lenlam+(ii-1)*lenveckk*lenlam,10) = info.time;
            result(qq+(jj-1)*lenlam+(ii-1)*lenveckk*lenlam,11) = obj;
            if OPTIONS.Warm_starting == 1
                result(qq+(jj-1)*lenlam+(ii-1)*lenveckk*lenlam,12) = info.warm_time;
            end
        end
    end
end
 

save result_MM_Gurobi.mat result

%profile viewer


 

