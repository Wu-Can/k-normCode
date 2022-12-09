%%=========================================================================================
%% Test the MM+N-PPA for the nonconvex truncated CVaR-based model with UCI data
%%=========================================================================================
clear all; clc;
rng('default');
%%
%profile on
%================================== INPUT =================================
prob = [10 11 5 8 9];
%==========================================================================
HOME = pwd;
addpath(genpath(HOME));
filepath = fileparts(HOME);
datadir = [filepath,filesep,'UCIdata'];

fname{1} = 'E2006.train';
fname{2} = 'log1p.E2006.train';
fname{3} = 'E2006.test';
fname{4} = 'log1p.E2006.test';
fname{5} = 'pyrim_scale_expanded5';
fname{6}= 'triazines_scale_expanded4';
fname{7} = 'abalone_scale_expanded7';
fname{8} = 'bodyfat_scale_expanded7';
fname{9} = 'housing_scale_expanded7';
fname{10} = 'mpg_scale_expanded7';
fname{11} = 'space_ga_scale_expanded9';

lamcvec = [1e-5 1e-6];
lenprob = length(prob);
lenlamc = length(lamcvec);
result = zeros(lenprob*6,11);
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
    OPTIONS.UCI = 1;
    OPTIONS.rho0 = 15;
    OPTIONS.etascale = 1.2;
    alpha = [0.1,0.9];
    veckk1 = [ceil((1-alpha(1))*m),ceil((1-alpha(2))*m),ceil((1-alpha(1))*m)];
    veckk2 = [veckk1(1)-1,veckk1(2)-1,veckk1(2)];
    lenk = length(veckk1);
    for jj = 1:lenk
        kk1 = veckk1(jj);
        kk2 = veckk2(jj);
        if (abs(kk1-kk2) > 0.5*m) || (( m < 300) && (n/m) > 400)
            OPTIONS.rhoscale = 0.2;
            OPTIONS.sigmascale = 0.2;
            OPTIONS.eta0 = 10;
        else
            OPTIONS.rhoscale = 0.7;
            OPTIONS.sigmascale = 0.6;
            OPTIONS.eta0 = 1;
        end
        
        for qq = 1:lenlamc
            lamc = lamcvec(qq);
            lambda = (kk1-kk2)*lamc*max(abs(b'*A));
            OPTIONS.lamc = (kk1-kk2)*lamc;
            if abs(kk1-kk2) == 1
                OPTIONS.sigma0 = -10*log10(lamc)-45;
            else
                OPTIONS.sigma0 = -5*log10(lamc)-10;
            end
            if qq == 1
                OPTIONS.Warm_starting = 0;
            else
                OPTIONS.Warm_starting = 1;
            end
            
            [~,obj,info,runhist] = MM_NPPA(A,b,kk1,kk2,lambda,OPTIONS);
            
            result(qq+(jj-1)*lenlamc+(ii-1)*(lenlamc*lenk),1) = kk1;
            result(qq+(jj-1)*lenlamc+(ii-1)*(lenlamc*lenk),2) = kk2;
            result(qq+(jj-1)*lenlamc+(ii-1)*(lenlamc*lenk),3) = lamc;
            result(qq+(jj-1)*lenlamc+(ii-1)*(lenlamc*lenk),4) = info.nnzeros_x;
            result(qq+(jj-1)*lenlamc+(ii-1)*(lenlamc*lenk),5) = info.obj_gap;
            result(qq+(jj-1)*lenlamc+(ii-1)*(lenlamc*lenk),6) = info.iter;
            result(qq+(jj-1)*lenlamc+(ii-1)*(lenlamc*lenk),7) = info.iterPPA;
            result(qq+(jj-1)*lenlamc+(ii-1)*(lenlamc*lenk),8) = info.iterSSN;
            result(qq+(jj-1)*lenlamc+(ii-1)*(lenlamc*lenk),9) = info.time;
            result(qq+(jj-1)*lenlamc+(ii-1)*(lenlamc*lenk),10) = obj;
            if  OPTIONS.Warm_starting == 1
                result(qq+(jj-1)*lenlamc+(ii-1)*(lenlamc*lenk),11) = info.warm_time;
            end
            
        end
    end
end
%save result_MM_NPPA_UCI.mat result
%profile viewer





