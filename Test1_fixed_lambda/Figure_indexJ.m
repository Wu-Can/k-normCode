%%=========================================================================
%% Generate the Figure 1 in the supplementary materials
%% Input:
%% iter_0 = the number of starting test iteration of the N-ALM 
%% iter_1 = the number of final test iteration of the N-ALM
%% m = the sample size
%% n = the feature size
%%=========================================================================
clc;clear
%============= INPUT ==========
iter_0 = 6;
iter_1 = 20;
m = 3000; 
n = 500000;
%==============================

load('Result_NALM_random_flagtol_0_tol_1e-06_flagJ_1.mat'); 
eval(['load result_',num2str(m),'_',num2str(n),'_J.mat']);
itermin = min(result(:,8));

eval(['lenalp=size(result_',num2str(m),'_',num2str(n),'_J,1);']); 
eval(['Test_Matrix_s = result_',num2str(m),'_',num2str(n),'_J(:,','1:itermin);']);
if (1 > iter_0) || (iter_0 >= iter_1) || (iter_1 > itermin)
    fprintf('Notice: iter_0 and iter_1 must satisfy 1 <= iter_0 < iter_1 <= %2.0d',...
        itermin);
    return;
end

figure
for ii = 1:lenalp
    s_vec = Test_Matrix_s(ii,iter_0:iter_1);
    switch ii
        case 1
            plot(iter_0:iter_1,s_vec,'-*b','Linewidth',1.5,'MarkerSize',8); hold on
        case 2
            plot(iter_0:iter_1,s_vec,'-dr','Linewidth',1.5,'MarkerSize',8); hold on
        case 3
            plot(iter_0:iter_1,s_vec,'-ok','Linewidth',1.5,'MarkerSize',8); hold on
        case 4
            plot(iter_0:iter_1,s_vec,'-pm','Linewidth',1.5,'MarkerSize',8);
    end
end
 
xlabel('iteration','FontSize',15);
ylabel('Average number of s','FontSize',15);

eval(['title(''The average number of s for the N-ALM on synthetic data with n = ',...
    num2str(m),' and d = ',num2str(n),''');']);
legend({'k = 60%n','k = 70%n','k = 80%n','k = 90%n'},'Location','northeast','FontSize',12);

set(gca, 'xTick', [iter_0:1:iter_1]);
set(gca, 'xTickLabel', [iter_0:1:iter_1]);
set(gca, 'FontSize', 15);

maxs = 2050;
set(gca, 'yTick', [0:250:maxs]);
set(gca, 'yTickLabel', [0:250:maxs]);
set(gca, 'FontSize', 15);
xlim([iter_0,iter_1]);
ylim([0,2050]);

saveas(gcf,'Figure_s_indexJ.png')

