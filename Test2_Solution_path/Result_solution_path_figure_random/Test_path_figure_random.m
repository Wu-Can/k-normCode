%% ======================================================================================
%% Generate figures for the solution path with equally divided grid points on random data
%% Input:
%% flag_time = 1, generate the figure of "Time comparison"; 
%%           = 0, does not generate the figure of "Time comparison";
%% flag_nnz = 1, generate the figure of "Selected active features in the AS";
%%          = 0, does not generate the figure "Selected active features in the AS".
%% ======================================================================================
clc;clear
%========= INPUT =========
flag_time = 1;
flag_nnz = 1;
%=========================
Prob = [1 2 3];
lenp = length(Prob);
for i = 1:lenp
    p = Prob(i);
    if p == 1
        m = 500; n = 100000;
        esgp = 45;
    elseif p == 2
        m = 1000; n = 200000;
        esgp = 45;
    elseif p == 3
        m = 3000; n = 500000;
        esgp = 90;
    end
    
    if flag_nnz == 1
        eval(['load new_Result_AS_NALM_path_nnz_',num2str(m),'_',num2str(n),'_',num2str(esgp)]);
        nnzx_new = result_xnnz_path;
        clear result_xnnz_path;
        eval(['load new_Result_AS_NALM_path_n_mean_',num2str(m),'_',num2str(n),'_',num2str(esgp)]);
        mean_n_new =  result_n_mean;
        clear result_n_mean;
    end
    if flag_time == 1
        eval(['load new_Result_AS_NALM_path_time_',num2str(m),'_',num2str(n),'_',num2str(esgp)]);
        time_AS_NALM_new = result_time_path;
        clear result_time_path;
        eval(['load new_Result_warm_NALM_time_',num2str(m),'_',num2str(n),'_',num2str(esgp)]);
        time_warm_NALM_new = result_ttime_path;
        clear result_ttime_path;
        eval(['load new_Result_NALM_path_time_',num2str(m),'_',num2str(n),'_',num2str(esgp)]);
        time_NALM_new = result_ttime_path;
        clear result_ttime_path;
    end
    
    if p == 1
        lambda_vec = [20:(8-20)/esgp:8;60:(45-60)/esgp:45;70:(55-70)/esgp:55];
    elseif p == 2
        lambda_vec = [40:(25-40)/esgp:25;80:(65-80)/esgp:65;110:(95-110)/esgp:95];
    elseif p == 3
        lambda_vec = [80:-30/esgp:50;160:-30/esgp:130;200:-30/esgp:170];
        time_AS_NALM_new = time_AS_NALM_new./3600;
        time_warm_NALM_new = time_warm_NALM_new./3600;
        time_NALM_new = time_NALM_new./3600;
    end
    [lenalp,lenlam] = size(lambda_vec);
    alpha_vec = [0.9,0.5,0.1];
    
   %% Time comparison
    if flag_time == 1
        figure(1)
        for jj = 1:lenalp
            alpha = alpha_vec(jj);
            subplot(lenalp,lenp,i+lenp*(jj-1))
            lambda = lambda_vec(jj,1:lenlam);
            maxlam = max(lambda); minlam = min(lambda);
            if (p == 1) || (p == 2)
                if mod(maxlam-minlam,2) == 0
                    lam_ed = (maxlam-minlam)/2;
                else
                    lam_ed = (maxlam-minlam)/3;
                end
            elseif (p == 3)
                lam_ed = 6;
            end
            plot(lambda,time_AS_NALM_new(jj,1:lenlam),'-pr');
            hold on
            plot(lambda,time_warm_NALM_new(jj,1:lenlam),'-ob');
            hold on
            plot(lambda,time_NALM_new(jj,1:lenlam),'-*g');
            xlabel('\lambda');
            if p == 3
                ylabel('Time(h)');
            else
                ylabel('Time(s)');
            end
            
            xlim([minlam,maxlam]);
            title({['Time comparison '];['(n=',num2str(m),', d=',num2str(n),', \alpha=',...
                num2str(alpha),', ','esgp=',num2str(esgp+1),')']},'FontSize',12);
            legend({'AS+N-ALM','Warm+N-ALM','N-ALM'},'Location','northwest','FontSize',9);
            set(gca,'XDir','reverse','XTick',[minlam:(maxlam-minlam)/lam_ed:maxlam])
        end
    end
     %% Selected active features in the AS
    if flag_nnz == 1
       figure(2)
       for jj = 1:lenalp
           alpha = alpha_vec(jj);
           subplot(lenalp,lenp,i+lenp*(jj-1))
           lambda = lambda_vec(jj,1:lenlam);
           maxlam = max(lambda); minlam = min(lambda);
           if (p == 1) || (p == 2)
               if mod(maxlam-minlam,2) == 0
                   lam_ed = (maxlam-minlam)/2;
               else
                   lam_ed = (maxlam-minlam)/3;
               end
           else
               lam_ed = 6;
           end
           plot(lambda,nnzx_new(jj,1:lenlam),'-*r');
           hold on
           plot(lambda,mean_n_new(jj,1:lenlam),'-db');
           xlabel('\lambda','FontSize',12);
           ylabel('Feature size','FontSize',12);
           xlim([min(lambda),max(lambda)]);
           title({['Selected active features'];[ '(n=',num2str(m),', d=',num2str(n),...
               ', \alpha=',num2str(alpha),', ','esgp=',num2str(esgp+1),')']},'FontSize',12);
           legend({'nnz(x)','mean(n)'},'Location','northwest','FontSize',10);
           set(gca,'XDir','reverse','XTick',[minlam:(maxlam-minlam)/lam_ed:maxlam])
       end
    end
end

saveas(figure(1),'Time_comparison.fig');
saveas(figure(2),'Selected_active_features.fig');


