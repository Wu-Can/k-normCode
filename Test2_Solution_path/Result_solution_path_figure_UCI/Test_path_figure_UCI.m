%% =====================================================================================
%% Generate figures for the solution path with equally divided grid points on UCI data
%% Input:
%% flag_time = 1, generate the figure of "Time comparison"; 
%%           = 0, does not generate the figure of "Time comparison";
%% flag_nnz = 1, generate the figure of "Selected active features in the AS";
%%          = 0, does not generate the figure "Selected active features in the AS".
%% =====================================================================================
clc;clear
%========= INPUT =========
flag_time = 1;
flag_nnz = 1;
%=========================

Prob = [9, 4, 2];
lenp = length(Prob);
for i = 1:lenp
    p = Prob(i);
    if (p == 9)
        Probname = 'housing7';
        m = 506; n = 77520; esgp = 40; lam_ed = 6; ymax = 350; time_ed = 1;
        lambda_vec = [4.5:(1.5-4.5)/esgp:1.5;4.5:(1.5-4.5)/esgp:1.5;4.5:(1.5-4.5)/esgp:1.5];
    elseif (p == 4)
        Probname = 'loglp.E2006.test';
        m = 3308; n = 1771946; esgp = 90; lam_ed = 5; ymax = 300; time_ed = 60;
        lambda_vec = [45:(20-45)/esgp:20;75:(45-75)/esgp:45;85:(60-85)/esgp:60];
    elseif (p == 2)
        Probname = 'loglp.E2006.train';
        m = 16087; n = 4265669; esgp = 300; lam_ed = 4; ymax = 350; time_ed = 3600;
        lambda_vec = [140:-100/esgp:40;300:-200/esgp:100;330:-200/esgp:130];
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
        time_AS_NALM_new = result_time_path./time_ed;
        clear result_time_path;
        eval(['load new_Result_warm_NALM_',num2str(m),'_',num2str(n),'_',num2str(esgp)]);
        time_warm_SNIPAL_new = result_ttime_path./time_ed;
        clear result_ttime_path;
        eval(['load new_Result_NALM_path_time_',num2str(m),'_',num2str(n),'_',num2str(esgp)]);
        time_SNIPAL_new = result_ttime_path./time_ed;
        clear result_ttime_path;
    end
    
    [lenalp,lenlam] = size(lambda_vec);
    alpha_vec = [0.9,0.5,0.1];
    
   %% Time comparison
    if flag_time == 1
        figure(1)
        for jj = 1:lenalp
            alpha = alpha_vec(jj);
            subplot(lenalp,lenp,i+lenp*(jj-1))
            lamvec = lambda_vec(jj,1:lenlam);
            minlam = min(lamvec); maxlam = max(lamvec);
            plot(lamvec,time_AS_NALM_new(jj,1:lenlam),'-pr');
            hold on
            plot(lamvec,time_warm_SNIPAL_new(jj,1:lenlam),'-ob');
            hold on
            plot(lamvec,time_SNIPAL_new(jj,1:lenlam),'-*g');
            xlabel('\lambda');
            if p == 2
                ylabel('Time(h)');
            elseif p == 4
                ylabel('Time(m)');
            else
                ylabel('Time(s)');
            end
            
            xlim([minlam,maxlam]);
            title( {[ 'Time comparison on ', Probname];[ ' (\alpha=', num2str(alpha),...
                ', ','esgp=',num2str(esgp+1),')' ]}, 'FontSize',12);
            legend({'AS+N-ALM','Warm+N-ALM','N-ALM'},'Location','northwest','FontSize',9);
            set(gca,'XDir','reverse','XTick',[minlam:(maxlam-minlam)/lam_ed:maxlam]);
        end
    end
     %% Selected active features in the AS
    if flag_nnz == 1
       figure(2)
       for jj = 1:lenalp
           alpha = alpha_vec(jj);
           subplot(lenalp,lenp,i+lenp*(jj-1))
           lamvec = lambda_vec(jj,1:lenlam);
           minlam = min(lamvec); maxlam = max(lamvec);
           plot(lamvec,nnzx_new(jj,1:lenlam),'-*r');
           hold on
           plot(lamvec,mean_n_new(jj,1:lenlam),'-db');
           xlabel('\lambda','FontSize',12);
           ylabel('Feature size','FontSize',12);
           xlim([minlam,maxlam]);
           ylim([0,ymax]);
           title({['Selected active features'];[' on ', Probname, ' (\alpha=', num2str(alpha),...
               ', ','esgp=', num2str(esgp+1),')']},'FontSize',12);
           legend({'nnz(x)','mean(n)'},'Location','northwest','FontSize',10);
           set(gca,'XDir','reverse','XTick',[minlam:(maxlam-minlam)/lam_ed:maxlam]);
       end
    end
end

saveas(figure(1),'Time_comparison_9_4_2.fig');
saveas(figure(2),'Selected_active_features_9_4_2.fig');



