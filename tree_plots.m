function [ LMdata_fin ] = tree_plots( LMdata, LM_list )
% This function will compare a variety of plant water status indicators for
% individual trees. Some comparisons will include:
    % SWP vs ACSI
    % SWP vs ISCI
    % ACSI vs ICSI (or ICSI vs ACSI)
    % CWSI vs ASCI
    % CWSI vs ISCI
    % and more
    
%% Import SWP data and Initialization for Plotting Features
year = LMdata.LM1.TimeStamps(1,1);
if year==2017
    SWP = xlsread('SWPdata_2017_amended.xls','MPa');
% elseif year==2016
%     
end

k=1; % counter for subplot number

%parameters for figure and panel size
plotheight=19.6;
plotwidth=16;
subplotsx=3;
subplotsy=3;   
leftedge=1.2;
rightedge=0.4;   
topedge=1;
bottomedge=1.5;
spacex=0.2;
spacey=0.6;
fontsize=5;  
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
%setting the Matlab figure      

%% Plotting 'For Loop'
counter = 1;
for g=[4,5,6,7,8,10,12,13,15]
    
    if k==1 || k==4 || k==7 
        i=1;
    elseif k==2 || k==5 || k==8 
        i=2;
    elseif k==3 || k==6 || k==9 
        i=3;
    end
    
    if k>=1 && k<=3
        ii=1;
    elseif k>=4 && k<=6
        ii=2;
    elseif k>=7 && k<=9
        ii=3;
    end

    if LM_list(g)==35 && year==2017 % skip Tree 35 for 2017
        continue
    end
   
    LMname = ['LM', num2str(LM_list(g))]; 
    % create a string to name the leaf monitor

    not_nan = isnan(SWP(:,g+1))==0; % get row numbers where SWP is not NaN
    LMdata.(LMname).SWP = SWP(not_nan,g+1); % trees are in numeric order in the 
    % Excel spreadsheet; this can be verified by opening the Excel 
    % spreadsheet. It's g+1 because we want to ignore the first column of
    % SWP, because it is the Julian Dates (we will use this column for
    % something else, just not this line of the code).
    
    LMdata.(LMname).SWPdays = SWP(not_nan,1); % the Julian dates that correspond 
    % to the SWP measurements
    
% We need to find the row & column #s of the stress indices that correspond 
    % to the days when we took SWP measurements (since we did not take SWP
    % measurements every day that we obtained data from the leaf monitors).
    n = ismember(LMdata.(LMname).days, LMdata.(LMname).SWPdays);
    
    % Create a table with results from the days SWP was measured
    T = table(LMdata.(LMname).SWPdays, LMdata.(LMname).SWP,...
        LMdata.(LMname).days(n), LMdata.(LMname).ratio_int(n),...
        LMdata.(LMname).ratio_daily_avg(n), LMdata.(LMname).CWSI(n));
    T.Properties.VariableNames = {'SWP_Days','SWP','LM_Days','ICSI',...
                                  'ACSI','CWSI'};
    disp(LMname)
    T    
    
    S = table(LMdata.(LMname).days, LMdata.(LMname).ratio_int,...
        LMdata.(LMname).ratio_daily_avg, LMdata.(LMname).CWSI);
    S.Properties.VariableNames = {'LM_Days','ICSI','ACSI','CWSI'};
    disp(LMname)
    S
    
%% SWP vs (Ta-Tl)average between 12 pm and 4 pm    
%      for f=1:length(LMdata.(LMname).days)
%          
%         p = LMdata.(LMname).JD==LMdata.(LMname).days(f) ...
%             & LMdata.(LMname).time>=12 & ...
%             LMdata.(LMname).time<=16; % 12 pm - 4 pm
%     
%         % Air Temperature
%         Tair_p_avg(f) = mean(LMdata.(LMname).Ta_fil(p));
%         
%         % Live Leaf Temperature
%         Tleaf_p_avg(f) = mean(LMdata.(LMname).Tl_fil(p));
%     
%         % Dry Leaf Temperature
%         Tdry_p_avg(f) = mean(LMdata.(LMname).Tdry_fil(p));
%                    
%     end
%     
%       % Temperature Differentials
%       dTal = Tair_p_avg - Tleaf_p_avg; dTal = dTal'; % Tair - Tleaf
%       dTdl = Tdry_p_avg - Tleaf_p_avg; dTdl = dTdl'; % Tdry - Tleaf
%       
%       x = dTal(n);
%       y = LMdata.(LMname).SWP;
%       
%        not_nan = isnan(x)==0; % get row numbers where it is not NaN
%         x = x(not_nan);
%         y = y(not_nan); 
%         
%         tbl = table(x,y);
%         mdl = fitlm(tbl,'linear');
% 
%         figure(299);
%         set(gcf, 'PaperUnits', 'centimeters');
%         set(gcf, 'PaperSize', [plotwidth plotheight]);
%         set(gcf, 'PaperPositionMode', 'manual');
%         set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
%         set(gcf, 'Position', get(0, 'Screensize'));
%         ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
%         plot(mdl)
%         ylim([0,3])
%         xlim([0,5])
%         xlabel('')
%         ylabel('')
%         title(['SWP vs Tair-Tleaf: ', LMname, ', R²=',num2str(mdl.Rsquared.Ordinary)]);
%        
%         if ii>1
%         set(ax,'xticklabel',[])
%         end
%         
%         if i>1
%         set(ax,'yticklabel',[])
%         end
% 
%         if i==1
%         ylabel('SWP (MPa)')
%         end
% 
%         if ii==1
%         xlabel('Tair-Tleaf (°C)')
%         end
                              
%% SWP vs ICSI
    x = LMdata.(LMname).ratio_int(n);
    y = LMdata.(LMname).SWP;    
    not_nan = isnan(x)==0; % get row numbers where it is not NaN
    x = x(not_nan);
    y = y(not_nan);    
    tbl = table(x,y);
    mdl = fitlm(tbl,'linear');

    figure(300);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [plotwidth plotheight]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
    set(gcf, 'Position', get(0, 'Screensize'));
    ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
    plot(mdl)
    ylim([0,3])
    xlim([0,30])
    xlabel('')
        ylabel('')
%     title(['SWP vs ICSI: ', LMname, ', R^2=',num2str(mdl.Rsquared.Adjusted)]);
    title(['SWP vs ICSI: ', LMname, ', R²=',num2str(mdl.Rsquared.Ordinary)]);

    if ii>1
    set(ax,'xticklabel',[])
    end

    if i>1
    set(ax,'yticklabel',[])
    end

    if i==1
    ylabel('SWP (MPa)')
    end

    if ii==1
    xlabel('ICSI (°C/kPa)')
    end
    
    set(gca,'fontsize',10)
    
    table_ICSI_SWP(counter,:) = table({LMname},mdl.Rsquared.Ordinary,...
        mdl.Rsquared.Adjusted,mdl.MSE,mdl.RMSE,mdl.Coefficients(1,1),...
        mdl.Coefficients(2,1),...
        mdl.Coefficients(1,4),mdl.Coefficients(2,4));
    table_ICSI_SWP.Properties.VariableNames = {'LM','RsquaredOrdinary',...
        'RsquaredAdjusted','MSE','RMSE','InterceptEstimate','xEstimate',...
        'pValueIntercept','pValuex'};

%% SWP vs ACSI
    figure(301);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [plotwidth plotheight]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
    set(gcf, 'Position', get(0, 'Screensize'));
    x = LMdata.(LMname).ratio_daily_avg(n);
    y = LMdata.(LMname).SWP;
    not_nan = isnan(x)==0; % get row numbers where it is not NaN
    x = x(not_nan);
    y = y(not_nan);  
    tbl = table(x,y);
    mdl = fitlm(tbl,'linear');
    ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
    plot(mdl)
    ylim([0,3])
    xlim([0,4])
    xlabel('')
        ylabel('')
    title(['SWP vs ACSI: ', LMname, ', R²=',num2str(mdl.Rsquared.Ordinary)]) 

    if ii>1
    set(ax,'xticklabel',[])
    end

    if i>1
    set(ax,'yticklabel',[])
    end

    if i==1
    ylabel('SWP (MPa)')
    end

    if ii==1
    xlabel('ACSI (°C/kPa)')
    end
    
    set(gca,'fontsize',10)
    
    table_ACSI_SWP(counter,:) = table({LMname},mdl.Rsquared.Ordinary,...
        mdl.Rsquared.Adjusted,mdl.MSE,mdl.RMSE,mdl.Coefficients(1,1),...
        mdl.Coefficients(2,1),...
        mdl.Coefficients(1,4),mdl.Coefficients(2,4));
    table_ACSI_SWP.Properties.VariableNames = {'LM','RsquaredOrdinary',...
        'RsquaredAdjusted','MSE','RMSE','InterceptEstimate','xEstimate',...
        'pValueIntercept','pValuex'};
    
%% ICSI vs ACSI
   figure(302);
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperSize', [plotwidth plotheight]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf, 'Position', get(0, 'Screensize'));
    x = LMdata.(LMname).ratio_daily_avg(n);
    y = LMdata.(LMname).ratio_int(n);
    not_nan = isnan(x)==0; % get row numbers where it is not NaN
    x = x(not_nan);
    y = y(not_nan);  
    tbl = table(x,y);
    mdl = fitlm(tbl,'linear');    
    ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
    plot(mdl)
    ylim([0,30])
    xlim([0,4])
    xlabel('')
        ylabel('')
    title(['ICSI vs ACSI: ', LMname, ', R²=',num2str(mdl.Rsquared.Ordinary)])   
        if ii>1
        set(ax,'xticklabel',[])
        end

        if i>1
        set(ax,'yticklabel',[])
        end

        if i==1
        ylabel('ICSI (°C/kPa)')
        end

        if ii==1
        xlabel('ACSI (°C/kPa)')
        end
        
        set(gca,'fontsize',10)
        
        table_ICSI_ACSI(counter,:) = table({LMname},mdl.Rsquared.Ordinary,...
        mdl.Rsquared.Adjusted,mdl.MSE,mdl.RMSE,mdl.Coefficients(1,1),...
        mdl.Coefficients(2,1),...
        mdl.Coefficients(1,4),mdl.Coefficients(2,4));
    table_ACSI_SWP.Properties.VariableNames = {'LM','RsquaredOrdinary',...
        'RsquaredAdjusted','MSE','RMSE','InterceptEstimate','xEstimate',...
        'pValueIntercept','pValuex'};
    
%% CWSI vs ICSI
    figure(303);
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperSize', [plotwidth plotheight]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
        set(gcf, 'Position', get(0, 'Screensize'));
    x = LMdata.(LMname).ratio_int;
    y = LMdata.(LMname).CWSI;
    not_nan = isnan(x)==0; % get row numbers where it is not NaN
    x = x(not_nan);
    y = y(not_nan);  
    tbl = table(x,y);
    mdl = fitlm(tbl,'linear');
    ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
    plot(mdl)
    ylim([0,1])
    xlim([0,30])
    xlabel('')
        ylabel('')
    title(['CWSI vs ICSI: ', LMname, ', R²=',num2str(mdl.Rsquared.Ordinary)])
        if ii>1
        set(ax,'xticklabel',[])
        end

        if i>1
        set(ax,'yticklabel',[])
        end

        if i==1
        ylabel('CWSI')
        end

        if ii==1
        xlabel('ICSI (°C/kPa)')
        end
        
        set(gca,'fontsize',10)
        
        table_CWSI_ICSI(counter,:) = table({LMname},mdl.Rsquared.Ordinary,...
        mdl.Rsquared.Adjusted,mdl.MSE,mdl.RMSE,mdl.Coefficients(1,1),...
        mdl.Coefficients(2,1),...
        mdl.Coefficients(1,4),mdl.Coefficients(2,4));
    table_ACSI_SWP.Properties.VariableNames = {'LM','RsquaredOrdinary',...
        'RsquaredAdjusted','MSE','RMSE','InterceptEstimate','xEstimate',...
        'pValueIntercept','pValuex'};
    
%% CWSI vs ACSI
    figure(304);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [plotwidth plotheight]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
    set(gcf, 'Position', get(0, 'Screensize'));
    x = LMdata.(LMname).ratio_daily_avg;
    y = LMdata.(LMname).CWSI;
    not_nan = isnan(x)==0; % get row numbers where it is not NaN
    x = x(not_nan);
    y = y(not_nan);  
    tbl = table(x,y);
    mdl = fitlm(tbl,'linear');
    ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
    plot(mdl)
    ylim([0,1])
    xlim([0,4])
    xlabel('')
        ylabel('')
    title(['CWSI vs ACSI: ', LMname, ', R²=',num2str(mdl.Rsquared.Ordinary)]) 
        if ii>1
        set(ax,'xticklabel',[])
        end

        if i>1
        set(ax,'yticklabel',[])
        end

        if i==1
        ylabel('CWSI')
        end

        if ii==1
        xlabel('ACSI (°C/kPa)')
        end
        
        set(gca,'fontsize',10)
        
        table_CWSI_ACSI(counter,:) = table({LMname},mdl.Rsquared.Ordinary,...
        mdl.Rsquared.Adjusted,mdl.MSE,mdl.RMSE,mdl.Coefficients(1,1),...
        mdl.Coefficients(2,1),...
        mdl.Coefficients(1,4),mdl.Coefficients(2,4));
    table_ACSI_SWP.Properties.VariableNames = {'LM','RsquaredOrdinary',...
        'RsquaredAdjusted','MSE','RMSE','InterceptEstimate','xEstimate',...
        'pValueIntercept','pValuex'};
    
%% SWP vs CWSI
        figure(305);
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperSize', [plotwidth plotheight]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
        set(gcf, 'Position', get(0, 'Screensize'));
        x = LMdata.(LMname).CWSI(n);
        y = LMdata.(LMname).SWP;
        not_nan = isnan(x)==0; % get row numbers where it is not NaN
        x = x(not_nan);
        y = y(not_nan);  
        tbl = table(x,y);
        mdl = fitlm(tbl,'linear');
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        plot(mdl)
        ylim([0,3])
        xlim([0,1])
        xlabel('')
        ylabel('')
        title(['SWP vs CWSI: ', LMname, ', R²=',num2str(mdl.Rsquared.Ordinary)])  
        if ii>1
        set(ax,'xticklabel',[])
        end

        if i>1
        set(ax,'yticklabel',[])
        end

        if i==1
        ylabel('SWP (MPa)')
        end

        if ii==1
        xlabel('CWSI (°C/kPa)')
        end
        
        set(gca,'fontsize',10)
    
        table_SWP_CWSI(counter,:) = table({LMname},mdl.Rsquared.Ordinary,...
        mdl.Rsquared.Adjusted,mdl.MSE,mdl.RMSE,mdl.Coefficients(1,1),...
        mdl.Coefficients(2,1),...
        mdl.Coefficients(1,4),mdl.Coefficients(2,4));
    table_ACSI_SWP.Properties.VariableNames = {'LM','RsquaredOrdinary',...
        'RsquaredAdjusted','MSE','RMSE','InterceptEstimate','xEstimate',...
        'pValueIntercept','pValuex'};
        
%% SWP vs IDANS
        figure(306);
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperSize', [plotwidth plotheight]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
        set(gcf, 'Position', get(0, 'Screensize'));
        x = LMdata.(LMname).IDANS(n);
        y = LMdata.(LMname).SWP;
        not_nan = isnan(x)==0; % get row numbers where it is not NaN
        x = x(not_nan);
        y = y(not_nan);  
        tbl = table(x,y);
        mdl = fitlm(tbl,'linear');
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        plot(mdl)
        ylim([0,3])
        xlim([0,225])
        xlabel('')
        ylabel('')
        title(['SWP vs IDANS: ', LMname, ', R²=',num2str(mdl.Rsquared.Ordinary)]) 
        if ii>1
        set(ax,'xticklabel',[])
        end

        if i>1
        set(ax,'yticklabel',[])
        end

        if i==1
        ylabel('SWP (MPa)')
        end

        if ii==1
        xlabel('IDANS (°C*hr)')
        end
        
        set(gca,'fontsize',10)
        
        table_SWP_IDANS(counter,:) = table({LMname},mdl.Rsquared.Ordinary,...
        mdl.Rsquared.Adjusted,mdl.MSE,mdl.RMSE,mdl.Coefficients(1,1),...
        mdl.Coefficients(2,1),...
        mdl.Coefficients(1,4),mdl.Coefficients(2,4));
    table_ACSI_SWP.Properties.VariableNames = {'LM','RsquaredOrdinary',...
        'RsquaredAdjusted','MSE','RMSE','InterceptEstimate','xEstimate',...
        'pValueIntercept','pValuex'};
        
%% ICSI vs IDANS
        figure(307);
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperSize', [plotwidth plotheight]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
        set(gcf, 'Position', get(0, 'Screensize'));
        x = LMdata.(LMname).IDANS;
        y = LMdata.(LMname).ratio_int;
        not_nan = isnan(x)==0; % get row numbers where it is not NaN
        x = x(not_nan);
        y = y(not_nan);  
        tbl = table(x,y);
        mdl = fitlm(tbl,'linear');
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        plot(mdl)
        ylim([0,30])
        xlim([0,225])
        xlabel('')
        ylabel('')
        title(['ICSI vs IDANS: ', LMname, ', R²=',num2str(mdl.Rsquared.Ordinary)])  
        if ii>1
        set(ax,'xticklabel',[])
        end

        if i>1
        set(ax,'yticklabel',[])
        end

        if i==1
        ylabel('ICSI (°C/kPa)')
        end

        if ii==1
        xlabel('IDANS (°C*hr)')
        end
        
        set(gca,'fontsize',10)
        
        table_ICSI_IDANS(counter,:) = table({LMname},mdl.Rsquared.Ordinary,...
        mdl.Rsquared.Adjusted,mdl.MSE,mdl.RMSE,mdl.Coefficients(1,1),...
        mdl.Coefficients(2,1),...
        mdl.Coefficients(1,4),mdl.Coefficients(2,4));
    table_ACSI_SWP.Properties.VariableNames = {'LM','RsquaredOrdinary',...
        'RsquaredAdjusted','MSE','RMSE','InterceptEstimate','xEstimate',...
        'pValueIntercept','pValuex'};
    
%% ACSI vs IDANS
        figure(308);
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperSize', [plotwidth plotheight]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
        set(gcf, 'Position', get(0, 'Screensize'));
        x = LMdata.(LMname).IDANS;
        y = LMdata.(LMname).ratio_daily_avg;
        not_nan = isnan(x)==0; % get row numbers where it is not NaN
        x = x(not_nan);
        y = y(not_nan);  
        tbl = table(x,y);
        mdl = fitlm(tbl,'linear');
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        plot(mdl)
        ylim([0,4])
        xlim([0,225])
        xlabel('')
        ylabel('')
        title(['ACSI vs IDANS: ', LMname, ', R²=',num2str(mdl.Rsquared.Ordinary)])  
        if ii>1
        set(ax,'xticklabel',[])
        end

        if i>1
        set(ax,'yticklabel',[])
        end

        if i==1
        ylabel('ACSI (°C/kPa)')
        end

        if ii==1
        xlabel('IDANS (°C*hr)')
        end
        
        set(gca,'fontsize',10)
        
        table_ACSI_IDANS(counter,:) = table({LMname},mdl.Rsquared.Ordinary,...
        mdl.Rsquared.Adjusted,mdl.MSE,mdl.RMSE,mdl.Coefficients(1,1),...
        mdl.Coefficients(2,1),...
        mdl.Coefficients(1,4),mdl.Coefficients(2,4));
    table_ACSI_SWP.Properties.VariableNames = {'LM','RsquaredOrdinary',...
        'RsquaredAdjusted','MSE','RMSE','InterceptEstimate','xEstimate',...
        'pValueIntercept','pValuex'};
    
%% CWSI vs IDANS
        figure(309);
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperSize', [plotwidth plotheight]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
        set(gcf, 'Position', get(0, 'Screensize'));
        x = LMdata.(LMname).IDANS;
        y = LMdata.(LMname).CWSI;
        not_nan = isnan(x)==0; % get row numbers where it is not NaN
        x = x(not_nan);
        y = y(not_nan);  
        tbl = table(x,y);
        mdl = fitlm(tbl,'linear');
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        plot(mdl)
        ylim([0,1])
        xlim([0,225])
        xlabel('')
        ylabel('')
        title(['CWSI vs IDANS: ', LMname, ', R²=',num2str(mdl.Rsquared.Ordinary)])
        if ii>1
        set(ax,'xticklabel',[])
        end

        if i>1
        set(ax,'yticklabel',[])
        end

        if i==1
        ylabel('CWSI')
        end

        if ii==1
        xlabel('IDANS (°C*hr)')
        end
        
        set(gca,'fontsize',10)
        
        table_CWSI_IDANS(counter,:) = table({LMname},mdl.Rsquared.Ordinary,...
        mdl.Rsquared.Adjusted,mdl.MSE,mdl.RMSE,mdl.Coefficients(1,1),...
        mdl.Coefficients(2,1),...
        mdl.Coefficients(1,4),mdl.Coefficients(2,4));
    table_ACSI_SWP.Properties.VariableNames = {'LM','RsquaredOrdinary',...
        'RsquaredAdjusted','MSE','RMSE','InterceptEstimate','xEstimate',...
        'pValueIntercept','pValuex'};
    
    k=k+1;
    counter = counter + 1;
end

disp(table_ICSI_SWP);
disp(table_ACSI_SWP);
disp(table_ICSI_ACSI);
disp(table_CWSI_ICSI);
disp(table_CWSI_ACSI);
disp(table_SWP_CWSI);
disp(table_SWP_IDANS);
disp(table_ICSI_IDANS);
disp(table_ACSI_IDANS);
disp(table_CWSI_IDANS);

% figure(300)
% print('fig0','-dpng')
% 
% figure(301)
% print('fig1','-dpng')
% 
% figure(302)
% print('fig2','-dpng')
% 
% figure(303)
% print('fig3','-dpng')
% 
% figure(304)
% print('fig4','-dpng')
% 
% figure(305)
% print('fig5','-dpng')
% 
% figure(306)
% print('fig6','-dpng')
% 
% figure(307)
% print('fig7','-dpng')
% 
% figure(308)
% print('fig8','-dpng')
% 
% figure(309)
% print('fig9','-dpng')

LMdata_fin = LMdata;

end

