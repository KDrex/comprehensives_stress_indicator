 %% Leaf Ratio Analysis Program
% Kelley Drechsler

clear
close all
clc

% Get files path to access program files
fileID = fopen('filesPath.txt','r');
A = fscanf(fileID,'%s'); % files path string
fclose(fileID);

% Get file path to access LM data
LMpath = get_LM_filepath( ); 

% Import LM Data
rawdata = importdata(LMpath); 

% Build LM Structure    
[ LMdata, LM_list ] = build_tree_struct( rawdata.data ); 

% Process LM Data
[ LMdata ] = LMdata_processing( LMdata, LM_list );

% Comprehensive Stress Ratio (CSR)
graphing_period = 196:202;
LMdata = Stom_Cond_Ratio( LMdata, graphing_period, LM_list, 'n' );

% Daily Ratio Integral & Daily Ratio Average
LMdata = Stom_Cond_Integral_and_Avg( LMdata, LM_list );

% % Write the ratio results to Excel
% write2Excel( LMdata, LM_list )

% CWSI
LMdata = CWSI(LMdata, LM_list);

% DANS & IDANS
graphing_period = 196:202;
LMdata = IDANS(LMdata, graphing_period, LM_list, 'n'); 

% Tree Plots
LMdata = tree_plots( LMdata, LM_list );

close Figure 1 % not sure where this blank figure is coming from

% % Treatment Zone Plots
% [TreatZone, LMdata] = treatment_zone_plots( LMdata, LM_list );

% Normalized Comprehensive Stress Ratio (NCSR)
% LMdata = NCSR( LMdata, LM_list );
% 
% % Least-Squares Linear Regression: Tair-Tleaf vs. VPD (Saturated Leaf Only)
% saturated_leaf_model( LMdata )

disp('Done!')

%% SWP vs Ta-Tl









%% verifying steady state assumption is okay (dT/dt=0 for dt=15 minutes)
dT=cell(1,length(LM_list));
max_dT=cell(1,length(LM_list));
mean_dT=cell(1,length(LM_list));
for f=1:length(LM_list)
    LMname = ['LM', num2str(LM_list(f))];
    n = LMdata.(LMname).TimeStamps(:,4)>=10 & LMdata.(LMname).TimeStamps(:,4)<=18; 
    % data from 10 AM to 6 PM
    Tl_n = LMdata.(LMname).Tl(n);
    dT{f} = Tl_n(2:end) - Tl_n(1:end-1);
    max_dT{f} = max(dT{f});
    mean_dT{f} = mean(dT{f});
    median_dT{f} = median(dT{f});
end
disp('program ran successfully')



