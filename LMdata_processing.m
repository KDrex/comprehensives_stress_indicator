function [ LMdata_processed ] = LMdata_processing( LMdata, LM_list )
% This function processes the leaf monitor data by:
% - Calibrating all RH and Tair sensors to the dry leaf monitor RH and Tair
% sensors (2016)
% - Adjusting for the Tleaf and RH offsets (both years)
% - smoothes the data with a moving average function



for g=1:length(LM_list)

    LMname = ['LM', num2str(LM_list(g))]; % create a string to name the leaf monitor

    % Calibrate the RH and Tair data for each LM to the Dry LM
    if LMdata.LM1.TimeStamps(1,1)==2016 % if it's 2016 (used LM1 as a sneaky 
        % way to access the year. Could have used any LM to get the year).

        if LM_list(g)~=40 % it's not meaningful to calibrate LM 40 against itself   
            LMdata.(LMname) = t_rh_cal(LMdata.(LMname), LMdata.LM40,...
                                           'n', 198, 199, g);     
        end

    end

% Offset Tleaf
graph='n';
if graph=='y'
    figure(200)
    subplot(4,4,g)
end
[ LMdata.(LMname).Tl, ~ ] = Tleaf_offset_daily( LMdata.(LMname).TimeStamps(:,4),...
                                    LMdata.(LMname).Ta,...
                                    LMdata.(LMname).Tl, ...
                                    LMdata.(LMname).JDD, ...
                                    LMdata.(LMname).JD,...
                                    LMdata.(LMname).sNum(1), graph,'Live Leaf IRT (\circC)' );    
hold on
graph='n';
if graph=='y'
    figure(201)
    subplot(4,4,g) 
end
 % Offset Tdry                              
 [ LMdata.(LMname).Tdry, ~ ] = Tleaf_offset_daily( LMdata.(LMname).TimeStamps(:,4),...
                                    LMdata.(LMname).Ta,...
                                    LMdata.(LMname).Tdry, ...
                                    LMdata.(LMname).JDD, ...
                                    LMdata.(LMname).JD,...
                                    LMdata.(LMname).sNum(1), graph, 'Dry Leaf IRT (\circC)' );                                    
 hold on



   if LM_list(g)~=4 % don't adjust the RH offset for the standard! 
       % (because that would be calibrating LM 4 RH sensor against itself!)
       graph='n';
       if graph=='y'
            figure(202)
            subplot(4,4,g)   
       end
        LMdata.(LMname) = RH_offset( LMdata.(LMname), LMdata.LM4, graph );
        hold on
   end


% Moving Average Filter

   % Air Temperature
   graph='n';
   if graph=='y'
       figure(203)
       subplot(4,4,g)   
   end
[LMdata.(LMname).Ta_fil ] = moving_average( LMdata.(LMname).Ta,...
                                                LMdata.(LMname).JDD, graph, ...
                                                LMdata.(LMname).sNum(1), ...
                                                'Air Temperature (\circC)' );
    hold on

   % Leaf Temperature
   graph='n';
   if graph=='y'
       figure(204)
       subplot(4,4,g)
   end
[ LMdata.(LMname).Tl_fil ] = moving_average( LMdata.(LMname).Tl,...
                                                LMdata.(LMname).JDD, graph, ...
                                                LMdata.(LMname).sNum(1), ...
                                                'Leaf Temperature (\circC)' );
   hold on
                                            
   % Dry Leaf Temperature         
   graph='n';
   if graph=='y'
       figure(205)
       subplot(4,4,g)
   end
[ LMdata.(LMname).Tdry_fil ] = moving_average( LMdata.(LMname).Tdry,...
                                                LMdata.(LMname).JDD, graph, ...
                                                LMdata.(LMname).sNum(1), ...
                                                'Dry Leaf Temperature (\circC)' );
   hold on
                                            
   % Relative Humidity 
   graph='n';
   if graph=='y'
       figure(206)
       subplot(4,4,g)
   end
[ LMdata.(LMname).RH_fil ] = moving_average( LMdata.(LMname).RH,...
                                                LMdata.(LMname).JDD, graph, ...
                                                LMdata.(LMname).sNum(1), ...
                                                'Relative Humidity (%)' );
    hold on

   % PAR         
   graph='n';
   if graph=='y'
       figure(207)
       subplot(4,4,g)
   end
[ LMdata.(LMname).PAR_fil ] = moving_average( LMdata.(LMname).PAR,...
                                                LMdata.(LMname).JDD, graph, ...
                                                LMdata.(LMname).sNum(1), ...
                                                'PAR (\mumol/m^2s)' );
    hold on

% Compute Tair-Tleaf and Tdry-Tleaf with smoothed data
LMdata.(LMname).Tdiff_fil = LMdata.(LMname).Ta_fil - LMdata.(LMname).Tl_fil;
LMdata.(LMname).Tdiffdry_fil = LMdata.(LMname).Tdry_fil - LMdata.(LMname).Tl_fil; 
    

% Compute es(Tl)-e(Ta) and VPD after the data has been processed                                            
[ LMdata.(LMname).E_fil, LMdata.(LMname).VPD_fil ] ...
                                = env_calcs(LMdata.(LMname).Tl_fil,...
                                LMdata.(LMname).Ta_fil,...
                                LMdata.(LMname).RH_fil);                                            

end
    

% This should be the last line (but before the nested functions)
LMdata_processed = LMdata;

%% Data Processing Functions

function [ data_i ] = t_rh_cal(data_i, data_s, p, fig1, fig2, g)
% Inputs:
% data_i: data of the ith leaf monitor 
% data_s: data of the standard (s for standard) leaf monitor to which the
% ith leaf monitor will be calibrated to
% p should be 'y' to display plots of calibrated vs uncalibrated data.

% Tair Calibration

   % Least Squares Regression
   x = data_i.Ta;
   y = data_s.Ta;
   X = [ones(length(x),1), x];
   b = X\y;
   yCalc = X*b;
   Rsq = 1 - sum((y - yCalc).^2)/sum((y - mean(y)).^2); % R-Squared
   data_i.Ta = yCalc; % calibrate Ta
   if p=='y'
       figure(fig1)
       subplot(4,4,g)
       scatter(x,y,'.')
       hold on
       scatter(x,yCalc,'.')
       title(['LM ', num2str(data_i.sNum(1)), ' Tair: R^2=',num2str(Rsq),...
           ', y=',num2str(b(2)),'*x+',num2str(b(1))])
       xlabel('Input')
       ylabel('Output')
       hold on
   end
 
% RH calibration

% Least Squares Regression
   x = data_i.RH;
   y = data_s.RH;
   X = [ones(length(x),1), x];
   b = X\y;
   yCalc = X*b;
   Rsq = 1 - sum((y - yCalc).^2)/sum((y - mean(y)).^2); % R-Squared
   data_i.RH = yCalc; % calibrate Ta
   if p=='y'
       figure(fig2)
       subplot(4,4,g)
       scatter(x,y,'.')
       hold on
       scatter(x,yCalc,'.')
       title(['LM ', num2str(data_i.sNum(1)), ' RH: R^2=',num2str(Rsq),...
           ', y=',num2str(b(2)),'*x+',num2str(b(1))])
       xlabel('Input')
       ylabel('Output')
       hold on
   end

end



function [ IRT_adj, offset_dates ] = Tleaf_offset_daily( hour, Tair, IRT,...
                                      JDD, JD, leafmon, p, description )
% This function adjusts for the offset in Tair-Tleaf. A unique offset is 
% calculated each day. We know that Tair
% should not be lower than Tleaf at night. We are trusting the number on
% the Tair sensor and blaming the Tleaf sensor for the negative offset
% in Tair-Tleaf during the night.
% Input: IRT data (must run this function for every IRT sensor you have on
% the leaf monitor. The 2017 version of the leaf monitor has two IRT
% sensors, so this function needs to be run for both IRT sensors. 

%Calculate the Tdiff offset by summing all of the Tdiff between the time
%period chosen below. 

JD_unique = unique(JD);
i=zeros(size(JD_unique,1),1);
sum=zeros(size(JD_unique,1),1);
offset=zeros(size(JD_unique,1),1);
IRT_adj=zeros(size(IRT,1),1); 
% pre-allocate array of zeros the size of the number of IRT data points

for k = 1:size(JD_unique,1) % for each unique day    
    n = find(JD == JD_unique(k)); 
    % find the indices of all the data points that correspond to day k
    
    for j=n(1):n(end) % for each index corresponding to day k
         if hour(j) >= 1 && hour(j) < 4 % if the hour is between 1 and 4 AM
             i(k) = i(k) + 1; % add one to the counter
             sum(k) = sum(k) + (Tair(j) - IRT(j)); 
             % add the offset between Tair and the IRT to the sum
         end
    end

    offset(k) = sum(k)/i(k); % calculate the offset by dividing the sum of 
    %the offsets by the number of data  points between 1 and 4 AM
    offset_dates = [JD_unique, offset];
    
    for j=n(1):n(end) % for each index corresponding to day k
        IRT_adj (j) = IRT(j) + offset(k); % + or - the offset to every 
        %data point of the IRT sensor. (Note: If the offset is negative, 
        %the addition will turn into a subtraction sign).
    end


end

if p=='y'
    
    plot(JDD,Tair-IRT,... 
        JDD,Tair-IRT_adj)
    % plot original IRT data versus offset adjusted IRT data
    xlabel('Julian Decimal Date')
    ylabel('Tair - IRT Sensor')
    title(['LM ' num2str(leafmon) ': ', description])
    legend('Original Data','Offset Adjusted')
    xlim([200,210])
    
    hold on
    
    n=find(hour >= 1 & hour < 4); % find all data points between 1 and 4 AM
    scatter(JDD(n),Tair(n)-IRT(n),'xb') % plot in 1-4 AM data points with a 
    %different symbol to show what data points were used for adjusting 
    %the offset
    scatter(JDD(n),Tair(n)-IRT_adj(n),'xr')
    
    grid on
    
end

end


function [ data_i_adj ] = RH_offset( data_i, data_s, p )
    % This function calibrates each relative humidity sensor for the ith leaf 
    % monitor to a standard relative humidity sensor. Each relative humdity 
    % sensor will be adjusted to match the relative humidity of the standard 
    % between 1 AM and 4 AM. The averages RH of the standard and the ith leaf
    % monitor will be computed from 1 AM to 4 AM and the difference between the
    % averages of the two sensors will be considered the RH offset. 


    % 1-4 AM Relative Humidity data for ith sensor
    morning_RH_i = data_i.RH( data_i.TimeStamps(:,4) >= 4 ...
        & data_i.TimeStamps(:,4) < 7);    

    % 1-4 AM Relative Humidity data for standard sensor
    morning_RH_s = data_s.RH( data_s.TimeStamps(:,4) >= 4 ...
        & data_s.TimeStamps(:,4) < 7); 
    
    % average RH of all 1-4 AM RH data throughout the season for ith sensor    
    avg_RH_i = mean(morning_RH_i,1);

    % average RH of all 1-4 AM RH data throughout the season for standard 
    %sensor
    avg_RH_s = mean(morning_RH_s,1);

    % difference between the standard and the ith sensor
    offset = avg_RH_s - avg_RH_i;

    % Add the offset to the raw relative humidity data of the ith sensor
    RH_adj = data_i.RH + offset;


        if p=='y'
            plot(data_i.JDD,data_i.RH,'b')
            hold on
            plot(data_i.JDD, RH_adj,'r')
            plot(data_s.JDD,data_s.RH,'g')
            title(['Relative Humidity LM ', num2str(data_i.sNum(1)),...
                ' versus Standard LM'])
            legend('Tested Sensor','Tested Sensor Adjusted','Standard Sensor')
            ylabel('Relative Humidity')
            xlabel('Julian Decimal Date')
            xlim([195 205])

            hold on

            n=find(data_i.TimeStamps(:,4) >= 4 & data_i.TimeStamps(:,4) < 7); 
            % find all data points between 1 and 4 AM
            scatter(data_i.JDD(n),data_i.RH(n),'xb') % plot in 1-4 AM data 
            %points with a different symbol to show what data points were 
            %used for adjusting the offset
            scatter(data_i.JDD(n),RH_adj(n),'xr')
            xlim([195 205])

            hold on

            m=find(data_s.TimeStamps(:,4) >= 4 & data_s.TimeStamps(:,4) <7);
            scatter(data_s.JDD(m),data_s.RH(m),'xg')
            xlim([195 205])
            
            grid on
            
           
            
        end
        
        data_i.RH = RH_adj;
        data_i_adj = data_i;

end


function [ FilteredData ] = moving_average( RawData, TimeStamps, fig, sNum, description )
% This is the function that performs a moving average on the data

    sampleRate=4; % four data points per hour
        numDataPointsToAvg = 5;
        coeffMovA = ones(1,numDataPointsToAvg)/numDataPointsToAvg;
        MovAvg = filter(coeffMovA, 1, RawData);

%     if fig=='y'
%         plot(TimeStamps,[RawData, MovAvg])
%         legend('Raw Data','Hourly Average (delayed)','location','best')
%         xlabel('Time')
%         ylabel('Raw Data')
%         title('Raw Data Versus Moving Average Filtered Data (delayed)')
%     end

    % Eliminating Moving Average Filter Delay
    MovAvgDelay = (length(coeffMovA)-1)/2; 
    noDelayTimeStamps = TimeStamps-MovAvgDelay/(24*sampleRate); 
    % adjust for X minutes ahead

    if fig=='y'
        plot(TimeStamps,RawData,... 
            noDelayTimeStamps,MovAvg)
        legend('Raw Data','Smoothed','location','best')
        xlabel('Julian Decimal Date')
        ylabel('Raw Data')
        title(['LM ',num2str(sNum), ' ', description])
        xlim([195 205])
    end

    FilteredData = MovAvg;

end



function [ E, VPD ] = env_calcs(Tl,Ta,RH)
%This function calculates es(Tl)-e(Ta) and the vapor pressure deficit. 
%   You must have the Tleaf, Tair, and Relative Humidity data. You must
%   have the same number of data points for each of these types of data,
%   otherwise matrix dimensions will not agree.

xl=(17.27.*Tl)./(Tl+237.3); %xl, the exponent part of es(Tl)

xa=(17.27.*Ta)./(Ta+237.3); %xa, the exponent part of es(Ta)

es_Tl=0.6108.*exp(xl); %es(Tl). Partial vapor pressure of the water within the stomatal cavity

es_Ta=0.6108.*exp(xa); %es(Ta). Partial vapor pressure of the water in the ambient air

e_Ta=es_Ta.*(RH./100); %e(Ta)

E=es_Tl-e_Ta; %es(Tl)-e(Ta)

VPD=es_Ta.*(1-RH./100); %VPD = vapor pressure deficit

end




end
