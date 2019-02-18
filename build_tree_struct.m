function [ LMdata, LM_list ] = build_tree_struct( alldata_pre )
% This function will build a structure of data pertaining to each tree. 
% alldata is the data from all the leaf monitor
% LMdata is the data from each leaf monitor #


normaldate = excel2dates(floor(alldata_pre(1,1)));
year = normaldate.Year;


if year==2017
    [alldata, LM_list] = leaf_ratio_2017_adjustments(alldata_pre);
      
elseif year==2016 % no LM changes in 2016
    alldata = alldata_pre( alldata_pre(:,1)>=42538 & alldata_pre(:,1)<42593,: ); % only use data from June 17 to August 10 (harvest date)
    LM_list = unique(alldata(:,12));

elseif year==2018 
    alldata = alldata_pre; 
    LM_list = unique(alldata(:,12));
end



for g=1:length(LM_list)

    LMname = ['LM', num2str(LM_list(g))]; % creat a string to name the leaf monitor

    LMdata.(LMname) = LM; % give each LMdata structure the properties of class LM

    n = find(LM_list(g)==alldata(:,12)); % find row numbers of alldata that pertain to LM g
    
    if year==2016        
        d = 40==alldata(:,12); % find row numbers of alldata that pertain to LM 40
        % For 2016 Data: Interpolation with respect to dry leaf monitor (LM 40) goes here
        alldata_i = interpolate_DryLM( alldata(d,:), alldata(n,:), 'n');    
    elseif year==2017            
        % Eliminate repeated time stamps. Important for LM 35.
        LM_i = alldata(n,:);
        [~,ind] = unique(LM_i(:,1),'rows');    
        alldata_i = LM_i(ind,:);
        
         if LM_i(1,12)==40 || LM_i(1,12)==41
             n = alldata_i(:,1)>=42937; % only include data on or after July 21 for
             % trees 40 & 41
             alldata_i = alldata_i(n,:);
         end
         
    else
        % Eliminate repeated time stamps. 
        LM_i = alldata(n,:);
        [~,ind] = unique(LM_i(:,1),'rows');    
        alldata_i = LM_i(ind,:);       
    end

   

    LMdata.(LMname).Tl = alldata_i(:,3); % leaf temperature

    LMdata.(LMname).RH = alldata_i(:,4); % relative humidity

    LMdata.(LMname).Ta = alldata_i(:,5); % air temperature 
    
    LMdata.(LMname).wind = alldata_i(:,6); % wind speed

    LMdata.(LMname).Tdiff = alldata_i(:,5) - alldata_i(:,3); % air temperature - leaf temperature 

    LMdata.(LMname).PAR = alldata_i(:,9); % photosynthetically active radiation (PAR)

    LMdata.(LMname).sNum = alldata_i(:,12); % serial number/leaf monitor number
    
    LMdata.(LMname).Tdry = alldata_i(:,7); % dry leaf temperature 

    LMdata.(LMname).Tdiffdry = alldata_i(:,7) - alldata_i(:,3); % dry leaf temperature - live leaf temperature 
    
    
    % Time Stamps
    serial = x2mdate(alldata_i(:,1),0);
    LMdata.(LMname).TimeStamps = datevec(datestr(serial));
    
    % Julian Days
    if year==2016
         LMdata.(LMname).JD = Date_to_Julian(LMdata.(LMname).TimeStamps(:,2),...
                             LMdata.(LMname).TimeStamps(:,3),'y');
    elseif year==2017 || year==2018
         LMdata.(LMname).JD = Date_to_Julian(LMdata.(LMname).TimeStamps(:,2),...
                             LMdata.(LMname).TimeStamps(:,3),'n');
    end
    LMdata.(LMname).days = unique(LMdata.(LMname).JD); % only the unique Julian days
    
    % Julian Decimal Days
    LMdata.(LMname).JDD = Julian_Decimal_Days( LMdata.(LMname).JD,...
                                        LMdata.(LMname).TimeStamps(:,4),...
                                        LMdata.(LMname).TimeStamps(:,5),...
                                        LMdata.(LMname).TimeStamps(:,6) );
                                    
    % Decimal Hours
    [ LMdata.(LMname).time ] = Decimal_Hours(LMdata.(LMname).TimeStamps(:,4),...
                                             LMdata.(LMname).TimeStamps(:,5));

                                    
      [ LMdata.(LMname).E,...
        LMdata.(LMname).VPD ] = env_calcs(LMdata.(LMname).Tl,...
                                LMdata.(LMname).Ta,...
                                LMdata.(LMname).RH);
    
    
end

 
%% Functions




function [ Julian_Day ] = Date_to_Julian(month,day,leap)
    %This is a simple function for converting a calender day to a Julian day
    % The inputs are (month,day,leap). Input the month for month, day for day,
    % and for leap put 'y' if a leap year, or 'n' if not.
    md=[31,59,90,120,151,181,212,243,273,304,334,365];  %Days of the year if it is not a leap year
    mdl=[31,60,91,121,152,182,213,244,274,305,335,366]; %Days of the year if it is a leap year
 
    Julian_Day = zeros(length(month),1);
    for f=1:length(month)
        if (leap=='y'|| leap=='Y')
            Julian_Day(f)=mdl(month(f)-1)+day(f);
        else
            Julian_Day(f) = md(month(f)-1)+day(f);
        end
    end
    
end


function [ JDD ] = Julian_Decimal_Days( JD, Hour, Minute, Second )
 
    JDD = JD + Hour/24+...
        Minute/(24*60) + Second/(24*60*60); 
    % Julian decimal dates

end


function [ DecHours ] = Decimal_Hours(Hour, Minute)
        
    DecHours=zeros(length(Hour),1); % placeholder for decimal hours to speed up program

    for g=1:length(Hour)
         DecHours(g) = Hour(g) + Minute(g)./60; %i.e. 8:30 a.m. would be 8.5 in 
         %decimal hours and 16.25 would be 4:15 pm in decimal hours (military
         %time plus fraction of minutes in an hour)
    end

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


function [ LM40_fixed ] = fix_LM40_2017( LM40 )
    % The dry IRT and the live IRT were switched on July 24. Prior
    % to July 24, the live IRT was looking at a live leaf and the dry IRT was
    % looking at a dry leaf. On July 24, they were switched.
    n = find( LM40.sNum==40 & LM40.JD>205 );

    % temporarily rename Tdry and Tleaf
    Tleaf_LM40_205 = LM40.Tl(n);
    Tdry_LM40_205 = LM40.Tdry(n);

    % switch them
    LM40.Tl(n) = Tdry_LM40_205;
    LM40.Tdry(n) = Tleaf_LM40_205;

    LM40_fixed = LM40;

end


function [ LM_i_int ] = interpolate_DryLM( DryLM, LM_i, p )
% We want to interpolate each leaf monitor to the timestamps of the dry
% leaf monitor, so that the timestamps for all leaf monitors match the dry
% leaf monitor.

    figure 

    % Eliminate repeated time stamps
    % LM_i
    [~,ind] = unique(LM_i(:,1),'rows');
    LM_i = LM_i(ind,:);
   
    % Dry_i
    [~,ind] = unique(DryLM(:,1),'rows');
    DryLM = DryLM(ind,:);

    % Initialize LM_i interpolated
    LM_i_int = zeros(size(DryLM));
    
    k=1;
 
    for f = [3,4,5,6,9] % skip the first row, which contains the time stamps 

        x = LM_i(:,1); % LM_i Excel time stamps
        y = LM_i(:,f); % data to be interpolated
        xq = DryLM(:,1); % x query data (Excel timestamps of Dry LM)

        s = interp1(x,y,xq,'spline',exp(1)); % make points outside the 
        % domain of x equal to e, which is a very unique number that can
        % easily be identified in the next if/elseif statement.
        
        % The point in this part is to identify if the last and/or the first
        % element of the data was outside the domain of x and was given a
        % default value of e. I used e because I know that no other number
        % in the dataset will be perfectly equal to e, so no other number
        % will unintentionally be 'picked up' by this if statement. 
        if s(end)==exp(1) % if e was put in the last column of s
            s(end)=s(end-1); % replace e with the nearest data point value, 
            % the second to last one
        elseif s(1)==exp(1) % if e was put in the first column of s
            s(1)=s(2); % replace e with the nearest data point value, 
            % the second one
        end 

        LM_i_int(:,f) = s; % set the output as the interpolated data 

        if p=='y'
            % Plot
            subplot(3,3,k)
            plot(x,y,'.')
            hold on
            plot(xq,s,'.')
            xlim([4.258e4  4.259e4])
            title(['Interpolated LM ', num2str(LM_i(1,12)),' data with dry leaf monitor as query'])
            xlabel('Time (Excel)')
            if f==3
                 ylabel('Leaf Temperature')
            elseif f==4
                ylabel('Relative Humidity')
            elseif f==5
                ylabel('Air Temperature')
            elseif f==6
                ylabel('Wind Speed')
            elseif f==9
                ylabel('Photosynthetically Active Radiation')
            end
            legend('Interpolated','Query')
            k=k+1; % counter to move up subplot index each time a data type as been plotted
            hold off
        end
           
    end

    LM_i_int(:,1) = DryLM(:,1); % LM_i now has the same time stamps as the Dry LM
    
    LM_i_int(:,8) = LM_i_int(:,5) - LM_i_int(:,3); % instead of interpolating Tdiff, 
    % recompute Tdiff with the interpolated Tair and Tleaf
    
    LM_i_int(:,12) = LM_i(1,12)*ones(length(LM_i_int),1); % serial number of LM
    
    % Entering dry leaf temperature into LM property (because the 2016 leaf
    % monitors did not have unique dry leaf sensors for each tree)
    LM_i_int(:,7) = DryLM(:,3); % put dry leaf temperature data into LM_i data matrix
    
 
    % remove rows with NaNs (using leaf temperature column to remove NaNs).
    % Anywhere there are a NaNs in the leaf temperature data, there should
    % also be NaNs in the other data types.
    LM_i_int(isnan(LM_i_int(:,3)),:)=[]; 
    
end





end

