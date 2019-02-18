function [ LMdata_with_results ] = Stom_Cond_Integral_and_Avg( LMdata, LM_list )

% Integral of Ratio
% calculating the integral of the ratio with respect to the time of day
% will contain the integral of the ratio for each day for each tree

year = LMdata.LM1.TimeStamps(1,1); 

for g=1:length(LM_list)

LMname = ['LM', num2str(LM_list(g))]; 
days = unique(LMdata.(LMname).JD); 
k=1;
% 
%     if year==2016 || (year==2017 && LM_list(g)~=35) % LM 35 had a problem in 2017 
%         % (still need to identify the problem)
%     
        for d=1:length(days)
            
            % Daily Integral of the Ratio from 10 AM to 6 PM
            n = (LMdata.(LMname).JD==days(d) & LMdata.(LMname).time>=10 ...
            & LMdata.(LMname).time<=18); 
        
%            m = length(LMdata.(LMname).time(n));
%            area = zeros(m,1);
%            for f=1:m-1     
%                dt = abs(LMdata.(LMname).time(f+1) - LMdata.(LMname).time(f));
%                area(m)  = (1/2).*(LMdata.(LMname).ratio(m+1)...
%                    + LMdata.(LMname).ratio(m)).*dt;          
%            end
%            Area = sum(area);
%            LMdata.(LMname).ratio_int(k,1) = Area;
           
%            if days(d)~=199 && LM_list(g)~=35
           if LM_list(g)~=35
               LMdata.(LMname).ratio_int(k,1) = trapz(LMdata.(LMname).time(n),...
               LMdata.(LMname).ratio(n)); 
           end
           
            % Daily Average of the Ratio from 1 PM to 5 PM
            m = (LMdata.(LMname).JD==days(d) & LMdata.(LMname).time>=13 ...
            & LMdata.(LMname).time<=15);

            LMdata.(LMname).ratio_daily_avg(k,1) = mean(LMdata.(LMname).ratio(m));

            k=k+1;

        end
        
%     end

end

LMdata_with_results = LMdata;

end

