function [ LMdata_with_CWSI ] = CWSI( LMdata, LM_list )
% This function computes the Crop Water Stress Index (CWSI) using data from
% the leaf monitor.

% Saturated Leaf Temperature
satdays = unique(LMdata.LM16.JD);    
Tsat_m_avg = zeros(length(satdays),1);
for f=1:length(satdays)
    m = LMdata.LM16.JD==satdays(f) & LMdata.LM16.time>=12 & ...
        LMdata.LM16.time<=16; % 12 pm - 4 pm data on day f
    Tsat_m_avg(f) = mean(LMdata.LM16.Tl_fil(m)); % average saturated leaf 
    % temperature between 12pm and 4 pm on day f
end
    
for g=1:length(LM_list)
    
    LMname = ['LM', num2str(LM_list(g))];
    
    % Initialize variables
    Tair_n_avg = zeros(length(LMdata.(LMname).days),1);
    Tleaf_n_avg = zeros(length(LMdata.(LMname).days),1);
    Tdry_n_avg = zeros(length(LMdata.(LMname).days),1);
    
    n = ismember(satdays, LMdata.(LMname).days); % find what days of the 
    % saturated leaf temp are also members of the days of the monitored leaf
    
    Tsat_n_avg = Tsat_m_avg(n); % only get saturated leaf temp from days 
    % when the monitored leaf had data
    
    for f=1:length(LMdata.(LMname).days)
         
        n = LMdata.(LMname).JD==LMdata.(LMname).days(f) ...
            & LMdata.(LMname).time>=12 & ...
            LMdata.(LMname).time<=16; % 12 pm - 4 pm
    
        % Air Temperature
        Tair_n_avg(f) = mean(LMdata.(LMname).Ta_fil(n));
        
        % Live Leaf Temperature
        Tleaf_n_avg(f) = mean(LMdata.(LMname).Tl_fil(n));
    
        % Dry Leaf Temperature
        Tdry_n_avg(f) = mean(LMdata.(LMname).Tdry_fil(n));
                   
    end
    
      % Temperature Differentials
      dTm = Tair_n_avg - Tleaf_n_avg;
      dTll = Tair_n_avg - Tsat_n_avg;
      dTul = Tair_n_avg - Tdry_n_avg;
      
       % CWSI
      CWSI_n = (dTm - dTll)./(dTul - dTll);

      LMdata.(LMname).CWSI = CWSI_n;
    
end

LMdata_with_CWSI = LMdata;

end

