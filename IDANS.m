function [ LMdata_fin ] = IDANS( LMdata, graphing_period, LM_list, graph )
% This function computes the Degrees Above Non-Stressed (DANS) and 
% Integrated Degrees Above Non-Stressed (IDANS) using the leaf monitor data.

% Non-stressed canopy temperature
Tc_ns = LMdata.LM16.Tl_fil;

F=1;

for g=1:length(LM_list)
    
   LMname = ['LM', num2str(LM_list(g))];
   
   days = unique(LMdata.(LMname).JD);
   
   IDANS_n = zeros(length(days),1); % initialize variable   
       
   % Interpolate the saturated leaf temperature so it has the same number
   % of data points as leaf monitor g.
   s = interp1(LMdata.LM16.JDD, Tc_ns, LMdata.(LMname).JDD,'spline',exp(1)); 

   if s(end)==exp(1) 
       s(end)=s(end-1); 
   elseif s(1)==exp(1) 
       s(1)=s(2); 
   end 
   
   Tc_ns_interp = s;
   
   LMdata.(LMname).DANS = LMdata.(LMname).Tl_fil - Tc_ns_interp; % DANS
   
   % IDANS
   for f=1:length(days) % integrate daily      
       n = LMdata.(LMname).JD == days(f); % find row numbers of day f
       IDANS_n(f) = sum(LMdata.(LMname).Tl_fil(n) - Tc_ns_interp(n)); % integral   
   end
   
   LMdata.(LMname).IDANS = IDANS_n;
  
   
   if graph=='y'
        % Plot the DANS vs time for a certain graphing period 
        n = find( LMdata.(LMname).JD(:)>=graphing_period(1)...
            & LMdata.(LMname).JD(:)<=graphing_period(end));

        set(gcf, 'Position', get(0, 'Screensize')); 

        figure(500)
        subplot(4,4,F)
        days=unique(LMdata.(LMname).JD(n));  % find the unique days in the analysis period 

                for k=1:numel(days(:,1)) % for each day in the selected 
                    % analysis period  

                    n_i = find(LMdata.(LMname).JD==days(k,1)); 
                    % determine the matrix indices of data points that fall on day 
                    % k in the analysis period        

                    x = LMdata.(LMname).time(n_i);
                    y = LMdata.(LMname).DANS(n_i) ;

                    plot(x,y);

                    grid on           

                    ylim([-1 6]);
                    ylabel('DANS')

                    xlim([0 24]);
                    xlabel('Time (Hr)')

                    title(LMname)
                                      
                    hold on

                end

   end
    
    F=F+1;   
    
    grph = cell(1,length(days));
    for p=1:length(days)
         grph{1,p} = num2str(days(p));  
    end
    h_legend = legend(grph);
    set(h_legend,'FontSize',6) 
    
end


LMdata_fin = LMdata; 

end

