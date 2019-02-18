function [ LMdata_with_ratio ] = Stom_Cond_Ratio( LMdata, graphing_period, LM_list, graph )
% This function computes a canopy temperature-based index that indicates
% the plant water status (PWS). 
% This function plots (Tleaf-Tdry)/(es(Tleaf)-e(Tair)) vs time 

f=1;
for g=1:length(LM_list)

    LMname = ['LM', num2str(LM_list(g))]; 
    % create a string to name the leaf monitor

    % Compute the ratio for the entire growing season
    LMdata.(LMname).ratio = (LMdata.(LMname).Tdry_fil -...
    LMdata.(LMname).Tl_fil)./LMdata.(LMname).E_fil;

    if graph=='y'
        % Plot the ratio vs time for a certain graphing period (it would be
        % silly to try to plot the entire season's ratio plots on one graph!
        % That would be WAY too much information on one plot, so let's just
        % look at a sample 3-7 days or so at a time, depending on how many days
        % you want to view at one time).
        n =  LMdata.(LMname).JD(:)>=graphing_period(1)...
            & LMdata.(LMname).JD(:)<=graphing_period(end);

        set(gcf, 'Position', get(0, 'Screensize')); 

        figure(600); 
        subplot(4,4,f)
        thedays=unique(LMdata.(LMname).JD(n));  % find the unique days in the analysis period 

                for k=1:numel(thedays(:,1)) % for each day in the selected 
                    % analysis period  

                    n_i = find(LMdata.(LMname).JD==thedays(k,1)); 
                    % determine the matrix indices of data points that fall on day 
                    % k in the analysis period        

                    x = LMdata.(LMname).time(n_i);
                    y = LMdata.(LMname).ratio(n_i) ;

                    plot(x,y);

                    grid on           

                    ylim([-1 3.5]);
                    ylabel('^{T_D-T_L}/_{e_s(T_L)-e(T_A)}')

                    xlim([0 24]);
                    xlabel('Time (Hr)')

                    title(LMname)

                    hold on

                end

    f=f+1;    

    grph = cell(1,length(thedays));
    for p=1:length(thedays)
         grph{1,p} = num2str(thedays(p));  
    end
    h_legend = legend(grph);
    set(h_legend,'FontSize',6) 
    
    end
    
end

LMdata_with_ratio = LMdata;

end

