function [alldata, LM_list] = leaf_ratio_2017_adjustments(alldata_pre)

% Fix 2017 LM Changes    
    serial = x2mdate(alldata_pre(:,1),0);
    timestamps = datevec(datestr(serial));
    JD = Date_to_Julian(timestamps(:,2), timestamps(:,3), 'n');
    alldata_pre = [alldata_pre,JD]; % concatenate JD to alldata_JD
    [ splitLM ] = split_leafmon(alldata_pre,alldata_pre(:,12));
    [ split_corrected, ~ ] = correct_LM_changes( splitLM );
    alldata = [ ];
    for k=1:length(split_corrected)    
       alldata = [alldata; split_corrected{1,k}];    
    end   
    LM_list = [1;2;4;15;16;18;21;31;33;34;35;36;37;38;40;41]; % tree numbers
%     alldata = alldata( alldata(:,1)>42853,:); % only include data after April 28 
    % (this was the first FULL day all leaf monitors were installed expect 
    % for LM 37, which was installed earlier).
%     alldata = alldata( alldata(:,1)>42888,:); % only include data after June 22

% DATA OMISSIONS
    alldata = alldata( alldata(:,13)>153 & alldata(:,13)~=159,:); 
    % only include data after June 2 (Julian Date 153) and
    % delete data from June 8 (JD 159) because it was raining and paper was used
    % instead of a dry leaf...the paper got wet because of the rain, so it
    % did not simulate the dry leaf condition at all.
    
    % Tree 1
    alldata( alldata(:,12)==1 & alldata(:,13)==184,:)=[]; % from field notes, 'yellowing leaf on 184, leaf changed on 186'
    alldata( alldata(:,12)==1 & alldata(:,13)==185,:)=[]; % from field notes, 'yellowing leaf on 184, leaf changed on 186'
    alldata( alldata(:,12)==1 & alldata(:,13)==186,:)=[]; % from field notes, 'yellowing leaf on 184, leaf changed on 186'
    
    % Tree 2
    alldata( alldata(:,12)==2 & alldata(:,13)==188,:)=[]; % from field notes, 'IRT not looking at live leaf'
    
    % Tree 4
    alldata( alldata(:,12)==4 & alldata(:,13)==170,:)=[]; % from field notes 'data bad'
    alldata( alldata(:,12)==4 & alldata(:,13)==205,:)=[]; % from field notes 'stem broke'
    
    % Tree 34
    alldata( alldata(:,12)==34 & alldata(:,13)==191,:)=[]; 
    % eliminate data from LM 34 on July 10 (JD 191) because the leaf
    % shifted out of its holder (from field notes)
    
    % Tree 38
    alldata( alldata(:,12)==38 & alldata(:,13)==181,:)=[]; % suspected bug in IRT sensor or dying leaf
    alldata( alldata(:,12)==38 & alldata(:,13)==182,:)=[];% suspected bug in IRT sensor or dying leaf
    alldata( alldata(:,12)==38 & alldata(:,13)==183,:)=[];% suspected bug in IRT sensor or dying leaf
    alldata( alldata(:,12)==38 & alldata(:,13)==184,:)=[];% suspected bug in IRT sensor or dying leaf
    alldata( alldata(:,12)==38 & alldata(:,13)==185,:)=[]; % from field notes 'dying leaf'
    alldata( alldata(:,12)==38 & alldata(:,13)==186,:)=[]; % from field notes 'dying leaf, leaf changed'
    alldata( alldata(:,12)==38 & alldata(:,13)==187,:)=[]; % negative ICSR and ACSR
    alldata( alldata(:,12)==38 & alldata(:,13)==188,:)=[];% suspected bug in IRT sensor
    alldata( alldata(:,12)==38 & alldata(:,13)==189,:)=[];% suspected bug in IRT sensor
    alldata( alldata(:,12)==38 & alldata(:,13)==190,:)=[];% suspected bug in IRT sensor
    alldata( alldata(:,12)==38 & alldata(:,13)==191,:)=[];% suspected bug in IRT sensor
    alldata( alldata(:,12)==38 & alldata(:,13)==192,:)=[];% suspected bug in IRT sensor
    alldata( alldata(:,12)==38 & alldata(:,13)==193,:)=[];% bug observed in IRT sensor
    
    
    
end

%% Functions

function [ splitLM ] = split_leafmon(data,serial)
%This function splits the almond data into separate cells for each leaf
%monitor
%data: the matrix with ALL your data
%serial: the column in the matrix that contains the leaf monitor serial
%numbers

    unique_LM=unique(serial);

    for u=1:length(unique_LM)

        ind_LM=serial==unique_LM(u);

        splitLM{u}=data(ind_LM,:);

    end

end


function [ split_corrected, TN ] = correct_LM_changes( split )
% There were several leaf monitor changes throughout the 2017 season in
% almonds at Nickels. It does not make sense to identify the trees being
% studied by leaf monitor numbers because there were so many leaf monitor
% changes. This function assigns each tree being studied to a tree number 
% rather than a leaf monitor number. 
% For the trees that had the same leaf monitor during the entire season,
% the tree number is the same as the leaf monitor number. For trees that
% had leaf monitor changes during the season, the tree number will be the
% same number as the leaf monitor number when they were originally
% installed. Since there were no leaf monitor changes in 2016, we will use 
% the 2016 season leaf monitor numbers as the tree numbers for both the
% 2016 and 2017 seasons. That way, the tree numbers are the same
% between the 2016 and 2017 seasons. 

   JD = 13; % Julian Date is in column 13 of split{i,j}
    
    split_corrected = cell(1,16);

    TN = zeros(1,16); % 16 trees

    for k=1:size(split,2)

        LM = split{1,k}(1,12); % leaf monitor number for k

        if LM==4 || LM==21 || LM==36 || LM==34 || LM==16 || LM==18 || LM==15 || LM==37 || LM==33 || LM==31 % for the leaf monitors that were never moved

            split_corrected{1,k} = split{1,k};

            TN(k) = LM; % tree number equals leaf monitor number for leaf monitors that were never moved

        end

    end

    % Tree 1/Row 4 Middle 
    split_corrected{:,1} = split{:,1} ( split{:,1}(:,JD) < 191,:); % Tree 1
    TN(1) = 1;
    split_corrected{:,1}(:,12) = 1.*ones(size(split_corrected{:,1}(:,12),1),1);

    % Tree 2/Row 4 North
    split_corrected{:,2} = split{:,2}( split{:,2}(:,JD) < 193,:); % Tree 2
    TN(2) = 2;
    split_corrected{:,2}(:,12) = 2.*ones(size(split_corrected{:,2}(:,12),1),1);

    % Tree 35/Row 2 North
    fromLM35a = split{:,11} (split{:,11}(:,JD)<=135,:);
    fromLM40a = split{:,15} (split{:,15}(:,JD)>135 & split{:,15}(:,JD)<192,:);
    fromLM35b = split{:,11} (split{:,11}(:,JD)>193 & split{:,11}(:,JD)<=200,:);
    fromLM40b = split{:,15} (split{:,15}(:,JD)>200,:);
    split_corrected{:,11} = [fromLM35a; fromLM40a; fromLM35b; fromLM40b];  
    TN(11)=35;
    split_corrected{:,11}(:,12) = 35.*ones(size(split_corrected{:,11}(:,12),1),1);
    
    % Tree 38/Row 3 North
    fromLM38 = split{1,14} (split{:,14}(:,JD)<=192,:);
    fromLM1 = split{1,1} (split{:,1}(:,JD)>193,:);
    split_corrected{:,14} = [fromLM38;fromLM1];
    TN(14)=38;
    split_corrected{:,14}(:,12) = 38.*ones(size(split_corrected{:,14}(:,12),1),1);

    % Tree 40/Row 3 Middle South
    fromLM40 =  split{:,15} (split{:,15}(:,JD)<=135,:);
    fromLM2 =   split{:,2} (split{:,2}(:,JD)>200,:);
    split_corrected{:,15} = [fromLM40;fromLM2];
    TN(15)=40;
    split_corrected{:,15}(:,12) = 40.*ones(size(split_corrected{:,15}(:,12),1),1);

    % Tree 41/Row 2 South
    split_corrected{:,16} = split{:,14} (split{:,14}(:,JD)>=193,:);
    TN(16)=41;
    split_corrected{:,16}(:,12) = 41.*ones(size(split_corrected{:,16}(:,12),1),1);

end


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
