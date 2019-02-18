function [LMfile] = get_LM_filepath( )
        [FileName,PathName,FilterIndex] = uigetfile('*.*',  'All Files (*.*)');
        LMfile = strcat(PathName,FileName);        
       
%         
%         fl=load('flag.mat');
%         switch fl.flag
%             case 1
%                 save(strcat(pwd,'\files\LMfile_amd.mat'),'LMfile');
%             case 2
%                 save(strcat(pwd,'\files\LMfile_grp.mat'),'LMfile');
%         end

    end


