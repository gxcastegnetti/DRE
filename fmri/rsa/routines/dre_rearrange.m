function arrData = dre_rearrange(dir,subs,taskOrd,bData)
%% function dre_rearrange(dirSub,sub,runType)
% ~~~
% INPUTS:
%   dir: struct of directories
%   subj: subject numbers
%   taskOrd: FBFB or BFBF
%   bData: behavioural data output of dre_extractData
% OUTPUTS:
%   bData: struct with trial onsets and behavioural measures
% ~~~
% GX Castegnetti --- 2018

fs = filesep;
n_sess = 4;

dir.data = [dir.dre,fs,'data'];

%% read objects
objs        = readtable([dir.beh,fs,'Objects.csv']);
objsIdxCorr = [(1:120)',find(table2array(objs(:,7)))]; clear objs

%% loop over subjects and sessions
for s = 1:length(subs)
    for r = 1:n_sess
        
        %% extract subject-specific task order
        if taskOrd(s) == 1
            day2Order = {'F','B','F','B'};
        elseif taskOrd(s) == 2
            day2Order = {'B','F','B','F'};
        end
        
        %% load behavioural data day 2
        dirPsyO = [dir.data,fs,'fmri',fs,'psychOut',fs,'SF',num2str(subs(s),'%03d')];
        Mday2 = csvread([dirPsyO,fs,'DRE_mri_S',num2str(subs(s),'%03d'),'_B',num2str(r),day2Order{r},'.csv']);
        
        %% object IDs
        % extract IDs (not idx) from day 2
        sessions(:,r) = Mday2(Mday2(:,3)==0,4);
        
    end
    
    %% arrange according to presentation order
    % replace obj ID with obj idx
    arrSessions = NaN(size(sessions));
    for i = 1:120
        if taskOrd(s) == 1
            arrSessions(sessions(:,1) == objsIdxCorr(i,2),1) = objsIdxCorr(i,1);
            arrSessions(sessions(:,2) == objsIdxCorr(i,2),2) = objsIdxCorr(i,1) + 120;
            arrSessions(sessions(:,3) == objsIdxCorr(i,2),3) = objsIdxCorr(i,1);
            arrSessions(sessions(:,4) == objsIdxCorr(i,2),4) = objsIdxCorr(i,1) + 120;
        else
            arrSessions(sessions(:,1) == objsIdxCorr(i,2),1) = objsIdxCorr(i,1) + 120;
            arrSessions(sessions(:,2) == objsIdxCorr(i,2),2) = objsIdxCorr(i,1);
            arrSessions(sessions(:,3) == objsIdxCorr(i,2),3) = objsIdxCorr(i,1) + 120;
            arrSessions(sessions(:,4) == objsIdxCorr(i,2),4) = objsIdxCorr(i,1);
        end
    end
    
    arrData(subs(s)).acc2sessions = arrSessions(:);
    
    %% arrange according to value (3 levels)
    bData(subs(s)).imagination(r).val
    valAll = 
    prctile
    
    
end