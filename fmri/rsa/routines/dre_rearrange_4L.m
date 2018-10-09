function arrData = dre_rearrange_4L(dir,subs,taskOrd,bData)
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
        sessions(:,r) = Mday2(Mday2(:,3)==0,4); %#ok<*AGROW>
        
    end
    
    %% arrange according to presentation order
    % replace obj ID with obj idx
    norm2sessions = NaN(size(sessions));
    for i = 1:120
        if taskOrd(s) == 1
            norm2sessions(:,1) = bData(subs(s)).imagination(1).objIdx;
            norm2sessions(:,2) = bData(subs(s)).imagination(2).objIdx + 120;
            norm2sessions(:,3) = bData(subs(s)).imagination(3).objIdx;
            norm2sessions(:,4) = bData(subs(s)).imagination(4).objIdx + 120;
        else
            norm2sessions(:,1) = bData(subs(s)).imagination(1).objIdx + 120;
            norm2sessions(:,2) = bData(subs(s)).imagination(2).objIdx;
            norm2sessions(:,3) = bData(subs(s)).imagination(3).objIdx + 120;
            norm2sessions(:,4) = bData(subs(s)).imagination(4).objIdx;
        end
    end
    
    % put in output struct
    arrData(subs(s)).norm2sessions = norm2sessions(:);
    
    %% arrange according to value (3 levels)
    
    sessIdx = norm2sessions(:);
    
    %%%%%%%%%
    % value %
    %%%%%%%%%
    
    valAll = [bData(subs(s)).imagination(1).val;
        bData(subs(s)).imagination(2).val;
        bData(subs(s)).imagination(3).val;
        bData(subs(s)).imagination(4).val];
    
    valAll(isnan(valAll)) = 1;
    
    % perturb values for univoque percentiles calculation
    valAll = valAll + 0.00000001*(1:length(valAll))';
    prctile25 = prctile(valAll,25);
    prctile50 = prctile(valAll,50);
    prctile75 = prctile(valAll,75);
    
    val_L = find(valAll < prctile25);
    val_m = find(valAll > prctile25 & valAll < prctile50);
    val_M = find(valAll > prctile50 & valAll < prctile75);
    val_H = find(valAll > prctile75);
    
    val_L = val_L(randperm(length(val_L)));
    val_m = val_m(randperm(length(val_m)));
    val_M = val_M(randperm(length(val_M)));
    val_H = val_H(randperm(length(val_H)));
    valAll_4L = [val_L; val_m; val_M; val_H];
    
    [~,valSort] = sort(valAll);
    
    % put in output struct
    arrData(subs(s)).norm2val = sessIdx(valAll_4L);
%     arrData(subs(s)).norm2val = sessIdx(valSort);
    
    %%%%%%%%%%%%%%%
    % familiarity %
    %%%%%%%%%%%%%%%
    
    famAll = [bData(subs(s)).imagination(1).fam;
        bData(subs(s)).imagination(2).fam;
        bData(subs(s)).imagination(3).fam;
        bData(subs(s)).imagination(4).fam];
    
    famAll(isnan(famAll)) = 1;
    
    % perturb familiarity for univoque percentiles calculation
    famAll = famAll + 0.00000001*(1:length(famAll))';
    prctile33 = prctile(famAll,100/3);
    prctile66 = prctile(famAll,200/3);
    fam_L = find(famAll < prctile33);
    fam_M = find(famAll > prctile33 & famAll < prctile66);
    fam_H = find(famAll > prctile66);
    fam_L = fam_L(randperm(length(fam_L)));
    fam_M = fam_M(randperm(length(fam_M)));
    fam_H = fam_H(randperm(length(fam_H)));
    famAll_3L = [fam_L; fam_M; fam_H];
    
    [~,famSort] = sort(famAll);
    
    % put in output struct
    arrData(subs(s)).norm2fam = sessIdx(famAll_3L);
%     arrData(subs(s)).norm2fam = sessIdx(famSort);
    
end