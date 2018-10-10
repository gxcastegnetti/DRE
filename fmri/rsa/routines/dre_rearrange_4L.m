function arrIdx = dre_rearrange_4L(dir,subs,taskOrd,bData)
%% function dre_rearrange_4L(dirSub,sub,runType)
% finds indices to arrange response patterns according to four levels of
% value and familiarity, or in order of presentation
% ~~~
% INPUTS:
%   dir: struct of directories
%   subj: subject numbers
%   taskOrd: FBFB or BFBF
%   bData: behavioural data output of dre_extractData
% OUTPUTS:
%   arrIdx: indices to rearrange response patterns
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
    arrIdx(subs(s)).norm2sessions = norm2sessions(:);
    
    %% arrange according to value (3 levels)
    
    sessIdx = norm2sessions(:);
    
    %%%%%%%%%
    % value %
    %%%%%%%%%
    
    % put all value in the same vector
    valAll = [bData(subs(s)).imagination(1).val;
        bData(subs(s)).imagination(2).val;
        bData(subs(s)).imagination(3).val;
        bData(subs(s)).imagination(4).val];
        
    valAll(isnan(valAll)) = floor(50*rand(1,sum(isnan(valAll))));
    
    % perturb value for univoque percentiles calculation
    valAll = valAll + 0.00000001*(1:length(valAll))';
    
    % computes percentiles
    prctile25 = prctile(valAll,25);
    prctile50 = prctile(valAll,50);
    prctile75 = prctile(valAll,75);
    
    % divide into four levels and stack them
    val_L = find(valAll < prctile25);
    val_m = find(valAll > prctile25 & valAll < prctile50);
    val_M = find(valAll > prctile50 & valAll < prctile75);
    val_H = find(valAll > prctile75);
    
    % repeatable random permutation
    rng(1)
    val_L = val_L(randperm(length(val_L)));
    val_m = val_m(randperm(length(val_m)));
    val_M = val_M(randperm(length(val_M)));
    val_H = val_H(randperm(length(val_H)));
    
    % put them together 
    valAll_4L = [val_L; val_m; val_M; val_H];
    
    % sort value for the continuous arrangement
    [~,valSort] = sort(valAll);
    
    % put in output struct
    arrIdx(subs(s)).norm2val_disc = sessIdx(valAll_4L);
    arrIdx(subs(s)).norm2val_cont = sessIdx(valSort);
    
    %%%%%%%%%%%%%%%
    % familiarity %
    %%%%%%%%%%%%%%%
    
    % put all familiarities in one vector
    famAll = [bData(subs(s)).imagination(1).fam;
        bData(subs(s)).imagination(2).fam;
        bData(subs(s)).imagination(3).fam;
        bData(subs(s)).imagination(4).fam];
    
    famAll(isnan(famAll)) = floor(50*rand(1,sum(isnan(famAll))));
    
    % perturb familiarity for univoque percentiles calculation
    famAll = famAll + 0.00000001*(1:length(famAll))';
    
    % find percentiles
    prctile25 = prctile(famAll,25);
    prctile50 = prctile(famAll,50);
    prctile75 = prctile(famAll,75);
    
    % divide into three levels and stack them
    fam_L = find(famAll < prctile25);
    fam_m = find(famAll > prctile25 & famAll < prctile50);
    fam_M = find(famAll > prctile50 & famAll < prctile75);
    fam_H = find(famAll > prctile75);
    
    % repeatable random permutation
    rng(1)
    fam_L = fam_L(randperm(length(fam_L)));
    fam_m = fam_m(randperm(length(fam_m)));    
    fam_M = fam_M(randperm(length(fam_M)));
    fam_H = fam_H(randperm(length(fam_H)));
    
    % put them together
    famAll_4L = [fam_L; fam_m; fam_M; fam_H];
    
    % sort familiarity for the continuous arrangement
    [~,famSort] = sort(famAll);
    
    % put in output struct
    arrIdx(subs(s)).norm2fam_disc = sessIdx(famAll_4L);
    arrIdx(subs(s)).norm2fam_cont = sessIdx(famSort);
    
end