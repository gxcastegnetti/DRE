%% dre_rsa_roi_see
% ~~~
% GX Castegnetti --- start ~ 17.07.18 --- last ~ 03.08.18

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_roi_test_1';

%% folders
dir.root = pwd;
fs       = filesep;
idcs     = strfind(dir.root,'/');
dir.dre  = dir.root(1:idcs(end-2)-1);
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.msk  = [dir.dre,fs,'out',fs,'fmri',fs,'masks'];
dir.beh  = [dir.dre,fs,'data',fs,'behaviour'];
dir.out  = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'roi'];
addpath([dir.root,fs,'routines'])
addpath([dir.dre,fs,'codes',fs,'fmri',fs,'stats',fs,'routines'])
addpath(genpath([dir.root,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))

%% subjects
subs = [5 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39];
taskOrd = [ones(1,8),2*ones(1,11),1,2,1];

%% extract behavioural data
bData = dre_extractData(dir,subs,taskOrd,0);

%% some options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

%% load response patters computed in dre_rsa_roi_run
load([dir.out,fs,analysisName,fs,'rsaPatterns_roi.mat'],'responsePatterns')

%% construct RDMs
RDMs_data = constructRDMs(responsePatterns, 'SPM', userOptions);
% RDM_average = averageRDMs_subjectSession(RDMs_data,'subject');

%% plot RDMs
% matrices
% figureRDMs(RDM_average,userOptions)

% dendrograms
% dendrogramConditions(RDM_average,userOptions)

% MDS
% MDSConditions(RDM_average,userOptions)

%% extract models of value, confidence, familiarity, price
RDMs_models = dre_extractRDMs(dir,subs,taskOrd);
mat_ID = [diag(ones(120,1)), diag(ones(120,1));
          diag(ones(120,1)), diag(ones(120,1))];
for s = 1:length(subs)
    RDMs_val{s}.name = 'value';
    RDMs_val{s}.RDM = RDMs_models{s}.val; %#ok<*SAGROW>
    RDMs_val{s}.color = [0 1 0];
    RDMs_con{s}.name = 'confidence';
    RDMs_con{s}.RDM = RDMs_models{s}.con;
    RDMs_con{s}.color = [0 1 0];
    RDMs_fam{s}.name = 'familiarity';
    RDMs_fam{s}.RDM = RDMs_models{s}.fam;
    RDMs_fam{s}.color = [0 1 0];
    RDMs_pri{s}.name = 'price';
    RDMs_pri{s}.RDM = RDMs_models{s}.pri;
    RDMs_pri{s}.color = [0 1 0];
    RDMs_oid{s}.name = 'obj ID';
    RDMs_oid{s}.RDM = mat_ID;
    RDMs_oid{s}.color = [0 1 0];
    
end

for m = 1:size(RDMs_data,1)
    for s = 1:length(subs)
        
        RDMs = concatenateRDMs(RDMs_data(m,s),RDMs_con{s});
        corrMat = RDMCorrMat(RDMs, 1);
        c(m,s) = corrMat(1,2);
        
%         RDMs = concatenateRDMs(RDMs_data(m,s),RDMs_val{s});
%         corrMat = RDMCorrMat(RDMs, 1);
%         c(m,s) = corrMat(1,2);
%         
%         RDMs = concatenateRDMs(RDMs_data(m,s),RDMs_val{s});
%         corrMat = RDMCorrMat(RDMs, 1);
%         c(m,s) = corrMat(1,2);
%         
    end
    [h,p(m),~,~] = ttest(c(m,:))
end

%% test goal
% find ROIs considered in this analysis
roiNames = fieldnames(responsePatterns);

for i = 1:length(roiNames)
    % for each subject take the vector of presented items
    for s = 1:length(subs)
        
        % find indices of items presented in each session. Add 120 to the boat
        % conditions because in responsePatters the 120 boat activity
        % patterns are stacked after the 120 fire patterns.
        if taskOrd(s) == 1
            idx_S1 = bData(subs(s)).imagination(1).fire.objIdx;
            idx_S2 = bData(subs(s)).imagination(2).boat.objIdx + 120;
            idx_S3 = bData(subs(s)).imagination(3).fire.objIdx;
            idx_S4 = bData(subs(s)).imagination(4).boat.objIdx + 120;
        else
            idx_S1 = bData(subs(s)).imagination(1).boat.objIdx + 120;
            idx_S2 = bData(subs(s)).imagination(2).fire.objIdx;
            idx_S3 = bData(subs(s)).imagination(3).boat.objIdx + 120;
            idx_S4 = bData(subs(s)).imagination(4).fire.objIdx;
        end
        
        subjName = ['SF',num2str(subs(s),'%03d')];
        
        % fetch activation patterns
        actPatt_S1 = responsePatterns.(roiNames{i}).(subjName)(:,idx_S1);
        actPatt_S2 = responsePatterns.(roiNames{i}).(subjName)(:,idx_S2);
        actPatt_S3 = responsePatterns.(roiNames{i}).(subjName)(:,idx_S3);
        actPatt_S4 = responsePatterns.(roiNames{i}).(subjName)(:,idx_S4);
        
        corr_S1S2 = corr(actPatt_S1,actPatt_S2);
        corr_S1S3 = corr(actPatt_S1,actPatt_S3);
        foo(s) = mean(mean(corr_S1S3)) - mean(mean(corr_S1S2));
        
    end
    
    roiNames{i}
    [h,p,~,~] = ttest(foo)
    
end

