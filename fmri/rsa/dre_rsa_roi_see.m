%% dre_rsa_roi_see
% ~~~
% GX Castegnetti --- start ~ 17.07.18 --- last ~ 03.08.18

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_roi_box';

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

%% load masks
roiNames = {'HPC','ACC','infFG','medFG','midFG','supFG'};

%% subjects
subs = [5 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39];
taskOrd = [ones(1,8),2*ones(1,11),1,2,1];

%% some options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

%% load response patters computed in dre_rsa_roi_run
load([dir.out,fs,analysisName,fs,'rsaPatterns_roi.mat'],'responsePatterns')

%% construct RDMs
RDMs_data = constructRDMs(responsePatterns, 'SPM', userOptions);
RDM_average = averageRDMs_subjectSession(RDMs_data,'subject');

%% plot RDMs
% matrices
figureRDMs(RDM_average,userOptions)

% dendrograms
% dendrogramConditions(RDM_average,userOptions)

% MDS
MDSConditions(RDM_average,userOptions)

%% extract models of value, confidence, familiarity, price
RDMs_models = dre_extractRDMs(dir,subs,taskOrd);

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
end

for m = 1:size(RDMs_data,1)
    for s = 1:length(subs)
        RDMs = concatenateRDMs(RDMs_data(m,s),RDMs_val{s});
        corrMat = RDMCorrMat(RDMs, 1);
        c(m,s) = corrMat(1,2);
    end
    [h,p,~,~] = ttest(c(m,:))
end
