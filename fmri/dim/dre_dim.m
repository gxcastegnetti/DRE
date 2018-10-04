%% dre_dim
% ~~~
% analyses dimensionality in brain activity
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'dim_test_1';

%% folders
fs      = filesep;
dir.dim = pwd;
idcs    = strfind(dir.dim,'/');
dir.dre = dir.dim(1:idcs(end-2)-1);
dir.sta = [dir.dre,fs,'codes',fs,'fmri',fs,'uni'];
dir.rsa = [dir.dre,fs,'codes',fs,'fmri',fs,'rsa'];
dir.beh  = [dir.dre,fs,'data',fs,'behaviour'];
dir.msk = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'atlas'];
dir.out = [dir.dre,fs,'out',fs,'fmri',fs,'dim',fs,'sl'];
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
addpath(genpath([dir.rsa,fs,'rsatoolbox']))
addpath(genpath([dir.sta,fs,'routines']))
addpath(genpath('functionalDimensionality'))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% subjects
subs = [5 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39 40 41 43 47:49];
taskOrd = [ones(1,8),2*ones(1,11),1,2,ones(1,4),2*ones(1,3)];

%% behaviour
bData = dre_extractData(dir,subs,taskOrd,0);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% native space masks and L1 are created in dre_rsa_roi_run %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load betas from RSA folder (because they are nicely extracted)
% in this initial analysis we use normalised and smoothed data, just to see
% if the toolbox works
load([dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'sl',fs,'rsa_sl_box_sw',fs,'rsaPatterns_sl.mat'],'responsePatterns')

subNames = fieldnames(responsePatterns);
for s = 1:length(subNames)
    
    % reorder objects according to presentation
    sessType1 = bData(subs(s)).sessType{1};
    idx_s1 = bData(subs(s)).imagination(1).(sessType1).objIdx;
    if strcmp(sessType1,'boat')
        idx_s1 = idx_s1 + 120;
    end
    
    sessType2 = bData(subs(s)).sessType{2};
    idx_s2 = bData(subs(s)).imagination(2).(sessType2).objIdx;
    if strcmp(sessType2,'boat')
        idx_s2 = idx_s2 + 120;
    end
    
    sessType3 = bData(subs(s)).sessType{3};
    idx_s3 = bData(subs(s)).imagination(3).(sessType3).objIdx;
    if strcmp(sessType3,'boat')
        idx_s3 = idx_s3 + 120;
    end
    
    sessType4 = bData(subs(s)).sessType{4};
    idx_s4 = bData(subs(s)).imagination(4).(sessType4).objIdx;
    if strcmp(sessType4,'boat')
        idx_s4 = idx_s4 + 120;
    end
    
    % for every subject, this should be nVox x nCond x nRuns
    respPatt_presentOrder{s}(:,:,1) = responsePatterns.(subNames{s})(:,idx_s1);
    respPatt_presentOrder{s}(:,:,2) = responsePatterns.(subNames{s})(:,idx_s2);
    respPatt_presentOrder{s}(:,:,3) = responsePatterns.(subNames{s})(:,idx_s3);
    respPatt_presentOrder{s}(:,:,4) = responsePatterns.(subNames{s})(:,idx_s4);
end

%% load gm mask
binaryMask = [dir.msk,fs,'rgm.nii'];

%% subjects
[bestn,r_outer, r_alter, test_tfce] = functional_dimensionality(respPatt_presentOrder, binaryMask, 'sphere',9);