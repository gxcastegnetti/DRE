%% dre_dim
% ~~~
% analyses dimensionality in brain activity
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% folders
fs      = filesep;
dir.rsa = pwd;
idcs    = strfind(dir.rsa,'/');
dir.dre = dir.rsa(1:idcs(end-2)-1);
dir.sta = [dir.dre,fs,'codes',fs,'fmri',fs,'stats'];
dir.msk = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'gm_subj'];
dir.beh = [dir.dre,fs,'data',fs,'behaviour'];
dir.out = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'sl'];
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
addpath([dir.rsa,fs,'routines'])
addpath([dir.sta,fs,'routines'])
addpath(genpath([dir.rsa,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% analysisName
analysisName = 'dim_sl_box';

%% load betas
nameBeta = ['level1',fs,'rsa_box_sw',fs,'none']; % <------------------------ set here which betas to look for
dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,nameBeta];
userOptions.betaPath = [dir.beta,filesep,'[[subjectName]]',filesep,'[[betaIdentifier]]'];
if ~exist([dir.out,fs,analysisName,fs,'rsaPatterns_sl.mat'],'file')
    [~, responsePatterns] = fMRIDataPreparation('SPM', userOptions);
    save([dir.out,fs,analysisName,fs,'rsaPatterns_sl.mat'],'responsePatterns','-v7.3')
else
    load([dir.out,fs,analysisName,fs,'rsaPatterns_sl.mat'],'responsePatterns')
end

%% subjects
subs = [5 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39];
taskOrd = [ones(1,8),2*ones(1,11),1,2,1];
[bestn,r_outer, r_alter, test_tfce] = functional_dimensionality(wholebrain_all, mask, varargin);