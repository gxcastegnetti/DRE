%% dre_dim
% ~~~
% analyses dimensionality in brain activity
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'dim_sl_box';

%% folders
fs      = filesep;
dir.dim = pwd;
idcs    = strfind(dir.dim,'/');
dir.dre = dir.dim(1:idcs(end-2)-1);
dir.sta = [dir.dre,fs,'codes',fs,'fmri',fs,'stats'];
dir.rsa = [dir.dre,fs,'codes',fs,'fmri',fs,'rsa'];
dir.msk = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'atlas'];
dir.out = [dir.dre,fs,'out',fs,'fmri',fs,'dim',fs,'sl'];
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
addpath(genpath([dir.rsa,fs,'rsatoolbox']))
addpath(genpath('functionalDimensionality'))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% load betas from RSA folder (because they are nicely extracted)
% in this initial analysis we use normalised and smoothed data, just to see
% if the toolbox works
load([dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'sl',fs,'rsa_sl_box_sw',fs,'rsaPatterns_sl.mat'],'responsePatterns')
subNames = fieldnames(responsePatterns);
for s = 1:length(subNames)
    % for every subject, this should be nVox x nCond x nRuns
    respPatt_2goals{s}(:,:,1) = responsePatterns.(subNames{s})(:,1:60);
    respPatt_2goals{s}(:,:,2) = responsePatterns.(subNames{s})(:,61:120);
    respPatt_2goals{s}(:,:,3) = responsePatterns.(subNames{s})(:,121:180);
    respPatt_2goals{s}(:,:,4) = responsePatterns.(subNames{s})(:,181:240); 
end

%% load gm mask
binaryMask = [dir.msk,fs,'rgm.nii'];

%% subjects
[bestn,r_outer, r_alter, test_tfce] = functional_dimensionality(respPatt_2goals, binaryMask, 'sphere',12);