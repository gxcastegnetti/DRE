%% dre_reprReinstatem
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'reprReinstatement_v0';

%% directories
fs         = filesep;
dir.geoCod = pwd;
idcs       = strfind(dir.geoCod,'/');
dir.dre    = dir.geoCod(1:idcs(end-2)-1); clear idcs
dir.uniCod = [dir.dre,fs,'codes',fs,'fmri',fs,'uni'];
dir.rsaCod = [dir.dre,fs,'codes',fs,'fmri',fs,'rsa'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks'];
dir.behDat = [dir.dre,fs,'data',fs,'behaviour'];
dir.out    = [dir.dre,fs,'out',fs,'fmri',fs,'class'];
dir.rsaOut = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];

% paths
addpath([pwd,fs,'routines'])
addpath([dir.uniCod,fs,'routines'])
addpath([dir.rsaCod,fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% subjects
subs = [4 5 7 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39 40 41 43 47:49 50];
taskOrd = [ones(1,10),2*ones(1,11),1,2,ones(1,4),2*ones(1,3) 1];

%% load response patterns during choice and during imagination
filePattIma = '/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/rsa/sl/_responsePatterns/rsa_pulse_ons0/rsaPatterns_sl.mat';
filePattCho = '/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/rsa/sl/_responsePatterns/rsa_pulse_choice/rsaPatterns_sl.mat';
load(filePattIma,'responsePatterns'), clear filePattIma
respPattIma_unmasked = responsePatterns; clear responsePatterns
load(filePattCho,'responsePatterns'), clear filePattCho
respPattCho_unmasked = responsePatterns; clear responsePatterns

%% ROI
roiNames = {'sphere_9--28_34_-19'};

%% apply grey matter and ROI masks
for r = 1:length(roiNames)
    for s = 1:length(subs)
        
        % mask files
        roiMaskFile = [dir.mskOut,fs,roiNames{r},'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',roiNames{r},'.nii'];
        gmMaskFile =  [dir.mskOut,fs,'gm_subj',fs,'gm_SF',num2str(subs(s),'%03d'),'.nii'];
        roiMask = spm_read_vols(spm_vol(roiMaskFile));
        gmMask = spm_read_vols(spm_vol(gmMaskFile));
        
        % vectorise it
        roiMask = reshape(roiMask, 1, []);
        gmMask = reshape(gmMask, 1, []);
        
        
        subjName = ['SF',num2str(subs(s),'%03d')]; % subject codename
        
        % apply masks
        respPattIma_foo = respPattIma_unmasked.(subjName)(logical(roiMask) & logical(gmMask),:);
        respPattIma.(['roi',num2str(r)]).(subjName) = respPattIma_foo(~isnan(respPattIma_foo(:,1)),:);
        respPattCho_foo = respPattCho_unmasked.(subjName)(logical(roiMask) & logical(gmMask),:);
        respPattCho.(['roi',num2str(r)]).(subjName) = respPattCho_foo(~isnan(respPattCho_foo(:,1)),:);
        
        clear subjName roiMaskFile gmMaskFile roiMask gmMask respPattIma_foo respPattCho_foo 
        
    end
end, clear r s respPattIma_unmasked respPattCho_unmasked

%% behavioural data
bData = dre_extractData(dir,subs,taskOrd,0);

% objID during choice
% take response pattern of the corresponding objects during imagination