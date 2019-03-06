%% dre_rsa_sl_compare
% ~~~
% GX Castegnetti --- 2019

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_sl_ima';

%% Subjects
subs = [4 5 7 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39:43 47:49];

%% directories
dir.rsaCod = pwd;
fs         = filesep;
idcs       = strfind(dir.rsaCod,'/');
dir.dre    = dir.rsaCod(1:idcs(end-2)-1); clear idcs
dir.datScn = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'atlas'];
dir.dimOut = [dir.dre,fs,'out',fs,'fmri',fs,'dim'];
dir.behDat = [dir.dre,fs,'data',fs,'behaviour'];
dir.analys = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'sl',fs,analysisName];

% paths
dir.spm  = '/Users/gcastegnetti/Desktop/tools/matlab/spm12';
addpath([dir.rsaCod,fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath(dir.spm))

%% specify model(s) to test
model2 = 'valPcon';
model1 = 'val';

dir_mod2 = [dir.analys,fs,model2];
dir_mod1 = [dir.analys,fs,model1];

for s = 1:length(subs)
    
    swrMapFile_mod2 = [dir_mod2,fs,'swrMap_',model2,'_SF',num2str(subs(s),'%03d'),'.nii'];
    swrMapFile_mod1 = [dir_mod1,fs,'swrMap_',model1,'_SF',num2str(subs(s),'%03d'),'.nii'];
    
    metadataStruct_mod2 = spm_vol(swrMapFile_mod2);
    metadataStruct_mod1 = spm_vol(swrMapFile_mod1);
    
    data_mod2 = spm_read_vols(metadataStruct_mod2);
    data_mod1 = spm_read_vols(metadataStruct_mod1);
    
    data_diff(:,:,:,s) = data_mod2 - data_mod1;
    
    metadataStruct_mod2.fname = [dir_mod2,fs,'swrMap_vPcMv_SF',num2str(subs(s),'%03d'),'.nii'];
    spm_write_vol(metadataStruct_mod2, data_diff(:,:,:,s));

    % update user
    disp(['sub#',num2str(s),' of ',num2str(length(subs))])
       
end