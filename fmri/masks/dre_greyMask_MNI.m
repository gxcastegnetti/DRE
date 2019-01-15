%% dre_greyMask_MNI
% Creates individual grey matter masks from thresholding tissue
% probability from the segmentation.
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% folders
fs      = filesep;
dir.msk = pwd;
idcs    = strfind(dir.msk,'/');
dir.dre = dir.msk(1:idcs(end-2)-1);
dir.out = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'gm_subj'];
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'sl',fs,'_nS_masks_gm'];
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))

%% subjects
subs = [4 5 7 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39 40 41 43 47:49 50];

for s = 1:length(subs)   
    
    %% take normalised grey matter masks
    maskFile = [dir.mskOut,fs,'wnS_gm_SF',num2str(subs(s),'%03d'),'.nii'];
    sS_maskMetaData = spm_vol(maskFile);
    sS_subMasks(:,:,:,s) = spm_read_vols(sS_maskMetaData);

end

sS_avgMaskFile = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'_useNow',fs,'gmAvg.nii'];
sS_avgMask = mean(sS_subMasks,4) > 0.75;
sS_maskMetaData.fname = sS_avgMaskFile;
spm_write_vol(sS_maskMetaData, sS_avgMask);