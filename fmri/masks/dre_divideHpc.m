%% dre_divideHpc
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
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'_useNow'];
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))

%% open hippocampal masks
hpc_l_maskFile = [dir.mskOut,fs,'l_hpc.nii'];
hpc_r_maskFile = [dir.mskOut,fs,'r_hpc.nii'];

hpc_l_maskStruct = spm_vol(hpc_l_maskFile);
hpc_r_maskStruct = spm_vol(hpc_r_maskFile);

hpc_l_maskData = spm_read_vols(hpc_l_maskStruct);
hpc_r_maskData = spm_read_vols(hpc_r_maskStruct);

% figure
% for i = 1:91
%     subplot(9,10,i)
%     imagesc(squeeze(hpc_r_maskData(35,:,:))')
% end

hpc_la_maskData = hpc_l_maskData;
hpc_lp_maskData = hpc_l_maskData;
hpc_la_maskData(:,1:51,:) = 0;
hpc_lp_maskData(:,52:end,:) = 0;

hpc_ra_maskData = hpc_r_maskData;
hpc_rp_maskData = hpc_r_maskData;
hpc_ra_maskData(:,1:51,:) = 0;
hpc_rp_maskData(:,52:end,:) = 0;

newMaskMetadataStruct_la = hpc_l_maskStruct;
newMaskMetadataStruct_la.fname = [dir.mskOut,fs,'la_hpc.nii'];
spm_write_vol(newMaskMetadataStruct_la, hpc_la_maskData);

newMaskMetadataStruct_lp = hpc_l_maskStruct;
newMaskMetadataStruct_lp.fname = [dir.mskOut,fs,'lp_hpc.nii'];
spm_write_vol(newMaskMetadataStruct_lp, hpc_lp_maskData);

newMaskMetadataStruct_ra = hpc_l_maskStruct;
newMaskMetadataStruct_ra.fname = [dir.mskOut,fs,'ra_hpc.nii'];
spm_write_vol(newMaskMetadataStruct_ra, hpc_ra_maskData);

newMaskMetadataStruct_rp = hpc_l_maskStruct;
newMaskMetadataStruct_rp.fname = [dir.mskOut,fs,'rp_hpc.nii'];
spm_write_vol(newMaskMetadataStruct_rp, hpc_rp_maskData);





