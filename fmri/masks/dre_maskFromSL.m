%% dre_maskFromSL
% ~~~
% GX Castegnetti --- 2018

clear
close all

%% analysis and contrast name
analysisName = 'rsa_sl_choice';
contrastName = 'chos';

%% directories
dir.mskCod = pwd;
fs         = filesep;
idcs       = strfind(dir.mskCod,'/');
dir.dre    = dir.mskCod(1:idcs(end-2)-1); clear idcs
dir.datScn = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks'];
dir.rsaOut = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'sl',fs,analysisName];

% paths
addpath([dir.dre,fs,'codes',fs,'fmri',fs,'uni',fs,'routines'])
addpath(genpath([dir.mskCod,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))

%% load clusters and write them into separate masks
dir.cluOut = [dir.rsaOut,fs,contrastName,fs,'cluster'];
cluFile = [dir.cluOut,fs,'_SnPM_filtered.nii'];

% load in spm
cluStruct = spm_vol(cluFile);
cluMatrix = spm_read_vols(cluStruct);

% find connected clusters
k = 1;
cluster = {};
L = bwlabeln(cluMatrix > 0);

% plot to see if things are right
% for i = 1:79
%     figure,imagesc(L(:,:,i))
% end

%% separate clusters
for i = 1:max(L(:))
    cluster_foo{i} = L == i; %#ok<*SAGROW>
    clusterSize(i) = sum(cluster_foo{i}(:));
    
    % remove very small clusters
    if clusterSize(i) >= 10
        cluster{k} = cluster_foo{i};
        
        % update index
        k = k+1;
    end
end, clear i k L

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% separate hippocampus from lingual gyrus %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reference normalised EPI
sw_epi = '/Users/gcastegnetti/Desktop/stds/DRE/data/fmri/scanner/SF050/fun/S4/swuafMQ04886-0010-00225-000225-01.nii';

% copy atlas mask to subject's folder
maskHpc_l = [dir.mskOut,fs,'_useNow',fs,'l_hpc.nii'];
maskHpc_r = [dir.mskOut,fs,'_useNow',fs,'r_hpc.nii'];

%% first coregister left HPC...
jobL{1}.spm.spatial.coreg.write.ref = {sw_epi};
jobL{1}.spm.spatial.coreg.write.source = {maskHpc_l};

% defaults
jobL{1}.spm.spatial.coreg.write.roptions.interp = 4;
jobL{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
jobL{1}.spm.spatial.coreg.write.roptions.mask = 0;
jobL{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

% run job
disp('Coregistering left HPC')
spm_jobman('run',jobL);
clear jobL

%% then right HPC...
jobR{1}.spm.spatial.coreg.write.ref = {sw_epi};
jobR{1}.spm.spatial.coreg.write.source = {maskHpc_r};

% defaults
jobR{1}.spm.spatial.coreg.write.roptions.interp = 4;
jobR{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
jobR{1}.spm.spatial.coreg.write.roptions.mask = 0;
jobR{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

% run job
disp('Coregistering left HPC')
spm_jobman('run',jobR);
clear jobR

%% then reload the masks back in
% se�ect coregistered masks
rmaskHpc_l = [dir.mskOut,fs,'_useNow',fs,'rl_hpc.nii'];
rmaskHpc_r = [dir.mskOut,fs,'_useNow',fs,'rr_hpc.nii'];
maskHpc_l_mat = spm_read_vols(spm_vol(rmaskHpc_l));
maskHpc_r_mat = spm_read_vols(spm_vol(rmaskHpc_r));

%% now set to zero HPC voxels in lingual mask and only retain HPC voxels in HPC mask
for i = 1:numel(cluster)
    
    % overlay clusters with HPC cluster
    clusMat = cluster{i};
    clusMat_lHpc = clusMat.*maskHpc_l_mat;
    clusMat_rHpc = clusMat.*maskHpc_r_mat;
    
    % consider only those in which there is something in HPC
    thereIsSmth_l = sum(clusMat_lHpc(:)) > 10;
    thereIsSmth_r = sum(clusMat_rHpc(:)) > 10;
    
    % if there is something, separate hippocampal cluster and append it at
    % the end of the cluster cell array
    if thereIsSmth_l
        cluster{end+1} = clusMat_lHpc;
    elseif thereIsSmth_r
        cluster{end+1} = clusMat_rHpc;
    end
    
    % now update the current cluster to exclude HPC
    cluster{i} = clusMat.*(1-maskHpc_l_mat).*(1-maskHpc_r_mat);
    
end, clear i

%% write masks
for i = 1:numel(cluster)
    maskMetadataStruct_sS = spm_vol(sw_epi);
    maskMetadataStruct_sS.fname = [dir.mskOut,fs,'fromRSA',fs,contrastName,fs,'mask_sl_',contrastName,'_',num2str(i),'.nii'];
    maskMetadataStruct_sS.dim = size(cluMatrix);
    spm_write_vol(maskMetadataStruct_sS, cluster{i});
end