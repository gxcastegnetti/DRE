%% dre_rsa
% ~~~
% GX Castegnetti --- start ~ 17.08.18 --- last ~ 20.08.18

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_sl_box';

%% Folders
dir.root = pwd;
fs       = filesep;
idcs     = strfind(dir.root,'/');
dir.dre  = dir.root(1:idcs(end-2)-1);
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.msk  = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'atlas'];
dir.beh  = [dir.dre,fs,'data',fs,'behaviour'];
dir.out  = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
dir.spm  = '/Users/gcastegnetti/Desktop/tools/matlab/spm12';
addpath([dir.root,fs,'routines'])
addpath(genpath([dir.root,fs,'rsatoolbox']))
addpath(genpath([dir.spm]))

%% Subjects
subs = [5 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39];

%% Set options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

%% mask <- should be 79x95x79 because we put correlation map in MNI

% reslice mask to the size of the correlation maps, if needed
if ~exist([dir.msk,fs,'rgm.nii'],'file')
    
    % load sample EPI from current subject for coregistration
    dirFun = [dir.data,fs,'SF039',fs,'fun',fs,'S4'];
    d = spm_select('List', dirFun, '^wuaf.*\.nii$');
    d = d(end-1,:);
    epi_file = [dirFun fs d];
    job{1}.spm.spatial.coreg.write.ref = {epi_file};
    
    % select mask to coregister
    job{1}.spm.spatial.coreg.write.source = {[dir.msk,fs,'gm.nii']};
    
    % defaults
    job{1}.spm.spatial.coreg.write.roptions.interp = 4;
    job{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    job{1}.spm.spatial.coreg.write.roptions.mask = 0;
    job{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    
    % run job
    disp('Coregistering grey matter mask')
    spm_jobman('run',job);
    clear job d epi_file
end

mask = niftiread([dir.msk,fs,'rgm.nii']);
mask = logical(mask);

%% what are we looking at?
scorP = 'gol';

%% load correlation maps and make them SPM-like
dirSl = [userOptions.rootPath,filesep,'sl',fs,analysisName,fs,scorP];
for s = 1:length(subs)
    
    % load correlations and save as nifti
    load([dirSl,fs,'sl_',scorP,'_SF',num2str(subs(s),'%03d'),'.mat']);
    niftiwrite(rs,[dirSl,fs,'rs_SF',num2str(subs(s),'%03d')])
    
    % read nifti and modify matrix
    V = spm_vol([dirSl,fs,'rs_SF',num2str(subs(s),'%03d'),'.nii']);
    mri
    % load sample EPI from current subject
    dirFun = [dir.data,fs,'SF',num2str(subs(s),'%03d'),fs,'fun',fs,'S4'];
    d = spm_select('List', dirFun, '^uaf.*\.nii$');
    d = d(end-1,:);
    epi_file = {[dirFun fs d]};
    V_epi = spm_vol(epi_file);
    
    % change header
    V.mat = V_epi{1}.mat;
    V.descrip = V_epi{1}.descrip;
    
    % save modified nifti
    spm_write_vol(V,rs);
    
    fprintf('Loading correlation maps for sub#%d \n',subs(s));
end

%% normalise to MNI
for s = 1:length(subs)
    % select forward deformation images from T1 segmentation step
    dirStruct = [dir.data,fs,'SF',num2str(subs(s),'%03d'),fs,'struct'];
    d = spm_select('List', dirStruct, '^y_.*\.nii$');
    y_file = {[dirStruct fs d]};
    s_file = {[dirSl,fs,'rs_SF',num2str(subs(s),'%03d'),'.nii']};
    job{1}.spatial{1}.normalise{1}.write.subj.def = y_file;
    job{1}.spatial{1}.normalise{1}.write.subj.resample = s_file;
    job{1}.spatial{1}.normalise{1}.write.woptions.bb = [-78 -112 -70; 78 76 85]; % bounding box of volume
    job{1}.spatial{1}.normalise{1}.write.woptions.vox = [2 2 2]; % voxel size of normalised images; DEFAULT = 2x2x2 (CHANGED to acquisition resolution)
    job{1}.spatial{1}.normalise{1}.write.woptions.interp = 4; % changed default to 7th degree B-spline
    
    % run job
    disp(['Normalising sub#', num2str(subs(s),'%03d')])
    d = spm_jobman('run',job);
    clear job
end

%% smooth
for s = 1:length(subs)
    P = [dirSl,fs,'wrs_SF',num2str(subs(s),'%03d'),'.nii'];
    Q = [dirSl,fs,'swrs_SF',num2str(subs(s),'%03d'),'.nii'];
    spm_smooth(P,Q,[6 6 6]);
end

%% concatenate across subjects
for s = 1:length(subs)
    thisRs = niftiread([dirSl,fs,'swrs_SF',num2str(subs(s),'%03d'),'.nii']);
    rMaps(:,:,:,s) = thisRs(:,:,:); %#ok<*SAGROW>
    %     figure,imagesc(squeeze(rMaps(:,:,40,s)))
end

% obtain a pMaps from applying a 1-sided signrank test and also t-test to
rMaps(isnan(rMaps)) = 0;
% the model similarities:
for x=1:size(rMaps,1)
    for y=1:size(rMaps,2)
        for z=1:size(rMaps,3)
            if mask(x,y,z) == 1
                [h p1(x,y,z)] = ttest(squeeze(rMaps(x,y,z,:)),0,0.05,'right');
                [p2(x,y,z)] = signrank_onesided(squeeze(rMaps(x,y,z,:)));
            else
                p1(x,y,z) = NaN;
                p2(x,y,z) = NaN;
            end
        end
    end
    disp(x);
end

% apply FDR correction
pThrsh_t  = FDRthreshold(p1,0.05,mask);
pThrsh_sr = FDRthreshold(p2,0.05,mask);

for i = 1:79
    figure,imagesc(p1(:,:,i));
end

% % mark the suprathreshold voxels in yellow
supraThreshMarked_t = zeros(size(p1));
supraThreshMarked_t(p1 <= pThrsh_t) = 1;
supraThreshMarked_sr = zeros(size(p2));
supraThreshMarked_sr(p2 <= pThrsh_sr) = 1;

niftiwrite(supraThreshMarked_sr,[dirSl,fs,'signVox_sr_',scorP])
aaa = niftiread([dirSl,fs,'signVox_sr_',scorP]);

%% coregister significance map

% read nifti and modify matrix
V = spm_vol([dirSl,fs,'signVox_sr_',scorP,'.nii']);

% load reference mask
V_mask = spm_vol([dir.msk,fs,'rgm.nii']);

% change header
V.mat = V_mask.mat;
V.descrip = V_mask.descrip;

% save modified nifti
spm_write_vol(V,supraThreshMarked_sr);

