%% dre_rsa
% ~~~
% GX Castegnetti --- start ~ 17.08.18 --- last ~ 20.08.18

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_pulse';

%% Folders
dir.root = pwd;
fs       = filesep;
idcs     = strfind(dir.root,'/');
dir.dre  = dir.root(1:idcs(end-2)-1);
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.beh  = [dir.dre,fs,'data',fs,'behaviour'];
dir.out  = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
addpath([dir.root,fs,'routines'])
addpath(genpath([dir.root,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))

%% Subjects
subs = [5 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39];

%% Set options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

%% mask
mask = niftiread([dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'atlas',fs,'rgm.nii']);
mask = logical(mask);

%% prepare correlation volumes
dirSl = [userOptions.rootPath,filesep,'sl',fs,analysisName];
for s = 1:length(subs)
    
    % load correlations and save as nifti
    load([dirSl,fs,'sl_val_SF',num2str(subs(s),'%03d'),'.mat']);
    niftiwrite(rs,[dirSl,fs,'rs_SF',num2str(subs(s),'%03d')])
    
    % read nifti and modify matrix
    V = spm_vol([dirSl,fs,'rs_SF',num2str(subs(s),'%03d'),'.nii']);
    
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
    spm_smooth(P,Q,[3 3 3]);
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

% % mark the suprathreshold voxels in yellow
supraThreshMarked_t = zeros(size(p1));
supraThreshMarked_t(p1 <= pThrsh_t) = 1;
supraThreshMarked_sr = zeros(size(p2));
supraThreshMarked_sr(p2 <= pThrsh_sr) = 1;

for i = 1:79
    figure
    subplot(1,2,1)
    imagesc(squeeze(p1(:,:,i)))
    subplot(1,2,2)
    imagesc(squeeze(supraThreshMarked_sr(:,:,i)))
    %     figure,imagesc(squeeze(thisRs(:,:,i)))
end
