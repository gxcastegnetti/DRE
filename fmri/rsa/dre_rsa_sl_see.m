%% dre_rsa
% ~~~
% GX Castegnetti --- start ~ 17.08.18 --- last ~ 20.08.18

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_sl_grey_1';

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
subs = [4:5 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39];

%% Set options
userOptions = DRE_RSA_userOptions(dir,subs);
userOptions.analysisName = 'RSA_prova';
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

%% mask
maskAnat = niftiread('/Users/gcastegnetti/Desktop/stds/DRE/data/fmri/scanner/SF008/fun/S3/wuafMQ04605-0008-00204-000204-01.nii');
mask = niftiread('/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/masks/rgrey_mni.nii');
mask = logical(mask);

%% Load previously computed rMaps and save correlation as nifti
for s = 1:length(subs)
    load([userOptions.rootPath,filesep,'searchLight',fs,analysisName,fs,'slOut_SF',num2str(subs(s),'%03d'),'.mat']);
    rMap{s} = rs; %#ok<*SAGROW>
    niftiwrite(rs,[userOptions.rootPath,filesep,'searchLight',fs,analysisName,fs,'rs_SF',num2str(subs(s),'%03d')])
    fprintf('Loading correlation maps for sub#%d \n',subs(s));
end

%% normalise to MNI
for s = 1:length(subs)
    % select forward deformation images from T1 segmentation step
    dirStruct = [dir.data,fs,'SF',num2str(subs(s),'%03d'),fs,'struct'];
    d = spm_select('List', dirStruct, '^y_.*\.nii$');
    y_file = {[dirStruct fs d]};
    s_file = {[userOptions.rootPath,filesep,'searchLight',fs,analysisName,fs,'rrs_SF',num2str(subs(s),'%03d'),'.nii']};
    job{1}.spatial{1}.normalise{1}.write.subj.def = y_file;
    job{1}.spatial{1}.normalise{1}.write.subj.resample = s_file;
    job{1}.spatial{1}.normalise{1}.write.woptions.bb = [-78 -112 -70; 78 76 85]; % bounding box of volume
    job{1}.spatial{1}.normalise{1}.write.woptions.vox = [3 3 3]; % voxel size of normalised images; DEFAULT = 2x2x2 (CHANGED to acquisition resolution)
    job{1}.spatial{1}.normalise{1}.write.woptions.interp = 4; % changed default to 7th degree B-spline
    
    % run job
    disp(['Normalising sub#', num2str(subs(s),'%03d')])
    d = spm_jobman('run',job);
    clear job
end

%% smooth
for s = 1:length(subs)
    P = [userOptions.rootPath,filesep,'searchLight',fs,analysisName,fs,'wrs_SF',num2str(subs(s),'%03d'),'.nii'];
    Q = [userOptions.rootPath,filesep,'searchLight',fs,analysisName,fs,'swrs_SF',num2str(subs(s),'%03d'),'.nii'];
    spm_smooth(P,Q,[6 6 6]);
end

%% concatenate across subjects
for s = 1:length(subs)
    thisRs = niftiread([userOptions.rootPath,filesep,'searchLight',fs,analysisName,fs,'wrs_SF',num2str(subs(s),'%03d'),'.nii']);
    rMaps(:,:,:,s) = thisRs(:,:,:);
    figure,imagesc(squeeze(rMaps(:,:,20,s))), colorbar
end

% obtain a pMaps from applying a 1-sided signrank test and also t-test to
rMaps(isnan(rMaps)) = 0;
% the model similarities:
for x=1:size(rMaps,1)
    for y=1:size(rMaps,2)
        for z=1:size(rMaps,3)
            if mask(x,y,z) == 1
                [h p1(x,y,z)] = ttest(squeeze(rMaps(x,y,z,:)),0,0.05,'right');
                [p2(x,y,z)] = rsa.util.signrank_onesided(squeeze(rMaps(x,y,z,:)));
            else
                p1(x,y,z) = NaN;
                p2(x,y,z) = NaN;
            end
        end
    end
    disp(x);
end

% apply FDR correction
% pThrsh_t  = rsa.stat.FDRthreshold(p1,0.05,mask);
% pThrsh_sr = rsa.stat.FDRthreshold(p2,0.05,mask);
%
% % mark the suprathreshold voxels in yellow
% supraThreshMarked_t = zeros(size(p1));
% supraThreshMarked_t(p1 <= pThrsh_t) = 1;
% supraThreshMarked_sr = zeros(size(p2));
% supraThreshMarked_sr(p2 <= pThrsh_sr) = 1;

for i = 1:48
    figure,imagesc(squeeze(p1(:,:,i)))
    %     figure,imagesc(squeeze(supraThreshMarked_sr(:,:,i)))
    %     figure,imagesc(squeeze(thisRs(:,:,i)))
end


% display the FDR-thresholded maps on a sample anatomy (signed rank test) :
brainVol = rsa.fmri.addRoiToVol(rsa.util.map2vol(maskAnat),rsa.util.mask2roi(mask),[1 0 0],1);
brainVol_sr = rsa.fmri.addBinaryMapToVol(brainVol,supraThreshMarked_sr.*mask,[1 1 0]);
rsa.fig.showVol(brainVol_sr,'signrank, E(FDR) < .05',3)
% rsa.fig.handleCurrentFigure([returnHere,filesep,'DEMO4',filesep,'results_DEMO4_signRank'],userOptions);

% display the FDR-thresholded maps on a sample anatomy (t-test) :
brainVol = rsa.fmri.addRoiToVol(rsa.util.map2vol(maskAnat),rsa.util.mask2roi(mask),[1 0 0],1);
brainVol_t = rsa.fmri.addBinaryMapToVol(brainVol,supraThreshMarked_t.*mask,[1 1 0]);
rsa.fig.showVol(brainVol_t,'t-test, E(FDR) < .05',4)
