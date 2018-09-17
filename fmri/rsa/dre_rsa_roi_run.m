%% dre_rsa_roi_run
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_roi_box_w';

%% folders
dir.root = pwd;
fs       = filesep;
idcs     = strfind(dir.root,'/');
dir.dre  = dir.root(1:idcs(end-2)-1);
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.msk  = [dir.dre,fs,'out',fs,'fmri',fs,'masks'];
dir.beh  = [dir.dre,fs,'data',fs,'behaviour'];
dir.out  = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'roi'];
dir.rsa  = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
addpath([dir.root,fs,'routines'])
addpath([dir.dre,fs,'codes',fs,'fmri',fs,'stats',fs,'routines'])
addpath(genpath([dir.root,fs,'rsatoolbox',fs,'rsa']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% load masks
% roiNames = {'hpc','ba8','ba9','ba10','ba11','ba25','ba32','ba33','ba34','ba44','ba45','ba46','ba47'};
% roiNames = {'gm'};
roiNames = {'HPC','ACC','infFG','medFG','midFG','supFG','infOcc','supOcc','PCC','AG','Ins','paraHPC'};

%% subjects
subs = [5 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39];
taskOrd = [ones(1,8),2*ones(1,11),1,2,1];

%% reverse normalise mask to subjective space and coregister
if false
    for i = 1:length(roiNames)
        for s = 1:length(subs)
            
            % create folder for subjective masks
            dirSub = [dir.msk,fs,roiNames{i},'_subj',fs,'SF',num2str(subs(s),'%03d')];
            if ~exist(dirSub,'dir'),mkdir(dirSub),end
            cd(dirSub)
            
            % copy atlas mask to subject's folder
            fileSource = [dir.msk,fs,'atlas',fs,roiNames{i},'.nii'];
            copyfile(fileSource)
            
            % select inverse deformation images from T1 segmentation step
            dirStruct = [dir.data,fs,'SF',num2str(subs(s),'%03d'),fs,'struct'];
            d = spm_select('List', dirStruct, '^iy_.*\.nii$');
            y_file = {[dirStruct fs d]};
            job1{1}.spatial{1}.normalise{1}.write.subj.def = y_file;
            
            % select mask
            msk_file = {[dirSub,fs,roiNames{i},'.nii']};
            job1{1}.spatial{1}.normalise{1}.write.subj.resample = msk_file;
            
            % defaults
            job1{1}.spatial{1}.normalise{1}.write.woptions.bb = [-78 -112 -70; 78 76 85]; % bounding box of volume
            job1{1}.spatial{1}.normalise{1}.write.woptions.vox = [2 2 2]; % voxel size of normalised images; DEFAULT = 2x2x2 (CHANGED to acquisition resolution)
            job1{1}.spatial{1}.normalise{1}.write.woptions.interp = 4; % changed default to 7th degree B-spline
            
            % run job
            disp(['Normalising sub#', num2str(subs(s),'%03d')])
            d = spm_jobman('run',job1);
            clear job1 d
            
            % delete unwarped file
            delete([dirSub,fs,roiNames{i},'.nii'])
            
            % load sample EPI from current subject for coregistration
            dirFun = [dir.data,fs,'SF',num2str(subs(s),'%03d'),fs,'fun',fs,'S4'];
            d = spm_select('List', dirFun, '^uaf.*\.nii$');
            d = d(end-1,:);
            epi_file = [dirFun fs d];
            job2{1}.spm.spatial.coreg.write.ref = {epi_file};
            
            % select mask to coregister 
            job2{1}.spm.spatial.coreg.write.source = {[dirSub,fs,'w',roiNames{i},'.nii']};
            
            % defaults
            job2{1}.spm.spatial.coreg.write.roptions.interp = 4;
            job2{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            job2{1}.spm.spatial.coreg.write.roptions.mask = 0;
            job2{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
            
            % run job
            disp(['Coregistering sub#', num2str(subs(s),'%03d')])
            d = spm_jobman('run',job2);
            clear job2 d
            
        end
    end
end

%% 1st level
roiNames = {'none'};
if true
    for i = 1:length(roiNames)
        nameBeta = ['level1',fs,'rsa_box',fs,roiNames{i}];
        bData = dre_extractData(dir,subs,taskOrd,0);
        timing.iOns = 0;
        timing.iDur = 5;
        dre_level1_rsa(dir,nameBeta,subs,bData,timing,roiNames{i});
    end
end

%% load betas and build response patterns
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

for i = 1:length(roiNames)
    nameBeta = ['level1',fs,'rsa_box',fs,roiNames{i}];
    dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,nameBeta];
    userOptions.betaPath = [dir.beta,filesep,'[[subjectName]]',filesep,'[[betaIdentifier]]'];
    [fullBrainVols, ~] = fMRIDataPreparation('SPM', userOptions);
    responsePatterns.(roiNames{i}) = fullBrainVols; clear fullBrainVols
end
save([dir.out,fs,analysisName,fs,'rsaPatterns_roi'],'responsePatterns','-v7.3')

