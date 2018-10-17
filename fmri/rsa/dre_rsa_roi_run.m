%% dre_rsa_roi_run
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_roi_pulse_cxt';
dirBeta = 'rsa_pulse_ons0';

%% directories
dir.rsaCod = pwd;
fs         = filesep;
idcs       = strfind(dir.rsaCod,'/');
dir.dre    = dir.rsaCod(1:idcs(end-2)-1); clear idcs
dir.datScn = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks'];
dir.behDat = [dir.dre,fs,'data',fs,'behaviour'];
dir.rsaOut = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
dir.out    = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'roi'];

% paths
addpath([dir.rsaCod,fs,'routines'])
addpath([dir.dre,fs,'codes',fs,'fmri',fs,'uni',fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox',fs,'rsa']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% load masks
% roiNames = {'none'};
% roiNames = {'HPC','mPFC_cS_pulse','verm_iV_pulse','rANG','paraHPC','l_midFC','insula','ACC','PCC','supOcc'};
roiNames = {'ba10','ba11','ba24','ba25','ba47','hpc'};
roiNames = {'sl_cxt'};

%% subjects
subs = [4 5 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39 40 41 43 47:49];
taskOrd = [ones(1,9),2*ones(1,10),1,2,ones(1,4),2*ones(1,3)];

%% reverse normalise mask to subjective space and coregister
if false
    for i = 1:length(roiNames)
        for s = 1:length(subs)
            
            % create folder for subjective masks
            dirSubjMask = [dir.mskOut,fs,roiNames{i},'_subj',fs,'SF',num2str(subs(s),'%03d')];
            if ~exist(dirSubjMask,'dir'),mkdir(dirSubjMask),end
            cd(dirSubjMask)
            
            % copy atlas mask to subject's folder
            fileSource = [dir.mskOut,fs,'_useNow',fs,roiNames{i},'.nii'];
            copyfile(fileSource)
            
            % select inverse deformation images from T1 segmentation step
            dirStruct = [dir.datScn,fs,'SF',num2str(subs(s),'%03d'),fs,'struct'];
            d = spm_select('List', dirStruct, '^iy_.*\.nii$');
            y_file = {[dirStruct fs d]}; clear d dirStruct
            job1{1}.spatial{1}.normalise{1}.write.subj.def = y_file;
            
            % select mask
            msk_file = {[dirSubjMask,fs,roiNames{i},'.nii']};
            job1{1}.spatial{1}.normalise{1}.write.subj.resample = msk_file;
            
            % defaults
            job1{1}.spatial{1}.normalise{1}.write.woptions.bb = [-78 -112 -70; 78 76 85]; % bounding box of volume
            job1{1}.spatial{1}.normalise{1}.write.woptions.vox = [2 2 2]; % voxel size of normalised images; DEFAULT = 2x2x2 (CHANGED to acquisition resolution)
            job1{1}.spatial{1}.normalise{1}.write.woptions.interp = 4; % changed default to 7th degree B-spline
            
            % run job
            disp(['Normalising sub#', num2str(subs(s),'%03d')])
            spm_jobman('run',job1);
            clear job1 msk_file y_file
            
            % delete unwarped file
            delete([dirSubjMask,fs,roiNames{i},'.nii'])
            
            % load sample EPI from current subject for coregistration
            dirFun = [dir.datScn,fs,'SF',num2str(subs(s),'%03d'),fs,'fun',fs,'S4'];
            d = spm_select('List', dirFun, '^uaf.*\.nii$');
            d = d(end-1,:);
            epi_file = [dirFun fs d];
            job2{1}.spm.spatial.coreg.write.ref = {epi_file};
            
            % select mask to coregister 
            job2{1}.spm.spatial.coreg.write.source = {[dirSubjMask,fs,'w',roiNames{i},'.nii']};
            
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
if false
    for i = 1:length(roiNames)
        nameBeta = ['level1',fs,dirBeta,fs,roiNames{i}];
        bData = dre_extractData(dir,subs,taskOrd,0);
        timing.iOns = 0;
        timing.iDur = 0;
        dre_level1_rsa(dir,nameBeta,subs,bData,timing,roiNames{i});
    end
end

%% load betas and build response patterns
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

for i = 1:length(roiNames)
    nameBeta = ['level1',fs,dirBeta,fs,roiNames{i}];
    dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,nameBeta];
    userOptions.betaPath = [dir.beta,filesep,'[[subjectName]]',filesep,'[[betaIdentifier]]'];
    [fullBrainVols, ~] = fMRIDataPreparation('SPM', userOptions);
    responsePatterns.(roiNames{i}) = fullBrainVols; clear fullBrainVols
end
save([dir.out,fs,analysisName,fs,'rsaPatterns_roi'],'responsePatterns','-v7.3')
