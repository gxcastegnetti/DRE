%% dre_rsa_roi
% ~~~
% GX Castegnetti --- start ~ 17.07.18 --- last ~ 03.08.18

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_roi_1';

%% folders
dir.root = pwd;
fs       = filesep;
idcs     = strfind(dir.root,'/');
dir.dre  = dir.root(1:idcs(end-2)-1);
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.msk  = [dir.dre,fs,'out',fs,'fmri',fs,'masks'];
dir.beh  = [dir.dre,fs,'data',fs,'behaviour'];
dir.out  = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'roi'];
addpath([dir.root,fs,'routines'])
addpath([dir.dre,fs,'codes',fs,'fmri',fs,'stats',fs,'routines'])
addpath(genpath([dir.root,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))

%% load masks
roiNames = {'ba8'};

%% subjects
subs = [4:5 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39];
taskOrd = [ones(1,9),2*ones(1,11),1,2,1];

%% reverse normalise mask to subjective space and coregister
for i = 1:length(roiNames)
    
    for s = 1:length(subs)
        
        % create folder for subjective masks
        dirSub = [dir.msk,fs,roiNames{i},fs,'SF',num2str(subs(s),'%03d')];
        mkdir(dirSub)
        cd(dirSub)
        
        % copy atlas mask to subject's folder
        fileSource = [dir.msk,fs,'atlas',fs,roiNames{i},'.nii'];
        copyfile(fileSource)
        
        % select inverse deformation images from T1 segmentation step
        dirStruct = [dir.data,fs,'SF',num2str(subs(s),'%03d'),fs,'struct'];
        d = spm_select('List', dirStruct, '^iy_.*\.nii$');
        y_file = {[dirStruct fs d]};
        msk_file = {[dirSub,fs,roiNames{i},'.nii']};
        job1{1}.spatial{1}.normalise{1}.write.subj.def = y_file;
        job1{1}.spatial{1}.normalise{1}.write.subj.resample = msk_file;
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
        job2{1}.spm.spatial.coreg.write.source = {[dirSub,fs,'w',roiNames{i},'.nii']};
        job2{1}.spm.spatial.coreg.write.roptions.interp = 4;
        job2{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
        job2{1}.spm.spatial.coreg.write.roptions.mask = 0;
        job2{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
        
        % run job
        disp(['Normalising sub#', num2str(subs(s),'%03d')])
        d = spm_jobman('run',job2);
        clear job2 d
        
    end
end

%% 1st level
if true
    for i = 1:length(roiNames)
        nameBeta = ['_level1_',roiNames{i}];
        bData = dre_extractData(dir,subs,taskOrd,0);
        dre_level1_rsa(dir,nameBeta,subs,bData,'ba8');
    end
end

%% load betas
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

for i = 1:length(roiNames)
    dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,nameBeta];
    userOptions.betaPath = [dir.beta,filesep,'[[subjectName]]',filesep,'[[betaIdentifier]]'];
    [fullBrainVols, ~] = rsa.fmri.fMRIDataPreparation('SPM', userOptions);
    responsePatterns.(roiNames{i}) = fullBrainVols; clear fullBrainVols
end

%% construct RDMs
RDMs = rsa.constructRDMs(responsePatterns, 'SPM', userOptions);
RDM_average = rsa.rdm.averageRDMs_subjectSession(RDMs,'subject');

%% plot RDMs
% matrices
rsa.figureRDMs(RDM_average,userOptions)

% dendrograms
rsa.dendrogramConditions(RDM_average,userOptions)

%% extract scores given on day 1
% behavioural data
bData = dre_extractData(dir,subs,taskOrd,0);

for s = 1:length(subs)
    for r = 1:4
        % find out whether this session is F or B
        sessType = bData(subs(s)).sessType{r};
        
        % extract value, confidence, familiarity and put them in a vector
        valSess = bData(subs(s)).imagination(r).(sessType).value;
        conSess = bData(subs(s)).imagination(r).(sessType).confidence;
        famSess = bData(subs(s)).imagination(r).(sessType).familiarity;
    end
end

goal = [zeros(120),ones(120); ones(120),zeros(120)];

