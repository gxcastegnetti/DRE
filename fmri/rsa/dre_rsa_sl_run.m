%% dre_rsa_sl_run
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_sl_pulse_ons-1';
betaid       = 'rsa_pulse_ons-1';

%% folders
fs      = filesep;
dir.here = pwd;
idcs    = strfind(dir.here,'/');
dir.dre = dir.here(1:idcs(end-2)-1);
dir.sta = [dir.dre,fs,'codes',fs,'fmri',fs,'uni'];
dir.msk = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'gm_subj'];
dir.beh = [dir.dre,fs,'data',fs,'behaviour'];
dir.out = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'sl'];
dir.rsa = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
addpath([dir.here,fs,'routines'])
addpath([dir.sta,fs,'routines'])
addpath(genpath([dir.here,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% subjects
subs = [4 5 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39 40 41 43 47:49];
taskOrd = [ones(1,9),2*ones(1,10),1,2,ones(1,4),2*ones(1,3)];
% subsBest = [23 18 5 3 21 11 10 20 17 28 24  1 15 22];
% subsWors = [2  14 6 4  7  9 27 26 12 16 19 13  8 25];

%% user options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';
userOptions.overwriteflag = 'r';

%% 1st level
roiNames = {'none'};
if true
    for i = 1:length(roiNames)
        nameBeta = ['level1',fs,betaid,fs,roiNames{i}];
        bData = dre_extractData(dir,subs,taskOrd,0);
        timing.iOns = -1;
        timing.iDur = 0;
        dre_level1_rsa(dir,nameBeta,subs,bData,timing,roiNames{i});
    end
end

%% load betas
dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'level1',fs,betaid,fs,'none'];
userOptions.betaPath = [dir.beta,filesep,'[[subjectName]]',filesep,'[[betaIdentifier]]'];
filePatterns = [dir.out,fs,'_responsePatterns',fs,betaid,fs,'rsaPatterns_sl.mat'];
if ~exist(filePatterns,'file')
    [~, responsePatterns] = fMRIDataPreparation('SPM', userOptions);
    save(filePatterns,'responsePatterns','-v7.3')
else
    load(filePatterns,'responsePatterns')
end

%% extract models of value, confidence, familiarity, price
RDMs = dre_extractRDMs(dir,subs,taskOrd);

%% searchlight options
userOptions.voxelSize = [3 3 3];
userOptions.searchlightRadius = 9;
searchlightOptions.monitor = false;
searchlightOptions.fisher = true;
searchlightOptions.nSessions = 1;
searchlightOptions.nConditions = 240;

%% run searchlight for conitnuous value, confidence, familiarity, price, object ID
% create matrix for object ID
mat_ID = [diag(ones(120,1)), diag(ones(120,1));
          diag(ones(120,1)), diag(ones(120,1))];
for s = 1:length(subs)
    disp(['Computing correlation for sub#',num2str(s),' of ',num2str(length(subs))])
    binaryMask = niftiread([dir.msk,fs,'gm_SF',num2str(subs(s),'%03d'),'.nii']);
    binaryMask = logical(binaryMask);
    thisSubject = userOptions.subjectNames{s};
    model(1).name = 'val';
    model(1).RDM = RDMs{s}.val;
    model(1).color = [0 1 0];
    model(2).name = 'fam';
    model(2).RDM = RDMs{s}.fam;
    model(2).color = [0 1 0];
    model(3).name = 'oid';
    model(3).RDM = 1-mat_ID;
    model(3).color = [0 1 0];
    model(4).name = 'cxt';
    model(4).RDM = RDMs{s}.cxt;
    model(4).color = [0 1 0];
    [rs,~,~,~] = searchlightMapping_fMRI(responsePatterns.(thisSubject), model, binaryMask, userOptions, searchlightOptions); %#ok<*ASGLU>
    save([dir.out,fs,analysisName,fs,'sl_SF',num2str(subs(s),'%03d')],'rs','model')
    clear model rs binaryMask
end

%% run searchlight for goal
% 
% % take behavioural data for reordering response patterns
% bData = dre_extractData(dir,subs,taskOrd,0);
% 
% subNames = fieldnames(responsePatterns);
% subNames(25) = [];
% for s = 1:length(subNames)
%     
%     % update user
%     disp(['Computing correlation (context) for sub#',num2str(s),' of ',num2str(length(subs))])
%     
%     % define (subjective) mask
%     binaryMask = niftiread([dir.msk,fs,'gm_SF',num2str(subs(s),'%03d'),'.nii']);
%     binaryMask = logical(binaryMask);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % reorder response patterns according to presentation %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % S1
%     sessType1 = bData(subs(s)).sessType{1};
%     idx_s1 = bData(subs(s)).imagination(1).objIdx;
%     if strcmp(sessType1,'B')
%         idx_s1 = idx_s1 + 120;
%     end
%     
%     % S2
%     sessType2 = bData(subs(s)).sessType{2};
%     idx_s2 = bData(subs(s)).imagination(2).objIdx;
%     if strcmp(sessType2,'B')
%         idx_s2 = idx_s2 + 120;
%     end
%     
%     % S3
%     sessType3 = bData(subs(s)).sessType{3};
%     idx_s3 = bData(subs(s)).imagination(3).objIdx;
%     if strcmp(sessType3,'B')
%         idx_s3 = idx_s3 + 120;
%     end
%     
%     % S4
%     sessType4 = bData(subs(s)).sessType{4};
%     idx_s4 = bData(subs(s)).imagination(4).objIdx;
%     if strcmp(sessType4,'B')
%         idx_s4 = idx_s4 + 120;
%     end
%     
%     % for every subject, this should be nVox x nCond x nRuns
%     respPatt_presentOrder{s}(:,:,1) = responsePatterns.(subNames{s})(:,idx_s1); %#ok<*SAGROW>
%     respPatt_presentOrder{s}(:,:,2) = responsePatterns.(subNames{s})(:,idx_s2);
%     respPatt_presentOrder{s}(:,:,3) = responsePatterns.(subNames{s})(:,idx_s3);
%     respPatt_presentOrder{s}(:,:,4) = responsePatterns.(subNames{s})(:,idx_s4);
%     
%     % run searchlight
%     rs = searchlightGoal(respPatt_presentOrder{s}, binaryMask, userOptions, searchlightOptions);
%     save([dir.out,fs,analysisName,fs,'sl_context_SF',num2str(subs(s),'%03d')],'rs')
% end
% 
% 
