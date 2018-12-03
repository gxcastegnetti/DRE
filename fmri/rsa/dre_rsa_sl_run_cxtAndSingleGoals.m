%% dre_rsa_sl_run_cxtAndSingleGoals
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_sl_pulse_ons0_cxtAndSingleGoals';
betaid       = 'rsa_pulse_ons0';

%% directories
fs         = filesep;
dir.rsaCod = pwd;
idcs       = strfind(dir.rsaCod,'/');
dir.dre    = dir.rsaCod(1:idcs(end-2)-1); clear idcs
dir.uniCod = [dir.dre,fs,'codes',fs,'fmri',fs,'uni'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'gm_subj'];
dir.behDat = [dir.dre,fs,'data',fs,'behaviour'];
dir.out    = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'sl'];
dir.rsaOut = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
dir.datScn = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];

% paths
addpath([dir.rsaCod,fs,'routines'])
addpath([dir.uniCod,fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% subjects
subs = [4 5 7 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39 40 41 43 47:50];
taskOrd = [ones(1,10),2*ones(1,11),1,2,ones(1,4),2*ones(1,3) 1];

%% user options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';
userOptions.overwriteflag = 'r';

%% load betas
dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'level1',fs,betaid,fs,'none'];
userOptions.betaPath = [dir.beta,filesep,'[[subjectName]]',filesep,'[[betaIdentifier]]'];
filePatterns = [dir.out,fs,'_responsePatterns',fs,betaid,fs,'rsaPatterns_sl.mat'];
if ~exist(filePatterns,'file')
    responsePatterns = fMRIDataPreparation('SPM', userOptions);
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

%% run separate searchlights for the two goals

% subNames = fieldnames(responsePatterns);
% for s = 1:length(subNames)
%     
%     % update user
%     disp(['Computing correlation (separate goals) for sub#',num2str(s),' of ',num2str(length(subs))])
%     
%     % define (subjective) mask
%     binaryMask = niftiread([dir.mskOut,fs,'gm_SF',num2str(subs(s),'%03d'),'.nii']);
%     binaryMask = logical(binaryMask);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % take two reduced RDMs for the two goals %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     RDM_F{s} = RDMs{s}.val(1:120,1:120);
%     RDM_B{s} = RDMs{s}.val(121:240,121:240);
%     
%     % run searchlight
%     rs_F = searchlightGoal(responsePatterns.(subNames(s))(:,1:120), binaryMask, userOptions, searchlightOptions);
%     save([dir.out,fs,analysisName,fs,'onlyFire',fs,'sl__SF',num2str(subs(ss),'%03d')],'rs_F')
%     
% end

%% run searchlight for goal

% take behavioural data for reordering response patterns
bData = dre_extractData(dir,subs,taskOrd,0);

subNames = fieldnames(responsePatterns);
for s = 1:length(subNames)
    
    % update user
    disp(['Computing correlation (goal) for sub#',num2str(s),' of ',num2str(length(subs))])
    
    % define (subjective) mask
    binaryMask = niftiread([dir.mskOut,fs,'gm_SF',num2str(subs(s),'%03d'),'.nii']);
    binaryMask = logical(binaryMask);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reorder response patterns according to presentation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % S1
    sessType1 = bData(subs(s)).sessType{1};
    idx_s1 = bData(subs(s)).imagination(1).objIdx;
    if strcmp(sessType1,'B')
        idx_s1 = idx_s1 + 120;
    end
    
    % S2
    sessType2 = bData(subs(s)).sessType{2};
    idx_s2 = bData(subs(s)).imagination(2).objIdx;
    if strcmp(sessType2,'B')
        idx_s2 = idx_s2 + 120;
    end
    
    % S3
    sessType3 = bData(subs(s)).sessType{3};
    idx_s3 = bData(subs(s)).imagination(3).objIdx;
    if strcmp(sessType3,'B')
        idx_s3 = idx_s3 + 120;
    end
    
    % S4
    sessType4 = bData(subs(s)).sessType{4};
    idx_s4 = bData(subs(s)).imagination(4).objIdx;
    if strcmp(sessType4,'B')
        idx_s4 = idx_s4 + 120;
    end
    
    % for every subject, this should be nVox x nCond x nRuns
    respPatt_presentOrder{s}(:,:,1) = responsePatterns.(subNames{s})(:,idx_s1); %#ok<*SAGROW>
    respPatt_presentOrder{s}(:,:,2) = responsePatterns.(subNames{s})(:,idx_s2);
    respPatt_presentOrder{s}(:,:,3) = responsePatterns.(subNames{s})(:,idx_s3);
    respPatt_presentOrder{s}(:,:,4) = responsePatterns.(subNames{s})(:,idx_s4);
    
    % run searchlight
    rs = searchlightGoal(respPatt_presentOrder{s}, binaryMask, userOptions, searchlightOptions);
    save([dir.out,fs,analysisName,fs,'goal',fs,'sl_goal_SF',num2str(subs(s),'%03d')],'rs')
    
end