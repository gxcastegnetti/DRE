%% dre_rsa_sl_run
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_sl_pulse_ons0_pw';
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
dir.betaid = betaid;

% paths
addpath([dir.rsaCod,fs,'routines'])
addpath([dir.uniCod,fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% subjects
% subs = [4 5 7 8 9 13:17 19 20 21 23 25:26 29:32 34 35 37 39 40 41 43 47:49 50];
% taskOrd = [ones(1,10),2*ones(1,11),1,2,ones(1,4),2*ones(1,3) 1];
subs = [4 5 7 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39 40 41 43 47:50];
taskOrd = [ones(1,10),2*ones(1,11),1,2,ones(1,4),2*ones(1,3) 1];

%% user options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';
userOptions.overwriteflag = 'r';

%% 1st level
roiName = 'none';
if false
    nameBeta = ['level1',fs,betaid,fs,roiName];
    bData = dre_extractData(dir,subs,taskOrd,0);
    timing.iOns = 0;
    timing.iDur = 0;
    dre_level1_rsa(dir,nameBeta,subs,bData,timing,roiName);
end

%% load betas
dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'level1',fs,betaid,fs,roiName];
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

%% run searchlight for imagination
% create matrix for object ID
mat_ID = [diag(ones(120,1)), diag(ones(120,1));
    diag(ones(120,1)), diag(ones(120,1))];
for s = 1:length(subs)
    
    % update user
    disp(['Computing correlation for sub#',num2str(s),' of ',num2str(length(subs))])
    
    % prepare mask
    fileMask = [dir.mskOut,fs,'gm_SF',num2str(subs(s),'%03d'),'.nii'];
    thisSubject = userOptions.subjectNames{s};
    
    %% models
    
    % value
    model(1).name = 'val';
    model(1).RDM = RDMs{s}.val;
    model(1).color = [0 1 0];
    
    % familiarity
    model(2).name = 'fam';
    model(2).RDM = RDMs{s}.fam;
    model(2).color = [0 1 0];
    
    % object identity
    model(3).name = 'oid';
    model(3).RDM = 1-mat_ID;
    model(3).color = [0 1 0];
    
    % context
    model(4).name = 'cxt';
    model(4).RDM = RDMs{s}.cxt;
    model(4).color = [0 1 0];
    
    %% run searchlight
    model = RDMs{s}.val;
    rs = searchlight_pw(dir,subs(s),analysisName,fileMask,model); %#ok<*ASGLU>
    save([dir.out,fs,analysisName,fs,'sl_val_SF',num2str(subs(s),'%03d')],'rs','model')
    clear model rs binaryMask
end
