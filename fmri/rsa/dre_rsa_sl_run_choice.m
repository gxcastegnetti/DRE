%% dre_rsa_sl_run_choice
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_sl_pulse_choice';
betaid       = 'rsa_pulse_choice';

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
subs = [4 5 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39 40 41 43 47:49];
taskOrd = [ones(1,9),2*ones(1,10),1,2,ones(1,4),2*ones(1,3)];
% subsBest = sort([23 18 5 3 21 11 10 20 17 28 24  1 15 22]);
% subsWors = sort([2  14 6 4  7  9 27 26 12 16 19 13  8 25]);
% subs = subs(subsWors);
% taskOrd = taskOrd(subsWors);

%% user options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';
userOptions.overwriteflag = 'r';

%% modify condition names
objsName = 1:48;
for i = 1:length(objsName)
    objsF{i} = ['F-',num2str(objsName(i))]; %#ok<*SAGROW>
    objsB{i} = ['B-',num2str(objsName(i))];
end
condLabels = [objsF';objsB'];

% text lables which may be attached to the conditions for MDS plots
userOptions.conditionLabels = condLabels;

% colours for the conditions
userOptions.conditionColours = kron([1 0 0; 0 0 1], ones(48,1));

%% 1st level
roiNames = {'none'};
if false
    for i = 1:length(roiNames)
        nameBeta = ['level1',fs,betaid,fs,roiNames{i}];
        bData = dre_extractData(dir,subs,taskOrd,0);
        dre_level1_rsa_choice(dir,nameBeta,subs,bData,roiNames{i});
    end
end

%% load betas
dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'level1',fs,betaid,fs,'none'];
userOptions.betaPath = [dir.beta,filesep,'[[subjectName]]',filesep,'[[betaIdentifier]]'];
filePatterns = [dir.out,fs,'_responsePatterns',fs,betaid,fs,'rsaPatterns_sl.mat'];
if ~exist(filePatterns,'file')
    [~, responsePatterns] = fMRIDataPreparation('SPM', userOptions);
    if ~exist([dir.out,fs,'_responsePatterns',fs,betaid],'dir'),mkdir([dir.out,fs,'_responsePatterns',fs,betaid]),end
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
searchlightOptions.nConditions = 96;

%% run searchlight for choice
for s = 1:length(subs)
    disp(['Computing correlation for sub#',num2str(s),' of ',num2str(length(subs))])
    binaryMask = niftiread([dir.mskOut,fs,'gm_SF',num2str(subs(s),'%03d'),'.nii']);
    binaryMask = logical(binaryMask);
    thisSubject = userOptions.subjectNames{s};
    model(1).name = 'dval';
    model(1).RDM = RDMs{s}.choice.dVal;
    model(1).color = [0 1 0];
    model(2).name = 'vCho';
    model(2).RDM = RDMs{s}.choice.Chos;
    model(2).color = [0 1 0];
    model(3).name = 'vUnc';
    model(3).RDM = RDMs{s}.choice.Unch;
    model(3).color = [0 1 0];
    model(4).name = 'cMun';
    model(4).RDM = RDMs{s}.choice.cMun;
    model(4).color = [0 1 0];
    model(5).name = 'ccxt';
    model(5).RDM = RDMs{s}.choice.ccxt;
    model(5).color = [0 1 0];
    [rs,~,~,~] = searchlightMapping_fMRI(responsePatterns.(thisSubject), model, binaryMask, userOptions, searchlightOptions); %#ok<*ASGLU>
    save([dir.out,fs,analysisName,fs,'sl_SF',num2str(subs(s),'%03d')],'rs','model')
    clear model rs binaryMask
end
