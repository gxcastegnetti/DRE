%% dre_rsa_roi
% ~~~
% GX Castegnetti --- start ~ 17.07.18 --- last ~ 03.08.18

clear
close all
restoredefaultpath

%% folders
roiNames = {'hpc','vmpfc','parietalSup'};
dir.root = pwd;
fs       = filesep;
idcs     = strfind(dir.root,'/');
dir.dre  = dir.root(1:idcs(end-2)-1);
dir.beh  = [dir.dre,fs,'data',fs,'behaviour'];
dir.out  = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
addpath([dir.root,fs,'routines'])
addpath([dir.dre,fs,'codes',fs,'fmri',fs,'stats',fs,'routines'])
addpath(genpath([dir.root,fs,'rsatoolbox']))
addpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12')

%% load masks and map them back to subjective space
roiNames = {'hpc','vmpfc','parietal_sup','parietal_inf'};

%% subjects
subss = [4:5 7:9 12:17 19:23 25:27 29:37 39];
taskOrd = [ones(1,12),2*ones(1,14),1,1,2,1];

%% load fMRI data
userOptions = DRE_RSA_userOptions(dir,subss);
userOptions.analysisName = 'RSA_prova';
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';

% load data with different masks
for i = 1:length(roiNames)
    dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'stats',fs,'for_RSA_pulse',fs,'r',roiNames{i}];
    userOptions.betaPath = [dir.beta,filesep,'[[subjectName]]',filesep,'[[betaIdentifier]]'];
    [fullBrainVols fullBrainVols_full] = rsa.fmri.fMRIDataPreparation('SPM', userOptions);
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
bData = dre_extractData(dir,subss,taskOrd,0);

for s = 1:length(subss)
    for r = 1:4
        % find out whether this session is F or B
        sessType = bData(subss(s)).sessType{r};
        
        % extract value, confidence, familiarity and put them in a vector
        valSess = bData(subss(s)).imagination(r).(sessType).value;
        conSess = bData(subss(s)).imagination(r).(sessType).confidence;
        famSess = bData(subss(s)).imagination(r).(sessType).familiarity;
    end
end

goal = [zeros(120),ones(120); ones(120),zeros(120)];

