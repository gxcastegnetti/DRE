%% dre_rsa_geom
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'geom_F';

%% directories
fs         = filesep;
dir.geoCod = pwd;
idcs       = strfind(dir.geoCod,'/');
dir.dre    = dir.geoCod(1:idcs(end-2)-1); clear idcs
dir.uniCod = [dir.dre,fs,'codes',fs,'fmri',fs,'uni'];
dir.rsaCod = [dir.dre,fs,'codes',fs,'fmri',fs,'rsa'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'gm_subj'];
dir.behDat = [dir.dre,fs,'data',fs,'behaviour'];
dir.out    = [dir.dre,fs,'out',fs,'fmri',fs,'geom',fs,'sl'];
dir.rsaOut = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];

% paths
addpath([pwd,fs,'routines'])
addpath([dir.uniCod,fs,'routines'])
addpath([dir.rsaCod,fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% subjects
subs = [4 5 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39 40 41 43 47:49];
taskOrd = [ones(1,9),2*ones(1,10),1,2,ones(1,4),2*ones(1,3)];

%% user options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';
userOptions.overwriteflag = 'r';

%% load betas
filePatterns = '/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/rsa/sl/_responsePatterns/rsa_pulse_ons0/rsaPatterns_sl.mat';
load(filePatterns,'responsePatterns')

%% searchlight options
userOptions.voxelSize = [3 3 3];
userOptions.searchlightRadius = 9;
searchlightOptions.monitor = false;
searchlightOptions.fisher = true;
searchlightOptions.nSessions = 1;
searchlightOptions.nConditions = 240;

%% extract behavioural data and rearrange
% bData = dre_extractData(dir,subs,taskOrd,0);
% ordData = dre_rearrange_4L(dir,subs,taskOrd,bData);
%
% % respPatt_acc2ses_vis = responsePatterns;
% for s = 1:length(subs)
%     subjName = ['SF',num2str(subs(s),'%03d')];
%     respPatt_acc2ses_vis.gm.(subjName) = responsePatterns.(subjName)(~isnan(responsePatterns.(subjName)(:,1)),ordData(subs(s)).norm2sessions);
%     respPatt_acc2ses{s} = responsePatterns.(subjName)(:,ordData(subs(s)).norm2sessions);
%
% end

%% construct and plot RDMs for visual check
% RDMs_data = constructRDMs(respPatt_acc2ses_vis, 'SPM', userOptions);
% RDM_average = averageRDMs_subjectSession(RDMs_data,'subject');
% figureRDMs(RDM_average,userOptions)

%% run searchlight for imagination
for s = 1:length(subs)
    
    % update user
    disp(['Computing correlation for sub#',num2str(s),' of ',num2str(length(subs))])
    
    % prepare mask
    binaryMask = niftiread([dir.mskOut,fs,'gm_SF',num2str(subs(s),'%03d'),'.nii']);
    binaryMask = logical(binaryMask);
    thisSubject = userOptions.subjectNames{s};
    
    % run searchlight
    geomDiff = searchlightMapping_geom(responsePatterns.(thisSubject), binaryMask, userOptions, searchlightOptions); %#ok<NASGU>
    save([dir.out,fs,analysisName,fs,'geom_SF',num2str(subs(s),'%03d')],'geomDiff','-v7.3')
    clear geomDiff binaryMask
end
