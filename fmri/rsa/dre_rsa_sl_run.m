%% dre_rsa_sl_run
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_sl_box_L1_unmasked';

%% folders
fs      = filesep;
dir.rsa = pwd;
idcs    = strfind(dir.rsa,'/');
dir.dre = dir.rsa(1:idcs(end-2)-1);
dir.sta = [dir.dre,fs,'codes',fs,'fmri',fs,'stats'];
dir.msk = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'gm_subj'];
dir.beh = [dir.dre,fs,'data',fs,'behaviour'];
dir.out = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'sl'];
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
addpath([dir.rsa,fs,'routines'])
addpath([dir.sta,fs,'routines'])
addpath(genpath([dir.rsa,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))
mkdir([dir.out,fs,analysisName])

%% subjects
subs = [5 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39];
taskOrd = [ones(1,8),2*ones(1,11),1,2,1];

%% user options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';
userOptions.overwriteflag = 'r';

%% load betas
nameBeta = ['level1',fs,'rsa_box',fs,'none']; % <------------------------ set here which betas to look for
dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,nameBeta];
userOptions.betaPath = [dir.beta,filesep,'[[subjectName]]',filesep,'[[betaIdentifier]]'];
if ~exist([dir.out,fs,analysisName,fs,'rsaPatterns_sl.mat'],'file')
    [~, responsePatterns] = fMRIDataPreparation('SPM', userOptions);
    save([dir.out,fs,analysisName,fs,'rsaPatterns_sl.mat'],'responsePatterns','-v7.3')
else
    load([dir.out,fs,analysisName,fs,'rsaPatterns_sl.mat'],'responsePatterns')
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

%% run searchlight for value, confidence, familiarity, price
for s = 1:length(subs)
    disp(['Computing correlation for sub#',num2str(s),' of ',num2str(length(subs))])
    binaryMask = niftiread([dir.msk,fs,'gm_SF',num2str(subs(s),'%03d'),'.nii']);
    binaryMask = logical(binaryMask);
    thisSubject = userOptions.subjectNames{s};
    model_val.name = 'value';
    model_val.RDM = RDMs{s}.val;
    model_val.color = [0 1 0];
    model_con.name = 'confidence';
    model_con.RDM = RDMs{s}.con;
    model_con.color = [0 1 0];
    model_fam.name = 'familiarity';
    model_fam.RDM = RDMs{s}.fam;
    model_fam.color = [0 1 0];
    model_pri.name = 'price';
    model_pri.RDM = RDMs{s}.pri;
    model_pri.color = [0 1 0];
    [rs_val{s},ps_val{s},ns_val{s},~] = searchlightMapping_fMRI(responsePatterns.(thisSubject), model_val, binaryMask, userOptions, searchlightOptions); %#ok<SAGROW>
    [rs_con{s},ps_con{s},ns_con{s},~] = searchlightMapping_fMRI(responsePatterns.(thisSubject), model_con, binaryMask, userOptions, searchlightOptions); %#ok<SAGROW>
    [rs_fam{s},ps_fam{s},ns_fam{s},~] = searchlightMapping_fMRI(responsePatterns.(thisSubject), model_fam, binaryMask, userOptions, searchlightOptions); %#ok<SAGROW>
    [rs_pri{s},ps_pri{s},ns_pri{s},~] = searchlightMapping_fMRI(responsePatterns.(thisSubject), model_pri, binaryMask, userOptions, searchlightOptions); %#ok<SAGROW>
end

% save stuff
for ss = 1:length(subs)
    rs = rs_val{ss}; %#ok<*NASGU>
    ps = ps_val{ss};
    ns = ns_val{ss};
    save([dir.out,fs,analysisName,fs,'val',fs,'sl_val_SF',num2str(subs(ss),'%03d')],'rs','ps','ns')
    rs = rs_con{ss};
    ps = ps_con{ss};
    ns = ns_con{ss};
    save([dir.out,fs,analysisName,fs,'con',fs,'sl_con_SF',num2str(subs(ss),'%03d')],'rs','ps','ns')
    rs = rs_fam{ss};
    ps = ps_fam{ss};
    ns = ns_fam{ss};
    save([dir.out,fs,analysisName,fs,'fam',fs,'sl_fam_SF',num2str(subs(ss),'%03d')],'rs','ps','ns')
    rs = rs_pri{ss}; 
    ps = ps_pri{ss};
    ns = ns_pri{ss};
    save([dir.out,fs,analysisName,fs,'pri',fs,'sl_pri_SF',num2str(subs(ss),'%03d')],'rs','ps','ns')
end

%% run searchlight for goal
% for s = 1:length(subs)
%     [rs_val{s},ps_val{s},ns_val{s},~] = searchlightMapping_fMRI(responsePatterns.(thisSubject), model_val, binaryMask, userOptions, searchlightOptions); %#ok<SAGROW>
% end