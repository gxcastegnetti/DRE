%% dre_rsa_sl_run
% ~~~
% GX Castegnetti --- start ~ 17.07.18 --- last ~ 24.08.18

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_sl_grey_1';

%% folders
fs      = filesep;
dir.rsa = pwd;
idcs    = strfind(dir.rsa,'/');
dir.dre = dir.rsa(1:idcs(end-2)-1);
dir.sta = [dir.dre,fs,'codes',fs,'fmri',fs,'stats'];
dir.msk = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'grey'];
dir.beh = [dir.dre,fs,'data',fs,'behaviour'];
dir.out = [dir.dre,fs,'out',fs,'fmri',fs,'rsa'];
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
addpath([dir.rsa,fs,'routines'])
addpath([dir.sta,fs,'routines'])
addpath(genpath([dir.rsa,fs,'rsatoolbox']))
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))

%% subjects
subs = [4:5 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39];
taskOrd = [ones(1,9),2*ones(1,11),1,2,1];

%% user options
userOptions = dre_rsa_userOptions(dir,subs);
userOptions.analysisName = analysisName;
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';
userOptions.overwriteflag = 'r';

%% 1st level for RSA
nameBeta = '_level1_gm';
if false
    bData = dre_extractData(dir,subs,taskOrd,0);
    dre_level1_rsa(dir,nameBeta,subs,bData);
end

%% load betas
dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,nameBeta];
userOptions.betaPath = [dir.beta,filesep,'[[subjectName]]',filesep,'[[betaIdentifier]]'];
if ~exist([dir.out,fs,'searchLight',fs,'fullBrainVols.mat'],'file')
    [~, fullBrainVols] = rsa.fmri.fMRIDataPreparation('SPM', userOptions);
    save([dir.out,fs,'searchLight',fs,'fullBrainVols.mat'],'fullBrainVols','-v7.3')
else
    load([dir.out,fs,'searchLight',fs,'fullBrainVols.mat'],'fullBrainVols')
end

%% extract models of value, confidence, familiarity, price
RDMs = dre_extractRDMs(dir,subs,taskOrd);
userOptions.voxelSize = [3 3 3];
userOptions.searchlightRadius = 9;
searchlightOptions.monitor = false;
searchlightOptions.fisher = true;
searchlightOptions.nSessions = 1;
searchlightOptions.nConditions = 240;

%% run searchlight
for s = 1:length(subs)
    disp(['Computing correlation for sub#',num2str(s),' of ',num2str(length(subs))])
    binaryMask = niftiread([dir.msk,fs,'grey_SF',num2str(subs(s),'%03d'),'.nii']);
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
    [rs_val{s},ps_val{s},ns_val{s},~] = rsa.fmri.searchlightMapping_fMRI(fullBrainVols.(thisSubject), model_val, binaryMask, userOptions, searchlightOptions); %#ok<SAGROW>
    [rs_con{s},ps_con{s},ns_con{s},~] = rsa.fmri.searchlightMapping_fMRI(fullBrainVols.(thisSubject), model_con, binaryMask, userOptions, searchlightOptions); %#ok<SAGROW>
    [rs_fam{s},ps_fam{s},ns_fam{s},~] = rsa.fmri.searchlightMapping_fMRI(fullBrainVols.(thisSubject), model_fam, binaryMask, userOptions, searchlightOptions); %#ok<SAGROW>
    [rs_pri{s},ps_pri{s},ns_pri{s},~] = rsa.fmri.searchlightMapping_fMRI(fullBrainVols.(thisSubject), model_pri, binaryMask, userOptions, searchlightOptions); %#ok<SAGROW>
end

% save stuff
for ss = 1:length(subs)
    rs = rs_val{ss}; %#ok<*NASGU>
    ps = ps_val{ss};
    ns = ns_val{ss};
    save([dir.out,fs,'searchLight',fs,analysisName,fs,'sl_val_SF',num2str(subs(ss),'%03d')],'rs','ps','ns')
    rs = rs_con{ss}; %#ok<*NASGU>
    ps = ps_con{ss};
    ns = ns_con{ss};
    save([dir.out,fs,'searchLight',fs,analysisName,fs,'sl_con_SF',num2str(subs(ss),'%03d')],'rs','ps','ns')
    rs = rs_fam{ss}; %#ok<*NASGU>
    ps = ps_fam{ss};
    ns = ns_fam{ss};
    save([dir.out,fs,'searchLight',fs,analysisName,fs,'sl_fam_SF',num2str(subs(ss),'%03d')],'rs','ps','ns')
    rs = rs_pri{ss}; %#ok<*NASGU>
    ps = ps_pri{ss};
    ns = ns_pri{ss};
    save([dir.out,fs,'searchLight',fs,analysisName,fs,'sl_pri_SF',num2str(subs(ss),'%03d')],'rs','ps','ns')
end
