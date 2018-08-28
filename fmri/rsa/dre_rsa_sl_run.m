%% dre_rsa_sl_run
% ~~~
% GX Castegnetti --- start ~ 17.07.18 --- last ~ 24.08.18

clear
close all
restoredefaultpath

%% folders
maskName = 'brain';
fs      = filesep;
dir.rsa = pwd;
idcs    = strfind(dir.rsa,'/');
dir.dre = dir.rsa(1:idcs(end-2)-1);
dir.sta = [dir.dre,fs,'codes',fs,'fmri',fs,'stats'];
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

%% 1st level for RSA
if false
    
    % extract behavioural data
    bData = dre_extractData(dir,subs,taskOrd,0);
        
    % reference EPI
    epiRef = [dir.data,fs,'SF039',fs,'fun',fs,'S4',fs,'wuafMQ04784-0008-00225-000225-01.nii,1'];
    maskUnwarp = [dir.dre,fs,'codes',fs,'fmri',fs,'stats',fs,'masks',fs,maskName,'.nii,1'];
    job{1}.spm.spatial.coreg.write.ref = {epiRef};
    job{1}.spm.spatial.coreg.write.source = {maskUnwarp};
    job{1}.spm.spatial.coreg.write.roptions.interp = 4;
    job{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    job{1}.spm.spatial.coreg.write.roptions.mask = 0;
    job{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    d = spm_jobman('run',job);
    dre_level1_rsa(dir,['r',maskName],subs,bData);
end


%% load fMRI data
userOptions = DRE_RSA_userOptions(dir,subs);
userOptions.analysisName = 'RSA_prova_rand';
userOptions.rootPath = dir.out;
userOptions.forcePromptReply = 'r';
userOptions.overwriteflag = 'r';

% load data with different masks
dir.beta = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'for_RSA_pulse',fs,'r',maskName];
userOptions.betaPath = [dir.beta,filesep,'[[subjectName]]',filesep,'[[betaIdentifier]]'];
[~, fullBrainVols_full] = rsa.fmri.fMRIDataPreparation('SPM', userOptions);


%% extract models of value, confidence, familiarity, price
RDMs = dre_extractRDMs(dir,subs,taskOrd);

userOptions.voxelSize = [3 3 3];
userOptions.searchlightRadius = 9;
searchlightOptions.monitor = false;
searchlightOptions.fisher = true;
searchlightOptions.nSessions = 1;
searchlightOptions.nConditions = 240;
binaryMask = niftiread([dir.sta,fs,'masks',fs,'rbrain.nii']);
binaryMask = logical(binaryMask);

for s = 1:length(subs)
    disp(['Computing correlation for sub#',num2str(s),' of ',num2str(length(subs))])
    thisSubject = userOptions.subjectNam es{s};
    model.name = 'value';
    model.RDM = RDMs{s}.val;
    model.color = [0 1 0];
    [rs{s},ps{s},ns{s},~] = rsa.fmri.searchlightMapping_fMRI(fullBrainVols_full.(thisSubject), model, binaryMask, userOptions, searchlightOptions); %#ok<ASGLU>
    
end

% save stuff
for ss = 1:length(subs)
    foo_rs = rs{ss};
    foo_ps = ps{ss};
    foo_ns = ns{ss};
    save([dir.out,fs,'searchLight',fs,'rs_rand_SF',num2str(subs(ss),'%03d')],'foo_rs')
    save([dir.out,fs,'searchLight',fs,'ps_rand_SF',num2str(subs(ss),'%03d')],'foo_ps')
    save([dir.out,fs,'searchLight',fs,'ns_rand_SF',num2str(subs(ss),'%03d')],'foo_ns')
end
clear rs ps ns
