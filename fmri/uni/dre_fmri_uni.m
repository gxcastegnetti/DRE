%% DRE_fMRI_stats
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% Folders
foodir  = pwd;
fs      = filesep;
idcs    = strfind(foodir,'/');
dir.dre = foodir(1:idcs(end-2)-1);
dir.datScn = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.behDat = [dir.dre,fs,'data',fs,'behaviour'];
dir.phy = [dir.dre,fs,'data',fs,'physio'];
dir.out = [dir.dre,fs,'out',fs,'fmri',fs,'uni'];
newdir4 = foodir(1:idcs(end-4));
dir.spm = [newdir4,'tools/matlab/spm12'];

addpath(genpath(dir.spm))
addpath(genpath([foodir,fs,'routines'])), clear foodir idcs newdir4 dirSPM

%% Subjects
subs = [4 5 7 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39 40 41 43 47:49 50];
taskOrd = [ones(1,10),2*ones(1,11),1,2,ones(1,4),2*ones(1,3) 1];

%% extract behavioural data
bData = dre_extractData(dir,subs,taskOrd,1);

%% pulse - value, confid, famil, pmod of imagination; value chosen - unchosen pmod of choice
if true
    analysisName = 'uni_pulse_iVCF_cS';
    context = 'F';
    
    % 1st level
    timing.iOns = 0; % onset for imagination
    timing.cOns = 0; % onset for choice
    timing.iDur = 0; % duration for imagination
    timing.cDur = 0; % duration for choice
%     dre_L1_iVCF_cS(dir,analysisName,subs,timing,bData);
    
    % contrasts
    dre_con_i3_c1(dir,analysisName,subs,taskOrd,context);
    
    % 2nd level
    dre_L2(dir,analysisName,[context,'_imagination_onset'],subs,1);
    dre_L2(dir,analysisName,[context,'_imagination_value'],subs,2);
    dre_L2(dir,analysisName,[context,'_imagination_confid'],subs,3);
    dre_L2(dir,analysisName,[context,'_imagination_famil'],subs,4);
    dre_L2(dir,analysisName,[context,'_choice_onset'],subs,5);
    dre_L2(dir,analysisName,[context,'_choice_vChosen'],subs,6);
end

%% FOR RFXPLOT
if false
    analysisName = 'uni_pulse_rfx_iV_cS';
    
    % 1st level
    timing.iOns = 0; % onset for imagination
    timing.cOns = 0; % onset for choice
    timing.iDur = 0; % duration for imagination
    timing.cDur = 0; % duration for choice
    dre_L1_rfx_iV_cS(dir,analysisName,subs,timing,bData);
    
    % contrasts
    dre_con_i1_c1(dir,analysisName,subs);
    
    % 2nd level
    dre_L2(dir,analysisName,'imagination_onset',subs,1);
    dre_L2(dir,analysisName,'imagination_value',subs,2);
    dre_L2(dir,analysisName,'choice_onset',subs,3);
    dre_L2(dir,analysisName,'choice_valueChosen',subs,4);
end

%% median splits
if false
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % imagin.: V; choice: dV %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    analysisName = 'uni_pulse_ikV_ckV';
    
    % 1st level
    timing.iOns = 0; % onset for imagination
    timing.cOns = 0; % onset for choice
    timing.iDur = 0; % duration for imagination
    timing.cDur = 0; % duration for choice
    dre_L1_ikV_ckV(dir,analysisName,subs,timing,bData);
    
    % contrasts
    dre_con_ik1_ck1(dir,analysisName,subs);
    
    % 2nd level
    dre_L2(dir,analysisName,'ima. val. H-L',subs,1);
    dre_L2(dir,analysisName,'dV H-L',subs,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % imagin.: C; choice: S %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    analysisName = 'uni_pulse_ikC_ckS';
    
    % 1st level
    timing.iOns = 0; % onset for imagination
    timing.cOns = 0; % onset for choice
    timing.iDur = 0; % duration for imagination
    timing.cDur = 0; % duration for choice
    dre_L1_ikC_ckS(dir,analysisName,subs,timing,bData);
    
    % contrasts
    dre_con_ik1_ck1(dir,analysisName,subs);
    
    % 2nd level
    dre_L2(dir,analysisName,'ima. conf. H-L',subs,1);
    dre_L2(dir,analysisName,'chosen val H-L',subs,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % imagin.: F; choice: U %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    analysisName = 'uni_pulse_ikF_ckU';
    
    % 1st level
    timing.iOns = 0; % onset for imagination
    timing.cOns = 0; % onset for choice
    timing.iDur = 0; % duration for imagination
    timing.cDur = 0; % duration for choice
    dre_L1_ikF_ckU(dir,analysisName,subs,timing,bData);
    
    % contrasts
    dre_con_ik1_ck1(dir,analysisName,subs);
    
    % 2nd level
    dre_L2(dir,analysisName,'ima. fam. H-L',subs,1);
    dre_L2(dir,analysisName,'cho. - unc. val H-L',subs,2);
    
end