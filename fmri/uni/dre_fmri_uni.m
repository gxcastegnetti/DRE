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
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.beh = [dir.dre,fs,'data',fs,'behaviour'];
dir.phy = [dir.dre,fs,'data',fs,'physio'];
dir.out = [dir.dre,fs,'out',fs,'fmri',fs,'uni'];
newdir4 = foodir(1:idcs(end-4));
dir.spm = [newdir4,'tools/matlab/spm12'];

addpath(genpath(dir.spm))
addpath(genpath([foodir,fs,'routines'])), clear foodir idcs newdir4 dirSPM

%% Subjects
subs = [4 5 8 9 13:17 19 21 23 25:26 29:32 34 35 37 39 40 41 43 47:49];
taskOrd = [ones(1,9),2*ones(1,10),1,2,ones(1,4),2*ones(1,3)];

%% extract behavioural data
bData = dre_extractData(dir,subs,taskOrd,1);

%% pulse - value, confid, famil, pmod of imagination; value selected pmod of choice
if false
    analysisName = 'uni_pulse_iVCF_cS';
    
    % 1st level
    timing.iOns = 0; % onset for imagination
    timing.cOns = 0; % onset for choice
    timing.iDur = 0; % duration for imagination
    timing.cDur = 0; % duration for choice
    dre_L1_iVCF_cS(dir,analysisName,subs,timing,bData);
    
    % contrasts
    dre_con_i3_c1(dir,analysisName,subs,taskOrd);
    
    % 2nd level
    dre_L2(dir,analysisName,'imagination_onset',subs,1);
    dre_L2(dir,analysisName,'imagination_value',subs,2);
    dre_L2(dir,analysisName,'imagination_confid',subs,3);
    dre_L2(dir,analysisName,'imagination_famil',subs,4);
    dre_L2(dir,analysisName,'choice_onset',subs,5);
    dre_L2(dir,analysisName,'choice_valueChosen',subs,6);
end

%% box - value, confid, famil, pmod of imagination; value selected pmod of choice
if false
    analysisName = 'uni_box_iVCF_cS';
    
    % 1st level
    timing.iOns = 0; % onset for imagination
    timing.cOns = 0; % onset for choice
    timing.iDur = 5; % duration for imagination
    timing.cDur = 3.5; % duration for choice
    dre_L1_iVCF_cS(dir,analysisName,subs,timing,bData);
    
    % contrasts
    dre_con_i3_c1(dir,analysisName,subs,taskOrd);
    
    % 2nd level
    dre_L2(dir,analysisName,'imagination_onset',subs,1);
    dre_L2(dir,analysisName,'imagination_value',subs,2);
    dre_L2(dir,analysisName,'imagination_confid',subs,3);
    dre_L2(dir,analysisName,'imagination_famil',subs,4);
    dre_L2(dir,analysisName,'choice_onset',subs,5);
    dre_L2(dir,analysisName,'choice_valueChosen',subs,6);
end

%% pulse - value, confid, famil, pmod of imagination; value selected pmod of choice
if false
    analysisName = 'uni_pulse_iV_cV';
    
    % 1st level
    timing.iOns = 0; % onset for imagination
    timing.cOns = 0; % onset for choice
    timing.iDur = 0; % duration for imagination
    timing.cDur = 0; % duration for choice
    dre_L1_iV_cV(dir,analysisName,subs,timing,bData);
    
    % contrasts
    dre_con_i1_c1(dir,analysisName,subs);
    
    % 2nd level
    dre_L2(dir,analysisName,'imagination_onset',subs,1);
    dre_L2(dir,analysisName,'imagination_value',subs,2);
    dre_L2(dir,analysisName,'choice_onset',subs,3);
    dre_L2(dir,analysisName,'choice_dV',subs,4);
end

%% pulse - value pmod of imagination; dV pmod of choice
if true
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

%% box - value, confid, famil, pmod of imagination; value selected pmod of choice
% analysisName = 'uni_box_iV_cV';
%
% % 1st level
% timing.iOns = 0; % onset for imagination
% timing.cOns = 0; % onset for choice
% timing.iDur = 5; % duration for imagination
% timing.cDur = 3.5; % duration for choice
% dre_L1_iV_cV(dir,analysisName,subs,timing,bData);
%
% % contrasts
% dre_con_i1_c1(dir,analysisName,subs);
%
% % 2nd level
% dre_L2(dir,analysisName,'imagination_onset',subs,1);
% dre_L2(dir,analysisName,'imagination_value',subs,2);
% dre_L2(dir,analysisName,'choice_onset',subs,3);
% dre_L2(dir,analysisName,'choice_valueDiff',subs,4);

%% pulse - value weighed by conf. pmod of imagination; value chosen - unchosen pmod of choice
if false
    analysisName = 'uni_pulse_iVw_cSs';
    
    % 1st level
    timing.iOns = 0; % onset for imagination
    timing.cOns = 0; % onset for choice
    timing.iDur = 0; % duration for imagination
    timing.cDur = 0; % duration for choice
    dre_L1_iVw_cU(dir,analysisName,subs,timing,bData);
    
    % contrasts
    dre_con_i1_c1(dir,analysisName,subs);
    
    % 2nd level
    dre_L2(dir,analysisName,'imag. onset',subs,1);
    dre_L2(dir,analysisName,'imag. value weighed by conf.',subs,2);
    dre_L2(dir,analysisName,'choice onset',subs,3);
    dre_L2(dir,analysisName,'choice value ch. - unch.',subs,4);
end

%% box - price pmod of imagination; value WRONG GOAL chosen - unchosen pmod of choice
% analysisName = 'uni_box_iVw_cSs_opp';
%
% % 1st level
% timing.iOns = 0; % onset for imagination
% timing.cOns = 0; % onset for choice
% timing.iDur = 5; % duration for imagination
% timing.cDur = 3.5; % duration for choice
% dre_L1_iVw_cSs_opp(dir,analysisName,subs,timing,bData);
%
% % contrasts
% dre_con_i1_c1(dir,analysisName,subs);
%
% % 2nd level
% dre_L2(dir,analysisName,'imag. onset',subs,1);
% dre_L2(dir,analysisName,'imag. price',subs,2);
% dre_L2(dir,analysisName,'choice onset',subs,3);
% dre_L2(dir,analysisName,'choice value ch. - unch. OPP',subs,4);


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