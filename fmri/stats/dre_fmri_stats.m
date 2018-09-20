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
dir.out = [dir.dre,fs,'out',fs,'fmri',fs,'stats'];
newdir4 = foodir(1:idcs(end-4));
dir.spm = [newdir4,'tools/matlab/spm12'];

addpath(genpath(dir.spm))
addpath(genpath([foodir,fs,'routines'])), clear foodir idcs newdir4 dirSPM

%% Subjects
subs = [4 5 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39];
taskOrd = [ones(1,9),2*ones(1,11),1,2,1];

%% extract behavioural data
bData = dre_extractData(dir,subs,taskOrd,1);

%% pulse - value, confid, famil, pmod of imagination; value selected pmod of choice
% analysisName = 'uni_pulse_iVCF_cS';
% 
% % 1st level
% timing.iOns = 0; % onset for imagination
% timing.cOns = 0; % onset for choice
% timing.iDur = 0; % duration for imagination
% timing.cDur = 0; % duration for choice
% dre_L1_iVCF_cS(dir,analysisName,subs,timing,bData);
% 
% % contrasts
% dre_con_i3_c1(dir,analysisName,subs);
% 
% % 2nd level
% dre_L2(dir,analysisName,'imagination_onset',subs,1);
% dre_L2(dir,analysisName,'imagination_value',subs,2);
% dre_L2(dir,analysisName,'imagination_confid',subs,3);
% dre_L2(dir,analysisName,'imagination_famil',subs,4);
% dre_L2(dir,analysisName,'choice_onset',subs,5);
% dre_L2(dir,analysisName,'choice_valueChosen',subs,6);

%% box - value, confid, famil, pmod of imagination; value selected pmod of choice
analysisName = 'uni_box_iVCF_cS_wPhy';

% 1st level
timing.iOns = 0; % onset for imagination
timing.cOns = 0; % onset for choice
timing.iDur = 5; % duration for imagination
timing.cDur = 3.5; % duration for choice
% dre_L1_iVCF_cS(dir,analysisName,subs,timing,bData);

% contrasts
dre_con_i3_c1(dir,analysisName,subs);

% 2nd level
dre_L2(dir,analysisName,'imagination_onset',subs,1);
dre_L2(dir,analysisName,'imagination_value',subs,2);
dre_L2(dir,analysisName,'imagination_confid',subs,3);
dre_L2(dir,analysisName,'imagination_famil',subs,4);
dre_L2(dir,analysisName,'choice_onset',subs,5);
dre_L2(dir,analysisName,'choice_valueChosen',subs,6);

%% pulse - value, confid, famil, pmod of imagination; value selected pmod of choice
% analysisName = 'uni_pulse_iV_cV';
% 
% % 1st level
% timing.iOns = 0; % onset for imagination
% timing.cOns = 0; % onset for choice
% timing.iDur = 0; % duration for imagination
% timing.cDur = 0; % duration for choice
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

%% box - value weighed by conf. pmod of imagination; value chosen - unchosen pmod of choice
% analysisName = 'uni_box_iVw_cSs';
% 
% % 1st level
% timing.iOns = 0; % onset for imagination
% timing.cOns = 0; % onset for choice
% timing.iDur = 5; % duration for imagination
% timing.cDur = 3.5; % duration for choice
% dre_L1_iVw_cSs(dir,analysisName,subs,timing,bData);
% 
% % contrasts
% dre_con_i1_c1(dir,analysisName,subs);
% 
% % 2nd level
% dre_L2(dir,analysisName,'imag. onset',subs,1);
% dre_L2(dir,analysisName,'imag. value weighed by conf.',subs,2);
% dre_L2(dir,analysisName,'choice onset',subs,3);
% dre_L2(dir,analysisName,'choice value ch. - unch.',subs,4);

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
