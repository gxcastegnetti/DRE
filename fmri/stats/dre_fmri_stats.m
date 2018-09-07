%% DRE_fMRI_stats
% ~~~
% GX Castegnetti --- start ~ 13.06.18 --- last ~ 17.08.18

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

%% pulse - value pmod of imagination; dValue pmod of choice
analysisName = 'uni_pulse_iV_cV';

% 1st level
timing.iDur = 0; % duration for imagination
timing.cDur = 0; % duration for choice
dre_L1_iV_cV(dir,analysisName,subs,timing,bData);

% contrasts
dre_con_i1_c1(dir,analysisName,subs);

% 2nd level
dre_L2(dir,analysisName,'imagination_onset',subs,1);
dre_L2(dir,analysisName,'imagination_value',subs,2);
dre_L2(dir,analysisName,'choice_onset',subs,3);
dre_L2(dir,analysisName,'choice_dValue',subs,4);

%% box - value pmod of imagination; dValue pmod of choice
analysisName = 'uni_box_iV_cV';

% 1st level
timing.iDur = 5;   % duration for imagination
timing.cDur = 3.5; % duration for choice
dre_L1_iV_cV(dir,analysisName,subs,timing,bData);

% contrasts
dre_con_i1_c1(dir,analysisName,subs);

% 2nd level
dre_L2(dir,analysisName,'imagination_onset',subs,1);
dre_L2(dir,analysisName,'imagination_value',subs,2);
dre_L2(dir,analysisName,'choice_onset',subs,3);
dre_L2(dir,analysisName,'choice_dValue',subs,4);


%% check movement
% if false
%     analysisName = 'uni_pulse_checkMov_3cond';
%
%     % extract behavioural data
%     bData = dre_extractData(dir,subs,taskOrd,1);
%
%     % extract behavioural data
%     bData = dre_extractData(dir,subs,taskOrd,1);
%
%     % 1st level
%     dre_level1_checkMov(dir,analysisName,subs,bData);
%
%     % contrasts
%     dre_contrasts_checkMov(dir,analysisName,subs);
%
%     % 2nd level
%     dre_level2(dir,analysisName,'I',subs,1);
%     dre_level2(dir,analysisName,'C',subs,2);
%     dre_level2(dir,analysisName,'C-I',subs,3);
%     dre_level2(dir,analysisName,'I',subs,1);
%     dre_level2(dir,analysisName,'C',subs,2);
%     dre_level2(dir,analysisName,'C-I',subs,3);
%
% end
