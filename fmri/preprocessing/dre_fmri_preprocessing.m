%% dre_fmri_preprocessing
% ~~~
% GX Castegnetti --- start ~ 09.05.18 --- last ~ 16.06.18

clear
close all
restoredefaultpath

%% folders
foodir  = pwd;
fs      = filesep;
idcs    = strfind(foodir,'/');
dir.dre = foodir(1:idcs(end-2));
dir.dsk = foodir(1:idcs(end-4));
dir.scn = [dir.dre,'data',fs,'fmri',fs,'scanner'];
dir.fmp = [dir.dre,'out',fs,'fmri',fs,'fmap'];
dir.spm = [dir.dsk,'tools',fs,'matlab',fs,'spm12'];
addpath([pwd,filesep,'routines'])

addpath(genpath(dir.spm))

%% subjects
subs = [4 5 8 9 13:17 19:21 23 25:26 29:32 34 35 37 39];

%% input run types (fun: functional; struct: structural; fmapm: fmap magnitude; fmapp: fmap phase; loc: localiser)
runType{1}  = {'loc','ignore','fun','fun','fmapm','fmapp','struct','fun'};
runType{2}  = {};
runType{3}  = {'loc','fun','fmapm','fmapp','fun','loc','struct','ignore','ignore','ignore','fun','fmapm','fmapp','fun'};
runType{4}  = {'loc','fun','fmapm','fmapp','fun','struct','fun','fun'};
runType{5}  = {'loc','ignore','fmapm','fmapp','ignore','fun','fun','fun','fun','struct'};
runType{6}  = {'loc','ignore','fmapm','fmapp','fun','fun','ignore','loc','ignore','fmapm','fmapp','fun','ignore','fun','struct'};
runType{7}  = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{8}  = {'loc','ignore','fmapm','fmapp','struct','fun','fun','fun','fun'};
runType{9}  = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{10} = {};
runType{11} = {};
runType{12} = {'loc','ignore','fmapm','fmapp','ignore','fun','loc','ignore','ignore','ignore','fmapm','fmapp','fun','loc','ignore','fmapm','fmapp','fun','fun','struct'};
runType{13} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{14} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{15} = {'loc','ignore','fmapm','fmapp','fun','ignore','fun','fun','fun','struct'};
runType{16} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{17} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{18} = {};
runType{19} = {'loc','ignore','fmapm','fmapp','ignore','fun','fun','fun','fun','struct'};
runType{20} = {'loc','ignore','fmapm','fmapp','ignore','fun','fun','fun','fun','struct'};
runType{21} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','loc','ignore','fmapm','fmapp','fun','struct'};
runType{22} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{23} = {'loc','ignore','fmapm','fmapp','ignore','fun','ignore','fun','fun','fun','struct'};
runType{24} = {};
runType{25} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{26} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{27} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{28} = {};
runType{29} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{30} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{31} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{32} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{33} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{34} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{35} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{36} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{37} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
runType{38} = {};
runType{39} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};
% runType{40} = {'loc','ignore','fmapm','fmapp','fun','fun','fun','fun','struct'};

%% Import images
DICOM_import(dir,subs,runType)

%% Field map preprocessing
% process new field maps
fieldmap_preproc(dir,subs,runType);

% load all field maps names
fmapFiles = cell(subs(end),1);
for s = 1:length(subs)
    load([dir.fmp,fs,'fmapNames_SF',num2str(subs(s),'%03d')],'fileNames')
    fmapFiles{subs(s)} = fileNames;
end

%% spatial preprocessing
if false
    spatial_preproc(dir,subs,fmapFiles);
end

%% physiological regressors
addpath([dir.spm,fs,'toolbox',fs,'physio']);
addpath([dir.spm,fs,'toolbox',fs,'physio',fs,'son']);
physioRegressors(dir,subs)

%% save movement data on a .pdf
% fmri_plotMovements(dir,subs)

