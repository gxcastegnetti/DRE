%% dre_rsa_sl_perm
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% analysisName
analysisName = 'rsa_sl_ima';

%% directories
dir.rsaCod = pwd;
fs         = filesep;
idcs       = strfind(dir.rsaCod,'/');
dir.dre    = dir.rsaCod(1:idcs(end-2)-1); clear idcs
dir.datScn = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
dir.mskOut = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'atlas'];
dir.dimOut = [dir.dre,fs,'out',fs,'fmri',fs,'dim'];
dir.behDat = [dir.dre,fs,'data',fs,'behaviour'];
dir.analys = [dir.dre,fs,'out',fs,'fmri',fs,'rsa',fs,'sl',fs,analysisName];

% paths
dir.spm  = '/Users/gcastegnetti/Desktop/tools/matlab/spm12';
addpath([dir.rsaCod,fs,'routines'])
addpath(genpath([dir.rsaCod,fs,'rsatoolbox']))
addpath(genpath(dir.spm))

%% specify model(s) to test
models = {'cxt'};

%% run batch for each model
for m = 1:length(models)
    
    model = models{m};
    
    %%%%%%%%%%
    % define %
    %%%%%%%%%%
    
    dirData = [dir.analys,fs,model];
    
    % output directory
    dirOut = [dir.analys,fs,model,fs,'cluster'];
    job{1}.spm.tools.snpm.des.OneSampT.dir = {dirOut};
    
    % files
    d = spm_select('List', dirData, '^swrMap_c.*\.nii$');
    files  = cellstr([repmat([dirData fs],size(d,1),1) d]);
    job{1}.spm.tools.snpm.des.OneSampT.P = files;
    
    % cluster inference
    job{1}.spm.tools.snpm.des.OneSampT.bVolm = 1;
    %     job{1}.spm.tools.snpm.des.OneSampT.ST.ST_later = -1;
    
    % defaults
    job{1}.spm.tools.snpm.des.OneSampT.DesignName = 'MultiSub: One Sample T test on diffs/contrasts';
    job{1}.spm.tools.snpm.des.OneSampT.DesignFile = 'snpm_bch_ui_OneSampT';
    job{1}.spm.tools.snpm.des.OneSampT.cov = struct('c', {}, 'cname', {});
    job{1}.spm.tools.snpm.des.OneSampT.nPerm = 1000;
    job{1}.spm.tools.snpm.des.OneSampT.vFWHM = [0 0 0];
    job{1}.spm.tools.snpm.des.OneSampT.masking.tm.tm_none = 1;
    job{1}.spm.tools.snpm.des.OneSampT.masking.im = 1;
    job{1}.spm.tools.snpm.des.OneSampT.ST.ST_U = 0.001;
%     job{1}.spm.tools.snpm.des.OneSampT.masking.em = {'/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/masks/atlas/frontalLobe.nii'};
%     job{1}.spm.tools.snpm.des.OneSampT.masking.em = {'/Users/gcastegnetti/Desktop/stds/DRE/out/fmri/masks/_useNow/gmAvg.nii'};
    job{1}.spm.tools.snpm.des.OneSampT.masking.em = {''};
    job{1}.spm.tools.snpm.des.OneSampT.globalc.g_omit = 1;
    job{1}.spm.tools.snpm.des.OneSampT.globalm.gmsca.gmsca_no = 1;
    job{1}.spm.tools.snpm.des.OneSampT.globalm.glonorm = 1;
    
    % run job
    spm_jobman('run',job)
    clear job
    
    %%%%%%%%%%%
    % compute %
    %%%%%%%%%%%
    
    job{1}.spm.tools.snpm.cp.snpmcfg = {[dirOut,fs,'SnPMcfg.mat']};
    
    % run job
    spm_jobman('run',job)
    clear job
    
    %%%%%%%%%%%%%
    % inference %
    %%%%%%%%%%%%%
    job{1}.spm.tools.snpm.inference.SnPMmat = cellstr([dirOut,fs,'SnPM.mat']);
    job{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = nan;
    job{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.FWEthC = 0.05;
%         job{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.PthC = 0.15;
    job{1}.spm.tools.snpm.inference.Tsign = 1;
    job{1}.spm.tools.snpm.inference.WriteFiltImg.name = '_SnPM_filtered';
    job{1}.spm.tools.snpm.inference.Report = 'MIPtable';
    
    % run job
    spm_jobman('run',job)
    clear job
    
end