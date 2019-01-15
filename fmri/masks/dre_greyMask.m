%% dre_greyMask
% Creates individual grey matter masks from thresholding tissue
% probability from the segmentation.
% ~~~
% GX Castegnetti --- 2018

clear
close all
restoredefaultpath

%% folders
fs      = filesep;
dir.msk = pwd;
idcs    = strfind(dir.msk,'/');
dir.dre = dir.msk(1:idcs(end-2)-1);
dir.out = [dir.dre,fs,'out',fs,'fmri',fs,'masks',fs,'gm_subj'];
dir.data = [dir.dre,fs,'data',fs,'fmri',fs,'scanner'];
addpath(genpath('/Users/gcastegnetti/Desktop/tools/matlab/spm12'))

%% subjects
subs = [4 5 7 8 9 13:17 19 20 21 23 25:26 29:32 34 35 37 39 40 41 43 47:49 50];
subs = [42 44 45];

for s = 1:length(subs)   
    
    % directories of functional and structural images
    dirStruct = [dir.data,fs,'SF',num2str(subs(s),'%03d'),fs,'struct'];
    dirFun = [dir.data,fs,'SF',num2str(subs(s),'%03d'),fs,'fun',fs,'S1'];
    
    %% coregister with EPI
    % take coregistered gm map from segmentation output
    d = spm_select('List', dirStruct, '^c1.*\.nii$');
    c1_file = cellstr([dirStruct,fs,d]);
    job{1}.spm.spatial.coreg.write.source = c1_file;
    
    % take mean coregistered EPI as a reference
    d = spm_select('List', dirFun, '^meanua.*\.nii$');  
    epiRef = cellstr([dirFun,fs,d]);
    job{1}.spm.spatial.coreg.write.ref = epiRef;
    
    % defaults
    job{1}.spm.spatial.coreg.write.roptions.interp = 4;
    job{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    job{1}.spm.spatial.coreg.write.roptions.mask = 0;
    job{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    
    % run job
    spm_jobman('run',job);
    clear job
    
    %% threshold map to extract mask
    % take coregistered gm map from segmentation output
    d = spm_select('List', dirStruct, '^rc1.*\.nii$');
    rc1_file = cellstr([dirStruct,fs,d]);
    job{1}.spm.util.imcalc.input = rc1_file;
    
    % output mask
    job{1}.spm.util.imcalc.output = ['gm_SF',num2str(subs(s),'%03d')];
    job{1}.spm.util.imcalc.outdir = {dir.out};
    
    % probability threshold
    job{1}.spm.util.imcalc.expression = 'i1>0.3';
    
    % defaults
    job{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    job{1}.spm.util.imcalc.options.dmtx = 0;
    job{1}.spm.util.imcalc.options.mask = 0;
    job{1}.spm.util.imcalc.options.interp = 1;
    job{1}.spm.util.imcalc.options.dtype = 4;
    
    % run job
    disp(['Creating grey matter mask for sub#',num2str(subs(s),'%03d'),'...']);
    spm_jobman('run',job);
    clear job
    
end