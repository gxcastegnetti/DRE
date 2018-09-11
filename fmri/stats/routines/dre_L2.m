function d = dre_L2(dir,analysisName,modName,subs,contrNum)
%% create 2nd level model
% ~~~
% GX Castegnetti --- 2018

fs = filesep;

% output folder
dirOut = [dir.out,fs,analysisName,fs,'2nd_level',fs,modName];
job{1}.spm.stats.factorial_design.dir = {dirOut};

% Select contrast images
files = cell(length(subs),1);
for s = 1:length(subs)
    dirSub = [dir.out,fs,analysisName,fs,'SF',num2str(subs(s),'%03d')];
    files{s} = [dirSub,fs,'con_00',num2str(contrNum,'%02d'),'.nii'];
end

% defaults
job{1}.spm.stats.factorial_design.des.t1.scans = files;
job{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
job{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
job{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
job{1}.spm.stats.factorial_design.masking.im = 1;
job{1}.spm.stats.factorial_design.masking.em = {''};
job{1}.spm.stats.factorial_design.globalc.g_omit = 1;
job{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
job{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% run job
d = spm_jobman('run',job);
clear job

%% estimate model
job{1}.spm.stats.fmri_est.spmmat = d{1}.spmmat;
job{1}.spm.stats.fmri_est.write_residuals = 0;
job{1}.spm.stats.fmri_est.method.Classical = 1;

% run job
spm_jobman('run',job);