function spmMat_name = dre_level1(dir,analysisName,subs,bData)
%% function dre_1stLevel(dirSub,sub,runType)
% ~~~
% INPUTS:
%   dir: directories
%   subj: subject number
%   analysisName: output folder
%   bdData: behavioural data struct created in dre_extractData.m
% ~~~
% GX Castegnetti --- start ~ 14.07.18 --- last ~ 18.08.18

fs = filesep;
n_sess = 4;

%% loop subjects
for s = 1:length(subs)
    
    % update user
    disp(['Preparing 1st level model for sub#', num2str(subs(s),'%03d'),'...']);
    
    %% parameters
    job1LM{1}.spm.stats.fmri_spec.timing.units = 'secs';
    job1LM{1}.spm.stats.fmri_spec.timing.RT = 3.36;
    job1LM{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    job1LM{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    %% folders
    dirSub = [dir.dre,fs,'data',fs,'fmri',fs,'scanner',fs,'SF',num2str(subs(s),'%03d')];
    dirOut = [dir.out,fs,analysisName,fs,'SF',num2str(subs(s),'%03d')];
    mkdir(dirOut)
    job1LM{1}.spm.stats.fmri_spec.dir = {dirOut};
    
    for r = 1:n_sess
        %% select EPI files
        dirFun = [dirSub,'/fun/S',num2str(r)];
        d = spm_select('List', dirFun, '^swuaf.*\.nii$');
        files = cellstr([repmat([dirFun fs],size(d,1),1) d]);
        job1LM{1}.spm.stats.fmri_spec.sess(r).scans = files; clear files
        
        %% extract session type
        sessType = bData(subs(s)).sessType{r};
        
        %% imagination
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).name = ['imagina_',sessType];
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).onset = bData(subs(s)).imagination(r).(sessType).onset;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).duration = 5;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).tmod = 0;
        
        % parametric modulations of value and confidence
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(1).name = 'value';
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(1).param = bData(subs(s)).imagination(r).(sessType).value;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(1).poly = 1;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(2).name = 'confidence';
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(2).param = bData(subs(s)).imagination(r).(sessType).confidence;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(2).poly = 1;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(3).name = 'familiarity';
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(3).param = bData(subs(s)).imagination(r).(sessType).familiarity;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(3).poly = 1;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).orth = 0;
        
        
        %% choice
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).name = ['choice_',sessType];
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).onset = bData(subs(s)).choice(r).(sessType).onset;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).duration = 3.5;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).tmod = 0;
        
        % parametric modulation of value and confidence difference
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(1).name = 'valueDiff';
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(1).param = bData(subs(s)).choice(r).(sessType).valueDiff;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(1).poly = 1;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(2).name = 'valueChosen';
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(2).param = bData(subs(s)).choice(r).(sessType).valueChosen;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(2).poly = 1;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).orth = 0;
        
        
        %% select movement regressors
        d = spm_select('List', dirFun, '^rp_af.*\.txt$');
        rp_file = cellstr([repmat([dirFun fs],size(d,1),1) d]);
        
        job1LM{1}.spm.stats.fmri_spec.sess(r).multi = {''};
        job1LM{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
        job1LM{1}.spm.stats.fmri_spec.sess(r).multi_reg = rp_file;
        job1LM{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
        
    end
    
    %% other stuff
    job1LM{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    job1LM{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    job1LM{1}.spm.stats.fmri_spec.volt = 1;
    job1LM{1}.spm.stats.fmri_spec.global = 'None';
    job1LM{1}.spm.stats.fmri_spec.mthresh = 0.8;
    job1LM{1}.spm.stats.fmri_spec.mask = {''};
    job1LM{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    %% run job and save file names
    d = spm_jobman('run',job1LM);
    spmMat_name{s} = d{1}.spmmat{1}; %#ok<AGROW>
    clear job1LM
    
    %% estimate model
    jobEst{1}.spm.stats.fmri_est.spmmat = {spmMat_name{s}};
    jobEst{1}.spm.stats.fmri_est.write_residuals = 0;
    jobEst{1}.spm.stats.fmri_est.method.Classical = 1;
    mkdir([dirSub,fs,'out',fs,'uni'])
    cd([dirSub,fs,'out',fs,'uni'])
    d = spm_jobman('run',jobEst);
    clear jobEst
    
end