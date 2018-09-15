function spmMat_name = dre_L1_iVw_cSs_opp(dir,analysisName,subs,timing,bData)
%% function dre_L1_iV_cV(dirSub,sub,runType)
% ~~~
% First level analysis with conditions:
%   * imagination
%       - pmod: value weighed by confidence
%   * choice
%       - pmod: value of chosen item minus value of unchosen item
% ~~~
% GX Castegnetti --- start ~ 07.09.18 --- last ~ 07.09.18

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
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).name = ['imagination_',sessType];
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).onset = bData(subs(s)).imagination(r).(sessType).onset + timing.iOns;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).duration = timing.iDur;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).tmod = 0;
        
        % parametric modulations by value weighed by confidence
        confNorm = bData(subs(s)).imagination(r).(sessType).confidence/50; % extract confidence and map it onto [0,1]
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(1).name = 'price';
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(1).param = bData(subs(s)).imagination(r).(sessType).price;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(1).poly = 1;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).orth = 0;
        
        
        %% choice
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).name = ['choice_',sessType];
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).onset = bData(subs(s)).choice(r).(sessType).onset + timing.cOns;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).duration = timing.cDur;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).tmod = 0;
        
        % parametric modulation by value of the chosen item minus the
        % value of the unchosed item
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(1).name = 'value chosen - unchosen OPP';
        choVal = bData(subs(s)).choice(r).(sessTypeWRONG).valueChosen; % value of chosen item
        uncVal = bData(subs(s)).choice(r).(sessTypeWRONG).valueUnchosen; % value of unchosen item
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(1).param = choVal - uncVal;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(1).poly = 1;
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