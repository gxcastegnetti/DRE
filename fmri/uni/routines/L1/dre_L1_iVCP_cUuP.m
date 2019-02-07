function spmMat_name = dre_L1_iVCP_cUuP(dir,analysisName,subs,timing,bData)
%% function dre_L1_iVCF_cV(dirSub,sub,runType)
% ~~~
% First level analysis with conditions:
%   * imagination
%       - pmod: value, familiarity
%   * choice
%       - pmod: value chosen - value unchosen
% ~~~
% GX Castegnetti --- 2018

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
    %     dirPhy = [dir.phy,fs,num2str(subs(s),'%03d')];
    %     dirPhyOut = [dir.dre,fs,'out',fs,'fmri',fs,'preprocessing',fs,'body_regr',fs,'SF',num2str(subs(s),'%03d')];
    %     if ~exist(dirPhyOut,'dir'), mkdir(dirPhyOut), end
    if ~exist(dirOut,'dir'), mkdir(dirOut), end
    
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
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).name = ['ima. ',sessType];
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).onset = bData(subs(s)).imagination(r).onset + timing.iOns;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).duration = timing.iDur;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).tmod = 0;
        
        % parametric modulations of value
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(1).name = 'imagin. value';
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(1).param = bData(subs(s)).imagination(r).val;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(1).poly = 1;
        
        % parametric modulations of confidence
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(2).name = 'imagin. conf.';
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(2).param = bData(subs(s)).imagination(r).con;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(2).poly = 1;
        
        % parametric modulations of familiarity
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(3).name = 'imagin. price';
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(3).param = bData(subs(s)).imagination(r).pri;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod(3).poly = 1;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(1).orth = 0;
        
        
        %% choice
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).name = ['choice_',sessType];
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).onset = bData(subs(s)).choice(r).onset + timing.cOns;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).duration = timing.cDur;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).tmod = 0;
        
        % parametric modulation by value of the chosen minus unchosed item
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(1).name = 'chosen CONGRUENT';
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(1).param = bData(subs(s)).choice(r).valCho;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(1).poly = 1;
        
        % parametric modulation by value of the chosen minus unchosed item
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(2).name = 'chosen INCONGRUENT';
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(2).param = bData(subs(s)).choice(r).valChoWrong;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(2).poly = 1;
        
        % parametric modulation by value of the chosen minus unchosed item
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(3).name = 'chosen PRICE';
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(3).param = bData(subs(s)).choice(r).priCho;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod(3).poly = 1;
        job1LM{1}.spm.stats.fmri_spec.sess(r).cond(2).orth = 0;      
        
        %% movement (and physio) regressors
        
        d_mov = spm_select('List', dirFun, '^rp_af.*\.txt$');
        rp_file = cellstr([repmat([dirFun fs],size(d_mov,1),1) d_mov]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % use this if only movement regressors required %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        job1LM{1}.spm.stats.fmri_spec.sess(r).multi = {''};
        job1LM{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
        job1LM{1}.spm.stats.fmri_spec.sess(r).multi_reg = rp_file;
        job1LM{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % use this if also physiological regressors required %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %         % load rp file to concatenate with physio regressors
        %         rp = load(rp_file{1});
        %
        %         % see which file is correct
        %         fileF = [dirPhy,fs,'SF',num2str(subs(s),'%03d'),'_S',num2str(r),'F_R_S1.mat'];
        %         fileB = [dirPhy,fs,'SF',num2str(subs(s),'%03d'),'_S',num2str(r),'B_R_S1.mat'];
        %         try
        %             try
        %                 mov = load(fileF);
        %             catch
        %                 mov = load(fileB);
        %             end
        %             mov = mov.R;
        %         catch
        %             mov = [];
        %         end
        %
        %
        %         % concatenate movement and physiological regressors
        %         R = [rp, mov]; %#ok<NASGU>
        %
        %         % save
        %         allRegr_file = [dirPhyOut,fs,'body_r_SF',num2str(subs(s),'%03d'),'_S',num2str(r),'.mat'];
        %         save(allRegr_file,'R');
        %
        %         job1LM{1}.spm.stats.fmri_spec.sess(r).multi = {''};
        %         job1LM{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
        %         job1LM{1}.spm.stats.fmri_spec.sess(r).multi_reg = {allRegr_file};
        %         job1LM{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
        
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