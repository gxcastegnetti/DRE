function spmMat_name = dre_level1_rsa_choice(dir,anName,subs,bData,maskName)
%% function DRE_First(dirSub,sub,runType)
% ~~~
% INPUTS:
%   dir: structure with directories
%   anName: name of the analysis
%   subs: subjects
%   bData: behavioural data
% ~~~
% GX Castegnetti --- 2018

fs = filesep;

%% loop subjects
for s = 1:length(subs)
    
    % update user
    disp(['Preparing 1st level model for sub#', num2str(subs(s),'%03d'),'...']);
    
    %% parameters
    job{1}.spm.stats.fmri_spec.timing.units = 'secs';
    job{1}.spm.stats.fmri_spec.timing.RT = 3.36;
    job{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    job{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    %% folders
    dirSub = [dir.dre,fs,'data',fs,'fmri',fs,'scanner',fs,'SF',num2str(subs(s),'%03d')];
    dirOut = [dir.rsaOut,fs,anName,fs,'SF',num2str(subs(s),'%03d')];
    mkdir(dirOut)
    job{1}.spm.stats.fmri_spec.dir = {dirOut};
    
    %% loop over sessions
    presentObj_F = 0;
    presentObj_B = 0;
    for r = 1:4
        %% select EPI files
        dirFun = [dirSub,'/fun/S',num2str(r)];
        d = spm_select('List', dirFun, '^uaf.*\.nii$');
        files = cellstr([repmat([dirFun fs],size(d,1),1) d]);
        job{1}.spm.stats.fmri_spec.sess(r).scans = files;
        
        %% extract session type and length
        sessType = bData(subs(s)).sessType{r};
        numChoices = 24;
        
        %% loop over trials in the session
        for choNum = 1:numChoices
            if strcmp(sessType,'F')
                presentObj_F = presentObj_F + 1;
                job{1}.spm.stats.fmri_spec.sess(r).cond(choNum).name = ['F-',num2str(presentObj_F)];
            elseif strcmp(sessType,'B')
                presentObj_B = presentObj_B + 1;
                job{1}.spm.stats.fmri_spec.sess(r).cond(choNum).name = ['B-',num2str(presentObj_B)];
            end
            job{1}.spm.stats.fmri_spec.sess(r).cond(choNum).onset = bData(subs(s)).choice(r).onset(choNum);
            job{1}.spm.stats.fmri_spec.sess(r).cond(choNum).duration = 0;
            job{1}.spm.stats.fmri_spec.sess(r).cond(choNum).tmod = 0;
            job{1}.spm.stats.fmri_spec.sess(r).cond(choNum).orth = 0;
        end
        
        %% select movement regressors
        d = spm_select('List', dirFun, '^rp_af.*\.txt$');
        rp_file = cellstr([repmat([dirFun fs],size(d,1),1) d]);
        
        job{1}.spm.stats.fmri_spec.sess(r).multi = {''};
        job{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
        job{1}.spm.stats.fmri_spec.sess(r).multi_reg = rp_file;
        job{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
        
    end
    
    %% other stuff
    job{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    job{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    job{1}.spm.stats.fmri_spec.volt = 1;
    job{1}.spm.stats.fmri_spec.global = 'None';
    job{1}.spm.stats.fmri_spec.mthresh = 0;
    
    % if there is a mask, apply it
    if strcmp(maskName,'none')
        mask = '';
    elseif strcmp(maskName,'gm')
        mask = [dir.mskOut,fs,'gm_subj',fs,'gm_SF',num2str(subs(s),'%03d'),'.nii'];
    else
        mask = [dir.mskOut,fs,maskName,'_subj',fs,'SF',num2str(subs(s),'%03d'),fs,'rw',maskName,'.nii'];
    end
    job{1}.spm.stats.fmri_spec.mask = {mask};
    job{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    %% run job and save file names
    d = spm_jobman('run',job);
    spmMat_name{s} = d{1}.spmmat{1}; %#ok<AGROW>
    clear job1LM
    
    %% estimate model
    jobEst{1}.spm.stats.fmri_est.spmmat = {spmMat_name{s}};
    jobEst{1}.spm.stats.fmri_est.write_residuals = 0;
    jobEst{1}.spm.stats.fmri_est.method.Classical = 1;
    d = spm_jobman('run',jobEst);
    clear jobEst
    
end