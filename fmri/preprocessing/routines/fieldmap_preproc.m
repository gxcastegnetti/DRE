function fileNames = fieldmap_preproc(dir,subs,runType)
%% function fieldmap_preprocess(dirSPM,dirMain,subs,runType)
% ~~~
% INPUTS:
%   dirSPM: SPM directory
%   dirMain: root directory
%   subs: subjects
%   runType: cell array of run types ('fun','fmap','struct','ignore')
% ~~~
% GX Castegnetti --- start ~ 04.06.18 --- last ~ 13.06.18

fs = filesep;

for s = 1:length(subs)
    
    % subject directory
    dirSub = [dir.dre,'data',fs,'fmri',fs,'scanner',fs,'SF',num2str(subs(s),'%03d')];
    
    if ~exist([dir.fmp,filesep,'fmapNames_SF',num2str(subs(s),'%03d'),'.mat'],'file')
        
        % find how many field maps and functional sessions were acquired
        nFmap = sum(strcmp(runType{subs(s)},'fmapm'));
        nFun = sum(strcmp(runType{subs(s)},'fun'));
        
        % if more than one field map, find which belongs to which
        if nFmap == 2
            runLoc = find(strcmp(runType{subs(s)},'loc'));
            runFun = find(strcmp(runType{subs(s)},'fun'));
            fun_fmap{1} = 1:sum(runFun < runLoc(2));
            fun_fmap{2} = (fun_fmap{1}(end)+1):nFun;
        else
            fun_fmap{1} = 1:nFun;
        end
        
        if subs(s) == 12
            fun_fmap{1} = 1;
            fun_fmap{2} = 2;
            fun_fmap{3} = [3 4];
        end
        
        
        %% loop over field map and apply correction to the associated functional scans
        l = 1;
        for i = 1:nFmap
            
            % select magnitude and phase difference files
            fmPhaseDir = [dirSub,fs,'fmap',fs,'S',num2str(i),fs,'phase'];
            fmMagniDir = [dirSub,fs,'fmap',fs,'S',num2str(i),fs,'magni'];
            filePhase = spm_select('List', fmPhaseDir,'^sM.*\.nii$');
            fileMagni = spm_select('List', fmMagniDir,'^sM.*\.nii$');
            fileMagni = fileMagni(1,:);
            job{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = {[fmPhaseDir,fs,filePhase]};
            job{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = {[fmMagniDir,fs,fileMagni]};
            
            % select defaults for FIL's Siemens Prisma 3.0T (Quattro)
            job{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsfile = {[dir.spm,fs,'toolbox',fs,'FieldMap',fs,'pm_defaults_Prisma.m']};
            
            % loop over runs
            for j = 1:length(fun_fmap{i})
                
                % select folder for current session
                dirSess = [dirSub,fs,'fun',fs,'S',num2str(fun_fmap{i}(j))];
                
                % select first EPI from each session
                fileEPIs = spm_select('List', dirSess,'^f.*\.nii$');
                fileEPI = fileEPIs(1,:); clear fileEPIs
                job{1}.spm.tools.fieldmap.calculatevdm.subj.session(j).epi = {[dirSess,fs,fileEPI]};
            end
            
            job{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
            job{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
            job{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
            job{1}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
            job{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
            
            % run job
            disp(['Preprocessing Field map for sub#',num2str(subs(s),'%03d'),'...']);
            d = spm_jobman('run',job);
            clear job
            
            for k = 1:length(fun_fmap{i})
                fileNames{l} = d{1}.vdmfile{k}{1};  %#ok<AGROW>
                l = l+1;
            end
            
            % save fieldmap
            save([dir.fmp,fs,'fmapNames_SF',num2str(subs(s),'%03d')],'fileNames')
        end
    end
end