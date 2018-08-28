function dre_contrasts(dir,anName,subs)
%% function fieldmap_preprocess(dirSub,sub,runType)
% ~~~
% INPUTS:
%   subs: subjects
%   spmMat_name: name of the SPM.mat
% ~~~
% GX Castegnetti --- start ~ 11.07.18 --- last ~ 18.08.18

movNull = zeros(1,6);

for s = 1:length(subs)
    
    % update user
    disp(['Writing contrasts for sub#', num2str(subs(s),'%03d'),'...']);
    
    spmFile = [dir.out,filesep,anName,filesep,'SF',num2str(subs(s),'%03d'),filesep,'SPM.mat'];
    job{1}.spm.stats.con.spmmat = {spmFile};
    
    %% imagination
    job{1}.spm.stats.con.consess{1}.tcon.name = 'imaOns';
    job{1}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0 0 0 0 movNull 1 0 0 0 0 0 0 movNull 1 0 0 0 0 0 0 movNull 1 0 0 0 0 0 0 movNull];
    job{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    job{1}.spm.stats.con.consess{2}.tcon.name = 'imaVal';
    job{1}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 0 0 0 0 movNull 0 1 0 0 0 0 0 movNull 0 1 0 0 0 0 0 movNull 0 1 0 0 0 0 0 movNull];
    job{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    job{1}.spm.stats.con.consess{3}.tcon.name = 'imaCon';
    job{1}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 0 0 0 0 movNull 0 0 1 0 0 0 0 movNull 0 0 1 0 0 0 0 movNull 0 0 1 0 0 0 0 movNull];
    job{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
    job{1}.spm.stats.con.consess{4}.tcon.name = 'imaFam';
    job{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1 0 0 0 movNull 0 0 0 1 0 0 0 movNull 0 0 0 1 0 0 0 movNull 0 0 0 1 0 0 0 movNull];
    job{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
    %% choice
    job{1}.spm.stats.con.consess{5}.tcon.name = 'choOns';
    job{1}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 0 1 0 0 movNull 0 0 0 0 1 0 0 movNull 0 0 0 0 1 0 0 movNull 0 0 0 0 1 0 0 movNull];
    job{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    
    job{1}.spm.stats.con.consess{6}.tcon.name = 'choVal';
    job{1}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 0 0 1 0 movNull 0 0 0 0 0 1 0 movNull 0 0 0 0 0 1 0 movNull 0 0 0 0 0 1 0 movNull];
    job{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    
    job{1}.spm.stats.con.consess{7}.tcon.name = 'choCho';
    job{1}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0 0 0 1 movNull 0 0 0 0 0 0 1 movNull 0 0 0 0 0 0 1 movNull 0 0 0 0 0 0 1 movNull];
    job{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
    
    %% choice - imagination
    
    job{1}.spm.stats.con.consess{8}.tcon.name = 'cho-ima';
    job{1}.spm.stats.con.consess{8}.tcon.weights = [-1 0 0 0 1 0 0 movNull -1 0 0 0 1 0 0 movNull -1 0 0 0 1 0 0 movNull -1 0 0 0 1 0 0 movNull];
    job{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    
    job{1}.spm.stats.con.delete = 1;
    spm_jobman('run',job)
    
end