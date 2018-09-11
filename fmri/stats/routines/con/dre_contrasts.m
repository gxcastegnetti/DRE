function dre_contrasts(dir,anName,subs)
%% function fieldmap_preprocess(dirSub,sub,runType)
% ~~~
% INPUTS:
%   dir: struct with directories
%   anName: analysis name
%   subs: vector of subject numbers
% ~~~
% GX Castegnetti --- start ~ 11.07.18 --- last ~ 18.08.18

movNull = zeros(1,6);

for s = 1:length(subs)
    
    % update user
    disp(['Writing contrasts for sub#', num2str(subs(s),'%03d'),'...']);
    
    spmFile = [dir.out,filesep,anName,filesep,'SF',num2str(subs(s),'%03d'),filesep,'SPM.mat'];
    job{1}.spm.stats.con.spmmat = {spmFile};
    
    %% imagination
    job{1}.spm.stats.con.consess{1}.tcon.name = 'imagination_onset';
    job{1}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0 0 0 movNull 1 0 0 0 0 0 movNull 1 0 0 0 0 0 movNull 1 0 0 0 0 0 movNull];
    job{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    job{1}.spm.stats.con.consess{2}.tcon.name = 'imagination_value';
                                                   [0 1 0 0 0 0 movNull 0 1 0 0 0 0 movNull 0 1 0 0 0 0 movNull 0 1 0 0 0 0 movNull];
    job{1}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 0 0 0 movNull 0 1 0 0 0 0 movNull 0 1 0 0 0 0 movNull 0 1 0 0 0 0 movNull];
    job{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    job{1}.spm.stats.con.consess{3}.tcon.name = 'imagination_confid';
    job{1}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 0 0 0 movNull 0 0 1 0 0 0 movNull 0 0 1 0 0 0 movNull 0 0 1 0 0 0 movNull];
    job{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
    job{1}.spm.stats.con.consess{4}.tcon.name = 'imagination_famil';
    job{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1 0 0 movNull 0 0 0 1 0 0 movNull 0 0 0 1 0 0 movNull 0 0 0 1 0 0 movNull];
    job{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
    %% choice
    job{1}.spm.stats.con.consess{5}.tcon.name = 'choice_onset';
    job{1}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 0 1 0 movNull 0 0 0 0 1 0 movNull 0 0 0 0 1 0 movNull 0 0 0 0 1 0 movNull];
    job{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    
    job{1}.spm.stats.con.consess{6}.tcon.name = 'choice_dValue';
    job{1}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 0 0 1 movNull 0 0 0 0 0 1 movNull 0 0 0 0 0 1 movNull 0 0 0 0 0 1 movNull];
    job{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    
    job{1}.spm.stats.con.delete = 1;
    spm_jobman('run',job)
    
end