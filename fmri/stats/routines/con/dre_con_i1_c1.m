function dre_con_i1_c1(dir,anName,subs)
%% function dre_con_i1_c1(dir,anName,subs)
% computes contrasts of parametric modulators in the case of ONE parametric
% modulator for imagination and one for choice
% ~~~
% GX Castegnetti --- start ~ 11.07.18 --- last ~ 18.08.18

movNull = zeros(1,6);

for s = 1:length(subs)
    
    % update user
    disp(['Writing contrasts for sub#', num2str(subs(s),'%03d'),'...']);
    
    spmFile = [dir.out,filesep,anName,filesep,'SF',num2str(subs(s),'%03d'),filesep,'SPM.mat'];
    job{1}.spm.stats.con.spmmat = {spmFile};
    
    %% imagination
    
    % onset 
    job{1}.spm.stats.con.consess{1}.tcon.name = 'imagination_onset';
    job{1}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0 movNull 1 0 0 0 movNull 1 0 0 0 movNull 1 0 0 0 movNull];
    job{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    % first parametric modulator
    job{1}.spm.stats.con.consess{2}.tcon.name = 'imagination_1';
    job{1}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 0 movNull 0 1 0 0 movNull 0 1 0 0 movNull 0 1 0 0 movNull];
    job{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    %% choice
    
    % onset
    job{1}.spm.stats.con.consess{3}.tcon.name = 'choice_onset';
    job{1}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 0 movNull 0 0 1 0 movNull 0 0 1 0 movNull 0 0 1 0 movNull];
    job{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
    % first parametric modulator
    job{1}.spm.stats.con.consess{4}.tcon.name = 'choice_1';
    job{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1 movNull 0 0 0 1 movNull 0 0 0 1 movNull 0 0 0 1 movNull];
    job{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
    job{1}.spm.stats.con.delete = 1;
    spm_jobman('run',job)
    
end