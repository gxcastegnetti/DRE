function dre_con_ik1_ck1(dir,anName,subs)
%% function dre_con_ik1_ck1(dir,anName,subs)
% computes contrasts of median split subjective evaluations
% ~~~
% GX Castegnetti --- 2018

movNull = zeros(1,6);

for s = 1:length(subs)
    
    % update user
    disp(['Writing contrasts for sub#', num2str(subs(s),'%03d'),'...']);
    
    spmFile = [dir.out,filesep,anName,filesep,'SF',num2str(subs(s),'%03d'),filesep,'SPM.mat'];
    job{1}.spm.stats.con.spmmat = {spmFile};
    
    %% imagination high - low whatever 
    job{1}.spm.stats.con.consess{1}.tcon.name = 'ima. H-L';
    job{1}.spm.stats.con.consess{1}.tcon.weights = [-1 1 0 0 movNull -1 1 0 0 movNull -1 1 0 0 movNull -1 1 0 0 movNull];
    job{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    %% choice high - low whatever 
    job{1}.spm.stats.con.consess{2}.tcon.name = 'cho. H-L';
    job{1}.spm.stats.con.consess{2}.tcon.weights = [0 0 -1 1 movNull 0 0 -1 1 movNull 0 0 -1 1 movNull 0 0 -1 1 movNull];
    job{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    %% run job
    job{1}.spm.stats.con.delete = 1;
    spm_jobman('run',job)
    
end