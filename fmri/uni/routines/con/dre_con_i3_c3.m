function dre_con_full(dir,anName,subs,taskOrd,context)
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
    
    %% find which sessions to consider
    if strcmp(context,'F')
        if taskOrd(s) == 1
            whichSessions = [1 0 1 0];
        else
            whichSessions = [0 1 0 1];
        end
    elseif strcmp(context,'B')
        if taskOrd(s) == 1
            whichSessions = [0 1 0 1];
        else
            whichSessions = [1 0 1 0];
        end
    elseif strcmp(context,'all')
        whichSessions = [1 1 1 1];
    end
    wS = whichSessions; % just to make it shorter
    
    %% imagination
    
    % onset
    job{1}.spm.stats.con.consess{1}.tcon.name = 'imagination_onset';
    job{1}.spm.stats.con.consess{1}.tcon.weights = [wS(1) 0 0 0 0 0 0 0 movNull wS(2) 0 0 0 0 0 0 0 movNull wS(3) 0 0 0 0 0 0 0 movNull wS(4) 0 0 0 0 0 0 0 movNull];
    job{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    % value
    job{1}.spm.stats.con.consess{2}.tcon.name = 'value';
    job{1}.spm.stats.con.consess{2}.tcon.weights = [0 wS(1) 0 0 0 0 0 0 movNull 0 wS(2) 0 0 0 0 0 0 movNull 0 wS(3) 0 0 0 0 0 0 movNull 0 wS(4) 0 0 0 0 0 0 movNull];
    job{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    % confidence
    job{1}.spm.stats.con.consess{3}.tcon.name = 'confidence';
    job{1}.spm.stats.con.consess{3}.tcon.weights = [0 0 wS(1) 0 0 0 0 0 movNull 0 0 wS(2) 0 0 0 0 0 movNull 0 0 wS(3) 0 0 0 0 0 movNull 0 0 wS(4) 0 0 0 0 0 movNull];
    job{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
    % price
    job{1}.spm.stats.con.consess{4}.tcon.name = 'price';
    job{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 wS(1) 0 0 0 0 movNull 0 0 0 wS(2) 0 0 0 0 movNull 0 0 0 wS(3) 0 0 0 0 movNull 0 0 0 wS(4) 0 0 0 0 movNull];
    job{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
    
    %% choice
    
    % onset
    job{1}.spm.stats.con.consess{5}.tcon.name = 'choice_onset';
    job{1}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 0 wS(1) 0 0 0 movNull 0 0 0 0 wS(2) 0 0 0 movNull 0 0 0 0 wS(3) 0 0 0 movNull 0 0 0 0 wS(4) 0 0 0 movNull];
    job{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    
    % first parametric modulator
    job{1}.spm.stats.con.consess{6}.tcon.name = 'chosen-unchosen CONGRUENT';
    job{1}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 0 0 wS(1) 0 0 movNull 0 0 0 0 0 wS(2) 0 0 movNull 0 0 0 0 0 wS(3) 0 0 movNull 0 0 0 0 0 wS(4) 0 0 movNull];
    job{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    
    % second parametric modulator
    job{1}.spm.stats.con.consess{7}.tcon.name = 'chosen-unchosen INCONGRUENT';
    job{1}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0 0 0 wS(1) 0 movNull 0 0 0 0 0 0 wS(2) 0 movNull 0 0 0 0 0 0 wS(3) 0 movNull 0 0 0 0 0 0 wS(4) 0 movNull];
    job{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
    
    % second parametric modulator
    job{1}.spm.stats.con.consess{8}.tcon.name = 'chosen-unchosen PRICE';
    job{1}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 0 0 0 0 wS(1) movNull 0 0 0 0 0 0 0 wS(2) movNull 0 0 0 0 0 0 0 wS(3) movNull 0 0 0 0 0 0 0 wS(4) movNull];
    job{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    
    %% run job
    job{1}.spm.stats.con.delete = 1;
    spm_jobman('run',job)
    
end