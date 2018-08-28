function fmri_plotMovements(dir,subs)
%% fmri_plotMovements(dir,subs)
% This function takes movement file rp* from realignment from every
% subject, concatenates sessions, and plots them in a pdf
% ~~~
% GX Castegnetti 2018

fs = filesep;

addpath('/Users/gcastegnetti/Desktop/tools/matlab/export_fig')

for s = 1:length(subs)
    for r = 1:4
        dirSess = [dir.scn,fs,'SF',num2str(subs(s),'%03d'),fs,'fun',fs,'S',num2str(r)];
        d = spm_select('List', dirSess, '^rp.*\.txt$');
        fileID = fopen([dirSess,fs,d],'r');
        formatSpec = '%f %f %f %f %f %f';
        sizeA = [6 inf];
        movSess = fscanf(fileID,formatSpec,sizeA);
        transSess{r} = movSess(1:3,:)';
        rotatSess{r} = movSess(4:6,:)';
    end
    movSub = [transSess{1}; transSess{2}; transSess{3}; transSess{4}];
    trans = figure('color',[1 1 1]); plot(movSub);
    xlabel('Volume (4 sessions concat.)'), ylabel('\Deltar'),title(['sub#',num2str(subs(s),'%03d')])
    
    rotSub = [rotatSess{1}; rotatSess{2}; rotatSess{3}; rotatSess{4}];
    rotat = figure('color',[1 1 1]); plot(rotSub);
    xlabel('Volume (4 sessions concat.)'), ylabel('\Delta\theta'),title(['sub#',num2str(subs(s),'%03d')])
    
    export_fig('Translations','-pdf','-nocrop','-append',trans)
    export_fig('Rotations','-pdf','-nocrop','-append',rotat)
end