function RDMs = dre_extractRDMs(dir,subs,taskOrd)
% ~~~
% INPUTS:
%   dirSub: directory of subject's data
%   subj: subject number
%   taskOrd: FBFB or BFBF
% OUTPUTS:
%   RDMs: representation dissimilarity matrices of behavioural judgements
% ~~~
% GX Castegnetti --- start ~ 24.08.18 --- last ~ 24.08.18

fs = filesep;
n_sess = 4;

dir.data = [dir.dre,fs,'data'];

%% read objects
objs        = readtable([dir.beh,fs,'Objects.csv']);
objsName    = table2cell(objs(:,2)); clear foo objs

%% loop over subjects and sessions
for s = 1:length(subs)
    for r = 1:n_sess
        
        %% extract subject-specific task order
        if taskOrd(s) == 1
            day2Order = {'F','B','F','B'};
        elseif taskOrd(s) == 2
            day2Order = {'B','F','B','F'};
        end
        
        %% subject directories
        dirPsyO = [dir.data,fs,'fmri',fs,'psychOut',fs,'SF',num2str(subs(s),'%03d')];
        dirBeha = [dir.data,fs,'behaviour',fs,'SF',num2str(subs(s),'%03d')];
        
        %% trial onsets
        
        % load behavioural matrix from day 2
        if ~exist([dirPsyO,fs,'DRE_mri_S',num2str(subs(s),'%03d'),'_B',num2str(r),day2Order{r},'.csv']), keyboard, end
        Mday2 = csvread([dirPsyO,fs,'DRE_mri_S',num2str(subs(s),'%03d'),'_B',num2str(r),day2Order{r},'.csv']);
        
        %% object IDs
        % extract IDs from day 2
        id(:,r) = Mday2(Mday2(:,3)==0,4);
        id_vec = id(:);
        
    end
    
    %% separate F and B
    if taskOrd(s) == 1
        id_F = [id(:,1); id(:,3)];
        id_B = [id(:,2); id(:,4)];
    else
        id_F = [id(:,2); id(:,4)];
        id_B = [id(:,1); id(:,3)];
    end
    
    % sort them
    [fooF,sort_F] = sort(id_F);
    [fooB,sort_B] = sort(id_B);
    
    %% value and confidence
    Mday1_F = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_B1_DRE.csv']);
    Mday1_B = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_B2_DRE.csv']);
    
    % object names
    name_F = objsName(id_F);
    name_B = objsName(id_B);
    
    %% find corresponding indices on day 1
    idx_day1_F = NaN(length(id_F),1);
    idx_day1_B = NaN(length(id_B),1);
    for i = 1:length(id_F)
        idx_day1_F(i) = find(id_F(i) == Mday1_F(:,2));
        idx_day1_B(i) = find(id_B(i) == Mday1_B(:,2));
    end
    
    % value
    val_F = Mday1_F(idx_day1_F,3);
    val_B = Mday1_B(idx_day1_B,3);
    val_all = [val_F(sort_F); val_B(sort_B)];
    
    % confidence
    con_F = Mday1_F(idx_day1_F,4);
    con_B = Mday1_B(idx_day1_B,4);
    con_all = [con_F(sort_F); con_B(sort_B)];
    
    % familiarity and value
    Mday1_FP = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_PE_DRE.csv']);
    
    for i = 1:length(id_F)
        idx_day1_F_FP(i) = find(id_F(i) == Mday1_FP(:,2)); %#ok<*AGROW>
        idx_day1_B_FP(i) = find(id_B(i) == Mday1_FP(:,2));
    end
    
    % familiarity
    fam_F = Mday1_FP(idx_day1_F_FP,3);
    fam_B = Mday1_FP(idx_day1_B_FP,3);
    fam_all = [fam_F(sort_F); fam_B(sort_B)];
    
    % monetary value
    pri_F = Mday1_FP(idx_day1_F_FP,4);
    pri_B = Mday1_FP(idx_day1_B_FP,4);
    pri_all = [pri_F(sort_F); pri_B(sort_B)];
    
    % create continuous RDMs
    for i = 1:length(val_all)
        RDMs{s}.val(:,i) = abs(val_all(i) - val_all)/50;
        RDMs{s}.con(:,i) = abs(con_all(i) - con_all)/50;
        RDMs{s}.fam(:,i) = abs(fam_all(i) - fam_all)/50;
        RDMs{s}.pri(:,i) = abs(pri_all(i) - pri_all);
    end
    
    % median splits
    val_all_pert = val_all + 0.00001*[1:length(val_all)]';
    medVal = median(val_all_pert);
    idxVal_H = val_all_pert > medVal;
    
    con_all_pert = con_all + 0.00001*[1:length(con_all)]';
    medCon = median(con_all_pert);
    idxCon_H = con_all_pert > medCon;
    
    fam_all_pert = fam_all + 0.00001*[1:length(fam_all)]';
    medFam = median(fam_all_pert);
    idxFam_H = fam_all_pert > medFam;
    
    pri_all_pert = pri_all + 0.00001*[1:length(pri_all)]';
    medPri = median(pri_all_pert);
    idxPri_H = pri_all_pert > medPri;
    
    % create medain split RDMs
    for i = 1:length(val_all)
        RDMs{s}.valMed(:,i) = abs(idxVal_H(i) - idxVal_H);
        RDMs{s}.conMed(:,i) = abs(idxCon_H(i) - idxCon_H);
        RDMs{s}.famMed(:,i) = abs(idxFam_H(i) - idxFam_H);
        RDMs{s}.priMed(:,i) = abs(idxPri_H(i) - idxPri_H);
    end
    
end