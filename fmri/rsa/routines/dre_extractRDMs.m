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

dirData = [dir.dre,fs,'data'];

%% read objects
objs        = readtable([dir.behDat,fs,'Objects.csv']);
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
        
        %% load behavioural data day 2
        dirPsyO = [dirData,fs,'fmri',fs,'psychOut',fs,'SF',num2str(subs(s),'%03d')];
        Mday2 = csvread([dirPsyO,fs,'DRE_mri_S',num2str(subs(s),'%03d'),'_B',num2str(r),day2Order{r},'.csv']);
        
        %% object IDs
        % extract IDs from day 2
        sessions(:,r) = Mday2(Mday2(:,3)==0,4);
        id_vec = sessions(:);
        
    end
    
    %% separate F and B
    if taskOrd(s) == 1
        obj_F = [sessions(:,1); sessions(:,3)];
        obj_B = [sessions(:,2); sessions(:,4)];
    else
        obj_F = [sessions(:,2); sessions(:,4)];
        obj_B = [sessions(:,1); sessions(:,3)];
    end
    
    % sort them
    [obj_F_sort,idx_sort_F] = sort(obj_F);
    [obj_B_sort,idx_sort_B] = sort(obj_B);
    
    %% day 1 - load behavioural data
    dirBeha = [dirData,fs,'behaviour',fs,'SF',num2str(subs(s),'%03d')];
    Mday1_F = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_B1_DRE.csv']);
    Mday1_B = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_B2_DRE.csv']);
    Mday1_p = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_PE_DRE.csv']);
    
    %% day 2 - compute scores assigned to presented objects
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % value and confidence %
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % for each object presented on day 2, find corresponding idx on day 1
    idx_day1_F = NaN(length(obj_F),1);
    idx_day1_B = NaN(length(obj_B),1);
    for i = 1:length(obj_F)
        idx_day1_F(i) = find(obj_F(i) == Mday1_F(:,2));
        idx_day1_B(i) = find(obj_B(i) == Mday1_B(:,2));
    end
    
    % extract value, confidence
    val_F = Mday1_F(idx_day1_F,3);
    val_B = Mday1_B(idx_day1_B,3);
    val_all = [val_F(idx_sort_F); val_B(idx_sort_B)];
    con_F = Mday1_F(idx_day1_F,4);
    con_B = Mday1_B(idx_day1_B,4);
    con_all = [con_F(idx_sort_F); con_B(idx_sort_B)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % familiarity and price %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:length(obj_F)
        idx_day1_F_FP(i) = find(obj_F(i) == Mday1_p(:,2)); %#ok<*AGROW>
        idx_day1_B_FP(i) = find(obj_B(i) == Mday1_p(:,2));
    end
    
    % extract familiarity, price
    fam_F = Mday1_p(idx_day1_F_FP,3);
    fam_B = Mday1_p(idx_day1_B_FP,3);
    fam_all = [fam_F(idx_sort_F); fam_B(idx_sort_B)];
    pri_F = Mday1_p(idx_day1_F_FP,4);
    pri_B = Mday1_p(idx_day1_B_FP,4);
    pri_all = [pri_F(idx_sort_F); pri_B(idx_sort_B)];
    
    %% create continuous RDMs
    for i = 1:length(val_all)
        RDMs{s}.val(:,i) = abs(val_all(i) - val_all)/50;
        RDMs{s}.con(:,i) = abs(con_all(i) - con_all)/50;
        RDMs{s}.fam(:,i) = abs(fam_all(i) - fam_all)/50;
        RDMs{s}.pri(:,i) = abs(pri_all(i) - pri_all);
    end
    
    %% create median split RDMs
    % the 0.00001*idx is to make median univocally defined
    
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
    
    % create median split RDMs
    for i = 1:length(val_all)
        RDMs{s}.valMed(:,i) = abs(idxVal_H(i) - idxVal_H);
        RDMs{s}.conMed(:,i) = abs(idxCon_H(i) - idxCon_H);
        RDMs{s}.famMed(:,i) = abs(idxFam_H(i) - idxFam_H);
        RDMs{s}.priMed(:,i) = abs(idxPri_H(i) - idxPri_H);
    end
    
    %% divide into three percentiles
    
    %% create context model
    % to avoid overestimation of correlation due to time correlation within
    % sessions, we set correlation within session to NaN. Dissimilarity
    % across sessions is zero and one for same and different contect, respectively.
    if taskOrd(s) == 1
        sessions_F = sessions(:,[1 3]);
        sessions_B = sessions(:,[2 4]);
    else
        sessions_F = sessions(:,[2 4]);
        sessions_B = sessions(:,[1 3]);
    end
    cxt_model = nan(240);
    
    for i = 1:120
        for j = 1:120
            
            % fire
            [~,sess_F_i] = find(sessions_F == obj_F_sort(i));
            [~,sess_F_j] = find(sessions_F == obj_F_sort(j));
            
            % boat
            [~,sess_B_i] = find(sessions_B == obj_B_sort(i));
            [~,sess_B_j] = find(sessions_B == obj_B_sort(j));
            
            if sess_F_i ~= sess_F_j
                cxt_model(i,j) = 0;
            end
            
            if sess_F_i ~= sess_B_j
                cxt_model(i,120+j) = 1;
            end
            
            if sess_F_j ~= sess_B_i
                cxt_model(120+i,j) = 1;
            end
            
            if sess_B_i ~= sess_B_j
                cxt_model(120+i,120+j) = 0;
            end
            
        end
    end
    RDMs{s}.cxt = cxt_model;
    
%     cxt_model(isnan(cxt_model)) = -1;
%     figure,imagesc(cxt_model)
    
end