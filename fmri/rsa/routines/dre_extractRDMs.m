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
clear foo objs

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
        
        %% imagination
        % extract IDs
        sessIma(:,r) = Mday2(Mday2(:,3) == 0,4);
        
        %% choice
        % extract IDs (for both objects)
        sessCho_ID{r} = Mday2(Mday2(:,3) == 1,4:5);
        
        % extract choice (-1=L; 1=R)
        sessCho_LR(:,r) = Mday2(Mday2(:,3) == 1,7);
        
    end
    
    %% separate F and B
    if taskOrd(s) == 1
        sessIma_F_vec = [sessIma(:,1); sessIma(:,3)];
        sessIma_B_vec = [sessIma(:,2); sessIma(:,4)];
        sessIma_F = sessIma(:,[1 3]);
        sessIma_B = sessIma(:,[2 4]);
        sessCho_F = [sessCho_ID{1},sessCho_LR(:,1);sessCho_ID{3},sessCho_LR(:,3)];
        sessCho_B = [sessCho_ID{2},sessCho_LR(:,2);sessCho_ID{4},sessCho_LR(:,4)];
    else
        sessIma_F_vec = [sessIma(:,2); sessIma(:,4)];
        sessIma_B_vec = [sessIma(:,1); sessIma(:,3)];
        sessIma_F = sessIma(:,[2 4]);
        sessIma_B = sessIma(:,[1 3]);
        sessCho_F = [sessCho_ID{2},sessCho_LR(:,2);sessCho_ID{4},sessCho_LR(:,4)];
        sessCho_B = [sessCho_ID{1},sessCho_LR(:,1);sessCho_ID{3},sessCho_LR(:,3)];
    end
    
    % sort them
    [obj_F_sort,idx_sort_F] = sort(sessIma_F_vec);
    [obj_B_sort,idx_sort_B] = sort(sessIma_B_vec);
    
    %% load behavioural data from day 1
    dirBeha = [dirData,fs,'behaviour',fs,'SF',num2str(subs(s),'%03d')];
    Mday1_F = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_B1_DRE.csv']);
    Mday1_B = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_B2_DRE.csv']);
    Mday1_p = csvread([dirBeha,fs,'SF',num2str(subs(s),'%03d'),'_PE_DRE.csv']);
    
    %% assign scores to imagination trials
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % value and confidence %
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % for each object presented on day 2, find corresponding idx on day 1
    idx_day1_F = NaN(length(sessIma_F_vec),1);
    idx_day1_B = NaN(length(sessIma_B_vec),1);
    for i = 1:length(sessIma_F_vec)
        idx_day1_F(i) = find(sessIma_F_vec(i) == Mday1_F(:,2));
        idx_day1_B(i) = find(sessIma_B_vec(i) == Mday1_B(:,2));
    end
    
    % extract value, confidence
    val_F = Mday1_F(idx_day1_F,3);
    val_B = Mday1_B(idx_day1_B,3);
    val_F = val_F(idx_sort_F);
    val_B = val_B(idx_sort_B);
    val_all = [val_F(idx_sort_F); val_B(idx_sort_B)];
    con_F = Mday1_F(idx_day1_F,4);
    con_B = Mday1_B(idx_day1_B,4);
    con_all = [con_F(idx_sort_F); con_B(idx_sort_B)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % familiarity and price %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:length(sessIma_F_vec)
        idx_day1_F_FP(i) = find(sessIma_F_vec(i) == Mday1_p(:,2)); %#ok<*AGROW>
        idx_day1_B_FP(i) = find(sessIma_B_vec(i) == Mday1_p(:,2));
    end
    
    % extract familiarity, price
    fam_F = Mday1_p(idx_day1_F_FP,3);
    fam_B = Mday1_p(idx_day1_B_FP,3);
    fam_all = [fam_F(idx_sort_F); fam_B(idx_sort_B)];
    pri_F = Mday1_p(idx_day1_F_FP,4);
    pri_B = Mday1_p(idx_day1_B_FP,4);
    pri_all = [pri_F(idx_sort_F); pri_B(idx_sort_B)];
    
    %% create continuous RDMs
    % all trials
    for i = 1:length(val_all)
        RDMs{s}.val(:,i) = abs(val_all(i) - val_all);
        RDMs{s}.valF(:,i) = abs(val_all(i) - val_all);
        RDMs{s}.valB(:,i) = abs(val_all(i) - val_all);
        RDMs{s}.con(:,i) = abs(con_all(i) - con_all);
        RDMs{s}.fam(:,i) = abs(fam_all(i) - fam_all);
        RDMs{s}.pri(:,i) = abs(pri_all(i) - pri_all);
    end
    
    for i = 1:length(val_F)
        RDMs{s}.val_F(:,i) = abs(val_F(i) - val_F);
        RDMs{s}.val_B(:,i) = abs(val_B(i) - val_B);
    end
    
    zval_all = val_all;
    zval_all(~isnan(val_all)) = zscore(val_all(~isnan(val_all)));
    zcalQuad_all = zval_all.^2;
    for i = 1:length(val_all)
        RDMs{s}.valQuad(:,i) = abs(zcalQuad_all(i) - zcalQuad_all);
    end
    
    %% create median split RDMs
    % the 0.00001*idx is to make median univocally defined
    
    %     val_all_pert = val_all + 0.00001*(1:length(val_all))';
    %     medVal = nanmedian(val_all_pert);
    %     idxVal_H = val_all_pert > medVal;
    %
    %     con_all_pert = con_all + 0.00001*(1:length(con_all))';
    %     medCon = median(con_all_pert);
    %     idxCon_H = con_all_pert > medCon;
    %
    %     fam_all_pert = fam_all + 0.00001*(1:length(fam_all))';
    %     medFam = median(fam_all_pert);
    %     idxFam_H = fam_all_pert > medFam;
    %
    %     pri_all_pert = pri_all + 0.00001*(1:length(pri_all))';
    %     medPri = median(pri_all_pert);
    %     idxPri_H = pri_all_pert > medPri;
    %
    %     % create median split RDMs
    %     for i = 1:length(val_all)
    %         RDMs{s}.valMed(:,i) = abs(idxVal_H(i) - idxVal_H);
    %         RDMs{s}.conMed(:,i) = abs(idxCon_H(i) - idxCon_H);
    %         RDMs{s}.famMed(:,i) = abs(idxFam_H(i) - idxFam_H);
    %         RDMs{s}.priMed(:,i) = abs(idxPri_H(i) - idxPri_H);
    %     end
    
    %% now divide scores into high, low and create RDMs
    %     valHigh = val_all;
    %     valHigh(val_all_pert < medVal) = nan;
    %     valLow = val_all;
    %     valLow(val_all_pert > medVal) = nan;
    %     conHigh = con_all;
    %     conHigh(con_all_pert < medCon) = nan;
    %     conLow = con_all;
    %     conLow(con_all_pert > medCon) = nan;
    %     famHigh = fam_all;
    %     famHigh(fam_all_pert < medFam) = nan;
    %     famLow = fam_all;
    %     famLow(fam_all_pert > medFam) = nan;
    %
    %     for i = 1:length(val_all)
    %         RDMs{s}.valHigh(:,i) = abs(valHigh(i) - valHigh);
    %         RDMs{s}.valLow(:,i) = abs(valLow(i) - valLow);
    %         RDMs{s}.conHigh(:,i) = abs(conHigh(i) - conHigh);
    %         RDMs{s}.conLow(:,i) = abs(conLow(i) - conLow);
    %         RDMs{s}.famHigh(:,i) = abs(famHigh(i) - famHigh);
    %         RDMs{s}.famLow(:,i) = abs(famLow(i) - famLow);
    %     end
    
    %% create context model
    cxt_model = nan(240);
    
    for i = 1:120
        for j = 1:120
            
            % fire
            [~,sess_F_i] = find(sessIma_F == obj_F_sort(i));
            [~,sess_F_j] = find(sessIma_F == obj_F_sort(j));
            
            % boat
            [~,sess_B_i] = find(sessIma_B == obj_B_sort(i));
            [~,sess_B_j] = find(sessIma_B == obj_B_sort(j));
            
            % dissimilarity is 0 between different sessions with same context
            if sess_F_i ~= sess_F_j
                cxt_model(i,j) = 0;
            end
            
            if sess_B_i ~= sess_B_j
                cxt_model(120+i,120+j) = 0;
            end
            
            % dissimilarity is 1 across contexts (between S1,S4 and S2,S3 only)
            if sess_F_i ~= sess_B_j
                cxt_model(i,120+j) = 1;
            end
            
            if sess_F_j ~= sess_B_i
                cxt_model(120+i,j) = 1;
            end
        end
    end
    RDMs{s}.cxt = cxt_model;
    
    %% choice
    
    % assign scores to choice trials
    for r = 1:4
        
        numChoices = 48;
        
        % extract value of chosen and unchosen item
        valCho_F = NaN(numChoices,1);
        valUnc_F = NaN(numChoices,1);
        valCho_B = NaN(numChoices,1);
        valUnc_B = NaN(numChoices,1);
        for i = 1:numChoices
            idxCho_day1_F_L(i) = find(Mday1_F(:,2) == sessCho_F(i,1));
            idxCho_day1_F_R(i) = find(Mday1_F(:,2) == sessCho_F(i,2));
            idxCho_day1_B_L(i) = find(Mday1_B(:,2) == sessCho_B(i,1));
            idxCho_day1_B_R(i) = find(Mday1_B(:,2) == sessCho_B(i,2));
            if sessCho_F(i,3) == -1
                valCho_F(i) = Mday1_F(idxCho_day1_F_L(i),3);
                valUnc_F(i) = Mday1_F(idxCho_day1_F_R(i),3);
            elseif sessCho_F(i,3) == 1
                valCho_F(i) = Mday1_F(idxCho_day1_F_R(i),3);
                valUnc_F(i) = Mday1_F(idxCho_day1_F_L(i),3);
            end
            if sessCho_B(i,3) == -1
                valCho_B(i) = Mday1_B(idxCho_day1_B_L(i),3);
                valUnc_B(i) = Mday1_B(idxCho_day1_B_R(i),3);
            elseif sessCho_B(i,3) == 1
                valCho_B(i) = Mday1_B(idxCho_day1_B_R(i),3);
                valUnc_B(i) = Mday1_B(idxCho_day1_B_L(i),3);
            end
        end
        
        % put F and B together
        valCho = [valCho_F; valCho_B];
        valUnc = [valUnc_F; valUnc_B];
        
        % compute (signed and unsigned) value difference
        chMunc = valCho - valUnc;
        dVal = abs(valCho - valUnc);
        
        
        % movement onset and side
        movCho = Mday2(Mday2(:,3)==1,8);
        sidCho = Mday2(Mday2(:,3)==1,7);
        
    end
    
    % create context model for choice
    cxt_model_choice = [nan(24),zeros(24),nan(24),ones(24)
        zeros(24),nan(24),ones(24),nan(24)
        nan(24),ones(24),nan(24),zeros(24)
        ones(24),nan(24),zeros(24),nan(24)];
    
    %% create RDMs for choice trials
    % unlike imagination trials, I keep in task order within sessions (FFBB),
    % because there's no natural ordering and choice trials differed across subjects
    for i = 1:length(valCho)
        RDMs{s}.choice.dVal(:,i) = abs(dVal(i) - dVal)/50;
        RDMs{s}.choice.cMun(:,i) = abs(chMunc(i) - chMunc)/50;
        RDMs{s}.choice.Chos(:,i) = abs(valCho(i) - valCho)/50;
        RDMs{s}.choice.Unch(:,i) = abs(valUnc(i) - valUnc)/50;
        RDMs{s}.choice.ccxt = cxt_model_choice;
    end
    
end