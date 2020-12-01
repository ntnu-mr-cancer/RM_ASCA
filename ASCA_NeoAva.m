%% Importing data
addpath(genpath('/mnt/work/RM_ASCA+'))

opts = spreadsheetImportOptions("NumVariables", 20);
opts.Sheet = "Supplementary Table2";
opts.DataRange = "A15:T284";
opts.VariableNames = ["PatientNumber", "TimePoint", "Bevacizumabrandomizationgroup", "Pathologicalminimalresidualdiseaseresponse", "Glucose", "Ascorbate", "Lactate", "Tyrosine", "Myoinositol", "Glycine", "Taurine", "Glycerophosphocholine", "Phosphocholine", "Choline", "Creatine", "Glutathione", "Glutamine", "Succinate", "Glutamate", "Alanine"];
opts.VariableTypes = ["double", "categorical", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, ["TimePoint", "Bevacizumabrandomizationgroup", "Pathologicalminimalresidualdiseaseresponse"], "EmptyFieldRule", "auto");
MOESM3ESM = readtable("/mnt/work/RM_ASCA+/NeoAva/Data/11306_2017_1168_MOESM3_ESM.xlsx", opts, "UseExcel", false);
clear opts
data = MOESM3ESM; clear MOESM3ESM

data.Properties.VariableNames(12) = "GPC";
data.Properties.VariableNames(13) = "PC";

% Renaming variables and extracting metabolite names
data.Properties.VariableNames = ["ID", "timepoint", "treatment", "response", data.Properties.VariableNames(5:end)];
abbreviations = data.Properties.VariableNames(5:end);

% Sort data
data = sortrows(data,'ID','ascend');

%% Make design matrix X
T = dummyvar(data.timepoint);
T(:,1) = []; % Baseline as reference

D = dummyvar(data.treatment);
D(:,2) = []; % Control as reference

I = [T(:,1).*D, T(:,2).*D]; % Make interaction columns

X = [ones(size(data,1),1), T, I];
Z = [ones(size(data,1),1)]; % Random intercept
Y = data{:,5:end}; % Response matrix
G = data.ID; % Subject ID-vector

scalingfactor = std(Y(data.timepoint=="TP1",:)); % Scale to baseline standard deviation

A = [2:3]; % Column numbers for the time effect
B = [4:size(X,2)]; % Column numbers for the treatment+time*treatment effect
C = [2:size(X,2)]; % Column numbers for all fixed effects except intercept

%% Making effect matrices

clear b e rEffects
for i = 1:size(Y,2)
    lme = fitlmematrix(X,Y(:,i)./scalingfactor(i), Z, G);
    b(:,i) = lme.Coefficients.Estimate;
    rEffects(:,i) = randomEffects(lme);
    E(:,i) = residuals(lme);
end

M_A = X(:,A)*b(A,:); % Making effect matrix for time
M_B = X(:,B)*b(B,:); % Making effect matrix for treatment-time interactions
M_C = X(:,C)*b(C,:); % Making combined effect matrix

%% Extract scores for time, treatment and time+treatment effect matrices
[loadings_time, scores_time, eigen_time] = pca(M_A);
[loadings_treatment, scores_treatment, eigen_treatment] = pca(M_B);
[loadings_time_treatment, scores_time_treatment, eigen_time_treatment] = pca(M_C);

dummy_z = full(designMatrix(lme, 'random'));
randomeffects = dummy_z*rEffects;

aug_scores_time_treatment = M_C + randomeffects + E;
aug_scores_time_treatment = (aug_scores_time_treatment-mean(aug_scores_time_treatment))/loadings_time_treatment';

aug_scores_time_treatment2 = M_C + randomeffects;
aug_scores_time_treatment2 = (aug_scores_time_treatment2-mean(aug_scores_time_treatment2))/loadings_time_treatment';

% Finding scores to plot for Time
timepoints = sort(unique(data.timepoint), 'ascend');
for i = 1:length(unique(timepoints))
    PC1_scores_time(:,i) = scores_time(find(data.timepoint==timepoints(i),1,'first'), 1);
    PC2_scores_time(:,i) = scores_time(find(data.timepoint==timepoints(i),1,'first'), 2);
end

% Finding score-values to plot for Treatment*Time
treatments = sort(unique(data.treatment), 'ascend');
for e = 1:length(treatments)
    for i = 1:length(timepoints)
       PC1_scores_treatments(e,i) =  scores_treatment(find(data.timepoint==timepoints(i) & data.treatment==treatments(e), 1, 'first'), 1);
       PC2_scores_treatments(e,i) =  scores_treatment(find(data.timepoint==timepoints(i) & data.treatment==treatments(e), 1, 'first'), 2);
    end
end

% Finding score values to plot for Time + Time*Treatment
for e = 1:length(treatments)
    for i = 1:length(timepoints)
       PC1_scores_time_treatments(e,i) =  scores_time_treatment(find(data.timepoint==timepoints(i) & data.treatment==treatments(e), 1, 'first'), 1);
       PC2_scores_time_treatments(e,i) =  scores_time_treatment(find(data.timepoint==timepoints(i) & data.treatment==treatments(e), 1, 'first'), 2);
    end
end

%% Jackknife-validation
% Leave n patients out cross validation
allpatients = unique(data.ID, 'stable'); % Create vector containing all unique patient IDs
for i = 1:size(allpatients,1)
    groups_allpatients(i) = unique(data.treatment(find(data.ID==allpatients(i)))); % Find group membership for each patient
end

clear b lme
iterations = 50;
for i = 1:iterations;
    patients = cvpartition(groups_allpatients, "kfold", 7); % Partition patients into k subsets at random, keeping relative group sizes constant
    idx = training(patients, 1); % Create vector containing the indexes of the patients to include in training set. k-1 subsets are used.
    patients_to_include = allpatients(idx); % Make vector containing patient IDs of included patients
    idx = ismember(data.ID, patients_to_include); % Create selection vector for rows in data belonging to the patients in the training set.
    
    % Making training set
    X2 = X(idx,:);
    Y2 = Y(idx,:);
    Z2 = Z(idx,:);
    G2 = G(idx,:);
    scalingfactor = std(Y2(data.timepoint(idx)=="TP1", :)); % Making scaling factor for metabolites
    data2 = data(idx,:);
    
    for e = 1:size(Y2,2) % Calculate new effects based on training set
        lme = fitlmematrix(X2, Y2(:,e)./scalingfactor(e),Z2,G2);
        b(:,e) = lme.Coefficients.Estimate;
    end
    
     % Making new effect matrix for Time
    M_time2 = X2(:,A)*b(A,:); 
    [loadings2_time, scores2_time, eigen2_time] = pca(M_time2);
    
    % Rotate new loadings toward original loadings, and collect the values
    % in the variable cvloadings
    [rotated_loadings_time, rotation_matrix_time] = rotatefactors(loadings2_time(:,1:2),'Method','procrustes','type', 'orthogonal', 'Target',loadings_time(:,1:2)); 
    cvloadings_PC1_time(i,:) = rotated_loadings_time(:,1)';
    cvloadings_PC2_time(i,:) = rotated_loadings_time(:,2)';
    
    % Rotating the scores using the rotation matrix obtained from loading
    % rotation
    rotatedscores = scores2_time(:,1:2)*inv(rotation_matrix_time);
    
    % Finding the score value for each timepoint, and collect them in
    % vectors.
    for a = 1:length(unique(timepoints))
        PC1_scores2_time(:,a) = rotatedscores(find(data2.timepoint==timepoints(a), 1, 'first'), 1);
        PC2_scores2_time(:,a) = rotatedscores(find(data2.timepoint==timepoints(a), 1, 'first'), 2);
    end
    
    % Collect new scores in matrices. Every row contains the estimates from
    % one iteration.
    cvscores_PC1_time(i,:) = PC1_scores2_time;
    cvscores_PC2_time(i,:) = PC2_scores2_time;
    
    % Making new effect matrix and PCA for time-treatment interaction
    M_treatment2 = X2(:,B)*b(B,:); 
    [loadings2_treatment, scores2_treatment, eigen2_treatment] = pca(M_treatment2);
    
    % Rotate new loadings toward original loadings, and calculate the difference 
    [rotated_loadings_treatment, rotation_matrix_treatment] = rotatefactors(loadings2_treatment(:,1:2),'Method','procrustes','type', 'orthogonal', 'Target',loadings_treatment(:,1:2)); 
    cvloadings_PC1_treatment(i,:) = rotated_loadings_treatment(:,1)';
    cvloadings_PC2_treatment(i,:) = rotated_loadings_treatment(:,2)';
    
    % Rotate scores
    rotatedscores = scores2_treatment(:,1:2)*inv(rotation_matrix_treatment);
    
    for f = 1:length(treatments) % Loop for collecting scores to plot in a vector
        for g = 1:length(timepoints)
            PC1_scores2_treatments(f,g) =  rotatedscores(find(data2.timepoint==timepoints(g) & data2.treatment==treatments(f), 1, 'first'), 1);
            PC2_scores2_treatments(f,g) =  rotatedscores(find(data2.timepoint==timepoints(g) & data2.treatment==treatments(f), 1, 'first'), 2);
        end
    end
    
    % Collect new scores in matrices. Every row contains the estimates
    % from one iteration.
    cvscores_PC1_treatment(i,:) = PC1_scores2_treatments(1,:);
    cvscores_PC1_control(i,:) = PC1_scores2_treatments(2,:);
    
    cvscores_PC2_treatment(i,:) = PC2_scores2_treatments(1,:);
    cvscores_PC2_control(i,:) = PC2_scores2_treatments(2,:);
    
    % Making new effect matrix and PCA for time+time-treatment interaction
    M_time_treatment2 = X2(:,C)*b(C,:); 
    [loadings2_time_treatment, scores2_time_treatment, eigen2_time_treatment] = pca(M_time_treatment2);
    
    % Rotate new loadings toward original loadings, and calculate the difference 
    [rotated_loadings_time_treatment, rotation_matrix_time_treatment] = rotatefactors(loadings2_time_treatment(:,1:2),'Method','procrustes','type', 'orthogonal', 'Target',loadings_time_treatment(:,1:2)); 
    cvloadings_PC1_time_treatment(i,:) = rotated_loadings_time_treatment(:,1)';
    cvloadings_PC2_time_treatment(i,:) = rotated_loadings_time_treatment(:,2)';
    
    % Rotate scores
    rotatedscores = scores2_time_treatment(:,1:2)*inv(rotation_matrix_time_treatment);
    
    for f = 1:length(treatments) % Loop for collecting scores to plot in a vector
        for g = 1:length(timepoints)
            PC1_scores2_time_treatments(f,g) =  rotatedscores(find(data2.timepoint==timepoints(g) & data2.treatment==treatments(f), 1, 'first'), 1);
            PC2_scores2_time_treatments(f,g) =  rotatedscores(find(data2.timepoint==timepoints(g) & data2.treatment==treatments(f), 1, 'first'), 2);
        end
    end
    
    % Collect new scores in matrices. Every row contains the estimates
    % from one iteration.
    cvscores_PC1_time_treatment(i,:) = PC1_scores2_time_treatments(1,:);
    cvscores_PC1_time_control(i,:) = PC1_scores2_time_treatments(2,:);
    
    cvscores_PC2_time_treatment(i,:) = PC2_scores2_time_treatments(1,:);
    cvscores_PC2_time_control(i,:) = PC2_scores2_time_treatments(2,:);
end

% Find 2.5th and 97.5th percentiles of scores for the time effect
percentiles_cvscores_PC1_time = prctile(cvscores_PC1_time, [2.5, 97.5]); % First row contains 2.5th percentiles, second row contains 97.5th percentiles.
percentiles_cvscores_PC2_time = prctile(cvscores_PC2_time, [2.5, 97.5]);

% Find 2.5th and 97.5th percentiles of the loadings for the time effect
percentiles_cvloadings_PC1_time = prctile(cvloadings_PC1_time, [2.5, 97.5]);
percentiles_cvloadings_PC2_time = prctile(cvloadings_PC2_time, [2.5, 97.5]);

% Find 2.5th and 97.5th percentiles of scores for time*treatment effect for
% each of the groups.
lower_PC1_treatment = [prctile(cvscores_PC1_treatment, 2.5); prctile(cvscores_PC1_control, 2.5)];
upper_PC1_treatment = [prctile(cvscores_PC1_treatment, 97.5); prctile(cvscores_PC1_control, 97.5)];

lower_PC2_treatment = [prctile(cvscores_PC2_treatment, 2.5); prctile(cvscores_PC2_control, 2.5)];
upper_PC2_treatment = [prctile(cvscores_PC2_treatment, 97.5); prctile(cvscores_PC2_control, 97.5)];

% Find 2.5th and 97.5th percentiles of the loadings for the time*treatment
% effect
percentiles_cvloadings_PC1_treatment = prctile(cvloadings_PC1_treatment, [2.5, 97.5]);
percentiles_cvloadings_PC2_treatment = prctile(cvloadings_PC2_treatment, [2.5, 97.5]);

% Find 2.5th and 97.5th percentiles of scores for time+time*treatment effect for
% each of the groups.
lower_PC1_time_treatment = [prctile(cvscores_PC1_time_treatment, 2.5); prctile(cvscores_PC1_time_control, 2.5)];
upper_PC1_time_treatment = [prctile(cvscores_PC1_time_treatment, 97.5); prctile(cvscores_PC1_time_control, 97.5)];

lower_PC2_time_treatment = [prctile(cvscores_PC2_time_treatment, 2.5); prctile(cvscores_PC2_time_control, 2.5)];
upper_PC2_time_treatment = [prctile(cvscores_PC2_time_treatment, 97.5); prctile(cvscores_PC2_time_control, 97.5)];

% Find 2.5th and 97.5th percentiles of the loadings for the time*treatment
% effect
percentiles_cvloadings_PC1_time_treatment = prctile(cvloadings_PC1_time_treatment, [2.5, 97.5]);
percentiles_cvloadings_PC2_time_treatment = prctile(cvloadings_PC2_time_treatment, [2.5, 97.5]);

data.timepoint = string(data.timepoint);
data.timepoint(data.timepoint=="TP1") = "1";
data.timepoint(data.timepoint=="TP2") = "2";
data.timepoint(data.timepoint=="TP3") = "3";

%% Time effect
% Scree plot
figure
tiledlayout(2,3, 'TileSpacing',"none")
nexttile([2 1])
bar(eigen_time(1:5)./sum(eigen_time)*100)
ylim([0,100])
xlabel("Principal components", 'FontSize',8)
ylabel("Explained variance (%)", 'FontSize',8)
title("Scree plot", 'FontSize',8)
grid on

%  Plotting PC1-scores
nexttile
hold on 
errorbar([0, 1, 2], PC1_scores_time, PC1_scores_time-percentiles_cvscores_PC1_time(1,:), percentiles_cvscores_PC1_time(2,:)-PC1_scores_time)
set(gca,'XTick',[0 1 2]);
ylabel("PC1 (" + num2str((eigen_time(1)/sum(eigen_time))*100, '%.2f') + "%)", "FontSize",8)
title("Scores", 'FontSize',8)
grid on
hold off

% Loading plot
nexttile
bar(1:size(Y,2), loadings_time(:,1)')
hold on
errorbar([1:size(Y,2)], loadings_time(:,1)', loadings_time(:,1)' - percentiles_cvloadings_PC1_time(1,:), percentiles_cvloadings_PC1_time(2,:) - loadings_time(:,1)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(Y,2)]);
set(gca,'XTickLabel',abbreviations, 'XTickLabelRotation', 90, 'FontSize', 4);
title("Loadings", 'FontSize',8)
grid on
hold off

%  Plotting PC2-scores
nexttile
hold on 
errorbar([0, 1, 2], PC2_scores_time, PC2_scores_time-percentiles_cvscores_PC2_time(1,:), percentiles_cvscores_PC2_time(2,:)-PC2_scores_time)
set(gca,'XTick',[0 1 2]);
xlabel("Months", 'FontSize', 8)
ylabel("PC2 (" + num2str((eigen_time(2)/sum(eigen_time))*100, '%.2f') + "%)","FontSize",8)
grid on
hold off

% Loading plot
nexttile
bar(1:size(Y,2), loadings_time(:,2)')
hold on
errorbar([1:size(Y,2)], loadings_time(:,2)', loadings_time(:,2)'-percentiles_cvloadings_PC2_time(1,:), percentiles_cvloadings_PC2_time(2,:)-loadings_time(:,2)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(Y,2)]);
set(gca,'XTickLabel',abbreviations, 'XTickLabelRotation', 90, 'FontSize', 4);
grid on
hold off

%% Time*Treatment interaction
% Making figure
figure
tiledlayout(2,3, 'TileSpacing',"none")
nexttile([2 1])
bar(eigen_treatment(1:5)./sum(eigen_treatment)*100)
ylim([0,100])
xlabel("Principal components", 'FontSize',8)
ylabel("Explained variance (%)", 'FontSize',8)
title("Scree plot", 'FontSize',8)
grid on

nexttile
hold on
title("Scores", "FontSize", 8)
ylabel("PC1 (" + num2str((eigen_treatment(1)/sum(eigen_treatment))*100, '%.2f') + " %)", 'FontSize', 8);
for i = 1:size(PC1_scores_treatments,1)
    errorbar([0 1 2], PC1_scores_treatments(i,:), PC1_scores_treatments(i,:) - lower_PC1_treatment(i,:), upper_PC1_treatment(i,:) - PC1_scores_treatments(i,:))
end
set(gca,'XTick',[0 1 2]);
legend({'CTX + B', 'CTX'}, 'FontSize', 4, 'Position', [0.5,2,1,1])
grid on
hold off

nexttile
bar(1:size(Y,2), loadings_treatment(:,1)')
hold on
errorbar([1:size(Y,2)], loadings_treatment(:,1)', loadings_treatment(:,1)' - percentiles_cvloadings_PC1_treatment(1,:), percentiles_cvloadings_PC1_treatment(2,:) - loadings_treatment(:,1)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(Y,2)]);
set(gca,'XTickLabel',abbreviations, 'XTickLabelRotation', 90, 'FontSize', 4);
title("Loadings", "FontSize", 8)
grid on
hold off

% Plotting results for PC2
nexttile
hold on
ylabel("PC2 (" + num2str((eigen_treatment(2)/sum(eigen_treatment))*100, '%.2f') + " %)", 'FontSize', 8);
for i = 1:size(PC2_scores_treatments,1)
    errorbar([0 1 2], PC2_scores_treatments(i,:), PC2_scores_treatments(i,:) - lower_PC2_treatment(i,:), upper_PC2_treatment(i,:) - PC2_scores_treatments(i,:))
end
set(gca,'XTick',[0 1 2]);
grid on
xlabel("Timepoint", 'FontSize',8)
hold off

nexttile
bar(1:size(Y,2), loadings_treatment(:,2)')
hold on
errorbar([1:size(Y,2)], loadings_treatment(:,2)', loadings_treatment(:,2)' - percentiles_cvloadings_PC2_treatment(1,:), percentiles_cvloadings_PC2_treatment(2,:) - loadings_treatment(:,2)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(Y,2)]);
set(gca,'XTickLabel', abbreviations, 'XTickLabelRotation', 90, 'FontSize', 4);
grid on
hold off

%% Time + Time*Treatment interaction
figure
tiledlayout(2,3, 'TileSpacing',"none")
nexttile([2 1])
bar(eigen_time_treatment(1:5)./sum(eigen_time_treatment)*100)
ylim([0,100])
xlabel("Principal components", 'FontSize',8)
ylabel("Explained variance (%)", 'FontSize',8)
title("Scree plot", 'FontSize',8)
grid on

nexttile
hold on
title("Scores", "FontSize", 8)
ylabel("PC1 (" + num2str((eigen_time_treatment(1)/sum(eigen_time_treatment))*100, '%.2f') + " %)", 'FontSize', 8);
for i = 1:size(PC1_scores_time_treatments,1)
    errorbar([0 1 2], PC1_scores_time_treatments(i,:), PC1_scores_time_treatments(i,:) - lower_PC1_time_treatment(i,:), upper_PC1_time_treatment(i,:) - PC1_scores_time_treatments(i,:))
end
set(gca,'XTick',[0 1 2]);
legend({'CTX + B', 'CTX'}, 'FontSize', 4, 'Position', [0.5,2,1,1])
grid on
hold off

nexttile
bar(1:size(Y,2), loadings_time_treatment(:,1))
hold on
errorbar([1:size(Y,2)], loadings_time_treatment(:,1), loadings_time_treatment(:,1)' - percentiles_cvloadings_PC1_time_treatment(1,:), percentiles_cvloadings_PC1_time_treatment(2,:) - loadings_time_treatment(:,1)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(Y,2)]);
set(gca,'XTickLabel',abbreviations, 'XTickLabelRotation', 90, 'FontSize', 4);
title("Loadings", "FontSize", 8)
grid on
hold off

% Plotting results for PC2
nexttile
hold on
ylabel("PC2 (" + num2str((eigen_time_treatment(2)/sum(eigen_time_treatment))*100, '%.2f') + " %)", 'FontSize', 8);
for i = 1:size(PC2_scores_time_treatments,1)
    errorbar([0 1 2], PC2_scores_time_treatments(i,:), PC2_scores_time_treatments(i,:) - lower_PC2_time_treatment(i,:), upper_PC2_time_treatment(i,:) - PC2_scores_time_treatments(i,:))
end
set(gca,'XTick',[0 1 2]);
grid on
xlabel("Timepoint", 'FontSize',8)
hold off

nexttile
bar(1:size(Y,2), loadings_time_treatment(:,2)')
hold on
errorbar([1:size(Y,2)], loadings_time_treatment(:,2)', loadings_time_treatment(:,2)' - percentiles_cvloadings_PC2_time_treatment(1,:), percentiles_cvloadings_PC2_time_treatment(2,:) - loadings_time_treatment(:,2)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(Y,2)]);
set(gca,'XTickLabel', abbreviations, 'XTickLabelRotation', 90, 'FontSize', 4);
grid on
hold off

%% Augmented effect matrix figure

% Find rows of patients with all three measurements available
for d = 1:size(allpatients)
    idx(d) = length(data.timepoint(data.ID==allpatients(d)));
end
idx = find(idx==3);

% Randomly select half of them
% idx = idx(randperm(46, 23));
% idx = sort(idx);
idx = [7     9    19    28    37    41    54    55    58    60    63    67    68    72    74    75    76    81    86    89    91    93   100];

figure
subplot(1,2,1)
hold on
axis('square')
title("PC1", "FontSize", 8)
ylabel("PC1 (" + num2str((eigen_time_treatment(1)/sum(eigen_time_treatment))*100, '%.2f') + " %)", 'FontSize', 8);

colors = ['r', 'c'];
for d = 1:size(idx, 2)
    plot(str2double(data.timepoint(data.ID==idx(d))), aug_scores_time_treatment2(data.ID==idx(d), 1), 'color', colors(double(unique(data.treatment(data.ID==idx(d)))=="Bevacizumab treated")+1))
end

for d = 1:size(idx, 2)
    plot(str2double(data.timepoint(data.ID==idx(d))), aug_scores_time_treatment(data.ID==idx(d), 1), '--', 'LineWidth', 0.2, 'color', colors(double(unique(data.treatment(data.ID==allpatients(d)))=="Bevacizumab treated")+1))
end

set(gca,'XTick',[1 2 3]);
grid on
xlabel("Timepoint", 'FontSize',8)
legend([plot(str2double(data.timepoint(data.ID==allpatients(3))), aug_scores_time_treatment2(data.ID==allpatients(3), 1), 'color', colors(double(unique(data.treatment(data.ID==allpatients(3)))=="Bevacizumab treated")+1)),...
    plot(str2double(data.timepoint(data.ID==allpatients(1))), aug_scores_time_treatment2(data.ID==allpatients(1), 1), 'color', colors(double(unique(data.treatment(data.ID==allpatients(1)))=="Bevacizumab treated")+1))],...
    {"CTX + B", "CTX"})
hold off

subplot(1,2,2)
hold on
axis('square')
title("PC2", "FontSize", 8)
ylabel("PC2 (" + num2str((eigen_time_treatment(2)/sum(eigen_time_treatment))*100, '%.2f') + " %)", 'FontSize', 8);

for d = 1:size(idx, 2)
    plot(str2double(data.timepoint(data.ID==idx(d))), aug_scores_time_treatment2(data.ID==idx(d), 2), 'color', colors(double(unique(data.treatment(data.ID==idx(d)))=="Bevacizumab treated")+1))
end

for d = 1:size(idx, 2)
    plot(str2double(data.timepoint(data.ID==idx(d))), aug_scores_time_treatment(data.ID==idx(d), 2), '--', 'LineWidth', 0.2, 'color', colors(double(unique(data.treatment(data.ID==allpatients(d)))=="Bevacizumab treated")+1))
end
set(gca,'XTick',[1 2 3]);
grid on
xlabel("Timepoint", 'FontSize',8)
hold off
