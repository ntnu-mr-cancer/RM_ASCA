%% Importing data
addpath(genpath('/mnt/work/RM_ASCA+'))

opts = spreadsheetImportOptions("NumVariables", 9);
opts.Sheet = "Samples";
opts.DataRange = "A2:I108";
opts.VariableNames = ["Indiv", "NMRExperiment", "VarName3", "SURGERY", "Pre_operation", "months", "months1", "months2", "months3"];
opts.VariableTypes = ["double", "string", "string", "categorical", "string", "string", "string", "string", "string"];
opts = setvaropts(opts, ["NMRExperiment", "VarName3", "Pre_operation", "months", "months1", "months2", "months3"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["NMRExperiment", "VarName3", "SURGERY", "Pre_operation", "months", "months1", "months2", "months3"], "EmptyFieldRule", "auto");
samples = readtable("/mnt/work/EBBA-II/Data/NMR_data_bariatric.xlsx", opts, "UseExcel", false);
clear opts

% Importing Metabolite levels
opts = spreadsheetImportOptions("NumVariables", 22);
opts.Sheet = "Metabolite levels";
opts.DataRange = "A2:V466";
opts.VariableNames = ["VarName1", "hydroxybutyrate", "acetate", "acetoacetate", "alanine", "CHEBI16414", "CHEBI17368", "CID15012407", "creatinine", "glutamine", "glycine", "histidine", "isoleucine", "isopropylalcohol", "lactate", "leucine", "lipoproteins", "methanol", "Methylsulfonylmethane", "phenylalanine", "pyruvate", "tyrosine"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, "VarName1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "VarName1", "EmptyFieldRule", "auto");
metabolites = readtable("/mnt/work/EBBA-II/Data/NMR_data_bariatric.xlsx", opts, "UseExcel", false);
clear opts

% Make data table
samples.NMRExperiment = erase(samples.NMRExperiment, ["Obs", "0_", "s"]);
metabolites.VarName1 = char(metabolites.VarName1);
timepoint = [];
for i = 1:size(metabolites,1);
    if metabolites.VarName1(i) == '0'
        timepoint(i) = 0;
    elseif metabolites.VarName1(i) == '1'
        timepoint(i) = 1;
    elseif metabolites.VarName1(i) == '2'
        timepoint(i) = 2;
    elseif metabolites.VarName1(i) == '3'
        timepoint(i) = 3;
    elseif metabolites.VarName1(i) == '4'
        timepoint(i) = 4;
    end
end
timepoint = table(timepoint');
timepoint.Properties.VariableNames = "timepoint";
metabolites = [metabolites(:,1), timepoint, metabolites(:,2:end)];
metabolites.VarName1 = erase(string(metabolites.VarName1), ["_S", "-S", "_s"]);
metabolites.VarName1 = erase(string(metabolites.VarName1), ["0_", "1_", "2_", "3_", "4_", "0-", "1-", "2-", "3-", "4-"]);

% Select out patients with all timepoints available
select = metabolites(ismember(metabolites.VarName1, samples.NMRExperiment),:);

% Find patient IDs
for i = 1:size(select,1)
    ID(i) = samples.Indiv(samples.NMRExperiment == select.VarName1(i));
end

data = [table(ID'), select];
data.Properties.VariableNames(1) = "ID";
treatment =  [];
for i = 1:size(data,1);
    treatment(i) = samples.SURGERY(samples.Indiv == data.ID(i));
end

data = [data(:,1:2), table(treatment'), data(:,3:end)];
data.Properties.VariableNames(3) = "treatment";
data.VarName1 = [];

clear select metabolites ID i timepoint treatment i
clear baseline_matrix i samples baseline

%% Modifying variable names and types
% Classifying time and treatment variables as nominal
data.timepoint = nominal(data.timepoint);
data.treatment = nominal(data.treatment);
data.ID = nominal(data.ID);

% Changing variable names for three metabolites
data.Properties.VariableNames(8) = "valine";
data.Properties.VariableNames(9) = "hypoxanthine";
data.Properties.VariableNames(10) = "citrate";

% Making vector containing metabolite names
abbreviations = ["Hydroxybutyrate", "Acetate", "Acetoacetate", "Alanine", "Valine", "Hypoxanthine", "Citrate", "Creatine", "Glutamine", "Glycine", "Histidine", "Isocitrate", "Isopropylalcohol", "Lactate", "Leucine", "Lipoproteins", "Methanol", "MSM", "Phenylalanine", "Pyruvate", "Tyrosine"];

%% Creating design matrix for the unconstrained longitudinal data analysis model
T = dummyvar(data.timepoint);
T(:,1) = []; % Baseline as reference

D = [(data.treatment == "1") - (data.treatment == "3"), (data.treatment=="2") - (data.treatment=="3")]; % Sum coded treatment variable

I = [T(:,1).*D(:,1), T(:,2).*D(:,1), T(:,3).*D(:,1), T(:,4).*D(:,1), T(:,1).*D(:,2), T(:,2).*D(:,2), T(:,3).*D(:,2), T(:,4).*D(:,2)]; % Make interaction column

X = [ones(size(data,1),1), T, D, I];
Z = [ones(size(data,1),1)]; % Random intercept
Y = sqrt(data{:,4:end}); % Square root transformation for data
G = data.ID; % Grouping factor for random effects

%% Run mixed models on each response variable and create effect matrices
scalingfactor = std(Y(data.timepoint=="0",:));

% Estimate coefficients b for each metabolite
clear b e 
for i = 1:size(Y,2)
    lme = fitlmematrix(X,Y(:,i)./scalingfactor(i), Z, G);
    b(:,i) = lme.Coefficients.Estimate;
    rEffects(:,i) = randomEffects(lme);
    E(:,i) = residuals(lme);
end

% Specify which factors to include in the different effect matrices
A = [2:5]; % Time factor
B = [6:size(X,2)]; % Treatment factor and Time*Treatment interactions
C = [2:size(X,2)]; % Time + Treatment + Time*Treatment interactions

M_time = X(:,A)*b(A,:); % Making effect matrix for time
M_treatment = X(:,B)*b(B,:); % Making effect matrix for treatment-time interactions
M_time_treatment = X(:,C)*b(C,:); % Making effect matrix for treatment-time interactions

%% PCA on effect matrices
[loadings_time, scores_time, eigen_time] = pca(M_time);
[loadings_treatment, scores_treatment, eigen_treatment] = pca(M_treatment);
[loadings_time_treatment, scores_time_treatment, eigen_time_treatment] = pca(M_time_treatment);

% Finding scores to plot for Time factor
timepoints = sort(unique(data.timepoint), 'ascend');
for i = 1:length(unique(timepoints))
    PC1_scores_time(:,i) = scores_time(find(data.timepoint==timepoints(i), 1, 'first'), 1);
    PC2_scores_time(:,i) = scores_time(find(data.timepoint==timepoints(i), 1, 'first'), 2);
end

% Finding score-values to plot for Treatment + Treatment*Time
treatments = sort(unique(data.treatment), 'ascend');
for e = 1:length(treatments)
    for i = 1:length(timepoints)
       PC1_scores_treatments(e,i) =  scores_treatment(find(data.timepoint==timepoints(i) & data.treatment==treatments(e), 1, 'first'), 1);
       PC2_scores_treatments(e,i) =  scores_treatment(find(data.timepoint==timepoints(i) & data.treatment==treatments(e), 1, 'first'), 2);
    end
end

% Finding score-values to plot for Time + Treatment + Treatment*Time
treatments = sort(unique(data.treatment), 'ascend');
for e = 1:length(treatments)
    for i = 1:length(timepoints)
       PC1_scores_time_treatments(e,i) =  scores_time_treatment(find(data.timepoint==timepoints(i) & data.treatment==treatments(e), 1, 'first'), 1);
       PC2_scores_time_treatments(e,i) =  scores_time_treatment(find(data.timepoint==timepoints(i) & data.treatment==treatments(e), 1, 'first'), 2);
    end
end

%% Jackknife-validation
% 7-fold random subset jackknifing, 100 iterations. Only the training
% set is used to estmate variability in model parameters.

allpatients = unique(data.ID); % Create vector containing all unique patient IDs
for i = 1:size(allpatients,1)
    groups_allpatients(i) = unique(data.treatment(find(data.ID==allpatients(i)))); % Find group membership for each patient
end

iterations = 100;
for i = 1:iterations;
    patients = cvpartition(groups_allpatients, "kfold", 7); % Partition patients into k subsets at random, keeping relative group sizes constant
    idx = training(patients, 1); % Create vector containing the indexes of the patients to include in training set. k-1 subsets are used.
    patients_to_include = allpatients(idx); % Make vector containing patient IDs of included patients
    idx = ismember(data.ID, patients_to_include); % Create selection vector for rows in data belonging to the patients in the training set.
    
    % Making training set
    X2 = X(idx,:);
    Y2 = Y(idx,:);
    Z2 = Z(idx);
    G2 = G(idx);
    scalingfactor = Y(idx,:);
    scalingfactor = std(scalingfactor(data.timepoint(idx)=="0", :)); % Making scaling factor for metbolites
    data2 = data(idx,:);
    
    for e = 1:size(Y,2) % Calculate new effects based on training set
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
    rotated_scores = scores2_time(:,1:2)*inv(rotation_matrix_time);
    
    % Finding the score value for each timepoint, and collect them in
    % vectors.
    for a = 1:length(unique(timepoints))
        PC1_scores2_time(:,a) = rotated_scores(find(data2.timepoint==timepoints(a), 1, 'first'), 1);
        PC2_scores2_time(:,a) = rotated_scores(find(data2.timepoint==timepoints(a), 1, 'first'), 2);
    end
    
    % Collect new scores in matrices. Every row contains the estimates from
    % one iteration.
    cvscores_PC1_time(i,:) = PC1_scores2_time;
    cvscores_PC2_time(i,:) = PC2_scores2_time;
    
    % Making new effect matrix and PCA for time-treatment interaction
    M_treatment2 = X2(:,B)*b(B,:); 
    [loadings2_treatment, scores2_treatment, eigen2_treatment] = pca(M_treatment2);
    
    % Rotate new loadings toward original loadings
    [rotated_loadings_treatment, rotation_matrix_treatment] = rotatefactors(loadings2_treatment(:,1:2),'Method','procrustes','type', 'orthogonal', 'Target',loadings_treatment(:,1:2)); 
    cvloadings_PC1_treatment(i,:) = rotated_loadings_treatment(:,1)';
    cvloadings_PC2_treatment(i,:) = rotated_loadings_treatment(:,2)';
    
    % Rotate scores
    rotated_scores = scores2_treatment(:,1:2)*inv(rotation_matrix_treatment);
    
    for f = 1:length(treatments) % Loop for collecting scores to plot in a vector
        for g = 1:length(timepoints)
            PC1_scores2_treatments(f,g) =  rotated_scores(find(data2.timepoint==timepoints(g) & data2.treatment==treatments(f), 1, 'first'), 1);
            PC2_scores2_treatments(f,g) =  rotated_scores(find(data2.timepoint==timepoints(g) & data2.treatment==treatments(f), 1, 'first'), 2);
        end
    end
    
    % Collect new scores in matrices. Every row contains the estimates
    % from one iteration.
    cvscores_PC1_treatment_distal(i,:) = PC1_scores2_treatments(1,:);
    cvscores_PC1_treatment_proximal(i,:) = PC1_scores2_treatments(2,:);
    cvscores_PC1_treatment_sleeve(i,:) = PC1_scores2_treatments(3,:);
    
    cvscores_PC2_treatment_distal(i,:) = PC2_scores2_treatments(1,:);
    cvscores_PC2_treatment_proximal(i,:) = PC2_scores2_treatments(2,:);
    cvscores_PC2_treatment_sleeve(i,:) = PC2_scores2_treatments(3,:);
    
    % Making new effect matrix and PCA for time + treatment + time-treatment interaction
    M_time_treatment2 = X2(:,C)*b(C,:); 
    [loadings2_time_treatment, scores2_time_treatment, eigen2_time_treatment] = pca(M_time_treatment2);
    
    % Rotate new loadings toward original loadings, and calculate the difference 
    [rotated_loadings_time_treatment, rotation_matrix_time_treatment] = rotatefactors(loadings2_time_treatment(:,1:2),'Method','procrustes','type', 'orthogonal', 'Target',loadings_time_treatment(:,1:2)); 
    cvloadings_PC1_time_treatment(i,:) = rotated_loadings_time_treatment(:,1)';
    cvloadings_PC2_time_treatment(i,:) = rotated_loadings_time_treatment(:,2)';
    
    % Rotate scores
    rotated_scores = scores2_time_treatment(:,1:2)*inv(rotation_matrix_time_treatment);
    
    for f = 1:length(treatments) % Loop for collecting scores to plot in a vector
        for g = 1:length(timepoints)
            PC1_scores2_time_treatments(f,g) =  rotated_scores(find(data2.timepoint==timepoints(g) & data2.treatment==treatments(f), 1, 'first'), 1);
            PC2_scores2_time_treatments(f,g) =  rotated_scores(find(data2.timepoint==timepoints(g) & data2.treatment==treatments(f), 1, 'first'), 2);
        end
    end
    
    % Collect new scores in matrices. Every row contains the estimates
    % from one iteration.
    cvscores_PC1_time_treatment_distal(i,:) = PC1_scores2_time_treatments(1,:);
    cvscores_PC1_time_treatment_proximal(i,:) = PC1_scores2_time_treatments(2,:);
    cvscores_PC1_time_treatment_sleeve(i,:) = PC1_scores2_time_treatments(3,:);
    
    cvscores_PC2_time_treatment_distal(i,:) = PC2_scores2_time_treatments(1,:);
    cvscores_PC2_time_treatment_proximal(i,:) = PC2_scores2_time_treatments(2,:);
    cvscores_PC2_time_treatment_sleeve(i,:) = PC2_scores2_time_treatments(3,:);
end

% Find 2.5th and 97.5th percentiles of scores for the time effect
percentiles_cvscores_PC1_time = prctile(cvscores_PC1_time, [2.5, 97.5]); % First row contains 2.5th percentiles, second row contains 97.5th percentiles.
percentiles_cvscores_PC2_time = prctile(cvscores_PC2_time, [2.5, 97.5]);

% Find 2.5th and 97.5th percentiles of the loadings for the time effect
percentiles_cvloadings_PC1_time = prctile(cvloadings_PC1_time, [2.5, 97.5]);
percentiles_cvloadings_PC2_time = prctile(cvloadings_PC2_time, [2.5, 97.5]);

% Find 2.5th and 97.5th percentiles of scores for time*treatment effect for
% each of the groups.
lower_PC1_treatment = [prctile(cvscores_PC1_treatment_distal, 2.5); prctile(cvscores_PC1_treatment_proximal, 2.5); prctile(cvscores_PC1_treatment_sleeve, 2.5)];
upper_PC1_treatment = [prctile(cvscores_PC1_treatment_distal, 97.5); prctile(cvscores_PC1_treatment_proximal, 97.5); prctile(cvscores_PC1_treatment_sleeve, 97.5)];

lower_PC2_treatment = [prctile(cvscores_PC2_treatment_distal, 2.5); prctile(cvscores_PC2_treatment_proximal, 2.5); prctile(cvscores_PC2_treatment_sleeve, 2.5)];
upper_PC2_treatment = [prctile(cvscores_PC2_treatment_distal, 97.5); prctile(cvscores_PC2_treatment_proximal, 97.5); prctile(cvscores_PC2_treatment_sleeve, 97.5)];

% Find 2.5th and 97.5th percentiles of the loadings for the time*treatment
% effect
percentiles_cvloadings_PC1_treatment = prctile(cvloadings_PC1_treatment, [2.5, 97.5]);
percentiles_cvloadings_PC2_treatment = prctile(cvloadings_PC2_treatment, [2.5, 97.5]);

% Find 2.5th and 97.5th percentiles of scores for time*treatment effect for
% each of the groups.
lower_PC1_time_treatment = [prctile(cvscores_PC1_time_treatment_distal, 2.5); prctile(cvscores_PC1_time_treatment_proximal, 2.5); prctile(cvscores_PC1_time_treatment_sleeve, 2.5)];
upper_PC1_time_treatment = [prctile(cvscores_PC1_time_treatment_distal, 97.5); prctile(cvscores_PC1_time_treatment_proximal, 97.5); prctile(cvscores_PC1_time_treatment_sleeve, 97.5)];

lower_PC2_time_treatment = [prctile(cvscores_PC2_time_treatment_distal, 2.5); prctile(cvscores_PC2_time_treatment_proximal, 2.5); prctile(cvscores_PC2_time_treatment_sleeve, 2.5)];
upper_PC2_time_treatment = [prctile(cvscores_PC2_time_treatment_distal, 97.5); prctile(cvscores_PC2_time_treatment_proximal, 97.5); prctile(cvscores_PC2_time_treatment_sleeve, 97.5)];

% Find 2.5th and 97.5th percentiles of the loadings for the time*treatment
% effect
percentiles_cvloadings_PC1_time_treatment = prctile(cvloadings_PC1_time_treatment, [2.5, 97.5]);
percentiles_cvloadings_PC2_time_treatment = prctile(cvloadings_PC2_time_treatment, [2.5, 97.5]);

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
errorbar([0, 3, 6, 9, 12], PC1_scores_time, PC1_scores_time-percentiles_cvscores_PC1_time(1,:), percentiles_cvscores_PC1_time(2,:)-PC1_scores_time)
set(gca,'XTick',[0 3 6 9 12]);
ylabel("PC1 (" + num2str((eigen_time(1)/sum(eigen_time))*100, '%.2f') + "%)", "FontSize",8)
title("Scores", 'FontSize',8)
grid on
hold off

% Loading plot
nexttile
bar(1:size(Y, 2), loadings_time(:,1)')
hold on
errorbar([1:size(Y, 2)], loadings_time(:,1)', loadings_time(:,1)' - percentiles_cvloadings_PC1_time(1,:), percentiles_cvloadings_PC1_time(2,:) - loadings_time(:,1)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(Y, 2)]);
set(gca,'XTickLabel',abbreviations, 'XTickLabelRotation', 90, 'FontSize', 4);
title("Loadings", 'FontSize',8)
grid on
hold off

%  Plotting PC2-scores
nexttile
hold on 
errorbar([0, 3, 6, 9, 12], PC2_scores_time, PC2_scores_time-percentiles_cvscores_PC2_time(1,:), percentiles_cvscores_PC2_time(2,:)-PC2_scores_time)
set(gca,'XTick',[0 3 6 9 12]);
xlabel("Months", 'FontSize', 8)
ylabel("PC2 (" + num2str((eigen_time(2)/sum(eigen_time))*100, '%.2f') + "%)","FontSize",8)
grid on
hold off

% Loading plot
nexttile
bar(1:size(Y, 2), loadings_time(:,2)')
hold on
errorbar([1:size(Y, 2)], loadings_time(:,2)', loadings_time(:,2)'-percentiles_cvloadings_PC2_time(1,:), percentiles_cvloadings_PC2_time(2,:)-loadings_time(:,2)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(Y, 2)]);
set(gca,'XTickLabel',abbreviations, 'XTickLabelRotation', 90, 'FontSize', 4);
grid on
hold off

%% Time*Treatment interaction
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
    errorbar([0, 3, 6, 9, 12], PC1_scores_treatments(i,:), PC1_scores_treatments(i,:) - lower_PC1_treatment(i,:), upper_PC1_treatment(i,:) - PC1_scores_treatments(i,:))
end
set(gca,'XTick',[0 3 6 9 12]);
legend({'Distal', 'Proximal', 'Sleeve'}, 'FontSize', 4, 'Location', 'NorthWest')
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
    errorbar([0, 3, 6, 9, 12], PC2_scores_treatments(i,:), PC2_scores_treatments(i,:) - lower_PC2_treatment(i,:), upper_PC2_treatment(i,:) - PC2_scores_treatments(i,:))
end
set(gca,'XTick',[0 3 6 9 12]);
grid on
xlabel("Months", 'FontSize',8)
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
    errorbar([0, 3, 6, 9, 12], PC1_scores_time_treatments(i,:), PC1_scores_time_treatments(i,:) - lower_PC1_time_treatment(i,:), upper_PC1_time_treatment(i,:) - PC1_scores_time_treatments(i,:))
end
set(gca,'XTick',[0 3 6 9 12]);
legend({'Distal', 'Proximal', 'Sleeve'}, 'FontSize', 4, 'Location', 'South')
grid on
hold off

nexttile
bar(1:size(Y,2), loadings_time_treatment(:,1)')
hold on
errorbar([1:size(Y,2)], loadings_time_treatment(:,1)', loadings_time_treatment(:,1)' - percentiles_cvloadings_PC1_time_treatment(1,:), percentiles_cvloadings_PC1_time_treatment(2,:) - loadings_time_treatment(:,1)', 'LineStyle', 'none')
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
    errorbar([0, 3, 6, 9, 12], PC2_scores_time_treatments(i,:), PC2_scores_time_treatments(i,:) - lower_PC2_time_treatment(i,:), upper_PC2_time_treatment(i,:) - PC2_scores_time_treatments(i,:))
end
set(gca,'XTick',[0 3 6 9 12]);
grid on
xlabel("Months", 'FontSize',8)
hold off

nexttile
bar(1:size(Y,2), loadings_time_treatment(:,2)')
hold on
errorbar([1:size(Y,2)], loadings_time_treatment(:,2)', loadings_time_treatment(:,2)' - percentiles_cvloadings_PC2_time_treatment(1,:), percentiles_cvloadings_PC2_time_treatment(2,:) - loadings_time_treatment(:,2)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(Y,2)]);
set(gca,'XTickLabel', abbreviations, 'XTickLabelRotation', 90, 'FontSize', 4);
grid on
hold off