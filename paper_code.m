%% RM-ASCA+ - NeoAva 
addpath(genpath('/mnt/work/RM_ASCA+/Functions'))
addpath(genpath('/mnt/work/RM_ASCA+/Paper_code'))

NeoAva = readtable('NeoAva.txt'); % Import NeoAva-data
NeoAva.timepoint = nominal(NeoAva.timepoint); % change from cell format to nominal

% Specify modeling options
options.iterations = 2000; % Setting number of bootstrap iterations
options.baseline = "cLDA"; % Constrain the baseline means
options.Y_vars = NeoAva.Properties.VariableNames(5:end); % Input list of response variable names
options.show_varnames = "yes"; % Show variable names on x-axis
options.coding = "PRC"; % Use reference coding for the treatment factor
options.pval = "yes"; % Calculate univariate p-values
options.plot = ""; % Turn off automatic plotting

% Run RM-ASCA+
[M_T1, M_TG1, M_TTG1, ZU1, E1, unipval_tab1] = RM_ASCA(NeoAva, options);

% Time effect - NeoAva
% Scree plot
figure
tiledlayout(2,5); nexttile([2 1])
bar(M_T1.eigen(1:rank(M_T1.scores,2))./sum(M_T1.eigen)*100)
ylim([0,100])
xlabel("Principal components", 'FontSize',7)
set(gca, 'xtick', []);
ylabel("Explained variance (%)", 'FontSize',8)
title("Scree plot", 'FontSize',8)
grid on

timepoints = unique(NeoAva.timepoint);
treatments = unique(NeoAva.treatment);

%  Plotting PC1-scores
for i = 1:length(unique(timepoints))
    id(i) = find(NeoAva.timepoint==timepoints(i), 1, 'first');
end

nexttile([1 2])
hold on 
errorbar([0:length(timepoints)-1], M_T1.scores(id,1)', M_T1.scores(id,1)'-prctile(M_T1.scores_boot{1}, 2.5), prctile(M_T1.scores_boot{1}, 97.5) - M_T1.scores(id,1)')
set(gca,'XTick',[0:length(timepoints)-1]);
ylabel("PC1 (" + num2str((M_T1.eigen(1)/sum(M_T1.eigen))*100, '%.2f') + "%)", "FontSize",8)
title("Scores", 'FontSize',8)
grid on
hold off

% Loading plot
nexttile([1 2]); hold on
bar(1:size(M_T1.loadings,1), M_T1.loadings(:,1)');
errorbar([1:size(M_T1.loadings,1)], M_T1.loadings(:,1)', M_T1.loadings(:,1)' - prctile(M_T1.loadings_boot{1}, 2.5), prctile(M_T1.loadings_boot{1}, 97.5) - M_T1.loadings(:,1)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(M_T1.loadings,1)]);
set(gca,'XTickLabel',options.Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
title("Loadings", 'FontSize',8)
grid on
hold off

% Plotting PC2-scores
nexttile([1 2])
hold on 
errorbar([0:length(timepoints)-1], M_T1.scores(id,2)', M_T1.scores(id,2)'-prctile(M_T1.scores_boot{2}, 2.5), prctile(M_T1.scores_boot{2}, 97.5) - M_T1.scores(id,2)')
set(gca,'XTick',[0:length(timepoints)-1]);
ylabel("PC2 (" + num2str((M_T1.eigen(2)/sum(M_T1.eigen))*100, '%.2f') + "%)", "FontSize",8)
xlabel("Timepoint", "FontSize", 8)
grid on
hold off

% Loading plot
nexttile([1 2]); hold on
bar(1:size(M_T1.loadings,1), M_T1.loadings(:,2)');
errorbar([1:size(M_T1.loadings,1)], M_T1.loadings(:,2)', M_T1.loadings(:,2)' - prctile(M_T1.loadings_boot{2}, 2.5), prctile(M_T1.loadings_boot{2}, 97.5) - M_T1.loadings(:,2)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(M_T1.loadings,1)]);
set(gca,'XTickLabel',options.Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
grid on
hold off

% Time*Treatment effect - NeoAva
figure('Position', [100 100 800 400])
tiledlayout(1,5)
nexttile([1 1])
bar(M_TG1.eigen(1:rank(M_TG1.scores))./sum(M_TG1.eigen)*100)
set(gca, 'xtick', []);
ylim([0,100])
xlabel("Principal components", 'FontSize',7)
ylabel("Explained variance (%)", 'FontSize',8)
title("Scree plot", 'FontSize',8)
grid on

nexttile([1 2])
hold on
title("Scores", "FontSize", 8)
xlabel("Timepoint", "FontSize", 8)
ylabel("PC1 (" + num2str((M_TG1.eigen(1)/sum(M_TG1.eigen))*100, '%.2f') + " %)", 'FontSize', 8);
for d = 1:length(treatments)
    for k = 1:length(timepoints)
        id(k) = find(NeoAva.timepoint == timepoints(k) & NeoAva.treatment == treatments(d), 1, 'first');
    end
    errorbar([0:length(unique(timepoints))-1], M_TG1.scores(id,1)', M_TG1.scores(id,1)' - prctile(M_TG1.scores_boot{1,d}, 2.5), prctile(M_TG1.scores_boot{1,d}, 97.5) - M_TG1.scores(id,1)')
end
set(gca,'XTick',[0:length(timepoints)-1]);
legend({'CTX', 'CTX + B'}, 'FontSize', 4) %, 'FontSize', 4, 'Position', [0.5,2,1,1])
grid on
hold off

nexttile([1 2])
bar(1:size(M_T1.loadings,1), M_TG1.loadings(:,1)');
hold on
errorbar([1:size(M_T1.loadings,1)], M_TG1.loadings(:,1)', M_TG1.loadings(:,1)' - prctile(M_TG1.loadings_boot{1}, 2.5), prctile(M_TG1.loadings_boot{1}, 97.5) - M_TG1.loadings(:,1)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(M_T1.loadings,1)]);
set(gca,'XTickLabel',options.Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
title("Loadings", 'FontSize',8)
grid on
hold off

% Time + Time*Treatment effect - NeoAva
figure('Position', [100 100 800 400])
tiledlayout(1,5)
nexttile([1 1])
bar(M_TTG1.eigen(1:rank(M_TTG1.scores))./sum(M_TTG1.eigen)*100)
set(gca, 'xtick', []);
ylim([0,100])
xlabel("Principal components", 'FontSize',7)
ylabel("Explained variance (%)", 'FontSize',8)
title("Scree plot", 'FontSize',8)
grid on

nexttile([1 2])
hold on
title("Scores", "FontSize", 8)
xlabel("Timepoint", "FontSize", 8)
ylabel("PC1 (" + num2str((M_TTG1.eigen(1)/sum(M_TTG1.eigen))*100, '%.2f') + " %)", 'FontSize', 8);
for d = 1:length(treatments)
    for k = 1:length(timepoints)
        id(k) = find(NeoAva.timepoint == timepoints(k) & NeoAva.treatment == treatments(d), 1, 'first');
    end
    errorbar([0:length(unique(timepoints))-1], M_TTG1.scores(id,1)', M_TTG1.scores(id,1)' - prctile(M_TTG1.scores_boot{1,d}, 2.5), prctile(M_TTG1.scores_boot{1,d}, 97.5) - M_TTG1.scores(id,1)')
end
set(gca,'XTick',[0:length(timepoints)-1]);
legend({'CTX', 'CTX + B'}, 'FontSize', 4) %, 'FontSize', 4, 'Position', [0.5,2,1,1])
grid on
hold off

nexttile([1 2])
bar(1:size(M_T1.loadings,1), M_TTG1.loadings(:,1)');
hold on
errorbar([1:size(M_T1.loadings,1)], M_TTG1.loadings(:,1)', M_TTG1.loadings(:,1)' - prctile(M_TTG1.loadings_boot{1}, 2.5), prctile(M_TTG1.loadings_boot{1}, 97.5) - M_TTG1.loadings(:,1)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(M_T1.loadings,1)]);
set(gca,'XTickLabel',options.Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
title("Loadings", 'FontSize',8)
grid on
hold off

%% Augmenting effect matrix to visualize response variability

M_TTG1.scores_aug1 = ((M_TTG1.M + ZU1) - mean(M_TTG1.M + ZU1 + E1))*M_TTG1.loadings;
M_TTG1.scores_aug2 = ((M_TTG1.M + ZU1 + E1) - mean(M_TTG1.M + ZU1 + E1))*M_TTG1.loadings;

timepoints = sort(unique(NeoAva.timepoint));
allpatients = sort(unique(NeoAva.ID));

for i = 1:size(allpatients,1)
    groups_allpatients(i) = unique(NeoAva.treatment(find(NeoAva.ID==allpatients(i)))); % Find group membership for each patient
end

% Randomly generated vector to select which patients to plot
idx = [7     9    19    28    37    41    54    55    58    60    63    67    68    72    74    75    76    81    86    89    91    93   100];

for i = 1:length(idx)
    for t = 1:length(timepoints)
        PC1_aug1(i,t) = M_TTG1.scores_aug1(NeoAva.ID==idx(i) & NeoAva.timepoint == timepoints(t),1);
        PC1_aug2(i,t) = M_TTG1.scores_aug2(NeoAva.ID==idx(i) & NeoAva.timepoint == timepoints(t),1);
    end
end

colors = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410];

figure('Position', [500 500 500 400])
hold on
for i = 1:length(idx)
    plot([0 1 2], PC1_aug1(i,:), 'color', colors(double(groups_allpatients(idx(i))==0)+1,:))
    plot([0 1 2], PC1_aug2(i,:), '--', 'LineWidth', 0.2, 'color', colors(double(groups_allpatients(idx(i))==0)+1,:))
end
set(gca, 'xtick', [0 1 2])
xlabel("Timepoint")
title("PC1, augmented")


%% RM-ASCA+ - Bariatric surgery
import_bariatricsurgery

options.baseline = "ucLDA"; % Unconstrained baseline means
options.Y_vars = data.Properties.VariableNames(4:end); % List of response variable names
options.coding = "ASCA"; % Use sum coding for the treatment factor

% Run RM-ASCA+
[M_T3, M_TG3, M_TTG3, ZU3, E3, unipval_tab2] = RM_ASCA(data, options);

% Time effect - Bariatric surgery
% Scree plot
figure
tiledlayout(2,5); nexttile([2 1])
bar(M_T3.eigen(1:rank(M_T3.scores,2))./sum(M_T3.eigen)*100)
ylim([0,100])
xlabel("Principal components", 'FontSize',7)
set(gca, 'xtick', []);
ylabel("Explained variance (%)", 'FontSize',8)
title("Scree plot", 'FontSize',8)
grid on

timepoints = unique(data.timepoint);
treatments = unique(data.treatment);

%  Plotting PC1-scores
for i = 1:length(unique(timepoints))
    id(i) = find(data.timepoint==timepoints(i), 1, 'first');
end

nexttile([1 2])
hold on 
errorbar([0:length(timepoints)-1], M_T3.scores(id,1)', M_T3.scores(id,1)'-prctile(M_T3.scores_boot{1}, 2.5), prctile(M_T3.scores_boot{1}, 97.5) - M_T3.scores(id,1)')
set(gca,'XTick',[0:length(timepoints)-1]);
ylabel("PC1 (" + num2str((M_T3.eigen(1)/sum(M_T3.eigen))*100, '%.2f') + "%)", "FontSize",8)
title("Scores", 'FontSize',8)
grid on
hold off

% Loading plot
nexttile([1 2]); hold on
bar(1:size(M_T3.loadings,1), M_T3.loadings(:,1)');
errorbar([1:size(M_T3.loadings,1)], M_T3.loadings(:,1)', M_T3.loadings(:,1)' - prctile(M_T3.loadings_boot{1}, 2.5), prctile(M_T3.loadings_boot{1}, 97.5) - M_T3.loadings(:,1)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(M_T3.loadings,1)]);
set(gca,'XTickLabel',options.Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
title("Loadings", 'FontSize',8)
grid on
hold off

% Plotting PC2-scores
nexttile([1 2])
hold on 
errorbar([0:length(timepoints)-1], M_T3.scores(id,2)', M_T3.scores(id,2)'-prctile(M_T3.scores_boot{2}, 2.5), prctile(M_T3.scores_boot{2}, 97.5) - M_T3.scores(id,2)')
set(gca,'XTick',[0:length(timepoints)-1]);
ylabel("PC2 (" + num2str((M_T3.eigen(2)/sum(M_T3.eigen))*100, '%.2f') + "%)", "FontSize",8)
xlabel("Timepoint", "FontSize", 8)
grid on
hold off

% Loading plot
nexttile([1 2]); hold on
bar(1:size(M_T3.loadings,1), M_T3.loadings(:,2)');
errorbar([1:size(M_T3.loadings,1)], M_T3.loadings(:,2)', M_T3.loadings(:,2)' - prctile(M_T3.loadings_boot{2}, 2.5), prctile(M_T3.loadings_boot{2}, 97.5) - M_T3.loadings(:,2)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(M_T3.loadings,1)]);
set(gca,'XTickLabel',options.Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
grid on
hold off

% Treatment + Time*Treatment effect - Bariatric surgery
figure('Position', [100 100 800 400])
tiledlayout(1,5)
nexttile([1 1])
bar(M_TG3.eigen(1:rank(M_TG3.scores))./sum(M_TG3.eigen)*100)
set(gca, 'xtick', []);
ylim([0,100])
xlabel("Principal components", 'FontSize',7)
ylabel("Explained variance (%)", 'FontSize',8)
title("Scree plot", 'FontSize',8)
grid on

nexttile([1 2])
hold on
title("Scores", "FontSize", 8)
ylabel("PC1 (" + num2str((M_TG3.eigen(1)/sum(M_TG3.eigen))*100, '%.2f') + " %)", 'FontSize', 8);
xlabel("Timepoint", "FontSize", 8)
for d = 1:length(treatments)
    for k = 1:length(timepoints)
        id(k) = find(data.timepoint == timepoints(k) & data.treatment == treatments(d), 1, 'first');
    end
    errorbar([0:length(unique(timepoints))-1], M_TG3.scores(id,1)', M_TG3.scores(id,1)' - prctile(M_TG3.scores_boot{1,d}, 2.5), prctile(M_TG3.scores_boot{1,d}, 97.5) - M_TG3.scores(id,1)')
end
set(gca,'XTick',[0:length(timepoints)-1]);
legend({'Distal', 'Proximal', 'Sleeve'}, 'FontSize', 4, 'Location', 'NorthWest') %, 'FontSize', 4, 'Position', [0.5,2,1,1])
grid on
hold off

nexttile([1 2])
bar(1:size(M_T3.loadings,1), M_TG3.loadings(:,1)');
hold on
errorbar([1:size(M_T3.loadings,1)], M_TG3.loadings(:,1)', M_TG3.loadings(:,1)' - prctile(M_TG3.loadings_boot{1}, 2.5), prctile(M_TG3.loadings_boot{1}, 97.5) - M_TG3.loadings(:,1)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(M_T3.loadings,1)]);
set(gca,'XTickLabel',options.Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
title("Loadings", 'FontSize',8)
grid on
hold off

% Time + Treatment + Time*Treatment effect - Bariatric surgery
figure
tiledlayout(2,5)
nexttile([2 1])
bar(M_TTG3.eigen(1:rank(M_TTG3.scores))./sum(M_TTG3.eigen)*100)
set(gca, 'xtick', []);
ylim([0,100])
xlabel("Principal components", 'FontSize',7)
ylabel("Explained variance (%)", 'FontSize',8)
title("Scree plot", 'FontSize',8)
grid on

nexttile([1 2])
hold on
title("Scores", "FontSize", 8)
ylabel("PC1 (" + num2str((M_TTG3.eigen(1)/sum(M_TTG3.eigen))*100, '%.2f') + " %)", 'FontSize', 8);
for d = 1:length(treatments)
    for k = 1:length(timepoints)
        id(k) = find(data.timepoint == timepoints(k) & data.treatment == treatments(d), 1, 'first');
    end
    errorbar([0:length(unique(timepoints))-1], M_TTG3.scores(id,1)', M_TTG3.scores(id,1)' - prctile(M_TTG3.scores_boot{1,d}, 2.5), prctile(M_TTG3.scores_boot{1,d}, 97.5) - M_TTG3.scores(id,1)')
end
set(gca,'XTick',[0:length(timepoints)-1]);
legend({'Distal', 'Proximal', 'Sleeve'}, 'FontSize', 4, 'Location', 'SouthEast') %, 'FontSize', 4, 'Position', [0.5,2,1,1])
grid on
hold off

nexttile([1 2])
bar(1:size(M_T3.loadings,1), M_TTG3.loadings(:,1)');
hold on
errorbar([1:size(M_T3.loadings,1)], M_TTG3.loadings(:,1)', M_TTG3.loadings(:,1)' - prctile(M_TTG3.loadings_boot{1}, 2.5), prctile(M_TTG3.loadings_boot{1}, 97.5) - M_TTG3.loadings(:,1)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(M_T3.loadings,1)]);
set(gca,'XTickLabel',options.Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
title("Loadings", 'FontSize',8)
grid on
hold off

nexttile([1 2])
hold on
xlabel("Timepoint", "FontSize", 8)
ylabel("PC2 (" + num2str((M_TTG3.eigen(2)/sum(M_TTG3.eigen))*100, '%.2f') + " %)", 'FontSize', 8);
for d = 1:length(treatments)
    for k = 1:length(timepoints)
        id(k) = find(data.timepoint == timepoints(k) & data.treatment == treatments(d), 1, 'first');
    end
    errorbar([0:length(unique(timepoints))-1], M_TTG3.scores(id,2)', M_TTG3.scores(id,2)' - prctile(M_TTG3.scores_boot{2,d}, 2.5), prctile(M_TTG3.scores_boot{2,d}, 97.5) - M_TTG3.scores(id,2)')
end
set(gca,'XTick',[0:length(timepoints)-1]);
grid on
hold off

nexttile([1 2])
bar(1:size(M_T3.loadings,1), M_TTG3.loadings(:,2)');
hold on
errorbar([1:size(M_T3.loadings,1)], M_TTG3.loadings(:,2)', M_TTG3.loadings(:,2)' - prctile(M_TTG3.loadings_boot{2}, 2.5), prctile(M_TTG3.loadings_boot{2}, 97.5) - M_TTG3.loadings(:,2)', 'LineStyle', 'none')
set(gca,'XTick',[1:size(M_T3.loadings,1)]);
set(gca,'XTickLabel',options.Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
grid on
hold off

