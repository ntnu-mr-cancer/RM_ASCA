function varargout = RM_ASCA(data, options)
% Function for doing repeated measures-ASCA.
% Requires the Statistics and Machine Learning toolbox.
% 
% [M_A, M_B, M_C, ZU, E, fig1, fig2, fig3] = RM_ASCA(data, options)
% The output elements M_A, M_B, and M_C describe the time effect, treatment + time*treatment effect, 
% and time + treatment + time*treatment effect, respectively. Each is a structure with the following
% fields:
% 
% M: The effect matrix for the factor(s)
% loadings: The loadings from PCA on the the effect matrix.
% scores: The scores from PCA on the effect matrices.
% eigen: Eigenvalues for the principal components for the effect matrices.
% scores_boot: A cell array containing the bootstrapped scores.
% loadings_boot: A cell array containing the bootstrapped loadings.
% 
% The outputs ZU and E are the random intercept and residual matrix, respectively. 
% [fig1, fig2, fig3] are figure handles for the figures for M_A, M_B, and M_C, respectively. 
%
% The function requires a data table (data), which in addition to the response variables also includes the following variables: 
% data.ID: Subject ID vector
% data.timepoint: Vector indicating timepoint
% data.treatment: Vector indicating treatment assignment
% 
% options is a structure with the following required fields:
% options.baseline: Set to either "ucLDA" or "cLDA".
% options.Y_vars: String array with response variable names
% 
% Optional fields:
% options.coding: Determines coding for the treatment variable. Default is sum coded treatment variable
% Set to "PRC" if treatment should be reference coded. If the treatment variable is a string, the group appearing first in the
% alphabet is the reference group. If it is a number, then lowest number is
% the reference group. If it is logical, then "false" is the reference category.
% options.iterations: Number of iterations to use when bootstrapping (1000 by default).
% options.CI: Confidence intervals are always calculated by default. To
% avoid this, set options.CI = [];
% options.plot: By default, plots are always produced. To avoid this, set to [];
% options.show_varnames: Set to "yes" if variable names should be displayed
% on x-axis for loadings.
% options.color: P x 3 matrix where P is the number of response variables,
% and every row is a color vector (e.g. [0 1 0]) indicating the color of
% the variable in the loading plot

%% Generate design- and response matrix
% Make handle for a function to create the design matrix
designmat = @(T, G) [ones(size(T,1),1), T, G, T(:,repmat([1:size(T,2)], 1, size(G,2))).*G(:,kron(1:size(G,2),ones(1,size(T,2))))];

% Set baseline as reference timepoint
dummy_time = dummyvar(data.timepoint);
dummy_time(:,1) =  [];

% Set sum coding as default for treatment variable
if isfield(options, 'coding') == 0
    options.coding = "ASCA"; 
end

if options.coding == "PRC"
    dummy_treatment = dummyvar(nominal(data.treatment));
    dummy_treatment(:,1) = [];
elseif options.coding == "ASCA"
    treatments = sort(unique(data.treatment));
    for i = 1:length(treatments)-1
        dummy_treatment(:,i) = (data.treatment == treatments(i)) - (data.treatment == treatments(length(treatments)));
    end
end

% Create fixed effect design matrix and specify column numbers for the effects
if options.baseline == "cLDA"
    X = designmat(dummy_time, dummy_treatment);
    X(:,(1+length(unique(data.timepoint))):(length(unique(data.timepoint))+size(dummy_treatment,2))) = [];
    effects = {[2:length(unique(data.timepoint))], ...
    [(1+length(unique(data.timepoint))):size(X,2)], ...
    [2:size(X,2)]};
elseif options.baseline == "ucLDA"
    X = designmat(dummy_time, dummy_treatment);
    effects = {[2:length(unique(data.timepoint))], ...
    [(1+length(unique(data.timepoint))):size(X,2)], ...
    [2:size(X,2)]};
end

% Define other matrices for mixed model analysis
Y_vars = options.Y_vars;
Y = data{:,Y_vars};
G = data.ID;
Z = ones(size(data,1),1);

% List all timepoints and treatments
timepoints = sort(unique(data.timepoint), 'ascend');
treatments = sort(unique(data.treatment), 'ascend');

% Scaling response variables to their baseline standard deviation
scalingfactor = nanstd(Y(data.timepoint==timepoints(1),:));
for d = 1:size(Y,2)
    Y(:,d) = Y(:,d)./scalingfactor(d);
end

%% Fit model

% Run LMMs and collect coefficients, random effects, and residuals
for i = 1:size(Y,2)
    lme = fitlmematrix(X, Y(:,i), Z, G);
    b(:,i) = lme.Coefficients.Estimate;
    E(:,i) = residuals(lme);
    U(:,i) = randomEffects(lme);
end

% Make effect matrices for each of the effect combinations specified in "effects"-input
for k = 1:length(effects)
    M_f = X(:,effects{k})*b(effects{k},:);
    [loadings_f, scores_f, eigen_f] = pca(M_f, 'NumComponents', rank(M_f));
    varargout{k} = struct('M', M_f, 'loadings', loadings_f, 'scores', scores_f, 'eigen', eigen_f);
end

% Get random effect matrix ZU and residual matrix E and add to output
varargout{length(varargout)+1} = full(designMatrix(lme, 'random'))*U;
varargout{length(varargout)+1} = E;

%% Nonparametric bootstrapping to calculate confidence intervals
if isfield(options, 'CI') == 0;
    options.CI = "yes"; % do bootstrapping-CI by default
end

if options.CI == "yes"
    allpatients = unique(G, 'stable'); % Create vector containing all unique patient IDs

    for i = 1:size(allpatients,1)
        groups_allpatients(i) = unique(data.treatment(find(G==allpatients(i)))); % Find group membership for each patient
    end

    for i = 1:options.iterations;
        
        % Make bootstrap sample
        X_b = [];
        Y_b = [];
        G_b = 0;  
        data2 = table();
        for t = 1:length(treatments)
            for k = 1:length(allpatients(groups_allpatients==treatments(t)))
                idx = find(groups_allpatients == treatments(t));
                idx = idx(randi(length(idx)));

                X_add = X(data.ID == allpatients(idx,:),:);
                X_b = [X_b; X_add];

                Y_add = data{data.ID == allpatients(idx),Y_vars};
                Y_b = [Y_b; Y_add];

                add = data(data.ID == allpatients(idx,:),:);
                data2 = [data2; add];
            end
        end
        G_b = data2.ID;

        clear Y_add X_add G_add
        Z_b = ones(size(X_b,1),1);
        
        % Scale to baseline SD
        scalingfactor = std(Y_b(data2.timepoint==timepoints(1),:)); % Scale to baseline standard deviation
        for d = 1:size(Y_b,2)
            Y_b(:,d) = Y_b(:,d)./scalingfactor(d);
        end
        
       % Run LMMs and collect coefficients and residuals
        for d = 1:size(Y_b,2)
            lme = fitlmematrix(X_b, Y_b(:,d), Z_b, G_b);
            b(:,d) = lme.Coefficients.Estimate;
        end

        % Make effect matrices for each of the effect combinations specified in "effects"-input
        M_Aj = X_b(:,effects{1})*b(effects{1},:);
        [loadings, scores, eigen] = pca(M_Aj, 'NumComponents', rank(M_Aj));
        M_Aj = struct('M', M_Aj, 'loadings', loadings, 'scores', scores, 'eigen', eigen);
        
        M_Bj = X_b(:,effects{2})*b(effects{2},:);
        [loadings, scores, eigen] = pca(M_Bj, 'NumComponents', rank(M_Bj));
        M_Bj = struct('M', M_Bj, 'loadings', loadings, 'scores', scores, 'eigen', eigen);
        
        M_Cj = X_b(:,effects{3})*b(effects{3},:);
        [loadings, scores, eigen] = pca(M_Cj, 'NumComponents', rank(M_Cj));
        M_Cj = struct('M', M_Cj, 'loadings', loadings, 'scores', scores, 'eigen', eigen);

        % Rotate scores and loadings for A
        [loadings_Aj_rot, rotationmatrix_Aj] = rotatefactors(M_Aj.loadings, 'Method', 'procrustes', 'type', 'orthogonal', 'Target', varargout{1}.loadings); 
        scores_Aj_rot = M_Aj.scores*inv(rotationmatrix_Aj);

        % Rotate scores and loadings for B
        [loadings_Bj_rot, rotationmatrix_Bj] = rotatefactors(M_Bj.loadings, 'Method', 'procrustes', 'type', 'orthogonal', 'Target', varargout{2}.loadings); 
        scores_Bj_rot = M_Bj.scores*inv(rotationmatrix_Bj);

         % Rotate scores and loadings for C
        [loadings_Cj_rot, rotationmatrix_Cj] = rotatefactors(M_Cj.loadings, 'Method', 'procrustes', 'type', 'orthogonal', 'Target', varargout{3}.loadings); 
        scores_Cj_rot = M_Cj.scores*inv(rotationmatrix_Cj); 

        % Finding scores to plot
        for r = 1:size(varargout{1}.scores, 2)
            for f = 1:length(unique(timepoints))
                scores_Aj{r,1}(i,f) = scores_Aj_rot(find(data2.timepoint == timepoints(f),1, 'first'), r);
            end
            loadings_Aj{r,1}(i,:) = loadings_Aj_rot(:,r)';
        end

        for r = 1:size(varargout{2}.scores, 2)
            for t = 1:length(treatments)
                for f = 1:length(timepoints)
                    scores_Bj{r,t}(i,f) = scores_Bj_rot(find(data2.timepoint == timepoints(f) & data2.treatment == treatments(t),1, 'first'), r);
                end
            end
            loadings_Bj{r,1}(i,:) = loadings_Bj_rot(:,r)';
        end

        for r = 1:size(varargout{3}.scores, 2)
            for t = 1:length(treatments)
                for f = 1:length(timepoints)
                    scores_Cj{r,t}(i,f) = scores_Cj_rot(find(data2.timepoint == timepoints(f) & data2.treatment == treatments(t),1, 'first'), r);
                end
            end
            loadings_Cj{r,1}(i,:) = loadings_Cj_rot(:,r)';
        end
    end
varargout{1}.scores_boot = scores_Aj;
varargout{1}.loadings_boot = loadings_Aj;

varargout{2}.scores_boot = scores_Bj;
varargout{2}.loadings_boot = loadings_Bj;

varargout{3}.scores_boot = scores_Cj;
varargout{3}.loadings_boot = loadings_Cj;
end

%% Create figures

if isfield(options, 'plot') == 0;
    options.plot = "yes"; % make plots by default
end

if options.plot == "yes"
    varargout{length(varargout)+1} = figure;
    % Scree plot
    tiledlayout(2,5); nexttile([2 1])
    bar(varargout{1}.eigen(1:rank(varargout{1}.scores,2))./sum(varargout{1}.eigen)*100)
    ylim([0,100])
    xlabel("Principal components", 'FontSize',7)
    set(gca, 'xtick', []);
    ylabel("Explained variance (%)", 'FontSize',8)
    title("Scree plot", 'FontSize',8)
    grid on

    %  Plotting PC1-scores
    for i = 1:length(unique(timepoints))
        id(i) = find(data.timepoint==timepoints(i), 1, 'first');
    end

    nexttile([1 2])
    hold on 
    errorbar([0:length(timepoints)-1], varargout{1}.scores(id,1)', varargout{1}.scores(id,1)'-prctile(scores_Aj{1}, 2.5), prctile(scores_Aj{1}, 97.5) - varargout{1}.scores(id,1)')
    set(gca,'XTick',[0:length(timepoints)-1]);
    ylabel("PC1 (" + num2str((varargout{1}.eigen(1)/sum(varargout{1}.eigen))*100, '%.2f') + "%)", "FontSize",8)
    title("Scores", 'FontSize',8)
    grid on
    hold off

    % Loading plot
    nexttile([1 2])
    h = bar(1:size(Y,2), varargout{1}.loadings(:,1)');
    h.FaceColor = 'flat';
    if isfield(options, 'color')
        h.CData = options.color;
    end
    hold on
    errorbar([1:size(Y,2)], varargout{1}.loadings(:,1)', varargout{1}.loadings(:,1)' - prctile(loadings_Aj{1}, 2.5), prctile(loadings_Aj{1}, 97.5) - varargout{1}.loadings(:,1)', 'LineStyle', 'none')
    if isfield(options, 'show_varnames') == 0
        options.show_varnames = "";
    end
    if options.show_varnames == "yes"
        set(gca,'XTick',[1:size(Y,2)]);
        set(gca,'XTickLabel',Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
    end
    title("Loadings", 'FontSize',8)
    grid on
    hold off

    %  Plotting PC2-scores
    nexttile([1 2])
    hold on 
    errorbar([0:length(timepoints)-1], varargout{1}.scores(id,2)', varargout{1}.scores(id,2)'-prctile(scores_Aj{2}, 2.5), prctile(scores_Aj{2}, 97.5) - varargout{1}.scores(id,2)')
    set(gca,'XTick',[0:length(timepoints)-1]);
    ylabel("PC1 (" + num2str((varargout{1}.eigen(2)/sum(varargout{1}.eigen))*100, '%.2f') + "%)", "FontSize",8)
    grid on
    xlabel("Timepoint")
    hold off

    % Loading plot
    nexttile([1 2])
    h = bar(1:size(Y,2), varargout{1}.loadings(:,2)');
    h.FaceColor = 'flat';
    if isfield(options, 'color')
        h.CData = options.color;
    end
    hold on
    errorbar([1:size(Y,2)], varargout{1}.loadings(:,2)', varargout{1}.loadings(:,2)' - prctile(loadings_Aj{2}, 2.5), prctile(loadings_Aj{2}, 97.5) - varargout{1}.loadings(:,2)', 'LineStyle', 'none')
    if options.show_varnames == "yes"
        set(gca,'XTick',[1:size(Y,2)]);
        set(gca,'XTickLabel',Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
    end
    grid on
    hold off

    % (Treatment +) Time*Treatment
    varargout{length(varargout)+1} = figure;
    tiledlayout(2,5)
    nexttile([2 1])
    bar(varargout{2}.eigen(1:rank(varargout{2}.scores))./sum(varargout{2}.eigen)*100)
    set(gca, 'xtick', []);
    ylim([0,100])
    xlabel("Principal components", 'FontSize',7)
    ylabel("Explained variance (%)", 'FontSize',8)
    title("Scree plot", 'FontSize',8)
    grid on

    nexttile([1 2])
    hold on
    title("Scores", "FontSize", 8)
    ylabel("PC1 (" + num2str((varargout{2}.eigen(1)/sum(varargout{2}.eigen))*100, '%.2f') + " %)", 'FontSize', 8);
    for d = 1:length(treatments)
        for k = 1:length(timepoints)
            id(k) = find(data.timepoint == timepoints(k) & data.treatment == treatments(d), 1, 'first');
        end
        errorbar([0:length(unique(timepoints))-1], varargout{2}.scores(id,1)', varargout{2}.scores(id,1)' - prctile(scores_Bj{1,d}, 2.5), prctile(scores_Bj{1,d}, 97.5) - varargout{2}.scores(id,1)')
    end
    set(gca,'XTick',[0:length(timepoints)-1]);
    legend(string(treatments), 'FontSize', 4) %, 'FontSize', 4, 'Position', [0.5,2,1,1])
    grid on
    hold off

    nexttile([1 2])
    h = bar(1:size(Y,2), varargout{2}.loadings(:,1)');
    h.FaceColor = 'flat';
    if isfield(options, 'color')
        h.CData = options.color;
    end
    hold on
    errorbar([1:size(Y,2)], varargout{2}.loadings(:,1)', varargout{2}.loadings(:,1)' - prctile(loadings_Bj{1}, 2.5), prctile(loadings_Bj{1}, 97.5) - varargout{2}.loadings(:,1)', 'LineStyle', 'none')
    if options.show_varnames == "yes"
        set(gca,'XTick',[1:size(Y,2)]);
        set(gca,'XTickLabel',Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
    end
    title("Loadings", 'FontSize',8)
    grid on
    hold off

    % Plotting results for PC2
    nexttile([1 2])
    hold on
    ylabel("PC2 (" + num2str((varargout{2}.eigen(2)/sum(varargout{2}.eigen))*100, '%.2f') + " %)", 'FontSize', 8);
    for d = 1:length(treatments)
        for k = 1:length(timepoints)
            id(k) = find(data.timepoint == timepoints(k) & data.treatment == treatments(d), 1, 'first');
        end
        errorbar([0:length(unique(timepoints))-1], varargout{2}.scores(id,2)', varargout{2}.scores(id,2)' - prctile(scores_Bj{2,d}, 2.5), prctile(scores_Bj{2,d}, 97.5) - varargout{2}.scores(id,2)')
    end
    set(gca,'XTick',[0:length(timepoints)-1]);
    xlabel("Timepoint")
    grid on
    hold off

    nexttile([1 2])
    h = bar(1:size(Y,2), varargout{2}.loadings(:,2)');
    h.FaceColor = 'flat';
    if isfield(options, 'color')
        h.CData = options.color;
    end
    hold on
    errorbar([1:size(Y,2)], varargout{2}.loadings(:,2)', varargout{2}.loadings(:,2)' - prctile(loadings_Bj{2}, 2.5), prctile(loadings_Bj{2}, 97.5) - varargout{2}.loadings(:,2)', 'LineStyle', 'none')
    if options.show_varnames == "yes"
        set(gca,'XTick',[1:size(Y,2)]);
        set(gca,'XTickLabel',Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
    end
    title("Loadings", 'FontSize',8)
    grid on
    hold off

    % Time + (Treatment +) Time*Treatment
    varargout{length(varargout)+1} = figure;
    tiledlayout(2,5)
    nexttile([2 1])
    bar(varargout{3}.eigen(1:rank(varargout{3}.scores))./sum(varargout{3}.eigen)*100)
    ylim([0,100])
    set(gca, 'xtick', []);
    xlabel("Principal components", 'FontSize',7)
    ylabel("Explained variance (%)", 'FontSize',8)
    title("Scree plot", 'FontSize',8)
    grid on

    nexttile([1 2])
    hold on
    title("Scores", "FontSize", 8)
    ylabel("PC1 (" + num2str((varargout{3}.eigen(1)/sum(varargout{3}.eigen))*100, '%.2f') + " %)", 'FontSize', 8);
    for d = 1:length(treatments)
        for k = 1:length(timepoints)
            id(k) = find(data.timepoint == timepoints(k) & data.treatment == treatments(d), 1, 'first');
        end
        errorbar([0:length(unique(timepoints))-1], varargout{3}.scores(id,1)', varargout{3}.scores(id,1)' - prctile(scores_Cj{1,d}, 2.5), prctile(scores_Cj{1,d}, 97.5) - varargout{3}.scores(id,1)')
    end
    set(gca,'XTick',[0:length(timepoints)-1]);
    legend(string(treatments), 'FontSize', 4) %, 'FontSize', 4, 'Position', [0.5,2,1,1])
    grid on
    hold off

    nexttile([1 2])
    h = bar(1:size(Y,2), varargout{3}.loadings(:,1)');
    h.FaceColor = 'flat';
    if isfield(options, 'color')
        h.CData = options.color;
    end
    hold on
    errorbar([1:size(Y,2)], varargout{3}.loadings(:,1)', varargout{3}.loadings(:,1)' - prctile(loadings_Cj{1}, 2.5), prctile(loadings_Cj{1}, 97.5) - varargout{3}.loadings(:,1)', 'LineStyle', 'none')
    if options.show_varnames == "yes"
        set(gca,'XTick',[1:size(Y,2)]);
        set(gca,'XTickLabel',Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
    end
    title("Loadings", 'FontSize',8)
    grid on
    hold off

    % Plotting results for PC2
    nexttile([1 2])
    hold on
    ylabel("PC2 (" + num2str((varargout{3}.eigen(2)/sum(varargout{3}.eigen))*100, '%.2f') + " %)", 'FontSize', 8);
    for d = 1:length(treatments)
        for k = 1:length(timepoints)
            id(k) = find(data.timepoint == timepoints(k) & data.treatment == treatments(d), 1, 'first');
        end
        errorbar([0:length(unique(timepoints))-1], varargout{3}.scores(id,2)', varargout{3}.scores(id,2)' - prctile(scores_Cj{2,d}, 2.5), prctile(scores_Cj{2,d}, 97.5) - varargout{3}.scores(id,2)')
    end
    set(gca,'XTick',[0:length(timepoints)-1]);
    xlabel("Timepoint")
    grid on
    hold off

    nexttile([1 2])
    h = bar(1:size(Y,2), varargout{3}.loadings(:,2)');
    h.FaceColor = 'flat';
    if isfield(options, 'color')
        h.CData = options.color;
    end
    hold on
    errorbar([1:size(Y,2)], varargout{3}.loadings(:,2)', varargout{3}.loadings(:,2)' - prctile(loadings_Cj{2}, 2.5), prctile(loadings_Cj{2}, 97.5) - varargout{3}.loadings(:,2)', 'LineStyle', 'none')
    if options.show_varnames == "yes"
        set(gca,'XTick',[1:size(Y,2)]);
        set(gca,'XTickLabel',Y_vars, 'XTickLabelRotation', 90, 'FontSize', 4);
    end
    title("Loadings", 'FontSize',8)
    grid on
    hold off
end

%% For univariate p-values

if isfield(options, 'pval') == 0;
    options.pval = ""; % no table by default
end

if options.pval == "yes"
    for i = 1:size(Y,2)
        lme = fitlmematrix(X, Y(:,i), Z, G);
        B(:,i) = lme.Coefficients.Estimate;
        result = anova(lme, 'DFMethod', 'satterthwaite');
        p(:,i) = result.pValue;
    end
    
    for i = 1:size(Y,2)
        [~, ~, q(:,i)] = fdr_bh(p(:,i));
    end
    
    pvaltab = table();
    for i = 1:size(X,2)
        tab_P = table(B(i,:)', p(i,:)', q(i,:)');
        tab_P.Properties.VariableNames = ["coef"+i, "p"+i, "q"+i];
        pvaltab = [pvaltab, tab_P];
    end
    pvaltab.VarNames = Y_vars';
    varargout{length(varargout)+1} = pvaltab(:, ["VarNames", pvaltab.Properties.VariableNames(1:end-1)]);
end

end