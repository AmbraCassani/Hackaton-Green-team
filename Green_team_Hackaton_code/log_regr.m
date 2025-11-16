function [] = log_regr(path)


%extract the excel file containing demographic data
X =[]; %pre allocate X
opt=detectImportOptions(path); 
shts=sheetnames(path); % usefull to not analyze the last excel sheet containing the legend
shts(end) = []; % deleting the "legend" sheet
for i=1:numel(shts)
    tab= readtable(path,opt,'Sheet',shts(i)); 
    l(i) = height(tab); % extract how many control and fallers are in the table
    X = [X; table2cell(tab)]; %transform the table into cell for easier analysi
end
X_temp = X;
X_temp(:,1:2) = [];
X_temp = cell2mat(X_temp); % Delete the columns containing the names and data of the recording and transform it into a matrix

[rows, ~] = find(isnan(X_temp)); % find any rows containing at least a NaN
row_to_remove = unique(rows); 

n_control = l(1) - length(find(row_to_remove<=l(1))); 
n_fallers = l(2) - length(find(row_to_remove>l(1))); %usefull to then form the label vector

X_withoutNaN = X_temp;
X_withoutNaN(row_to_remove,:) = [];


X_withoutNaN = zscore(X_withoutNaN);

y = [zeros(n_control,1); ones(n_fallers,1)]; %labels vector

[~, score, ~, ~, explained] = pca(X_withoutNaN);
how_much_variance = cumsum(explained);
idx_var = find(how_much_variance<=90);
idx_var = idx_var(end);
X = score(:,1:idx_var); % dimensionality reduction considering only principal components containing the 90% variance of the data

%%%%% Machine learning algorithm using stratified cross-validation and
%%%%% lasso regularization
numFolds = 5;
numRepeats = 5;

all_y_prob = [];
all_y_true = [];

rng('default'); % For reproducibility

for r = 1:numRepeats
    cvp = cvpartition(y, 'KFold', numFolds, 'Stratify', true);

    for i = 1:numFolds
        train_idx = training(cvp, i);
        test_idx = test(cvp, i);

        % Lasso Logistic Regression with built-in CV on training fold to select lambda
        [B, FitInfo] = lassoglm(X(train_idx,:), y(train_idx), 'binomial', 'CV', 5);

        % Get coefficients for best lambda
        idxLambda = FitInfo.IndexMinDeviance;
        coef = [FitInfo.Intercept(idxLambda); B(:, idxLambda)];

        % Add intercept term to test data
        X_test = [ones(sum(test_idx), 1), X(test_idx, :)];

        % Compute probabilites
        log_odds = X_test * coef;
        prob = 1 ./ (1 + exp(-log_odds));

        % Collect predictions
        all_y_prob = [all_y_prob; prob];
        all_y_true = [all_y_true; y(test_idx)];
    end
end

% Computation of ROC curve and Area Under the Curve (AUC)
[Xroc, Yroc, ~, AUC] = perfcurve(all_y_true, all_y_prob, 1);

% Plot ROC
figure;
plot(Xroc, Yroc, 'LineWidth', 2);
hold on;
grid on;
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title(sprintf('Repeated Stratified 5-Fold CV ROC Curve with L1 Regularization (AUC=%.2f)', AUC));




end
