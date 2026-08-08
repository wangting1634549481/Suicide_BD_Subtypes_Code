function model = train_preprocessed_svm(X_train, y_train, cov_train, C)
    % Training preprocessing:
    % 1. Residualize FC features for age, sex, education within training data.
    % 2. Z-score features using training mean and SD.
    % 3. Train linear SVM.

    y_train = y_train(:);

    % Covariate residualization
    D_train = [ones(size(cov_train, 1), 1), cov_train];
    beta_cov = D_train \ X_train;

    X_res_train = X_train - D_train * beta_cov;

    % Z-score normalization
    mu_train = mean(X_res_train, 1);
    sigma_train = std(X_res_train, 0, 1);
    sigma_train(sigma_train == 0) = 1;

    X_z_train = (X_res_train - mu_train) ./ sigma_train;

    % Train linear SVM
    svm_model = fitcsvm( ...
        X_z_train, y_train, ...
        'KernelFunction', 'linear', ...
        'BoxConstraint', C, ...
        'ClassNames', [0; 1], ...
        'Prior', 'uniform', ...
        'Standardize', false);

    model = struct();
    model.beta_cov = beta_cov;
    model.mu_train = mu_train;
    model.sigma_train = sigma_train;
    model.svm_model = svm_model;
    model.C = C;
end
