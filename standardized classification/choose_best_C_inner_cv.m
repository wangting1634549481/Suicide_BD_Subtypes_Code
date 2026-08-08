function best_C = choose_best_C_inner_cv(X, y, covariates, C_grid, innerK)
    % Select best SVM BoxConstraint C using inner cross-validation.

    cvp_inner = cvpartition(y, 'KFold', innerK);

    mean_bacc = zeros(length(C_grid), 1);

    for c = 1:length(C_grid)
        C = C_grid(c);
        bacc_folds = zeros(innerK, 1);

        for fold = 1:innerK
            train_idx = training(cvp_inner, fold);
            val_idx = test(cvp_inner, fold);

            X_train = X(train_idx, :);
            y_train = y(train_idx);
            cov_train = covariates(train_idx, :);

            X_val = X(val_idx, :);
            y_val = y(val_idx);
            cov_val = covariates(val_idx, :);

            model = train_preprocessed_svm(X_train, y_train, cov_train, C);

            result = evaluate_preprocessed_svm(model, X_val, y_val, cov_val, 0, ...
                'Inner validation');

            bacc_folds(fold) = result.balanced_accuracy;
        end

        mean_bacc(c) = mean(bacc_folds);
    end

    % If multiple C values tie, choose the smallest C for stronger regularization.
    max_bacc = max(mean_bacc);
    best_idx = find(mean_bacc == max_bacc, 1, 'first');
    best_C = C_grid(best_idx);
end
