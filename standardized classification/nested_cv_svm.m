function [cv_results, cv_summary] = nested_cv_svm(X, y, covariates, C_grid, outerK, innerK)
    % Nested cross-validation:
    % Outer loop evaluates model performance.
    % Inner loop selects SVM BoxConstraint C.
    % clear X y covariates;
    % X=X_disc;
    % y=y_disc;
    % covariates=cov_disc;


    cvp_outer = cvpartition(y, 'KFold', outerK);

    fold_id = zeros(outerK, 1);
    best_C_list = zeros(outerK, 1);
    acc_list = zeros(outerK, 1);
    bacc_list = zeros(outerK, 1);
    auc_list = zeros(outerK, 1);
    macro_f1_list = zeros(outerK, 1);
    sens0_list = zeros(outerK, 1);
    sens1_list = zeros(outerK, 1);

    for fold = 1:outerK
        train_idx = training(cvp_outer, fold);
        test_idx = test(cvp_outer, fold);

        X_train = X(train_idx, :);
        y_train = y(train_idx);
        cov_train = covariates(train_idx, :);

        X_test = X(test_idx, :);
        y_test = y(test_idx);
        cov_test = covariates(test_idx, :);

        best_C = choose_best_C_inner_cv(X_train, y_train, cov_train, C_grid, innerK);

        model = train_preprocessed_svm(X_train, y_train, cov_train, best_C);

        result = evaluate_preprocessed_svm(model, X_test, y_test, cov_test, 0, ...
            sprintf('Outer fold %d', fold));

        fold_id(fold) = fold;
        best_C_list(fold) = best_C;
        acc_list(fold) = result.accuracy;
        bacc_list(fold) = result.balanced_accuracy;
        auc_list(fold) = result.auc;
        macro_f1_list(fold) = result.macro_f1;
        sens0_list(fold) = result.sensitivity_class0;
        sens1_list(fold) = result.sensitivity_class1;
    end

    cv_results = table( ...
        fold_id, best_C_list, acc_list, bacc_list, auc_list, ...
        macro_f1_list, sens0_list, sens1_list, ...
        'VariableNames', {'Fold', 'Best_C', 'Accuracy', 'BalancedAccuracy', ...
        'AUC', 'MacroF1', 'Sensitivity_Visual', 'Sensitivity_DMNCEN'});

    cv_summary = table( ...
        mean(acc_list), std(acc_list), ...
        mean(bacc_list), std(bacc_list), ...
        mean(auc_list), std(auc_list), ...
        mean(macro_f1_list), std(macro_f1_list), ...
        mean(sens0_list), std(sens0_list), ...
        mean(sens1_list), std(sens1_list), ...
        'VariableNames', {'Accuracy_Mean', 'Accuracy_SD', ...
        'BalancedAccuracy_Mean', 'BalancedAccuracy_SD', ...
        'AUC_Mean', 'AUC_SD', ...
        'MacroF1_Mean', 'MacroF1_SD', ...
        'Sensitivity_Visual_Mean', 'Sensitivity_Visual_SD', ...
        'Sensitivity_DMNCEN_Mean', 'Sensitivity_DMNCEN_SD'});
end
