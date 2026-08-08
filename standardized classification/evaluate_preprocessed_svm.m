
function result = evaluate_preprocessed_svm(model, X, y, covariates, nPerm, dataset_name)
    % Evaluate trained SVM classifier.

    y = y(:);

    X_z = apply_preprocessing(model, X, covariates);

    [y_pred, score] = predict(model.svm_model, X_z);

    class_names = model.svm_model.ClassNames;
    idx_class1 = find(class_names == 1);

    score_class1 = score(:, idx_class1);

    cm = confusionmat(y, y_pred, 'Order', [0; 1]);

    % cm:
    % row = true class
    % column = predicted class
    %
    % class 0 = visual subtype
    % class 1 = DMN-CEN subtype

    n_total = sum(cm(:));

    accuracy = sum(diag(cm)) / n_total;

    sensitivity_class0 = cm(1, 1) / sum(cm(1, :));
    sensitivity_class1 = cm(2, 2) / sum(cm(2, :));

    balanced_accuracy = mean([sensitivity_class0, sensitivity_class1]);

    precision_class0 = cm(1, 1) / sum(cm(:, 1));
    precision_class1 = cm(2, 2) / sum(cm(:, 2));

    f1_class0 = 2 * precision_class0 * sensitivity_class0 / ...
        (precision_class0 + sensitivity_class0);

    f1_class1 = 2 * precision_class1 * sensitivity_class1 / ...
        (precision_class1 + sensitivity_class1);

    macro_f1 = mean([f1_class0, f1_class1], 'omitnan');

    % AUC
    try
        [~, ~, ~, auc] = perfcurve(y, score_class1, 1);
    catch
        auc = NaN;
    end

    % Permutation p value for balanced accuracy
    if nPerm > 0
        perm_bacc = zeros(nPerm, 1);

        for p = 1:nPerm
            y_perm = y(randperm(length(y)));
            cm_perm = confusionmat(y_perm, y_pred, 'Order', [0; 1]);

            sens0_perm = cm_perm(1, 1) / sum(cm_perm(1, :));
            sens1_perm = cm_perm(2, 2) / sum(cm_perm(2, :));

            perm_bacc(p) = mean([sens0_perm, sens1_perm]);
        end

        p_perm = (sum(perm_bacc >= balanced_accuracy) + 1) / (nPerm + 1);
    else
        p_perm = NaN;
    end

    result = struct();
    result.dataset = dataset_name;
    result.confusion_matrix = cm;
    result.accuracy = accuracy;
    result.balanced_accuracy = balanced_accuracy;
    result.auc = auc;
    result.macro_f1 = macro_f1;
    result.sensitivity_class0 = sensitivity_class0;
    result.sensitivity_class1 = sensitivity_class1;
    result.f1_class0 = f1_class0;
    result.f1_class1 = f1_class1;
    result.permutation_p = p_perm;

    fprintf('\n===== %s =====\n', dataset_name);
    fprintf('Accuracy: %.4f\n', accuracy);
    fprintf('Balanced accuracy: %.4f\n', balanced_accuracy);
    fprintf('AUC: %.4f\n', auc);
    fprintf('Macro-F1: %.4f\n', macro_f1);
    fprintf('Sensitivity visual subtype: %.4f\n', sensitivity_class0);
    fprintf('Sensitivity DMN-CEN subtype: %.4f\n', sensitivity_class1);

    if nPerm > 0
        fprintf('Permutation p for balanced accuracy: %.5f\n', p_perm);
    end

    fprintf('Confusion matrix:\n');
    disp(cm);
end