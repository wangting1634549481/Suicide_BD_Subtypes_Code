function X_z = apply_preprocessing(model, X, covariates)
    % Apply preprocessing learned from training data to new data.

    D = [ones(size(covariates, 1), 1), covariates];

    X_res = X - D * model.beta_cov;
    X_z = (X_res - model.mu_train) ./ model.sigma_train;
end

