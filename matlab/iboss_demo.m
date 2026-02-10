%% IBOSS Demo: Information-Based Optimal Subdata Selection
%  This script demonstrates how to use iboss_od() in MATLAB.
%  It replicates the example from the R IBOSS package documentation.

%% 1. Setup parameters
d = 50;                          % Number of covariates
beta_true = ones(d + 1, 1);     % True coefficients (intercept + 50 slopes)
n = 500000;                        % Full dataset size
k = 100;                         % Desired subdata size
noise_sd = 3;                    % Standard deviation of the error term

%% 2. Generate covariate correlation matrix
%  AR(1)-like structure with correlation 0.5
corr_val = 0.5;
Sigma = corr_val * ones(d, d) + (1 - corr_val) * eye(d);

%% 3. Generate data from a multivariate t-distribution (df = 2)
rng(0);  % Set seed for reproducibility

% Generate multivariate t with 2 degrees of freedom
% Method: X = Z / sqrt(chi2/df), where Z ~ N(0, Sigma)
df = 2;
Z_norm = mvnrnd(zeros(1, d), Sigma, n);
chi2_samples = chi2rnd(df, n, 1);
X = Z_norm ./ sqrt(chi2_samples / df);

% Generate response
mu = beta_true(1) + X * beta_true(2:end);
Y = mu + noise_sd * randn(n, 1);

%% 4. Run IBOSS
fprintf('Running IBOSS with n = %d, d = %d, k = %d ...\n', n, d, k);
fit = iboss_od(X, Y, k);

%% 5. Display results
fprintf('\n--- IBOSS Results ---\n');
fprintf('Number of observations selected: %d\n', length(fit.index));
fprintf('Estimated error variance (sigma): %.4f\n', fit.sigma);
fprintf('Adjusted intercept: %.4f\n\n', fit.beta0_adj);

% Compare estimated vs true coefficients
fprintf('Coefficient comparison (first 10 of %d):\n', d + 1);
fprintf('%-12s %-12s %-12s\n', 'Coeff', 'True', 'Estimated');
fprintf('%-12s %-12s %-12s\n', '-----', '----', '---------');
for i = 1:min(10, d + 1)
    if i == 1
        label = 'Intercept';
    else
        label = sprintf('beta_%d', i - 1);
    end
    fprintf('%-12s %-12.4f %-12.4f\n', label, beta_true(i), fit.beta(i));
end
if d + 1 > 10
    fprintf('... (%d more coefficients)\n', d + 1 - 10);
end

%% 6. Compare with full-data OLS
fprintf('\n--- Comparison with full-data OLS ---\n');
X_full = [ones(n, 1), X];
beta_full = X_full \ Y;

% Mean squared error of coefficients
mse_iboss = mean((fit.beta - beta_true).^2);
mse_full  = mean((beta_full - beta_true).^2);
fprintf('MSE of IBOSS coefficients:    %.6f\n', mse_iboss);
fprintf('MSE of full-data coefficients: %.6f\n', mse_full);
fprintf('Subdata fraction used:         %.1f%%\n', 100 * k / n);

%% 7. Visualize coefficient estimates
figure;
hold on;
plot(0:d, beta_true, 'k--', 'LineWidth', 1.5, 'DisplayName', 'True');
plot(0:d, fit.beta, 'ro', 'MarkerSize', 5, 'DisplayName', 'IBOSS');
plot(0:d, beta_full, 'b+', 'MarkerSize', 5, 'DisplayName', 'Full OLS');
hold off;
xlabel('Coefficient index');
ylabel('Value');
title('IBOSS vs Full OLS Coefficient Estimates');
legend('Location', 'best');
grid on;
