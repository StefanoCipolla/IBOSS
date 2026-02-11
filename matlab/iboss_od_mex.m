function result = iboss_od_mex(Z, Y, k, int_adj)
% IBOSS_OD_MEX  Information-Based Optimal Subdata Selection for linear regression.
%
%   result = iboss_od_mex(Z, Y, k)
%   result = iboss_od_mex(Z, Y, k, int_adj)
%
%   Same algorithm as iboss_od, but calls compiled MEX functions
%   (getIdx_mex, getIdxR_mex) for the index selection step, matching the
%   R package which calls the C++ routines in src/getIdx.cpp via .Call().
%
%   Before first use, compile the MEX files:
%       mex getIdx_mex.cpp
%       mex getIdxR_mex.cpp
%   (or run compile_mex.m)
%
%   INPUTS:
%       Z       - Covariate matrix (n x d) or covariate vector (n x 1).
%       Y       - Response vector (n x 1).
%       k       - Desired subdata size. For multivariate Z, k/(d*2) should
%                 be an integer; for univariate Z, k/2 should be an integer.
%       int_adj - (Optional) Logical flag for adjusting the intercept estimate
%                 using full-data means. Default is true.
%
%   OUTPUTS:
%       result  - A struct with fields:
%           beta      - Least squares coefficient estimates from the subdata.
%           se        - Standard errors of the coefficient estimates.
%           sigma     - Variance estimate for the error term.
%           index     - Indices of the selected subdata observations.
%           beta0_adj - Adjusted intercept estimate (if int_adj is true).
%
%   EXAMPLE:
%       See iboss_demo.m for a complete usage example.
%
%   REFERENCE:
%       Wang, H., Zhu, R., and Ma, P. (2018). Optimal Subsampling for Large
%       Sample Logistic Regression. Journal of the American Statistical
%       Association, 113(522), 829-844.

    if nargin < 4
        int_adj = true;
    end

    [n, d] = size(Z);

    % Handle univariate case (column vector)
    if d == 1
        if mod(k, 2) ~= 0
            warning('iboss_od_mex:nonIntegerR', ...
                'k/2 is not an integer; its floor is used as r.');
        end
        r = max(floor(k / 2), 1);
        idx_od = getIdx_mex(int32(r), Z);
        x_od = [ones(length(idx_od), 1), Z(idx_od)];
        if int_adj
            Zbar = mean(Z);
        end

    % Handle multivariate case (matrix)
    else
        if mod(floor(k / d), 2) ~= 0
            warning('iboss_od_mex:nonIntegerR', ...
                'k/d/2 is not an integer; its floor is used as r.');
        end
        r = floor(floor(k / d) / 2);
        if r < 1
            error('iboss_od_mex:rTooSmall', ...
                'Subdata size k is too small relative to the number of covariates d.');
        end

        % Select indices from first covariate
        idx_od = getIdx_mex(int32(r), Z(:, 1));

        % Select indices from remaining covariates, excluding already chosen
        for j = 2:d
            tmp = getIdxR_mex(int32(r), Z(:, j), idx_od);
            idx_od = [idx_od; tmp]; %#ok<AGROW>
        end

        x_od = [ones(length(idx_od), 1), Z(idx_od, :)];
        if int_adj
            Zbar = mean(Z, 1);
        end
    end

    % Fit least squares on subdata using eigendecomposition
    y_od = Y(idx_od);
    XtX = x_od' * x_od;
    [V, D] = eig(XtX);
    eigenvalues = diag(D);
    iI_od = V * diag(1 ./ eigenvalues) * V';
    beta_od = iI_od * (x_od' * y_od);

    % Adjusted intercept using full-data means
    beta0_adj = NaN;
    if int_adj
        Ybar = mean(Y);
        beta0_adj = Ybar - Zbar * beta_od(2:end);
    end

    % Residuals and variance estimate
    res_od = y_od - x_od * beta_od;
    sigma_od = sum(res_od.^2) / (n - d);

    % Standard errors
    se_od = diag(iI_od) * sigma_od;

    % Pack results into a struct
    result.beta      = beta_od;
    result.se        = se_od;
    result.sigma     = sigma_od;
    result.index     = idx_od;
    result.beta0_adj = beta0_adj;
end
