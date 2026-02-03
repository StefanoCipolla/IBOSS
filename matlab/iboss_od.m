function result = iboss_od(Z, Y, k, int_adj)
% IBOSS_OD  Information-Based Optimal Subdata Selection for linear regression.
%
%   result = iboss_od(Z, Y, k)
%   result = iboss_od(Z, Y, k, int_adj)
%
%   This function implements the IBOSS method (Wang, Zhu, and Ma, 2018) for
%   selecting an optimal subset of data for linear regression. It selects
%   observations with extreme covariate values to maximize the information
%   contained in the subdata.
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
            warning('iboss_od:nonIntegerR', ...
                'k/2 is not an integer; its floor is used as r.');
        end
        r = max(floor(k / 2), 1);
        idx_od = get_idx(r, Z);
        x_od = [ones(length(idx_od), 1), Z(idx_od)];
        if int_adj
            Zbar = mean(Z);
        end

    % Handle multivariate case (matrix)
    else
        if mod(floor(k / d), 2) ~= 0
            warning('iboss_od:nonIntegerR', ...
                'k/d/2 is not an integer; its floor is used as r.');
        end
        r = floor(floor(k / d) / 2);
        if r < 1
            error('iboss_od:rTooSmall', ...
                'Subdata size k is too small relative to the number of covariates d.');
        end

        % Select indices from first covariate
        idx_od = get_idx(r, Z(:, 1));

        % Select indices from remaining covariates, excluding already chosen
        for j = 2:d
            tmp = get_idx_restricted(r, Z(:, j), idx_od);
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


function idx = get_idx(r, z)
% GET_IDX  Find indices of the r smallest and r largest values in vector z.
%
%   idx = get_idx(r, z) returns a (2*r x 1) vector of indices:
%   the first r correspond to the smallest values, the last r to the largest.

    n = length(z);

    % Find the r-th smallest value (threshold for lower tail)
    z_sorted_asc = sort(z, 'ascend');
    yrl = z_sorted_asc(r);

    % Find the r-th largest value (threshold for upper tail)
    z_sorted_desc = sort(z, 'descend');
    yru = z_sorted_desc(r);

    % Collect indices by scanning (matches C++ behavior)
    locl = zeros(r, 1);
    locu = zeros(r, 1);
    jl = 0;
    ju = 0;
    for i = 1:n
        if z(i) <= yrl && jl < r
            jl = jl + 1;
            locl(jl) = i;
        end
        if z(i) >= yru && ju < r
            ju = ju + 1;
            locu(ju) = i;
        end
        if jl >= r && ju >= r
            break;
        end
    end
    idx = [locl; locu];
end


function idx = get_idx_restricted(r, z, del)
% GET_IDX_RESTRICTED  Find extreme-value indices, excluding already-selected ones.
%
%   idx = get_idx_restricted(r, z, del) returns (2*r x 1) indices of the
%   r smallest and r largest values in z, skipping any index present in del.

    n = length(z);
    del_sorted = sort(del);
    m = length(del_sorted);

    % Build filtered values (exclude del indices)
    y = zeros(n - m, 1);
    j = 1;
    k = 0;
    for i = 1:n
        if j <= m && del_sorted(j) == i
            j = j + 1;
        else
            k = k + 1;
            y(k) = z(i);
        end
    end
    y = y(1:k);

    % Find the r-th smallest threshold
    y_sorted_asc = sort(y, 'ascend');
    yrl = y_sorted_asc(r);

    % Find the r-th largest threshold
    y_sorted_desc = sort(y, 'descend');
    yru = y_sorted_desc(r);

    % Collect indices (excluding del)
    locl = zeros(r, 1);
    locu = zeros(r, 1);
    jl = 0;
    ju = 0;
    j = 1;
    for i = 1:n
        if j <= m && del_sorted(j) == i
            j = j + 1;
        else
            if z(i) <= yrl && jl < r
                jl = jl + 1;
                locl(jl) = i;
            end
            if z(i) >= yru && ju < r
                ju = ju + 1;
                locu(ju) = i;
            end
        end
        if jl >= r && ju >= r
            break;
        end
    end
    idx = [locl; locu];
end
