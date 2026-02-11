function overall_kur = kurtosis_summary(X, aggregate, do_log)
    % kurtosis_summary: Computes an aggregate kurtosis value across rows of X.
    %
    % Parameters:
    %   X               : A matrix of size (n_rows x n_columns)
    %   aggregate : Aggregation method: 'mean', 'min', 'max', or 'median'
    %   do_log      : Boolean flag to apply log10 transform before aggregation
    %
    % Returns:
    %   overall_kur : A scalar aggregate kurtosis value

    if nargin < 2
        aggregate = 'mean';
    end
    if nargin < 3
        do_log = true;
    end

    % Create symmetric version of X
    bigX = [ -X, X ];  % Horizontally concatenate -X and X

    % Compute row-wise kurtosis using Pearson definition (normal = 3)
    row_kurs = kurtosis(bigX, 0, 2);  % 0 => normal (not excess), 2 => operate along rows

    % Apply log transform if requested
    if do_log
        row_kurs_log = log10(row_kurs);
        switch lower(aggregate)
            case 'mean'
                agg_val = mean(row_kurs_log);
            case 'min'
                agg_val = min(row_kurs_log);
            case 'max'
                agg_val = max(row_kurs_log);
            case 'median'
                agg_val = median(row_kurs_log);
            otherwise
                error('%s is not a valid aggregation option', aggregate);
        end
        overall_kur = 10^agg_val;
    else
        switch lower(aggregate)
            case 'mean'
                overall_kur = mean(row_kurs);
            case 'min'
                overall_kur = min(row_kurs);
            case 'max'
                overall_kur = max(row_kurs);
            case 'median'
                overall_kur = median(row_kurs);
            otherwise
                error('%s is not a valid aggregation option', aggregate);
        end
    end
end
