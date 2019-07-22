function [F_average, F_survival] = Compute_Metrics(value, T, delta)
    %COMPUTE_METRICS - Computes the evaluation metrics.
    %
    % Syntax:  [F_average, F_survival] = Compute_Metrics(value, T, delta)
    %
    % Inputs:
    %    value - Function values. Matrix (n,n_t), where n is the number of trial points and n_t the number of time instants.
    %    T     - Lenght of the averaging window. Scalar.
    %    delta - Threshold above which the survival metric is computed. Scalar or vector (n_delta,1).
    %
    % Outputs:
    %    F_average  - Averaging each row over time window metrics. Vector (n,1) or matrix (n,n_t).
    %               - If T is scalar, the time window lenght is T.
    %               - If T == 'all', then it computes the averages over all columns.
    %    F_survival - Survival metric computing how long each long each row stays above delta. Vector (n,1) or matrix (n,n_delta).
    %
    % Example:
    %    [F_average, F_survival] = Compute_Metrics(value, T, delta)
    %    [F_average, F_survival] = Compute_Metrics(value, 'all', delta)
    %    [F_average, ~]          = Compute_Metrics(value, T)
    %    [~, F_survival]         = Compute_Metrics(value, [], delta)
    %
    % Author: Lukas Adam
    % Paper: L. Adam, X. Yao: A Simple Yet Effective Approach to Robust Optimization Over Time
    % Email: adam@utia.cas.cz
    % July 2019; Last revision: 17-Jul-2019
    
    % Compute the average metric
    if isempty(T)
        F_average = [];
    elseif strcmp(T, 'all')
        F_average = cumsum(value,2) ./ (1:size(value,2));
    elseif isscalar(T)
        F_average = sum(value(:,1:T),2) ./ T;
    else
        error('T does not have the proper argument');
    end
    
    % Compute the survival metric
    if nargin == 3
        
        n          = size(value,1);
        n_delta    = length(delta);
        F_survival = zeros(n, n_delta);
        for i_delta=1:n_delta
            [ok, idx]               = max(value < delta(i_delta), [], 2);
            F_survival(:,i_delta)   = idx-1;
            F_survival(~ok,i_delta) = size(value,2);
        end
    end
    
end


