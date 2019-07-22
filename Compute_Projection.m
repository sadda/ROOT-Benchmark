function x = Compute_Projection(x, lb, ub, boundary)
    %COMPUTE_PROJECTION - Computes the projection onto the bounds.
    %
    % Syntax:  x = Compute_Projection(x, lb, ub, boundary)
    %
    % Inputs:
    %    x        - Data to be projected. Matrix (n,m).
    %    lb       - Lower bound for x. Matrix (n,m), vector (n,1) or scalar.
    %    ub       - Upper bound for x. Matrix (n,m), vector (n,1) or scalar.
    %    boundary - Specifies how the projection is computed. 1: standard projection, 2: reflection, 3: ignore.
    %
    % Outputs:
    %    x - Projected data. Matrix (n,m).
    %
    % Example:
    %    x_proj = Compute_Projection(rand(10,1), zeros(10,1), ones(10,1), 1)
    %
    % Author: Lukas Adam
    % Paper: L. Adam, X. Yao: A Simple Yet Effective Approach to Robust Optimization Over Time
    % Email: adam@utia.cas.cz
    % July 2019; Last revision: 17-Jul-2019
    
    if boundary == 1 % Project onto the bounds
        x = min(max(x, lb), ub);
    elseif boundary == 2 % Reflect from the bounds
        if isscalar(lb)
            lb = lb*ones(size(x));
        end
        if isscalar(ub)
            ub = ub*ones(size(x));
        end
        if ~isvector(x)
            if size(lb,2) < size(x,2)
                lb = repmat(lb, 1, size(x,2));
            end
            if size(ub,2) < size(x,2)
                ub = repmat(ub, 1, size(x,2));
            end
        end
        
        % The actual reflection
        ii    = x<lb;
        x(ii) = 2*lb(ii)-x(ii);
        ii    = x>ub;
        x(ii) = 2*ub(ii)-x(ii);
        
        % In case multiple reflections are necessary, repeat.
        if ~all(x(:)>=lb(:) & x(:)<=ub(:))
            x = Compute_Projection(x, lb, ub, boundary);
        end
    else % Do nothing
        
    end
end


