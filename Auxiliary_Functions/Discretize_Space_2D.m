function [X, Y] = Discretize_Space_2D(n, x_min, x_max, y_min, y_max)
    %DISCRETIZE_SPACE_2D - Computes the uniform discretization of [x_min,x_max]*[y_min,y_max] with approximately n points.
    %
    % Syntax:  [X, Y] = Discretize_Space_2D(n, x_min, x_max, y_min, y_max)
    %
    % Inputs:
    %    n     - The approximate number of points. Scalar.
    %    x_min - Lower bound for coordinate x. Scalar.
    %    x_max - Upper bound for coordinate x. Scalar.
    %    y_min - Lower bound for coordinate y. Scalar.
    %    y_max - Upper bound for coordinate y. Scalar.
    %
    % Outputs:
    %    X - Coordinate x of the discretization. Vector (?,1).
    %    Y - Coordinate y of the discretization. Vector (?,1).
    %
    % Example:
    %    [X, Y] = Discretize_Space_2D(100, 0, 1, 0, 2)
    %
    % Author: Lukas Adam
    % Paper: L. Adam, X. Yao: A Simple Yet Effective Approach to Robust Optimization Over Time
    % Email: adam@utia.cas.cz
    % July 2019; Last revision: 17-Jul-2019
    
    n_x   = floor(sqrt(n));
    n_y   = floor(n/n_x);
    x     = linspace(x_min, x_max, n_x);
    y     = linspace(y_min, y_max, n_y);
    [X,Y] = meshgrid(x,y);
    X     = X(:);
    Y     = Y(:);
    
end