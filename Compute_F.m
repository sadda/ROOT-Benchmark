function fun = Compute_F(x, h, w, c, bench)
    %COMPUTE_F - Computes the objective function for the moving peak benchmark.
    %
    % Syntax:  fun = Compute_F(x, h, w, c, bench)
    %
    % Inputs:
    %    x     - Points at which the function is to be evaluated. Matrix (n,d), where n is the number of trial points and d is the dimension.
    %    h     - Peak heights. Vector (m,1) for Benchmark 1 and matrix (m,d) for Benchmark 2.
    %    w     - Peak widths. Vector (m,1) for Benchmark 1 and matrix (m,d) for Benchmark 2.
    %    c     - Peak centers. Matrix (m,d) for Benchmark 1 and matrix (m,d) for Benchmark 2.
    %    bench - Either 'bench1' for Benchmark 1 or 'bench2' for benchmark2.
    %
    % Outputs:
    %    fun - The objective value. Vector (n,1).
    %
    % Example:
    %    fun = Compute_F(rand(10,3), rand(5,1), rand(10,1), rand(5,3), 'bench1')
    %    fun = Compute_F(rand(10,3), rand(5,3), rand(10,3), rand(5,3), 'bench2')
    %
    % Author: Lukas Adam
    % Paper: L. Adam, X. Yao: A Simple Yet Effective Approach to Robust Optimization Over Time
    % Email: adam@utia.cas.cz
    % July 2019; Last revision: 17-Jul-2019
    
    m   = size(c,1);
    n   = size(x,1);
    
    % Evaluate the objective. The inner if only minimizes the lenght of the for loop (makes the code faster).
    fun = -Inf(n,1);
    if strcmp(bench, 'bench1')
        if m<n
            for i=1:m
                fun = max(fun, h(i) - w(i)*vecnorm(x-c(i,:), 2, 2));
            end
        else
            for j=1:n
                fun(j) = max(h - w.*vecnorm(x(j,:)-c, 2, 2));
            end
        end
    elseif strcmp(bench, 'bench2')
        if m<n
            for i_m=1:m
                fun = max(fun, h(i_m,:) - w(i_m,:).*abs(x-c(i_m,:)));
            end
            fun = mean(fun,2);
        else
            for i_n=1:n
                fun(i_n) = mean(max(h - w.*abs(x(i_n,:)-c)));
            end
        end
    end
    
end

