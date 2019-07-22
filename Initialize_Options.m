function opts = Initialize_Options(file_name)
    %INITIALIZE_OPTIONS - Initializes the options.
    %
    % Syntax:  opts = Initialize_Options(file_name)
    %
    % Inputs:
    %    file_name - Name of the paper from which should be called. String or empty.
    %
    % Outputs:
    %    opts - Structure of options. It has the following fields:
    %       n_eval      - Number of allowed function evalutions.
    %       boundary    - Described what happens when a parameter hits thte boundary. See Compute_Projection().
    %       bench       - Describes whether Benchmark 1 ('bench1') or Benchmark 2 ('bench2') is used.
    %       d           - Dimension of the problem.
    %       m           - Total number of peaks. Not along each dimension.
    %       T           - The time length on which the computation is run.
    %       T_add       - Number of additional time instances for evaluating performance.
    %       x_min       - Lower bound for x and peak centers c.
    %       x_max       - Upper bound for x and peak centers c.
    %       h_min       - Lower bound for peak heights h.
    %       h_max       - Upper bound for peak heights h.
    %       h_s_min     - Lower bound for peak height severities h_s.
    %       h_s_max     - Upper bound for peak height severities h_s.
    %       w_min       - Lower bound for peak widths w.
    %       w_max       - Upper bound for peak widths w.
    %       w_s_min     - Lower bound for peak width severities w_s.
    %       w_s_max     - Upper bound for peak width severities w_s.
    %            The following fields appear only for Benchmark 1.
    %       h_init      - Initial peak heights.
    %       w_init      - Initial peak widths.
    %       s           - Stepsizes for center movements.
    %       lambda      - Scaling parameter describing random and drift motion.
    %            The following fields appear only for Benchmark 2.
    %       theta_min   - Lower bound for angles theta.
    %       theta_max   - Upper bound for angles theta.
    %       theta_s_min - Lower bound for angle severities theta_s.
    %       theta_s_max - Upper bound for angle severities theta_s.
    %       gen_dim     - Described whether the peaks were generated uniformly along each dimension (with Cartesian product) or truly randomly.
    %       f_eval      - Evaluates the objective at time t for row vectors x. Call f_eval(t, x, h, w, c).
    %
    % Example:
    %    opts = Initialize_Options("Default1")
    %
    % Author: Lukas Adam
    % Paper: L. Adam, X. Yao: A Simple Yet Effective Approach to Robust Optimization Over Time
    % Email: adam@utia.cas.cz
    % July 2019; Last revision: 17-Jul-2019
    
    opts          = [];
    opts.n_eval   = 2500;
    opts.boundary = 1;
    
    if nargin < 1 || isempty(file_name)
        file_name = 'Default1';
    end
    
    if strcmp(file_name, 'Fu2012') || strcmp(file_name, 'Fu2013') || strcmp(file_name, 'Yazdani2017') || strcmp(file_name, 'Default1')
        
        opts.bench   = 'bench1';
        opts.d       = 2;
        opts.m       = 5;
        opts.x_min   = 0;
        opts.x_max   = 50;
        opts.h_min   = 30;
        opts.h_max   = 70;
        opts.h_init  = 50;
        opts.h_s_min = 1;
        opts.h_s_max = 10;
        opts.w_min   = 1;
        opts.w_max   = 12;
        opts.w_s_min = 0.1;
        opts.w_s_max = 1;
        opts.w_init  = 6;
        opts.F_eval  = @(t, x, h, w, c) Compute_F(x, h(:,t), w(:,t), squeeze(c(:,:,t)), 'bench1');
        
        if strcmp(file_name, 'Fu2012')
            opts.s      = [0; 0.5; 1; 1.5; 2];
            opts.lambda = 1;
            opts.T      = 20;
            opts.T_add  = 50;
        elseif strcmp(file_name, 'Fu2013') || strcmp(file_name, 'Guo2014')
            opts.s      = 1;
            opts.lambda = 1;
            opts.T      = 150;
            opts.T_add  = 50;
        elseif strcmp(file_name, 'Jin2013') || strcmp(file_name, 'Yazdani2017')
            opts.s      = 1;
            opts.lambda = 0;
            opts.T      = 150;
            opts.T_add  = 50;
        elseif strcmp(file_name, 'Default1')
            opts.s      = 1;
            opts.lambda = 0;
            opts.T      = 100;
            opts.T_add  = 50;
        end
        
    elseif strcmp(file_name, 'Fu2015') || strcmp(file_name, 'Hernandez2018') || strcmp(file_name, 'Default2')
        
        opts.bench       = 'bench2';
        opts.d           = 2;
        opts.m           = 25;
        opts.x_min       = -25;
        opts.x_max       = 25;
        opts.h_min       = 30;
        opts.h_max       = 70;
        opts.h_s_min     = 5;
        opts.h_s_max     = 5;
        opts.w_min       = 1;
        opts.w_max       = 13;
        opts.w_s_min     = 0.5;
        opts.w_s_max     = 0.5;
        opts.theta_min   = -pi;
        opts.theta_max   = pi;
        opts.theta_s_min = 1;
        opts.theta_s_max = 1;
        opts.F_eval      = @(t, x, h, w, c) Compute_F(x, squeeze(h(:,:,t)), squeeze(w(:,:,t)), squeeze(c(:,:,t)), 'bench2');
        
        if strcmp(file_name, 'Fu2015') || strcmp(file_name, 'Hernandez2018')
            opts.gen_dim = 1;
            opts.T       = 100;
            opts.T_add   = 100;
        elseif strcmp(file_name, 'Default2')
            opts.gen_dim = 0;
            opts.T       = 100;
            opts.T_add   = 100;
        end
        
    else
        
        error('Wrong file_name. Change it or run "opts = Initialize_Options(''Default1'');"');
        
    end
    
end


