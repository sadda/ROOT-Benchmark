function [h, w, c] = Compute_Evolution(opts)
    %COMPUTE_EVOLUTION - Computes the evolution of the heights, widths and centers for the moving peak benchmark.
    %
    % Syntax:  [h, w, c] = Compute_Evolution(opts)
    %
    % Inputs:
    %    opts - Options. Struct described in Initialize_Options().
    %
    % Outputs:
    %    h - Peak heights. Matrix (m,n_t) for Benchmark 1 and matrix (m,d,n_t) for Benchmark 2.
    %    w - Peak widths. Matrix (m,n_t) for Benchmark 1 and matrix (m,d,n_t) for Benchmark 2.
    %    c - Peak centers. Matrix (m,d,n_t) for Benchmark 1 and matrix (m,d,n_t) for Benchmark 2.
    %
    % Example:
    %    [h, w, c] = Compute_Evolution(opts)
    %
    % Author: Lukas Adam
    % Paper: L. Adam, X. Yao: A Simple Yet Effective Approach to Robust Optimization Over Time
    % Email: adam@utia.cas.cz
    % July 2019; Last revision: 17-Jul-2019
    
    d     = opts.d;
    m     = opts.m;
    T_enh = opts.T+opts.T_add;
    if strcmp(opts.bench, 'bench1')
        
        s      = opts.s;
        lambda = opts.lambda;
        
        %% Generate initial random variables
        
        % Propagate the step lenghts to handle that sometimes they are different for different peaks.
        if length(s) == 1
            s = s*ones(m, 1);
        end
        
        % Initialize the height and width severities
        h_s = opts.h_s_min + (opts.h_s_max-opts.h_s_min)*rand(m,1);
        w_s = opts.w_s_min + (opts.w_s_max-opts.w_s_min)*rand(m,1);
        
        % Get the trend vector
        v   = -0.5 + rand(m,d);
        v   = s.*v./vecnorm(v, 2, 2);
        
        %% Compute evolution of random variables
        
        h = zeros(m,T_enh);
        w = zeros(m,T_enh);
        c = zeros(m,d,T_enh);
        
        % Get the values at the initial time instant
        h(:,1)   = opts.h_init;
        w(:,1)   = opts.w_init;
        c(:,:,1) = opts.x_min + (opts.x_max-opts.x_min)*rand(m,d);
        
        % Get the value at the next time instant
        for t=1:T_enh-1
            
            r          = -0.5 + rand(m,d);
            r          = s.*r./vecnorm(r, 2, 2);
            v          = (1-lambda)*r+lambda*v;
            v          = s.*v./vecnorm(v, 2, 2);
            v(s==0,:)  = 0;
            
            h(:,t+1)   = Compute_Projection(h(:,t) + h_s.*randn(m,1), opts.h_min, opts.h_max, opts.boundary);
            w(:,t+1)   = Compute_Projection(w(:,t) + w_s.*randn(m,1), opts.w_min, opts.w_max, opts.boundary);
            c(:,:,t+1) = Compute_Projection(c(:,:,t) + v, opts.x_min, opts.x_max, opts.boundary);
            
        end
        
    elseif strcmp(opts.bench, 'bench2')
        
        %% Generate initial random variables
        
        % Initialize the height, width and angle severities
        h_s     = opts.h_s_min + (opts.h_s_max-opts.h_s_min)*rand(m,d);
        w_s     = opts.w_s_min + (opts.w_s_max-opts.w_s_min)*rand(m,d);
        theta_s = opts.theta_s_min + (opts.theta_s_max-opts.theta_s_min)*rand(1,1);
        
        %% Compute evolution of random variables
        
        h     = zeros(m,d,T_enh);
        w     = zeros(m,d,T_enh);
        c     = zeros(m,d,T_enh);
        theta = zeros(d-1,T_enh);
        
        % Get the values at the initial time instant
        h(:,:,1) = opts.h_min + (opts.h_max-opts.h_min)*rand(m,d);
        w(:,:,1) = opts.w_min + (opts.w_max-opts.w_min)*rand(m,d);
        if opts.gen_dim
            % This one is a bit more complicated. But it should be correct
            m_dim = round(m^(1/d));
            ii    = unique(nchoosek(repmat(1:m_dim, 1, m_dim), d), 'rows');
            c_aux = opts.x_min + (opts.x_max-opts.x_min)*rand(m_dim,d);
            for i_d = 1:d
                for i_m = 1:m_dim
                    jj          = ii(:,i_d) == i_m;
                    c(jj,i_d,1) = c_aux(i_m,i_d);
                end
            end
        else
            c(:,:,1) = opts.x_min + (opts.x_max-opts.x_min)*rand(m,d);
        end
        theta(:,1) = 0;
        
        % Get the value at the initial time instant
        for t=1:T_enh-1
            
            % Compute the rotation matrix
            R = eye(d);
            for i_d=1:d-1
                ii           = i_d:i_d+1;
                R_aux        = eye(d);
                R_aux(ii,ii) = [cos(theta(i_d,t)), sin(theta(i_d,t)); -sin(theta(i_d,t)), cos(theta(i_d,t))];
                R            = R*R_aux;
            end
            
            % Perform the rotation
            for i_m=1:m
                c(i_m,:,t+1) = R*squeeze(c(i_m,:,t))';
            end
            
            h(:,:,t+1)   = Compute_Projection(h(:,:,t) + h_s.*randn(m,d), opts.h_min, opts.h_max, opts.boundary);
            w(:,:,t+1)   = Compute_Projection(w(:,:,t) + w_s.*randn(m,d), opts.w_min, opts.w_max, opts.boundary);
            theta(:,t+1) = Compute_Projection(theta(:,t) + theta_s.*randn(d-1,1), opts.theta_min, opts.theta_max, opts.boundary);
            c(:,:,t+1)   = Compute_Projection(c(:,:,t+1), opts.x_min, opts.x_max, 1);
            
        end
        
    else
        error('Benchmark name not supported. Run "opts = Initialize_Options(''Default1'')"');
    end
    
end

