clear all;
close all;

addpath('./Auxiliary_Functions');

% Create the folder for results
folder_name = 'Results';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

% Set the parameters including the number of repetitions, thresholds delta and max function evaluations for the gradient method
n_rep      = 5000;
T_range    = 20;
delta_all  = [40 50];
opts_optim = optimoptions('fmincon', 'display', 'off', 'MaxFunctionEvaluations', 200, 'Algorithm', 'sqp');

% Set the parameters for evalution
eval_T_point = 30;
eval_T_int   = 20:100;
eval_S       = [2 6];

% Set the file names
file_names = {'Default1', 'Default2', 'Fu2013', 'Fu2015'};

% Loop over all datasets
for i_name = 1:length(file_names)
    
    % Initialize the options
    file_name = file_names{i_name};
    opts      = Initialize_Options(file_name);
    T         = opts.T;
    
    % Prepare data to save results to
    survival_mesh = zeros(T,length(delta_all),n_rep);
    survival_tmo  = zeros(T,length(delta_all),n_rep);
    survival_rob  = zeros(T,length(delta_all),n_rep);
    average_mesh  = zeros(T,T_range,n_rep);
    average_tmo   = zeros(T,T_range,n_rep);
    average_rob   = zeros(T,T_range,n_rep);
    center_mesh   = zeros(T,opts.m,n_rep);
    center_tmo    = zeros(T,opts.m,n_rep);
    center_rob    = zeros(T,opts.m,n_rep);
    
    % Loop over repetitions
    for i_rep=1:n_rep
        
        if i_rep == 1 || mod(i_rep, 50) == 0
            fprintf('Running %4d out of %4d. Setting = %s\n', i_rep, n_rep, file_name);
        end
        
        %% Prepare the data before optimization
        
        % Generate random evolution of the peaks
        [h, w, c] = Compute_Evolution(opts);
        
        % Discretize the space into mesh [X, Y]
        [X, Y] = Discretize_Space_2D(opts.n_eval, opts.x_min, opts.x_max, opts.x_min, opts.x_max);
        
        % Evaluate the function at the mesh for all times (during optimization no future values will be used)
        T_enh = opts.T+opts.T_add;
        Z_all = zeros(T_enh,length(X));
        for t=1:T_enh
            Z_all(t,:) = opts.F_eval(t, [X,Y], h, w, c);
        end
        
        %% Find the optimal value on the mesh to get the Mesh solution (note that at time t, we use only present values)
        
        i_mesh = zeros(T_enh,1);
        for t=1:T_enh
            [~,i_mesh(t)] = max(Z_all(t,:));
        end
        
        %% Improve the solution by using gradient descent to get the Time-optimal solution
        
        x_tmo = zeros(T_enh,2);
        for t=1:T_enh
            
            % Randomly sample opts.n_eval-opts_optim.MaxFunctionEvaluations points on the mesh
            Z_red     = Z_all(t,:)';
            ii        = randperm(opts.n_eval)';
            ii        = ii(1:(opts.n_eval-opts_optim.MaxFunctionEvaluations));
            
            % Select the best point on the smaller mesh
            [~,i_max] = max(Z_red(ii));
            i_max     = ii(i_max);
            
            % Improve the best point using the remaining opts_optim.MaxFunctionEvaluations function evalutions
            x0         = [X(i_max); Y(i_max)];
            objective  = @(xy) -opts.F_eval(t, xy', h, w, c);
            x          = fmincon(objective,x0,[],[],[],[],[opts.x_min;opts.x_min],[opts.x_max;opts.x_max],[],opts_optim);
            x_tmo(t,:) = x';
            
        end
        
        % Compute the function values for the newly found solution
        Z_tmo = zeros(T_enh,T_enh);
        for t=1:T_enh
            Z_tmo(t,:) = opts.F_eval(t, x_tmo, h, w, c);
        end
        
        %% Find the robust solution by averaging on the mesh
        
        i_rob = zeros(T_enh,1);
        for t=1:T
            
            % Average on the mesh with window size n_conv
            n_conv       = 7;
            Z_mod        = Z_all(t,:)';
            Z_mod        = reshape(Z_mod, sqrt(size(Z_mod,1)), []);
            Z_rob        = conv2(Z_mod, ones(n_conv), 'same') ./ conv2(ones(size(Z_mod)), ones(n_conv), 'same');
            Z_rob        = Z_rob(:);
            
            % Select the solution with the largest value
            [~,i_rob(t)] = max(Z_rob);
            
        end
        
        %% Evaluate the solutions
        
        for t=1:T
            
            % Evaluate the average and survival metrics
            [average_mesh(t,:,i_rep), survival_mesh(t,:,i_rep)] = Compute_Metrics(Z_all(t:t+T_range-1,i_mesh(t))', 'all', delta_all);
            [average_tmo(t,:,i_rep), survival_tmo(t,:,i_rep)]   = Compute_Metrics(Z_tmo(t:t+T_range-1,t)', 'all', delta_all);
            [average_rob(t,:,i_rep), survival_rob(t,:,i_rep)]   = Compute_Metrics(Z_all(t:t+T_range-1,i_rob(t))', 'all', delta_all);
            
            % Compute distance of the solutions to centers
            if strcmp(opts.bench, 'bench1')
                h_sorted = sort(h(:,t), 'descend');
            elseif strcmp(opts.bench, 'bench2')
                h_sorted = mean(sort(h(:,:,t), 'descend'),2);
            end
            center_mesh(t,:,i_rep) = abs(Z_all(t,i_mesh(t)) - h_sorted);
            center_tmo(t,:,i_rep)  = abs(Z_tmo(t,t) - h_sorted);
            center_rob(t,:,i_rep)  = abs(Z_all(t,i_rob(t)) - h_sorted);
            
        end
        
    end
    
    %% Save results
    
    % Average all results with respect to the repetions
    survival_mesh_mean = mean(survival_mesh,3)';
    survival_tmo_mean  = mean(survival_tmo,3)';
    survival_rob_mean  = mean(survival_rob,3)';
    average_mesh_mean  = mean(average_mesh,3)';
    average_tmo_mean   = mean(average_tmo,3)';
    average_rob_mean   = mean(average_rob,3)';
    center_mesh_mean   = mean(center_mesh,3);
    center_tmo_mean    = mean(center_tmo,3);
    center_rob_mean    = mean(center_rob,3);
    
    % Save the results while omitting some variables to make the files smaller
    variables_save = who;
    variables_excl = {'fig', 'survival_mesh', 'survival_tmo', 'survival_rob', 'average_mesh', 'average_tmo', 'average_rob', 'center_mesh', 'center_tmo', 'center_rob', 'Z_all', 'Z_tmo'};
    variables_save = variables_save(~ismember(variables_save, variables_excl));
    save(fullfile(folder_name, sprintf('Results_%s.mat', file_name)), variables_save{:});
    
    %% Analyze results
    
    % Plot survival time
    for i_delta=1:length(delta_all)
        fig = figure();
        plot([survival_mesh_mean(i_delta,:); survival_tmo_mean(i_delta,:); survival_rob_mean(i_delta,:)]');
        xlabel('Time');
        ylabel('Averaged survival (max 20)');
        ylim([-1 20]);
        legend({'Mesh', 'Time-optimal', 'Robust'});
        saveas(fig, fullfile(folder_name, sprintf('Results_%s_Survival_Delta=%d.jpg', file_name, delta_all(i_delta))));
    end
    
    % Plot average survival at time t_eval
    fig = figure();
    plot([average_mesh_mean(:,eval_T_point)'; average_tmo_mean(:,eval_T_point)'; average_rob_mean(:,eval_T_point)']');
    xlabel('Time');
    ylabel('Average f');
    ylim([30 70]);
    legend({'Mesh', 'Time-optimal', 'Robust'});
    saveas(fig, fullfile(folder_name, sprintf('Results_%s_Average.jpg', file_name)));
    
    %% Display the average and survival function for time interval T_int and distance from the first center
    
    res_aver = mean(average_tmo_mean(eval_S,eval_T_int),2);
    res_surv = mean(survival_tmo_mean(:,eval_T_int),2);
    res_cent = [mean(center_mesh_mean(eval_T_int,1)), mean(center_tmo_mean(eval_T_int,1)), mean(center_rob_mean(eval_T_int,1))];
    
    disp(res_aver);
    disp(res_surv);
    disp(res_cent);
    
    %% Save csv with data
    
    for i_delta=1:length(delta_all)
        csvwrite(fullfile(folder_name, sprintf('Results_%s_Survival_Delta=%d.csv', file_name, delta_all(i_delta))), [1:T; survival_mesh_mean(i_delta,:); survival_tmo_mean(i_delta,:); survival_rob_mean(i_delta,:)]');
    end
    csvwrite(fullfile(folder_name, sprintf('Results_%s_Average.csv', file_name)), [1:T_range; average_mesh_mean(:,eval_T_point)'; average_tmo_mean(:,eval_T_point)'; average_rob_mean(:,eval_T_point)']');
    
end




