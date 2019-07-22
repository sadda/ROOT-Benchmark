clear all;
close all;

%% Compute the maximal peak height over time

n      = 1e5;
T      = 20;
h_init = 50;
m_all  = [5, 25];

% Perform the evolution over time
n_m   = length(m_all);
h_max = zeros(T, n_m);
for i_m = 1:n_m
    
    m = m_all(i_m);
    h = h_init*ones(T,m,n);
    
    for t=1:T-1
        
        h(t+1,:,:) = squeeze(h(t,:,:)) + 5.5*randn(m,n);
        h(t+1,:,:) = Compute_Projection(h(t+1,:,:), 30, 70, 1);
        
    end
    
    h_max(:,i_m) = mean(squeeze(max(h, [], 2)), 2);
    
end

fig = figure();
plot(h_max);
xlabel('Time');
ylabel('Average max height');
ylim([50 70])
saveas(fig, 'Results_Max.jpg');

csvwrite('Results_Max.csv', [(1:T); h_max']');

%% Compare the uniform generation on the square and circle

n = 1e7;

% Geneate on the square
r1         = -0.5+rand(n,2);
r1         = r1./vecnorm(r1, 2, 2);
ii         = r1(:,2) <= 0;
theta1     = acos(r1(:,1));
theta1(ii) = 2*pi - theta1(ii);

% Generate on the circle
r2         = randn(n,2);
r2         = r2./vecnorm(r2, 2, 2);
ii         = r2(:,2) <= 0;
theta2     = acos(r2(:,1));
theta2(ii) = 2*pi - theta2(ii);

% Compute the histogram data
n_bins   = 100;
edges    = linspace(0,2*pi,n_bins+1);
bins_mid = edges(1:end-1) + diff(edges);
N1       = histcounts(theta1,edges);
N2       = histcounts(theta2,edges);
N1       = N1/sum(N1.*diff(edges));
N2       = N2/sum(N2.*diff(edges));

fig = figure();
plot(bins_mid, [N1; N2]);
saveas(fig, 'Results_Counts.jpg');

csvwrite('Results_Counts.csv', [bins_mid; N1; N2]');

%% Make video

addpath('./Auxiliary_Functions');

video_name = 'Results_Evolution.avi';
file_name  = 'Default1';
opts       = Initialize_Options(file_name);
[h, w, c]  = Compute_Evolution(opts);

% Compute the function values
T_range = 20;
[X, Y]  = Discretize_Space_2D(opts.n_eval, opts.x_min, opts.x_max, opts.x_min, opts.x_max);
Z_all   = zeros(T_range,length(X));

for t=1:T_range
    
    Z_all(t,:) = opts.F_eval(t, [X,Y], h, w, c);
    
end
Z = Z_all';

% Start figure
figure();
set(gca,'nextplot','replacechildren');

% Start video
v = VideoWriter(video_name);
v.FrameRate = 2;
open(v);

% Save all frames
for t = 1:T_range
    
    X_r = reshape(X, sqrt(size(X,1)), []);
    Y_r = reshape(Y, sqrt(size(X,1)), []);
    Z_r = reshape(Z(:,t), sqrt(size(X,1)), []);
    
    surf(X_r, Y_r, Z_r);
    view(2);
    shading interp;
    colormap jet;
    colorbar;
    xlim([min(X) max(X)]);
    ylim([min(Y) max(Y)]);
    axis off;
    caxis([25 75]);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
end

% Close the video
close(v);


