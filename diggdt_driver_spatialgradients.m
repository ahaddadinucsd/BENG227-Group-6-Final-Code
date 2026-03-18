% Runs the 5 node IgG transport model simulation. Plots results on: fetal
% IgG subclass concentration, stromal heatmap spatial gradient, subclass
% specific spatial gradient.

%% ================= DigGDT Driver: Erdogan 5-node Transport Model =================
clear; clc; close all;

%% Load parameters
parameters_Erdogan_transport_5nodes;

%% Initialize model
x0 = [p.igg1, p.igg2, p.igg3, p.igg4, zeros(1,47)]; % 51 states total

%% Define time span
t_end = 40; % weeks
tspan = 10:0.5:t_end;

%% Solve ODE system
solution = ode15s(@(t,x) diggdt_transport_5nodes(t,x,p), tspan, x0);

%% Common definitions
Nsub = 4;
distance = 0 : p.dxstr : p.Lstr;

% Correct stromal indexing helper
idx_str = @(i,j) (j==1) .* (13 + i) + (j>=2) .* (35 + (j-2)*4 + i);

%% ================= Plot fetal IgG concentrations =================
figure(1)
ColorOrder = [
    0.2 0.4 0.8;
    0.2 0.7 0.8;
    0.5 0.6 0.9;
    0.4 0.5 0.7
];

for i = 1:4
    subplot(1,4,i)
    plot(solution.x, solution.y(26+i,:), 'LineWidth', 2, 'Color', ColorOrder(i,:));
    title(['IgG', num2str(i)]);
    xlabel('Gestational Age (weeks)');
    ylabel(['Fetal IgG', num2str(i), ' (M)']);
    ylim([0 solution.y(26+i,end)])
    grid on;
end

%% ================= Heatmap: Total IgG in stroma (all subclasses combined) =================
IgG_total = zeros(p.Nstr, length(solution.x));

% Collect TRUE stromal IgG across all subclasses at each node
for node = 1:p.Nstr
    for sub = 1:Nsub
        IgG_total(node,:) = IgG_total(node,:) + solution.y(idx_str(sub,node), :);
    end
end

% Smooth over gestation for prettier figure
IgG_total_smooth = smoothdata(IgG_total, 2, 'movmean', 5);

% Plot heatmap
figure;
pcolor(distance, solution.x, IgG_total_smooth');
shading interp;
set(gca, 'YDir', 'normal');
colormap(parula);
colorbar;
xlabel('Distance Across Stroma (Normalized)');
ylabel('Gestational Age (weeks)');
title('Total IgG Concentration Across Stroma Over Gestation');

%% ===== Spatial Gradient for Each IgG Subclass at Selected Gestational Ages =====
figure

weeks_to_plot = [15 25 35 40];

for sub = 1:4
    subplot(2,2,sub)
    hold on

    % Extract this subclass across stromal nodes
    IgG_sub = zeros(p.Nstr, length(solution.x));
    for node = 1:p.Nstr
        IgG_sub(node,:) = solution.y(idx_str(sub,node), :);
    end

    % Smooth over gestation for cleaner curves
    IgG_sub_smooth = smoothdata(IgG_sub, 2, 'movmean', 5);

    % Plot selected gestational ages
    for i = 1:length(weeks_to_plot)
        [~, tidx] = min(abs(solution.x - weeks_to_plot(i)));
        plot(distance, IgG_sub_smooth(:,tidx), 'LineWidth', 2)
    end

    xlabel('Distance Across Stroma (Normalized)')
    ylabel(['IgG', num2str(sub), ' Concentration'])
    title(['Spatial Gradient of IgG', num2str(sub)])
    legend('15 weeks','25 weeks','35 weeks','40 weeks','Location','best')
    grid on
end
