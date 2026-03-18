%% run_preeclampsia_comparison.m
%% Runs the preeclampsia model under all PE severity patterns and
%% compares fetal IgG transfer.  Also overlays the healthy baseline
%% and endometriosis model for three-way disease comparison.
%%
%% Figures produced:
%%   Figure 1 - Stromal diffusivity profiles for each PE pattern
%%   Figure 2 - Fetal IgG subclass trajectories per PE pattern
%%   Figure 3 - F:M ratios at delivery: PE patterns + healthy + endo
%%   Figure 4 - Stromal IgG concentration profiles at week 28
%%
%% Dolatshahi Lab - Preeclampsia Extension

clear; clc; close all;

%% PE patterns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pe_patterns = {'mild_pe', 'severe_pe', 'early_onset', 'late_onset'};
pe_labels   = {'Mild PE', 'Severe PE', 'Early-Onset PE', 'Late-Onset PE'};
colors_pe   = [0.20 0.60 0.86;   % blue
               0.86 0.20 0.20;   % red
               0.60 0.10 0.10;   % dark red
               0.95 0.60 0.10];  % orange

%% Storage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results_pe = struct();

for pidx = 1:length(pe_patterns)
    p = struct();
    p.pe_pattern = pe_patterns{pidx};
    parameters_preeclampsia;         % populates p

    opts   = odeset('RelTol',1e-8,'AbsTol',1e-12,'NonNegative',1:51);
    odefun = @(t,x) diggdt_transport_5nodes_fibrosis(t, x, p);
    [T, X] = ode15s(odefun, tspan, x0, opts);

    results_pe(pidx).pattern  = pe_patterns{pidx};
    results_pe(pidx).label    = pe_labels{pidx};
    results_pe(pidx).T        = T;
    results_pe(pidx).X        = X;
    results_pe(pidx).Dstr     = p.Dstr;
    results_pe(pidx).ustr     = p.ustr;
    results_pe(pidx).p        = p;
    results_pe(pidx).IgG_fetal    = X(:, 27:30);
    results_pe(pidx).IgG_maternal = [p.igg1, p.igg2, p.igg3, p.igg4];
    results_pe(pidx).FM_ratio = X(end, 27:30) ./ [p.igg1,p.igg2,p.igg3,p.igg4];

    fprintf('Completed: %s\n', pe_labels{pidx});
end

%% Also run healthy baseline and endometriosis for comparison %%%%%%%%%%
p = struct();
parameters_improved;
p.Dstr = p.Dstr * ones(1, p.Nstr-1);   % convert scalar to uniform vector
p_h = p;
opts   = odeset('RelTol',1e-8,'AbsTol',1e-12,'NonNegative',1:51);
[T_h, X_h] = ode15s(@(t,x) diggdt_transport_5nodes_fibrosis(t,x,p_h), ...
                     tspan, x0, opts);
FM_healthy = X_h(end,27:30) ./ [p_h.igg1,p_h.igg2,p_h.igg3,p_h.igg4];

p_e = struct();
p_e.fibrosis_pattern = 'patchy_focal';
parameters_endometriosis_fibrosis;
p_e = p;
[T_e, X_e] = ode15s(@(t,x) diggdt_transport_5nodes_fibrosis(t,x,p_e), ...
                     tspan, x0, opts);
FM_endo = X_e(end,27:30) ./ [p_e.igg1,p_e.igg2,p_e.igg3,p_e.igg4];

subclass_names = {'IgG1','IgG2','IgG3','IgG4'};
Nstr = 5;
node_positions    = linspace(0,1,Nstr);
interval_midpoints = (node_positions(1:end-1)+node_positions(2:end))/2;
idx_str_fn = @(i,j) (j==1).*(13+i) + (j>=2).*(35+(j-2)*4+i);

%% =========================================================================
%% FIGURE 1: Stromal diffusivity profiles
%% =========================================================================
figure(1); clf; set(gcf,'Position',[100 100 720 420]);
hold on;
for pidx = 1:length(pe_patterns)
    plot(interval_midpoints, results_pe(pidx).Dstr, '-o', ...
        'Color', colors_pe(pidx,:), 'LineWidth', 2.5, ...
        'MarkerFaceColor', colors_pe(pidx,:), 'MarkerSize', 9, ...
        'DisplayName', sprintf('%s  (u=%.2f)', ...
            pe_labels{pidx}, results_pe(pidx).ustr));
end
plot(interval_midpoints, ones(1,4), '--k', 'LineWidth', 1.5, ...
    'DisplayName', 'Healthy baseline');
hold off;
xlabel('Stromal Position  (STB \leftarrow 0 ———— 1 \rightarrow EC)', 'FontSize',12);
ylabel('Effective Diffusivity  D', 'FontSize',12);
title('Preeclampsia: Stromal Diffusivity Profiles', 'FontSize',13,'FontWeight','bold');
legend('Location','northeast','FontSize',9);
ylim([0 1.2]); grid on; box on;
xticks(interval_midpoints);
xticklabels({'Node 1-2','Node 2-3','Node 3-4','Node 4-5'});

%% =========================================================================
%% FIGURE 2: Fetal IgG trajectories
%% =========================================================================
figure(2); clf; set(gcf,'Position',[100 600 900 680]);
for sc = 1:4
    subplot(2,2,sc); hold on;
    for pidx = 1:length(pe_patterns)
        T = results_pe(pidx).T;
        C = results_pe(pidx).IgG_fetal(:,sc);
        plot(T, C*1e6, '-', 'Color', colors_pe(pidx,:), 'LineWidth', 2, ...
            'DisplayName', pe_labels{pidx});
    end
    plot(T_h, X_h(:,26+sc)*1e6, '--k', 'LineWidth',1.5, 'DisplayName','Healthy');
    plot(T_e, X_e(:,26+sc)*1e6, ':',  'Color',[0.4 0.4 0.4], ...
         'LineWidth',1.5, 'DisplayName','Endometriosis');
    hold off;
    xlabel('Gestational Age (weeks)','FontSize',11);
    ylabel('Fetal IgG (\muM)','FontSize',11);
    title(subclass_names{sc},'FontSize',12,'FontWeight','bold');
    if sc==1; legend('Location','northwest','FontSize',8); end
    grid on; box on;
    xlim([10, max([results_pe(end).T(end), T_h(end), T_e(end)])]);
end
sgtitle('Fetal IgG Accumulation — Preeclampsia vs Comparators','FontSize',14);

%% =========================================================================
%% FIGURE 3: F:M ratios at delivery — grouped bar across all conditions
%% =========================================================================
figure(3); clf; set(gcf,'Position',[840 100 820 500]);

all_labels = [pe_labels, {'Healthy','Endometriosis'}];
all_FM     = vertcat( ...
    cell2mat(arrayfun(@(x) x.FM_ratio, results_pe, 'UniformOutput', false)'), ...
    FM_healthy, FM_endo );

bar_colors = [colors_pe; 0 0 0; 0.5 0.5 0.5];
b = bar(all_FM, 'grouped');
for k = 1:4
    b(k).FaceColor = [0.3+0.1*k, 0.5-0.05*k, 0.8-0.1*k];
    b(k).DisplayName = subclass_names{k};
end
set(gca,'XTickLabel', all_labels, 'FontSize',10);
xtickangle(20);
ylabel('Fetal:Maternal IgG Ratio at Delivery','FontSize',12);
title('F:M Ratios at Delivery — PE Patterns vs Healthy vs Endometriosis', ...
      'FontSize',13,'FontWeight','bold');
legend(subclass_names,'Location','northeast','FontSize',10);
yline(1.0,'--k','F:M = 1','LineWidth',1.5,'LabelHorizontalAlignment','left');
grid on; box on;

%% =========================================================================
%% FIGURE 4: Stromal concentration profiles at week 28
%% =========================================================================
target_week = 28;
figure(4); clf; set(gcf,'Position',[840 620 900 480]);
for sc = 1:4
    subplot(1,4,sc); hold on;
    for pidx = 1:length(pe_patterns)
        T = results_pe(pidx).T;
        X = results_pe(pidx).X;
        [~,tidx] = min(abs(T - target_week));
        C_nodes = zeros(Nstr,1);
        for j = 1:Nstr
            C_nodes(j) = X(tidx, idx_str_fn(sc,j));
        end
        plot(node_positions, C_nodes*1e6, '-o', ...
            'Color', colors_pe(pidx,:), 'LineWidth',2, ...
            'MarkerFaceColor', colors_pe(pidx,:), 'MarkerSize',7, ...
            'DisplayName', pe_labels{pidx});
    end
    hold off;
    xlabel('Stromal Position','FontSize',10);
    if sc==1; ylabel(['IgG (\muM) at wk ' num2str(target_week)],'FontSize',10); end
    title(subclass_names{sc},'FontSize',11,'FontWeight','bold');
    if sc==4; legend('Location','best','FontSize',7); end
    xticks(node_positions);
    xticklabels({'STB','','','','EC'});
    grid on; box on;
end
sgtitle(['Stromal IgG Profiles at Week ' num2str(target_week)], 'FontSize',13);

%% =========================================================================
%% PRINT SUMMARY TABLE
%% =========================================================================
fprintf('\n============================================================\n');
fprintf('  F:M Ratios at Delivery\n');
fprintf('============================================================\n');
fprintf('%-22s %8s %8s %8s %8s\n','Condition','IgG1','IgG2','IgG3','IgG4');
fprintf('------------------------------------------------------------\n');
fprintf('%-22s %8.3f %8.3f %8.3f %8.3f\n','Healthy', FM_healthy);
fprintf('%-22s %8.3f %8.3f %8.3f %8.3f\n','Endometriosis', FM_endo);
for pidx = 1:length(pe_patterns)
    fprintf('%-22s %8.3f %8.3f %8.3f %8.3f\n', ...
        pe_labels{pidx}, results_pe(pidx).FM_ratio);
end
fprintf('============================================================\n');

fprintf('\n  Percent reduction vs Healthy baseline\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-22s %8s %8s %8s %8s\n','Condition','IgG1','IgG2','IgG3','IgG4');
fprintf('------------------------------------------------------------\n');
fprintf('%-22s %7.1f%% %7.1f%% %7.1f%% %7.1f%%\n','Endometriosis', ...
    (FM_healthy-FM_endo)./FM_healthy*100);
for pidx = 1:length(pe_patterns)
    fm = results_pe(pidx).FM_ratio;
    pct = (FM_healthy - fm) ./ FM_healthy * 100;
    fprintf('%-22s %7.1f%% %7.1f%% %7.1f%% %7.1f%%\n', pe_labels{pidx}, pct);
end
fprintf('============================================================\n\n');
