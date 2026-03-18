% %% Initialize model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y0 = [p.igg1,p.igg2,p.igg3,p.igg4,zeros(1,31)];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Define color palette %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ColorOrder = [0.87, 0.443, 0.3569;
%               0.706, 0.87, 0.286;
%               0.302, 0.851, 1;
%               0.251, 0, 1];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% Solve and plot fetal IgG levels over time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solution = ode15s(@(t,x) diggdt_transport_5nodes(t,x,p), tspan, x0);
% % solution= ode15s(@(t,x) diggdt(t,x,p), tspan, y0);
% 
% t = solution.x(:);
% y = solution.y.';   % rows = time points, cols = states
% 
% %% Fetal:maternal ratio at 40 weeks gestation
% % fetal states in this model:
% % x(27) = fetal IgG1
% % x(28) = fetal IgG2
% % x(29) = fetal IgG3
% % x(30) = fetal IgG4
% 
% % find time point closest to 40 weeks
% [~, idx40] = min(abs(t - 40));
% 
% % fetal concentrations at 40 weeks
% fetal_40 = [y(idx40,27), y(idx40,28), y(idx40,29), y(idx40,30)];
% 
% % maternal concentrations (constant parameters)
% maternal = [p.igg1, p.igg2, p.igg3, p.igg4];
% 
% % fetal:maternal ratio
% FM_ratio = fetal_40 ./ maternal;
% 
% figure;
% b = bar(FM_ratio);
% set(gca, 'XTickLabel', {'IgG1','IgG2','IgG3','IgG4'});
% xlabel('IgG subclass');
% ylabel('Fetal:Maternal ratio at 40 weeks');
% title('F:M Ratio for Each IgG Subclass at 40 Weeks Gestation');
% grid on;
% ylim([0, max(FM_ratio)*1.2]);
% 
% for k = 1:length(FM_ratio)
%     text(k, FM_ratio(k) + 0.02*max(FM_ratio), sprintf('%.2f', FM_ratio(k)), ...
%         'HorizontalAlignment', 'center', 'FontSize', 10);
% end
% 
% % %% Optional: display values in Command Window
% disp('Fetal concentrations at 40 weeks:')
% disp(array2table(fetal_40, ...
%     'VariableNames', {'IgG1','IgG2','IgG3','IgG4'}))
% 
% disp('Fetal:Maternal ratios at 40 weeks:')
% disp(array2table(FM_ratio, ...
%     'VariableNames', {'IgG1','IgG2','IgG3','IgG4'}))
% 
% 
% 
% % figure
% 
% for k = 1:length(t)
% 
%     IgG1_nodes = [y(k,14), y(k,36), y(k,40), y(k,44), y(k,48)];
% 
%     plot(xpos,IgG1_nodes,'-o','LineWidth',2)
% 
%     xlabel('Position across stroma')
%     ylabel('IgG1 concentration')
% 
%     title(sprintf('IgG Transport Through Stroma (Week %.1f)',t(k)))
% 
%     ylim([0 max(y(:,27))*1.2])
%     drawnow
% end








% %% Spatial positions across stroma
% N = p.Nstr;
% xpos = linspace(0,p.Lstr,N);
% 
% %% Extract IgG1 stromal nodes from solution
% % node1 = x(14)
% % node2 = x(36)
% % node3 = x(40)
% % node4 = x(44)
% % node5 = x(48)
% 
% IgG1_nodes = [y(:,14), y(:,36), y(:,40), y(:,44), y(:,48)];
% 
% %% Create mesh
% [X,T] = meshgrid(xpos,t);
% 
% figure
% surf(X,T,IgG1_nodes,'EdgeColor','none')
% colormap(parula)
% colorbar
% 
% xlabel('Position across stroma')
% ylabel('Gestational age (weeks)')
% zlabel('IgG1 concentration')
% 
% title('3D Transport of IgG1 Across Placental Stroma')
% view(45,30)
% t = solution.x;
% 
% IgG = [
%     solution.y(27,:);
%     solution.y(28,:);
%     solution.y(29,:);
%     solution.y(30,:)
% ];


% figure(1)
% subplot(1,4,1)
% plot(solution.x,solution.y(27,:),'linewidth',2,'color',ColorOrder(1,:))
% title('IgG1'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG1 (M)')
% subplot(1,4,2)
% plot(solution.x,solution.y(28,:),'linewidth',2,'color',ColorOrder(1,:))
% title('IgG2'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG2 (M)')
% subplot(1,4,3)
% plot(solution.x,solution.y(29,:),'linewidth',2,'color',ColorOrder(1,:))
% title('IgG3'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG3 (M)')
% subplot(1,4,4)
% plot(solution.x,solution.y(30,:),'linewidth',2,'color',ColorOrder(1,:))
% title('IgG4'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG4 (M)')


% Heatmap of antibody levels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% 
% imagesc(t,1:4,IgG)
% 
% set(gca,'YTick',1:4,'YTickLabel',{'IgG1','IgG2','IgG3','IgG4'})
% xlabel('Gestational Age (weeks)')
% ylabel('Antibody Subclass')
% 
% title('Fetal IgG Concentration Heatmap')
% 
% colorbar
% colormap(parula)





% %% Initialize model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y0 = [p.igg1,p.igg2,p.igg3,p.igg4,zeros(1,31)];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Define color palette for plotting downstream figure %%%%%%%%%%%%%%%%%%%%
% ColorOrder = [0.87, 0.443, 0.3569; 0.706, 0.87, 0.286; 0.302, 0.851, 1; 0.251, 0, 1];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Solve and plot fetal IgG levels over time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solution = ode15s(@(t,x) diggdt_transport_5nodes(t,x,p), tspan, x0);
% % solution= ode15s(@(t,x) diggdt(t,x,p), tspan, y0);
% 
% 
% figure(1)
% subplot(1,4,1)
% plot(solution.x,solution.y(27,:),'linewidth',2,'color',ColorOrder(1,:))
% title('IgG1'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG1 (M)')
% subplot(1,4,2)
% plot(solution.x,solution.y(28,:),'linewidth',2,'color',ColorOrder(1,:))
% title('IgG2'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG2 (M)')
% subplot(1,4,3)
% plot(solution.x,solution.y(29,:),'linewidth',2,'color',ColorOrder(1,:))
% title('IgG3'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG3 (M)')
% subplot(1,4,4)
% plot(solution.x,solution.y(30,:),'linewidth',2,'color',ColorOrder(1,:))
% title('IgG4'); xlabel('Gestational Age (weeks)'); ylabel('Fetal IgG4 (M)')



%  Comparison bar chart %%%

% % ============================================================
% Compare F:M ratios at 40 weeks:
%   1) Original model
%   2) New transport model
%   3) Clinical data from paper (Fig. 1C)
% =============================================================

clear; clc; close all;

%% ---------------- ORIGINAL MODEL ----------------
parameters_Erdogan;   % should create p, x0, tspan

sol_orig = ode15s(@(t,x) diggdt(t,x,p), tspan, x0);

t_orig = sol_orig.x(:);
y_orig = sol_orig.y.';   % rows = time points

% index closest to 40 weeks
[~, idx40_orig] = min(abs(t_orig - 40));

% fetal states in original model
% x(27)=IgG1_fetus, x(28)=IgG2_fetus, x(29)=IgG3_fetus, x(30)=IgG4_fetus
fetal_orig_40 = [y_orig(idx40_orig,27), ...
                 y_orig(idx40_orig,28), ...
                 y_orig(idx40_orig,29), ...
                 y_orig(idx40_orig,30)];

maternal_orig = [p.igg1, p.igg2, p.igg3, p.igg4];
FM_orig = fetal_orig_40 ./ maternal_orig;


%% ---------------- NEW TRANSPORT MODEL ----------------
parameters_Erdogan_transport_5nodes;   % should create p, x0, tspan

sol_new = ode15s(@(t,x) diggdt_transport_5nodes(t,x,p), tspan, x0);

t_new = sol_new.x(:);
y_new = sol_new.y.';

[~, idx40_new] = min(abs(t_new - 40));

% fetal states remain 27:30 in the transport version
fetal_new_40 = [y_new(idx40_new,27), ...
                y_new(idx40_new,28), ...
                y_new(idx40_new,29), ...
                y_new(idx40_new,30)];

maternal_new = [p.igg1, p.igg2, p.igg3, p.igg4];
FM_new = fetal_new_40 ./ maternal_new;


%% ---------------- CLINICAL DATA FROM PAPER ----------------
clinical_ratio = [1.4601, 1.0453, 1.2274, 1.24];

%% ---------------- GROUPED BAR GRAPH ----------------
labels = {'IgG1','IgG2','IgG3','IgG4'};

compare_mat = [FM_orig(:), FM_new(:), clinical_ratio(:)];

figure;
b = bar(compare_mat, 'grouped');
set(gca, 'XTickLabel', labels, 'FontSize', 11);
xlabel('IgG subclass');
ylabel('Fetal:Maternal ratio at 40 weeks');
title('Comparison of IgG F:M Ratios at 40 Weeks');
legend({'Original model','New transport model','Clinical data'}, ...
       'Location','northeast');
grid on;

% optional horizontal line at F:M = 1
hold on;
yline(1, '--k', 'LineWidth', 1);
hold off;


%% ---------------- OPTIONAL: ADD VALUE LABELS ----------------
for k = 1:size(compare_mat,2)
    xend = b(k).XEndPoints;
    yend = b(k).YEndPoints;
    labels_txt = strings(size(yend));

    for j = 1:numel(yend)
        if isnan(yend(j))
            labels_txt(j) = "";
        else
            labels_txt(j) = sprintf('%.2f', yend(j));
        end
    end

    text(xend, yend, labels_txt, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 9);
end


%% ---------------- PRINT VALUES TO COMMAND WINDOW ----------------
disp('Original model F:M ratios at 40 weeks:')
disp(array2table(FM_orig, 'VariableNames', labels))

disp('New transport model F:M ratios at 40 weeks:')
disp(array2table(FM_new, 'VariableNames', labels))

disp('Clinical data F:M ratios at 40 weeks:')
disp(array2table(clinical_ratio, 'VariableNames', labels))



figure;

% IgG1
subplot(2,2,1)
plot(t_orig, y_orig(:,27), 'LineWidth', 2); hold on
plot(t_new,  y_new(:,27), '--', 'LineWidth', 2);
title('IgG1')
xlabel('Weeks'); ylabel('Conc')
legend('Dolatshahi & Wessel','Improved'); grid on

% IgG2
subplot(2,2,2)
plot(t_orig, y_orig(:,28), 'LineWidth', 2); hold on
plot(t_new,  y_new(:,28), '--', 'LineWidth', 2);
title('IgG2')
xlabel('Weeks'); ylabel('Conc')
legend('Dolatshahi & Wessel','Improved'); grid on

% IgG3
subplot(2,2,3)
plot(t_orig, y_orig(:,29), 'LineWidth', 2); hold on
plot(t_new,  y_new(:,29), '--', 'LineWidth', 2);
title('IgG3')
xlabel('Weeks'); ylabel('Conc')
legend('Dolatshahi & Wessel','Improved'); grid on

% IgG4
subplot(2,2,4)
plot(t_orig, y_orig(:,30), 'LineWidth', 2); hold on
plot(t_new,  y_new(:,30), '--', 'LineWidth', 2);
title('IgG4')
xlabel('Weeks'); ylabel('Conc')
legend('Dolatshahi & Wessel','Improved'); grid on

sgtitle('Fetal IgG Subclasses Over Gestation: Original vs Improved Model')