%% parameters_endometriosis.m
%% Endometriosis disease model - focused parameter set
%% Dolatshahi Lab
%%
%% PARAMETERS CHANGED FROM BASELINE (5 total):
%%
%%   1. p.fcrn    - STB FcRn expression     (-30%)
%%                  #1 driver of total IgG transfer (highest VIP score)
%%                  Endometriosis impairs trophoblast differentiation,
%%                  directly reducing FcRn expression in STBs
%%
%%   2. p.fcgr2b  - EC FcgRIIb expression   (-30%)
%%                  #1 driver of subclass selectivity (highest VIP score)
%%                  Chronic inflammation alters EC phenotype and
%%                  downregulates FcgRIIb
%%
%%   3. p.k_t     - Transcytosis rate        (-25%)
%%                  #1 kinetic parameter by VIP score
%%                  Inflammatory cytokine environment (TNF-a, IL-6)
%%                  disrupts vesicular trafficking in both STBs and ECs
%%
%%   4. p.Dstr    - Stromal diffusivity      (patchy reduction)
%%                  Fibrotic stromal remodeling unique to endometriosis
%%                  Modeled as a vector to capture spatially heterogeneous
%%                  fibrotic patches between stromal nodes
%%
%%   5. p.t_end   - Gestational length       (40 -> 37 weeks)
%%                  Endometriosis increases preterm birth risk
%%                  Truncates the transport window, especially the steep
%%                  third-trimester accumulation phase
%%
%% ALL OTHER PARAMETERS ARE IDENTICAL TO BASELINE

%% =========================================================================
%% START FROM BASELINE - load all original parameters first
%% =========================================================================
parameters_Erdogan_transport_5nodes;   % loads all baseline values into p

%% =========================================================================
%% DISEASE-SPECIFIC CHANGES
%% =========================================================================

%% CHANGE 1: STB FcRn expression (-30%)
%% Rationale: Highest VIP score for total IgG transfer in sensitivity
%% analysis. Endometriosis impairs trophoblast invasion and differentiation,
%% directly reducing FcRn (FCGRT) transcription in STBs.
p.fcrn = mean([3.33e-5, 4.94e-5]) * 0.70;   % -30% vs baseline

%% CHANGE 2: EC FcgRIIb expression (-30%)
%% Rationale: Highest VIP score for subclass selectivity. Chronic
%% inflammation from endometriosis alters endothelial cell phenotype.
%% Inflammatory cytokines (TNF-a, IL-6) downregulate FCGR2B expression.
p.fcgr2b = mean([2.73e-5, 3.33e-5]) * 0.70;   % -30% vs baseline

%% CHANGE 3: Transcytosis rate (-25%)
%% Rationale: Highest-ranked kinetic parameter by VIP score. Acts in an
%% FcR-dependent manner to amplify selectivity effects. Disrupted
%% vesicular trafficking from the inflammatory microenvironment reduces
%% the rate at which receptor-bound IgG is shuttled across cell layers.
p.k_t = mean([0.0569, 0.0743]) * 0.75;   % -25% vs baseline

%% CHANGE 4: Stromal diffusivity - patchy fibrosis
%% Rationale: Unique pathological feature of endometriosis. Fibrotic
%% stromal remodeling creates spatially heterogeneous transport resistance.
%% Using a vector rather than a scalar captures the patchy nature of
%% fibrotic lesions between stromal nodes (STB side -> EC side).
%%   Node 1->2: healthy tissue (near STB)
%%   Node 2->3: moderate fibrosis (central lesion)
%%   Node 3->4: severe fibrosis (core of lesion)
%%   Node 4->5: partial recovery (near EC)
p.Dstr = [1.0, 0.30, 0.10, 0.40];   % patchy fibrosis pattern

%% CHANGE 5: Gestational length (37 weeks)
%% Rationale: Endometriosis increases preterm birth risk approximately
%% 2-fold. Delivering at 37 weeks truncates the steep third-trimester
%% IgG accumulation phase, reducing total fetal antibody acquisition.
p.t_end = 37;
tspan       = (10:1:p.t_end);
time_points = linspace(10, p.t_end-1, 100);

%% =========================================================================
%% UPDATE DEPENDENT PARAMETERS
%% Re-fit the FcR expression curves with the new peak values
%% (must be done after changing p.fcrn and p.fcgr2b)
%% =========================================================================
days_in_rat   = [1 15 21]' * 40/21;
fcrn_vals     = [0, p.fcrn/2,    p.fcrn]';
fcgr2b_vals   = [0, p.fcgr2b/2,  p.fcgr2b]';
fcrn_ec_vals  = [0, p.fcrn_ec/2, p.fcrn_ec]';   % fcrn_ec unchanged

p.fcrn_STB_curve  = fit(days_in_rat, fcrn_vals,   'poly2');
p.fcgr2b_EC_curve = fit(days_in_rat, fcgr2b_vals, 'poly2');
p.fcrn_EC_curve   = fit(days_in_rat, fcrn_ec_vals,'poly2');

%% =========================================================================
%% SUMMARY PRINTOUT
%% =========================================================================
fprintf('\n--- Endometriosis Model Parameters ---\n');
fprintf('  p.fcrn    = %.3e M  (baseline * 0.70, -30%%)\n', p.fcrn);
fprintf('  p.fcgr2b  = %.3e M  (baseline * 0.70, -30%%)\n', p.fcgr2b);
fprintf('  p.k_t     = %.4f L/wk  (baseline * 0.75, -25%%)\n', p.k_t);
fprintf('  p.Dstr    = [%.2f, %.2f, %.2f, %.2f]  (patchy fibrosis)\n', p.Dstr);
fprintf('  p.t_end   = %d weeks  (preterm risk)\n', p.t_end);
fprintf('  All other parameters: baseline values\n\n');
