%% parameters_preeclampsia.m
%% Preeclampsia disease model - focused parameter set
%% Dolatshahi Lab

%% Only set default if caller did not already specify it

%%
%% PARAMETERS CHANGED FROM BASELINE (5 total):
%%
%%   1. p.fcrn    - STB FcRn expression     (-55%)
%%                  #1 driver of total IgG transfer (highest VIP score)
%%                  PE is the most well-documented cause of placental FcRn
%%                  loss. Hypoxia directly downregulates FCGRT transcription
%%                  in STBs. More severe than endometriosis.
%%                  (Atyeo et al 2021; Borghi et al 2020)
%%
%%   2. p.fcgr2b  - EC FcgRIIb expression   (-50%)
%%                  #1 driver of subclass selectivity (highest VIP score)
%%                  Elevated sFlt-1 and reduced PlGF in PE cause systemic
%%                  endothelial dysfunction, severely downregulating
%%                  FcgRIIb on fetal ECs. More severe than endometriosis.
%%
%%   3. p.k_t     - Transcytosis rate        (-45%)
%%                  #1 kinetic parameter by VIP score
%%                  Hypoxia-reperfusion injury in PE severely damages
%%                  vesicular trafficking machinery in STBs and ECs.
%%                  More severe impairment than endometriosis.
%%
%%   4. p.ustr    - Stromal advection velocity  (0 -> -0.20)
%%                  PE-specific parameter: impaired spiral artery remodeling
%%                  reduces uteroplacental perfusion, introducing an opposing
%%                  convective flow that counteracts STB->EC diffusion.
%%                  Not present in endometriosis or healthy placenta.
%%
%%   5. p.t_end   - Gestational length       (40 -> 34 weeks)
%%                  PE typically requires delivery at 34-37 weeks.
%%                  Early-onset severe PE often delivers at 34 weeks,
%%                  truncating more of the transport window than endometriosis.
%%
%% ALL OTHER PARAMETERS ARE IDENTICAL TO BASELINE

%% =========================================================================
%% START FROM BASELINE - load all original parameters first
%% =========================================================================
parameters_improved;   % loads all baseline values into p

%% =========================================================================
%% DISEASE-SPECIFIC CHANGES
%% =========================================================================


%% CHANGE 1: STB FcRn expression (-55%)
%% Rationale: Highest VIP score for total IgG transfer. PE causes more
%% severe FcRn loss than endometriosis because hypoxia directly suppresses
%% FCGRT gene expression in trophoblasts. Multiple clinical studies
%% document significantly lower FcRn in preeclamptic placentas compared
%% to both healthy controls and other pregnancy complications.
p.fcrn = mean([3.33e-5, 4.94e-5]) * 0.45;   % -55% vs baseline

%% CHANGE 2: EC FcgRIIb expression (-50%)
%% Rationale: Highest VIP score for subclass selectivity. The sFlt-1/PlGF
%% imbalance characteristic of PE causes systemic endothelial dysfunction,
%% more severely reducing FcgRIIb than the local inflammatory effects of
%% endometriosis. This also disproportionately impairs IgG3/IgG4 transfer
%% since FcgRIIb has highest affinity for these subclasses.
p.fcgr2b = mean([2.73e-5, 3.33e-5]) * 0.50;   % -50% vs baseline

%% CHANGE 3: Transcytosis rate (-45%)
%% Rationale: Highest-ranked kinetic parameter by VIP score. Hypoxia-
%% reperfusion injury in PE causes oxidative damage to trophoblast and EC
%% cytoskeletal machinery, more severely disrupting vesicular trafficking
%% than the inflammatory environment of endometriosis alone.
p.k_t = mean([0.0569, 0.0743]) * 0.55;   % -45% vs baseline

%% CHANGE 4: Stromal advection velocity (-0.20, opposing flow)
%% Rationale: PE-specific parameter not present in endometriosis or
%% baseline. Impaired spiral artery remodeling in PE reduces
%% uteroplacental blood flow, introducing an opposing convective term
%% that works against the STB->EC diffusion direction. This reduces
%% the effective flux of IgG through the stroma beyond what diffusivity
%% reduction alone would predict. The negative sign indicates flow
%% opposing the STB-to-EC transport direction.
%% Note: p.Dstr kept as uniform scalar (ischemic edema rather than
%% patchy fibrosis - diffusivity uniformly reduced across stroma)
%% Only set default if caller did not already specify it
if ~isfield(p, 'pe_pattern')
    p.pe_pattern = 'severe_pe';
end

switch p.pe_pattern

    case 'mild_pe'
        p.Dstr = p.k_up * [0.55, 0.55, 0.55, 0.55];
        p.ustr = -0.05;

    case 'severe_pe'
        p.Dstr = p.k_up * [0.10, 0.25, 0.40, 0.30];
        p.ustr = -0.20;

    case 'early_onset'
        p.Dstr = p.k_up * [0.05, 0.15, 0.20, 0.25];
        p.ustr = -0.35;

    case 'late_onset'
        p.Dstr = p.k_up * [0.70, 0.55, 0.35, 0.20];
        p.ustr = -0.08;

    case 'custom'
        p.Dstr = p.k_up * [0.30, 0.30, 0.30, 0.30];
        p.ustr = -0.10;

    otherwise
        error('Unknown pe_pattern: %s', p.pe_pattern);
end
%% CHANGE 5: Gestational length (34 weeks)
%% Rationale: PE requires earlier delivery than endometriosis. Severe
%% and early-onset PE typically mandates delivery at 34 weeks or earlier
%% to protect maternal health. This cuts off 6 weeks of the third-
%% trimester IgG accumulation phase versus the healthy baseline,
%% and 3 weeks more than the endometriosis model.
p.t_end = 34;
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
fprintf('\n--- Preeclampsia Model Parameters ---\n');
fprintf('  p.fcrn    = %.3e M  (baseline * 0.45, -55%%)\n', p.fcrn);
fprintf('  p.fcgr2b  = %.3e M  (baseline * 0.50, -50%%)\n', p.fcgr2b);
fprintf('  p.k_t     = %.4f L/wk  (baseline * 0.55, -45%%)\n', p.k_t);
fprintf('  p.ustr    = %.2f  (opposing advection, PE-specific)\n', p.ustr);
fprintf('  p.Dstr    = %.2f  (uniform edema reduction)\n', p.Dstr);
fprintf('  p.t_end   = %d weeks  (earlier delivery)\n', p.t_end);
fprintf('  All other parameters: baseline values\n\n');
