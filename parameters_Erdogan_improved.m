%% KD binding parameters for IgG subclasses binding to FcRn and FcgRIIb %%%
%% k_d_fcrn and k_d_fcgr2b are given as a vector with KD parameters corresponding to
%% [IgG1, IgG2, IgG3, IgG4].  Parameters derived from Bruhns, P. et al (2009).
p.k_d_fcrn = 1./[8e7 5e7 3e7 2e7]; % Subclass dissociation constants for FcRn (M)
p.k_d_fcgr2b = 1./[1e5,2e4,2e5,2e5]; % Subclass dissociation constants for FcgRIIb (M)
p.mass_igg = 1.5e5; %molar mass of IgG (Da)
p.igg1 = 5.67/p.mass_igg; %maternal IgG1 concentration (M)
p.igg2 = 2.71/p.mass_igg; %maternal IgG2 concentration (M)
p.igg3 = 0.4/(p.mass_igg+2e4); % maternal IgG3 concentration (M) - note the higher molar mass.
p.igg4 = 0.34/p.mass_igg; % maternal IgG4 concentration (M)
%% Initialize parameters - Optimized parameter ranges coming from CaliPro adj for  k_on/k_off
   p.k_up = mean([0.0872,0.0982]); %transcytosis rate (L/week)
   p.k_t = mean([0.0569,0.0743]); %IgG uptake rate into cells (L/week)
   p.fcrn = mean([3.33e-5,4.94e-5]); %FcRn expression level on STBs (M)
   p.fcgr2b = mean([2.73e-5,3.33e-5]); %FcgRIIb expression level on ECs (M)
   p.v_m = mean([5.552,5.619]); %volume of plasma in mom (L)
   p.v_f = mean([0.254,0.271]); %volume of plasma in fetus (L) (Prior estimation based on - Lin and Tran, 2013)
   p.v_str = mean([0.165,0.179]); %total volume of endothelial cell endosomes (L)
% --- Stroma transport (method-of-lines) parameters ---
p.Nstr  = 5;        % number of stromal spatial nodes
p.Lstr  = 1.0;      % stromal thickness (model units)
p.Dstr  = 1;     % effective stromal diffusivity (model units) able to be changed
p.ustr  = 0.0;      % advection velocity where 0 = diffusion-only
p.dxstr = p.Lstr/(p.Nstr-1);
   p.v_stb = mean([0.224,0.258]); %total volume of syncytiotrophoblast endosomes (L)
   p.v_ec = mean([0.247,0.275]); %volume of villous stroma (L) (Prior estimation based on - Abdalla et al 2016)
   p.k_deg = mean([8.4733,9.0809]); %lysosomal degradation rate in STBs (L/week)
   p.fcrn_ec = mean([2.36e-6,2.86e-6]); %FcgRIIb expression level on ECs (M)
   p.d_ab = log(2)/(31/7)*p.v_f; %antibody decay rate (L/week)
%% Optional transport decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% These parameters split the original lumped transport terms into
%% compartment-to-compartment fluxes. The defaults recover the original
%% model. Set p.k_diff_str > 0 to enable passive, gradient-driven exchange
%% between villous stroma and endothelial endosomes.
p.k_m_to_stb = p.k_up; %maternal blood -> STB uptake (L/week)
p.k_stb_to_str = p.k_t; %STB complex export -> stroma (L/week)
p.k_str_to_ec = p.k_up; %stroma -> EC uptake (L/week)
p.k_ec_to_f = p.k_t; %EC export -> fetal blood (L/week)
p.k_diff_str = 0; %passive stroma <-> EC diffusion coefficient (L/week)
%% Establish time span of simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This can be adjusted to consider shorter gestational lengths
%% (i.e., premature birth).
t_end=40; % length of the simulation (weeks)
tspan=(10:1:t_end);   % time points where the output is calculated
time_points=linspace(10,t_end-1,100); % time points of interest
%% Model initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [p.igg1,p.igg2,p.igg3,p.igg4,zeros(1,47)]; % 51 states (adds 16 stromal-node states)
%% Extract kon (rate of forward reaction) and koff (rate of reverse reaction)
%% from KD values. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FcRn kinetic binding parameters.
%https://www.tandfonline.com/doi/epdf/10.1080/19420862.2015.1008353?needAccess=true&role=button
%IgG1
% p.kon1 = 5.4e6; %(1/Ms)
p.koff1 = 1.5e-1; %(1.4 1/s)
p.kon1 = p.koff1/p.k_d_fcrn(1); %(1/Ms)
p.koff2=p.koff1; p.koff3=p.koff1; p.koff4=p.koff1;
%IgG2
p.kon2 = p.koff2/p.k_d_fcrn(2); %(1/Ms)
%IgG3
p.kon3 = p.koff3/p.k_d_fcrn(3); %(1/Ms)
%IgG4
p.kon4 = p.koff4/p.k_d_fcrn(4); %(1/Ms)
%FcgRIIb kinetic binding parameters
%IgG1
p.koff1b = 1.5e-1; %(1.4 1/s)
p.kon1b = p.koff1b/p.k_d_fcgr2b(1); %(1/Ms)
p.koff2b=p.koff1b; p.koff3b=p.koff1b; p.koff4b=p.koff1b;
%IgG2
p.kon2b = p.koff2b/p.k_d_fcgr2b(2); %(1/Ms)
%IgG3
p.kon3b = p.koff3b/p.k_d_fcgr2b(3); %(1/Ms)
%IgG4
p.kon4b = p.koff4b/p.k_d_fcgr2b(4); %(1/Ms)
%% This code turns the expression level of FcRn and FcgRIIb into a dynamic
%% (parabolic) curve, which is sampled in the model to simulate increasing
%% receptor expression over time. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FcRn expression data from Wang et al (AJRI) was used to calibrate the trend
%% in human placenta.
days_in_rat = [1 15 21]'*40/21; %Scale the time data from rats in Wang et al to human gestational length.
fcgr2b = [0 p.fcgr2b/2 p.fcgr2b]'; %relative FcgRIIb expression in rat placenta
fcrn = [0 p.fcrn/2 p.fcrn]'; %relative FcRn expression in rat placenta
fcrn_ec = [0 p.fcrn_ec/2 p.fcrn_ec]';
% Curve fit a parabola to the data scaled from Wang et al.
p.fcrn_STB_curve = fit(days_in_rat,fcrn,'poly2');
p.fcgr2b_EC_curve = fit(days_in_rat,fcgr2b,'poly2');
p.fcrn_EC_curve = fit(days_in_rat,fcrn_ec,'poly2');
%% Vaccine simulation parameters (used in module 04) %%%%%%%%%%%%%%%%%%%%%%
%% (Some model parameters were extracted from White et al (2014))
p.rho = 0.96; %fraction of S-ASCs (95%) ***fixed***
p.k_asc = 0.6*7*p.v_m; %rate of ASC generation
p.k_ab = 4.82e-11*7*p.v_m; %rate of antibody secretion (scaled from IU/ml to M)
p.delta_l = log(2)/1050*7*p.v_m; %decay rate of L-ASCs (1/weeks) ***fixed***
p.delta_s = log(2)/4*7*p.v_m; %decay rate of S-ASCs (1/weeks) ***fixed***
p.delta_ab = log(2)/21*7*p.v_m; %decay rate of antibody (1/weeks) ***fixed***
p.delta_ag = 0.15*7*p.v_m;
