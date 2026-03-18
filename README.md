# BENG227-Group-6-Final-Code
These files have been altered from Remziye E Wessel, Sepideh Dolatshahi. Quantitative mechanistic model reveals key determinants of placental IgG transfer and informs prenatal immunization strategies. PLoS Comput Biol 2023;19:. https://doi.org/10.1371/journal.pcbi.1011109. Github repository.

diggdt_transport_improved.m: Defines the system of ODEs and implementation of method of lines for the PDE that solves the stromal placental transport across a spatial domain.

parameters_improved.m: Initializes all model parameters for all non-disease modeling including biological constants, binding kinetics, transport rates, and diffusion properties. Parameter p.Dstr here can be altered to represent different diffusion coefficients for the spatial gradient modeling.

diggdt_driver_spatialgradients.m: Runs the transport model simulation and generates plots of original fetal IgG subclass concentrations, heatmap of spatial concentration gradients, and IgG subclass spatial concentration gradients.

diggdt_driver_FM.m: Driver function that computes model dynamics for fetal and maternal system, creating a fetal:maternal plot of IgG transport.

parameters_endometriosis.m: Defines parameters specific to modeling endometriosis conditions.

parameters_preeclampsia.m: Defines parameters specific to modeling preeclampsia conditions.

run_preeclampsia_comparison.m: Script that runs simulations to compare healthy vs. diseased model outputs.
