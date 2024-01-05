 # Fetal_MEG_Entropy
Code for the analyses of fetal and neonatal MEG in Frohlich et al. 2024 Nature Mental Health

Last updated 2024-01-05

This project uses archived data from Moser et al. 2020 (https://zenodo.org/records/4018827) and Moser et al. 2021 (https://zenodo.org/records/4541463). Please visit those archives on Zenodo for the corresponding data. 

NOTE: I have tried to include as many external dependencies here as possible (and please note that the repo license only applies to my code and not these dependencies), but if you run into problems, try also adding Fieldtrip to your path before running my code (https://www.fieldtriptoolbox.org/). 

# Principal scripts

## Entropy Computation

- **ana_MEG_fetal_entropy_frohlichetal2024.m:** Computes fetal entropy measures and runs stochasticity tests for fetal data.
- **ana_MEG_neonate_entropy_frohlichetal2024.m:** Similar to the above, but for newborn data.

## Analysis, Statistics, and Plotting

- **ana_fMEG_stats_fetal_frohlichetal2024.m:** Analysis script for fetal data, runs statistics, and generates figures.
- **ana_fMEG_stats_neonates_frohlichetal2024.m:** Similar to the above, but for newborn data.
- **ana_fetal_entropy_decomp_frohlichetal2024.m:** Runs entropy decomposition for fetal data (younger vs. older fetuses).
- **ana_neonate_entropy_decomp_frohlichetal2024.m:** Runs entropy decomposition for neonatal data (younger vs. older newborns).
- **FetalMediationTest_frohlichetal2024.R:** Runs the mediation analysis to determine if effects of sex on entropy are mediated by amplitude.
-  **Big_FDR_correction_BH_frohlichetal2024.m:** Run an false discover rate correction on all results.
  
# entropy_decomp_code_adapted_from_PMediano
## (These scripts are called by the principal scripts.) 

- **CTWDecomposition_fMEG.m:** Run the entropy decomposition on contrext-tree weighted complexity.
- **LZDecomposition_fMEG.m:** Run the entropy decomposition on Lempel-Ziv complexity.
- **PermEnDecomposition_fMEG.m:** Run the entropy decomposition on permutation entropy (with tau = 32 ms and 64 ms).
- **mMSEDecomposition_fMEG.m:** Run the entropy decomposition on modified multiscale entropy.
- **mSampEnDecomposition_fMEG.m:** Run the entropy decomposition on modified sample entropy.

# Other
## (You probably won't need this code, but it's here just in case ...)

- **GetAverageERFFetal.m:** Generate an average ERF template from fetal data (using older fetuses based on results of Moser et al. 2021).
- **GetAverageERFNeonatal.m:** Generate an average ERF template from neonatal data. 

# Dependencies
## (The prinicipal scripts have a lot of dependencies! I've tried my best to list and describe as many as possible below ... )

### Scripts (These include helper functions, my custom mods of native MATLAB functions, external code, etc.)

- **cohens_d.m:** Compute the Cohen's d effect size.
- **customcolormap.m:** Creates alternative colormaps (downloaded from MathWorks File Exchange).
- **customcolormap_preset.m:** Works with the above function.
- **CI_highlight.m:** Wrapper function for easier use of the `patch()` function to plot confidence interval highlights.
- **CTWEntropRate.m:** Pedro's code for CTW entropy.
- **data_range.m:** Computes the improvement (smaller range across recordings) after normalizing fetal data.
- **filter_effects.m:** Shows effects of different filter lowpass filter settings on neonatal data. 
- **GetAverageERFNeonatal.m:** Generates ERFs for the newborn ERF template.
- **GetAverageERFFetal.m:** Similar to the above, but for fetuses >= 35 weeks.
- **LZ76.m:** Pedro's Lempel-Ziv code.
- **makefighandsome.m:** Makes aesthetic changes to figures.
- **makefigpretty.m:** Another one of my scripts to make aesthetic changes to figures.
- **myICC.m:** Modification of external code to compute the intraclass correlation (ICC).
- **mycolorbar.m:** Nicer-looking colorbar.
- **myfigure.m:** Bigger window for the figure.
- **myfigure2.m:** Even bigger window for the figure.
- **mylsline.m:** Nicer-looking least squares line.
- **mypcolor.m:** Improvement over MATLAB's pcolor.
- **myviolin.m:** Generates violin plots (modified from code on MathWorks File Exchange).
- **NewbornSubjectIDDecoder.m:** Fixes errors in decoding subject IDs.
- **PermEn.m:** Computes permutation entropy (modified from code downloaded [here](https://de.mathworks.com/matlabcentral/fileexchange/37289-permutation-entropy)).
- **PermEnDecomp_fMEG.m:** Permutation entropy decomposition, adapted from code by Pedro Mediano.
- **prop_test.m:** Computes the chi-squared test (external code downloaded from eleswhere). 
- **ro_cluster_timefreq.m:** Needed for permutation clusters statistics.
- **ro_freq_wavelet.m:** Joerg Hipp's code to get power spectrum (averaged over time) with Morlet wavelets.
- **ro_freq_wavelet_TFT.m:** Computes the Morlet wavelet transform, written by Joerg Hipp with slight modifications.
- **ro_mse.m:** Joerg's multiscale entropy code.
- **sortref.m:** Sorts one list with respect to another.
- **surrogate_JF.m:** Modification of external code to compute surrogate data.
- **stochastic_test_JF.m:** Tests whether a signal is stochastic or deterministic, modified from Toker et al. 2020.
- **SubjectTable.m:** Generates a table of all subjects.
- **swtest.m:** Shapiro Wilk test (downloaded from elsewhere, actually unsure if this is even still used but including here to be safe).
- **timeresentropy_preissl.m:** My code to compute time-resolved entropy (with CTW) for this project. 
- **TFCE.m:** Function for threshold-free cluster enhancement statistics.

### Spreadsheets

- **HypothesisTesting.csv:** Neonatal data from Julia (public repository).
- **MeasurementQualityTable.csv:** Table with age information on newborns (shared privately).
- **NeonatalEntropyPower.csv:** Entropy and power values for neonates (fetal file too large for GitHub).
- **NeonatalTable.csv:** Entropy values from neonates with variations (ERF template   subtracted, truncation, etc.).

### .MAT Files

- **Data_traces.mat:** Fetal data.
- **grand_average_ERFs.m:** ERF templates from fetuses.
- **Newborn_ERFs.mat:** ERF templates from newborns.


