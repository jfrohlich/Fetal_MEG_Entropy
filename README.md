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
## The prinicipal scripts have a lot of dependencies! I've tried my best to list and describe as many as possible below ... 

### Scripts (These include helper functions, my custom mods of native MATLAB functions, external code, etc.) 
### Note that many of these are publicly available external scripts! See code comments for proper attribution. 

- **betacdf.m:** Octave function, CDF of the beta distribution. **external code**
- **betainv.m:** Octave function, quantile function of the beta distribution. **external code**
- **betapdf.m:** Octave function, PDF of the beta distribution. **external code**
- **cohens_d.m:** Compute the Cohen's d effect size.
- **customcolormap.m:** Creates alternative colormaps (downloaded from [MathWorks File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/69470-custom-colormap)). **external code**
- **customcolormap_preset.m:** Works with the above function. **external code**
- **CI_highlight.m:** Wrapper function for easier use of the `patch()` function to plot confidence interval highlights.
- **CTWEntropRate.m:** Pedro Mediano's code for CTW entropy. **external code**
- **data_range.m:** Computes the improvement (smaller range across recordings) after normalizing fetal data.
- **filter_effects.m:** Shows effects of different filter lowpass filter settings on neonatal data. 
- **GetAverageERFNeonatal.m:** Generates ERFs for the newborn ERF template.
- **GetAverageERFFetal.m:** Similar to the above, but for fetuses (specifically, older fetuses).
- **LZ76.m:** Pedro Mediano's Lempel-Ziv code. **external code**
- **makefighandsome.m:** Makes aesthetic changes to figures.
- **makefigpretty.m:** Another one of my scripts to make aesthetic changes to figures.
- **myICC.m:** Modification of external code to compute the intraclass correlation (ICC).  **external code**
- **mycolorbar.m:** Nicer-looking colorbar.
- **myfigure.m:** Bigger window for the figure.
- **myfigure2.m:** Even bigger window for the figure.
- **mylsline.m:** Nicer-looking least squares line.
- **mypcolor.m:** Improvement over MATLAB's pcolor.
- **myviolin.m:** Generates violin plots (modified from code on [MathWorks File Exchange](https://de.mathworks.com/matlabcentral/fileexchange/45134-violin-plot)).  **external code**
- **NewbornSubjectIDDecoder.m:** Fixes errors in decoding subject IDs.
- **PermEn.m:** Computes permutation entropy (modified from code downloaded from [MathWorks File Exchange](https://de.mathworks.com/matlabcentral/fileexchange/37289-permutation-entropy)).  **external code**
- **prop_test.m:** Computes the chi-squared test (downloaded from [Mathworks File Exchange](https://de.mathworks.com/matlabcentral/fileexchange/45966-compare-two-proportions-chi-square)).  **external code**
- **ro_freq_wavelet.m:** Joerg Hipp's code to get power spectrum (averaged over time) with Morlet wavelets.  **external code**
- **ro_freq_wavelet_TFT.m:** Computes the Morlet wavelet transform, written by Joerg Hipp with slight modifications. **external code**
- **ro_mse.m:** Joerg Hipp's multiscale entropy code. **external code**
- **sortref.m:** Sorts one list with respect to another.
- **surrogate_JF.m:** Modification of *external code* to compute surrogate data.
- **stochastic_test_JF.m:** Tests whether a signal is stochastic or deterministic, modified from *external code* by Toker et al. 2020.
- **SubjectTable.m:** Generates a table of all subjects.
- **table2latex.m:** External code for writing tables in LaTeX format. 
- **timeresentropy_preissl.m:** My code to compute time-resolved entropy (with CTW) for this project. 
- **TFCE.m:** Function for threshold-free cluster enhancement statistics.

### Spreadsheets

- **FetalModel.csv:** Test statistics from linear mixed models predicting fetal entropy measures.
- **FetalModelDynamics.csv:** Test statistics from linear mixed model predicting fetal cortical dynamics (stochastic versus deterministic).
- **FetalSurrogate.csv:** Test statistics predicting surrogacy for each fetal entropy measure. 
- **NeonatalModel.csv:** Test statistics from linear mixed models predicting neonatal entropy measures.
- **NeonatalModelDynamics.csv:** Test statistics from linear mixed model predicting neonatal cortical dynamics (stochastic versus deterministic).
- **NeonatalSurrogate.csv:** Test statistics predicting surrogacy for each neonatal entropy measure. 
- **NeonatalTable.csv:** MEG entropy values from neonatal recordings.
- **NeonatalEntropyPower.csv:** This is an updated version of the table from the Zenodo repository with added values for entropy, spectral power, dynamics, etc.
- **fMEGEntropyDecomp.csv:** Test statics from models testing for differences between amplitude and non-amplitude components in fetal entropy measures.
- **neonate_decomp.csv:** Test statics from models testing for differences between amplitude and non-amplitude components in neonatal entropy measures.
  
### .MAT Files

- **Alltable.mat:** Huge table with most relevant variables for all subjects
- **CTWDecompfMEG2023_Nsur=250.mat:** Output of fetal CTW entropy decomposition
- **CTWDecompneonate2023_Nsur=250.mat:** Output of neonatal CTW entropy decomposition
- **LZCDecompfMEG2023_Nsur=250.mat:** Output of fetal LZC entropy decomposition
- **LZCDecompneonate2023_Nsur=250.mat:** Output of neonatal LZC entropy decomposition
- **MSEDecompfMEG2023_Nsur=250.mat:** Output of fetal mMSE entropy decomposition
- **MSEDecompneonate2023_Nsur=250.mat:** Output of neonatal mMSE entropy decomposition
- **Newborn_ERFs.mat:** Contains a "template" of the average neonatal ERF response -- isn't strictly necessary anymore but some functions will still ask for it
- **PermEnDecompfMEG2023_Nsur=250.mat:** Output of fetal PermEn entropy decomposition (both 32 and 64 ms values for tau)
- **PermEnDecompneonate2023_Nsur=250.mat:** Output of neonatal PermEn entropy decomposition (both 32 and 64 ms values for tau)
- **grand_average_ERFs.mat:**  Contains a "template" of the average fetal ERF response -- isn't strictly necessary anymore but some functions will still ask for it
- **mSampEnDecompfMEG2023_Nsur=250.mat:** Output of fetal mSampEn entropy decomposition
- **mSampEnDecompneonate2023_Nsur=250.mat:** Output of neonatal mSampEn entropy decomposition
- **newbornIDs.mat:** Needed at some point to understand neonatal subject IDs


