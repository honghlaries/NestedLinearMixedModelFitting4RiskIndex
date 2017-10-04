# Brief Introduction
These codes were provided as examples for manuscript entitled "Risk assessment for Artificial Topography on Element Heterogeneity: case study at Jiangsu, China"


# Package Dependencies
The following R packages are required:

*dplyr*, *tidyr*, *ggplot2*, *lme4*, *merTools*, *lmerTest*, *lsmeans*, *multcompView*, *car*, *lattice*, *psych*, *MASS*, *gridExtra*, *maptools*, *grid*, *ggmap*, *maptools*

All of these packages are available on CRAN.

I have regrouped some codes in the *merTools* package to fitting the linear mixed models in parallel. Now the model selection (based on AIC), LMM fitting and ploting (fixed effects and random effects) can be finished automatically and some log files are generated for review. 

# Summary of the scripts
## Igeo.R (for geo-accumulation index) and EF.R (for enrichment factor)
These two scripts are used for calculation and visualization of geo-accumulation index and enrichment factor.

### datareadln()
Read the data files (.csv) in the \\data.

### meanseCal()
Calculate and output the means and se as .csv file.

### shapiroTest()
Perform Shapiro Test on the selected data. 

### modelCompare()
Compare different models based on AIC, output the best model in the model list. Log files (optional, output as default) are generated to tracing the model comparing process.

### plotFEsim2()
Regrouped from *merTools* package. Add new features to emphasis the targeted factors. Changes made from *merTools* packages could be found in the source fils.

### plotFEsim2facet()
Similar to the plotFEsim2() but output a grouped faceted plot. Supproting by the facet_wrap() commend of *ggplot2* package. Changes made from *merTools* packages could be found in the source fils.

### multiElementMod()
This is the main function of these two script files. It finishes tasks as fellow: 

1. Perform a Shapiro Test on the dataset using shapiroTest();
2. Select a best fit models based on AIC using modelCompare();
3. Run ANOVA and post hoc analysis (Tukey test) on the selected model;
4. Extract the fixed effects and random effects of every target and output the single-target plot for fixed effect and random effect;
5. Extract the residuals of LMMs and perform Shapiro Tests;
6. Output the log files, including Shapiro Test results, Model choices, Fixed effects, Random effects, ANOVA results, and Tukey Test results.

### denPlot()
Prepare a density plot for selected target.

## Pca.R
Calculate and output the PCA result, with selected factors emphasized by a different color (red as default). 

### datareadln()
Read the data files (.csv) in the \\data.

### pcaLoadingCal()
Calculate the PCA loading and output as a .csv file. Scree Plot generated as review file.

### pcaLoadingPlot()
visualize the PCA loading, with target farotrs (emphtag) emphasized by a different color (red as default).

### pairPlot()
Provide a x-y plot in different groups, perform a linear regression to confirm the relation between targeted factors.  


## samplemap.R
Script to generate map components representing the station location.



