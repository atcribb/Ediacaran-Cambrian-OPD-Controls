# Ediacaran-Cambrian bioturbation did not extensively oxygenate sediments in shallow marine ecosystems

Alison T. Cribb(1*+), Sebastiaan J. van de Velde(2,3+), William M. Berelson(1), David J. Bottjer(1), Frank A. Corsetti(1)

1: Department of Earth Sciences, University of Southern California 

2: Department of Geosciences, Environment and Society, Université Libre de Bruxelles

3: Operational Directorate Natural Environment, Royal Belgian Institute of Natural Sciences

*Corresponding author: cribb@usc.edu 

+The R scripts in this repository were created by A.T.C. and S.J.V.

<a href="https://doi.org/10.5281/zenodo.6129380"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.6129380.svg" alt="DOI"></a>

The purpose of this model is to simulate the oxygen penetration depth (OPD) given our paramaterized Ediacaran-Cambrian bioturbation values. This repository contains the biogeochemical reactive-tranpsort model, sensitivity experiment simulation scripts, figure plotting scripts, and output files needed to reproduce the results in "Ediacaran-Cambrian bioturbation did not oxygenate sediments in shallow marine ecosystems" by Cribb et al. The repository contains three directories: <b>1) Datasets</b> - the ichnotaxa dataset with biomixing and bioirrigation assignments, reformatted from Buatois et al. (2020) Science Advances, and a small .csv file of ogranized biodiffusion coefficients from BI/ii and mixed layer depths (Table S1 in the paper); <b>2) Model</b> - the carbon, oxygen, nitrogen, and sulfur reactive-transport model; <b>3) Simulations</b> - the eight scripts used to produce Figures 2-4 and supplementary material 6; and <b>4) Figures</b> - the scripts used to create Figures 2, 3, 4, and 5, and supplementary material figures S1-3. Note that some figure editing was done in Illustrator after the figures were created in R, but you should be able to broadly reproduce all figures in the paper.

No data input is needed, but details of how we paramaterized bioturbation in these models is explained in the manuscript text. Before running this model, set your working directory to a folder containing the model (Model_Function.R) and the simulation scripts. I recommend an empty folder for this, as each script will save an .RData output to the working directory.  You must run the Model_Function.R script to have to function in your environment in order for the simulations scripts to work. To recreate the figures, the output .RData files must be in the same working directory as the figure scripts. Each simulation script will also write the simulation output to an excel file. If you do not wish to save your results in an excel tfile, I recommend you comment out that line of code, which is at the end of each script.

# Datasets
There are two .csv files in <b>Datasets</b>:
* Biodiffusion_MixedLayerDepths.csv - the data displayed in Table S1 of the supplemental
* Ichnotaxa_dataset.csv - the entire Ediacaran and Terreneuvian trace fossil dataset taken from Buatois et al. (2020) Sciences Advances, reformatted into a workable format for the code here.

# Model
The following package is required to run the model:
* ReacTran

This file contains the model function used for all simulations:

* Model_Function.R - This script contains the function needed to run the OPD simulations. It is a reactive-transport model, which will solve differential equations for the following seven carbon, oxygen, nitrogen, and sulfur cycling reactions:
  * R1: CH2O + O2 -> HCO3 + H
  * R2: CH2O + (4/5)NO3 -> HCO3 + (2/5)N2 + (2/5)H2O + (1/5)H
  * R3: CH2O + (1/5)SO4 -> HCO3 + (1/2)HS + (1/2)H
  * R4: HS + 2O2 -> SO4 + H
  * R5: NH4 + 2O2 -> NO3 + H2O + 2H
  * R6: FeS + (9/4)O2 + (3/2)H2O -> FeOOH + SO4 + 2H
  * R7: HS + (8/5)NO3 + (3/5)H -> (4/5)N2 + SO4 + (4/5)H2O

IMPORTANT: Do not change the model function script. If you would like to edit any of the model for your use, do so in the simulation files.

# Simulations
The following packages are required to run the OPD simiulation scripts:
* ReacTran
* rootSolve
* marelac
* writexl

To use these simulations for your own use, you simply need to change the biodiffusion and bioirrigation variables under "Bioturbation Options" at the top of each script. 

This will produce the data in Figure 2 (biodiffusion sensitivity analysis) and Figure S1 (bioirrigation sensitivity analysis):
* Sim_Dba0sensitivitytests.R - file ultimately will produce three .RData files, one for the biodiffusion sensitivity analysis, and two for the bioirrigation sensitivity analysis (one Ediacaran and one Terreneuvian)

The following files were used to produce the data in Figures 3. Each simulation will produce an .RData output file that is used in the associated plotting output figure script. Each lists which figure it will produce the output for, and what the bioturbation parameters (biodiffusion(Db,0) and bioirrigation(a0) are:
* Sim1_NoBioturbation - No bioturbation; Db,0=0, a0=0
* Sim2_EdiacaranBiomixing.R - Ediacaran biomixing; Db,0=0.1 cm2 yr-1, a0=0 yr-1
* Sim3_EdiacaranBioirrigation.R - Ediacaran bioirrigation; Db,0=0 cm2 yr-1, a0=3.65 yr-1
* Sim4_EdiacaranBioturbation.R - Ediacaran biomixing + bioirrigation; Db,0=0.1 cm2 yr-1, a0=3.65 yr-1
* Sim5_TerreneuvianBiomixing.R - Terreneuvian biomixing; Db,0=0.98 cm2 yr-1, a0=0 yr-1
* Sim6_TerreneuvianBioirrigation.R - Terreneuvian bioirrigation; Db,0=0 cm2 yr-1, a0=35.77 yr-1
* Sim7_TerreneuvianBioturbation.R - Terreneuvian biomixing + bioirrigation; Db,0=0.98 cm2 yr-1, a0=35.77 yr-1
* Sim8_ModernBiomixing.R - Modern biomixing; Db,0=10 cm2 yr-1, a0=0 yr-1
* Sim9_ModernBioirrigation.R - Modern bioirrigation; Db,0=0 cm2 yr-1, a0=365 yr-1
* Sim10_ModernBioturbation.R - Modern biomixing + bioirrigation; Db,0=10 cm2 yr-1, a0=365 yr-1

This will produce the data in Figure 4:
* Sim_Figure4_O2EndMembers.R - this script will produce both the high bottom-water oxygen endmember and low bottom-water oxygen endmember results in Figure 4. This will produce two .RData output files that are used in the plotting output figure script.

This will produce the oxygen consumption profiles in Figure 5:
* Sim_Figure5_O2ConsumptionProfiles.R - this script will produce all 10 .RData output files of oxygen consumption curves used to produce Figure 6. 

Thus will produce the oxygen flux data in Figure S4:
* Sim_O2irrflux.R

Theses will produce the burial rates in Figure S5:
* Sim_BurialRates_nobioturbation.R
* Sim_BurialRates_Ediacaran.R
* Sim_BurialRates_Terreneuvian.R

If for any of these simulations you receive the error "Error in steady.1D(y = yini, func = model, parms = PL, nspec = PL$N.var,  : 
  object 'model' not found", you need to load the model function from Model_Function.R into your environment.

# Figures
The following packages are required to recreate the figures in this paper:
* ggplot2
* tidyverse
* ggmap
* MASS
* metR
* egg
* tayloRswift

Before running any of the scripts that reproduce the figures, download both of the files in <b>Datasets</b> and put them in your working directory. 

This script will produce the biomixing vs. bioirrigation analysis in Figure 1:
* Figure1_Biomixing_vs_Bioirrigation.R - This script contains the analysis and plotting output to recreate Figure 1.

This script will produce Figure 2 and Figure S2:
* PlottingOutput_Figure2.R - This script will reproduce the biodiffusion sensitivity analysis for different oxygen levels in Figure 2.
* PlottingOutput_FigureS2.R - This script will reproduce the bioirrgation sensitivity analysis for different oxygen levels at Ediacaran and Terreneuvian biodiffusion coefficients for Figure S1.

This script will produce the figures in Figure 3 and Figure S3:
* PlottingOutput_Figure3_FigureS2.R - This script will reproduce the fan diagrams/contour plots of Figure 3 and Figure S2. You will need to select which output file you would like to plot from the simulation files by commenting/un-commenting the dataset you want. For example, the default set here is to plot the no bioturbation output for Figure 2. Note that the simulation files produce OPD outputs for a larger range of organic matter flux and bottom-water oxygen concentrations, and these scripts will subset the data for the Corg=100-450 umol cm-2 yr-1 and O2=0.014-0.28 mM shown in the figures and discussed in the text of the paper. Lower bottom-water oxygen concentrations and organic matter fluxes are not applicable to the Ediacaran and early Cambrian.

This script will produce Figure 4:
* PlottingOutput_Figure4.R - This script will reproduce Figure 4 and uses the two end member .RData files produced by Sim_Figure4_O2EndMembers.R

This will produce the oxygen consumption profiles in Figure 5 and Figure S3:
* PlottingOutput_Figure5_FigureS3 - This script will reproduce the oxygen consumption profiles in the paper. It combines the consumption rates of all reoxidation reactions into a summated reoxidation rate, and plots consumption of oxygen via aerobic respiration and all reoxidation reactions. This script will use all ten .RData files produced from Sim_Figure5_O2ConsumptionProfiles.R and create ten different oxygen consumption profiles for the upper 1 cm of the sediment, which are all given in Figure S3.

This will produce the oxygen flux figure in Figure S4:
* PlottingOUtput_FigureS4.R - This needs the .RData files from Sim_O2irrflux.R

This will produce the burial rates profile in Figure S5:
* PlottingOutput_FigureS5.R - This needs the .RData files from Sim_BurialRates_[nobioturbation/Ediacaran/Terreneuvian].R

# Contact
If you have any questions, please feel free to contact me at cribb@usc.edu
