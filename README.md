# Ediacaran-Cambrian bioturbation did not oxygenate sediments in shallow marine ecosystems

Alison T. Cribb(1*+), Sebastiaan J. van de Velde(2,3+), William M. Berelson(1), David J. Bottjer(1), Frank A. Corsetti(1)

1: Department of Earth Sciences, University of Southern California 

2: Department of Geosciences, Environment and Society, Université Libre de Bruxelles

3: Operational Directorate Natural Environment, Royal Belgian Institute of Natural Sciences

*Corresponding author: cribb@usc.edu 

+The R scripts in this repository were created by A.T.C. and S.J.V.

Current release (initial submission):
<a href="https://zenodo.org/badge/latestdoi/386430928"><img src="https://zenodo.org/badge/386430928.svg" alt="DOI"></a>

The purpose of this model is to simulate the oxygen penetration depth (OPD) given our paramaterized Ediacaran-Cambrian bioturbation values. This repository contains the biogeochemical reactive-tranpsort model, sensitivity experiment simulation scripts, figure plotting scripts, and output files needed to reproduce the results in "Ediacaran-Cambrian bioturbation did not oxygenate sediments in shallow marine ecosystems" by Cribb et al. The repository contains three directories: <b>1) Model</b> - the carbon, oxygen, nitrogen, and sulfur reactive-transport model; <b>2) Simulations</b> - the eight scripts used to produce Figures 2-4 and supplementary material 6; and <b>3) Figures</b> - the scripts used to create Figures 2, 3, 4, and supplementary material 6. Note that some figure editing was done in Illustrator after the figures were created in R, but you should be able to broadly reproduce all figures in the paper.

No data input is needed, but details of how we paramaterized bioturbation in these models is explained in the manuscript text. Before running this model, set your working directory to a folder containing the model (Model_Function.R) and the simulation scripts. I recommend an empty folder for this, as each script will save an .RData output to the working directory.  You must run the Model_Function.R script to have to function in your environment in order for the simulations scripts to work. To recreate the figures, the output .RData files must be in the same working directory as the figure scripts. Each simulation script will also write the simulation output to an excel file. If you do not wish to save your results in an excel tfile, I recommend you comment out that line of code, which is at the end of each script.

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

The following files were used to produce the data in Figure 2 and Figure 3. Each lists which figure it will produce the output for, and what the bioturbation parameters (biodiffusion(Db,0) and bioirrigation(a0) are:
* Sim1_NoBioturbation.R - Produces output for Figure 2. No bioturbation; Db,0=0, a0=0
* Sim2_WeakBiomixing.R - Produces output for Figure 3A. Weak biomixing; Db,0=1 cm2 yr-1, a0=0 yr-1
* Sim3_WeakBioirrigation.R - Produces output for Figure 3B. Weak bioirrigation; Db,0=0 cm2 yr-1, a0=36.5 yr-1
* Sim4_WeakBioturbation.R - Produces output for Figure 3C. Weak biomixing + bioirrigation; Db,0=1 cm2 yr-1, a0=36.5 yr-1
* Sim5_StrongBiomixing.R - Produces output for Figure 3D. Strong biomixing; Db,0=10 cm2 yr-1, a0=0 yr-1
* Sim6_StrongBioirrigation.R - Produces output for Figure 3E. Strong bioirrigation; Db,0=0 cm2 yr-1, a0=365 yr-1
* Sim7_StrongBoioturbation.R - Produces output for Figure 3C. Strong biomixing + bioirrigation; Db,0=10 cm2 yr-1, a0=365 yr-1

This will produce the data in Figure 4:
* Figure4_simulations.R - this script will produce both the high bottom-water oxygen endmember and low bottom-water oxygen endmember results in Figure 4. This will produce two .RData output files that are used in the figure script.

This will produce the oxygen consumption profiles in Figure 6:
* SM6_OxygenConsumption.R - this script will produce all 7 .RData output files of oxygen consumption curves used to produce Figure 6. 

# Figures
The following packages are required to recreate the figures in this paper:
* ggplot2
* tidyverse
* ggmap
* MASS
* metR
* egg

This script will produce the figures in Figures 2 and 3:
* PlottingOutput_Figures23.R - This script will reproduce the fan diagrams/contour plots of Figure 2 and 3. You will need to select which output file you would like to plot from the simulation files by commenting/un-commenting the dataset you want. For example, the default set here is to plot the no bioturbation output for Figure 2. Note that the simulation files produce OPD outputs for a larger range of organic matter flux and bottom-water oxygen concentrations, and these scripts will subset the data for the Corg=100-450 umol cm-2 yr-1 and O2=0.014-0.28 mM shown in the figures and discussed in the text of the paper. Lower bottom-water oxygen concentrations and organic matter fluxes are not applicable to the Ediacaran and early Cambrian.

This script will produce Figure 4:
* PlottingOutput_Figure4.R - This script will reproduce Figure 4 and uses the two end member .RData files produced by Figure4_simulations.R

This will produce the oxygen consumption profiles in Supplementary Material 6:
* PlottingOutput_SM6.R - This script will reproduce the oxygen consumption profiles in Supplementary Material 6. It combines the consumption rates of all reoxidation reactions into a summated reoxidation rate, and plots consumption of oxygen via aerobic respiration and all reoxidation reactions. This script will use all seven .RData files produced from SM6_OxygenConsumption.R and create seven different oxygen consumption profiles for the upper 1 cm of the sediment.
