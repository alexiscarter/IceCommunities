## IceCommunities

This repository contains the computer code and data to replicate the results of the submitted manuscript of Ficetola et al.

For this study, environmental DNA metabarcoding and measurements of microhabitat characteristics were used to assess the development of terrestrial ecosystems emerging after glacier retreats across 46 glacier forelands worldwide.

The scripts available in the `R` folder are ordered in this way:
- `1.manage_data.R` for accessing the data and doing the first manipulation.
- `2.GLMM.R` for the generalized linear mixed-effect models for the biodiversity and soil features
- `3.SEM_alpha_diversity.R` for the structural equation modeling for the alpha diversity
- `4.SEM_beta_diversity.R` for the structural equation modeling for the beta diversity

Files available in the `data` folder:
- The filtered sequence data from all the markers in compressed format
- The phyloseq objects for all the markers
- `div.clim.chem.csv` Table containing biodiversity values, plot information, time since glacier retreat, soil chemistry, temperature, productivity, wetness and plot coordinates. 
- `IC_sampled_points.csv` Raw field data
- `labels.csv` Sample names
- `beta.biotic.plot.csv` Data transformed for beta diversity analyses
- `terrain_swrad.csv` Data on terrain (solar radiation and aspect)
- `tables_fig_3.xlsx` Values used to build the figure on relative mportance of biotic relationships, habitat and time for biodiversity development

Raw sequences are deposited at https://doi.org/10.5281/zenodo.6620359
If you used any of the data, please cite as: Alessia Guerrieri, Aurélie Bonin, Ludovic Gielly, & Gentile Francesco Ficetola. (2022). Raw sequencing data for studying the colonization of soil communities after glacier retreat [Data set]. Zenodo. https://doi.org/10.5281/zenodo.6620359

![](img/IceCom.png)