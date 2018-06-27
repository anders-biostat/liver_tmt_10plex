# liver_tmt_10plex

Annika Weber, 06/2018

TMT 10plex proteomics project from the Klingmueller group (DKFZ), 80 mice were fed with Standard and Western Diet for 0-26 weeks, liver proteomes were then analysed in order to find the point in time were fibrosis occurs 

### Files
- `Notebook-TMT10plex.Rmd` -- R Notebook, shows how to load the data (MaxQuant peptides table for the 90 samples (peptides.txt), and TMT_10plex randomisation-1.xlsx required), prepare the data for working with it and some first analysis 
- `Notebook-TMT10plex.nb` -- HMTL Notebook version of the above mentioned R Notebook
- `MSnbase.Rmd` -- R Notebook, quantile normalization with the MSn package and following analysis
- `MSnbase.nb` -- HMTL Notebook version of the above mentioned R Notebook
- `calc_density.R` -- R Script, function to calculate neighbourhood density (with Rcpp) and example code to execute the function on

You can find the peptides.txt table, protein tables and the TMT randomisation scheme on the bioinform5: /mnt/Data/Anders_group/u/annika/ZMBH