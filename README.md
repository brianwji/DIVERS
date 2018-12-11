# DIVERS: Decomposition of Variance Using Replicate Sampling
MATLAB code for DIVERS (Decomposition of Variance Using Replicate Sampling), including absolute abundance estimation from spike-in sequencing and variance/covariance decompostion of absolute bacterial abundances. We have tested this code for MATLAB_R2016a.

## MATLAB notebooks

### Absolute abundance estimation from spike-in sequencing
[_matlab/matlab_notebooks/preprocess_OTU_relative_abundance_table.ipynb_](matlab/matlab_notebooks/preprocess_OTU_relative_abundance_table.ipynb)

### DIVERS 3 component variance/covariance decomposition 
[_matlab/matlab_notebooks/DIVERS_decomposition.ipynb_](matlab/matlab_notebooks/DIVERS_decomposition.ipynb)

### DIVERS 2 component variance/covariance decomposition 
[_matlab/matlab_notebooks/DIVERS_2comp_decomposition.ipynb_](matlab/matlab_notebooks/DIVERS_2comp_decomposition.ipynb)

### Data analysis and figures from the paper 
[_matlab/matlab_notebooks/plot_main_figures.ipynb_](matlab/matlab_notebooks/plot_main_figures.ipynb)

## MATLAB scripts

### DIVERS.m 

```
usage: ./DIVERS.m 

 *User required to specify input and saving directories

 INPUTS: 1) Three absolute abundance tables (data_X, data_Y, data_Z) of
            equal size
              a) data_X and data_Y are technical replicates of the same
              biological samples (measured at every time point of a
              longtitudinal microbiome study)
              b) data_Z is a second replicate (from a second spatial
              location at every time point of a longitudinal microbiome
              study)
   
          *Assumes taxon (OTU) identifiers are provided in the first
          column and full taxonomies are provided in the last column

  OUTPUTS: 1) Variance decomposition of each taxon (DIVERS_variances.txt)
              a) Average abundances of each taxon
              b) Total abundances variances of each taxon
              c) Temporal, spatial, technical variances of each taxon

           2) Covariance decomposition for all pairs of taxa 
              a) Total correlation matrix between all pairs of taxa
              (DIVERS_cormat_total.txt)
              b) Temporal correlation matrix between all pairs of taxa
              (DIVERS_cormat_T.txt)
              c) Spatial correlation matrix between all pairs of taxa
              (DIVERS_cormat_S.txt)
              d) Technical correlation matrix between all pairs of taxa
              (DIVERS_cormat_N.txt)

          *Covariance decomposition output reflects abundant OTUs (log10
          mean absolute abundance > -4). This value was informed by the
          variance decomposition results. 

          *For large data sets, filtering of abundant OTUs may be
          required before covariance decomposition analysis
```
## DIVERS_dual.m

```

DIVERS: Biological (non-technical) and technical variance and covariance
  contribution estimates from absolute bacterial abundance data

  *User required to specify input and saving directories

 INPUTS: 1) Two absolute abundance tables (data_X, data_Y) of
            equal size
              a) data_X and data_Y are technical replicates of the same
              biological samples (measured at every time point of a
              longtitudinal microbiome study)
   
          *Assumes taxon (OTU) identifiers are provided in the first
          column and full taxonomies are provided in the last column

  OUTPUTS: 1) Variance decomposition of each taxon (DIVERS_dual_variances.txt)
              a) Average abundances of each taxon
              b) Total abundances variances of each taxon
              c) Biological and technical variances of each taxon

           2) Covariance decomposition for all pairs of taxa 
              a) Total correlation matrix between all pairs of taxa
              (DIVERS_dual_cormat_total.txt)
              b) Biological correlation matrix between all pairs of taxa
              (DIVERS_dual_cormat_B.txt)
              c) Technical correlation matrix between all pairs of taxa
              (DIVERS_dual_cormat_N.txt)

          *Covariance decomposition output reflects abundant OTUs (log10
          mean absolute abundance > -4). This value was informed by the
          variance decomposition results. 

          *For large data sets, filtering of abundant OTUs may be
          required before covariance decomposition analysis
```



