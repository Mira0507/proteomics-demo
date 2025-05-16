# proteomics-demo

This workflow demonstrates differential expression analysis of data-independent
acquisition (DIA) proteomics data.  Differential testing is conducted using 
the [Limma-trend](https://academic.oup.com/nar/article/43/7/e47/2414268) statistical
framework, employing linear models and empirical Bayes estimation for individual
proteins.

## Repository structure

```
$ tree
.
├── env.yaml
├── README.md
└── scripts
    ├── helper.R
    └── proteomics_de.Rmd

```

### Supporting Files

- `env.yaml` conda environment
- `scripts/helper.R`: helper functions
- `scripts/proteomics_de.Rmd`
    - Dataset overview
    - Quality control
        - sample similarity heatmap
        - PCA
        - Mean-variance relationship
        - Distribution of residual standard error (RSE) 
        - Pre-filtering
    - Differential testing
    - Data visualization
        - MA plot
        - Volcano plot
        - P-value distribution
    - Result table
