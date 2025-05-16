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
    ├── functional_enrichment.Rmd
    ├── helper.R
    └── proteomics_de.Rmd

```

### Supporting Files

- `env.yaml` conda environment
- `scripts/helper.R`: helper functions


### Analysis scripts

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
    - Result tables

- `scripts/functional_enrichment.Rmd`
    - Enriched gene sets calculated using [clusterProfiler](https://pubmed.ncbi.nlm.nih.gov/22455463/)
        - [GO](https://pubmed.ncbi.nlm.nih.gov/10802651/) Cellular Component (CC)
        - GO Biological Process (BP)
        - GO Molecular Function (MF)
        - [kyoto encyclopedia of genes and genomes (KEGG)](https://pubmed.ncbi.nlm.nih.gov/10592173/)
        - [Disease Ontology (DO)](https://pubmed.ncbi.nlm.nih.gov/22080554/)
    - Result tables
    - Plots
        - [Dotplot](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#dot-plot)
        - [Emapplot](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#enrichment-map)
        - [Cnetplot](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#cnetplot)
