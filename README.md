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

## Setup 

1. Clone the repository

Clone the repository to your local working directory. If your authentication 
method for GitHub is [SSH (Secure Shell Protocol)](https://www.ssh.com/academy/ssh-keys), 
run the following:

```
git clone git@github.com:Mira0507/proteomics-demo.git 
cd proteomics-demo
```

Otherwise, clone using the web URL below:

```
git clone https://github.com/Mira0507/proteomics-demo.git
cd proteomics-demo
```

2. Create Conda environment

The package management of the current workflow relies on conda and mamba. 
Ensure you have conda and mamba ready in your terminal. For more information, 
refer to the following pages:

- [Conda documentation](https://docs.conda.io/projects/conda/en/stable/)
- [Mamba user guide](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html)

If you are ready to set up your main conda environment, follow the command below:

```
# Assume you are in the proteomics-demo directory
$ mamba env create --prefix ./env --file env.yaml
```

This will create a new conda environment named `env` in the current directory.


## Analysis

The analysis is conducted by following the order listed under 
[Analysis scripts](###analysis-scripts).
