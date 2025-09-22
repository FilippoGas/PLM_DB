# HapScoreDB

HapScoreDB is a database providing protein language model–derived fitness scores for haplotype-resolved protein-coding sequences across all human transcript isoforms.

Built using GENCODE and Ensembl annotations combined with phased variant data from the 1000 Genomes Project, HapScoreDB catalogs:
- 130,000+ distinct protein haplotypes
- From 18,000+ genes and 78,000 transcripts
- Covering 94,000+ coding variants
  
Fitness scores for each haplotype were computed with state-of-the-art protein language models.

HapScoreDB includes an intuitive web interface for interactive exploration, visualization, and download of complete or customized datasets.
This repository contains the source code for the RShiny application powering the web interface.

Explore HapScoreDB here: https://bcglab.cibio.unitn.it/hapscoredb
