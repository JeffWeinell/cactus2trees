# cactus2trees

Pipeline for generating gene trees from a whole genome alignment in Hierarchical Alignment (HAL) Format and produced using Progressive Cactus

- Inputs:
  - HAL file with aligned genomes
  - BED file with gene features in reference genome
  - BED file with sliding window features in reference genome (steps for creating this are shown below)
  - (*Often also:*): BED files with masking features (regions of unaligned/gapless genomes to hardmask post alignment)

- Outputs:
  - locus alignment (one or more alignment blocks in separate fasta file) for each gene/window locus (subregion of full alignment intersecting gene/window features)
    - if masking features supplied, a copy of each locus alignment is generated with sequence-specific masking
  - maximum-likelihood gene trees inferred for locus alignments with at least four sequences from different genomes
    - if locus alignment has more than one sequence from the same genome (potentially resulting from gene duplications/losses) only one sequence (the first) is retained for tree inference

##### Steps

  **1.1.** 




intervals intersecting 

2.1