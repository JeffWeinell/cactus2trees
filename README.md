# cactus2trees

Pipeline for generating gene trees from a whole genome alignment in Hierarchical Alignment (HAL) Format and produced using Progressive Cactus

- Inputs:
  - HAL file with aligned genomes
  - BED file with gene features in reference genome
  - BED file with sliding window features in reference genome (steps for creating this are shown below)
  - (*Possibly also:*): Separate BED file for each genome with masking features (indicating sites to hardmask)

- Outputs: two trees datasets produced:
  - gene trees for annotated features (protein coding genes)
  - gene trees for loci in sliding windows along chromosomes/contigs/scaffolds of reference genome

##### Steps

  **1.1.** 





2.1