# cactus2trees

Pipeline for generating locus alignments and gene trees from a whole genome alignment in Hierarchical Alignment (HAL) Format and produced using Progressive Cactus

- Inputs:
  - HAL file with aligned genomes
  - BED file with gene features in reference genome
  - BED file with sliding window features in reference genome (steps for creating this are shown below)
  - optionally, BED files with mask features for post hoc genome-specific masking
    - I used this step to hardmask sites in chimera genomes assembled using pseudo-it. After masking, these genomes are no longer chimeric.

- Outputs:
  - locus alignment (one or more alignment blocks in separate fasta file) for each gene/window locus (subregion of full alignment intersecting gene/window features)
    - if masking features supplied, a copy of each locus alignment is generated with sequence-specific masking
  - maximum-likelihood gene trees inferred for locus alignments with at least four sequences from different genomes
    - if locus alignment has more than one sequence from the same genome (potentially resulting from gene duplications/losses) only one sequence (the first) is retained for tree inference

##### Steps

Convert .hal to taf.gz using cactus-hal2maf

```
# HAL alignment (input file)
HAL=lamps_genomes58_alignment_final.hal

# TAF alignment (output file)
TAF=lamps_genomes58_alignment_final.taf.gz

# Name of genome to use as the reference genome \\
# Pick an annotated genome that has gene features in a gff/gft/bed feature table file
REFERENCE_GENOME_NAME=Pantherophis-alleghaniensis-PanAll1

# Convert hal to taf.gz with cactus-hal2maf (installed with cactus)
# you could convert to maf instead but the maf file is huge
# I ran cactus-hal2maf on the NRP cluster and all subsquent steps on my laptop or mendel cluster
cactus-hal2maf ~/jobStore/ $HAL $TAF --refGenome $REFERENCE_GENOME_NAME --chunkSize 500000 --batchCores 20  --noAncestors --dupeMode "all"

```

<!--

Convert taf.gz to Multiple Alignment Format (.maf) using taffy

```
### Use the program taffy for taf indexing, conversion to maf, maf indexing, and extracting subregions from maf

### input/output files

# taf.gz alignment (input file)
TAF=lamps_genomes58_alignment_final.taf.gz

# MAF alignment (output file)
MAF=lamps_genomes58_alignment_final.maf

##### do stuff

# create taf index file (taf.gz.tai)
taffy index -i $TAF

# convert taf.gz to maf
taffy view --inputFile $TAFGZ_PATH_LOCAL --outputFile $MAF_PATH_LOCAL --maf

# create maf index file (maf.tai)
taffy index -i $MAF

```
-->

Extract locus alignments from taf.gz at genes in reference genome and save each as a maf

```
# taf alignment (input file)
TAF=lamps_genomes58_alignment_final.taf.gz

# reference genome used during conversion from .hal
REFERENCE_GENOME_NAME=Pantherophis-alleghaniensis-PanAll1

# bed format feature table with gene and possibly other features annotated in your reference genome (input file)
BED=Pantherophis-alleghaniensis-PanAll1.bed

# bed file with gene features exctracted from full bed (output file)
BED_GENES=Pantherophis-alleghaniensis-PanAll1_genes.bed

# where to save maf locus (genes) alignments
OUTPUT_DIR=~/genes_scaffolds-softmasked_repeats-softmasked_MAF/

# make bed with gene features
awk -F'\t' '$8=="gene"{print}' $BED > $BED_GENES

# create taf index file (taf.gz.tai)
taffy index -i $TAF

# Extract locus alignments and save each as a separate maf file with filename indicating \\
# genomic coordinates of locus relative to the reference genome

NUMLOCI=$(cat "$BED_GENES" | wc -l)
for i in $(seq 1 $NUMLOCI);
do
   BEDi=$(sed "${i}q;d" $BED_GENES)
   REGIONi=$(echo "$BEDi" | awk '{print $1":"$2"-"$3}')
   GENEi=$(echo "$BEDi" | awk '{print $10}')
   echo "i:"$i "region:"$REGIONi "gene:"$GENEi
   GENEi_MAF_PATH=$OUTPUT_DIR"/gene-"$GENEi"-"$(echo $REGIONi | sed 's|:|-|g')".maf"
   [[ ! -f "$GENEi_MAF_PATH" ]] && taffy view --inputFile $TAF --outputFile $GENEi_MAF_PATH --maf --region $REFERENCE_GENOME_NAME"."$REGIONi
done
```

Extract locus alignments in 10kb sliding windows along reference genome seqs

```
# HAL alignment (input file)
HAL=lamps_genomes58_alignment_final.hal

# TAF alignment (input file)
TAF=lamps_genomes58_alignment_final.taf.gz

# Reference genome used during conversion from hal to taf.gz
REFERENCE_GENOME_NAME=Pantherophis-alleghaniensis-PanAll1

# Where to save bed file with window features in reference genome (output file)
WINDOWS_BED_PATH=Pantherophis-alleghaniensis-PanAll1_10kb-windows.bed

# directory where window alignments should be saved
OUTPUT_DIR=~/windows-10kb_scaffolds-softmasked_repeats-softmasked_MAF/

# window widths in base pairs
WINDOW_SIZE=10000

# two column tab-separated table with sequence names and lengths for the reference genome
REFCHROMLENGTHS=$(halStats --chromSizes $REFERENCE_GENOME_NAME $HAL | sort -k2 -n -r )

# sequence names for reference genome (=first column of $REFCHROMLENGTHS )
CHROMS=$(echo "$REFCHROMLENGTHS" | awk '{print $1}')

# Create BED file with window features
NUMCHROMS=$(echo "$CHROMS" | wc -l)
for i in $(seq 1 $NUMCHROMS);
do
   CHROMi=$(echo "$CHROMS" | awk -v i=$i 'NR==i{print}')
   CHROMi_LENGTH=$(echo "$REFCHROMLENGTHS" | awk -v CHROMi=$CHROMi '$1==CHROMi{print $2}')
   CHROMi_NUMWINDOWS=$(($CHROMi_LENGTH/$WINDOW_SIZE))
   for j in $(seq 0 $CHROMi_NUMWINDOWS);
   do
      R1=$(($j*$WINDOW_SIZE))
      R2=$(($R1+$WINDOW_SIZE))
      [[ "$R2" -gt "$CHROMi_LENGTH" ]] && R2=$CHROMi_LENGTH
      REGIONij=$CHROMi"\t"$R1"\t"$R2
      [[ -f "$WINDOWS_BED_PATH" ]] && echo "$REGIONij" >> $WINDOWS_BED_PATH
      [[ ! -f "$WINDOWS_BED_PATH" ]] && echo "$REGIONij" > $WINDOWS_BED_PATH
   done
done

# create TAF index file if it doesnt already exist
[[ ! -f ${TAF}".tai" ]] && taffy index -i $TAF

# Extract and save alignment for each window
NUMLOCI=$(wc -l $WINDOWS_BED_PATH | awk '{print $1}')
for i in $(seq 1 $NUMLOCI);
do
   WINDOWS_BED=$(sed "${i}q;d" $WINDOWS_BED_PATH)
   REGIONi=$(echo "$WINDOWS_BED" | awk '{print $1":"$2"-"$3}')
   echo "i:"$i "region:"$REGIONi
   WINDOWi_MAF_PATH=${OUTPUT_DIR}$(echo $REGIONi | sed 's|:|-|g')".maf"
   [[ ! -f "$WINDOWi_MAF_PATH" ]] && taffy view --inputFile $TAF --outputFile $WINDOWi_MAF_PATH --maf --region $REFERENCE_GENOME_NAME"."$REGIONi
done
```

Locus alignments can contain one or more discontiguous subregion alignments (alignment blocks).

The next step extracts and writes each block to a separate MAF in a locus-specific folder.

Example input/output file structure

- /loci_MAF/
  - locus1.maf (with 3 discontiguous alignment blocks)
  - locus2.maf (with 1 alignment block)
  - locus3.maf (with 2 discontiguous alignment blocks)

- /loci_MAF_region-blocks/
  - locus1/
   - locus1_block1.maf
   - locus1_block2.maf
   - locus1_block3.maf
  - locus2/
   - locus2_block1.maf
  - locus3/
   - locus3_block1.maf
   - locus3_block2.maf

```
# directory with locus maf alignments (input files)
MAF_DIR_IN=/windows-10kb_scaffolds-softmasked_repeats-softmasked_MAF/

# base output directory where subdirectories should be created holding alignment blocks for input locus alignments
MAF_BASEDIR_OUT=/windows-10kb_scaffolds-softmasked_repeats-softmasked_MAF_region-blocks/

# create temporary file with paths to input files in $MAF_DIR_IN
MAFs=$(mktemp 2>&1)
find $MAF_DIR_IN -type f > $MAFs

# number of loci
NLOCI=$(wc -l $MAFs | awk '{print $1}')

# set imax to a small number to test on a few loci
imin=1 
imax=$NLOCI

for i in $(seq $imin $imax);
do
   # input path to ith locus alignment
   MAFi=$(sed "${i}q;d" $MAFs)
   
   # current locus name (filename)
   REGIONi=$(basename $MAFi | sed -e 's|.maf$||g')

   # output directory for $MAFi alignment blocks
   OUTDIRi=$(echo $MAF_BASEDIR_OUT"/"$REGIONi)"/"
   mkdir -p $OUTDIRi

   # two-column table with $MAFi line ranges for alignment blocks
   BLOCK_STARTS=$(grep -n ^a $MAFi | cut -f1 -d:)
   FILE_END=$(wc -l $MAFi | awk '{print $1}')
   BLOCK_ENDS=$(cat <(echo "$BLOCK_STARTS" | awk 'NR>1{print $1 - 1}') <(echo "$FILE_END"))
   
   NBLOCKS_MAFi=$(echo "$BLOCK_STARTS" | wc -l)
   for j in $(seq 1 $NBLOCKS_MAFi);
   do
      BLOCK_STARTij=$(echo "$BLOCK_STARTS" | sed "${j}q;d")
      BLOCK_ENDij=$(echo "$BLOCK_ENDS" | sed "${j}q;d")
      BLOCKij_PATH=$OUTDIRi"/"$REGIONi"_block"$j".maf"
      echo $i" "$j
      [[ ! -f "$BLOCKij_PATH" ]] && sed -n "${BLOCK_STARTij},${BLOCK_ENDij}p;${BLOCK_ENDij}q" "$MAFi" > $BLOCKij_PATH
   done
done

```

Sequence-specific alignment masking (if necessary)

```


```

Convert mafs to fasta

```


```

Infer gene tree for each locus with each alignment block as a separate partition

```


```

