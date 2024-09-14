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

### Steps

Convert .hal to taf.gz using the [cactus](https://github.com/ComparativeGenomicsToolkit/cactus) cactus-hal2maf function.

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
- halStats function requires [HAL](https://github.com/ComparativeGenomicsToolkit/hal), which is included in [cactus](https://github.com/ComparativeGenomicsToolkit/cactus).

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
MAF_DIR_IN=/genes_scaffolds-softmasked_repeats-softmasked_MAF/

# base output directory where subdirectories should be created holding alignment blocks for input locus alignments
MAF_BASEDIR_OUT=/genes_scaffolds-softmasked_repeats-softmasked_MAF_region-blocks/

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

# remove temporary file
rm $MAFs

```

Sequence-specific alignment masking (if necessary)

First, create a huge fasta file with masked versions of your unaligned genomes
- requires HAL (v2.1 tested) and bedtools (v2.29.2 tested)

```
# path to HAL alignment (input file)
HAL=lamps_genomes58_alignment_final.hal

# path to a two-column tab-separated file with each row containing genome name (column 1) and filepath to a BED file with intervals to mask in the genome (column 2).
BED_CONFIG_PATH=

# where to save masked-versions of unaligned genomes extracted from HAL (combined in a single fasta file)
GENOMES_PATH=lamps_genomes58_masked-nogaps.fa

# names of all genomes in the HAL
GENOME_NAMES=$(halStats --genomes "$HAL" | sed 's| |\n|g' | sort)

# genome names excluding ancestral genomes
GENOME_TIPNAMES=$(echo "$GENOME_NAMES" | grep -v '^Anc[0-9]*')

# number of tip genomes
NUMTIPGENOMES=$(echo "$GENOME_TIPNAMES" | wc -l)
 
# (1) get bed intervals to mask for each genome and add to $BED_PATH
# (2) get sequence lengths for each genome and append to $CHROMLENGTHS_PATH
# (3) extract each genome and save to $GENOMES_PATH as fasta with UCSC format sequence names
for i in $(seq 1 $NUMTIPGENOMES);
do
   GENOMEi=$(echo "$GENOME_TIPNAMES" | sed "${i}q;d")
   echo $i $GENOMEi

   # test if any intervals to mask in $GENOMEi (0=false, 1=true)
   TEST_BEDPATHi=$(awk -v gen=$GENOMEi '$1==gen{print $1}' $BED_CONFIG_PATH | wc -w)
   [[ "$TEST_BEDPATHi" -gt 1 ]] && echo "genome names in must be unique in "$BED_CONFIG_PATH && exit 1
   
   # (1) add $GENOMEi BED intervals (if any) to $BED_PATH
   [[ "$TEST_BEDPATHi" -eq 0 ]] && BEDPATHi=""
   [[ "$TEST_BEDPATHi" -eq 1 ]] && BEDPATHi=$(awk -v gen=$GENOMEi '$1==gen{print $2}' $BED_CONFIG_PATH)
   [[ "$TEST_BEDPATHi" -eq 1 ]] && [[ -f "$BEDPATHi" ]] && [[ -f "$BED_CONFIG_PATH" ]] && awk -v gen=$GENOMEi '{print gen"."$0}' $BEDPATHi > $BED_PATH
   [[ "$TEST_BEDPATHi" -eq 1 ]] && [[ -f "$BEDPATHi" ]] && [[ ! -f "$BED_CONFIG_PATH" ]] && awk -v gen=$GENOMEi '{print gen"."$0}' $BEDPATHi > $BED_PATH

   # (2) add $GENOMEi sequence names (genomeName.chromosomeName) and lengths to $CHROMLENGTHS_PATH
   [[ $i -eq 1 ]] && halStats --chromSizes $GENOMEi $HAL | awk -v gen=$GENOMEi '{print gen"."$0}' > $CHROMLENGTHS_PATH
   [[ $i -gt 1 ]] && halStats --chromSizes $GENOMEi $HAL | awk -v gen=$GENOMEi '{print gen"."$0}' >> $CHROMLENGTHS_PATH

   # (3) mask $GENOMEi sequences and add to $GENOMES_PATH
   [[ "$TEST_BEDPATHi" -eq 0 ]] && [[ $i -eq 1 ]] && hal2fasta --ucscSequenceNames $HAL $GENOMEi > $GENOMES_PATH
   [[ "$TEST_BEDPATHi" -eq 0 ]] && [[ $i -gt 1 ]] && hal2fasta --ucscSequenceNames $HAL $GENOMEi >> $GENOMES_PATH

   [[ "$TEST_BEDPATHi" -eq 1 ]] && [[ $i -eq 1 ]] && GENOMEi_TEMP_PATH=$(mktemp 2>&1) && 
              bedtools maskfasta -fi <(hal2fasta --ucscSequenceNames $HAL $GENOMEi) -fo $GENOMEi_TEMP_PATH -bed $BED_PATH && 
              cp "$GENOMEi_TEMP_PATH" $GENOMES_PATH && rm $GENOMEi_TEMP_PATH

   [[ "$TEST_BEDPATHi" -eq 1 ]] && [[ $i -gt 1 ]] && GENOMEi_TEMP_PATH=$(mktemp 2>&1) && 
              bedtools maskfasta -fi <(hal2fasta --ucscSequenceNames $HAL $GENOMEi) -fo $GENOMEi_TEMP_PATH -bed $BED_PATH && 
              cat "$GENOMEi_TEMP_PATH" >> $GENOMES_PATH && rm $GENOMEi_TEMP_PATH
done
```


Use the hardmasked unaligned genomes fasta (created in the previous step) to mask sites in your MAF locus alignments. Output alignments = fasta format.
Requires:
- bedtools (tested using version 2.29.2)
- R (tested using version 4.0.2)
- R packages dplyr, stringr, Biostrings, and GenomicRanges
- this R script: [mask-alignment.R](https://raw.githubusercontent.com/JeffWeinell/mask-alignment/refs/heads/main/current/mask-alignment.R)


```
# NOTE: zero-length intervals in input MAF are always filtered

### In the future this script will use inputs \\
# (1) HAL file (instead of fasta with all genomes) \\
# (2) MAF file generated from the HAL \\
# (3) BED with intervals to mask in output alignment fasta

# module load R/R-4.0.2
# module load Bedtools/bedtools-2.29.2

# path to fasta with all hardmasked, unaligned genomes (input file)
GENOMES_PATH=lamps_genomes58_masked-nogaps.fa

# where your R packages are installed
R_PACKAGES_DIR=

# path to your copy of R script mask-alignment.R
MASK_ALIGNMENT_RSCRIPT=mask-alignment.R

# directory containing input MAF files (including any in subdirectories)
MAF_DIR_IN=/genes_scaffolds-softmasked_repeats-softmasked_MAF_region-blocks/

# create temporary file with paths to input files in $MAF_DIR_IN
MAFs=$(mktemp 2>&1)
find $MAF_DIR_IN -type f > $MAFs

# output directory where masked fasta alignments should be saved
ALN_FA_DIR_OUT=/genes_scaffolds-softmasked_repeats-softmasked_fasta_region-blocks/

# first alignment to process
imin=1

# last alignment to process (default 0 to process from $imin to last alignment)
imax=0

# Number of alignments in input directory
NUMLOCI=$(wc -l $MAFs | awk '{print $1}')

# updates $imax
[[ "$imax" -eq 0]] && imax=$NUMLOCI
[[ "$imax" -gt "$NUMLOCI" ]] && imax=$NUMLOCI

for i in $(seq $imin $imax);
do
   MAFi=$(sed "${i}q;d" $MAFs)
   FAi=$ALN_FA_DIR_OUT"/"$(basename "$MAFi" | sed 's|.maf$|.fa|g')
   

   ### exit if $FAi already exists
   [[ -f "$FAi" ]] && echo "exiting because "$FAi" already exists" && exit 0
   
   ### temporary files
   ALN_FA_UNMASKED_PATH=$(mktemp 2>&1)
   REGIONS_TABLE_PATH_A=$(mktemp 2>&1)
   REGIONS_TABLE_PATH_B=$(mktemp 2>&1)
   REGIONS_TABLE_PATH=$(mktemp 2>&1)
   BEDPATH=$(mktemp 2>&1)
   GAPLESS_FA_MASKED_PATH=$(mktemp 2>&1)
   
   # make a table with genomic coordinates of each sequence regions and some other info
   ALNFILEi=$(basename "$MAFi")
   zcat -f $MAFi | awk -F'\t' -v x=$ALNFILEi '$1=="s"{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"x}' | awk '{print $0"\t"$2+$3}' | sed 's|\(^[^.]*\)[.]\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*$\)|\1.\2\t\1\t\2\t\3\t\4\t\5\t\6\t\7\t\8\t\1-\2|g' | awk '{print $0"("$6")/"$4+1"-"$9+1"\tgene\t.\t0"}' | awk '{print $1"\t"$13"\t"$7"\t"$1"\t"$12"\t"$6"\t"$4"\t"$9"\t"$5"\t"$10"\t"$2"\t"$3"\t"$8"\t"$11}' > $REGIONS_TABLE_PATH_A
   zcat -f $MAFi | awk '$1=="s"{print $2"("$5")/"$3+1"-"$3+$4+1}' > $REGIONS_TABLE_PATH_B
   paste $REGIONS_TABLE_PATH_A $REGIONS_TABLE_PATH_B > $REGIONS_TABLE_PATH
   
   # make a BED-format version of the table created in the previous step
   awk '{print $1"\t"$3"\t"$6"\t"$7"\t"$8"\t"$15}' $REGIONS_TABLE_PATH | 
       awk '$3=="-" { $7 = $2-$5 }1' | 
       awk '$3=="-" { $8 = $2-$4 }1' | 
       awk '$3=="+" { $7=$4 }1' | 
       awk '$3=="+" { $8=$5 }1' |
       awk '$7!=$8{print $1"\t"$7"\t"$8"\t"$6"\t.\t"$3}' > $BEDPATH
   
   # extract BED intervals from $GENOMES_PATH containing masked-versions of sequences spanning the same regions as in $MAFi (but without gaps)
   bedtools getfasta -s -nameOnly -fi $GENOMES_PATH -bed $BEDPATH | sed 's|[(][+-][)]$||g' > $GAPLESS_FA_MASKED_PATH
   
   # convert MAFi to fasta alignment
   awk '$1=="s"{print ">"$2"("$5")/"$3+1"-"$3+$4+1"\n"$7}' $MAFi > $ALN_FA_UNMASKED_PATH
   
   ### Idea for future faster version of code
   # Instead of extracting BED intervals from $GENOMES_PATH, it may be able to generate $GAPLESS_FA_MASKED_PATH this way:
   # (1) remove gaps from $ALN_FA_UNMASKED_PATH and save as $GAPLESS_FA_UNMASKED_PATH
   # (2) transform $BED_MASK to $BED_MASKi:
   #       (2.1) Get intersection between $BED_MASK and subsequence features and save as $BEDxi_MASK
   #       (2.2) transform $BEDxi_MASK from chromosomal coordinates to subsequence coordinates and save as $BED_MASKi
   # (3) apply masking to ungapped subsequences:
   #     bedtools maskfasta -fi $GAPLESS_FA_UNMASKED_PATH -fo $GAPLESS_FA_MASKED_PATH -bed $BED_MASKi

   # For each sequence in $ALN_FA_UNMASKED_PATH, replace non-gap sites with corresponding replacement sequence in $GAPLESS_FA_MASKED_PATH
   Rscript $MASK_ALIGNMENT_RSCRIPT $ALN_FA_UNMASKED_PATH $GAPLESS_FA_MASKED_PATH $FAi $R_PACKAGES_DIR
   echo "masked alignment written to: "$FAi
   
   ### remove temporary files
   rm $ALN_FA_UNMASKED_PATH
   rm $REGIONS_TABLE_PATH_A
   rm $REGIONS_TABLE_PATH_B
   rm $REGIONS_TABLE_PATH
   rm $BEDPATH
   rm $GAPLESS_FA_MASKED_PATH
   
   # module unload R/R-4.0.2
   # module unload Bedtools/bedtools-2.29.2
done
```



<!--
Convert mafs to fasta
```
-->


<!--
```
Infer gene tree for each locus with each alignment block as a separate partition

```
-->


<!--

```
-->










