# cactus2trees

Pipeline for generating gene trees from a whole genome alignment in Hierarchical Alignment (HAL) Format and produced using Progressive Cactus

- Inputs:
  - HAL file with aligned genomes
  - BED file with gene features in reference genome
  - BED file with sliding window features in reference genome (steps for creating this are shown below)
  - *optionally* BED files with mask features for post hoc genome-specific masking
    - I used this step to hardmask sites in chimera genomes assembled using pseudo-it. After masking, these genomes are no longer chimeric.

- Outputs:
  - locus alignment (one or more alignment blocks in separate fasta file) for each gene/window locus (subregion of full alignment intersecting gene/window features)
    - if masking features supplied, a copy of each locus alignment is generated with sequence-specific masking
  - maximum-likelihood gene trees inferred for locus alignments with at least four sequences from different genomes
    - if locus alignment has more than one sequence from the same genome (potentially resulting from gene duplications/losses) only one sequence (the first) is retained for tree inference

##### Steps

**Convert HAL to TAF**

```
#### convert .hal --> taf.gz

# HAL alignment (input file)
HAL=lamps_genomes58_alignment_final.hal

# TAF alignment (intermediate output file)
TAF=lamps_genomes58_alignment_final.taf.gz

# MAF alignment (output file)
MAF=lamps_genomes58_alignment_final.maf

# Name of genome to use as the reference genome \\
# Pick an annotated genome that has gene features in a gff/gft/bed feature table file
REFERENCE_GENOME_NAME=Pantherophis-alleghaniensis-PanAll1

# Convert hal to taf.gz with cactus-hal2maf (installed with cactus)
# you can go straight from hal to maf with this command but the maf file is huge
# I converted to .taf.gz first because it was easier for me to run cactus-hal2maf on \\
# the NRP cluster than on the mendel cluster or my laptop, where I run the other steps, and transfering the taf.gz is fast
cactus-hal2maf ~/jobStore/ $HAL $TAF --refGenome $REFERENCE_GENOME_NAME --chunkSize 500000 --batchCores 20  --noAncestors --dupeMode "all"

```

**Convert TAF to Multiple Alignment Format (MAF) and extract subregion alignments**

```
### Use the program taffy for taf indexing, conversion to maf, maf indexing, and extracting subregions from maf

### input/output files

# TAF alignment (input file)
TAF=lamps_genomes58_alignment_final.taf.gz

# MAF alignment (output file)
MAF=lamps_genomes58_alignment_final.maf

# bed format feature table with gene and possibly other features annotated in your reference genome (input file)
BED=Pantherophis-alleghaniensis-PanAll1.bed

# gene feature rows from $BED (output file)
BED_GENES=Pantherophis-alleghaniensis-PanAll1_genes.bed

##### do stuff

# create taf index file (taf.gz.tai)
taffy index -i $TAF

# convert taf.gz --> maf
taffy view --inputFile $TAFGZ_PATH_LOCAL --outputFile $MAF_PATH_LOCAL --maf

# create maf index file (maf.tai)
taffy index -i $MAF

# get gene feature rows from reference genome annotations bed file
awk -F'\t' '$8=="gene"{print}' $BED > $BED_GENES

```

Extract locus alignments from whole genome alignment at genes of reference genome

```
# full alignment (input file)
MAF=lamps_genomes58_alignment_final.maf

# reference genome used during converion from HAL
REFERENCE_GENOME_NAME=Pantherophis-alleghaniensis-PanAll1

# bed file with gene features in reference genome (input file)
BED_GENES=Pantherophis-alleghaniensis-PanAll1_genes.bed

# Where to save locus alignments
OUTPUT_DIR=~/genes_scaffolds-softmasked_repeats-softmasked_MAF/

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
	[[ ! -f "$GENEi_MAF_PATH" ]] && taffy view --inputFile $MAF  --outputFile $GENEi_MAF_PATH --maf --region $REFERENCE_GENOME_NAME"."$REGIONi
done
```

Extract locus alignments from whole genome alignment at 10kb sliding windows along reference genome seqs

```


# Create BED file defining window features (zero based index)



REFERENCE_CHROM_BED_PATH=/Users/jeffreyweinell/Documents/research-projects/lamps/TwoBit/Pantherophis-alleghaniensis-PanAll1.bed

BED_PATH_OUT="/Users/jeffreyweinell/Documents/lamps/MAF/Pantherophis-alleghaniensis-PanAll1_10kb-windows.bed"

CHROMS=$(awk '{print $1}' "$REFERENCE_CHROM_BED_PATH")

NUMCHROMS=$(echo "$CHROMS" | wc -l)

for x in $(seq 3 $NUMCHROMS);
do
	CHROMx=$(echo "$CHROMS" | awk -v x=$x 'NR==x{print}')
	CHROMx_LENGTH=$(awk -v CHROMx=$CHROMx '$1==CHROMx{print $3}' $REFERENCE_CHROM_BED_PATH)
	zmax=$((CHROMx_LENGTH/10000)) # = number of 10kb windows for CHROMx
	for z in $(seq 0 $zmax);
	do
		R1=$(($z*10000))
		R2=$(($R1+10000))
		[[ "$R2" -gt "$CHROMx_LENGTH" ]] && R2=$CHROMx_LENGTH
		REGIONxz=$CHROMx"\t"$R1"\t"$R2
		[[ -f "$BED_PATH_OUT" ]] && echo "$REGIONxz" >> $BED_PATH_OUT
		[[ ! -f "$BED_PATH_OUT" ]] && echo "$REGIONxz" > $BED_PATH_OUT
	done
done

# Extract alignment blocks in windows
WINDOWS_BED_PATH="/Users/jeffreyweinell/Documents/lamps/MAF/Pantherophis-alleghaniensis-PanAll1_10kb-windows.bed"

# path to MAF alignment
# MAF_PATH_LOCAL=/Users/jeffreyweinell/Documents/lamps/cactus/lamps_genomes58_alignment_final.maf
TAFGZ_PATH_LOCAL=/Users/jeffreyweinell/Documents/lamps/cactus/lamps_genomes58_alignment_final.taf.gz

# create alignment index file if it doesnt exist
[[ ! -f ${TAFGZ_PATH_LOCAL}".tai" ]] && taffy index -i ${TAFGZ_PATH_LOCAL}

# reference species, which must be the first species in each alignment block and the species in which regions defined in the bed file correspond to
REFSPECIES=Pantherophis-alleghaniensis-PanAll1
OUTPUT_DIR=/Users/jeffreyweinell/Documents/lamps/MAF/windows-10kb_scaffolds-softmasked_repeats-softmasked/
NUMWINDOWS=$(wc -l $WINDOWS_BED_PATH | awk '{print $1}')
zmax=$(($NUMWINDOWS/1000))
for z in $(seq 0 $zmax);
do
	R1=$(($z*1000))
	[[ "$R1" -eq 0 ]] && R1=1
	R2=$(($R1+999))
	WINDOWS_BED=$(cat $WINDOWS_BED_PATH | awk -v R1=$R1 -v R2=$R2 'NR>=R1 && NR<=R2 {print}')
	NUMLOCI=$(echo "$WINDOWS_BED" | wc -l)
	for i in $(seq 1 $NUMLOCI);
	do
		REGIONzi=$(echo "$WINDOWS_BED" | awk -v i=$i 'NR==i{print $1":"$2"-"$3}')
		echo "z:"$z "i:"$i "region:"$REGIONzi
		WINDOW_MAF_PATH_LOCAL=${OUTPUT_DIR}$(echo $REGIONzi | sed 's|:|-|g')".maf"
		[[ ! -f "$WINDOW_MAF_PATH_LOCAL" ]] && taffy view --inputFile $TAFGZ_PATH_LOCAL  --outputFile $WINDOW_MAF_PATH_LOCAL --maf --region $REFSPECIES"."$REGIONzi
	done
done
```




