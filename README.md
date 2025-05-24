# Scripts used for the RADSeq data of Atlantic horseshoe crab （Limulus polyphemus）

In this project, we collected 55 horseshoe crabs across twelve American and two Mexican states. We found two genetically distinct populations on the Atlantic horseshoe crab, reconstructed the evolutionary history, and projected their viability under impending anthropogenic climate change scenarios.

Manuscript title: Population structure of a living fossil, the Atlantic horseshoe crab Limulus polyphemus, and genomic offset in the face of imminent climate change.


## Data availability

Sequencing reads are archived in NCBI under Bioproject XXX. 

To download the data, you may use the command below (may need to install relevant NCBI tools accordingly):

```bash
project='PRJNAXXX'
esearch -db sra -query $project | efetch -format runinfo > runinfo.csv
cat runinfo.csv | cut -d "," -f 1 > SRR.numbers
sed '1d' SRR.numbers | parallel -j 12 fastq-dump --split-files --origfmt --gzip {}
```


## Step 1. NGS reads alignment (Align.txt)

Programs used: 

STACKS (process_radtags): https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php 

BWA (MEM): https://github.com/lh3/bwa

SAMTOOLS: http://www.htslib.org/


## Step 2a. SNP calling and filtering (SNPs.txt)

Programs used: 

STACKS (ref_map.pl): https://catchenlab.life.illinois.edu/stacks/comp/ref_map.php 

PLINK: https://www.cog-genomics.org/plink/


## Step 2b. Site frequency spectrum calculation (SFS.txt)

Program used: 

ANGSD (angsd & realSFS): http://www.popgen.dk/angsd/index.php/SFS_Estimation


## Step 3a. Genetic diversity analyses (Nucleotide diversity.txt)

Packages used: 

Pixy: https://pixy.readthedocs.io/en/latest/


## Step 3b. Reconstruction of demographical history

Scripts for simulating demographical events (Folder "FastSimCoal")

Packages used: 

FastSimCoal2: http://cmpg.unibe.ch/software/fastsimcoal27/


## Step 3c. Reconstruction of species distribution models for the present and the Last Glacial Maximum (LGM) 

Scripts for two species distribution models (Folder "SDM_LGM")

dismo (R package): https://cran.r-project.org/web/packages/dismo/index.html

## Step 3d. Genotype-environment association and viability prediction

Scripts for estimating genome niche index and resistance to dispersal (Folder "Genotype-climate")

Packages used: 

gradientForest (R package): https://rdrr.io/rforge/gradientForest/

resGF (R package): https://github.com/MVan35/resGF/ 
