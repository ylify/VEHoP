# VEHoP (version 1.4)

A **V**ersatile and **E**asy-to-use **Ho**mology-based **P**hylogenomic (VEHoP) pipeline accommodating multiple data types (DNA, RNA, protein sequences, or raw reads).

With the advent of affordable sequencing technologies, massive amounts of data have become publicly available. However, the annotation quality of many gene models remains suboptimal. Institutions such as the Wellcome Sanger Institute and IRADIAN GENOMICS, among others, are working to improve these annotations.

Typically, phylogenetic relationships are inferred from amino acid sequences of multiple genes. However, predicting eukaryotic genes directly from genomes is often time-consuming and complex. VEHoP aims to simplify this process by providing a unified pipeline that accommodates various input types and delivers reliable phylogenetic results.

Our benchmarking demonstrates that VEHoP's genome-based phylogenies are consistent with those based on high-quality protein datasets, offering robust support across all nodes (test datasets have been deposited at Figshare, https://doi.org/10.6084/m9.figshare.26370955.v1 ).

Supported Input Types:

1. Raw reads (e.g., from next-generation or third-generation DNA/RNA sequencing)
2. Assembled genomes or transcriptomes (drafts are acceptable, even those based on short reads)
3. Proteomes from high-quality genomes, accompanied by gene feature files (recommended as mapping databases)
4. Combinations of any of the above sources can be used to construct a single tree, expanding taxonomic coverage. Note that requirements for high-quality genomes or transcriptomes are strict (e.g., sample preservation in liquid nitrogen).

This pipeline leverages ethanol-preserved samples and extensive NGS data from mitochondrial and genome-survey projects.

## Workflow

![VEHoP Workflow](https://github.com/ylify/VEHoP/blob/main/Figs/Fig.1_VEHoP_pipeline_1.jpeg)

## Accuracy Benchmarks

Benchmarks based on oysters, catfish, and insects, compared with two other methods (Fig.2), evaluated with IQ-TREE2 (MFP model; node supports not shown):

![Benchmarks](https://github.com/ylify/VEHoP/blob/main/Figs/Fig.2.png)

## Dependencies

VEHoP requires the following software:
- hifiasm, Megahit, Trinity, Shasta, sra-tools, Trimmomatic, Java, miniprot, Python, CD-HIT, TransDecoder, OrthoFinder, FastTree, IQ-TREE2, MAFFT, BMGE, HmmCleaner (optional), BioPerl, uniqHaplo, AlignmentCompare

## Applicability

1. If all inputs are proteins, VEHoP works for any organism, including prokaryotic or eukaryotic.
2. For transcript inputs, adjust the genetic code in TransDecoder (`-g Universal`). VEHoP automatically adopts TransDecoder for transcript coding potential prediction if needed.
3. When working with genomic sequences or miniprot-derived transcripts, you must provide a protein database for alignment (`-d database`). Otherwise, the pipeline will raise an exception and exit.

## Installation

**Prerequisites:**
- [Mamba](https://github.com/conda-forge/miniforge#mambaforge) (strongly recommended) or Conda.

```shell
git clone https://github.com/ylify/VEHoP.git # or download from Releases
cd VEHoP
mamba env create --name phylogenomics -f environment.yml
# A new environment named 'phylogenomics' will be created with most dependencies installed
mamba activate phylogenomics
# If Mamba is not available, use Conda
```

Most dependencies have been integrated, but `HmmCleaner.pl` cannot be configured via Conda/Mamba.

**Installing HmmCleaner.pl (from VEHoP directory):**

```shell
chmod +x ./dependencies/cpanm 
cpanm Bio::MUST::Apps::HmmCleaner 
# This step might take ~20 minutes; failures during installation are common and not critical
./dependencies/cpanm Bio::MUST::Apps::HmmCleaner --force 
# Test HmmCleaner.pl for executability; errors prevent results, but its usage is optional
```

If you wish to install HmmCleaner.pl, follow the [guidelines here](https://metacpan.org/release/ARNODF/Bio-MUST-Apps-HmmCleaner-0.180750/source/INSTALL).

**Docker Image:**
```shell
docker pull agnostidae/vehop:1.3
```

## Dependency Check

- For local installations:
  ```shell
  python3 VEHoP.py -h
  ```
  It will check the dependencies, show the missing ones, and then install them via mamba automatically.
  
- For Docker (ensure all input files are in `host_input_working_dir`):
  ```shell
  docker run --name vehop -v host_input_working_dir:/container_working_dir -it agnostidae/vehop:1.3 /bin/bash # interactive container
  python /root/app/VEHoP/VEHoP.py -h  
  ```
  Output should display the help information below if dependencies are properly resolved.

## Usage
-
    chmod +x VEHoP.py
      #(if you don't want to call python3 every run)
    python3 VEHoP.py (with absolute path) [-h] [-p PREFIX] [-t THREADS] [-i INPUT] [-m MIN_TAXA] [-l LENGTH_CUTOFF] [-g GENETIC_CODE] [-d DATABASE]
    
    
    options:
          -h, --help
                  show this help message and exit
          -p PREFIX, --prefix PREFIX
                  The prefix used in the output (Required)
          -t THREADS, --threads THREADS
                  Threads used in running (Required, default: 40)
          -i INPUT, --input INPUT
                  Files containing sequences for tree construction (Required, must be in the working directory, default: raw)        
          -m MIN_TAXA, --min_taxa MIN_TAXA
                  The taxon threshold in partition (Required, default: 2/3 of the total inputs)
          -l LENGTH_CUTOFF, --length_cutoff LENGTH_CUTOFF
                  The length threshold in partition (Required, default: 100)
          -g GENETIC_CODE, --genetic_code GENETIC_CODE
                  Genetic code for protein prediction from transcripts, which might be different for the phylum, please check by 
                  "TransDecoder.LongOrfs -h" (If the parameter is given, it will adopt TransDecoder to predict coding potential 
                  in transcripts. Optional if only proteins and genomic sequences as inputs; Required if transcripts existed in inputs, default: Universal) 
          -d DATABASE, --database DATABASE
                  Protein sequences for homolog prediction from genomic sequences, it is suggested that proteins from their close 
                  relatives (three organisms from the same genus, family, order, class, or phylum are suggested, from public data) 
                  (Optional if proteins or transcripts as inputs; Required if genomic sequences existed in inputs; 
                  It must be provided with the absolute path.)
                  (Database will not be included in the matrix and tree)
   
## Configuration

See details in `example.config`.

## Input

Input directory (default: `raw` in working directory) containing sequence files. Use the following suffix conventions:

1. Raw reads or SRA accessions (see `example.reads.txt`)
2. Genomic fasta: `species_name.genomic.fasta`
3. Transcript fasta: `species_name.transcript.fasta`
4. Protein fasta: `species_name.pep.fasta`

**Note:** The `species_name` must be consistent across files for proper tree visualization. Naming conventions:

- `genus_species.genomic/transcript/pep.fasta`
- For multiple inputs from the same species:  
  `genus_species_1.genomic/transcript/pep.fasta`, `genus_species_2.genomic/transcript/pep.fasta`
- To distinguish assembly method/source:  
  `genus_species_megahit.genomic.fasta` (assembled with Megahit),  
  `genus_species_trinity.transcript.fasta` (assembled with Trinity)

## Output

- `$PREFIX.$NUM_INPUTS__$OCCUPANCY.FastTree.full.tre` (FastTreeMP -slow -gamma)
- `$PREFIX.$NUM_INPUTS__$OCCUPANCY.IQTREE2.full.tre` (IQ-TREE2 -m MFP)
- `homolog-phylogenomics.$PREFIX.$NUM_INPUTS__$OCCUPANC.RUN-Day.log` (VEHoP running log)
- `miniprot/`: Results of homolog inference via miniprot (GFF/GFF3 and predicted amino acid sequences)
- `transdecoder/`: Homolog inference results via TransDecoder (predicted amino acid sequences)
- `reads/`: Processed reads (e.g., `_RNA.transcript.fasta`, `_NGS.genomic.fasta`, etc.), organized into folders per input
- `cd-hit/`: Non-redundant amino acid sequences (miniprot/TransDecoder output, cutoff: 0.85)
- `$PREFIX.$NUM_INPUTS.orthofinder`: OrthoFinder results, PhyloPyPruner input files, and logs
- `$PREFIX.$NUM_INPUTS__$OCCUPANCY.Phylogenomics/`: Complete results, including:
  - OG*.fa and OG*.tre (PhyloPyPruner inputs)
  - Backup, checkpoint, rejected taxa, and filtered alignment folders
  - PhyloPyPruner run logs, concatenated matrices, partition info
  - Final trees, shell scripts for ASTRAL and PhyloBayes gene tree analysis
  - AlignmentCompare, BMGE, trimal checkpoints
  - Occupancy heatmap plots

## Tips for Running

1. The script checks for existing intermediate files. If files from previous runs exist with the same input/database, those steps are skipped.
2. More databases result in more predicted proteins, increasing runtime for the miniprot step. Two or three high-quality proteomes are recommended for practical performance.
3. To test different taxonomic occupancies using the same input/database, run in the same working directory. OrthoFinder steps will be skipped if completed previously.

## Example

- Working directory: `/home/yunlongli/Software/VEHoP/test`
- Database: `/home/yunlongli/mollusca_three.pep.fasta`
- Date: 2023-12-25
- Command:
  ```shell
  python3 /home/yunlongli/Software/VEHoP/VEHoP.py -i test -t 40 -m 10 -p mollusca -d /home/yunlongli/mollusca_three.pep.fasta
  ```
- Log file: `/home/yunlongli/Software/VEHoP/test/homolog-phylogenomics.mollusca.40__0.25.2023-12-25.log`
- Results: `/home/yunlongli/VEHoP/test/mollusca.40__0.25.Phylogenomics/phylopypruner/`
  - `/home/yunlongli/Software/VEHoP/test/mollusca.40__0.25.IQTREE2.full.tre`
  - `/home/yunlongli/Software/VEHoP/test/mollusca.40__0.25.FastTree.full.tre`

## Publication

VEHoP: A Versatile, Easy-to-use, Homology-based Phylogenomic Pipeline Accommodating Diverse Sequences  
Yunlong Li, Xu Liu, Chong Chen, Jian-Wen Qiu, Kevin Kocot, Jin Sun  
bioRxiv 2024.07.24.604968; doi: https://doi.org/10.1101/2024.07.24.604968

## Remarks

If you have questions, feel free to open an issue or email ylify@connect.ust.hk.

Please cite any integrated software you use from this pipeline, using the provided DOIs or websites listed below. Since not all dependencies may be included in your analysis, you can check the actual usage in the log file.

- Bioconda: https://doi.org/10.1038/s41592-018-0046-7  
- General shell pipeline: https://doi.org/10.1093/sysbio/syw079  
- AlignmentCompare: https://github.com/DamienWaits/Alignment_Compare.git  
- BMGE: https://doi.org/10.1186/1471-2148-10-210  
- cd-hit: https://doi.org/10.1093/bioinformatics/bts565  
- FastTree: https://doi.org/10.1371/journal.pone.0009490  
- HmmCleaner: https://doi.org/10.1186/s12862-019-1350-2  
- IQ-TREE 2: https://doi.org/10.1093/molbev/msaa015  
- miniprot: https://doi.org/10.1093/bioinformatics/btad014  
- OrthoFinder: https://doi.org/10.1186/s13059-019-1832-y  
- TransDecoder: https://github.com/TransDecoder/TransDecoder.git  
- uniqHaplo: http://raven.wrrb.uaf.edu/~ntakebay/teaching/programming/perl-scripts/uniqHaplo.pl
