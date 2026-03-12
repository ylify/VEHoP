# VEHoP (Versatile and Easy-to-use Homology-based Phylogenomic Pipeline)

## Introduction

VEHoP (**V**ersatile and **E**asy-to-use **Ho**mology-based **P**hylogenomic) is a robust pipeline for phylogenetic tree construction from genomic, transcriptomic, proteomic data, or raw reads. By facilitating reliable tree inference even from low-quality or incomplete inputs, VEHoP enables phylogenomic analysis in a broad range of organisms.

### Key Advantages:
- Supports diverse input types: raw sequencing reads, genomic/transcriptomic assemblies, and proteomes.
- Compatible with draft quality genomes and transcriptomes.
- Customizable workflows for homology and phylogenomics, while maintaining reproducibility.
- Benchmark consistency with high-quality protein-based phylogenies.

### Why VEHoP?
While many phylogenomics tools rely on manually curated gene models, this pipeline alleviates those constraints by generating trees directly from raw or minimally processed datasets, expanding accessibility to phylogenomic research.

---

## Workflow Overview

The pipeline follows these major steps:

1. **Data Input**:
   - Accepts raw reads (e.g., NGS, HiFi, ONT), assembled fasta (genomes/transcripts), and high-quality proteomes.

2. **Homology Prediction**:
   - Utilizes tools like `miniprot` for protein inference or `TransDecoder` for coding potential identification.

3. **Sequence Alignment**:
   - Automatically performs alignments with `MAFFT` and conducts quality trimming with `BMGE` and `trimal`.

4. **Pruning and Matrix Creation**:
   - Applies `PhyloPyPruner` to eliminate redundancy, ensure species occupancy, and generate high-quality supermatrices.

5. **Phylogenetic Inference**:
   - Supports maximum likelihood tree construction with `FastTree` and `IQ-TREE2`, and species tree analysis with `ASTRAL`.

7. **Output Results**:
   - Produces gene-specific phylogenies, supermatrices, and concatenated maximum-likelihood trees.

![Workflow Diagram](https://github.com/ylify/VEHoP/blob/main/Figs/Fig.1_VEHoP_pipeline_1.jpeg)

---

## Installation

### Prerequisites:
- **Mamba** (recommended for dependency installation) or **Conda**.
- Software dependencies (detailed in `environment.yml`).

### Step 1: Clone Repository
```bash
git clone https://github.com/ylify/VEHoP.git
cd VEHoP
```

### Step 2: Install Dependencies
```bash
mamba env create --name phylogenomics -f environment.yml
mamba activate phylogenomics
```

### Optional: Install HmmCleaner
```bash
chmod +x ./dependencies/cpanm
./dependencies/cpanm Bio::MUST::Apps::HmmCleaner --force
```

**Other Options:**
- [`environment.yml`](environment.yml) for full dependency management.
- The script will check and install the missing dependencies automatically via mamba.
- [Docker image](https://hub.docker.com): `docker pull agnostidae/vehop:1.3`.
- [Singularity image](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html): `singularity pull library://xuliuouc/collection/vehop`.

---

## Usage

```bash ($PATH should be the absolute path)
python3 $PATH/VEHoP.py -p PREFIX -t THREADS -i INPUT -g GENETIC_CODE -d DATABASE
```

For docker (**host_input_working_dir** should contain all the input files, the absolute path):
```bash
docker run --name vehop -v host_input_working_dir:/container_working_dir -it agnostidae/vehop:1.3 /bin/bash （the interactive image)
```

For singularity (**host_input_working_dir** should contain all the input files, the absolute path):
```bash
singularity run --contain --home host_input_working_dir --bind host_input_working_dir:/data vehop.sif #vehop.sif should be provided with absolute path
python /root/app/VEHoP/VEHoP.py -h
```

For Apptainer (It should be run under the directory with input_working_dir):
```bash
apptainer shell VEHop.sif
source /root/miniforge3/etc/profile.d/conda.sh && conda activate phylogenomics
python /root/app/VEHoP/VEHoP.py -h
``` 

Help information (shown below)
```bash
python3 VEHoP.py -h 
```

Running test (example data)    
The example data (demo.tar.gz and demo_db.fasta.gz) have been deposited at figshare (https://doi.org/10.6084/m9.figshare.28189616).    
Please copy the two files to the directory that consists of VEHoP.py.
```bash
tar -zxvf demo.tar.gz
gunzip demo_db.fasta.gz
python VEHoP.py -p test -t 60 -i raw -d demo_db.fasta
```

### Required Parameters:
| Option                 | Description                                                                            |
|------------------------|----------------------------------------------------------------------------------------|
| `-p PREFIX`            | Prefix for output files (required).                                                    |
| `-t THREADS`           | Number of CPU threads (default: all).                                                   |
| `-i INPUT`             | Directory containing sequence files (`raw/` by default).                               |
| `-g GENETIC_CODE`      | Genetic code for `TransDecoder` (e.g., Universal, required for transcript datasets when de novo prediction).    |
| `-d DATABASE`          | Protein database for homolog prediction (required if input includes genomes).          |
| `-r READS`             | Raw reads as inputs (tab-delimited), cases shown in 'example.reads.txt'.         |
| `-c CONFIGS`           | The customized commands or parameters for integrated software, cases shown in 'example.config'.     |
| `-l LENGTH`            | The length threshold in parition.       |
| `-min MINIMUM_TAXA`    | The taxonic sampling threshold in each parition (default: 2/3 of the total inputs).        |

---

## Accuracy Benchmarks

VEHoP's fidelity has been tested across taxonomically diverse datasets (e.g., oysters, insects, catfish), yielding consistent phylogenetic trees. Comparisons made with alternative approaches demonstrate enhanced robustness and better support for difficult nodes.

![VEHoP Benchmark Results](https://github.com/ylify/VEHoP/blob/main/Figs/Fig.2.png)

---

## Input Guidelines

### Accepted File Types:
1. **Raw Reads**:
    - Tab-delimited `.txt`, e.g. (shown in example.reads.txt):
    ```
    species_1    NGS    /path/to/read_r1.fq    /path/to/read_r2.fq
    species_2    HiFi   /path/to/reads.fq
    ```
2. **Genome/Transcriptome**:  
    - `species_name.genomic.fasta`, `species_name.transcript.fasta`.
3. **Proteomes**:
    - `species_name.pep.fasta` (recommended for alignment databases).

---

## Output Summary

### Key Results:
- Phylogenetic Trees:
  - `{PREFIX}.IQTREE2.full.tre`: Maximum Likelihood tree by **IQ-TREE2**.
  - `{PREFIX}.FastTree.full.tre`: Tree constructed by **FastTree**.
- Alignment results are stored in subdirectories:
  - `miniprot/`: Miniprot-predicted homologs.
  - `transdecoder/`: TransDecoder coding potential predictions.
  - `cd-hit/`: Non-redundant protein clusters.

### Full Logging:
Detailed logs of runtime steps, errors, and parameter settings are stored in:
```bash
homolog-phylogenomics.{PREFIX}.log
```

---

## Tips for Effective Analysis

1. Ensure all sequence files follow uniform naming:
   - Use consistent prefixes (e.g., `genus_species`).
   - Use the other labels if more than one input from the same taxonomy (e.g., `genus_species_1`, `genus_species_2`).
2. Choose databases wisely:
   - High-quality proteomes from related taxa yield the best results.
3. Adjust occupancy (`-m MIN_TAXA`) and length thresholds depending on dataset completeness, under the same working directory:
   - The intermediate files will be checked while running. Some steps will be skipped if the same inputs are detected, reducing the running time.
4. Configure the parameters in each dependency (listed in example.config), only for senior users.


---

### Parameter reference (as in example.config)

#### Assembly
- `shasta_cmd = --config Nanopore-May2022.conf`
  - Description: Options passed to the Shasta assembler.
  - Example flag used here: `--config <file>` — specify a configuration file for Shasta (e.g., tuned for Nanopore reads, such as `Nanopore-May2022.conf`).
  - Note: Provide the correct config file available to Shasta or alter options inline as needed.

- `hifiasm_cmd = --hom-cov auto`
  - Description: Options passed to hifiasm (a HiFi assembler).
  - `--hom-cov auto` tells hifiasm to automatically estimate homozygous coverage. You can instead set an explicit numeric coverage if desired.

- `megahit_cmd= --k-list 21,29,39,59,79,99,119,141 -m 0.8`
  - Description: Options passed to MEGAHIT (short-read assembler).
  - `--k-list 21,29,...` sets the k-mer sizes MEGAHIT will iterate through.
  - `-m 0.8` sets the memory limit (fraction of available memory or memory parameter depending on MEGAHIT version) — commonly used to limit RAM usage; check your installed MEGAHIT version for exact semantics.

- `Trinity_cmd = --seqType fq --max_memory 100G --CPU 6 --full_cleanup`
  - Description: Options passed to Trinity (RNA-seq assembler).
  - `--seqType fq` input sequencing type (fastq).
  - `--max_memory 100G` maximum memory Trinity is allowed to use.
  - `--CPU 6` number of CPU threads to use.
  - `--full_cleanup` remove intermediate files after assembly to save disk.

---

#### Protein-coding region calling
- `miniprot_cmd = -L 30 -j 1 -G 200k`
  - Description: Options for miniprot (protein-to-genome alignment).
  - Common meanings:
    - `-L 30` — minimum alignment length to report (e.g., 30 aa) or similar threshold; helps discard very short hits.
    - `-j 1` — number of threads (here set to 1). Pipelines typically replace this with a parallelism variable.
    - `-G 200k` — maximum intron/gap length (e.g., 200k bases). Adjust according to expected gene structure and genome size.
  - Note: Confirm exact flag semantics with the miniprot version you use.

- `transdecoder_cmd = -m 100`
  - Description: Options for TransDecoder (predicting coding regions/ORFs).
  - `-m 100` — minimum protein length (in amino acids) for retained ORFs; here set to 100 AA.

- `cd-hit_cmd = -c 0.85 -M 50000`
  - Description: Options for CD-HIT (sequence clustering/deduplication).
  - `-c 0.85` — sequence identity threshold (85%). Sequences with identity >= 0.85 are clustered together.
  - `-M 50000` — memory limit in MB (e.g., 50000 MB = 50 GB). Adjust for your environment.

---

#### Ortholog inference
- `orthofinder_cmd = -S diamond`
  - Description: Options for OrthoFinder.
  - `-S diamond` — use DIAMOND for sequence similarity searches (fast alternative to BLAST). OrthoFinder will run all-vs-all searches with DIAMOND.

---

#### Supermatrix construction (alignment processing)
- `uniqHaplo_cmd =`
  - Description: Placeholder for a tool/command that collapses or deduplicates highly similar haplotypes. Blank in example — add options for the tool you use.

- `mafft_cmd = --auto --localpair --quiet --maxiterate 1000`
  - Description: Options for MAFFT multiple sequence alignment.
  - `--auto` — let MAFFT select an appropriate algorithm based on input size.
  - `--localpair` — use iterative refinement method suitable for high accuracy (L-INS-i).
  - `--quiet` — suppress verbose output.
  - `--maxiterate 1000` — run up to 1000 refinement iterations.

- `HmmCleaner_cmd = --specificity`
  - Description: Options for HmmCleaner (or similar HMM-based alignment cleaner).
  - `--specificity` — favor specificity when cleaning (more aggressive removal of suspect columns/regions). Check tool docs for exact behavior.

- `trimal_cmd = -automated1`
  - Description: Options for trimAl (alignment trimming tool).
  - `-automated1` — automated method to choose a trimming strategy appropriate for the alignment.

- `BMGE_cmd = -t AA -g 0.2`
  - Description: Options for BMGE (Block Mapping and Gathering with Entropy).
  - `-t AA` — indicate amino-acid alignment (AA).
  - `-g 0.2` — gap threshold; for example, remove alignment columns with >20% gaps (confirm exact interpretation per BMGE version).

- `AlignmentCompare_cmd =`
  - Description: Placeholder for a tool/command that compares alternative alignments (empty in the example). Add arguments for whichever alignment comparison tool you use.

- `phylopypruner_cmd = --min-support 0.75 --mask pdist --trim-divergent 0.75 --min-pdist 0.01 --prune MI`
  - Description: Options for PhyloPyPruner (tree-based ortholog/pruning tool).
  - Typical meanings:
    - `--min-support 0.75` — minimum node support (e.g., bootstrap or other support metric) to consider clades reliable. Here expressed as fraction (75%).
    - `--mask pdist` — mask sequences or sites based on pairwise distance criteria (mode `pdist`).
    - `--trim-divergent 0.75` — trim sequences diverging beyond this threshold (relative measure).
    - `--min-pdist 0.01` — minimum pairwise distance to retain (avoid near-identical sequences causing artifacts).
    - `--prune MI` — pruning strategy (e.g., MI stands for a particular criterion; check tool docs).
  - Note: These controls tune how strict the pruning is when removing paralogs and divergent sequences from gene trees.

---

#### Phylogenetic relationship
- `FastTreeMP_cmd = -slow -gamma`
  - Description: Options for FastTree (parallel/MP version).
  - `-slow` — run a slower but more thorough tree search.
  - `-gamma` — use the Gamma model of rate heterogeneity across sites.

- `iqtree2_cmd = -B 1000 -m MFP`
  - Description: Options for IQ-TREE 2 (phylogenetic inference).
  - `-B 1000` — run 1000 ultrafast bootstrap replicates (UFBoot).
  - `-m MFP` — run ModelFinder Plus (automatic model selection across many candidate models).

---

### General notes & best practices

- Always confirm the exact meaning of each flag with the documentation for the specific tool version installed on your system; option semantics can change between versions.
- When editing the config:
  - Keep each variable as a single line mapping to the tool options.
  - Do not include the program binary name in the value (the pipeline is expected to prepend the executable name).
  - For parallelism and memory, prefer using pipeline-level variables (if available) to avoid hardcoding thread counts in each tool value. The example shows `-j 1` in `miniprot_cmd` — replace it with your thread variable if the pipeline supports it.
- For blank placeholders (e.g., `uniqHaplo_cmd`, `AlignmentCompare_cmd`), fill in options for the specific tool you plan to use or leave blank to use pipeline defaults.

---

### Example: customizing the config

Suppose you want to:
- Run Trinity with 24 CPUs and 200 GB memory:
  - Edit `Trinity_cmd` to:
    - `Trinity_cmd = --seqType fq --max_memory 200G --CPU 24 --full_cleanup`

- Run CD-HIT clustering at 95% identity and allocate 64 GB RAM:
  - Edit `cd-hit_cmd` to:
    - `cd-hit_cmd = -c 0.95 -M 64000`

- Use IQ-TREE with 1000 standard non-parametric bootstraps instead of ultrafast bootstraps:
  - Modify:
    - `iqtree2_cmd = -b 1000 -m MFP`
  - (Note: `-b` triggers standard non-parametric bootstrap; check IQ-TREE docs.)

---

### Troubleshooting

- If a tool fails with an unknown option, remove or adjust the offending flag and re-run. Version mismatches are a common cause.
- If memory-related crashes occur, lower `-m` (MEGAHIT), `-M` (CD-HIT), or `--max_memory` (Trinity) in the config before re-submitting jobs.
- If assemblies or alignments look poor, consider changing alignment/trimming options (MAFFT, trimAl, BMGE) to be less aggressive; conversely, increase strictness to remove artifactual regions.

---

## References

### Recommended Citation:
Yunlong Li, Xu Liu, Chong Chen, Jian-Wen Qiu, Kevin Kocot, Jin Sun (2026).  
**Reliable Inference of Phylogenomic Relationship via Assembly‐Based Strategy Accommodating Raw Reads and Proteins**.  
Molecular Ecology Resources. [doi:10.1111/1755-0998.70116](https://doi.org/10.1111/1755-0998.70116).    

### Integrated tools:
- **Bioconda**: [doi:10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7)
- **BMGE**: [doi:10.1186/1471-2148-10-210](https://doi.org/10.1186/1471-2148-10-210)
- **cd-hit**: [doi:0.1093/bioinformatics/bts565](https://doi.org/10.1093/bioinformatics/bts565)
- **IQ-TREE2**: [doi:10.1093/molbev/msaa015](https://doi.org/10.1093/molbev/msaa015) 
- **FastTree**: [doi:10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490)  
- **OrthoFinder**: [doi:10.1186/s13059-019-1832-y](https://doi.org/10.1186/s13059-019-1832-y)  
- **miniprot**: [doi:10.1093/bioinformatics/btad014](https://doi.org/10.1093/bioinformatics/btad014)
- **HmmCleaner**: [doi:10.1186/s12862-019-1350-2](https://doi.org/10.1186/s12862-019-1350-2)
- **TransDecoder**: [GitHub Repository](https://github.com/TransDecoder/TransDecoder)
- **AlignmentCompare**: [GitHub Repository](https://github.com/DamienWaits/Alignment_Compare)
- **General shell pipeline**: [doi:10.1093/sysbio/syw079](https://doi.org/10.1093/sysbio/syw079)
