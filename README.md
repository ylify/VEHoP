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
- [Docker image](https://hub.docker.com): `docker pull agnostidae/vehop:1.3`.

---

## Usage

```bash
python3 VEHoP.py -p PREFIX -t THREADS -i INPUT -g GENETIC_CODE -d DATABASE
```

For docker (**host_input_working_dir** should contain all the input files):
```bash
docker run --name vehop -v host_input_working_dir:/container_working_dir -it agnostidae/vehop:1.3 /bin/bash （the interactive image)
```

Help information (shown below)
```bash
python3 VEHoP.py -h 
```

### Required Parameters:
| Option                 | Description                                                                            |
|------------------------|----------------------------------------------------------------------------------------|
| `-p PREFIX`            | Prefix for output files (required).                                                    |
| `-t THREADS`           | Number of CPU threads (default: all).                                                   |
| `-i INPUT`             | Directory containing sequence files (`raw/` by default).                               |
| `-g GENETIC_CODE`      | Genetic code for `TransDecoder` (e.g., Universal, required for transcript datasets when de novo prediction).    |
| `-d DATABASE`          | Protein database for homolog prediction (required if input includes genomes).          |
| `-r READS`             | Raw reads as inputs (tab-delimited).         |
| `-c CONFIGS`           | The customized commands or parameters for integrated software.        |
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
    - Tab-delimited `.txt`, e.g.:
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
   - The intermediate files will be checked while running. Some steps will be skipped if the same inputs are detected.

---

## References

### Recommended Citation:
Yunlong Li, Xu Liu, Chong Chen, Jian-Wen Qiu, Kevin Kocot, Jin Sun.  
**VEHoP: A Versatile, Easy-to-use, Homology-based Phylogenomic Pipeline**.  
bioRxiv 2024.07.24.604968; [DOI:10.1101/2024.07.24.604968](https://doi.org/10.1101/2024.07.24.604968).

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
