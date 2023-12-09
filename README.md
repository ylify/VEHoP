# Homolog-phylogeny
A pipeline to construct a phylogenomic tree from genomic sequences，transcripts, and proteins.

In most cases, phylogenetic relationships are based on amino acids sequences of multi-genes, but it is time-consuming and complicated for researchers to predict eukaryotic genes from the genome. The commonly used methods include MAKER and EVidenceModeler. Thanks to the newly published Protein-to-genome aligner by Dr. Heng Li, miniprot, it really benefits us to obtain the homolog-based prediction in few minutes. The other point is that Wellcome Sanger Institute is working on expanding the genomic resources with thousands of organisms. Here, we develop a Python-based and one-step pipeline, which could be used to construct a phylogenomic tree from different sources, including genomic sequences, transcripts, and proteins.

Dependencies: miniprot, python3, CD-HIT, TransDecoder, OrthoFinder, FastTree, IQ-TREE2, Mafft, BMGE, HmmCleaner (optional), uniqHaplo, AlignmentCompare.

Usage: python3 homology-phylogeny.py prefix database (with a folder "raw" containing genomic sequences)

Output: prefix.FastTree.full.tre, prefix.IQTREE2.full.tre

#Installation
We integrate lots of software and packages in the pipeline. Of them, only HmmCleaner.pl could not be configured by conda/mamba, which needs administation access and cpanm.

HmmCleaner.pl

Pre-installation: mamba (highly suggested) or conda. Link: https://github.com/conda-forge/miniforge#mambaforge

Command: mamba env create --name phylogenomics -f environment.yml (if mamba is not installed in your system, use conda)
