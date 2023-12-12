# Homolog-phylogeny
A pipeline to construct a phylogenomic tree from genomic sequences, transcripts, and proteins.

In most cases, phylogenetic relationships are based on amino acids sequences of multi-genes, but it is time-consuming and complicated for researchers to predict eukaryotic genes from the genome. The commonly used methods include MAKER and EVidenceModeler. Thanks to the newly published Protein-to-genome aligner by Dr. Heng Li, miniprot, it really benefits us to obtain the homolog-based prediction in a few minutes. The other point is that Wellcome Sanger Institute is working on expanding the genomic resources with thousands of organisms. In addition, constructing a good tree is complicated and not user-friendly for most biologists, especially for the mushrooming data. Here, we develop a Python-based and one-step pipeline, which could be used to construct a phylogenomic tree from different sources, including genomic sequences, transcripts, and proteins.

Dependencies: java, miniprot, python, CD-HIT, TransDecoder, OrthoFinder, FastTree, IQ-TREE2, Mafft, BMGE, HmmCleaner (optional), BioPerl, uniqHaplo, AlignmentCompare.


#Installation

##Pre-installation: mamba (highly suggested) or conda. Link: https://github.com/conda-forge/miniforge#mambaforge 

Command: mamba env create --name phylogenomics -f environment.yml (if mamba is not installed in your system, use conda)

##note: We integrate many software and packages into the pipeline. Of them, only HmmCleaner.pl could not be configured by conda/mamba. 

##Installation of HmmCleaner.pl: run two commands in HmmCleaner_install.sh (it might take ~20 minutes; be patient；try HmmCleaner.pl to check whether it was executable.)

##As HmmCleaner.pl is unnecessary, and the installation cannot always be finished properly, I write it as the optional step in this pipeline. If such a file is not found or is not executable in your system or environment, it will automatically skip. You don't have to do anything. 

##If you insist on installing it, please see the guidelines at https://metacpan.org/release/ARNODF/Bio-MUST-Apps-HmmCleaner-0.180750/source/INSTALL


#Usage: python3 homology-phylogeny.py prefix database (optional, needed if genomic files in raw) (with a folder "raw" containing genomic sequences)

Output: prefix.FastTree.full.tre, prefix.IQTREE2.full.tre

Publication:

If any questions, feel free to post an issue or email to ylify@connenct.ust.hk

#Remark: Please cite the integrated software (below) in this pipeline if you will include this pipeline, with doi or website listed.
Bioconda: https://doi.org/10.1038/s41592-018-0046-7
general shell pipeline: https://doi.org/10.1093/sysbio/syw079
BMGE: https://doi.org/10.1186/1471-2148-10-210
cd-hit: https://doi.org/10.1093/bioinformatics/bts565
miniprot: https://doi.org/10.1093/bioinformatics/btad014
OrthoFinder: https://doi.org/10.1186/s13059-019-1832-y
IQ-TREE 2: https://doi.org/10.1093/molbev/msaa015
FastTree: https://doi.org/10.1371/journal.pone.0009490
uniqHaplo: http://raven.wrrb.uaf.edu/~ntakebay/teaching/programming/perl-scripts/uniqHaplo.pl
AlignmentCompare: https://github.com/DamienWaits/Alignment_Compare.git
TransDecoder: https://github.com/TransDecoder/TransDecoder.git
MAFFT: https://doi.org/10.1093/molbev/mst010
