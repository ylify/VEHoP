# homology-phylogeny
A pipeline to construct a phylogenomic tree from genomic sequences

In most cases, phylogenetic relationships are based on amino acids sequences if multi-genes, but it is time-consuming and complicated for researchers to predict eukaryotic genes from genome. The commonly used methods include MAKER and EVidenceModeler. Thanks to the newly published Protein-to-genome aligner by Dr. Heng Li, miniprot, it really benefit us to obtain the homolog-based prediction in few minutes. The other point is that Wellcome Sanger Institute is working on draw the genomic resources with thousands of organisms. Here, we develop a python-based and one-step pipeline, which could be used to construct a phylogenomic tree from genomic sequences.

Dependencies: miniprot, python3, cd-hit, TransDecoder, OrthoFinder, FastTree, IQ-TREE2, Mafft, BMGE, Hmmcleaner, uniqHaplo, AlignmentCompare.

Usage: python3 homology-phylogeny.py prefix database (with a folder "raw_genome" containing genomic sequences)
