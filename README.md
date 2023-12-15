# Homolog-phylogeny
A pipeline to construct a phylogenomic tree from genomic sequences, transcripts, and proteins.

In most cases, phylogenetic relationships are based on amino acid sequences of multi-genes, but it is time-consuming and complicated for researchers to predict eukaryotic genes from the genome. The commonly used methods include MAKER and EVidenceModeler. Thanks to the newly published Protein-to-genome aligner by Dr. Heng Li, miniprot, it really benefits us to obtain the homolog-based prediction in a few minutes. The other point is that Wellcome Sanger Institute is working on expanding the genomic resources with thousands of organisms. In addition, constructing a good tree is complicated and not user-friendly for most biologists, especially for the mushrooming data. Here, we develop a Python-based and one-step pipeline, which could be used to construct a phylogenomic tree from different sources, including genomic sequences, transcripts, and proteins.

We did the benchmark in the genomic-based phylogeny. 1) the assembled genome (draft genome shall be ok, even based on short-reads), 2) the newly sequenced samples (10X depth in short-reads is recommended if most of the samples in tree construction are draft genomes).

Dependencies: java, miniprot, python, CD-HIT, TransDecoder, OrthoFinder, FastTree, IQ-TREE2, Mafft, BMGE, HmmCleaner (optional), BioPerl, uniqHaplo, AlignmentCompare.


Situation: 
1) if all inputs are proteins, it should work in any organisms, including prokaryotes and eukaryotes.
2) if some of the inputs are transcripts. the genetic_code in TransDecoder should be adjusted. In this pipeline, we adopt the Universal (default).
3) if some of the inputs are genomic sequences, you should add the database for alignment.
          

#Installation

##Pre-installation: mamba (highly suggested) or conda. Link: https://github.com/conda-forge/miniforge#mambaforge 

Command: mamba env create --name phylogenomics -f environment.yml  #Once finished, a new environment named phylogenomics will be created, with most dependencies installed. 
Command: mamba activate phylogenomics (if mamba is not installed in your system, use conda)

##note: We integrate many software and packages into the pipeline. Of them, only HmmCleaner.pl could not be configured by conda/mamba. 

##Installation of HmmCleaner.pl: 
1) chmod +x ./dependencies/cpanm 
2) cpan Bio::MUST::Apps::HmmCleaner (This step might take ~20 minutes; be patient; this installation always fails; no worried about that)
3) ./dependencies/cpanm Bio::MUST::Apps::HmmCleaner --force (Try HmmCleaner.pl to check whether it was executable without errors. If errors, it will not produce results)

##As HmmCleaner.pl is unnecessary, and the installation cannot always be finished properly, I write it as the optional step in this pipeline. If such a file is not found or is not executable in your system or environment, it will automatically skip. You don't have to do anything. 

##If you insist on installing it, please see the guidelines at https://metacpan.org/release/ARNODF/Bio-MUST-Apps-HmmCleaner-0.180750/source/INSTALL

#Usage: usage: homolog-phylogenomics.py [-h] [-p PREFIX] [-t THREADS] [-i INPUT] [-m MIN_TAXA] [-l LENGTH_CUTOFF] [-g GENETIC_CODE] [-d DATABASE]

options:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        The prefix used in the output (Required)
  -t THREADS, --threads THREADS
                        Threads used in running (Required, defualt: 40).
  -i INPUT, --input INPUT
                        Files containing sequences for tree construction (Required, must be in the working directory, default: raw).
  -m MIN_TAXA, --min_taxa MIN_TAXA
                        The taxon threshold in parition (Required, defualt: 2/3 of the total inputs).
  -l LENGTH_CUTOFF, --length_cutoff LENGTH_CUTOFF
                        The length threshold in parition (Required, defualt: 100).
  -g GENETIC_CODE, --genetic_code GENETIC_CODE
                        Genetic code for proteins prediction from transcripts, which might be different based on phylum, please check by "TransDecoder.LongOrfs -h" (Optional if only
                        proteins and genomic sequences as inputs; Required if transcripts existed in inputs, default: Universal)
  -d DATABASE, --database DATABASE
                        Proteins sequences for homolog prediction from genomic sequences, it is suggested as proteins from its/their close relatives (three organisms from the the same
                        genus, family, order, class or phylum are suggested, from public data) (Optional if proteins or transcripts as inputs; Required if genomic sequences existed in
                        inputs; It must be provided with the absolute path).

Input: a folder (default 'raw' if not specified by -i) containing sequences. 
We define the rule of three sources with specific suffixes. 1) genomic fasta: species_name.genomic.fasta. 2) transcripts: species_name.transcript.fasta. 3) proteins: species_name.pep.fasta. 
Note: species_name should be identical to others, otherwise it will fail in the tree visualization. We recommend naming the input files as the below rules.
1) genus_species.genomic/transcript/pep.fasta
2) If more than one input from the same species, try genus_species_1.genomic/transcript/pep.fasta and genus_species_2.genomic/transcript/pep.fasta
3) To distinguish from assembly method or source, try genus_species_megahit.genomic/transcript/pep.fasta and genus_species_trinity.genomic/transcript/pep.fasta

Output: prefix.FastTree.full.tre, prefix.IQTREE2.full.tre


Publication:

If any questions, feel free to post an issue or email to ylify@connenct.ust.hk

#Remark: Please cite the integrated software (below) in this pipeline if you will include this pipeline, with doi or website listed.

Bioconda: https://doi.org/10.1038/s41592-018-0046-7

General shell pipeline: https://doi.org/10.1093/sysbio/syw079

AlignmentCompare: https://github.com/DamienWaits/Alignment_Compare.git

BMGE: https://doi.org/10.1186/1471-2148-10-210

cd-hit: https://doi.org/10.1093/bioinformatics/bts565

FastTree: https://doi.org/10.1371/journal.pone.0009490

IQ-TREE 2: https://doi.org/10.1093/molbev/msaa015

miniprot: https://doi.org/10.1093/bioinformatics/btad014

OrthoFinder: https://doi.org/10.1186/s13059-019-1832-y

TransDecoder: https://github.com/TransDecoder/TransDecoder.git

uniqHaplo: http://raven.wrrb.uaf.edu/~ntakebay/teaching/programming/perl-scripts/uniqHaplo.pl
