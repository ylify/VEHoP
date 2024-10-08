# VEHoP (version 1)
A **V**ersatile and **E**asy-to-use **Ho**mology-based **P**hylogenomic (VEHoP) pipeline accommodating multiple types (DNA, RNA, and protein sequences)

With affordable sequencing, mushrooming data is available in the public database. Most of them are not well annotated in the gene model. Wellcome Sanger Institute, IRADIAN GENOMICS, and others are working on expanding the genomic resources with thousands of organisms. How to use these data will be a valuable question to answer.

In most cases, phylogenetic relationships are based on amino acid sequences of multi-genes, but it is time-consuming and complicated for researchers to predict eukaryotic genes from the genome. The commonly used methods include MAKER and EVidenceModeler. Besides, the transcripts from RNA-seq are also widely used. Thanks to the newly published Protein-to-genome aligner by Dr. Heng Li, miniprot, it really benefits us to obtain the homolog-based prediction in a few minutes. In addition, constructing a good tree is complicated and not user-friendly for most biologists, especially for the mushrooming data. Here, we develop a Python-based and one-step pipeline, which could be used to construct a phylogenomic tree from different sources, including genomic sequences, transcripts, and proteins (the first figure).

We did the benchmark in the genomic-based phylogeny (the second figure). The topology is consistent with the result from high-quality proteins, with robust support in all nodes. (test dataset has been deposited in Figshare, https://doi.org/10.6084/m9.figshare.26370955.v1 )

    1) the assembled genome (draft genome shall be ok, even based on short-reads)
    2) the newly sequenced samples (10X depth in short-reads and megahit assembler are recommended if most of the samples in tree construction are draft genomes)
    3) three sources could be used in a single tree, which expands the coverage of taxonomy. The sample requirement for a high-quality genome or transcriptome is strict, for example, liquid nitrogen or RNALATER
This pipeline will take advantage of ethanol-preserved samples and also the massive NGS data from the mitochondrial and genome-survey projects.    
    
Workflow (Fig.1):
-
![image](https://github.com/ylify/VEHoP/blob/main/Pipeline.jpeg)    
Accuracy benchmark based on Ostreida (Fig.2, evaluated by IQ-TREE2 with MFP model, full support in nodes not shown):
-
![image](https://github.com/ylify/VEHoP/blob/main/Fig.2.jpg)    
    
Dependencies: 
-
java, miniprot, python, CD-HIT, TransDecoder, OrthoFinder, FastTree, IQ-TREE2, Mafft, BMGE, HmmCleaner (optional), BioPerl, uniqHaplo, AlignmentCompare, PhyloPyPruner.  

Applicability: 
-
1) If all inputs are proteins, it should work in any organisms, including prokaryotes and eukaryotes.  
2) If some of the inputs are transcripts, the genetic_code in TransDecoder should be adjusted (-g Universal). If the parameter is given, it will adopt TransDecoder in predicting coding potential in transcripts.  
3) If some of the inputs are genomic sequences or miniprot-based transcripts, you should add the database for alignment (-d database). Otherwise, it will raise exceptions and exit.    
          

Installation 
-
Pre-installation:

mamba (highly suggested) or conda. Link: https://github.com/conda-forge/miniforge#mambaforge

    git clone https://github.com/ylify/VEHoP.git #or download via release 
    cd VEHoP
    mamba env create --name phylogenomics -f environment.yml  
      #Once finished, a new environment named phylogenomics will be created, with most dependencies installed. 
    mamba activate phylogenomics
      #(if mamba is not installed in your system, use conda)  
Note: I have integrated most of many software and packages into the pipeline. Of them, only HmmCleaner.pl could not be configured by conda/mamba.   

Installation of HmmCleaner.pl (cd VEHoP): 
    
    chmod +x ./dependencies/cpanm 
    cpanm Bio::MUST::Apps::HmmCleaner 
      #(This step might take ~20 minutes; be patient; this installation always fails; no worried about that)
    ./dependencies/cpanm Bio::MUST::Apps::HmmCleaner --force 
      #(Try HmmCleaner.pl to check whether it was executable without errors. If errors, it will not produce results) 
As HmmCleaner.pl is unnecessary, and the installation cannot always be finished properly, I write it as the optional step in this pipeline. If such a file is not found or is not executable in your system or environment, it will automatically skip. You don't have to do anything.  
If you insist on installing it, please see the guidelines at https://metacpan.org/release/ARNODF/Bio-MUST-Apps-HmmCleaner-0.180750/source/INSTALL  

Dependency check:
-
    python3 VEHoP.py -h     
    It should be well resolved in dependency if the output includes help information.

Usage
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
                  Genetic code for protein prediction from transcripts, which might be different with phylum, please check by 
                  "TransDecoder.LongOrfs -h" (If the parameter is given, it will adopt TransDecoder to predict coding potential 
                  in transcripts. Optional if only proteins and genomic sequences as inputs; Required if transcripts existed in inputs, default: Universal) 
          -d DATABASE, --database DATABASE
                  Proteins sequences for homolog prediction from genomic sequences, it is suggested as proteins from its/their close
                  relatives (three organisms from the same genus, family, order, class, or phylum are suggested, from public data) 
                  (Optional if proteins or transcripts as inputs; Required if genomic sequences existed in inputs; 
                  It must be provided with the absolute path)
                  (Database will not be included in the matrix and tree)
  
Input
-
a folder (must be in the working directory, default: raw) containing sequences. We define the rule of three sources with specific suffixes. 

    1) genomic fasta: species_name.genomic.fasta
    2) transcript: species_name.transcript.fasta
    3) proteins: species_name.pep.fasta  
    
Note: species_name should be identical to others, otherwise it will fail in the tree visualization. We would recommend that you name the input files according to the rules below.  

    1) genus_species.genomic/transcript/pep.fasta  
    2) If more than one input from the same species, try:
      genus_species_1.genomic/transcript/pep.fasta and genus_species_2.genomic/transcript/pep.fasta  
    3) To distinguish from assembly method or source, try:
      genus_species_megahit.genomic.fasta for genomic assembly via megahit (the purge duplicated contig is not required), 
      genus_species_trinity.transcript.fasta for transcript assembly via Trinity (the selection of the longest isoform in gene is not required.)
      
Output
-
-   $PREFIX.$NUMBER-OF-INPUTS__$OCCCUPANCY.FastTree.full.tre (FastTreeMP -slow -gamma)
-   $PREFIX.$NUMBER-OF-INPUTS__$OCCCUPANCY.IQTREE2.full.tre (iqtree2 -m MFP)
-   homolog-phylogenomics.$PREFIX.$NUMBER-OF-INPUTS__$OCCCUPANC.$RUN-Day.log (running log of VEHoP)
-   miniprot/: the result of homolog-inference via miniprot, including gene feature files (gff and gff3) and predicted amino-acid sequences (pep.fasta)
-   transdecoder/: the result of homolog-inference via TransDecoder (default output), including the predicted amino-acid sequences (pep.fasta)
-   cd-hit/: the result of the non-redundant amino-acid sequences (from miniprot or TransDecoder) via CD-Hit (cut-off: 0.85)    
-   $PREFIX.$NUMBER-OF-INPUTS.orthofinder (the result of name-formatted amino-acid sequences (required in Phylopypruner) and the corresponding change log, and OrthoFinder. It also contains a Fullname_abbr.txt that records the formatted name and the original species name.)
-   $PREFIX.$NUMBER-OF-INPUTS__$OCCCUPANCY.Phylogenomics/ (the result of phylogenomic processes, including taxonomy occupancy, alignment, trimming, PhyloPyPruner, etc.)
    -   OG*.fa and OG*.tre (input files in PhyloPyPruner)
    -   01.backup_all_OGs/ (the top 40000 OGs)
    -   rejected_few_taxa_1/ (OGs with low taxonomy sampling)
    -   check_occupancy_1st.checkpoint.ok
    -   02.backup_preUniqHaplo/ (OGs before UniqHaplo)
    -   uniqHaplo.checkpoint.ok
    -   03.backup_alignments/ (OGs before mafft aligning)
    -   Mafft.checkpoint.ok
    -   04.back_pre_HmmCleaner/ (OGs before HmmCleaner processing)
    -   HmmCleaner.pl.checkpoint.ok
    -   05.backup_pre-trimal/ (OGs before trimal trimming)
    -   trimal.checkpoint.ok
    -   06.backup_pre-BMGE/ (OGs before BMGE processing)
    -   BMGE.checkpoint.ok
    -   07.back_pre_AlignmentCompare/ (OGs before AlignmentCompare processing)
    -   AlignmentCompare.checkpoint.ok
    -   08.backup_check_occupancy_2nd/ (OGs before occupancy check)
    -   rejected_few_taxa_2/ (OGs with low taxonomy sampling)
    -   check_occupancy_2nd.checkpoint.ok
    -   phylopypruner_output/
        -    All partitions (Folder: filtered)
        -    supermatrix.new.fas (concatenated matrix)
        -    partition_data.new.txt
        -    $PREFIX.$NUMBER-OF-INPUTS__$OCCCUPANCY.FastTree.full.tre (FastTreeMP -slow -gamma)
        -    $PREFIX.$NUMBER-OF-INPUTS__$OCCCUPANCY.IQTREE2.full.tre (iqtree2 -m MFP)
        -    run_ASTRAL_FastTreeMP.sh (gene trees implimented by FastTreeMP)
        -    run_ASTRAL_IQTREE2.sh (gene trees implimented by IQTREE2 with the best model)
        -    run_phylobayes.2500000.sh and $PREFIX.$NUMBER-OF-INPUTS__$OCCCUPANCY.2500000.fa
        -    run_phylobayes.5000000.sh and $PREFIX.$NUMBER-OF-INPUTS__$OCCCUPANCY.5000000.fa
        -    plots/ (heatmap of occupancy)
        -    output_alignments/ (each OG that was pseudo single-copy via PhyloPyPruner, but there might be some duplicated OGs that will raise errors in iqtree2)
        -    filtered/ (the deduplicated OGs from output_alignments)
        -    phylopypruner.log (PhyloPyPruner running log)
        -    partition_data.new.txt* (iqtree2 results)  

Tips in running
-
    
1) The script will check the existence of intermedia files in homolog-inference from genomic or transcriptional profiles. It will skip the step if the same input and database in the same working directory as before. If you want to run in a new directory with the same input, you could make soft links or copy them so that the pipeline will take less time. (We are considering to deposite this file in the package directory in the next version.)
2) Of course, more predicted proteins if more database. The consuming time in miniprot step will increase with the size of databases. We would recommend two or three high-quality proteomes based on our tests.
3) Please run in the same working directory if you would like to test the performance of difference taxonomy occupancy with the same input and database. It will skip the OrthoFinder step. 

Example
-
    working_directory: /home/yunlongli/Software/VEHoP/test
    database: /home/yunlongli/mollusca_three.pep.fasta 
    Time: 2023-12-25
    Command: python3 /home/yunlongli/Software/VEHoP/VEHoP.py -i test -t 40 -m 10 -p mollusca -d /home/yunlongli/mollusca_three.pep.fasta
    log_file: /home/yunlongli/Software/VEHoP/test/homolog-phylogenomics.mollusca.40__0.25.2023-12-25.log
    Result_directory: /home/yunlongli/VEHoP/test/mollusca.40__0.25.Phylogenomics/phylopyruner/
      /home/yunlongli/Software/VEHoP/test/mollusca.40__0.25.IQTREE2.full.tre
      /home/yunlongli/Software/VEHoP/test/mollusca.40__0.25.FastTree.full.tre
  
Publication
-
VEHoP: A Versatile, Easy-to-use, and Homology-based Phylogenomic pipeline accommodating diverse sequences    
Yunlong Li, Xu Liu, Chong Chen, Jian-Wen Qiu, Kevin Kocot, Jin Sun    
bioRxiv 2024.07.24.604968; doi: https://doi.org/10.1101/2024.07.24.604968    
     

Remark
-
If you have any questions, feel free to post an issue or email to ylify@connenct.ust.hk  
  
Please cite the integrated software (below) in this pipeline if you will include this pipeline, with doi or website listed. 
Since not all dependencies are included in your study, you could check the used ones in the log_file.  
Bioconda: https://doi.org/10.1038/s41592-018-0046-7  
General shell pipeline: https://doi.org/10.1093/sysbio/syw079  
AlignmentCompare: https://github.com/DamienWaits/Alignment_Compare.git  
BMGE: https://doi.org/10.1186/1471-2148-10-210  
cd-hit: https://doi.org/10.1093/bioinformatics/bts565  
FastTree: https://doi.org/10.1371/journal.pone.0009490  
HmmCleaner: https://doi.org/10.1186/s12862-019-1350-2  
IQ-TREE 2: https://doi.org/10.1093/molbev/msaa015  
miniprot: https://doi.org/10.1093/bioinformatics/btad014   
OrthoFinder: https://doi.org/10.1186/s13059-019-1832-y   
TransDecoder: https://github.com/TransDecoder/TransDecoder.git  
uniqHaplo: http://raven.wrrb.uaf.edu/~ntakebay/teaching/programming/perl-scripts/uniqHaplo.pl
