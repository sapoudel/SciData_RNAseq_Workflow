# RNAseq_workflow
This pipeline will guide you through the analysis of RNAseq data, from FASTQ files to finding differentially expressed genes. Please read the [RNASeq Wiki](https://github.com/SBRG/SBRG_wiki/wiki/RNA-seq-Processing) before starting your analysis.

What you need:
1. [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
1. [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) or [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
1. [Samtools](http://www.htslib.org/)
1. [HTSeq-count](http://www-huber.embl.de/HTSeq/doc/install.html#install) (if using count_reads_python)
1. FASTA and GFF/Genbank file for your organism
1. [Jupyter](http://jupyter.org/install.html) with [R kernel](https://irkernel.github.io/)

Pipeline:
1. Align reads (align_reads.ipynb)
2. Count reads (count_reads.ipynb)
3. Analysis (analysis.ipynb)

Note: If installing Bowtie on Mac, you may need to install [brew](https://brew.sh/) and run `brew install tbb`
