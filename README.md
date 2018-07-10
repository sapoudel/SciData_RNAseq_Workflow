# RNAseq_workflow
This pipeline will guide you through the analysis of RNAseq data, from FASTQ files to finding differentially expressed genes. Please read the [RNASeq Wiki](https://github.com/SBRG/SBRG_wiki/wiki/RNA-seq-Processing) before starting your analysis.

## Requirements:
1. [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
1. [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) or [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
1. [Samtools](http://www.htslib.org/)
1. [HTSeq-count](http://www-huber.embl.de/HTSeq/doc/install.html#install) (if using count_reads_python)
1. FASTA and GFF/Genbank file for your organism
1. [Jupyter](http://jupyter.org/install.html) with [R kernel](https://irkernel.github.io/)

Note: If installing Bowtie on Mac, you may need to install [brew](https://brew.sh/) and run `brew install tbb`

Alternatively, see Anand for access to a temporary prebuilt virtual machine with everything pre-installed.

## Pipeline:
1. Setup Organism ([0_setup_organism.ipynb](https://github.com/SBRG/RNAseq_workflow/blob/master/0_setup_organism.ipynb))
    * Creates Bowtie index and GFF files in the `ref/` directory
    * This step is not required for _E. coli_ K-12 MG1655
1. Align reads ([1_align_reads.ipynb](https://github.com/SBRG/RNAseq_workflow/blob/master/1_align_reads.ipynb))
    * Prepare your reads by creating `raw_files.csv`
    * `raw_files.csv` does not need absolute paths, as long as all the fastq files are in a common directory
    * Supports alignment to multiple genomes within a run (using the organism column)
    * Organism column must contain the same ID used in 0_setup_organism
    
| sample_id | R1  | R2  | organism |
|:----------:|:---:|:---:|:--------:|
| wt_fe2_1	| WT-Fe2-1_S1_L001_R1_001.fastq.trunc.gz | WT-Fe2-1_S1_L001_R2_001.fastq.trunc.gz | MG1655 |
|	wt_fe2_2	| WT-FE2-2_S2_L001_R1_001.fastq.trunc.gz | WT-FE2-2_S2_L001_R2_001.fastq.trunc.gz | MG1655 |

2. Count reads ([2_count_reads.ipynb](https://github.com/SBRG/RNAseq_workflow/blob/master/2_count_reads.ipynb))
    * This notebook requires the R kernel
3. Analysis ([3_analysis.ipynb](https://github.com/SBRG/RNAseq_workflow/blob/master/3_analysis.ipynb))
    * Contains workflows for:
        * Generating Transcripts Per Million (TPM)
        * Computing Principal Component Analysis (PCA)
        * Identifying differentially expressed genes (DEGs)

