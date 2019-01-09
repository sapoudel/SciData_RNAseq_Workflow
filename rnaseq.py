import os, subprocess, shutil, re


##########################
## Make Reference Files ##
##########################

def initialize_reference(org_id, fasta, gb, aligner):
    full_fasta = os.path.expanduser(fasta)
    full_gb = os.path.expanduser(gb)

    # Check that fasta and gb_full exist
    if not os.path.isfile(full_fasta):
        raise ValueError('File does not exist: %s' % full_fasta)
    if not os.path.isfile(full_gb):
        raise ValueError('File does not exist: %s' % full_gb)

    # Initialize reference directory
    org_dir = os.path.join('ref', org_id)
    if not os.path.isdir(org_dir):
        os.makedirs(org_dir)

    # Copy fasta and genbank file to organism directory
    new_fasta = os.path.join(org_dir, org_id + '.fasta')
    new_gb = os.path.join(org_dir, org_id + '.gb')
    shutil.copy(full_fasta, new_fasta)
    shutil.copy(full_gb, new_gb)

    # Build Bowtie Index
    bt_index = os.path.join(org_dir, org_id)
    build_index(new_fasta, bt_index, aligner=aligner)

    # Make GFF file
    print('Making GFF file...')
    gff = gb2gff(new_fasta, new_gb)

    # Add aligner info
    with open(os.path.join(org_dir, 'info.txt'), 'w') as f:
        f.write('\n'.join([aligner, bt_index, gff]))

    # Print verbose output
    print('')
    print('Created Directory: ' + org_dir)
    for f in sorted(os.listdir(org_dir)):
        print('-> ' + f)


def build_index(sequence, bt_index, aligner='bowtie'):
    # Check that aligner is properly installed
    try:
        subprocess.check_output([aligner, '--version'])
    except:
        raise ValueError('Aligner not installed correctly: ' + aligner)

    # Check that sequence exists
    if not os.path.isfile(sequence):
        raise ValueError('File does not exist: %s' % sequence)

    if aligner == 'bowtie':
        print('Building bowtie index: bowtie-build -f ' + sequence + ' ' + bt_index)
        cmd = ['bowtie-build', '-f', sequence, bt_index]
    elif aligner == 'bowtie2':
        print('Building bowtie index: bowtie2-build -f ' + sequence + ' ' + bt_index)
        cmd = ['bowtie2-build', '-f', sequence, bt_index]
    else:
        raise ValueError('Aligner must be "bowtie" or "bowtie2"')

    subprocess.call(cmd)

    return
        
def gb2gff(sequence,genbank,id_tag='locus_tag'):

    from Bio import SeqIO
    import pandas as pd
    import csv

    # Check that files exist
    if not os.path.isfile(sequence):
        raise ValueError('File does not exist: %s' % sequence)
    if not os.path.isfile(genbank):
        raise ValueError('File does not exist: %s' % genbank)

    out_dir = os.path.split(genbank)[0]

    out_file = os.path.splitext(genbank)[0] + '.gff'

    fasta_ids = [record.id for record in SeqIO.parse(sequence, 'fasta')]

    lines = []

    # Open genbank file
    with open(genbank, 'r') as gb_handle:
        # Open record in genbank file
        for rec in SeqIO.parse(gb_handle, "genbank"):
            # Check that record is also in fasta
            if rec.id not in fasta_ids:
                raise ValueError('Fasta file is missing record: ' + rec.id)
            # Parse through each feature in genbank record
            for feature in rec.features:
                # Only grab info if feature is a CDS
                if feature.type == 'CDS':
                    if len(feature.location.parts) == 1:
                        gene_id = feature.qualifiers[id_tag][0]
                        start = feature.location.start.position
                        end = feature.location.end.position
                        if feature.location.strand == 1:
                            strand = '+'
                        else:
                            strand = '-'

                        # Now append this line to the growing GFF list
                        attr = 'gene_id "%s"; transcript_id "%s"' % (gene_id, gene_id)
                        lines.append([rec.id, 'feature', 'exon', start, end,
                                      '.', strand, '.', attr])

                    # If gene is split between two parts
                    else:
                        i = 1
                        for part in feature.location.parts:
                            # Get gene info
                            gene_id = feature.qualifiers[id_tag][0] + '_' + str(i)
                            start = part.start.position
                            end = part.end.position
                            if part.strand == 1:
                                strand = '+'
                            else:
                                strand = '-'

                            # Now append this line to the growing GFF list
                            attr = 'gene_id "%s"; transcript_id "%s"' % (gene_id, gene_id)
                            lines.append([rec.id, 'feature', 'exon', start, end,
                                          '.', strand, '.', attr])
                            i += 1

    DF_gff = pd.DataFrame(lines).sort_values(by=[0, 3], ascending=[1, 1])

    DF_gff.to_csv(out_file, sep='\t', quoting=csv.QUOTE_NONE,
                  index=False, header=False)
    return out_file

#################
## Align Reads ##
#################

def gunzip(gz, out_dir):
    basename = os.path.split(gz)[1][:-3]
    result = os.path.join(out_dir, basename)
    subprocess.call('gunzip < ' + gz + ' > ' + result, shell=True)
    return result


def get_alignment_score(out_dir, name, aligner):
    filename = os.path.join(out_dir, name + '_bowtie_output.txt')
    with open(filename, 'r') as f:
        if aligner == 'bowtie':
            result = f.readlines()[1]
        else:
            result = f.readlines()[-1]
    if aligner == 'bowtie':
        match = re.search('\([\d\.]*%\)', result)
        return float(result[match.start() + 1:match.end() - 2])
    else:
        match = re.search('[\d{2,3}.]*%', result)
        return float(result[match.start():match.end() - 1])


def align_reads(name, R1, R2, organism, in_dir, out_dir, cores=8,
                insertsize=1000, force=False, verbose=False):
    '''
    Align RNAseq FASTQ files using bowtie or bowtie2 to a reference genome.
    Returns the final BAM file location ('/<out_dir>/<name>.bam'), and other
    quality scores.

    name: str
        Unique name of sample
    R1: str
        Absolute location of R1 fastq file
    R2: str
        Absolute location of R2 fastq file
    organism: str
        Organism ID for alignment 
    in_dir: str
        Directory with raw files
    out_dir: str
        Output directory  
    insertsize: int (default 1000)
        Maximum distance between paired ends
    cores: int (default 1)
        Number of cores
    force: boolean (default False)
        Re-runs alignment even if BAM file already exists
    verbose: boolean (default False)
        Update user with current process

    '''

    if verbose:
        print 'Processing %s' % name

    # Set-up and quality check
    if not os.path.isdir(out_dir):
        if verbose:
            print('Creating output directory %s' % out_dir)
        os.makedirs(out_dir)

    # Check that organism directory exists
    org_dir = os.path.join('ref', organism)
    if not os.path.isdir(org_dir):
        raise ValueError('Reference not created for organism. See 0_setup_organism')

    # Get bowtie index and aligner
    with open(os.path.join(org_dir, 'info.txt'), 'r') as f:
        aligner, bt_index, gff = [line[:-1] for line in f.readlines()]

    # Quit if output file already exists
    if os.path.isfile(os.path.join(out_dir, name + '.bam')) and not force:
        score = get_alignment_score(out_dir, name, aligner)
        bam_stats = get_bam_stats(out_dir, name)
        return os.path.join(out_dir, name + '.bam'), score, bam_stats

    # Split R1/R2 files if necessary and get filename
    r1_files = [os.path.join(in_dir, f) for f in R1.split(';')]
    r2_files = [os.path.join(in_dir, f) for f in R2.split(';')]

    # Check that all files exist
    for f in r1_files + r2_files:
        if not os.path.isfile(f):
            raise ValueError('File does not exist: %s' % f)
    if R1 == R2:
        raise ValueError('R1 and R2 files are identical')

    # Unzip fastq files
    tmp_dir = os.path.join(out_dir, 'tmp')
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)

    final_r1 = []
    final_r2 = []
    # In case the fastq files are split
    for fastq in r1_files:
        if fastq.endswith('.gz') and aligner == 'bowtie':
            if verbose:
                print 'Unzipping file: %s' % fastq

            final_r1.append(gunzip(fastq, tmp_dir))
        else:
            final_r1.append(fastq)

    
    for fastq in r2_files:
        if fastq.endswith('.gz') and aligner == 'bowtie':
            if verbose:
                print 'Unzipping file: %s' % fastq
            final_r2.append(gunzip(fastq, tmp_dir))
        else:
            final_r2.append(fastq)

    # Run Bowtie Aligner

    bowtie_out = os.path.join(tmp_dir, name + '.sam')
    bowtie_err = os.path.join(out_dir, name + '_bowtie_output.txt')

    if aligner == 'bowtie':
        cmd = [aligner, '-X', str(insertsize), '-n', '2', '-p', str(cores), '-3', '3',
               '-S', '-1', ','.join(final_r1), '-2', ','.join(final_r2), bt_index]
    elif aligner == 'bowtie2':
        cmd = [aligner, '-X', str(insertsize), '-N', '1', '-p', str(cores), '-3', '3',
               '-1', ','.join(final_r1), '-2', ','.join(final_r2),
               '-x', bt_index]
    else:
        raise ValueError('Aligner must be "bowtie" or "bowtie2"')


    if verbose:
        print 'Running bowtie: ' + ' '.join(cmd)

    with open(bowtie_out, 'w') as out:
        with open(bowtie_err, 'w') as err:
            subprocess.call(cmd, stdout=out, stderr=err)

    # Post-process files

    unsorted_bam = os.path.join(tmp_dir, name + '.unsorted.bam')
    sorted_bam = os.path.join(out_dir, name + '.bam')

    sam2bam = ['samtools', 'view', '-b', bowtie_out, '-@', str(cores), '-o', unsorted_bam]
    samsort = ['samtools', 'sort', unsorted_bam, '-@', str(cores), '-o', sorted_bam]
    if verbose:
        print 'Converting to BAM: ' + ' '.join(sam2bam)
    subprocess.call(sam2bam)
    if verbose:
        print 'Sorting BAM file: ' + ' '.join(samsort)
    subprocess.call(samsort)
    
    # Move final BAM and bowtie output files to output directory
    shutil.move(sorted_bam,out_dir)
    shutil.move(bowtie_err,out_dir)

    # Clear all non-bam files ###
    if verbose:
        print 'Cleaning up...'
    shutil.rmtree(tmp_dir)

    # Find alignment score
    score = get_alignment_score(out_dir, name, aligner)
    bam_stats = get_bam_stats(out_dir, name)

    # Add BAM file and alignment value to replicate ###
    return sorted_bam, score, bam_stats


def get_bam_stats(out_dir, name):
    """
    :param bamfile:
    :return: mean phred quality, total reads, mapped reads
    """
    bamfile = os.path.join(out_dir, name + '.bam')
    call = ['samtools', 'stats', '-@', '2', bamfile]
    output = subprocess.check_output(call)

    qsearch = re.search('average quality:(.*)', output)
    qual = float(qsearch.group().split('\t')[-1])

    tsearch = re.search('raw total sequences:(.*)', output)
    total_seq = int(tsearch.group().split('\t')[-1])

    msearch = re.search('reads mapped:(.*)', output)
    mapped_seq = int(msearch.group().split('\t')[-1])

    return qual, total_seq, mapped_seq
