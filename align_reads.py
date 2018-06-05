import os,warnings,subprocess,gzip,shutil,re

def gunzip(gz,out_dir):
    basename = os.path.split(gz)[1][:-3]
    result = os.path.join(out_dir,basename)
    with gzip.open(gz,'rb') as f:
        with open(result,'w') as f2:
            f2.write(f.read())
    return result

def get_alignment_score(out_dir,name,aligner):
    filename = os.path.join(out_dir,name+'_bowtie_output.txt')
    with open(filename,'r') as f:
        if aligner == 'bowtie':
            result = f.readlines()[1]
        else:
            result = f.readlines()[-1]
    if aligner == 'bowtie':
        match = re.search('\([\d\.]*%\)',result)
        return float(result[match.start()+1:match.end()-2])
    else:
        match = re.search('[\d{2,3}.]*%',result)
        return float(result[match.start():match.end()-1])

def align_reads(name,R1,R2,bt_index,out_dir,aligner='bowtie',cores=1,
                insertsize=1000,force=False,verbose=False):
    '''
    Align RNAseq FASTQ files using bowtie or bowtie2 to a reference genome.
    Returns the final BAM file location ('/<out_dir>/<name>.bam') and the
    bowtie alignment score

    name: str
        Unique name of sample
    R1: str
        Absolute location of R1 fastq file
    R2: str
        Absolute location of R2 fastq file
    bt_index: str
        Location of bowtie index to use for alignment
    out_dir: str
        Output directory
    aligner: 'bowtie' or 'bowtie2' (default 'bowtie')
        Designate which aligner to use
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
        print 'Processing %s'%name

    if not os.path.isdir(out_dir):
        warnings.warn('Creating output directory %s'%out_dir)
        os.makedirs(out_dir)

    # Quit if output file already exists
    if os.path.isfile(os.path.join(out_dir,name+'.bam')) and not force:
        score = get_alignment_score(out_dir,name,aligner)
        return os.path.join(out_dir,name+'.bam'),score

    # Check that all files exist
    if not os.path.isfile(R1):
        raise ValueError('File does not exist: %s'%R1)
    if not os.path.isfile(R2):
        raise ValueError('File does not exist: %s'%R2)
    if R1 == R2:
        raise ValueError('R1 and R2 files are identical')
    if not os.path.isfile(bt_index+'.1.ebwt') and \
       not os.path.isfile(bt_index+'.1.bt2'):
        raise ValueError('Bowtie index does not exist: %s'%bt_index)

    ### Unzip fastq files ###
    tmp_dir = os.path.join(out_dir,'tmp')
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)

    r1_files = []
    r2_files = []
    # In case the fastq files are split
    for fastq in R1.split(','):
        if fastq.endswith('.gz') and aligner=='bowtie':
            if verbose:
                print 'Unzipping file: %s'%fastq
            
            r1_files.append(gunzip(fastq,tmp_dir))
        else:
            r1_files.append(fastq)

    for fastq in R2.split(','):
        if fastq.endswith('.gz') and aligner=='bowtie':
            if verbose:
                print 'Unzipping file: %s'%fastq
            r2_files.append(gunzip(fastq,tmp_dir))
        else:
            r2_files.append(fastq)

    ### Run Bowtie Aligner ###

    
    bowtie_out = os.path.join(tmp_dir,name+'.sam')
    bowtie_err = os.path.join(out_dir,name+'_bowtie_output.txt')

    if aligner == 'bowtie':
        cmd = [aligner,'-X',str(insertsize),'-n','2','-p',str(cores),'-3','3',
               '-S','-1',','.join(r1_files),'-2',','.join(r2_files),bt_index]
    elif aligner == 'bowtie2':
        cmd = [aligner,'-X',str(insertsize),'-N','1','-p',str(cores),'-3','3',
               '-1',','.join(r1_files),'-2',','.join(r2_files),
               '-x',bt_index]
    else:
        raise ValueError('Aligner must be "bowtie" or "bowtie2"')


    if verbose:
        print 'Running bowtie: ' + ' '.join(cmd)

    with open(bowtie_out,'w') as out:
        with open(bowtie_err,'w') as err:
            subprocess.call(cmd,stdout=out,stderr=err)

    ### Post-process files ###

    unsorted_bam = os.path.join(tmp_dir,name+'.unsorted.bam')
    sorted_bam = os.path.join(out_dir,name+'.bam')

    sam2bam = ['samtools','view','-bS',bowtie_out,'-o',unsorted_bam]
    samsort = ['samtools','sort',unsorted_bam,'-o',sorted_bam]
    if verbose:
        print 'Converting to BAM: ' + ' '.join(sam2bam)
    subprocess.call(sam2bam)
    if verbose:
        print 'Sorting BAM file: ' + ' '.join(samsort)
    subprocess.call(samsort)

    ### Clear all non-bam files ###
    if verbose:
        print 'Cleaning up...'
    shutil.rmtree(tmp_dir)

    ### Find alignment score ###
    score = get_alignment_score(out_dir,name,aligner)

    ### Add BAM file and alignment value to replicate ###
    return sorted_bam,score

def build_index(sequence,bt_index,aligner='bowtie'):

    # Check that sequence exists
    if not os.path.isfile(sequence):
        raise ValueError('File does not exist: %s'%sequence)

    if aligner == 'bowtie':
        cmd = ['bowtie-build','-f',sequence,bt_index]
    elif aligner == 'bowtie2':
        cmd = ['bowtie2-build','-f',sequence,bt_index]
    else:
        raise ValueError('Aligner must be "bowtie" or "bowtie2"')

    subprocess.call(cmd)

    return

def gb2gff(sequence,genbank):

    from BCBio import GFF
    from Bio import SeqIO
    import pandas as pd
    import csv

    # Check that files exist
    if not os.path.isfile(sequence):
        raise ValueError('File does not exist: %s'%sequence)
    if not os.path.isfile(genbank):
        raise ValueError('File does not exist: %s'%genbank)

    out_dir = os.path.split(genbank)[0]

    out_file = os.path.splitext(genbank)[0]+'.gff'
        
    with open(sequence,'r') as f:
        header = f.readline()

    seqname = header[1:re.search('\s',header).start()]
    
    lines = []

    # Open genbank file
    with open(genbank,'r') as gb_handle:
        # Open record in genbank file
        for rec in SeqIO.parse(gb_handle, "genbank"):
            # Parse through each feature in genbank record
            for feature in rec.features:
                # Only grab info if feature is a CDS
                if feature.type == 'CDS':
                    if len(feature.location.parts) == 1:

                        # Get gene info
                        gene = feature.qualifiers['gene'][0]
                        try:
                            locus_tag = feature.qualifiers['locus_tag'][0]
                        except:
                            locus_tag = gene
                        start = feature.location.start.position
                        end = feature.location.end.position
                        if feature.location.strand == 1:
                            strand = '+'
                        else:
                            strand = '-'

                        # Now append this line to the growing GFF list
                        attr = 'gene_id "%s"; transcript_id "%s"; gene_name "%s";'%(locus_tag,locus_tag,gene)
                        lines.append([seqname,'feature','exon',start,end,
                                      '.',strand,'.',attr])

                    # If gene is split between two parts
                    else:
                        i = 1
                        for part in feature.location.parts:
                            # Get gene info
                            gene = feature.qualifiers['gene'][0]
                            try:
                                locus_tag = feature.qualifiers['locus_tag'][0] + '_' + str(i)
                            except:
                                locus_tag = gene + '_' + str(i)
                            start = part.start.position
                            end = part.end.position
                            if part.strand == 1:
                                strand = '+'
                            else:
                                strand = '-'

                            # Now append this line to the growing GFF list
                            attr = 'gene_id "%s"; transcript_id "%s"; gene_name "%s";'%(locus_tag,locus_tag,gene)
                            lines.append([seqname,'feature','exon',start,end,
                                          '.',strand,'.',attr])
                            i+=1

    DF_gff = pd.DataFrame(lines).sort_values(by=3,ascending=1)
    
    DF_gff.to_csv(out_file,sep='\t',quoting=csv.QUOTE_NONE,
                  index=False,header=False)
    return

    
