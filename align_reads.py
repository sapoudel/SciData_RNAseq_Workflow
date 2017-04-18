import os,warnings,subprocess,gzip,shutil,re

def gunzip(gz,out_dir):
    basename = os.path.split(gz)[1][:-3]
    result = os.path.join(out_dir,basename)
    with gzip.open(gz,'rb') as f:
        with open(result,'w') as f2:
            f2.write(f.read())
    return result

def get_alignment_score(out_dir,name):
    filename = os.path.join(out_dir,name+'_bowtie_output.txt')
    with open(filename,'r') as f:
        result = f.readlines()[1]
    match = re.search('\([\d\.]*%\)',result)
    return float(result[match.start()+1:match.end()-2])

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
		score = get_alignment_score(out_dir,name)
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
	os.makedirs(tmp_dir)

	r1_files = []
	r2_files = []
	# In case the fastq files are split
	for fastq in R1.split(','):
		if fastq.endswith('.gz'):
			if verbose:
				print 'Unzipping file: %s'%fastq
			
			r1_files.append(gunzip(fastq,tmp_dir))
		else:
			r1_files.append(fastq)

	for fastq in R2.split(','):
		if fastq.endswith('.gz'):
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
	score = get_alignment_score(out_dir,name)

    ### Add BAM file and alignment value to replicate ###
	return sorted_bam,score
