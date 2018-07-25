#/usr/bin/python
#PBS -l nodes=1:ppn=16
#PBS -l walltime=60:00:00
#PBS -l mem=600gb
#cd $PBS_O_WORKDIR

module load conda

#This pipeline aligne and reads against the reference file
#and generate the high quality variant call at the end. 
#In order to run this pipeline on Iridis you need to specify the
#PATH to bwa, GATK, PICHARD and SAMTOOLS.

#The whole genome data usually arives as multiple FastQ file;
#You can either concatanetate all the forward and all the reverse
#reads and run the pair of the concataneted FASTQ file through this pipeline
#or run each pair from each flowcell lane individually and then concatante the 
#bam file for the downstream processing.



import subprocess
import tempfile

LIB_PATH = '/data'

BWA_PATH = '{}/bwa-0.7.15'.format(LIB_PATH)
bwa = '{}/bwa'.format(BWA_PATH)

SAMTOOLS_PATH = '{}/samtools-1.2'.format(LIB_PATH)
samtools = '{}/samtools'.format(SAMTOOLS_PATH)

PICARD_PATH = '{}/picard-tools-1.119'.format(LIB_PATH)
create_seq_dict = '{}/CreateSequenceDictionary.jar'.format(PICARD_PATH)
mark_duplicates = '{}/MarkDuplicates.jar'.format(PICARD_PATH)
add_or_replace_read_groups = '{}/AddOrReplaceReadGroups.jar'.format(PICARD_PATH)

GATK_PATH = '{}/gatk-3.7'.format(LIB_PATH)
genome_analysis_tk = '{}/GenomeAnalysisTK.jar'.format(GATK_PATH)


DATA_PATH = '{}/input'.format(LIB_PATH)

FILE_PREFIX = 'chr19'
FA_FILE = '{}/{}.fa'.format(DATA_PATH, FILE_PREFIX)


# ------------ BWA ------------------------------------------------
subprocess.call([bwa, 'index', '-a', 'bwtsw', FA_FILE])
# chr19.fa.bwt
# chr19.fa.pac
# chr19.fa.ann
# chr19.fa.amb
# chr19.fa.sa

# ------------ SAMTOOLS -------------------------------------------
subprocess.call([samtools, 'faidx', FA_FILE])
# chr19.fa.fai

# ------------ CREATE SEQUENCE DICTIONARY -------------------------
DICT_FILE = '{}/{}.dict'.format(DATA_PATH, FILE_PREFIX)    # output file
subprocess.call([
    'java', '-jar',
    create_seq_dict,
    'R={}'.format(FA_FILE),
    'O={}'.format(DICT_FILE)])

# ------------ BWA-MEM --------------------------------------------
# NOTE: command for directing output of .sam to produce .bam
#     bwa-0.7.15/bwa mem -M -t 4 -aM -R
#     "@RG\tID:abc1\tSM:bwa2\tPL:illumia\tLB:unitabc1\tPU:Illuninate"
#     input/chr19.fa input/reads.fastq  | samtools view -bS - >
#     direct_output_using_term.bam

FASTQ_FILE = '{}/reads.fastq'.format(DATA_PATH)
dot_sam_file = '{}/reads.sam'.format(DATA_PATH)
dot_bam_file = '{}/reads.bam'.format(DATA_PATH)
sam_generator = subprocess.Popen([
    bwa,
    'mem',
    '-M',
    '-t',
    '4',
    '-aM',
    '-R',
    "@RG\tID:abc1\tSM:bwa2\tPL:illumia\tLB:unitabc1\tPU:Illuninate",
    FA_FILE,
    FASTQ_FILE],
    stdout=subprocess.PIPE)

# NOTE: The output should be directly into samtools to produce dot_bam_file.
# BUT samtools only accepts file path as input

# ------------ SAMTOOLS -------------------------------------------
# samtools view -s -h -b -t chr19.fa reads.sam -o reads.bam
bam_generator = subprocess.Popen([
    samtools,
    'view',
    '-bS', '-'],
    stdin=sam_generator.stdout,
    stdout=open(dot_bam_file, 'wb'))
stdout, _ = bam_generator.communicate()

# samtools sort reads.bam read.sorted.bam
dot_sorted_bam_file = '{}/reads.sorted'.format(DATA_PATH)
subprocess.call([
    samtools,
    'sort',
    dot_bam_file,
    dot_sorted_bam_file],
    stdout=subprocess.PIPE)

# extension .bam will be automaticall added to the output file
dot_sorted_bam_file = '{}.bam'.format(dot_sorted_bam_file)

# ------------ PICARD -------------------------------------------
sorted_de_dup_bam_file = '{}/reads.sortedDeDup.bam'.format(DATA_PATH)
subprocess.call([
    'java', '-jar',
    mark_duplicates,
    'INPUT={}'.format(dot_sorted_bam_file),
    'MAX_RECORDS_IN_RAM=2000',
    'REMOVE_DUPLICATES=false',
    'VALIDATION_STRINGENCY=SILENT',
    'ASSUME_SORTED=true',
    'METRICS_FILE={}/output.dups'.format(DATA_PATH),
    'OUTPUT={}'.format(sorted_de_dup_bam_file)],
    stdout=subprocess.PIPE)

sorted_de_dup_add_group_bam_file = '{}/reads.sortedDeDup_addGroup.bam'.format(DATA_PATH)
subprocess.call([
    'java', '-jar',
    add_or_replace_read_groups,
    'I={}/reads.sortedDeDup.bam'.format(DATA_PATH),
    'O={}'.format(sorted_de_dup_add_group_bam_file),
    'SORT_ORDER=coordinate',
    'CREATE_INDEX=true',
    'RGPL=illumina',
    'RGID=reads1',
    'RGSM=sample1',
    'RGLB=S1bar',
    'RGPU=pu1',
    'VALIDATION_STRINGENCY=LENIENT'],
    stdout=subprocess.PIPE)

# ------------ GATK -------------------------------------------
# ## Realigned
vcf_file = '{}/dbsnp_138.hg19.excluding_sites_after_129.vcf'.format(DATA_PATH)
sites_vcf_file = '{}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'.format(DATA_PATH)
output_interval_log = '{}/output.intervals.log'.format(DATA_PATH)
output_realigner_interval = '{}/output.ForIndelRealigner_reads.intervals'.format(DATA_PATH)
subprocess.call([
    'java', '-jar',
    genome_analysis_tk,
    '-T', 'RealignerTargetCreator',
    '-R', '{}'.format(FA_FILE),
    '-I', '{}'.format(sorted_de_dup_add_group_bam_file),
    '--known', '{}'.format(vcf_file),
    '--known', '{}'.format(sites_vcf_file),
    '-log', '{}'.format(output_interval_log),
    '-o', '{}'.format(output_realigner_interval)],
    stdout=subprocess.PIPE)

sorted_de_dup_add_group_realigned_bam_file = '{}/reads.sortedDeDup_addGroup_realigned.bam'.format(DATA_PATH)
subprocess.call([
    'java', '-jar',
    genome_analysis_tk,
    '-T',
    'IndelRealigner',
    '-R',
    FA_FILE,
    '-I',
    sorted_de_dup_add_group_bam_file,
    '-targetIntervals',
    output_realigner_interval,
    '-known',
    vcf_file,
    '-known',
    sites_vcf_file,
    '-o',
    sorted_de_dup_add_group_realigned_bam_file,
    '--filter_bases_not_stored'],
    stdout=subprocess.PIPE)

# ## Recalibration
recal_data_reads_table = '{}/recal_data_reads.table'.format(DATA_PATH)
subprocess.call([
    'java', '-jar',
    genome_analysis_tk,
    '-T',
    'BaseRecalibrator',
    '-R',
    FA_FILE,
    '-I',
    sorted_de_dup_add_group_realigned_bam_file,
    '-knownSites',
    vcf_file,
    '-o',
    recal_data_reads_table],
    stdout=subprocess.PIPE)

sorted_de_dup_add_group_re_aligned_recalib_bam_file = '{}/reads.sortedDeDup_addGroup_realigned_recalibrated.bam'.format(DATA_PATH)
subprocess.call([
    'sudo', 'java', '-jar',
    genome_analysis_tk,
    '-T',
    'PrintReads',
    '-R',
    FA_FILE,
    '-I',
    sorted_de_dup_add_group_realigned_bam_file,
    '-BQSR',
    recal_data_reads_table,
    '-o',
    sorted_de_dup_add_group_re_aligned_recalib_bam_file],
    stdout=subprocess.PIPE)

# ## Variant caller
sorted_de_dup_add_group_re_aligned_recalib_htc_vcf_file = '{}/reads.sortedDeDup_addGroup_realigned_recalibrated_haplotypecaller.vcf'.format(DATA_PATH)
subprocess.call([
    'java', '-jar',
    genome_analysis_tk,
    '-T',
    'HaplotypeCaller',
    '-R',
    FA_FILE,
    '-I',
    sorted_de_dup_add_group_re_aligned_recalib_bam_file,
    '--dbsnp',
    vcf_file,
    '-stand_call_conf',
    '30',
    '-o',
    sorted_de_dup_add_group_re_aligned_recalib_htc_vcf_file])
