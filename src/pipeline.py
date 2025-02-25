#!/usr/bin/python
'''
author: Chao-Jung Wu <wu.chaojung@gmail.com>
2020-09-18
2021-04-05 refactor
2021-04-28 update
2022-04-10 update
2022-04-30 major refactor, app name = ArcticFox
2022-07-19 tag 0.4.21.03, this version complete igv snapshots.
2022-07-20 tag 0.4.21.04, this version complete results collection to an upper folder.
'''
from __future__ import print_function
import os
import arg
import re
import subprocess
import datetime
import time

Version_ArcticFox='0.4.21.04'
Author='Chao-Jung Wu <wu.chaojung@gmail.com>'




def reverse_complement(seq):
  """Returns a reversed string"""
  seq = seq[::-1]
  """Returns a complement DNA sequence"""
  complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
  seq_list = list(seq)
  seq_list = [complement_dict[base] for base in seq_list]
  return ''.join(seq_list)
  
args = arg.init_params ()
hpc = args.hpc
libname = args.libname
rep_input = (args.rep_input).rstrip('/') + '/'
ref_fasta_key = args.genome_keyword + '_stdchr'
ref_fasta = '../dbs/' + ref_fasta_key + '.fa'#args.ref_fasta
cargo_fasta = args.cargo_fasta #'ITR2ITR_bsg168.fa'
backbone_fasta = args.backbone_fasta #'backbone_bsg168.fa'
plasmids_fasta = args.plasmids_fasta
#####
barcode = args.barcode
AdapterYtop = args.AdapterYtop
AdapterYbottom = args.AdapterYbottom
linker5p = args.linker5p
linker3p = args.linker3p
aligner = args.aligner
#####
trim_poorSequencing_Nanopore = args.trim_poorSequencing_Nanopore
trimAdapters = True #args.trimAdapters
trim_poorAlignment = args.trim_poorAlignment
min_LambdaCt = args.min_LambdaCt
#min_lenRead = args.min_lenRead
#min_lenAlnChr = args.min_lenAlnChr

getHaploFasta = args.getHaploFasta
#####

barcodeRC = reverse_complement(barcode)
adapter5prime = AdapterYtop + linker5p + barcode
adapter3prime = barcodeRC   + linker3p + AdapterYbottom

rep_output = '../output/'+ libname + '/' ## force the output to be generate within ArcticFox. In the end, I should zip a copy and ship it to its intended storage destination.
logfile = rep_output + 'CJWAAV.log.txt'
#if not os.path.exists(rep_output): os.makedirs(rep_output)

print(args)


inlineLoad = ''
inlineLoad_cutadapt = ''
inlineLoad_samtool = ''
inlineLoad_fastqc = ''
inlineLoad_bwa = ''
inlineLoad_mm2 = ''
  
if hpc == "edge":
  #=Edge login node does not need any inlineLoad hpc modules. But scheduler version needs inlineLoad.
  inlineLoad = 'module load anaconda3/4.9.2.compchem;'
  inlineLoad_cutadapt = 'module load anaconda3/4.9.2.compchem;'
  inlineLoad_samtool = 'module load anaconda3/4.9.2.compchem;module load samtools/1.12-gcc-11.1.0;'
  inlineLoad_fastqc = 'module load fastqc/0.11.9-gcc-11.1.0;'
  inlineLoad_bwa = 'module load bwa/0.7.17-gcc-11.1.0;'
  inlineLoad_mm2 = 'module load minimap2/2.14-gcc-11.1.0;' # omit testing this.

if hpc == "cambridge":
  inlineLoad = 'module load miniconda3/4.3.21;'
  inlineLoad_cutadapt = 'module load pymods/cutadapt/1.18;'
  inlineLoad_samtool = 'module load samtools/1.10;'
  inlineLoad_fastqc = 'module load java/jre/1.8.0_66 fastqc/0.11.5;'
  inlineLoad_bwa = 'module load bwa/0.7.15;'
  inlineLoad_mm2 = 'module load minimap2/2.17;' # omit testing this.





def logArg (): #(args, rep_output):
  if not os.path.exists(rep_output): os.makedirs(rep_output)

  time_a = datetime.datetime.now()
  NowStr = time_a.strftime("%a %b %d, %Y. %H:%M:%S")

  with open (logfile, 'w') as fh:
    print(NowStr, file=fh)
    print('', file=fh)
    print('ArcticFox Version:', Version_ArcticFox, file=fh)
    print('Author:', Author, file=fh)
    print('', file=fh)
    for arg in sorted(vars(args)):
      print (arg + ':\t\t\t', getattr(args, arg), file=fh)
    print('', file=fh)
    
  return time_a
      
def logcmd (cmd):
  with open (logfile, 'a') as fh:
    print(cmd, file=fh)
    

      
def subprocess_tryexcept (cmd):
  logcmd (cmd)
  try:
    subprocess.check_call(cmd, shell=True)
  except Exception as e:
    print("An exception occurred")
    print(e)
    logcmd ("An exception occurred")
    logcmd (e)

def cat_ref ():
  #cargo_fasta = 'ITR2ITR_bsg168.fa' ###
  #backbone_fasta = 'backbone_bsg168.fa' ###
  rep = '../dbs/'
  key = ' ' + rep
  L = [cargo_fasta] + '02_AmpR.fa,03_KanR.fa,04_Rep2.fa,05_Cap9.fa,06_HPV18.fa,07_Ad5.fa,09_lambda.fa'.split(',') + [backbone_fasta]
  files = rep + key.join(L)
  if not os.path.exists(rep + cargo_fasta): 
    print('error: file does not exist:', rep + cargo_fasta)
    exit()
  if not os.path.exists(rep + backbone_fasta): 
    print('error: file does not exist:', rep + backbone_fasta)
    exit()
  cmd = 'cat ' + files + ' > ' + ref_fasta
  print(cmd);subprocess_tryexcept (cmd)
  
def cat_ref_hasBug ():
  #cargo_fasta = 'ITR2ITR_bsg168.fa' ###
  #backbone_fasta = 'backbone_bsg168.fa' ###
  rep = '../dbs/'
  key = ' ' + rep
  L = [cargo_fasta] + '02_AmpR.fa,03_KanR.fa,04_Rep2.fa,05_Cap9.fa,06_HPV18.fa,07_Ad5.fa,09_lambda.fa'.split(',') + [backbone_fasta]

  if not os.path.exists(rep + cargo_fasta): 
    print('error: file does not exist:', rep + cargo_fasta)
    exit()
  if not os.path.exists(rep + backbone_fasta): 
    print('error: file does not exist:', rep + backbone_fasta)
    backbone_fasta = '10_F1pUCoris.fa'
    print('Note: use generic cis backbone instead:', rep + backbone_fasta)
    #exit()
  files = rep + key.join(L)
  cmd = 'cat ' + files + ' > ' + ref_fasta
  print(cmd);subprocess_tryexcept (cmd) 
 
def merge_fastq_compress (rep_input, outfile):
  from os import walk
  filenames = next(walk(rep_input), (None, None, []))[2]
  _, file_extension = os.path.splitext(filenames[0])
  cat = 'zcat' if file_extension == '.gz' else 'cat'
  gz = '.gz' if file_extension == '.gz' else ''

  cmd = cat + ' ' + rep_input + '*.fastq' + gz + '  | gzip > ' + outfile
  print(cmd)
  subprocess_tryexcept (cmd)    
    
def NanoFilt (infile_fastq_gz):
  #= syntax: gunzip -c test.fastq.gz | NanoFilt -q 10 -l 100 | gzip > test.seqQ10.fastq.gz
  Q, L = '10', '100'
  outfile = re.sub('.fastq.gz', '.seqQ' + Q + '.fastq.gz', infile_fastq_gz)
  cmd = inlineLoad
  cmd += 'gunzip -c ' + infile_fastq_gz
  cmd += ' | NanoFilt -q ' + Q  + ' -l ' + L #' | NanoFilt -q 10 -l 100'
  cmd += ' | gzip > ' + outfile
  print(cmd)
  subprocess_tryexcept (cmd)
  head, tail = os.path.split(infile_fastq_gz)
  cmd = 'mv NanoFilt.log ' + head;os.system(cmd)
  return outfile  
    
def cutadapt_5prime (infile, adapter, outfile, stdout_cutadapt_file):
  #cutadapt -g "XGGCGTCTGCTTG;max_error_rate=0.25" -o test.delAdapter5p.fastq.gz test.fastq.gz > cutadapt_stdout.txt
  max_error_rate = 0.25 ## 0.1 to 0.5
  cmd = inlineLoad_cutadapt
  cmd += 'cutadapt -g "X' + adapter + ';max_error_rate=' + str(max_error_rate) + '" -o ' + outfile + ' ' + infile 
  cmd += ' > ' + stdout_cutadapt_file
  print(cmd)
  subprocess_tryexcept (cmd)

def cutadapt_3prime (infile, adapter, outfile, stdout_cutadapt_file):
  #cutadapt -a "TTCGGATTCTATCGX;max_error_rate=0.25" -o test.delAdapter3p.fastq.gz test.fastq.gz > cutadapt_stdout.txt
  max_error_rate = 0.25 ## 0.1 to 0.5
  cmd = inlineLoad_cutadapt
  cmd += 'cutadapt -a "' + adapter + 'X;max_error_rate=' + str(max_error_rate) + '" -o ' + outfile + ' ' + infile 
  cmd += ' > ' + stdout_cutadapt_file
  print(cmd)
  subprocess_tryexcept (cmd) 

def fastqc (infile, rep_output):
  cmd = inlineLoad_fastqc
  cmd += 'fastqc ' + infile + ' --outdir ' + rep_output
  cmd += ' >/dev/null 2>&1'
  print(cmd)
  subprocess_tryexcept (cmd)

def NanoPlot (infile, rep_output):
  ''' Since 210410, nanoplot's dependencies are not supported in miniconda3/4.3.21's imports 
      However, the outputs of NanoPlot are still generated. Able to locate them in the output folder.
  '''
  #= infile fastq file can be fastq.gz or fastq
  #= infile_type = ['fastq', 'bam']
  if infile.endswith('fastq') or infile.endswith('fq') or infile.endswith('fastq.gz'): infile_type = 'fastq'
  elif infile.endswith('bam'): infile_type = 'bam'
  else:
    msg = 'error: wrong filetype, can not generate nanoplot', infile
    print(msg)
    logcmd (msg)
    return
  inBasename = os.path.basename(infile)
  cmd = inlineLoad
  cmd += 'NanoPlot -t 8 --prefix ' + rep_output + inBasename + '.' 
  cmd += ' --' + infile_type + ' '
  cmd += infile
  cmd += ' >/dev/null 2>&1'
  print(cmd)
  subprocess_tryexcept (cmd)

def bwa_build_index (ref_fasta):
  tmpref_fasta = ref_fasta + '.fai'
  if os.path.isfile(tmpref_fasta): return
  print('=====================================')
  print('== building reference genome index ==')
  cmd = inlineLoad_bwa
  cmd += 'bwa index -a bwtsw ' + ref_fasta
  cmd += ' >/dev/null 2>&1'
  print(cmd)
  subprocess_tryexcept (cmd)
  make_fai (ref_fasta)
  print('=====================================')

def minimap2_build_index (ref_fasta):
  #=syntax: minimap2 -d target.mmi target.fa
  prefix = os.path.splitext(ref_fasta)[0] #= "/path/to/file.txt.txt.fa" => "/path/to/file.txt.txt"
  outfile = prefix + '.mmi' #= "/path/to/file.txt.txt" => "/path/to/file.txt.txt.mmi"
  if os.path.isfile(outfile): return
  print('=====================================')
  print('== building reference genome index ==')
  cmd = inlineLoad_mm2
  cmd += 'minimap2 -d ' + outfile + ' ' + ref_fasta
  print(cmd)
  subprocess_tryexcept (cmd)
  print('=====================================')

def minimap2 (libname, ref_fasta, library_file, outfile_sam):
  #* Long-read alignment with CIGAR: (manual)
  # minimap2 -a [-x preset] target.mmi query.fa > output.sam
  # minimap2 -c [-H] [-k kmer] [-w miniWinSize] [...] target.fa query.fa > output.paf
  ##
  # ddl => minimap2 -a -x map-ont --eqx -L -O 5,56 -E 4,1 -B 5 -z 400,50 -r 2k -Y  -R "@RG\tID:[libname]\tSM:[libname]" [ref_fasta] [nanofiltered library file] > [sam file]
  minimap2_build_index (ref_fasta)
  cmd = inlineLoad_mm2
  cmd += 'minimap2 -a -x map-ont --eqx -L -O 5,56 -E 4,1 -B 5 -z 400,50 -r 2k -Y  -R '
  cmd += '"@RG\tID:' + libname + '\tSM:' + libname + '" '
  cmd += ref_fasta + ' ' + library_file + ' > ' + outfile_sam
  #cmd += ' >/dev/null 2>&1' #= do not use this line, it causes error for the next steps for yet unknown reasons.
  print(cmd)
  subprocess_tryexcept (cmd)

def bwa_mem (infile, ref_fasta, outfile):
  """ bwa mem ref.fa in.fastq > out.sam """
  bwa_build_index (ref_fasta)
  cmd = inlineLoad_bwa
  cmd += 'bwa mem ' + ref_fasta + ' ' + infile
  cmd += ' > ' + outfile
  print(cmd)
  subprocess_tryexcept (cmd)

def index_sam_into_bam (infile):
  cmd = 'samtools index ' + infile + ' >/dev/null 2>&1'
  print(cmd)
  subprocess_tryexcept (cmd)

def sam2bam (infile):
  bam_out = re.sub('.sam', '.bam', infile)
  cmd = inlineLoad_samtool
  cmd += 'samtools view -bS ' + infile + ' > ' + bam_out
  print(cmd)
  subprocess_tryexcept (cmd)
  return bam_out

def bam2sam (infile):
  #= samtools view -h -o out.sam in.bam
  sam_out = re.sub('.bam', '.sam', infile)
  cmd = inlineLoad_samtool
  cmd += 'samtools view -h -o ' + sam_out + ' ' + infile
  print(cmd)
  subprocess_tryexcept (cmd)
  return sam_out

def sort_bam (infile):
  sortedBam_out = re.sub('.bam', '.sorted.bam', infile)
  cmd = inlineLoad_samtool
  cmd += 'samtools sort -o ' + sortedBam_out + ' ' + infile
  print(cmd)
  subprocess_tryexcept (cmd)
  cmd = 'mv ' + sortedBam_out + ' ' + infile
  subprocess_tryexcept (cmd)
  return infile

def sam2SortedBam (alignmentSam):
  alignmentBam = sam2bam (alignmentSam)
  alignmentSortedBam = sort_bam (alignmentBam)
  return alignmentSortedBam

def make_fai (ref_fasta):
  outfile = ref_fasta + '.fai'
  if os.path.exists(outfile): return outfile
  cmd = inlineLoad_samtool
  cmd += 'samtools faidx ' + ref_fasta
  print(cmd)
  subprocess_tryexcept (cmd)
  return outfile

def make_bai (sortedBam):
  cmd = inlineLoad_samtool
  cmd += 'samtools index -b ' + sortedBam + ' ' + sortedBam + '.bai'
  print(cmd)
  subprocess_tryexcept (cmd)
  return sortedBam + '.bai'

def idxstats (infile):
  outfile = re.sub('.bam', '.idxstats.tsv', infile)
  bai = make_bai (infile) #= idxstats needs bai file
  cmd = inlineLoad_samtool
  cmd += 'samtools idxstats ' + infile
  cmd += ' > ' + outfile
  print(cmd)
  subprocess_tryexcept (cmd)
  return outfile

def flagstat (infile):
  outfile = re.sub('.bam', '.flagstat.txt', infile)
  cmd = inlineLoad_samtool
  cmd += 'samtools flagstat ' + infile
  cmd += ' > ' + outfile
  print(cmd)
  subprocess_tryexcept (cmd)

def basic_stats (infile):
  """
  I don't think there is much to see in this report, but it shows error rate
  error rate = # mismatches / bases mapped 
  """
  outfile = re.sub('.bam', '.stats.txt', infile)
  cmd = inlineLoad_samtool
  cmd += 'samtools stats ' + infile
  cmd += ' > ' + outfile
  print(cmd)
  subprocess_tryexcept (cmd)

  infile = outfile
  out_prefix = re.sub('.stats.txt', '.stats-', infile)
  cmd = inlineLoad_samtool
  cmd += 'plot-bamstats -p ' + out_prefix + ' ' + infile
  cmd += ' >/dev/null 2>&1'
  print(cmd)
  subprocess_tryexcept (cmd)

def filter_mappingQuality (infile):
  """ 
  remove reads with mapping quality < 10 
  infile = bam_sorted_aln
  """
  Q = '10'
  outHQ = re.sub('.bam', '.mapQ' + Q + '.bam', infile)
  cmd = inlineLoad_samtool
  cmd += 'samtools view -bh -q ' + Q + ' ' + infile + ' > ' + outHQ
  print(cmd)
  subprocess_tryexcept (cmd)
  return outHQ

def bam2haplotype (infile, rep_haplo, fastqmerged, f_idx, infile2):
  relocated_infile = rep_haplo + os.path.basename(infile)
  cmd = 'cp ' + infile + ' ' + relocated_infile;print(cmd);subprocess_tryexcept (cmd)
  infile = relocated_infile
  cmd = inlineLoad
  cmd += 'python decodeSAMpositions_pysam.py ' + infile + ' ' + f_idx
  cmd += ' >/dev/null 2>&1'
  print(cmd)
  subprocess_tryexcept (cmd)
  infileRH = re.sub('.bam', '.read_haplotype.tsv', infile)
  cmd = inlineLoad
  cmd += 'python haplotype_breakdown.py ' + infileRH + ' ' + f_idx + ' ' + infile2
  print(cmd)
  subprocess_tryexcept (cmd)
  cmd = 'rm ' + relocated_infile;print(cmd)
  subprocess_tryexcept (cmd)
  return infileRH

def getLambdaCt (f_idx):
  lambdaCt = 0
  with open (f_idx, 'r') as fh: DATA = [x.rstrip('\n') for x in fh.readlines()]
  for i in DATA:
    if i.startswith ('lambda'):
      data = i.split('\t')
      lambdaCt = data[2]
  return lambdaCt

def getFastaFiles_for_allHaplotypes (infileRH, fastqmerged):
  cmd = inlineLoad
  cmd += 'python extractReadsByHaplotypes.py ' + infileRH + ' ' + fastqmerged 
  print(cmd)
  subprocess_tryexcept (cmd)

def mass_pct (df, filenameRef):
  outfile1 = re.sub('.describe_reads.tsv', '.describe_reads.sumtmp.tsv', filenameRef)
  df.drop('readID', axis=1, inplace=True)
  df.drop('haplotype', axis=1, inplace=True)
  df.rename(columns={'read_length': 'library_withLambda'}, inplace=True)
  df.loc['total_basecount'] = df.sum(numeric_only=True, axis=0) ###
  info = df.loc['total_basecount']
  info.to_csv(outfile1, sep='\t')
  
  infile = outfile1 #"mass.txt"
  outfile2 = re.sub('.describe_reads.sumtmp.tsv', '.describe_reads.mass_pct.tsv', outfile1)
  with open (infile, "r") as fh: 
    DATA = [x.rstrip('\n').split("\t") for x in fh.readlines()]
  num_library_withLambda, num_lambda = 0, 0
  for i in DATA:
    if i[0].startswith("library_withLambda"):
      num_library_withLambda = int(float(i[1]))
    if i[0].startswith("lambda"):
      num_lambda = int(float(i[1]))
  num_library_noLambda = num_library_withLambda - num_lambda
  DATA[0] = DATA[0] + ["%_noLambda"]
  for i in range(1, len(DATA)):
    DATA[i][1] = int(float(DATA[i][1]))
    num = round((    DATA[i][1] / num_library_noLambda ) * 100, 2)
    DATA[i].append(num)
  data = ['library_noLambda', int(num_library_noLambda), 100.00]
  DATA.append(data)
  
  with open (outfile2, 'w') as fh:
    for i in DATA:
      print('\t'.join([str(x) for x in i]), file=fh)
    line = '# For cell host module, percentage is based on num_library_withLambda because lambda is not analyzed.'
    print(line, file=fh)
  
  cmd = 'rm ' + outfile1; print(cmd)
  subprocess_tryexcept (cmd)
  
  
  
  
  

def describe (infile):
  ''' This function generates two tables in one txt file.
      First table is the percentiles statistics for each chromosome.
      Second talbe is the sum of bases aligned to each chromosome.
      
      @ The second table is not formatted properly. TO DO.
  ''' 
  import pandas as pd
  import numpy as np
  pd.options.display.float_format = '{:.1f}'.format #= suppress scientific notation
  outfile = re.sub('.read_haplotype.tsv', '.describe_reads.tsv', infile)
  df = pd.read_csv(infile, sep='\t')
  df.drop('aln_baseCount', axis=1, inplace=True) #= remove a col
  #if 'lambda' in df.columns.values.tolist(): df = df[df['lambda'] == 0]
  df.replace(0, np.nan, inplace=True)
  #= now describe() chr aln
  info = df.describe().apply(lambda s: s.apply('{0:.2f}'.format))
  info.to_csv(outfile, sep='\t')
  with open (outfile, 'a') as fho:
    print('# Count is readcount. Other lines are basecount. ie: count (readcount), mean (bp), std (bp), min (bp), 25% (bp), 50% (bp), 75% (bp), max (bp)', file=fho)
    print('# % is percentile.', file=fho)
  #= now sum() chr aln, and put this table in the same txt file.
  mass_pct (df, outfile)

  
  
  

def isLargerThan10k (infile):
    line1 = ['zcat', infile]
    line2 = ['wc', '-l']
    p1 = subprocess.Popen(line1, stdout=subprocess.PIPE, env=os.environ)
    p2 = subprocess.Popen(line2, stdin=p1.stdout, stdout=subprocess.PIPE, env=os.environ)
    p1.stdout.close()  #= Allow p1 to receive a SIGPIPE if p2 exits.
    output = p2.communicate()[0].decode("utf-8").rstrip('\n') #= example: b'34416\n'
    print(output, 'communicate (infile)', infile)
    if int(output)/4 > 11000: return True
    return False 

    
def gzfastq_sampling (inGzFastq):
  if not isLargerThan10k (inGzFastq): return inGzFastq
  numSeq=str(10000)
  output = re.sub('.fastq.gz', '.10k.fastq', inGzFastq)
  outGz = re.sub('.10k.fastq', '.10k.fastq.gz', output)
  #= cmd
  cmd = "zcat " + inGzFastq + " | awk '{ printf(\"%s\",$0); n++; if(n%4==0) {"
  cmd += "printf(\"\\n\");} else { printf(\"\\t\");} }' |"
  cmd += "awk -v k=" + numSeq + " 'BEGIN{srand(systime() + PROCINFO[\"pid\"]);}{s=x++<k?x1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' |"
  cmd += "awk -F\"\\t\" '{print $1\"\\n\"$2\"\\n\"$3\"\\n\"$4 > \"sub.fastq\"}'"
  print(cmd);subprocess_tryexcept (cmd)
  cmd = "sleep 10"
  print(cmd);subprocess_tryexcept (cmd)
  cmd = "mv sub.fastq " + output
  print(cmd);subprocess_tryexcept (cmd)
  #= compress
  cmd = 'gzip -c ' + output + ' > ' + outGz
  print(cmd);subprocess_tryexcept (cmd)
  return outGz
  
def prepare_IGVready_files (ref_fasta, alignmentSortedBam, rep_igv):
  fai = make_fai (ref_fasta)
  bai = make_bai (alignmentSortedBam)
  cmd = 'cp ' + ref_fasta + ' ' + rep_igv + ';'
  cmd += 'cp ' + fai + ' ' + rep_igv + ';'
  cmd += 'cp ' + alignmentSortedBam + ' ' + rep_igv + ';'
  cmd += 'cp ' + bai + ' ' + rep_igv + ';'
  print(cmd)
  subprocess_tryexcept (cmd) 
  
def run_igv_demo (): 
  #= this will generate demo 2 png files in src/snapshot_dir/
  cmd = 'module load IGVgui/2.13.1;'
  cmd += 'igv.sh -b test_edge2.bat'
  print(cmd);subprocess_tryexcept (cmd) 







    
    
def pre_align_processing ():
  if not os.path.exists(rep_output): os.makedirs(rep_output)
  outfile = fastqmerged = rep_output + libname + '.fastq.gz'
  merge_fastq_compress (rep_input, outfile)
  #= Remove reads by read quality
  infile = fastqmerged
  if trim_poorSequencing_Nanopore:
    outfile = nanofiltered = NanoFilt (infile)
    infile = nanofiltered
  #= cutadapt
  if trimAdapters:
    rep_cutadapt = rep_output + 'cutadapt/'
    if not os.path.exists(rep_cutadapt): os.makedirs(rep_cutadapt)
    #= Remove 5prime adapter (must locate at extremity)
    inBasename = os.path.basename(infile)
    outfile = delADAPTER5p = rep_cutadapt + re.sub('.fastq.gz', '.delAdapter5p.fastq.gz', inBasename)
    stdout_cutadapt_file = rep_cutadapt + re.sub('.fastq.gz', '.cutadapt5p.stdout', inBasename)
    cutadapt_5prime (infile, adapter5prime, outfile, stdout_cutadapt_file)
    #= Remove 3prime adapter (must locate at extremity)
    infile = delADAPTER5p
    inBasename = os.path.basename(infile)
    outfile = delADAPTER5p3p = rep_cutadapt + re.sub('.delAdapter5p.fastq.gz', '.delAdapter5p3p.fastq.gz', inBasename) 
    stdout_cutadapt_file = rep_cutadapt + re.sub('.fastq.gz', '.cutadapt3p.stdout', inBasename)
    cutadapt_3prime (infile, adapter3prime, outfile, stdout_cutadapt_file)
    cmd = 'rm ' + delADAPTER5p;print(cmd);subprocess_tryexcept (cmd)
    infile = delADAPTER5p3p
  #= view read quality
  rep_qc = rep_output + 'preAlignQC/'
  if not os.path.exists(rep_qc): os.makedirs(rep_qc)
  fastqc (infile, rep_qc)
  NanoPlot (infile, rep_qc) #= not fully supported, no vis report
  return infile, fastqmerged
  
def map_haplotype (infile):
  #= mapping to genes: Cargo, AmpR, KanR, Rep2, Cap9, HPV18, Ad5, Lambda, plasmidBackbone
  cat_ref ()
  rep_mapping = rep_output + 'mapping_' + ref_fasta_key + '/'
  if not os.path.exists(rep_mapping): os.makedirs(rep_mapping)
  inBasename = os.path.basename(infile)
  if aligner == 'minimap2':
    outfile = alignmentSam = rep_mapping + re.sub('.fastq.gz', '.minimap.sam', inBasename)
    minimap2 (libname, ref_fasta, infile, outfile)
  elif aligner == 'bwamem':
    outfile = alignmentSam = rep_mapping + re.sub('.fastq.gz', '.bwamem.sam', inBasename)
    bwa_mem (infile, ref_fasta, outfile)
  alignmentSortedBam = sam2SortedBam (alignmentSam)
  infile = alignmentSortedBam
  #= filter low quality alignments
  if trim_poorAlignment:
    outfile = alignmentSortedBamHQ = filter_mappingQuality (infile)
    infile = alignmentSortedBamHQ
  #= basic stats on the alignment results (stats.txt, stats.html, idxstats.tsv, flagstat.txt)
  basic_stats (infile)
  flagstat (infile)
  f_idx = idxstats (infile) #= f_idx is required by bam2haplotype
  #= annotation of each read: convert bam to read_haplotype.tsv and read_coors.tsv
  #= and make haplotype distribution file, raw only
  rep_haplo = rep_mapping + 'haplotype/'
  if not os.path.exists(rep_haplo): os.makedirs(rep_haplo)
  file_read_haplotype_RAWDATA = bam2haplotype (infile, rep_haplo, fastqmerged, f_idx, '_')
  #time.sleep(60)
  describe (file_read_haplotype_RAWDATA) #= describe() dataset for alignment basecounts
  #= get reads in fasta format for each of all haplotypes
  infile = file_read_haplotype_RAWDATA
  if getHaploFasta: getFastaFiles_for_allHaplotypes (infile, fastqmerged)
  ##= draw heatmap representing haplotype for each read
  cmd = inlineLoad
  cmd += 'python heatmap_haplotype.py ' + infile
  print(cmd);subprocess_tryexcept (cmd)
  #= draw hist2d for residual chromosomes
  coors = re.sub('.read_haplotype.tsv', '.aln_coors.tsv', infile)
  cmd = inlineLoad
  cmd += 'python multiple_hist2d.py ' + coors
  print(cmd);subprocess_tryexcept (cmd)
  return f_idx, file_read_haplotype_RAWDATA

def map_lambda_fragments (infile, f_idx):
  #= mapping to 14 lambda fragments
  LambdaCt = getLambdaCt (f_idx)
  if int(LambdaCt) > min_LambdaCt: 
    print('================ Begin align lambda fragments =============================')
    ref_fasta_key2='lambda_fragments'
    ref_fasta='../dbs/lambda_fragments.fa'
    #= mapping
    rep_mapping = rep_output + 'mapping_' + ref_fasta_key2 + '/'
    if not os.path.exists(rep_mapping): os.makedirs(rep_mapping)
    inBasename = os.path.basename(infile)
    if aligner == 'minimap2':
      outfile = alignmentSam = rep_mapping + re.sub('.fastq.gz', '.minimap.sam', inBasename)
      minimap2 (libname, ref_fasta, infile, outfile)
    elif aligner == 'bwamem':
      outfile = alignmentSam = rep_mapping + re.sub('.fastq.gz', '.bwamem.sam', inBasename)
      bwa_mem (infile, ref_fasta, outfile)
    alignmentSortedBam = sam2SortedBam (alignmentSam)
    infile = alignmentSortedBam
    #= filter low quality alignments
    if trim_poorAlignment:
      outfile = alignmentSortedBamHQ = filter_mappingQuality (infile)
      infile = alignmentSortedBamHQ
    #= basic stats on the alignment results (stats.txt, stats.html, idxstats.tsv, flagstat.txt)
    basic_stats (infile)
    flagstat (infile)
    f_idx = idxstats (infile) #= f_idx is required by bam2haplotype
    #= annotation of each read: convert bam to read_haplotype.tsv and read_coors.tsv
    rep_haplo = rep_mapping + 'haplotype/'
    if not os.path.exists(rep_haplo): os.makedirs(rep_haplo)
    file_lambda_haplotype = bam2haplotype (infile, rep_haplo, fastqmerged, f_idx, file_read_haplotype_RAWDATA)
    #time.sleep(60)
    print('================ Make haplotype distribution file, raw and adjusted ====================')
    #= make haplotype distribution file by previous line: bam2haplotype using haplotype_breakdown.py
    infile = file_read_haplotype_RAWDATA
    file_read_haplotype_ADJDATA = re.sub('.read_haplotype.tsv', '.read_haplotype.adj.tsv', infile)
    return file_lambda_haplotype
   
def write_html_report ():
  cmd = inlineLoad
  cmd += 'python write_html.py ' + rep_output + ' ' + libname + '.fastq.gz ' + ref_fasta_key 
  print(cmd);subprocess_tryexcept (cmd)




def write_bat (ref_fasta, infile, rep_igv, chr_regions):
  '''
  new
  snapshotDirectory snapshot_dir
  load my.bam
  genome myplasmid.fa
  maxPanelHeight 10000
  colorBy NONE
  goto plasmid1:1-14047
  sort position
  squish
  snapshot myplasmid1.igvsnapshot.png
  exit 
  '''
  bam = rep_igv + os.path.basename(infile)
  fas = rep_igv + os.path.basename(ref_fasta)
 
  rep_igv_snapdir = rep_igv + 'snapshot_dir/'
  if not os.path.exists(rep_igv_snapdir): os.makedirs(rep_igv_snapdir)
      
  batfile = rep_igv + 'my.bat'
  with open (batfile, 'w') as fh:
    print('new', file=fh)
    print('snapshotDirectory ' + rep_igv_snapdir, file=fh)
    print('load ' + bam, file=fh)
    print('genome ' + fas, file=fh)
    print('maxPanelHeight 10000', file=fh)
    for chr_region in chr_regions: #= bsg168:1-14987
      png = 'igvsnapshot.' + re.sub(':', 'x', chr_region) + '.png'
      print('goto ' + chr_region, file=fh)
      print('sort position', file=fh)
      print('squish', file=fh)
      print('snapshot ' + png, file=fh)
    print('exit', file=fh)
  return

def getFastaStat (infile):
  D = {}
  with open (infile, 'r') as fh:
    for L in fh:
      L = L.rstrip('\n')
      if L.startswith('>'): 
        chrname = L[1:].split(' ')[0]
        D[chrname] = ''
      else: D[chrname] += L
  return D

def define_regions (ref_fasta):
  #chr_regions = ['bsg168:1-14987', ...]
  chr_regions = []
  D = getFastaStat (ref_fasta)
  for k, v in D.items():
    #print(k, len(v))
    region = k + ':1-' + str(len(v)-1)
    chr_regions.append(region)
  return chr_regions
  
def run_igv_realdeal (bat): 
  cmd = 'module load anaconda3/4.9.2.compchem;module load IGVgui/2.13.1;'
  cmd += 'igv.sh -b ' + bat # Can't connect to X11 window server 
  #cmd += 'xvfb-run --auto-servernum igv.sh -b ' + bat 
  print(cmd);subprocess_tryexcept (cmd) 
  
  
def map_plasmid (infile, option=0):
  rep_mapping = rep_output + 'mapping_plasmid' + '/'
  ref_fasta = '../dbs/' + plasmids_fasta
  if not os.path.exists(rep_mapping): os.makedirs(rep_mapping)
  inBasename = os.path.basename(infile)
  if aligner == 'minimap2':
    outfile = alignmentSam = rep_mapping + re.sub('.fastq.gz', '.minimap.sam', inBasename)
    minimap2 (libname, ref_fasta, infile, outfile)
  elif aligner == 'bwamem':
    outfile = alignmentSam = rep_mapping + re.sub('.fastq.gz', '.bwamem.sam', inBasename)
    bwa_mem (infile, ref_fasta, outfile)
  alignmentSortedBam = sam2SortedBam (alignmentSam)
  infile = alignmentSortedBam
  #= filter low quality alignments
  if trim_poorAlignment:
    outfile = alignmentSortedBamHQ = filter_mappingQuality (infile)
    infile = alignmentSortedBamHQ
  #= basic stats on the alignment results (stats.txt, stats.html, idxstats.tsv, flagstat.txt)
  basic_stats (infile)
  flagstat (infile)
  f_idx = idxstats (infile)
  #= igv viewer for the alignments
  if option == 1: #= when the input is sampled.
    rep_igv = rep_mapping + 'igv/'
    if not os.path.exists(rep_igv): os.makedirs(rep_igv)
    prepare_IGVready_files (ref_fasta, infile, rep_igv)
    if hpc == 'edge': 
      #run_igv_demo () #= generate demo snapshots in src/snapshots
      chr_regions = define_regions (ref_fasta)
      write_bat (ref_fasta, infile, rep_igv, chr_regions)
      bat = rep_igv + 'my.bat'
      run_igv_realdeal (bat)
  return



###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
def map_human_genome (infile):
  ''' TO DO '''
  rep_mapping = rep_output + 'mapping_humangenome' + '/'
  ref_fasta = '../dbs/08_hg38.clean.fa'
  if not os.path.exists(rep_mapping): os.makedirs(rep_mapping)
  inBasename = os.path.basename(infile)
  if aligner == 'minimap2':
    outfile = alignmentSam = rep_mapping + re.sub('.fastq.gz', '.minimap.sam', inBasename)
    minimap2 (libname, ref_fasta, infile, outfile)
  elif aligner == 'bwamem':
    outfile = alignmentSam = rep_mapping + re.sub('.fastq.gz', '.bwamem.sam', inBasename)
    bwa_mem (infile, ref_fasta, outfile)
  alignmentSortedBam = sam2SortedBam (alignmentSam)
  infile = alignmentSortedBam
  #= filter low quality alignments
  if trim_poorAlignment:
    outfile = alignmentSortedBamHQ = filter_mappingQuality (infile)
    infile = alignmentSortedBamHQ
  #= basic stats on the alignment results (stats.txt, stats.html, idxstats.tsv, flagstat.txt)
  basic_stats (infile)
  flagstat (infile)
  f_idx = idxstats (infile) #= f_idx is required by bam2haplotype
  #= annotation of each read: convert bam to read_haplotype.tsv and read_coors.tsv
  #= and make haplotype distribution file, raw only
  rep_haplo = rep_mapping + 'haplotype/'
  if not os.path.exists(rep_haplo): os.makedirs(rep_haplo)
  file_read_haplotype_RAWDATA = bam2haplotype (infile, rep_haplo, fastqmerged, f_idx, '_')
  #time.sleep(60)
  describe (file_read_haplotype_RAWDATA) #= describe() dataset for alignment basecounts
  #= get reads in fasta format for each of all haplotypes
  infile = file_read_haplotype_RAWDATA
  #= hg analysis by default extract reads (for now)
  getFastaFiles_for_allHaplotypes (infile, fastqmerged)
  return file_read_haplotype_RAWDATA
  #return
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################







  

  


  
if __name__ == '__main__':
  ActivateNanoporePreprocesing = True #args.n  
  ActivateLambda = args.l
  ActivateGenotyping = args.g
  ActivatePlasmid = args.p
  ActivateHG38 = args.c
  
  #= Lambda analysis requires output from genotyping analysis
  if ActivateLambda: ActivateGenotyping = True

  time_a = logArg ()
  if ActivateNanoporePreprocesing:
    preAlignFile, fastqmerged = pre_align_processing ()
  else:
    #==BUG BUG BUG
    #= if no preprocessing, rep_input must be a single file of fastq.gz
    #= some other bug still persists. Identify!   
    preAlignFile, fastqmerged = rep_input, rep_input

  rep_results = rep_output + 'results/'
  if not os.path.exists(rep_results): os.makedirs(rep_results)
  
  reporttxt = rep_results + 'report.txt'
  if ActivateGenotyping:
    f_idx, file_read_haplotype_RAWDATA = map_haplotype (preAlignFile)
    file1 = re.sub('.read_haplotype.tsv', '.chr_positivity_intact_rates.tsv', file_read_haplotype_RAWDATA)
    file2 = re.sub('.read_haplotype.tsv', '.describe_reads.tsv', file_read_haplotype_RAWDATA)
    file2_2 = re.sub('.read_haplotype.tsv', '.describe_reads.mass_pct.tsv', file_read_haplotype_RAWDATA)
    file3 = re.sub('.read_haplotype.tsv', '.read_haplotype_stats2.tsv', file_read_haplotype_RAWDATA)
    
    cmd = 'echo "################## table: Gene positivity and intact rates" > ' + reporttxt + ';'
    cmd += 'cat ' + file1  + ' >> ' + reporttxt + ';'
    cmd += 'echo >> ' + reporttxt + ';'
    cmd += 'echo "################## table: Gene described" >> ' + reporttxt + ';'
    cmd += 'cat ' + file2  + ' >> ' + reporttxt + ';'
    cmd += 'echo >> ' + reporttxt + ';'
    cmd += 'echo "################## table: Gene mass percentage" >> ' + reporttxt + ';'
    cmd += 'cat ' + file2_2  + ' >> ' + reporttxt + ';'
    cmd += 'echo >> ' + reporttxt + ';'
    if not ActivateLambda:
      cmd += 'echo "################## table: Gene haplotype distribution" >> ' + reporttxt + ';'
      cmd += 'cat ' + file3  + ' >> ' + reporttxt + ';'
      cmd += 'echo >> ' + reporttxt
    print(cmd);subprocess_tryexcept (cmd)
    
    rep_mapping = rep_output + 'mapping_' + ref_fasta_key + '/'
    cmd = 'cp ' + rep_mapping + 'haplotype/*.png ' + rep_results + ';'
    cmd += 'rm ' + rep_results + '*.target_breakdown.png'
    print(cmd);subprocess_tryexcept (cmd)
 
  if ActivateLambda: 
    fileLambda = map_lambda_fragments (preAlignFile, f_idx)
    file4 = re.sub('.read_haplotype.tsv', '.haplotype_distribution.adj.tsv', file_read_haplotype_RAWDATA)
    file5 = re.sub('.read_haplotype.tsv', '.chr_positivity_intact_rates.tsv', fileLambda)
    cmd = 'echo "################## table: Gene haplotype distribution (adj)" >> ' + reporttxt + ';'
    cmd += 'cat ' + file4  + ' >> ' + reporttxt + ';'
    cmd += 'echo >> ' + reporttxt + ';'
    cmd += 'echo "################## table: lambda fragment positivity and intact rates" >> ' + reporttxt + ';'
    cmd += 'cat ' + file5  + ' >> ' + reporttxt + ';'
    cmd += 'echo >> ' + reporttxt
    print(cmd);subprocess_tryexcept (cmd)
       
    #rep_mapping = rep_output + 'mapping_lambda_fragments/'
    #cmd = 'cp ' + rep_mapping + 'haplotype/*.png ' + rep_results
    #print(cmd);subprocess_tryexcept (cmd)
  
  
  if ActivateNanoporePreprocesing:
    file6 = rep_output + 'preAlignQC/' +libname + '.delAdapter5p3p.fastq.gz.NanoStats.txt'
    cmd = 'echo "################## table: Nanopore sequencing quality statistics" >> ' + reporttxt + ';'
    cmd += 'cat ' + file6  + ' >> ' + reporttxt + ';'
    cmd += 'echo >> ' + reporttxt
    print(cmd);subprocess_tryexcept (cmd)
    
  #write_html_report ()
  
 
  if ActivateHG38:
    fileHG = map_human_genome (preAlignFile)
    file7 = re.sub('.read_haplotype.tsv', '.describe_reads.tsv', fileHG)
    file7_2 = re.sub('.read_haplotype.tsv', '.describe_reads.mass_pct.tsv', fileHG)
    file8 = re.sub('.read_haplotype.tsv', '.read_haplotype_stats2.tsv', fileHG)
    cmd = 'echo "################## table: cell host chromosome described" >> ' + reporttxt + ';'
    cmd += 'cat ' + file7  + ' >> ' + reporttxt + ';'
    cmd += 'echo >> ' + reporttxt + ';'
    cmd += 'echo "################## table: cell host chromosome mass percentage" >> ' + reporttxt + ';'
    cmd += 'cat ' + file7_2  + ' >> ' + reporttxt + ';'
    cmd += 'echo >> ' + reporttxt + ';'
    cmd += 'echo "################## table: cell host haplotype distribution" >> ' + reporttxt + ';'
    cmd += 'cat ' + file8  + ' >> ' + reporttxt + ';'
    cmd += 'echo >> ' + reporttxt
    print(cmd);subprocess_tryexcept (cmd)
    
  if ActivatePlasmid:
    map_plasmid (preAlignFile)
    #= sampled and igv infiles generated.
    sampled = gzfastq_sampling (preAlignFile)
    map_plasmid (sampled, 1) 
    
    rep_snapshot = rep_output + 'mapping_plasmid/igv/snapshot_dir/'
    cmd = 'cp ' + rep_snapshot + '*.png ' + rep_results
    print(cmd);subprocess_tryexcept (cmd)

    
    
    


  




    
    
  #= log ending time
  time_b = datetime.datetime.now()
  NowStr = time_b.strftime("%a %b %d, %Y. %H:%M:%S")
  print('total running time: ', time_b - time_a)
  logcmd (' ')
  logcmd (NowStr)
  logcmd ('total running time: ')
  logcmd (time_b - time_a)




  
##= BUG: getHaploFasta not working
#time python longreadAAVpipelinecjw.py --rep_input ../input_lib_100reads/ --trim_poorSequencing_Nanopore --trimAdapters yes --trim_poorAlignment --getHaploFasta yes
#= stdout related to this bug
''' 
module load miniconda3/4.3.21;python extractReadsByHaplotypes.py ../output/myLibrary/mapping_universalAAVplatform/haplotype/myLibrary.seqQ10.delAdapter5p3p.bwamem.mapQ1.read_haplotype.tsv ../output/myLibrary/myLibrary.fastq.gz
Traceback (most recent call last):
  File "extractReadsByHaplotypes.py", line 64, in <module>
    RUN (infile, fqgzFile)
  File "extractReadsByHaplotypes.py", line 57, in RUN
    extarct_fastq_byID (outfile_listID, fqgzFile)
  File "extractReadsByHaplotypes.py", line 31, in extarct_fastq_byID
    subprocess.check_call(cmd, shell=True)
  File "/camhpc/pkg/miniconda3/4.3.21/centos6/lib/python3.6/subprocess.py", line 311, in check_call
    raise CalledProcessError(retcode, cmd)
subprocess.CalledProcessError: Command 'zcat ../output/myLibrary/myLibrary.fastq.gz | grep -F -A 3 -f ../output/myLibrary/mapping_universalAAVplatform/haplotype/myLibrary.seqQ10.delAdapter5p3p.bwamem.mapQ1.readIDs.tsv | grep -v '^--$' > ../output/myLibrary/mapping_universalAAVplatform/haplotype/myLibrary.seqQ10.delAdapter5p3p.bwamem.mapQ1..fastq' returned non-zero exit status 1.
An exception occurred
Command 'module load miniconda3/4.3.21;python extractReadsByHaplotypes.py ../output/myLibrary/mapping_universalAAVplatform/haplotype/myLibrary.seqQ10.delAdapter5p3p.bwamem.mapQ1.read_haplotype.tsv ../output/myLibrary/myLibrary.fastq.gz' returned non-zero exit status 1
'''


