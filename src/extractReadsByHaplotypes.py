'''
author: Chao-Jung Wu <wu.chaojung@gmail.com>
2020-10-07

extractReadsByHaplotypes.py
'''
import sys
import subprocess
import re
import os
#
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
#


def fastq2fasta (fq):
  fa = re.sub('.fastq', '.fasta', fq)
  cmd = "paste - - - - < " + fq + " | cut -f 1,2 | sed 's/^@/>/' | tr '\t' '\n' > " + fa
  #cmd += '; rm ' + fq
  #print(cmd)
  subprocess.check_call(cmd, shell=True)

def extarct_fastq_byID (listIDfile, fqgzFile):
  '''
  zcat ../output/test_new_2files/test_new_2files.fastq.gz | grep -F -A 3 -f listID.txt | grep -v '^--$' > extractOverpackage.fastq
  '''
  outfile_fq = re.sub('.readIDs.tsv', '.fastq', listIDfile)
  cmd = 'zcat ' + fqgzFile + " | grep -F -A 3 -f " + listIDfile + " | grep -v '^--$' > " + outfile_fq
  subprocess.check_call(cmd, shell=True)
  fq = outfile_fq
  fastq2fasta (fq)

def get_selected_ids (haplotypeF, somehaplotypes, outfile):
  '''
  infile = '../output/test_new_2files/test_new_2files.fastq.gz'
  haplotypeF = 'test_new_2files.delAdapter5p3p.bwamem.sorted.read_haplotype.tsv'
  somehaplotypes = ['h010001']
  '''
  df = pd.read_csv(haplotypeF, sep='\t')
  op = df.loc[df['haplotype'].isin(somehaplotypes)]
  rid = op['readID']
  rid.to_csv(outfile, sep='\t', index=False)

def RUN (infile, fqgzFile):
  '''
  infile = '../output/test_new_2files/test_new_2files.fastq.gz'
  haplotypeF = '../output/test_new_2files/mapping_modifiedRef/test_new_2files.delAdapter5p3p.bwamem.sorted.read_haplotype.tsv'
  somehaplotypes = ['h100010']
  '''
  stat_file = re.sub('.read_haplotype.tsv', '.read_haplotype_stats.tsv', infile)
  with open (stat_file, 'r') as fh: haplotypes = [x.rstrip('\n').split('\t')[0] for x in fh.readlines()]
  for haplotype in haplotypes:
    val = sum([int(x) for x in list(haplotype[1:])])
    if val == 0: continue
    outfile_listID = re.sub('.read_haplotype.tsv', '.' + haplotype + 'readIDs.tsv', infile)
    get_selected_ids (infile, [haplotype], outfile_listID)
    extarct_fastq_byID (outfile_listID, fqgzFile)

if __name__ == '__main__' :
  if not len(sys.argv) == 3:
    msg = 'usage: script.py in.read_haplotype.tsv in.fastq.gz'
    sys.exit(msg)
  infile, fqgzFile = sys.argv[1], sys.argv[2]
  RUN (infile, fqgzFile)


'''
module load miniconda3/4.3.21;python extractReadsByHaplotypes.py ../output/test_new_2files/mapping_modifiedRef/test_new_2files.delAdapter5p3p.bwamem.sorted.read_haplotype.tsv ../output/test_new_2files/test_new_2files.fastq.gz


220603
python extractReadsByHaplotypes.py "/camhpc/home/cwu3/gitRepo/ArcticFox/output/TST11802_bc01_sub10k/mapping_humangenome/haplotype/TST11802_bc01_sub10k.delAdapter5p3p.bwamem.read_haplotype.tsv" "/camhpc/home/cwu3/gitRepo/ArcticFox/output/TST11802_bc01_sub10k/TST11802_bc01_sub10k.fastq.gz"
'''
