'''
author: Chao-Jung Wu <wu.chaojung@gmail.com>
2020-10-07

decodeSAMpositions_pysam.py

pysam to decode q_start, q_end, etc
require: module load miniconda3/4.3.21
'''
from __future__ import print_function
import sys, re, pysam
import pandas as pd
#import os.path
#from os import listdir

def tsv2list (tsv):
  with open (tsv, 'r') as fh: L = [i.rstrip('\n').split('\t') for i in fh.readlines()]
  return L

def rank_haplotype (infile):
  in_haplotype = re.sub('.bam', '.read_haplotype.tsv', infile)
  df = pd.read_csv(in_haplotype, sep='\t')
  k = df['haplotype'].value_counts()
  outfile = re.sub('.bam', '.read_haplotype_stats.tsv', infile)
  k.to_csv(outfile, sep='\t')

def RUN (infile):
  ''' infile=BAM '''
  #outfile_1 = re.sub('.bam', '.read_coors.tsv', infile) #= original outfilename
  outfile_1 = re.sub('.bam', '.aln_coors.tsv', infile)
  outfile_2 = re.sub('.bam', '.read_haplotype.tsv', infile)
  samfile = pysam.AlignmentFile(infile, "rb")
  chrnames = samfile.references

  L, D, D_qlen = [], {}, {}
  ## L = [read_data] #= read_data = [q, rl, qs, qe, ql, sid, ss, se, sl]
  ## D = {read:[chr1_aln_pos, chr2_aln_pos, ...], ...}
  ## D_qlen = {read: len(read)}
  samfile = pysam.AlignmentFile(infile, "rb")
  for i in samfile:
    #f  = i.flag
    q  = i.query_name
    rl = len(i.query_sequence) ##rl = i.query_length #(This method is not reliable.)
    qs = i.query_alignment_start
    qe = i.query_alignment_end
    ql = i.query_alignment_length
    s  = i.reference_id
    if s == -1: sid = 'upmapped'
    else: sid = chrnames[s]
    ss = i.reference_start
    se = i.reference_end
    sl = i.reference_length
    data = [q, rl, qs, qe, ql, sid, ss, se, sl]
    L.append(data)
    if q not in D.keys():
      x = []
      for i in range(len(chrnames)+1): x.append([])
      D[q] = x
      D_qlen[q] = rl
    else: D_qlen[q] = max(rl, D_qlen[q]) #= bwa seems to have error on length so needs to select the max among reads of the same readID.
    for i in range(qs, qe): D[q][s].append(i) #### may result in bugs !!!!!!!!! Accessing on 201130 with Dongdong
  samfile.close()

  fho = open (outfile_1, 'w')
  title = 'readID, readLength, q_start, q_end, q_aln_len, chrName, s_start, s_end, s_aln_len'.split(', ')
  print('\t'.join(title), file=fho)
  for i in sorted(L):
    line = '\t'.join( str(j) for j in i )
    print(line, file=fho)
  fho.close()

  fho = open (outfile_2, 'w')
  title = ['readID'] + [x for x in chrnames] + ['unmapped', 'read_length', 'aln_baseCount', 'haplotype']
  print('\t'.join(title), file=fho)
  for k, v in sorted(D.items()):
    v = [len(set(x)) for x in v]
    val = '\t'.join([str(x) for x in v])
    rl = D_qlen[k]
    qpos = sum (v[:-1])
    haplotype = 'h' + ''.join( ['1' if x > 0 else '0' for x in v[:-1]] )
    results = [k, val, str(rl), str(qpos), haplotype]
    print('\t'.join(results), file=fho)
  fho.close()

#def get_f_idx (infile):
#  DIR = '/'.join(os.path.dirname(infile).split('/')[:-1])
#  infiles = [f for f in listdir(DIR) if os.path.isfile(os.path.join(DIR, f))]
#  for i in infiles:
#    if i.endswith('.idxstats.tsv'):
#      return DIR + '/' + i

def expanded_hapAnnotation (infile, f_idx):
  f_hstat = re.sub('.bam', '.read_haplotype_stats.tsv', infile)
  outfile = re.sub('.bam', '.read_haplotype_stats2.tsv', infile)
  #f_idx = get_f_idx (infile) 
  L = tsv2list (f_hstat)[1:] ## 210421 weird bug occured, it was tsv2list (f_hstat), now I have to read from [1:]
  chrs = [i[0] for i in tsv2list (f_idx)][:-1]
  fho = open (outfile, 'w')
  title = '#haplotype,numRead,%,description'
  for i in chrs: title += ',' + i
  print('\t'.join(title.split(',')), file=fho)
  mySum = sum([int(i[1]) for i in L])
  for i in L:
    h = i[0][1:]
    if not len(h) == len(chrs): print('error')
    desc = ''
    for j in range(len(h)):
      i.append(h[j])
      if h[j] == '1': desc += ' + ' + chrs[j]
    desc = desc.strip(' + ')
    i.insert(2, desc)
    percent = str(round(float(i[1]) / mySum, 4)*100)
    i.insert(2, percent)
    print('\t'.join(i), file=fho)
  fho.close()
  #cmd = 'mv ' + outfile + ' ' + f_hstat; os.system(cmd) ## can not remove this file yet, it causes bug because of the title line

if __name__ == '__main__' :
  if not len(sys.argv) == 3:
      msg = 'usage: script.py in.bam idxstats.tsv'
      sys.exit(msg)
  infile, f_idx = sys.argv[1], sys.argv[2]
  RUN (infile)
  rank_haplotype (infile)
  
  expanded_hapAnnotation (infile, f_idx)

##
## module load miniconda3/4.3.21; python test.py x
##
