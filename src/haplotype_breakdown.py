'''
author: Chao-Jung Wu <wu.chaojung@gmail.com>
2020-10-14

haplotype_breakdown.py

Process_breakdown_and_draw () in camhpc because of X server issue. Hopefully I can reactivate it after transfer to Edge.
'''
from __future__ import print_function
import re, os, sys
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mpl
import matplotlib.pyplot as plt
#
import lambda_spline as ls
import adjReadcounts as arc
import subprocess

rangeClasses = [-0.1, .50, .70, .90, 1.0]
fragment_bp = [8454,7242,6369,5686,4822,4324,3675,2323,1929,1371,1264,702,224,117]

def subprocess_tryexcept (cmd):
  try:
    subprocess.check_call(cmd, shell=True)
  except Exception as e:
    print("An exception occurred")
    print(e)
  #os.system(cmd) ==> os.system(cmd)
    
    
    

def tsv2list (tsv):
  with open (tsv, 'r') as fh: L = [i.rstrip('\n').split('\t') for i in fh.readlines()]
  return L

def tsv2dict (tsv):
  L = tsv2list (tsv)
  D = {}
  for i in L: D[i[0]] = i[1]
  return D

def __repeat (lower, upper, df, key='%target', adj=False):
  aslice = df[(lower <= df['read_length']) & (df['read_length'] < upper)]
  bins = pd.cut(aslice[key], rangeClasses)
  if adj: aslice_grouped = aslice.groupby(bins)['adj_readcount_basedon_readlength'].agg(['sum'])
  else:   aslice_grouped = aslice.groupby(bins)[key].agg(['count'])
  aslice_grouped.columns = [str(lower) + '-' + str(upper) + ' nt']
  return aslice_grouped, aslice

def RUN_chromosome (infile, outfile, chrName, chrLen, option, printStat=True):
  df = pd.read_csv(infile, sep='\t')
  #= remove lambda reads from breakdown statistics calculation
  len_NoLambda = len(df)
  if 'lambda' in df.columns.values.tolist():
    df_NoLambda = df[df['lambda'] == 0]
    len_NoLambda = len(df_NoLambda)
  #= calculate %target or %chr
  if option == '1': 
    key = '%target'
    df[key] = round (df[chrName] / df['read_length'], 2)
  if option == '2':
    key = '%chromosome'
    df[key] = round (df[chrName] / int(chrLen), 2)
  df_sub = df[df[chrName] > 0] #= any read that contains chr fragment is included in this analysis.
  ranges_definition = [(0, 500), (500, 1000), (1000, 1500), (1500, 2000), (2000, 2500), (2500, 3000), 
                       (3000, 3500), (3500, 4000), (4000, 4500), (4500, 5000), (5000, 5500), 
                       (5500, 6000), (6000, 6500), (6500, 7000), (7000, 7500), (7500, 8000), 
                       (8000, 8500), (8500, 9000), (9000, 9500), (9500, 10000), (10000, 10500), 
                       (10500, 10900), (10900, 11700), (11700, 12000),
                       (12000, 12500), (12500, 13000), (13000, 13500), (13500, 14000), (14000, 14500),
                       (14500, df_sub['read_length'].max() + 1)]
  flag = 0
  for (lower, upper) in ranges_definition:
    if flag == 0:
      aslice_grouped0, aslice = __repeat (lower, upper, df_sub, key)
      flag = 1
    else:
      aslice_grouped, aslice = __repeat (lower, upper, df_sub, key)
      merge = pd.merge(aslice_grouped0,aslice_grouped,on=key,how='left')
      aslice_grouped0 = merge
  if printStat == False: #= when processing lambda fragments
    cumsum = merge.cumsum();cumsum.to_csv(outfile, sep='\t')
    #return
  if printStat == True: merge.to_csv(outfile, sep='\t')

  #= print footnote
  with open (outfile, 'r') as fh:
    lines = fh.readlines()
  lastline = lines[-1].split(']\t')[1].rstrip('\n').split('\t')
  read90plus = sum([int(x) for x in lastline])
  seclastline = lines[-2].split(']\t')[1].rstrip('\n').split('\t')
  read70plus = read90plus + sum([int(x) for x in seclastline])
  
  fho = open (outfile, 'a')
  if len_NoLambda > 0:
      print('#The total number of reads mapped to ' + chrName + ': \tNB=\t' + str(len(df_sub)) + '\t%=\t' + str(round( (len(df_sub)*100) / float(len_NoLambda), 2)), file=fho)
  else: print('#The total number of reads mapped to ' + chrName + ': \tNB=\t' + str(len(df_sub)) + '\t%=\t0', file=fho)     
  print('#The total number of reads in the library without lambda: ' + str(len_NoLambda), file=fho)
  print('#Length of ' + chrName + ' = ' + chrLen + ' nt\n', file=fho)
  fho.close()

  #= calculates various rates for report
  positivityRate = round(len(df_sub)*100 / float(len_NoLambda), 2)
  intactRate = 0.00
  if len(df_sub) > 0: intactRate = round(read90plus*100 / float(len(df_sub)), 2)
  compositeRate = round (positivityRate * intactRate / 100.00, 2)
  data = [chrName, chrLen, str(len_NoLambda), str(len(df_sub)), str(read70plus), str(read90plus), str(positivityRate), str(intactRate), str(compositeRate)]
  return data


def draw_breakdown_stacked_barplot (infile, chrName, chrLen, ax, outpng):
  '''
  import pandas as pd
  import numpy as np
  import matplotlib.pyplot as plt
  '''
  df = pd.read_csv(infile, comment='#', sep='\t')
  Cols = df.columns.values.tolist()
  #= manual transpose
  index = Cols[1:] #['%target', '0-500 nt', ..., '8500-9653 nt']
  percent1 = df.loc[0] #['(-0.1, 0.5]']
  percent2 = df.loc[1] #['(0.5, 0.7]']
  percent3 = df.loc[2] #['(0.7, 0.9]']
  percent4 = df.loc[3] #['(0.9, 1.0]']
  df = pd.DataFrame({'0%-50%'  : percent1,
                     '50%-70%' : percent2,
                     '70%-90%' : percent3,
                     '90%-100%': percent4}, index=index)
  #= plot stacked bar
  df.plot.bar(ax=ax, stacked=True, colormap='Paired') #YlGnBu
  ax.set_title('Reads containing ' + chrName + ' (' + chrLen + ' nt)', fontsize=20)
  ax.set_xlabel('Read length')
  ax.set_ylabel('Read count')
  plt.savefig(outpng)

def Process_breakdown_and_draw (infile, word, D, option):
  #if option == '1': key = '%target'
  #elif option == '2': key = '%chromosome'
  #= draw stacked barplots for each chromosome in one figure
  plt.close('all')
  fig, axarr = plt.subplots(nrows=3, ncols=3, figsize=(24, 16))
  fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.1, hspace=0.3)
  if option == '1':
    fig.suptitle('%Target = len_aln_chromosome / len_read', fontsize=25, y=0.98)
  elif option == '2':
    fig.suptitle('%Chromosome = len_aln_chromosome / len_chromosome', fontsize=25, y=0.98)
  fig.set_tight_layout(True) #= this is for mac
  axes = [axarr[0][0], axarr[0][1], axarr[0][2], 
            axarr[1][0], axarr[1][1], axarr[1][2], 
            axarr[2][0], axarr[2][1], axarr[2][2]]
  count = 0
  outpng = re.sub('.read_haplotype.tsv', word + 'png', infile)

  ############################################################
  # save positivity and intact rates in a text file: xxx.chr_positivity_intact_rates.tsv
  ############################################################
  outrates_file = re.sub('.read_haplotype.tsv', '.chr_positivity_intact_rates.tsv', infile)
  fho = open (outrates_file, 'w')
  print('\t'.join('chr, chrlength, num_totalReads, num_chrReads, read70plus, read90plus, positivityRate(%), intactRate(%), compositeRate(%)'.split(', ')), file=fho)
  #= looping chromosomes
  for chrName, chrLen in D.items():
    outfile = re.sub('.read_haplotype.tsv', word + chrName + '.txt', infile)
    data = RUN_chromosome (infile, outfile, chrName, chrLen, option)
    if chrName != 'lambda': print('\t'.join(data), file=fho)
    if chrName.startswith('chr'): continue
    if len(chrName) == 3 and chrName.startswith('f'): continue
    draw_breakdown_stacked_barplot (outfile, chrName, chrLen, axes[count], outpng); count += 1
  fho.close()


def process_ref_genome (infile, D):
    head, tail = os.path.split(infile)
    word1, word2 = '.target_breakdown.', '.chr_breakdown.'
    Process_breakdown_and_draw  (infile, word1, D, '1')
    Process_breakdown_and_draw  (infile, word2, D, '2')
    #####Process_coverage_statistics (infile, word2, D) #= implementing
    #= merge out txt files into one
    outfile = re.sub('.read_haplotype.tsv', '.target_breakdowns.txt', infile)
    cmd = 'cat ' + head + '/*' + word1 + '*.txt > ' + outfile;print(cmd);os.system(cmd)
    cmd = 'rm ' + head + '/*' + word1 + '*.txt';print(cmd);os.system(cmd)
    outfile = re.sub('.read_haplotype.tsv', '.chr_breakdowns.txt', infile)
    cmd = 'cat ' + head + '/*' + word2 + '*.txt > ' + outfile;print(cmd);os.system(cmd)
    cmd = 'rm ' + head + '/*' + word2 + '*.txt';print(cmd);os.system(cmd)
 
def process_spikeins_genome (infile, D, outfile):
    head, tail = os.path.split(infile)
    word3 = '.lambdafragments.read_target_breakdown.'

    #= write positivity/intact rates report
    Process_breakdown_and_draw  (infile, '.lambdafragments.chr_breakdown.', D, '2')
    '''
    outrates_file = re.sub('.read_haplotype.tsv', '.chr_positivity_intact_rates.tsv', infile)
    count, data = 1, []
    with open (outrates_file, 'r') as fh:
      line = fh.readline().rstrip('\n').split('\t')
      if count == 1: toApp = 'Ratio_basedOnPositivity'
      if count == 2: toApp = baseline = float(line[6])
      if count > 2: toApp = myRatio = round(float(line[6]) / baseline, 2)
      line.append(str(toApp))
      count += 1
      data.append(line)
    with open (outrates_file, 'w') as fh:
      for i in data:
        print('\t'.join(i), file=fh)
    '''
    
    #= process spike-ins
    for chrName, chrLen in D.items():
      outfile3 = re.sub('.read_haplotype.tsv', word3 + chrName + '.txt', infile)
      RUN_chromosome (infile, outfile3, chrName, chrLen, '1', False)
      cmd = 'cat ' + head + '/*' + word3 + '*.txt > ' + outfile;print(cmd);os.system(cmd)
    cmd = 'rm ' + head + '/*' + word3 + '*.txt';print(cmd);os.system(cmd)


def Process_coverage_statistics (infile, word2, D):
  #infile = re.sub('.read_haplotype.tsv', word + chrName + '.txt', infile)

  #= draw stacked barplots for each chromosome in one figure
  plt.close('all')
  fig, axarr = plt.subplots(nrows=3, ncols=3, figsize=(24, 16))
  fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.1, hspace=0.3)
  fig.suptitle('Chromosome coverage level', fontsize=25, y=0.98)
  fig.set_tight_layout(True) #= this is for mac
  axes = [axarr[0][0], axarr[0][1], axarr[0][2], 
            axarr[1][0], axarr[1][1], axarr[1][2], 
            axarr[2][0], axarr[2][1], axarr[2][2]]
  count = 1
  chrName = 'ITR2ITR'

  infile = 'infile3_2.ITR2ITR.bwamem.target_breakdowns.txt'
  coverage_statistics (infile, axes[count], chrName)

  '''
  outpng = re.sub('.read_haplotype.tsv', word + 'png', infile)
  #= looping chromosomes
  for chrName, _chrLen in D.items():
      if chrName.startswith('chr'): continue
      infile = re.sub('.read_haplotype.tsv', word + chrName + '.txt', infile)
      #draw_breakdown_stacked_barplot (outfile, chrName, axes[count], outpng); count += 1
      coverage_statistics (infile, axes[count], chrName); count += 1
  #= properly close the plt
  plt.close()
  '''

def coverage_statistics (infile, ax, chrName):
  outpng = re.sub('.txt', '.stat.png', infile)
  with open (infile, 'r') as fh: 
    L = [x.rstrip('\n').rstrip('\t') for x in fh.readlines() if x.startswith ('#')]
  nbLibNoLambda = L[1].split(': ')[1]
  readsContainingThis = L[0].split('NB=')[1].split('%=')[0].strip('\t')

  df = pd.read_csv(infile, comment='#', sep='\t')
  df['#Read'] = df.sum(axis=1)
  df['%subLib'] = round ( df['#Read']*100 / float(readsContainingThis), 2)
  df['%Lib'] = round (df['#Read']*100 / float(nbLibNoLambda), 2)
  df['Coverage'] = 'Very Low (0%-50%),Low (50%-70%),Medium (70%-90%),High (90%-100%)'.split(',')

  dfSub = df['#Read,%subLib,%Lib,Coverage'.split(',')]
  #= manual transpose
  index = ['%Lib', '%subLib']
  list1, list2, list3, list4 = dfSub.loc[0], dfSub.loc[1], dfSub.loc[2], dfSub.loc[3]
  df2 = pd.DataFrame({'High coverage (90%-100%)'  : list4,
                      'Medium (70%-90%)' : list3,
                      'Low (50%-70%)'    : list2,
                      'Very Low (0%-50%)': list1}, index=index)
  #= plot stacked bar
  order = 'High coverage (90%-100%),Medium (70%-90%),Low (50%-70%),Very Low (0%-50%)'.split(',')
  df2[order].plot.bar(ax=ax, stacked=True, colormap='Paired_r', ylim=(0,100))
  #= annotate value in bar
  for p in ax.patches:
    width, height = p.get_width(), p.get_height()
    x, y = p.get_xy() 
    ax.text(x+width/2, 
            y+height/2, 
            '{:.0f} %'.format(height), 
            horizontalalignment='center', 
            verticalalignment='center')
  ax.set_title(chrName, fontsize=20)
  ax.set_ylabel('%')
  #plt.tight_layout() #(pad=0.4, w_pad=0.5, h_pad=1.0)
  plt.savefig(outpng)
  return

#Process_coverage_statistics ('', '', '')

if __name__ == '__main__' :
  if not len(sys.argv) == 4:
    msg = 'usage: script.py 1.in.lambda.read_haplotype.tsv 2.idxstats.tsv 3.in.library.read_haplotype.tsv'
    msg += '\nor, usage: script.py 1.in.library.read_haplotype.tsv 2.idxstats.tsv 3.--'
    sys.exit(msg)
  infile, f_idx, infile2 = sys.argv[1], sys.argv[2], sys.argv[3]

  D = tsv2dict (f_idx)
  D.pop('*', None)
  if 'AmpR' in D.keys() or 'KanR' in D.keys():
    process_ref_genome (infile, D)
    arc.noAdjReadcounts(infile)
  
  elif 'bsg098' in D.keys() or 'bsg116' in D.keys() or 'perg00328' in D.keys():
    process_ref_genome (infile, D)

  elif 'f14' in D.keys() or 'f05' in D.keys():
    outfile = re.sub('.read_haplotype.tsv', '.read_lambda_breakdown.txt', infile)
    process_spikeins_genome (infile, D, outfile) 
    spl = ls.run_spline(infile)
    arc.adjReadcounts(infile2, spl)
