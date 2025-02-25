'''
2021-03-09
author: Chao-Jung Wu <wu.chaojung@gmail.com>

to avoid macOS warning on tight layout: 
https://stackoverflow.com/questions/37309559/using-matplotlib-giving-me-the-following-warning-userwarning-tight-layout
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
fig.set_tight_layout(True)
fig.savefig('asd.pdf')  # No warning now
'''
# ../input/tai18_78/mapping_universalAAVplatform/haplotype/
# [x] tai18_78.bwamem.read_coors.tsv
# [o] tai18_78.bwamem.read_haplotype.tsv
# [?] tai18_78.bwamem.read_haplotype_stats2.tsv
#################################
# plan
#################################
# take read_haplotype.tsv => read_haplotype.adj.tsv
# [v] remove read length < 200 nt
# [v] combine hg38 chrs into one annotation
# [v] change alnLen = 0 if alnLen < 200
# [v] register genotype g12345678 in this file
# [v] register adj_readcount in this file
#
# make another file => haplotype_stats.adj.tsv
# based on read_haplotype.adj.tsv
# record unAdj haplotype distribution
# record Adj haplotype distribution
from __future__ import print_function
import re, os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#= cargo name will be added by: Cols_toCheck.insert(0, cargoName)
#Cols_toCheck = ['AmpR', 'KanR', 'Rep2', 'Cap9', 'HPV18', 'Ad5', 'hg38', 'lambda']
#Cols_toCheck = ['AmpR', 'KanR', 'Rep2', 'Cap9', 'HPV18', 'Ad5', 'hg38', 'backbone'] 
Cols_toCheck = ['AmpR', 'KanR', 'Rep2', 'Cap9', 'HPV18', 'Ad5', 'backbone']

def assign_genotype (aList):
  ''' aList = [12, 0, 0, 39, 0, 5] ==> 100101 '''
  genotypeBinary = 'g'
  for i in aList:
    if int(i) == 0: genotypeBinary += '0'
    else: genotypeBinary += '1'
  return genotypeBinary

def assign_genotype_description (aList):
  ''' aList = [12, 0, 0, 39, 0, 5] ==> 100101 '''
  genotypeDesc = ''
  count = 0
  for i in aList:
    if int(i) > 0: genotypeDesc += Cols_toCheck[count] + ' + '
    count += 1
  if genotypeDesc == '': genotypeDesc = 'no_annotation'
  return genotypeDesc.rstrip(' + ')
  
def process_reads_to_readsAdj (infile, spl, min_lenRead, min_lenAlnChr, option='spl'):
  #= infile = ../input/tai20_03/mapping_universalAAVplatform/haplotype/tai20_03.bwamem.read_haplotype.tsv
  #= outfile = xxx/yyy.delAdapter5p3p.bwamem.read_haplotype.adj.tsv
  df = pd.read_csv(infile, sep='\t')
  df.drop('haplotype', axis=1, inplace=True)
  Cols = df.columns.values.tolist()
  cargoName = [x for x in Cols if x.startswith('ITR2ITR')][0]
  if cargoName != 'ITR2ITR': print('Warning: cargo name is not ITR2ITR, this may cause error during making plot.')
  Cols_toCheck.insert(0, cargoName)
  #= remove reads if read_length < 200
  #df.drop(df[ (df['read_length'] < min_lenRead) ].index , inplace=True)
  #= combine chrs of hg38 into one chr
  #df['hg38'] = df[ [Chr for Chr in Cols if Chr.startswith('chr')] ].sum(axis=1)
  #= change alnLen = 0 if alnLen < 200
  #for col in Cols_toCheck: df.loc[df[col] < min_lenAlnChr, col] = 0
  #= register genotype g12345678
  df['genotype'] = df[Cols_toCheck].apply(assign_genotype, axis = 1)
  df['genotypeDescription'] = df[Cols_toCheck].apply(assign_genotype_description, axis = 1)
  #= adjust readcount based on read length
  if option == 'noSpl': df['adj_readcount_basedon_readlength'] = df['read_length'].apply(lambda x: 1)
  else: df['adj_readcount_basedon_readlength'] = df['read_length'].apply(lambda x: spl(x))
  #= SAVE df to file
  outfile = re.sub('.read_haplotype.tsv', '.read_haplotype.adj.tsv', infile)
  df.to_csv(outfile, sep='\t', index=False)
  return outfile

def process_readsAdj_to_haplotype_distribution (infile):
  #=infile=../input/tai18_78/mapping_universalAAVplatform/haplotype/tai18_78.bwamem.read_haplotype.adj.tsv
  outfile = re.sub('.read_haplotype.adj.tsv', '.haplotype_distribution.adj.tsv', infile)
  df = pd.read_csv(infile, sep='\t')
  noLambda = df.copy()
  noLambda.drop(df[ (df['lambda'] != 0) ].index, inplace=True)
  total_noLambda = len(noLambda)
  totalAdj_noLambda = noLambda['adj_readcount_basedon_readlength'].sum()
  #= groupby
  grouped_rawcount = df.groupby(['genotype', 'genotypeDescription'], as_index=False)['readID'].count()
  grouped_rawcount.rename(columns={'readID':'raw_readcount'}, inplace=True)
  grouped_adjcount = df.groupby('genotype', as_index=False)['adj_readcount_basedon_readlength'].sum()
  df = pd.merge(grouped_rawcount, grouped_adjcount, on='genotype',how='left')
  #= calculate %
  df['%LIB_noLambda'] = df['raw_readcount'].apply(lambda x: round(x*100/float(total_noLambda), 2) )
  df['%adjLIB_noLambda'] = df['adj_readcount_basedon_readlength'].apply(lambda x: round(x*100/float(totalAdj_noLambda), 2) )

  cols = 'genotype,genotypeDescription,raw_readcount,%LIB_noLambda,adj_readcount_basedon_readlength,%adjLIB_noLambda'.split(',')
  df = df[cols]

  #= remove the line of no_annotation
  indexNames = df[ df['genotypeDescription'] == 'no_annotation'].index
  df.drop(indexNames, inplace=True)

  #= SAVE df to file
  df.to_csv(outfile, sep='\t', index=False)
  #= print footnotes
  fho = open (outfile, 'a')
  print('# NB total reads after removing lambda:', total_noLambda, file=fho)
  print('# NB adjusted total reads after removing lambda: %.2f' % totalAdj_noLambda, file=fho)
  fho.close()
  print(outfile)
  return outfile

#########################################################
#########################################################
#########################################################
def draw_hapDist_barplot (infile):
  #import pandas as pd
  #import numpy as np
  #import matplotlib.pyplot as plt
  outpng = re.sub('.haplotype_distribution.adj.tsv', '.haplotype_distribution.adj.png', infile)
  df = pd.read_csv(infile, comment='#', sep='\t')
  #= remove genotypes no_annotation, lambda, ITR2ITR
  #= it is necessary to remove them, for (1) they are large, for (2) sometimes lambda is negative -> need to verify again and why
  indexNames1 = df[ df['genotypeDescription'] == 'no_annotation'].index
  indexNames2 = df[ df['genotypeDescription'] == 'lambda'].index
  indexNames3 = df[ df['genotypeDescription'] == 'ITR2ITR'].index
  df.drop(indexNames1, inplace=True)
  df.drop(indexNames2, inplace=True)
  df.drop(indexNames3, inplace=True)
  #= plot, set up to avoid warning on tight layout for MacOS
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  fig.set_tight_layout(True)
  # set width of bars
  barWidth = 0.25
  # set heights of bars
  bars1 = df['%LIB_noLambda'].tolist()
  bars2 = df['%adjLIB_noLambda'].tolist()
  # Set position of bar on X axis
  r1 = np.arange(len(bars1))
  r2 = [x + barWidth for x in r1]
  # Make the plot #=colors #7f6d5f #557f2d #rgbkymc
  plt.bar(r1, bars1, color='c', width=barWidth, edgecolor='white', label='%LIB')
  plt.bar(r2, bars2, color='g', width=barWidth, edgecolor='white', label='%LIB_adj')
  # Add xticks on the middle of the group bars
  plt.xlabel('Classification of read', fontweight='bold')
  plt.ylabel('% library without lambda reads', fontweight='bold')
  plt.xticks([r + barWidth for r in range(len(bars1))], df['genotypeDescription'].tolist(), rotation=90) ######
  # Create legend & save graphic
  plt.legend()
  plt.title('Distribution of residual haplotypes', fontsize = 15)
  fig.savefig(outpng) #plt.savefig(outpng)
  print('Haplotype distribution is drawn:', outpng)

#########################################################
#########################################################
#########################################################
#########################################################

def test ():
  infile1 = 'infile2.read_haplotype.tsv'
  infile2 = 'test_100.bwamem.read_haplotype.tsv'

  import lambda_spline as ls
  spl = ls.run_spline (infile1)
  outfile1 = process_reads_to_readsAdj (infile2, spl, min_lenRead, min_lenAlnChr)
  process_readsAdj_to_haplotype_distribution (outfile1)

def adjReadcounts (infile, spl):
  min_lenRead, min_lenAlnChr = 200, 200
  outfile1 = process_reads_to_readsAdj (infile, spl, min_lenRead, min_lenAlnChr)
  outfile2 = process_readsAdj_to_haplotype_distribution (outfile1)
  draw_hapDist_barplot (outfile2)

def noAdjReadcounts(infile):
  min_lenRead, min_lenAlnChr, option = 200, 200, 'noSpl'
  outfile1 = process_reads_to_readsAdj (infile, '-', min_lenRead, min_lenAlnChr, option)
  outfile2 = process_readsAdj_to_haplotype_distribution (outfile1)
  draw_hapDist_barplot (outfile2)


##############################
infile = 'test_100.bwamem.haplotype_distribution.adj.tsv'
#draw_hapDist_barplot (infile)
