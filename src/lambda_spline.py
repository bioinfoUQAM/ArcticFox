'''
author: Chao-Jung Wu <wu.chaojung@gmail.com>

lambda_spline

https://docs.scipy.org/doc/scipy/reference/interpolate.html
np.linspace (start, end, count) => Return evenly spaced numbers over a specified interval.

210417
spline smoothing factor causes the adj-readcount sometimes very negative. I tested several lambda libraries, and so far s=0 is the best
'''
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline

#module load miniconda3/4.3.21
fragment_bp = [8454,7242,6369,5686,4822,4324,3675,2323,1929,1371,1264,702,224,117]
fragment_name = 'f01,f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14'.split(',')
hap_name = 'h10000000000000,h01000000000000,h00100000000000,h00010000000000,h00001000000000,h00000100000000,h00000010000000,h00000001000000,h00000000100000,h00000000010000,h00000000001000,h00000000000100,h00000000000010,h00000000000001'.split(',')

s = 0 #= spline smoothing factor


def fit_UnivariateSpline (x, y, outpng='UnivariateSpline.png'):
  ''' 
  1-D smoothing spline fit to a given set of data points. 
  from scipy.interpolate import UnivariateSpline
  '''
  #= plot raw lambda points
  plt.plot(x, y, 'ro-', ms=5)
  #= plot regression curve
  spl = UnivariateSpline(x, y); spl.set_smoothing_factor(s) #= goodness: 0 > 0.1 > None > 1
  xs = np.linspace(min(x), max(x), 1000)
  plt.plot(xs, spl(xs), 'b--', lw=1)
  #= plot settings
  plt.legend(['raw_ratio', 'PolynomialSplineRegression'])
  plt.ylabel('Ratio', fontsize = 15)
  plt.xlabel('Length (nt)', fontsize = 15)
  plt.title('Lambda fragment template sizes vs read ratios (length, ratio)', fontsize = 15)
  plt.savefig(outpng)
  return spl

def plot_ratio_points (x, y, outpng):
  plt.plot(x, y, 'gd-', ms=5) #= adj ratio
  plt.legend(['raw_ratio', 'PolynomialSplineRegression', 'adj_ratio'])
  plt.savefig(outpng)
  
def take_spikeReads_file (infile):
  outfile = re.sub('.read_haplotype.tsv', '.read_lambda_breakdown.summaryTable_adj.tsv', infile)
  outpng = re.sub('.read_haplotype.tsv', '.spikeins_polynomial_spline.png', infile)
  df = pd.read_csv(infile, sep='\t')
  indexNames = df[ df['haplotype'] == 'h00000000000000'].index
  df.drop(indexNames, inplace=True)
  #= cal raw ratio
  grc = grouped_rawcount = df.groupby('haplotype', as_index=False)['readID'].count()
  grc.rename(columns={'readID': 'raw_readcount'}, inplace=True)
  #== to avoid bug, certain haplotypes not in the df => ValueError: 'h00000000001100' is not in list
  grc['bp'] = grc['haplotype'].apply(lambda x: fragment_bp[hap_name.index(x)] if x in hap_name else 0)
  indexNames_bp0 = grc[ grc['bp'] == 0].index
  grc.drop(indexNames_bp0, inplace=True)
  MAX = max (grc['raw_readcount'])
  grc['raw_ratio'] = grc['raw_readcount'].apply(lambda x: round(MAX/(x+0.000001), 5))
  #= plot all lambda points
  #x = grc['bp'].to_numpy(); y = grc['raw_ratio'].to_numpy(); plt.plot(x, y, 'ro-', ms=5)
  #== force outliers not to be considered for spline fitting
  rmOutliers = grc.copy()
  indexNames_outliers = grc[ grc['raw_ratio'] > 100].index
  rmOutliers.drop(indexNames_outliers , inplace=True)
  #== force f13 and f14 not to be considered for spline fitting (I think removing outliers is sufficient.)
  #f13, f14 = 'h00000000000001', 'h00000000000010'
  #indexNames_f13 = grc[ grc['haplotype'] == f13].index; rmOutliers.drop(indexNames_f13, inplace=True)
  #indexNames_f14 = grc[ grc['haplotype'] == f14].index; rmOutliers.drop(indexNames_f14, inplace=True)
  #= spline
  x = rmOutliers['bp'].to_numpy(); y = rmOutliers['raw_ratio'].to_numpy()
  spl = fit_UnivariateSpline (x, y, outpng)
  #= cal adj readcount and ratio, then merge adj table with raw table
  df['adj_readcount'] = df['read_length'].apply(lambda x: spl(x))
  grouped_adjcount = df.groupby('haplotype', as_index=False)['adj_readcount'].sum()
  df = pd.merge(grouped_rawcount, grouped_adjcount, on='haplotype',how='left')
  MAX = max (df['adj_readcount'])
  df['adj_ratio'] = df['adj_readcount'].apply(lambda x: round(MAX/(x+0.000001), 5))
  #= organize table, then save to file
  df = df.sort_values(by=['haplotype'], ascending=True) #= from shortest to longest
  df['fragment'] = df['haplotype'].replace(hap_name,fragment_name)
  df['bp'] = df['fragment'].apply(lambda x: fragment_bp[fragment_name.index(x)])
  df = df['fragment,bp,raw_readcount,adj_readcount,raw_ratio,adj_ratio'.split(',')]
  df.to_csv(outfile, sep='\t', index=False)
  #= plot 3rd line (the lambda adj_ratio line) in the same png file
  y = df['adj_ratio'].to_numpy()
  x = df['bp'].to_numpy()
  plot_ratio_points (x, y, outpng)
  plt.close() #= to close plt properly that was open in this script
  return spl

def run_spline (infile):
  #= infile= xxx.read_haplotype.tsv'
  spl = take_spikeReads_file (infile)
  return spl

########################################
#infile = 'infile2.read_haplotype.tsv'; run_spline (infile) #s=0
#infile = 'infile2.tai20_10.read_haplotype.tsv'; run_spline (infile) #s=0
#infile = 'infile2.tai20_07.read_haplotype.tsv'; run_spline (infile) #= s=0
#infile = 'pacbio_spikeins.bwamem.read_haplotype.tsv'; run_spline (infile) #= s=0
#infile = 'Tai20_04.bwamem.read_haplotype.tsv'; run_spline (infile) #s=0
#infile = 'Tai20_10.bwamem.read_haplotype.tsv'; run_spline (infile) #s=0

#run_spline (infile)
