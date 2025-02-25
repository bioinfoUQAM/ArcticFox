#= module load miniconda3/4.3.21
from __future__ import print_function
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import re, sys

'''
author: Chao-Jung Wu <wu.chaojung@gmail.com>
'''

def make_heatmap_several_png (infile):
  outpng = re.sub('.tsv', '.chrHeatmaps.png', infile)
  if '.adj.chrHeatmaps.png' in outpng: outpng = re.sub('.adj.chrHeatmaps.png', '.chrHeatmaps.png', outpng)
  sns.set(font_scale=3)
  fig, axarr = plt.subplots(nrows=3, ncols=3, figsize=(12*6, 8*6))
  fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.1, hspace=0.3)
  fig.suptitle("\n".join([infile, 'Chromosome (s-) start and end alignment positions (hist2d)']), fontsize=25, y=0.98)
  axes = [axarr[0][0], axarr[0][1], axarr[0][2], 
          axarr[1][0], axarr[1][1], axarr[1][2], 
          axarr[2][0], axarr[2][1], axarr[2][2]]
  df = pd.read_csv(infile, sep='\t', low_memory=False)
  Cols = df.columns.values.tolist()#; print(Cols)
  newCols, candidates = [], ['AmpR', 'KanR', 'Rep2', 'Cap9', 'HPV18', 'Ad5']
  for i in candidates:
    if i in Cols: newCols.append(i)
  if 'chr' in Cols:
    hg38 = [x for x in Cols if x.startswith('chr')]
    df['hg38'] = df[hg38].sum(axis=1)
    newCols.append('hg38')
  if 'ITR2ITR' in Cols:
    cargoName = [x for x in Cols if x.startswith('ITR2ITR')][0]
    newCols.insert(0, cargoName)
  if 'bsg098' in Cols: newCols.insert(0, 'bsg098') ###### specialized line 210504
  if 'bsg116' in Cols: newCols.insert(0, 'bsg116') ###### specialized line 210504
  if 'perg00328' in Cols: newCols.insert(0, 'perg00328') ###### specialized line 210504
  Cols = newCols #= [cargoName, 'AmpR', 'KanR', 'Rep2', 'Cap9', 'HPV18', 'Ad5', 'hg38'] ## 210407 added new line
  data = df[Cols]
  count = 0
  for i in Cols:
    rule = df[i] > 0
    dfExtracted = df[rule]
    result = dfExtracted[Cols]
    if len(result) < 1 or len(result) > 50000: continue
    g = sns.heatmap(result, cmap='Blues', yticklabels=False, vmin=30, ax=axes[count], cbar=False)
    g.set_title('Reads containing ' + i + ' (n=' + str(len(result)) + ', N=' + str(len(data))+ ')', fontsize =50)
    g.set_xticklabels(g.get_xmajorticklabels(), fontsize = 45)
    count += 1
  plt.tight_layout()
  plt.savefig(outpng)

if __name__ == '__main__' :
  infile = sys.argv[1]
  make_heatmap_several_png (infile)
