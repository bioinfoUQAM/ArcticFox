'''
author: Chao-Jung Wu <wu.chaojung@gmail.com>
'''
import os, re, sys
import matplotlib.pyplot as plt 
import pandas as pd

def plot_multiple_hist2d (infile):
  outpng = re.sub('.aln_coors.tsv', '.aln_coors.hist2d.png', infile)
  df = pd.read_csv(infile, sep='\t', low_memory=False)
  #Cols = df.columns.values.tolist(); print('infile cols: ', Cols) #= ['readID', 'readLength', 'q_start', 'q_end', 'q_aln_len', 'chrName', 's_start', 's_end', 's_aln_len']
  chrNames = sorted(list(set([x for x in df['chrName'] if not x.startswith('chr')])))
  for i in 'upmapped,hg38'.split(','):#ITR2ITR_PX602 
    if i in chrNames: chrNames.remove(i)
  #= print the sliced data to a tmp file, this step is necessary to avoid error msg.
  tmp = re.sub('.aln_coors.tsv', '.aln_coors.cleaned.tsv', infile) #  tmp = 'cleaned.tsv'
  result = df[ df['chrName'].isin(chrNames) ]
  result.to_csv(tmp, sep='\t', index=False)
  #= reading coors from the cleaned file
  df = pd.read_csv(tmp, sep='\t', low_memory=False)
  os.system('rm ' + tmp)
  if len(df) < 2000 and 'Rep2' in chrNames: 
    chrNames.remove('Rep2') ##= problematic?
    print('line24: Rep2 data were removed to avoid error.')
  #= draw hist2d plots for each chromosome in one figure
  print('Drawing hist2d subplots for each chromosome in one figure.\nchromosomes:', ','.join(chrNames))
  fig, axarr = plt.subplots(nrows=3, ncols=3, figsize=(24, 16))
  fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.1, hspace=0.3)
  fig.suptitle("\n".join([infile, 'Chromosome (s-) start and end alignment positions (hist2d)']), fontsize=25, y=0.98)
  #fig.set_tight_layout(True) #= this is for mac, but subplot hist2d does not work
  axes = [axarr[0][0], axarr[0][1], axarr[0][2], 
          axarr[1][0], axarr[1][1], axarr[1][2], 
          axarr[2][0], axarr[2][1], axarr[2][2]]
  count = 0
  for i in chrNames:
    result = df[ df['chrName'] == i ]
    ax, x, y = axes[count], 's_start', 's_end'; count += 1
    ax.set_title(i + ' (n=' + str(len(result)) + ', N=' + str(len(df))+ ')', fontsize =20)
    ax.set_xlabel(x)  
    ax.set_ylabel(y)
    _, _, _, im = ax.hist2d(result[x], result[y], bins=20, cmap=plt.cm.Blues) #=counts, xedges, yedges, im
    fig.colorbar(im, ax=ax)
  plt.savefig(outpng)
  print('outfigure:', outpng)


#infile = '../myData/SMN1_new.delAdapter5p3p.bwamem.aln_coors.tsv'

if __name__ == '__main__' :
  coors = sys.argv[1]
  plot_multiple_hist2d (coors)
