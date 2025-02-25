'''
program: arg.py
author: Chao-Jung Wu <wu.chaojung@gmail.com>
date: 2020-11-03
updated: 2021-04-26
'''
import argparse


def str2bool(v):
  if isinstance(v, bool): 
    return v
  if v.lower() in ('yes', 'true', 't', 'y', '1'): 
    return True
  elif v.lower() in ('no', 'false', 'f', 'n', '0'): 
    return False
  else: 
    raise argparse.ArgumentTypeError('Boolean value expected.')



def init_params ():
  parser = argparse.ArgumentParser()

  parser.add_argument('--hpc', default='edge', choices=['cambridge', 'edge'])
  #parser.add_argument('-n', type=str2bool, nargs='?', const=True, default=False, help='Activate Nanopore Preprocessing Analysis')  
  parser.add_argument('-l', type=str2bool, nargs='?', const=True, default=False, help='Activate Lambda Analysis')
  parser.add_argument('-g', type=str2bool, nargs='?', const=True, default=False, help='Activate Gene Analysis')
  parser.add_argument('-p', type=str2bool, nargs='?', const=True, default=False, help='Activate Plasmid Analysis')
  parser.add_argument('-c', type=str2bool, nargs='?', const=True, default=False, help='Activate host cell Human Genome Analysis')
  
  parser.add_argument('--libname', default='myLibrary', help='default=myLibrary')
  parser.add_argument('--rep_input', help='/path/to/input/folder/')
  parser.add_argument('--aligner', choices=['minimap2', 'bwamem'], default='bwamem', help='default=bwamem')
  parser.add_argument('--barcode', default='AAAAAAAAAAAAAAA')

  parser.add_argument('--genome_keyword', default='universalAAVplatform', help='Provide a keyword for the standardized genome for genotyping with customize ITR2ITR and backbone chromosomes, included chromosomes: Cargo, AmpR, KanR, Rep2, Cap9, HPV18, Ad5, Lambda, plasmidBackbone')
  parser.add_argument('--cargo_fasta', default='01b_ITR2ITR_PX602.fa', help='fasta file name of cargo ITR2ITR, must name it ITR2ITR (ie. the first line in the file is ">ITR2ITR"), must place this file in ../dbs/')
  parser.add_argument('--backbone_fasta', default='10_F1pUCoris.fa', help='fasta file name of cargo vactor minus cargo minus KanR/AmpR, must name it backbone (ie. the first line in the file is ">backbone"), must place this file in ../dbs/')
  parser.add_argument('--plasmids_fasta', default='pX602_AAV_TBG.fa', help='fasta file name of all cis- and trans- plasmids, must place this file in ../dbs/')
  
  parser.add_argument('--AdapterYtop', default='GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT', \
                         help='default=GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT')
  parser.add_argument('--AdapterYbottom', default='GCAATACGTAACTGAACGAAGT', help='default=GCAATACGTAACTGAACGAAGT') 
  parser.add_argument('--linker5p', default='AAGGTTAA', help='default=AAGGTTAA')
  parser.add_argument('--linker3p', default='TTAACCTTA', help='default=TTAACCTTA')
  parser.add_argument('--trim_poorSequencing_Nanopore', type=str2bool, nargs='?', const=True, default=False)
  ##parser.add_argument('--trimAdapters', type=str2bool, nargs='?', const=True, default=False)
  parser.add_argument('--trim_poorAlignment', type=str2bool, nargs='?', const=True, default=False)
  parser.add_argument('--min_LambdaCt', type=int, default=10)
  #parser.add_argument('--min_lenRead', type=int, default=200)
  #parser.add_argument('--min_lenAlnChr', type=int, default=200)
  
  parser.add_argument('--getHaploFasta', type=str2bool, nargs='?', const=True, default=False, help='extract reads of each haplotype in fasta format') 
  return parser.parse_args()



if __name__ == '__main__' :
  args = init_params ()
  #print(args)
  #print(args.linker3p)

  for arg in sorted(vars(args)):
    print (arg + ':\t\t\t', getattr(args, arg))
