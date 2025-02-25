# ArcticFox
Long read sequencing residual AAV analysis pipeline

__version__: v0.4.21.04

__date__: 2022-April for v0.4+

__author__: Chao-Jung Wu <wu.chaojung@gmail.com>


## Dependencies

#### edge.hpc.biogen.com
__Environment__: Login node and scheduler (sbatch to slurm in Conda environment `conda activate compchem`) are supported.

__Packages__:

`anaconda3/4.9.2.compchem`

`samtools/1.12-gcc-11.1.0`

`fastqc/0.11.9-gcc-11.1.0`

`bwa/0.7.17-gcc-11.1.0`

`minimap2/2.14-gcc-11.1.0`

`IGVgui/2.13.1`

The following packages are installed in anaconda3/4.9.2.compchem:

`NanoPlot`, `NanoFilt`

`Python 3.8.6` with libraries `pysam`, `pandas`, `numpy`, `sicpy`, `matplotlib`, `seaborn`

The following setting and package(s) are for compiling report using rmarkdown:

module use /usr/prog/modules/all

module load R/4.2.0-foss-2021b

`R` with libraries `rmarkdown` 

__Access to Biogen sequencing data__:
Be a member to `compbio` and `ngs` groups.
DNANexus directories are under `/edgehpc/dept/compbio/ngs` tree.


#### cambridge.hpc.biogen.com
__Environments__: Login node is supported, except the plasmid module can only generate indexed files but not the igv snapshots. Scheduler (qsub to Sun Grid Engine) is not supported.

__Packages__:

`java/jre/1.8.0_66`

`fastqc/0.11.5`

`pymods/cutadapt/1.18`

`bwa/0.7.15`

`samtools/1.10`

`minimap2/2.17`

The following packages are installed in miniconda3/4.3.21:

`NanoPlot`, `NanoFilt`

`Python 3.6.8` with libraries `pysam`, `pandas`, `numpy`, `sicpy`, `matplotlib`, `seaborn`


## Modules in ArcticFox
preprocessing (n): merge fastq files into fastq.gz, pre-align quality report, remove nanopore adapters. Always active.

genotyping (g): analyze residual DNA based on several pre-defined chromosomes, reporting positivity rate and intact rate. Optional.

lambda analysis (l): analyze lambda fragments and report positivity rate (picking up by template length) and intact rate (reading through by template length). This module requires genotyping module to be activated. Optional.

plasmids (p): analyze over-packaging and facilitate igv viewing on reduced samples. Optional.

cell host (c): analyze residual DNA coming from host cell line and extract reads containing human genome alignments. Optional.


### Note
One must not contain 'sam' in the input path, for example 'sampling'. A filename confusion will cause premature termination of the program.

Cargo name must be exactly "ITR2ITR", or distribution plot may not be optimal.
Vector backbone name must be exactly "backbone" to be considered in distribution analysis.

Lambda analysis is optional, and requires genotyping to be activated.

Must run dos2unix on submit.sh file, in order for the bash file to execute correctly in hpc.

backbone_fasta: fasta file name of cargo vactor minus cargo minus KanR/AmpR, must name it backbone 


## Julie's Debug run 1k (3 mins in login node)
### Must use a new login: exit, login, then run.
### Note that backslash (\\) does not show in markdown, do not copy paste.
### interactive mode
### srun -p interactive --time=0:30:0 --pty bash

exit

cd path/to/ArcticFox/src/



LIBNAME=myTraining

module use /usr/prog/modules/all
module load anaconda3/4.9.2.compchem samtools fastqc bwa IGVgui/2.13.1 

###   -l -g -p -c \
conda activate compchem
time xvfb-run --auto-servernum \
   python pipeline.py \
   -l -g -p -c \
   --hpc edge \
   --libname $LIBNAME \
   --rep_input /home/cwu3/input_home/libsub1k/ \
   --genome_keyword bsg168 \
   --cargo_fasta ITR2ITR_bsg168.fa \
   --backbone_fasta backbone_bsg168.fa \
   --barcode CACAAAGACACCGACAACTTTCTT \
   --plasmids_fasta bsg168.fa

mypath=../output/${LIBNAME}/mapping_plasmid/igv/
BAT=${mypath}my.bat
if test -f "$BAT"; then
	xvfb-run --auto-servernum igv.sh -b $BAT
	cp ${mypath}snapshot_dir/*.png ../output/${LIBNAME}/results/
fi


## Julie's Demo run 10k (6 mins in login node )
### Must use a new login.
### Note that backslash (\\) does not show in markdown, do not copy paste.
### igv snapshot is not generated correclty for larger dataset. Need to run an extra igv line, such as: xvfb-run --auto-servernum igv.sh -b bat, to simplfy, I will run this extra line anyway, regardless the size of dataset.
exit

cd path/to/ArcticFox/src/

LIBNAME=Demo_Edge_ArcticFox_10klib_0719_1452

module load anaconda3/4.9.2.compchem samtools fastqc bwa IGVgui/2.13.1 

conda activate compchem

time xvfb-run --auto-servernum \
   python pipeline.py \
   -l -g -p -c \
   --libname $LIBNAME \
   --rep_input /home/cwu3/input_home/libsub10k/ \
   --genome_keyword bsg168 \
   --cargo_fasta ITR2ITR_bsg168.fa \
   --backbone_fasta backbone_bsg168.fa \
   --barcode CACAAAGACACCGACAACTTTCTT \
   --plasmids_fasta bsg168.fa

mypath=../output/${LIBNAME}/mapping_plasmid/igv/
BAT=${mypath}my.bat
if test -f "$BAT"; then
	xvfb-run --auto-servernum igv.sh -b $BAT
	cp ${mypath}snapshot_dir/*.png ../output/${LIBNAME}/results/
fi

  
## Tasks to do
- [x] debug: haplotype_breakdown.py line 45: NOTE THAT if I change the architecture of genome, the definition of lambda will change too.
- [x] add read70plus in addition to read90plus in report
- [x] testsuite: change ref genome (biogenAAVplatform_bsg168) to be without hg38, facilitate code development.
      refgenome: bsg168_stdchr
- [x] flexible input cargo and backbone fasta
- [x] abandoned --> synch: sync haplotype naming to h, remove g
- [x] droph: remove h haplotype col from file=xxx/yyy.delAdapter5p3p.bwamem.read_haplotype.adj.tsv
- [x] develop plasmid module
- [x] debug getHaploFasta
- [x] automate cat_fasta
- [x] rename cisplasmid to plasmids, because I want to analyze cis- and trans-, not just cis.
- [x] igv input files must be sampled to 10K.
- [x] BUG gzip: ../output/TST11802_bc01_sub10k_2206031527/cutadapt/TST11802_bc01_sub10k_2206031527.delAdapter5p3p.10k.fastq.gz already exists; do you wish to overwrite (y or n)?
    "gzip -c" solves.
- [x] before sampling, make sure dataset is larger than 10000. If not larger, do not sample, it causes bug.
- [x] df.describe() debug
- [ ] test suite camhpc, partially done. If run directly in login node, the table of positivityRate(%) and intactRate(%) is generated correctly. However, if run using qsub, the table is not generated. Two of such tables are missing, one for genotyping, one for lambda fragments. In report.txt, there should be 7 tables, but using qsub, there are only 4 tables. The third missing table is the genotype distribution with adjustment by lambda.
- [x] keep or remove preloads inside pipeline.py? Abandon scheduler in camhpc.
      (1) if preload, still need to call by this: module load miniconda3; python pipeline.py, and the outputs are correct.
	  (2) if no preload for all, there is problem. Should test one by one to locate which proload is required. 
- [x] preprocessing (n) is always active, to avoid bug. It has little impact on non-nanopore datasets.
- [x] successfully migrate to Edge. Can run in login and sheduler with correct outputs.
- [x] See how many reads have adapter if default is longer AAAAAAAAAA (A x 15). Ax3: 6 percent reads have adpaters, and Ax15: 3 percent have adapters. Using nanopore demo dataset to evaluate.
- [x] integrate igv snapshots
- [x] 1) test if I can run a sudo demo snapshot within pipeline --> yes
- [x] 2) test if I can write a dynamic bat within pipeline --> yes
- [x] 3) test if I can execute igv bat and obtain snapshot dymanically
- [x] development note: with login node, the following WORKS. for libsub10k, output png=91kb.
module load samtools fastqc bwa anaconda3/4.9.2.compchem
conda activate compchem
time xvfb-run --auto-servernum python pipeline.py -p --hpc edge ...
Meanwhile inside pipeline.py for calling igv:
module load anaconda3/4.9.2.compchem IGVgui/2.13.1;igv.sh -b bat

With submit.sh
Inside submit, the extra line works.
xvfb-run --auto-servernum igv.sh -b bat
- [x] Edge submit, cell host module did not finish normally. Might be because it is after igv snapshot. It is fine in login node so the problem only exists in submit. --> mem insufficient, increas from 5G to 30G, it works.
- [x] add author name in "every" .py file


- [ ] how to add author name in "every" .sh file


- [ ] spikeins_polynomial_spline.png is drawn in a cell of 3x3 grid. Should be in a 1x1 grid.



- [x] cjwaav.html is broken. --> remove html function all together
- [ ] consider rmarkdown instead

- [ ] merge hg reads into one file
- [ ] If infile is already pre-processed from PacBio and other platforms, evaluate how much the preprocessing impacts the results.
- [ ] Cargo name must be exactly "ITR2ITR". To handle this during automation of ref genome making.



- [ ]  mass percentage has not been reported properly. Need to format that!
- [ ]  They want to make the gentyping chromsomes name and content flexible. Let me see how to do it.











- [ ] make a new public demo run. Get 1000 reads for the demo input library.
time python longreadAAVpipelinecjw.py \
  -g -l -p -c \
  --libname test_100reads_220616_1507 \
  --rep_input ../input_100reads/ \
  --barcode TAGGGAAACACGATAGAATCCGAA \
  --genome_keyword pX602 \
  --cargo_fasta 01b_ITR2ITR_PX602.fa \
  --backbone_fasta backbone_bsg168.fa \
  --plasmids_fasta pX602_AAV_TBG.fa


###################
# 220622 generate new outputs based on v0.4.21
###################

## sample01_nanopore
module load miniconda3
#qsub -V -b y -cwd -l h_rt=02:00:00 -m e -N ngl.ArcticFox.testsuite -M chaojung.wu@biogen.com \
time python longreadAAVpipelinecjw.py \
  -n -g -l -p \
  --libname TST11802_bc01_220622_2123 \
  --rep_input /home/cwu3/input_home/220406_TST11802/barcode01/ \
  --genome_keyword bsg168 \
  --cargo_fasta ITR2ITR_bsg168.fa \
  --backbone_fasta backbone_bsg168.fa \
  --barcode CACAAAGACACCGACAACTTTCTT \
  --plasmids_fasta bsg168.fa
  
## sample04_nanopore 
module load miniconda3
time python longreadAAVpipelinecjw.py \
  -n -g -l -p \
  --libname TST11802_bc04_220622_2123 \
  --rep_input /home/cwu3/input_home/220406_TST11802/barcode04/ \
  --genome_keyword bsg168 \
  --cargo_fasta ITR2ITR_bsg168.fa \
  --backbone_fasta backbone_bsg168.fa \
  --barcode TAGGGAAACACGATAGAATCCGAA \
  --plasmids_fasta bsg168.fa
  
## sample05_nanopore 
module load miniconda3
time python longreadAAVpipelinecjw.py \
  -n -g -l -p \
  --libname TST11802_bc05_220622_2123 \
  --rep_input /home/cwu3/input_home/220406_TST11802/barcode05/ \
  --genome_keyword bsg133 \
  --cargo_fasta ITR2ITR_bsg168.fa \
  --backbone_fasta backbone_bsg133.fa \
  --barcode AAGGTTACACAAACCCTGGACAAG \
  --plasmids_fasta bsg133.fa


## sample01_pacbio
module load miniconda3
time python longreadAAVpipelinecjw.py \
  -n -g -l -p \
  --libname TST11801_ds01_220622_2123 \
  --rep_input /home/cwu3/input_home/220512_TST11801/bc1001/ \
  --genome_keyword bsg168 \
  --cargo_fasta ITR2ITR_bsg168.fa \
  --backbone_fasta backbone_bsg168.fa \
  --barcode CACATATCAGAGTGCGT \
  --plasmids_fasta bsg168.fa
 
## sample04_pacbio
module load miniconda3
time python longreadAAVpipelinecjw.py \
  -n -g -l -p \
  --libname TST11801_ds04_220622_2123 \
  --rep_input /home/cwu3/input_home/220512_TST11801/bc1008/ \
  --genome_keyword bsg168 \
  --cargo_fasta ITR2ITR_bsg168.fa \
  --backbone_fasta backbone_bsg168.fa \
  --barcode ACAGTCGAGCGCTGCGT \
  --plasmids_fasta bsg168.fa
  
## sample05_pacbio
module load miniconda3
time python longreadAAVpipelinecjw.py \
  -n -g -l -p \
  --libname TST11801_ds05_220622_2123 \
  --rep_input /home/cwu3/input_home/220512_TST11801/bc1009/ \
  --genome_keyword bsg133 \
  --cargo_fasta ITR2ITR_bsg168.fa \
  --backbone_fasta backbone_bsg133.fa \
  --barcode ACACACGCGAGACAGAT \
  --plasmids_fasta bsg133.fa
  
