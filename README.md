# ArcticFox
Long read sequencing residual AAV analysis pipeline

__version__: v0.4.21.04

__update__: 2022-April

__author__: Chao-Jung Wu <wu.chaojung@gmail.com>


## Dependencies

__Packages__:

`anaconda3`

`samtools/1.12-gcc-11.1.0`

`fastqc/0.11.9-gcc-11.1.0`

`bwa/0.7.17-gcc-11.1.0`

`minimap2/2.14-gcc-11.1.0`

`IGVgui/2.13.1`

`NanoPlot`, `NanoFilt`

`Python 3.6.8` with libraries `pysam`, `pandas`, `numpy`, `sicpy`, `matplotlib`, `seaborn`


## Modules in ArcticFox
preprocessing (n): merge fastq files into fastq.gz, pre-align quality report, remove nanopore adapters. Always active.

genotyping (g): analyze residual DNA based on several pre-defined chromosomes, reporting positivity rate and intact rate. Optional.

lambda analysis (l): analyze lambda fragments and report positivity rate (picking up by template length) and intact rate (reading through by template length). This module requires genotyping module to be activated. Optional.

plasmids (p): analyze over-packaging and facilitate igv viewing on reduced samples. Optional.

cell host (c): analyze residual DNA coming from host cell line and extract reads containing human genome alignments. Optional.


### Note
Cargo name must be exactly "ITR2ITR", or distribution plot may not be optimal.

Vector backbone name must be exactly "backbone" to be considered in distribution analysis.

Lambda analysis is optional, and requires genotyping to be activated.

backbone_fasta: fasta file name of cargo vactor minus cargo minus KanR/AmpR, must name it backbone 