# Plots and analysis extracted from IGenotyper  

The original code is from https://github.com/oscarlr/IGenotyper  

## Usage  

usage: plot_sv_gene_cov.py [-h] sample_name input_bam phased_reads phased_snps phased_blocks gene_coords sv_coords target_regions  

Create IGenotyper reports and statistics  

positional arguments:  
  sample_name     - Sample name  
  input_bam       - Input bam file (ccs.bam)  
  phased_reads    - Phased reads file (ccs_to_ref_phased_sorted.bam)  
  phased_snps     - SNPs annotated with phase (1/1, 0|1, 1|0) (phased_snps.vcf)  
  phased_blocks   - Phased blocks (phased_blocks.txt)  
  gene_coords     - Gene coordinates (gene_coords.bed)  
  sv_coords       - SV coords (sv_coords.bed)  
  target_regions  - Target regions (target_regions.bed)  

options:  
  -h, --help      - show this help message and exit  

Various files will be created in the current directory - see below for notes on input and output files.  

## Notes on input files

input_bam - BAM file of reads produced by the PacBio sequencer.  
phased_reads - BAM file of reads phased by IGenotyper, aligned against the reference sequence.  
phased_snps - VCF file of phased SNPs in the reads, compared to the reference sequence. By default, IGenotyper produces this via 'whatshap phase'.  
phased_blocks - phased block list in TSV format, as produced by 'whatshapp stats --block-list'.  
gene_coords - BED file of core coding regions of V(D)J genes.  
sv_coords - BED file of structural variants included in the reference. Presence or absence of these in the sample will be highlighted in the plots.  
target_regions - BED file of capture-probe target regions. Coverage of these will be highlighted in the plots.  

## Output files (in the example directory)

primary.bam - BAM file of reads that are primary alignments in input_bam.  
stats.json - selected stats (calculated stats are pushed into this file and report.html via templates)  

### Plots

report.html - composite report bringing together the plots and statistics  
ccs_read_lengths.png - histogram of read lengths in input_bam  
ccs_read_quals.png - histogram of CCS quality scores in input_bam  
gene_cov.png - read coverage in primary.bam of all genes in gene_coords  
sv_gene_cov.png - read coverage in primary.bam of all genes in SVs  
phasing.png - genome tracks showing phasing, coverage, SVs etc  

### Intermediate files (not intended for further use)

hap 0,1,2 .bam - haplotyped sets of primary reads, used to produce bigwig coverage plots  
phased_snps.bed - SNPs explicitly phased to 0|1 or 1|0, used for coverage plots  
unphased_snps.bed - SNPs not explicitly phased to 0|1 or 1|0 (1/1), used for coverage plots  
hap 0,1,2,all .bw - coverage plots  
gene_cov.txt - file produced for rplot_gene_cov.R listing gene coverage by haplotype  
gene_coords.bed, sv_coords.bed - copied to the working directory for use by pygenometracks  
Rplots.pdf - intermediate file produced by rplot_gene_cov.R  

### Notes

rplot_gene_cov.R has an explicit list of igh SVs and the genes they contain  
There are a number of things that are specific to igh in the plots and statistics - searching for 'igh' finds these, I think  
