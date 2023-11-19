# Read utilities from IGenotyper

import subprocess
from collections import namedtuple

import pysam
from pybedtools import BedTool


def execute(cmd):
    print(f"\n-- executing {cmd}\n")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(result.stdout)
    print(result.stderr)


def snps_from_reads(phased_snps):
    snps = {}
    with open(phased_snps, 'r') as vcf_fh:
        for line in vcf_fh:
            if line.startswith("#"):
                continue
            line = line.rstrip().split("\t")
            chrom = line[0]
            position = line[1]
            genotype = line[9].split(":")[0]
            if chrom not in snps:
                snps[chrom] = {}
            snps[chrom][int(position) - 1] = genotype
    return snps


def bam_to_bigwig(bam, bigwig):
    execute(f'bamCoverage -b {bam} -o {bigwig}')


def select_hap_sequence(bam, hap, outbam):
    execute(f'samtools view -Sbh -F 3884 {bam} -r {hap} > {outbam}')
    execute(f'samtools index {outbam}')


def hap_bam_to_bigwig(bam, hap, bigwig):
    outbam = f"hap{hap}.bam"
    select_hap_sequence(bam, hap, outbam)
    bam_to_bigwig(outbam, bigwig)


def primary_alignments(inbam, outbam):
    execute(f'samtools view -Sbh -F 3884 {inbam} > {outbam}')
    execute(f'samtools index {outbam}')


def read_lengths(bam_file):
    lengths = []
    samfile = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
    for read in samfile:
        lengths.append(read.query_length)
    return lengths


def read_quality(bam_file):
    quality = []
    samfile = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
    for read in samfile:
        score = sum(map(float, read.query_qualities))/read.query_length
        quality.append(score)
    return quality


def load_bed_regions(bedfile, add_fourth=False):
    bed_regions = []
    with open(bedfile, 'r') as bedfh:
        for line in bedfh:
            line = line.rstrip().split('\t')
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            if add_fourth:
                annotation = True
                if len(line) == 4:
                    annotation = line[3]
                bed_regions.append([chrom, start, end, annotation])
            else:
                bed_regions.append([chrom, start, end])
    return bed_regions


def total_reads(bam):
    sam = pysam.AlignmentFile(bam, check_sq=False)
    # assert int(sam.mapped) == 0      # TODO - (WL) reinstate when we have a real input file to test with
    reads = int(sam.unmapped) if sam.unmapped else sam.mapped   # TODO - (WL) remove if statement when we have a real input file to test with
    return reads


def skip_read(read):
    skip = False
    if read.is_unmapped:
        skip = True
    if read.is_supplementary:
        skip = True
    if read.is_secondary:
        skip = True
    return skip


def num_target_reads(bam):
    sam = pysam.AlignmentFile(bam)
    count = 0
    for read in sam:
        if skip_read(read):
            continue
        count += 1
    return count


def target_region_coverage(primary_alignment_bam, target_bed, min_cov=10):
    regions = []
    sam = pysam.AlignmentFile(primary_alignment_bam)
    capture_regions = load_bed_regions(target_bed)
    for chrom, start, end in capture_regions:
        bases = 0
        for pileupcolumn in sam.pileup(chrom, start, end):
            if pileupcolumn.n >= min_cov:
                bases += 1
        regions.append([chrom, start, end, bases])
    return regions


def feat_coverage(features, bam):
    sam = pysam.AlignmentFile(bam)
    features = BedTool(features)
    coverage = []
    for row in features:
        chrom = str(row[0])
        start = int(row[1])
        end = int(row[2])
        feature = str(row[3])
        row_coverage = [chrom, start, end, feature, 0, 0, 0, 0]
        for read in sam.fetch(reference=chrom, start=start, end=end):
            if skip_read(read):
                continue
            if read.reference_start > start:
                continue
            if read.reference_end < end:
                continue
            read_hap = int(read.get_tag("RG", True)[0])
            row_coverage[4] += 1
            row_coverage[5 + read_hap] += 1
        coverage.append(row_coverage)
    return coverage


def get_phased_regions(phased_blocks, min_length=500, min_variants=2):
    blocks = []
    Block = namedtuple('Block', ['sample', 'chrom', 'start_1', 'start', 'end', 'num_variants'])
    with open(phased_blocks, 'r') as fh:
        _ = fh.readline()
        for line in fh:
            line = line.rstrip().split('\t')
            block = Block._make(line)
            if int(block.num_variants) < min_variants:
                continue
            if (int(block.end) - int(block.start)) < min_length:
                continue
            blocks.append([block.chrom, int(block.start), int(block.end)])
    return sorted(blocks, key=lambda x: x[1])


def add_haplotype_to_blocks(phased_blocks, regions, haplotype):
    for region in regions:
        block = [
            region.chrom,
            region.start,
            region.end,
            haplotype
            ]
        phased_blocks.append(block)
    return phased_blocks
