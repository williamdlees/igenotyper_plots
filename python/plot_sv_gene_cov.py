#!/bin/env python
# Create IGenotyper reports and statistics (code extracted from IGenotyper)

import argparse
import os
from pybedtools import BedTool
import numpy as np

from string import Template

from plotutils import plot_histogram
from readutils import execute, snps_from_reads, bam_to_bigwig, hap_bam_to_bigwig, read_lengths, read_quality, target_region_coverage, num_target_reads, total_reads, feat_coverage, primary_alignments, get_phased_regions, add_haplotype_to_blocks


parser = argparse.ArgumentParser(description='Create IGenotyper reports and statistics')
parser.add_argument('sample_name', help='Sample name')
parser.add_argument('input_bam', help='Input bam file (ccs.bam)')
parser.add_argument('phased_reads', help='Phased reads file (ccs_to_ref_phased_sorted.bam)')
parser.add_argument('phased_snps', help='SNPs vsf (phased_snps.vcf)')
parser.add_argument('phased_blocks', help='Phased blocks (phased_blocks.txt)')
parser.add_argument('gene_coords', help='Gene coordinates (gene_coords.bed)')
parser.add_argument('sv_coords', help='SV coords (sv_coords.bed)')
parser.add_argument('target_regions', help='Target regions (target_regions.bed)')
args = parser.parse_args()


def get_igh_region():
    igh_chrom = "igh"
    igh_start = None
    igh_end = None
    target_regions = BedTool(args.target_regions)
    for target_region in target_regions:
        if str(target_region.chrom) == igh_chrom:
            igh_start = int(target_region.start)
            igh_end = int(target_region.end)
    assert igh_start is not None
    assert igh_end is not None
    return (igh_chrom, igh_start, igh_end)


def coverage_stats(bamfn):
    bam = BedTool(bamfn)
    bedgraph = bam.genomecov(bg=True)
    targets = BedTool(args.target_regions)
    bedgraph = bedgraph.intersect(targets)
    coverage = []
    is_igh_present = False  # remove; Added so that IG can work for Rhesus

    for row in bedgraph:
        if row[0] == "igh":
            is_igh_present = True
            start = int(row[1])
            end = int(row[2])
            cov = int(row[3])
            for _ in range(start, end):
                coverage.append(cov)

    if is_igh_present:   # remove if statement only; Added so that IG can work for Rhesus
        igh_chrom, igh_start, igh_end = get_igh_region()
        igh = BedTool([(igh_chrom, igh_start, igh_end)])
        igh = igh.subtract(bedgraph)
        for row in igh:
            for _ in range(row.start, row.end):
                coverage.append(cov)

    mean = round(np.mean(coverage), 2)
    std = round(np.std(coverage), 2)
    cv = round(std/mean, 2)
    return (mean, std, cv)


def input_stats():
    ccs_read_lengths = read_lengths(args.input_bam)
    ccs_read_length_avg = sum(ccs_read_lengths)/len(ccs_read_lengths)
    plot_histogram(ccs_read_lengths, "CCS read lengths", 'ccs_read_lengths.png')

    ccs_quals = read_quality(args.input_bam)
    ccs_quals_avg = sum(ccs_quals)/len(ccs_quals)
    plot_histogram(ccs_quals, "CCS quality scores", 'ccs_read_quals.png')

    regions_covered = target_region_coverage('primary.bam', args.target_regions)
    igh_bases = None

    for chrom, start, end, num_bases in regions_covered:
        if chrom == "igh":
            igh_bases = num_bases

    target_read_count = num_target_reads('primary.bam')
    all_read_count = total_reads(args.input_bam)
    mean_cov, std_cov, cv_cov = coverage_stats('primary.bam')
    stats = {
        "num_subreads": total_reads(args.input_bam),
        "num_ccs": all_read_count,
        "on_target_count": target_read_count,
        "percent_on_target": round((float(target_read_count)/all_read_count), 2) * 100,
        "ref_bases_ccs_coverage": igh_bases,
        "ccs_read_length_avg": ccs_read_length_avg,
        "ccs_quals_avg": round(ccs_quals_avg, 2),
        "mean_coverage": mean_cov,
        "std_coverage": std_cov,
        "coeff_of_var": cv_cov
    }
    return stats


def phased_snps():
    num_phased_snps = {}
    num_unphased_snps = {}
    snps = snps_from_reads(args.phased_snps)
    for chrom in snps:
        if chrom not in num_phased_snps:
            num_phased_snps[chrom] = 0
        if chrom not in num_unphased_snps:
            num_unphased_snps[chrom] = 0
        for pos in snps[chrom]:
            if "/" in snps[chrom][pos]:
                num_unphased_snps[chrom] += 1
            if "|" in snps[chrom][pos]:
                num_phased_snps[chrom] += 1
    return (num_phased_snps, num_unphased_snps)


def get_phased_blocks():
    phased_blocks = []
    target_regions = BedTool(args.target_regions)
    phased_regions = target_regions.intersect(BedTool(get_phased_regions(args.phased_blocks)))
    unphased_regions = target_regions.subtract(phased_regions)

    for haplotype in ["1", "2"]:
        phased_blocks = add_haplotype_to_blocks(phased_blocks, phased_regions, haplotype)

    phased_blocks = add_haplotype_to_blocks(phased_blocks, unphased_regions, "0")
    return phased_blocks


def min_coverage(feature, cov=10):
    return int(feature.name) >= 10


def covered_regions(bam):
    a = BedTool(bam)
    b = a.genomecov(bg=True)
    min_cov = b.filter(min_coverage)
    return min_cov


def phased_bases_per_chrom():
    phased_bases = {}
    phased_regions = BedTool(get_phased_blocks())
    regions_with_coverage = covered_regions('primary.bam')
    phased_regions_with_cov = phased_regions.intersect(regions_with_coverage)
    for region in phased_regions_with_cov:
        hap = region.name
        if hap == "2":
            continue
        if hap == "0":
            continue
        chrom = region.chrom
        start = int(region.start)
        end = int(region.end)
        if chrom not in phased_bases:
            phased_bases[chrom] = 0
        phased_bases[chrom] += (end - start)
    return phased_bases


def phased_stats(stats):
    num_phased_snps, num_unphased_snps = phased_snps()
    phased_bases = phased_bases_per_chrom()

    if "igh" in phased_bases:
        stats["phased_bases"] = phased_bases["igh"]
        stats["num_phased_snps"] = num_phased_snps["igh"]
        stats["num_unphased_snps"] = num_unphased_snps["igh"]
    else:
        stats["phased_bases"] = 0
        stats["num_phased_snps"] = None
        stats["num_unphased_snps"] = None

    gene_cov = feat_coverage(args.gene_coords, 'primary.bam')
    gene_cov.sort(key=lambda x: x[-1])

    with open('gene_cov.txt', 'w') as fh:
        for row in gene_cov:
            fh.write("%s\n" % "\t".join(map(str, row)))

    execute("Rscript ../scripts/rplot_gene_cov.R gene_cov.txt gene_cov.png sv_gene_cov.png")
    stats["num_genes_no_cov"] = len([i for i in gene_cov if i[-1] == 0])
    return stats


def write_to_bashfile(template, outfile, params):
    with open(template, 'r') as fi, open(outfile, 'w') as fh:
        src = Template(fi.read())
        fh.write(src.safe_substitute(params))


def make_track_panel():
    # make haplo bigwigs
    if not os.path.isfile('hapall.bw'):
        bam_to_bigwig(args.phased_reads, 'hapall.bw')

    if not os.path.isfile('hapall.bw'):
        print('hapall.bw not created - quitting')
        exit(1)

    # Make phased snps
    if not os.path.isfile('phased_snps.bed') or not os.path.isfile('unphased_snps.bed'):
        snps = snps_from_reads(args.phased_snps)

        with open('phased_snps.bed', 'w') as fp, open('unphased_snps.bed', 'w') as fu:
            for chrom in snps:
                for pos in snps[chrom]:
                    gt = snps[chrom][pos]
                    out = [chrom, pos, pos+1]
                    if gt in ["1/1"]:
                        fu.write("%s\n" % "\t".join(map(str, out)))
                    if gt in ["0|1", "1|0"]:
                        fp.write("%s\n" % "\t".join(map(str, out)))

    # Make phased_blocks.bed
    if not os.path.isfile('phased_blocks.bed'):
        with open('phased_blocks.bed', 'w') as fh:
            with open(args.phased_blocks, 'r') as fofh:
                for line in fofh:
                    line = line.rstrip().split('\t')
                    if "#" in line[0]:
                        continue
                    if line[3] == line[4]:
                        continue
                    region = [line[1]] + line[3:]
                    fh.write("%s\n" % "\t".join(map(str, region)))

    # Make haplo bigwigs
    for i in ["0", "1", "2"]:
        hap_bw = f"hap{i}.bw"

        if not os.path.isfile(hap_bw):
            hap_bam_to_bigwig(args.phased_reads, i, hap_bw)

        if not os.path.isfile(hap_bw):
            print(f'{hap_bw} not created - quitting')
            exit(1)

    # Copy other required files
    if not os.path.isfile('gene_coords.bed'):
        execute(f'cp {args.gene_coords} gene_coords.bed')

    if not os.path.isfile('sv_coords.bed'):
        execute(f'cp {args.sv_coords} sv_coords.bed')

    # Run pygenometracks
    if not os.path.isfile('sv_gene_cov.png'):
        execute("pyGenomeTracks --tracks ../templates/track_panels.ini --region igh:1-1193129 -o phasing.png")


def phasing_stats():
    if not os.path.isfile('primary.bam'):
        primary_alignments(args.phased_reads, 'primary.bam')

    stats = input_stats()
    stats = phased_stats(stats)
    stats["sample"] = args.sample_name

    write_to_bashfile("../templates/report.html", 'report.html', stats)
    write_to_bashfile("../templates/stats.json", 'stats.json', stats)


make_track_panel()
phasing_stats()
