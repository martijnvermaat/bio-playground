#!/usr/bin/env python
"""
Write BAM file coverage info to a WIG file and BED file.

Run with no arguments for usage info. The script requires pysam [1] and is
partly inspired by [2].

Todo: Implement the window_size argument.
Todo: Use default filenames for coverage_file and summary_file based on
      bam_file.

[1] http://code.google.com/p/pysam/
[2] https://github.com/chapmanb/bcbb/blob/master/nextgen/scripts/bam_to_wiggle.py

Copyright (c) 2011 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2011 Martijn Vermaat <m.vermaat.hg@lumc.nl>
"""


from __future__ import division

import os
from contextlib import contextmanager
from itertools import repeat

import argparse
import pysam


def main(bam_file, coverage_file, summary_file, regions_file=None,
         window_size=1):
    #if not coverage_file:
    #    coverage_file = '%s.wig' % os.path.splitext(bam_file)[0]
    with open(coverage_file, 'w') as coverage:
        regions = write_coverage(bam_file, coverage, regions_file)
    with open(summary_file, 'w') as summary:
        write_summary(regions, summary)


@contextmanager
def indexed_bam(bam_file):
    if not os.path.exists(bam_file + '.bai'):
        pysam.index(bam_file)
    bam = pysam.Samfile(bam_file, 'rb')
    yield bam
    bam.close()


def read_regions(regions_file):
    with open(regions_file) as regions:
        for region in regions:
            if region.startswith('track'):
                continue
            name, start, end = region.strip().split('\t')[:3]
            yield name, int(start), int(end)


def write_coverage(bam_file, coverage, regions_file=None):
    coverage.write('track %s\n' % ' '.join(['type=wiggle_0',
        'name=%s' % os.path.splitext(os.path.split(bam_file)[-1])[0],
        'visibility=full']))
    with indexed_bam(bam_file) as bam:
        regions = []
        if regions_file is not None:
            guide = read_regions(regions_file)
        else:
            guide = zip(bam.references, repeat(0), bam.lengths)
        for name, start, end in guide:
            if not regions or name != regions[-1][0]:
                coverage.write('variableStep chrom=%s\n' % name)
            summed_coverage = 0
            for column in bam.pileup(name, start, end):
                if start <= column.pos < end:
                    summed_coverage += column.n
                    coverage.write('%s %.1f\n' % (column.pos + 1, column.n))
            regions.append( (name, start, end, summed_coverage / (end - start)) )
    return regions


def write_summary(regions, summary):
    for name, start, end, average in regions:
        summary.write('%s\t%i\t%i\t-\t%i\n' % (name, start, end, average))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    group = parser.add_argument_group()
    group.add_argument('bam_file', metavar='BAM_FILE',
                       help='file in BAM format to determine coverage for')
    group.add_argument('-c', dest='coverage_file', required=True,
                       help='write coverage in WIG format')
    group.add_argument('-s', dest='summary_file', required=True,
                       help='write coverage per region in BED format')
    parser.add_argument('-r', dest='regions_file',
                        help='regions to calculate coverage for in BED format')
    parser.add_argument('-w', dest='window_size', default=1, type=int,
                        help='window size for COVERAGE_FILE (default: 1)')
    args = parser.parse_args()
    main(args.bam_file, args.coverage_file, args.summary_file,
         args.regions_file, args.window_size)
